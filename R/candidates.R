#' @title Generate candidate mutations and consequences in peak regions
#'
#' @name generateCandidates
#'
#' @usage Follows the \code{\link{peakRefinement}} step and produces a
#' \code{\linkS4class{md}} object ready for
#' \code{\link{outputmd}}.
#'
#' @param md The \code{\linkS4class{md}} object to be analyzed.
#'
#' @return A \code{\linkS4class{md}} object with the \code{candidates}
#'   slot filled with a \code{\link[GenomicRanges]{GRanges}} object for each
#'   peak chromosome containing variants and predicted consequences from
#'   Ensembl's Variant Effect Predictor.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#' 
#' \dontrun{
#' md <- mmapprData(mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#'
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' }
#' 
NULL

generateCandidates <- function(md) {
  #get GRanges representation of peak
  .messageAndLog("Getting Variants in Peak", 
                 outputFolder(param(md)))
  peakGRanges <- lapply(md@peaks, .getPeakRange)
  
  #call variants in peak
  md@candidates$snps <- lapply(peakGRanges,
                               FUN=.getVariantsForRange,
                               param=md@param)
  
  #predict effects of variants
  .messageAndLog("Predicting Variant Effects", 
                 outputFolder(param(md)))
  md@candidates$effects <- lapply(md@candidates$snps,
                                  FUN=.predictEffects,
                                  param=md@param)
  
  #add diff expressed genes
  .messageAndLog("Identifying Differentially Expressed Genes in Peak", 
                 outputFolder(param(md)))
  md@candidates$diff <- lapply(peakGRanges,
                               FUN=.addDiff,
                               param=md@param)
  
  #density score and order variants
  .messageAndLog("Ordering Candidates by Position", 
                 outputFolder(param(md)))
  md@candidates <- .scoreVariants(md@candidates, md@peaks)
  
  return(md)
}


.getPeakRange <- function(peakList) {
  ir <- IRanges::IRanges(start=as.numeric(peakList$start),
                         end=as.numeric(peakList$end),
                         names=peakList$seqname)
  
  gr <- GRanges(seqnames=names(ir),
                ranges=ir)
  return(gr)
}


.getVariantsForRange <- function(inputRange, param) {
  # merge files in desired region if there are multiple
  mergedBam <- file.path(outputFolder(param), 'merged.tmp.bam')
  if (length(param@mutFiles) < 2) mutBam <- param@mutFiles[[1]]
  else{
    mutBam <- mergeBam(param@mutFiles,
                       destination=mergedBam,
                       region=inputRange)
  }
  
  # create param for variant calling
  tallyParam <- TallyVariantsParam(genome=param@refGenome,
                                   which=inputRange,
                                   indels=TRUE)
  
  resultVr <- callVariants(mutBam, tally.param=tallyParam)
  resultVr <- resultVr[altDepth(resultVr)/totalDepth(resultVr) > 0.8]
  if (file.exists(mergedBam)) file.remove(mergedBam)
  
  if (length(resultVr) > 0) {
    # need sampleNames to convert to VCF; using mutant file names
    Biobase::sampleNames(resultVr) <-
      paste0(names(param@mutFiles),
             collapse = " -- ")
    S4Vectors::mcols(resultVr) <- NULL
    return(resultVr)
  }
  else return(NULL)
}

.predictEffects <- function(inputVariants, param){
  suppressMessages(txdb <- makeTxDbFromGFF(gtf(param), 
                                           chrominfo = seqinfo(inputVariants)))    # This prevents incompatibility between the seqinfo for the txdb and inputVariants objects
  effects <- predictCoding(query = inputVariants,
                           subject = txdb,
                           seqSource = FaFile(refFasta(param)),
                           varAllele = DNAStringSet(alt(inputVariants)))
}

.addDiff <- function(peakGRange, param) {
  #prep data
  suppressMessages(genes <- fread(cmd=paste("gzcat", gtf(param)), 
                                  showProgress = FALSE))
  genes <- genes[V3 == "gene"
  ][, gene_id := gsub(".*gene_id \"(.*?)\";.*", "\\1", V9)
  ][, gene_name := gsub(".*gene_name \"(.*?)\";.*", "\\1", V9)
  ][, .(seqnames = V1, start = V4, end = V5, strand = V7, 
        gene_id, gene_name)]
  genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
  genes <- subsetByOverlaps(x=genes, ranges = peakGRange)
  
  #get and process counts
  readfiles <- c(param@wtFiles, param@mutFiles)
  counts <- GenomicAlignments::summarizeOverlaps(features=genes,
                                                 reads=readfiles,
                                                 param=ScanBamParam(which = peakGRange))
  num_wt <- length(param@wtFiles)
  countDF <- as.data.table(assays(counts)$counts)
  countDF[, ave_wt := rowMeans(countDF[, 1:(num_wt + 1)])    
  ][, ave_mt := rowMeans(countDF[, (num_wt + 1):length(countDF)])
  ][, log2FC := round(log2(ave_mt/ave_wt), 3)]
  
  #merge counts with genes
  mcols(genes) <- cbind(mcols(genes), countDF)
  
  #filter for differentially expressed genes
  return(genes[(abs(genes$log2FC) > 1 | is.na(genes$log2FC)) 
               & (genes$ave_wt > 10 | genes$ave_mt > 10)])
}


.scoreVariants <- function(candList, peaks) {
  for (GRname in names(candList)) {
    GR <- candList[[GRname]]
    for (seqname in names(GR)) {
      print(seqname)
      #density calculation
      densityFunc <- peaks[[seqname]]$densityFunction
      stopifnot(!is.null(densityFunc))
      positions <- BiocGenerics::start(GR[[seqname]]) +
        ((BiocGenerics::width(GR[[seqname]]) - 1) / 2)
      densityCol <- vapply(positions, densityFunc, FUN.VALUE=numeric(1))
      mcols(candList[[GRname]][[seqname]])$peakDensity <- densityCol
      
      #re-order
      candList[[GRname]][[seqname]] <-
        .orderVariants(candList[[GRname]][[seqname]])
    }
  }
  return(candList)
}


.orderVariants <- function(candidateGRanges) {
  if(!(class(candidateGRanges) %in% c("VRanges", "GRanges"))) { 
    .messageAndLog("invalid data type for sorting. (orderVariants)", 
                   outputFolder(param(md)))
    return(candidateGRanges)
  }
  if(!is.null(candidateGRanges$CONSEQUENCE)) {
    impactLevels <- c("synonymous", "nonsynonymous", "frameshift", "nonsense")
    orderVec <- order(match(candidateGRanges$CONSEQUENCE, impactLevels), 
                      candidateGRanges$peakDensity, decreasing=TRUE)
  } else {
    orderVec <- order(candidateGRanges$peakDensity, decreasing=TRUE)
  }
  candidateGRanges <- candidateGRanges[orderVec]
  return(candidateGRanges)
}

