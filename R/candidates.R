#' Generate potential causative mutations and consequences in peak regions
#'
#' Follows the \code{\link{peakRefinement}} step and produces a
#' \code{\linkS4class{MmapprData}} object ready for
#' \code{\link{outputMmapprData}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{candidates}
#'   slot filled with a \code{\link[GenomicRanges]{GRanges}} object for each
#'   peak chromosome containing variants and predicted consequences from
#'   Ensembl's Variant Effect Predictor.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)
#'         & all(Sys.which(c("samtools", "vep")) != "")) {
#'     mmappr_param <- MmapprParam(refFasta = MMAPPR2data::goldenFasta(),
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio",
#'                                outputFolder = tempOutputFolder())
#' }
#' \dontrun{
#' md <- new('MmapprData', param = mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#'
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' }

generateCandidates <- function(mmapprData) {
  
  #get GRanges representation of peak
  peakGRanges <- lapply(mmapprData@peaks, .getPeakRange)
  
  #call variants in peak
  mmapprData@candidates$snps <- lapply(peakGRanges,
                                       FUN=.getVariantsForRange,
                                       param=mmapprData@param)
  
  #run VEP
  mmapprData@candidates$effects <- lapply(mmapprData@candidates$snps,
                                          FUN=.runVEPForVariants,
                                          param=mmapprData@param)
  
  #filter out low impact effects
  mmapprData@candidates$effects <- lapply(mmapprData@candidates$effects, 
                                          FUN=.filterVariants, 
                                          impact = mmapprData@param@vep_impact)
  
  #add diff expressed genes
  mmapprData@candidates$diff <- lapply(peakGRanges,
                                       FUN=.addDiff,
                                       param=mmapprData@param)
  
  #density score and order variants
  mmapprData@candidates <- lapply(mmapprData@candidates,
                                  FUN=.scoreVariants,
                                  mmapprData@peaks)
  
  return(mmapprData)
}


.getPeakRange <- function(peakList) {
  ir <- IRanges::IRanges(start=as.numeric(peakList$start),
                         end=as.numeric(peakList$end),
                         names=peakList$seqname)
  
  gr <- GenomicRanges::GRanges(seqnames=names(ir),
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
                                   indels=TRUE
  )
  
  resultVr <- callVariants(mutBam, tally.param=tallyParam)
  resultVr <- resultVr[VariantAnnotation::altDepth(resultVr)/
                         VariantAnnotation::totalDepth(resultVr) > 0.8]
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


.peakVcf <- function(param) file.path(outputFolder(param), 'peak.vcf')


.runVEPForVariants <- function(inputVariants, param){
  vepFlags <- vepFlags(param)
  stopifnot(is(vepFlags, "VEPFlags"))
  stopifnot(is(inputVariants, 'VRanges'))
  
  vcf <- .peakVcf(param)
  tryCatch({
    VariantAnnotation::writeVcf(inputVariants, vcf)
    resultGRanges <- ensemblVEP::ensemblVEP(vcf, vepFlags)
  }, error=function(e) {
    stop(e)
  },finally={
    # if (file.exists(vcf)) file.remove(vcf)
  })
  
  return(resultGRanges)
}

.filterVariants <- function(candidateGRanges, impact) {
  impact_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  filter <-
    S4Vectors::mcols(candidateGRanges)$IMPACT %in% impact_levels[1:impact]
  filter[is.na(filter)] <- TRUE
  return(candidateGRanges[filter])
}


.addDiff <- function(peakGRange, param) {
  #prep data
  suppressMessages(genes <- data.table::fread(cmd=paste("zcat", param@gtf), 
                                              showProgress = FALSE))
  genes <- genes[V3 == "gene"
               ][, gene_id := gsub(".*gene_id \"(.*?)\";.*", "\\1", V9)
               ][, gene_name := gsub(".*gene_name \"(.*?)\";.*", "\\1", V9)
               ][, .(seqnames = V1, start = V4, end = V5, strand = V7, 
                     gene_id, gene_name)]
  genes <- GenomicRanges::makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
  genes <- IRanges::subsetByOverlaps(x=genes, ranges = peakGRange)
  
  #get and process counts
  readfiles <- c(param@wtFiles, param@mutFiles)
  counts <- GenomicAlignments::summarizeOverlaps(features=genes,
                                                 reads=readfiles,
                                                 param=Rsamtools::ScanBamParam(which = peakGRange))
  num_wt <- length(param@wtFiles)
  countDF <- data.table::as.data.table(
    SummarizedExperiment::assays(counts)$counts)
  countDF[, ave_wt := rowMeans(countDF[, 1:(num_wt + 1)])    
        ][, ave_mt := rowMeans(countDF[, (num_wt + 1):length(countDF)])
        ][, log2FC := round(log2(ave_mt/ave_wt), 3)]
  
  #merge counts with genes
  S4Vectors::mcols(genes) <- cbind(S4Vectors::mcols(genes), countDF)
  
  #filter for differentially expressed genes
  return(genes[(abs(genes$log2FC) > 1 | is.na(genes$log2FC)) 
               & (genes$ave_wt > 10 | genes$ave_mt > 10)])
}


.scoreVariants <- function(candList, peaks) {
  returnData = list()
  for (seqname in names(candList)) {
    #density calculation
    densityFunction <- peaks[[seqname]]$densityFunction
    stopifnot(!is.null(densityFunction))
    positions <- BiocGenerics::start(candList[[seqname]]) +
      ((BiocGenerics::width(candList[[seqname]]) - 1) / 2)
    densityCol <- vapply(positions, densityFunction, FUN.VALUE=numeric(1))
    S4Vectors::mcols(candList[[seqname]])$peakDensity <- densityCol
    
    #re-order
    returnData[[seqname]] <-
      .orderVariants(candList[[seqname]], densityCol)
  }
  return(returnData)
}


.orderVariants <- function(candidateGRanges, densityCol) {
  if(!(class(candidateGRanges) %in% c("VRanges", "GRanges"))) { 
    return(candidateGRanges)
  }
  if(!is.null(candidateGRanges$IMPACT)) {
    impact_levels <- c("MODIFIER", "LOW", "MODERATE", "HIGH")
    orderVec <- order(match(candidateGRanges$IMPACT, impact_levels), 
                      densityCol, decreasing=TRUE)
  } else {
    orderVec <- order(densityCol, decreasing=TRUE)
  }
  candidateGRanges <- candidateGRanges[orderVec]
  return(candidateGRanges)
}

