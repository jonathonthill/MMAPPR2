#' Read BAM files and generate Euclidean distance data
#'
#' First step in the MMAPPR2 pipeline. Precedes the \code{\link{loessFit}}
#' step.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{distance}
#'   slot filled.
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
#'
#'     md <- new('MmapprData', param = mmappr_param)
#'     postCalcDistMD <- calculateDistance(md)
#' }

calculateDistance <- function(mmapprData) {
  .indexBamFileList(wtFiles(param(mmapprData)), 
                    oF = outputFolder(param(mmapprData)))
  .indexBamFileList(mutFiles(param(mmapprData)), 
                    oF = outputFolder(param(mmapprData)))
  
  chrList <- suppressWarnings(.getFileReadChrList(mmapprData))
  
  mmapprData@distance <-
    BiocParallel::bplapply(chrList, .calcDistForChr, param = mmapprData@param)
  
  return(mmapprData)
}

.indexBamFileList <- function(bfl, oF) {
  for (i in seq_along(bfl)) {
    bamFile <- bfl[[i]]
    if (is.na(Rsamtools::index(bamFile))) {
      .messageAndLog("No index found. Indexing Bam now.", oF)
      Rsamtools::indexBam(bamFile)
    } else {
      .messageAndLog("Index found. Skipping index step.", oF)
    }
  }
}


.calcDistForChr <- function(chrRange, param){
  tryCatch({
    #parameter check
    stopifnot(length(GenomeInfoDb::seqnames(chrRange)) == 1)
    stopifnot(is(param, "MmapprParam"))
    
    stopifnot(!is.null(chrRange))
    
    #apply functions to wild type pool
    #CAUTION: functions must be applied in this order to work right
    tryCatch({
      pileupWT <- .samtoolsPileup(files = wtFiles(param), param, chrRange)       # Read and format files
      pileupWT <- .avgFiles(pileupWT, fileAggregation = fileAggregation(param))  # Average files

    }, error = function(e) {
      msg <- 'Insufficient data in wild-type file(s)'
      print(e)
      stop(msg)
    }
    )
    
    tryCatch({
      pileupMut <- .samtoolsPileup(files = mutFiles(param), param, chrRange)
      pileupMut <- .avgFiles(pileupMut, fileAggregation = fileAggregation(param))
    }, error = function(e) {
      msg <- 'Insufficient data in mutant file(s)'
      stop(msg)
    }
    )
    
    #inner_join (also removes rows without a match)
    setkey(pileupWT, CHROM, POS, REF, ALT)
    setkey(pileupMut, CHROM, POS, REF, ALT)
    distanceDf <- merge(pileupWT, pileupMut, suffixes = c(".WT", ".MT"))
    
    #filter out uninformative snps (both homozygous and identical)
    distanceDf <- distanceDf[!((AVE_REF_FREQ.WT > homozygoteCutoff(param) &
                                 AVE_REF_FREQ.MT > homozygoteCutoff(param)) |
                               (AVE_ALT_FREQ.WT > homozygoteCutoff(param) &
                                 AVE_ALT_FREQ.MT > homozygoteCutoff(param)))]
    
    if (length(distanceDf) == 0)
      stop('Empty dataframe after joining WT and Mut count tables')
    
    #calculate Euclidian distance
    distanceDf <- distanceDf[, DISTANCE := sqrt((AVE_REF_FREQ.WT - AVE_REF_FREQ.MT)^2 + 
                                                  (AVE_ALT_FREQ.WT - AVE_ALT_FREQ.MT)^2) ^ 
                               distancePower(param)]
    
    stopifnot(nrow(distanceDf) > 0)
    
    resultList <- list(wtCounts = pileupWT, mutCounts = pileupMut,
                       distanceDf = distanceDf)
    resultList$seqname <- as.character(GenomeInfoDb::seqnames(chrRange))
    
    .messageAndLog(paste("Finished with ", chrRange), 
                   outputFolder = outputFolder(param))
    return(resultList)
  },
  
  error = function(e) {
    msg <- paste(toString(GenomeInfoDb::seqnames(chrRange)),
                 e$message,
                 sep = ': ')
    return(msg)
  }
  )
}

.getFileReadChrList <- function(mmapprData) {suppressWarnings({
  bams <- Rsamtools::BamFileList(c(mmapprData@param@wtFiles,
                                   mmapprData@param@mutFiles))
  
  bamInfo <- Rsamtools::seqinfo(bams)
  chrRanges <- as(bamInfo, "GRanges")
  #cut to standard chromosomes
  chrRanges <- GenomeInfoDb::keepStandardChromosomes(chrRanges,
                                                     pruning.mode = 'coarse')
  chrRanges <- GenomeInfoDb::dropSeqlevels(chrRanges, 'chrM',
                                           pruning.mode = 'coarse')
  chrRanges <- GenomeInfoDb::dropSeqlevels(chrRanges, 'MT',
                                           pruning.mode = 'coarse')
  
  chrList <- list()
  # store range for each chromosome as list item
  for (i in suppressWarnings(
    GenomeInfoDb::orderSeqlevels(
      as.character(GenomeInfoDb::seqnames(chrRanges))))) {
    
    chrList[[toString(GenomeInfoDb::seqnames(chrRanges[i]))]] <-
      chrRanges[i]
  }
  
  
  return(chrList)
})}

###Get averages between files and divides by coverage
# function takes data.table with pos, nuc (ACGTcvg--with counts),
# file1, file2...returns it with files combined into one mean column
# Can weight by cvg in each file or simple average frequencies
.avgFiles <- function(chrDf, fileAggregation) {
  stopifnot(fileAggregation %in% c('simple', 'weighted'))
  setkey(chrDf, CHROM, POS, REF, ALT)
  #throw away non-mean columns and return
  if (fileAggregation == 'weighted') {
    chrDf <- chrDf[, .(AVE_CVG = mean(CVG), 
                       AVE_REF_FREQ = sum((REF_FREQ*CVG))/sum(CVG), 
                       AVE_ALT_FREQ = sum(ALT_FREQ*CVG)/sum(CVG)), 
                     .(CHROM, POS, REF, ALT)]
  }
  else if (fileAggregation == 'simple') {
    chrDf <- chrDf[, .(AVE_CVG = mean(CVG), 
                       AVE_REF_FREQ = mean(REF_FREQ), 
                       AVE_ALT_FREQ = mean(ALT_FREQ)), 
                     .(CHROM, POS, REF, ALT)]
  }
  return(chrDf)
}