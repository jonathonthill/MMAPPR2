#' @title Read BAM files and generate Euclidean distance data
#'
#' @name calculateDistance
#'
#' @usage First step in the MMAPPR2 pipeline. Precedes the \code{\link{loessFit}}
#' step.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{distance}
#'   slot filled.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#'
#'     md <- mmapprData(mmappr_param)
#'     postCalcDistMD <- calculateDistance(md)
#' }
#' 
NULL

calculateDistance <- function(mmapprData) {
  chrList <- suppressWarnings(.getFileReadChrList(param(mmapprData)))
  
  mmapprData@snpDistance <-
    bplapply(chrList, .calcDistForChr, param = mmapprData@param)
  
  return(mmapprData)
}


.getFileReadChrList <- function(param) {suppressWarnings({
  # get data for gtf 
  suppressMessages(gtfData <- fread(cmd=paste("gunzip -c", param@gtf), 
                                    showProgress = FALSE))
  genes <- gtfData[V3 == "gene"
  ][, .(seqnames = V1, seqnames2 = V1, start = V4, end = V5)]
  genes <- makeGRangesListFromDataFrame(genes, 
                                        split.field = "seqnames2")
  
  if (param@includeScaffolds != TRUE) {
    genes <- keepStandardChromosomes(genes,
                                     pruning.mode = 'coarse')
  } 
  
  # remove mitochondiral chromosome b/c not mendelian inheritance
  genes <- dropSeqlevels(genes, 'chrM',
                         pruning.mode = 'coarse')
  genes <- dropSeqlevels(genes, 'MT',
                         pruning.mode = 'coarse')
  
  return(genes)
})}


.calcDistForChr <- function(chrRange, param){
  tryCatch({
    #parameter check
    stopifnot(length(unique(seqnames(chrRange))) == 1)
    stopifnot(is(param, "MmapprParam"))
    
    stopifnot(!is.null(chrRange))
    
    #apply functions to wild type pool
    #CAUTION: functions must be applied in this order to work right
    tryCatch({
      
      pileupWT <- lapply(wtFiles(param), .getPileup, param, chrRange) # Read and format files
      pileupWT <- rbindlist(pileupWT)
      pileupWT <- .avgFiles(pileupWT, 
                            fileAggregation = fileAggregation(param)[1]) # Average files
      pileupWT <- pileupWT[AVE.CVG > param@minDepth]
      
    }, error = function(e) {
      msg <- 'Insufficient data in wild-type file(s)'
      print(e)
      stop(msg)
    }
    )
    
    tryCatch({
      
      pileupMut <- lapply(mutFiles(param), .getPileup, param, chrRange)
      pileupMut <- rbindlist(pileupMut)
      pileupMut <- .avgFiles(pileupMut, 
                             fileAggregation = fileAggregation(param)[1])
      pileupMut <- pileupMut[AVE.CVG > param@minDepth]
      
    }, error = function(e) {
      msg <- 'Insufficient data in mutant file(s)'
      stop(msg)
    }
    )
    
    #inner_join (also removes rows without a match)
    setkey(pileupWT, CHROM, POS)
    setkey(pileupMut, CHROM, POS)
    distanceDf <- merge(pileupWT, pileupMut, suffixes = c(".WT", ".MT"))
    
    if (length(distanceDf) == 0)
      stop('Empty dataframe after joining WT and Mut count tables')
    
    #calculate Euclidian distance
    distanceDf[, DISTANCE := sqrt((AVE.A.FREQ.WT - AVE.A.FREQ.MT)^2 + 
                                  (AVE.C.FREQ.WT - AVE.C.FREQ.MT)^2 +
                                  (AVE.G.FREQ.WT - AVE.G.FREQ.MT)^2 +
                                  (AVE.T.FREQ.WT - AVE.T.FREQ.MT)^2) ^
                                  distancePower(param)]
    
    #filter out uninformative snps (both homozygous and identical)
    homozygousWT <- apply(distanceDf[, AVE.A.FREQ.WT:AVE.T.FREQ.WT], 1, max) 
    homozygousWT <- homozygousWT > homozygoteCutoff(param)
    distanceDf <- distanceDf[!(homozygousWT)]
    
    stopifnot(nrow(distanceDf) > 0)
    
    resultList <- list(wtCounts = pileupWT, mutCounts = pileupMut,
                       distanceDf = distanceDf)
    
    .messageAndLog(paste("Finished", seqnames(chrRange)[1]), 
                   outputFolder = outputFolder(param))
    return(resultList)
  },
  
  error = function(e) {
    msg <- paste(toString(seqnames(chrRange)),
                 e$message,
                 sep = ': ')
    return(msg)
  }
  )
}


# Imports data from BAM files
.getPileup <- function(file, param, chrRange) {
  stopifnot(length(file) == 1)
  
  scanParam <- ScanBamParam(simpleCigar = TRUE,   # should try with FALSE
                            which = chrRange, 
                            mapqFilter=param@minMapQuality)
  
  pParam <- PileupParam(max_depth = 1000,
                        min_mapq = param@minMapQuality, 
                        min_base_quality = param@minBaseQuality,
                        distinguish_strands = FALSE, 
                        include_deletions = FALSE,
                        include_insertions = FALSE)
  
  pData <- data.table(pileup(file, 
                             scanBamParam = scanParam, 
                             pileupParam = pParam))
  
  pData <- dcast(pData, seqnames + pos ~ nucleotide, 
                 value.var = "count", 
                 fun.aggregate = sum)
  
  setnames(pData, colnames(pData), 
           c("CHROM", "POS", "A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ"))
  
  pData[, CVG := A.FREQ + C.FREQ + G.FREQ + T.FREQ
  ][, c("A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ") := 
      list(A.FREQ/CVG, C.FREQ/CVG, G.FREQ/CVG, T.FREQ/CVG)]
  
  return(pData)
}


###Get averages between files and divides by coverage
# function takes data.table with pos, nuc (ACGTcvg--with counts),
# file1, file2...returns it with files combined into one mean column
# Can weight by cvg in each file or simple average frequencies
.avgFiles <- function(chrDf, fileAggregation) {
  stopifnot(fileAggregation %in% c('simple', 'weighted'))
  setkey(chrDf, CHROM, POS)
  #throw away non-mean columns and return
  if (fileAggregation == 'weighted') {
    chrDf <- chrDf[, .(AVE.CVG = mean(CVG), 
                       AVE.A.FREQ = sum(A.FREQ*CVG)/sum(CVG), 
                       AVE.C.FREQ = sum(C.FREQ*CVG)/sum(CVG),
                       AVE.G.FREQ = sum(G.FREQ*CVG)/sum(CVG), 
                       AVE.T.FREQ = sum(T.FREQ*CVG)/sum(CVG)), 
                     .(CHROM, POS)]
  }
  else if (fileAggregation == 'simple') {
    chrDf <- chrDf[, .(AVE.CVG = mean(CVG), 
                       AVE.A.FREQ = mean(A.FREQ), 
                       AVE.C.FREQ = mean(C.FREQ),
                       AVE.G.FREQ = mean(G.FREQ), 
                       AVE.T.FREQ = mean(T.FREQ)), 
                     .(CHROM, POS)]
  }
  return(chrDf)
}