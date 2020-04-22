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
