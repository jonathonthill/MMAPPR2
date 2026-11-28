#' MMAPPR2: Mutation Mapping Analysis Pipeline for Pooled RNA-seq
#'
#' The main functionality of this package is described in the \code{\link{mmappr}}
#' function.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".", ".I", "..bases", "..wtFreqCols", ".row_id",
  "A", "C", "G", "T", "A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ",
  "AVE.A.FREQ.MT", "AVE.A.FREQ.WT", "AVE.C.FREQ.MT", "AVE.C.FREQ.WT",
  "AVE.G.FREQ.MT", "AVE.G.FREQ.WT", "AVE.T.FREQ.MT", "AVE.T.FREQ.WT",
  "AVE.CVG", "CHROM", "CVG", "DISTANCE", "FILE_ID", "POS",
  "attributes", "count", "end", "max_end", "nucleotide", "pos", "ref",
  "refDepth", "alt", "altDepth", "seqname", "seqnames", "start", "strand",
  "totalDepth", "type"
))

#' @import BiocParallel
#' @import data.table
#' @import GenomeInfoDb
#' @import Rsamtools
#' @import txdbmaker
#' @importFrom Biostrings DNAStringSet
#' @rawNamespace import(GenomicRanges, except = c(shift))
#' @importFrom graphics abline mtext par plot polygon axis legend
#' @importFrom grDevices dev.off dev.list pdf
#' @importFrom VariantAnnotation altDepth totalDepth refDepth predictCoding alt ref VRanges
#' @importFrom SummarizedExperiment assays
#' 
#' @importFrom stats approxfun density loess median sd var
#' @importFrom utils capture.output head object.size tail sessionInfo
#' @importFrom IRanges subsetByOverlaps
NULL
