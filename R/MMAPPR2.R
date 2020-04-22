#' Mutation Mapping Analysis Pipeline for Pooled RNA-seq
#'
#' The main functionality of this package is described in the \code{\link{mmappr}}
#' function.
#'
#' @docType package
#' @name MMAPPR2
#'
NULL

#' @import BiocParallel
#' @import data.table
#' @import GenomeInfoDb
#' @import Rsamtools
#' @rawNamespace import(GenomicRanges, except = c(shift))
#' @importFrom graphics abline mtext par plot polygon axis legend
#' @importFrom grDevices dev.off dev.list pdf rgb graphics.off
#' @importFrom ensemblVEP VEPFlags ensemblVEP flags
#' @importFrom VariantTools TallyVariantsParam callVariants
#' @importFrom VariantAnnotation altDepth totalDepth writeVcf
#' @importFrom SummarizedExperiment assays
#' @importFrom gmapR GmapGenome
#' @importFrom stats approxfun density.default loess median var
#' @importFrom utils capture.output object.size str write.table tail sessionInfo
NULL
