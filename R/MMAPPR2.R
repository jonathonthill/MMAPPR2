#' Mutation Mapping Analysis Pipeline for Pooled RNA-seq
#' 
#' The main functionality of this package is described in the \code{\link{mmappr}}
#' function.
#' 
#' @docType package
#' @name MMAPPR2
#' 
NULL

#' @import methods
#' @importFrom magrittr %>%
#' @importFrom graphics abline mtext par plot polygon axis legend
#' @importFrom grDevices dev.off pdf rgb graphics.off
#' @importFrom stats approxfun density.default loess median var
#' @importFrom utils capture.output object.size str write.table tail sessionInfo
#' @importFrom gmapR GmapGenome
#' @importFrom ensemblVEP VEPFlags
#' @importFrom Rsamtools BamFileList mergeBam
#' @importFrom VariantTools TallyVariantsParam callVariants
NULL
