#' Mutation Mapping Analysis Pipeline for Pooled RNA-seq
#' 
#' The main functionality of this package is described in the \code{\link{mmappr}}
#' function.
#' 
#' @inherit mmappr examples
#' @docType package
#' @name MMAPPR2
#' 
NULL

#' @import methods
#' @importFrom magrittr %>%
#' @importFrom graphics abline mtext par plot polygon
#' @importFrom grDevices dev.off pdf rgb
#' @importFrom stats approxfun density.default loess median
#' @importFrom utils capture.output object.size str write.table
#' @importFrom gmapR GmapGenome
#' @importFrom ensemblVEP VEPFlags
#' @importFrom Rsamtools BamFileList
NULL
