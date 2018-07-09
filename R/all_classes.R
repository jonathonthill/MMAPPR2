
#' MmapprParam object
#' 
#' \code{MmapprParam} stores options for running MMAPPR
#'
#' @slot refGenome GmapGenome. 
#' @slot wtFiles BamFileList. 
#' @slot mutFiles BamFileList. 
#' @slot vepFlags VEPFlags. 
#' @slot distancePower numeric. 
#' @slot peakIntervalWidth numeric. 
#' @slot minDepth numeric. 
#' @slot homozygoteCutoff numeric. 
#' @slot minBaseQuality numeric. 
#' @slot minMapQuality numeric. 
#' @slot loessOptResolution numeric. 
#' @slot loessOptCutFactor numeric. 
#' @slot naCutoff numeric. 
#' @slot outputFolder character. 
#' @slot fileAggregation A length-one character vector determining strategy for
#'   aggregating base calls when multiple wild-type or multiple mutant files are provided.
#'   When 'sum', average base call proportions are computed using
#'   the read counts from
#'   each file, effectively weighting files
#'   with higher counts at a given position. When equal to 'mean', the
#'   base call proportions as well as read depths, rather than the absolute count,
#'   are averaged across files, which is useful when you want to weight each
#'   replicate evenly without
#'   regards to differing depth.
#'
#' @return
#' @export
#'
#' @examples
setClass("MmapprParam",
         representation(
             refGenome = "GmapGenome",
             wtFiles = "BamFileList",
             mutFiles = "BamFileList",
             vepFlags = "VEPFlags",
             distancePower = "numeric",
             peakIntervalWidth = "numeric",
             minDepth = "numeric",
             homozygoteCutoff = "numeric", #the maximum WT allele frequency we'll accept for candidates
             minBaseQuality = "numeric",
             minMapQuality = "numeric",
             #resolution at which AICc will be calculated to find optimum Loess fit span
             loessOptResolution = "numeric",
             #factor between rounds of Loess fit optimization (e.g., factor of 0.1 results in spans of 0.1 apart, then 0.01 apart, etc.)
             loessOptCutFactor = "numeric",
             naCutoff = "numeric", # the most NAs we'll accept, that is, the number of files without data for that position
             outputFolder = "character",
             fileAggregation='character'
         )
)



#' MmapprData object
#' 
#' Stores data from each step of the MMAPPR pipeline.
#'
#' @slot param MmapprParam. TODO
#' @slot distance list. 
#' @slot peaks list. 
#' @slot candidates list. 
#'
#' @return
#' @export
#'
#' @examples
setClass("MmapprData",
         representation(
             param="MmapprParam",
             distance="list",
             peaks="list",
             candidates="list"
         )
)

