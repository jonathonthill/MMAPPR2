
#' MmapprParam Class
#' 
#' \code{MmapprParam} stores parameters running \code{\link{mmappr}}
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
#' @slot fileAggregation character.
#'
#' @rdname MmapprParam
#' @export
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
#' @export
setClass("MmapprData",
         representation(
             param="MmapprParam",
             distance="list",
             peaks="list",
             candidates="list"
         )
)

