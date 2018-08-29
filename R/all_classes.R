
#' MmapprParam Class and Constructor
#' 
#' \code{MmapprParam} stores parameters for running \code{\link{mmappr}}.
#'
#' @rdname MmapprParam
#' @export
setClass("MmapprParam",
         representation(
             refGenome = "GmapGenome",
             wtFiles = "BamFileList",
             mutFiles = "BamFileList",
             species = "character",
             vepFlags = "VEPFlags",
             distancePower = "numeric",
             peakIntervalWidth = "numeric",
             minDepth = "integer",
             homozygoteCutoff = "numeric",
             minBaseQuality = "numeric",
             minMapQuality = "numeric",
             loessOptResolution = "numeric",
             loessOptCutFactor = "numeric",
             naCutoff = "integer",
             outputFolder = "character",
             fileAggregation ='character'
         )
)



#' MmapprData Class
#' 
#' Stores data from each step of the MMAPPR2 pipeline.
#'
#' @slot param \code{\linkS4class{MmapprParam}} object storing parameters
#'   used in analysis.
#' @slot distance List containing raw counts and Euclidean distance data for
#'   each chromosome. After \code{\link{calculateDistance}}, chromosomes with
#'   sufficient data should have \code{$wtCounts},
#'   \code{$mutCounts}, and \code{$distanceDf} populated. After
#'   \code{\link{loessFit}}, the \code{$distanceDf} element for each chromosome
#'   list is replaced with a \code{$loess} element.
#' @slot peaks List of chromosomes containing peak regions. Initialized after
#'   \code{\link{prePeak}} and populated with density function after
#'   \code{\link{peakRefinement}}.
#' @slot candidates List containing \code{\link[GenomicRanges]{GRanges}} object
#'   for each peak, resulting from \code{\link{generateCandidates}} function.
#'   VEP data, including gene symbol and consequence for each variant, are
#'   included in metacolumns.
#'
#' @aliases MmapprData
#' @export
#' @seealso \link{MmapprDataGetters}
setClass("MmapprData",
         representation(
             param="MmapprParam",
             distance="list",
             peaks="list",
             candidates="list"
         )
)

