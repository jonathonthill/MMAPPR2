#' @title MmapprData Class
#' 
#' @name MmapprData-Class
#' @description
#' Stores data from each step of the MMAPPR2 pipeline.
#'
#' @slot param \code{\linkS4class{MmapprParam}} object storing parameters
#'   used in analysis.
#' @slot snpDistance List containing raw counts and Euclidean distance data for
#'   each chromosome. After \code{\link{calculateDistance}}, chromosomes with
#'   sufficient data should have \code{$wtCounts},
#'   \code{$mutCounts}, and \code{$distanceDf} populated. After
#'   \code{\link{loessFit}}, the \code{$distanceDf} element for each chromosome
#'   list is replaced with a \code{$loess} element.
#' @slot peaks List of chromosomes containing peak regions. Initialized after
#'   \code{\link{prePeak}} and populated with density function after
#'   \code{\link{peakRefinement}}.
#' @slot candidates List containing \code{\link[GenomicRanges]{GRanges}} objects
#'   with candidate variants, locally predicted coding effects, and descriptive expression data
#'   for each peak, all resulting from \code{\link{generateCandidates}} 
#'   function.
#'
#' @aliases MmapprData
#' @seealso \code{\link{mmappr}}, \link{MmapprData-functions}
#' @include param.R
#' 
NULL

setClass("MmapprData",
         representation(
           param="MmapprParam",
           snpDistance="list",
           peaks="list",
           candidates="list"
         )
)

#' @title MmapprData Constructor
#' 
#' @name mmapprData
#' @description
#' Creates a \code{MmapprData} class Object. This object stores data from each
#' step of the MMAPPR2 pipeline.
#'
#' @param param \code{\linkS4class{MmapprParam}} object storing parameters
#'   used in analysis.
#'   
#' @return A \code{MmapprData} object.
#' @export
#'
#' @seealso \code{\linkS4class{MmapprData}} \code{\linkS4class{MmapprParam}}
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#'     mmapprData(mmappr_param)
#' }
#' 
NULL


mmapprData <- function(param) {
  new("MmapprData", param = param)
}


#' @title MmapprData Getters and Setters
#' @description
#' Access slots of \code{\linkS4class{MmapprData}} object. The only slot with 
#' a setter is "param," as the others must be generated using one of the 
#' functions in the pipeline.
#'
#' @param obj Desired \code{\linkS4class{MmapprData}} object.
#'
#' @return Desired attribute.
#'
#' @name MmapprData-functions
#' @aliases param snpDistance peaks candidates
#' @seealso \code{\linkS4class{MmapprData}}
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#'
#'     md <- mmapprData(mmappr_param)
#'
#'     param(md)
#'     snpDistance(md)
#'     peaks(md)
#'     candidates(md)
#' }
#' 
NULL


#' @rdname MmapprData-functions
#' @export
setMethod("param", "MmapprData", function(obj) obj@param)
#' @rdname MmapprData-functions
#' @export
setMethod("snpDistance", "MmapprData", function(obj) obj@snpDistance)
#' @rdname MmapprData-functions
#' @export
setMethod("peaks", "MmapprData", function(obj) obj@peaks)
#' @rdname MmapprData-functions
#' @export
setMethod("candidates", "MmapprData", function(obj) obj@candidates)
#' @rdname MmapprData-functions
#' @export
setMethod("param<-", "MmapprData",
          function(obj, value) {
            if (!is(value, "MmapprParam")) stop("param must be a MmapprParam object")
            # Every derived slot depends on the analysis parameters.
            # Replacing param while retaining old distances/peaks/candidates creates
            # a silently inconsistent object, so invalidate downstream results.
            hadResults <- length(obj@snpDistance) > 0L || length(obj@peaks) > 0L ||
                          length(obj@candidates) > 0L
            obj@param <- value
            if (hadResults) {
              warning("Replacing param invalidates existing derived results; snpDistance, peaks, and candidates were cleared")
              obj@snpDistance <- list()
              obj@peaks <- list()
              obj@candidates <- list()
            }
            obj
          })


setMethod("show", "MmapprData", function(object) {
  margin <- "  "
  cat("MmapprData object with following slots:\n")
  cat("param:\n")
  .customPrint(object@param, margin)
  
  cat("snpDistance:\n")
  classes <- vapply(object@snpDistance, class, character(1))
  successes <- classes == "list"
  cat(margin, sprintf(
    "Contains Euclidean distance data for %i sequence(s)\n",
    sum(successes)), sep = "")
  loessFits <- 0
  try({loessFits <- sum(vapply(object@snpDistance[successes],
                               FUN.VALUE = logical(1),
                               FUN = function(seq) {
                                 if (!is.null(seq$loess)) return(TRUE)
                                 else return(FALSE)
                               }
  ))
  }, silent = TRUE
  )
  cat(margin, sprintf(
    "and Loess regression data for %i of those\n", loessFits
  ), sep = "")
  distanceSize <- object.size(object@snpDistance)
  cat(margin, sprintf("Memory usage = %.0f MB\n",
                      distanceSize/1000000), sep = "")
  
  cat("peaks:\n")
  for (peak in object@peaks){
    cat(margin, sprintf('%s: ', peak$seqname), sep = '')
    # this will fail and skip if start and end aren't calculated
    cat(sprintf('start = %.0f, end = %.0f',
                peak$start, peak$end), sep = "")
    cat('\n')
    if (!is.null(peak$densityFunction))
      cat(margin, margin, "Density function calculated\n", sep = "")
  }
  
  cat("candidates:\n")
  for (i in seq_along(object@candidates))
    .customPrint(object@candidates[i], margin, lineMax = 5)
})


