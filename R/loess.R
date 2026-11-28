#' Perform optimized Loess regression for each chromosome
#'
#' Called after the \code{\link{calculateDistance}} step and before
#' \code{\link{prePeak}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{$loess} 
#'   element of the \code{snpDistance} slot list filled.
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#' 
#' \dontrun{
#' md <- mmapprData(mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#'
#' postLoessMD <- loessFit(postCalcDistMD)
#' }
#' 
#' @export

loessFit <- function(mmapprData) {
    loessOptResolution <- mmapprData@param@loessOptResolution
    loessOptCutFactor <- mmapprData@param@loessOptCutFactor

    # each item (chr) of distance list has mutCounts, wtCounts,
    # distanceDf going in
    mmapprData@snpDistance <-
        BiocParallel::bplapply(mmapprData@snpDistance,
                 FUN = .loessFitForChr,   
                 loessOptResolution = loessOptResolution,
                 loessOptCutFactor = loessOptCutFactor)

    return(mmapprData)
}


.getLoess <- function(s, pos, eucDist, ...){
    suppressWarnings(try(stats::loess(eucDist ~ pos,
                                      span = s,
                                      degree = 1,
                                      family = "symmetric"),
                         silent = TRUE))
}


.localResolution <- function(spans, span) {
    spans <- unique(sort(spans[is.finite(spans)]))
    index <- match(span, spans)
    if (is.na(index)) return(Inf)
    neighbors <- numeric()
    if (index > 1L) neighbors <- c(neighbors, span - spans[index - 1L])
    if (index < length(spans)) neighbors <- c(neighbors, spans[index + 1L] - span)
    if (length(neighbors) == 0L) Inf else max(neighbors)
}


.aiccOpt <- function(distanceDf, spans, resolution, cutFactor) {
    evaluated <- data.frame(spans = numeric(), aiccValues = numeric())
    frontier <- unique(round(spans[spans > 0 & spans <= 1],
                             digits = .numDecimals(resolution)))

    converged <- FALSE
    for (iter in seq_len(50L)) {
        already <- evaluated$spans
        frontier <- setdiff(frontier, already)
        if (length(frontier) > 0L) {
            vals <- vapply(frontier, .aicc, FUN.VALUE = numeric(1),
                           eucDist = distanceDf$DISTANCE,
                           pos = distanceDf$POS)
            evaluated <- rbind(evaluated,
                               data.frame(spans = frontier, aiccValues = vals))
            evaluated <- evaluated[!duplicated(evaluated$spans), , drop = FALSE]
        }

        finite <- evaluated[is.finite(evaluated$aiccValues), , drop = FALSE]
        if (nrow(finite) == 0L)
            stop("All candidate LOESS/AICc fits failed for this chromosome")
        finite <- finite[order(finite$spans), , drop = FALSE]

        minimaIdx <- .localMinIndices(finite$aiccValues)
        if (length(minimaIdx) == 0L) minimaIdx <- which.min(finite$aiccValues)
        minimaIdx <- minimaIdx[order(finite$aiccValues[minimaIdx])]
        minimaIdx <- head(minimaIdx, 2L)
        minSpans <- finite$spans[minimaIdx]

        newSpans <- numeric()
        allSpans <- sort(unique(evaluated$spans))
        for (minSpan in minSpans) {
            localRes <- .localResolution(allSpans, minSpan)
            if (is.finite(localRes) && localRes > resolution) {
                step <- localRes * cutFactor
                addVector <- step * seq_len(9L)
                proposed <- c(minSpan - addVector, minSpan + addVector)
                proposed <- proposed[proposed > 0 & proposed <= 1]
                proposed <- round(proposed, digits = .numDecimals(resolution))
                newSpans <- c(newSpans, proposed)
            }
        }
        newSpans <- setdiff(unique(newSpans), evaluated$spans)
        if (length(newSpans) == 0L) {
            converged <- TRUE
            break
        }
        frontier <- newSpans
    }

    # Do not silently pretend the coarse-to-fine search converged if it
    # hit its defensive iteration cap. Returning the evaluated table is still useful,
    # but the warning makes an unexpectedly difficult objective visible to the user.
    if (!converged)
        warning("LOESS span optimization reached its 50-iteration safety cap before convergence")

    evaluated[order(evaluated$spans), , drop = FALSE]
}


.localMinIndices <- function(x) {
    if (!length(x)) return(integer())
    out <- integer()
    for (i in seq_along(x)) {
        if (!is.finite(x[i])) next
        left <- if (i == 1L) Inf else x[i - 1L]
        right <- if (i == length(x)) Inf else x[i + 1L]
        if ((!is.finite(left) || x[i] <= left) &&
            (!is.finite(right) || x[i] <= right)) out <- c(out, i)
    }
    out
}


.numDecimals <- function(x) {
    stopifnot(is.numeric(x), length(x) == 1L, is.finite(x), x > 0)
    txt <- format(x, scientific = FALSE, trim = TRUE, digits = 15)
    txt <- sub("0+$", "", txt)
    if (!grepl("\\.", txt)) 0L else nchar(sub("^.*\\.", "", txt))
}


.aicc <- function(s, eucDist, pos) {
    x <- .getLoess(s, pos, eucDist)
    if (inherits(x, "try-error")) return(NA_real_)
    n <- x$n
    traceL <- x$trace.hat
    denom <- n - traceL - 2
    if (!is.finite(denom) || denom <= 0 || n <= 1L) return(NA_real_)
    sigma2 <- sum(x$residuals^2, na.rm = TRUE) / (n - 1)
    if (!is.finite(sigma2) || sigma2 <= 0) return(NA_real_)
    value <- log(sigma2) + 1 + 2 * (2 * (traceL + 1)) / denom
    if (is.finite(value)) value else NA_real_
}


.chooseBestAiccSpan <- function(aiccTable, resolution) {
    finite <- aiccTable[is.finite(aiccTable$aiccValues), , drop = FALSE]
    if (nrow(finite) == 0L) stop("No finite AICc values were produced")
    bestRows <- finite$aiccValues == min(finite$aiccValues)
    tiedSpans <- sort(unique(finite$spans[bestRows]))
    tiedMedian <- stats::median(tiedSpans)
    # Return an ACTUALLY EVALUATED optimum. Averaging tied spans can create
    # an unevaluated value whose AICc is unknown. Ties equidistant from the median
    # resolve to the smaller evaluated span for deterministic behavior.
    distances <- abs(tiedSpans - tiedMedian)
    minDistance <- min(distances)
    # Binary floating-point can make mathematically symmetric
    # ties (for example 0.1 and 0.3 around 0.2) differ by ~1e-17. which.min()
    # could therefore choose the larger span despite the documented smaller-span
    # tie rule. Treat machine-scale differences as equal, then choose min().
    tol <- sqrt(.Machine$double.eps) * max(1, abs(tiedMedian), minDistance)
    tiedClosest <- tiedSpans[abs(distances - minDistance) <= tol]
    bestSpan <- min(tiedClosest)
    round(bestSpan, digits = .numDecimals(resolution))
}


.loessFitForChr <- function(resultList, loessOptResolution, loessOptCutFactor){
    startTime <- proc.time()
    tryCatch({
        if (is(resultList, "character")) stop("-- Loess fit failed: ", resultList)
        if (is.null(resultList$distanceDf) || nrow(resultList$distanceDf) < 5L)
            stop("Too few informative positions for LOESS fitting")

        startSpans <- c(seq(.01, .16, .01), seq(.21, .91, .10))
        resultList$aicc <- .aiccOpt(distanceDf = resultList$distanceDf,
                                    spans = startSpans,
                                    resolution = loessOptResolution,
                                    cutFactor = loessOptCutFactor)

        bestSpan <- .chooseBestAiccSpan(resultList$aicc, loessOptResolution)
        resultList$bestSpan <- bestSpan

        resultList$loess <- .getLoess(bestSpan,
                                      resultList$distanceDf$POS,
                                      resultList$distanceDf$DISTANCE)
        if (inherits(resultList$loess, "try-error") ||
            !inherits(resultList$loess, "loess"))
            stop("Final LOESS fit failed at optimized span ", bestSpan)

        resultList$distanceDf <- NULL
        resultList$seqname <- NULL
        resultList$loessTime <- proc.time() - startTime
        resultList
    }, error = function(e) {
        if (is(resultList, "character")) paste0(resultList, ": ", e$message) else e$message
    })
}
