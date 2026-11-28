#' 
#' @title Characterize Euclidean distance peaks using resampling simulation
#'
#' @name peakRefinement
#' @description
#' Follows the \code{\link{prePeak}} step and precedes
#' \code{\link{generateCandidates}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{peaks}
#'   slot filled and populated.
#' @export
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
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#'
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' }
#' 
#' @import BiocParallel
#' 
NULL


peakRefinement <- function(mmapprData){
    mmapprData@peaks <- BiocParallel::bplapply(
        mmapprData@peaks,
        .peakRefinementChr,
        mmapprData = mmapprData
    )
    mmapprData
}


.peakApexPosition <- function(x, fitted) {
    valid <- is.finite(x) & is.finite(fitted)
    if (!any(valid)) return(NA_real_)
    x <- as.numeric(x[valid])
    fitted <- as.numeric(fitted[valid])
    maxFit <- max(fitted)
    # Flat-topped fits should not be biased toward the first/left-most point.
    # The median coordinate is stable and symmetric across a plateau.
    stats::median(x[fitted == maxFit])
}


.getSubsampleLoessMax <- function(rawData, loessSpan) {
    n <- nrow(rawData)
    if (n < 5L) return(NA_real_)
    sampleSize <- max(5L, floor(n * 0.5))
    idx <- sample.int(n, size = sampleSize, replace = FALSE)
    tempData <- rawData[idx, , drop = FALSE]
    tempData <- tempData[order(tempData$pos), , drop = FALSE]

    fit <- suppressWarnings(try(
        stats::loess(euclideanDistance ~ pos,
                     data = tempData,
                     span = loessSpan,
                     degree = 1,
                     family = "symmetric"),
        silent = TRUE
    ))
    # One bad half-sample should not abort all 1,000 resamples.
    if (inherits(fit, "try-error") || !inherits(fit, "loess")) return(NA_real_)
    .peakApexPosition(as.numeric(fit$x), fit$fitted)
}


.densityPointMass <- function(x, y) {
    stopifnot(length(x) == length(y), length(x) >= 1L)
    if (length(x) == 1L) return(1)
    dx <- diff(x)
    if (any(!is.finite(dx)) || any(dx <= 0)) stop("Density x coordinates must be strictly increasing")
    # Approximate the integration cell represented by each KDE grid point. This
    # makes the mass calculation correct even if a future density grid is not
    # perfectly equally spaced.
    widths <- c(dx[1] / 2, (head(dx, -1L) + tail(dx, -1L)) / 2, tail(dx, 1L) / 2)
    rawMass <- pmax(y, 0) * widths
    if (!is.finite(sum(rawMass)) || sum(rawMass) <= 0) stop("Peak-density mass is empty")
    rawMass / sum(rawMass)
}


.getPeakFromTopP <- function(data, topP, method = c("hpd_span", "shortest_contiguous")) {
    stopifnot(ncol(data) == 2L, topP > 0, topP <= 1)
    method <- match.arg(method)
    names(data) <- c("x", "y")
    data <- data[is.finite(data$x) & is.finite(data$y) & data$y >= 0, , drop = FALSE]
    data <- data[order(data$x), , drop = FALSE]
    if (nrow(data) == 0L || sum(data$y) <= 0) stop("Peak-density data are empty")

    mass <- .densityPointMass(data$x, data$y)
    peakPos <- .peakApexPosition(data$x, data$y)

    if (method == "hpd_span") {
        ranked <- order(data$y, decreasing = TRUE)
        cumulative <- cumsum(mass[ranked])
        firstEnough <- which(cumulative >= topP)[1]
        cutoffDensity <- data$y[ranked[firstEnough]]
        selected <- data$y >= cutoffDensity
        return(list(minPos = min(data$x[selected]),
                    maxPos = max(data$x[selected]),
                    peakPos = peakPos))
    }

    cs <- c(0, cumsum(mass))
    bestLeft <- 1L
    bestRight <- nrow(data)
    bestWidth <- Inf
    for (left in seq_len(nrow(data))) {
        target <- cs[left] + topP
        right <- which(cs[-1L] >= target)[1]
        if (is.na(right) || right < left) next
        width <- data$x[right] - data$x[left]
        if (width < bestWidth) {
            bestWidth <- width
            bestLeft <- left
            bestRight <- right
        }
    }
    list(minPos = data$x[bestLeft],
         maxPos = data$x[bestRight],
         peakPos = peakPos)
}


.derivedSeed <- function(baseSeed, seqname) {
    # Deterministic per chromosome and stable across serial/parallel execution.
    chars <- utf8ToInt(as.character(seqname))
    offset <- if (length(chars)) sum(chars * seq_along(chars)) else 0
    as.integer((as.double(baseSeed) + offset) %% (.Machine$integer.max - 1L) + 1L)
}


.peakRefinementChr <- function(inputList, mmapprData) {
    stopifnot("seqname" %in% names(inputList))
    seqname <- inputList$seqname

    loessObj <- mmapprData@snpDistance[[seqname]]$loess
    loessSpan <- loessObj$pars$span
    # loess stores predictors in a matrix-like `x` component on some R
    # versions. Normalize the single genomic predictor to a plain numeric vector so
    # downstream data.frame/subsetting behavior is unambiguous.
    rawData <- data.frame(pos = as.numeric(loessObj$x),
                          euclideanDistance = as.numeric(loessObj$y))
    rawData <- rawData[is.finite(rawData$pos) & is.finite(rawData$euclideanDistance), , drop = FALSE]
    if (nrow(rawData) < 10L) stop("Too few finite markers to refine peak on ", seqname)

    # Peak refinement used to depend on the ambient RNG state, so identical
    # inputs could yield slightly different intervals. Use a deterministic seed
    # derived from the user-recorded base seed and chromosome name.
    seed <- .derivedSeed(randomSeed(mmapprData@param), seqname)
    oldSeedExists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (oldSeedExists) oldSeed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
        if (oldSeedExists) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)

    nIter <- as.integer(peakResampleIterations(mmapprData@param))
    maxValues <- replicate(nIter,
                           .getSubsampleLoessMax(rawData = rawData,
                                                loessSpan = loessSpan))
    success <- is.finite(maxValues)
    successRate <- mean(success)
    maxValues <- maxValues[success]
    if (length(maxValues) < max(5L, ceiling(0.5 * nIter)))
        stop("Too many failed peak-resampling LOESS fits on ", seqname,
             " (successful: ", length(maxValues), "/", nIter, ")")

    # A degenerate distribution can occur when every resample chooses exactly the
    # same marker. In that case density() cannot estimate bandwidth sensibly.
    if (length(unique(maxValues)) == 1L) {
        densityData <- list(x = maxValues[[1]], y = 1)
        peak <- list(minPos = maxValues[[1]], maxPos = maxValues[[1]],
                     peakPos = maxValues[[1]])
        densityFunction <- function(x) as.numeric(x == maxValues[[1]])
    } else {
        densityData <- stats::density(maxValues)
        densityDf <- data.frame(x = densityData$x, y = densityData$y)

        # KDE support can extend beyond the chromosome. Clip it before using
        # the density to define a physical genomic interval.
        fa_si <- .faSeqinfo(mmapprData@param@refGenome)
        chrLen <- GenomeInfoDb::seqlengths(fa_si)[seqname]
        if (length(chrLen) == 1L && is.finite(chrLen)) {
            densityDf <- densityDf[densityDf$x >= 1 & densityDf$x <= chrLen, , drop = FALSE]
        }
        if (nrow(densityDf) == 0L) stop("No valid KDE support remained on ", seqname)

        peak <- .getPeakFromTopP(densityDf, peakIntervalWidth(mmapprData@param),
                                 method = peakIntervalMethod(mmapprData@param))
        densityData$x <- densityDf$x
        densityData$y <- densityDf$y
        densityFunction <- stats::approxfun(x = densityDf$x, y = densityDf$y,
                                            yleft = 0, yright = 0, rule = 1)
    }

    # Keep integer, valid coordinates for downstream GRanges queries.
    fa_si <- .faSeqinfo(mmapprData@param@refGenome)
    chrLen <- GenomeInfoDb::seqlengths(fa_si)[seqname]
    lower <- max(1L, floor(peak$minPos))
    upper <- ceiling(peak$maxPos)
    if (length(chrLen) == 1L && is.finite(chrLen)) upper <- min(upper, as.integer(chrLen))

    loessX <- as.numeric(loessObj$x)
    loessPeak <- .peakApexPosition(loessX, loessObj$fitted)

    # This data model still represents one refined locus per chromosome.
    # Proper polygenic / multiple-same-chromosome support requires changing the
    # peaks structure and downstream output, not simply choosing a second maximum.
    list(seqname = seqname,
         start = lower,
         end = upper,
         densityFunction = densityFunction,
         peakPosition = peak$peakPos,       # resampling-density apex; preserves the old field meaning
         densityPeakPosition = peak$peakPos,
         loessPeakPosition = loessPeak,
         densityData = densityData,
         resampleSuccessRate = successRate,
         resampleSeed = seed,
         cutoff = inputList$cutoff,
         cutoffCenter = inputList$cutoffCenter,
         cutoffSpread = inputList$cutoffSpread,
         cutoffMethod = inputList$cutoffMethod,
         intervalMethod = peakIntervalMethod(mmapprData@param))
}


#' @title Identify chromosomes containing peaks
#'
#' @name prePeak
#' @description
#' Follows the \code{\link{loessFit}} step and precedes
#' \code{\link{peakRefinement}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \linkS4class{MmapprData} object with the \code{peaks}
#'   slot initialized.
#' @export
NULL


.calculatePeakCutoff <- function(snpDistance, method = c("legacy_current", "global_sd"), k = 3) {
    method <- match.arg(method)
    valid <- lapply(snpDistance, function(chr) {
        if (!is.list(chr) || !inherits(chr$loess, "loess")) return(numeric())
        chr$loess$fitted[is.finite(chr$loess$fitted)]
    })
    valid <- valid[vapply(valid, length, integer(1)) > 0L]
    if (length(valid) == 0L) stop("No finite LOESS values are available for peak thresholding")

    if (method == "legacy_current") {
        centers <- vapply(valid, mean, numeric(1))
        varianceTerms <- vapply(valid, function(x) if (length(x) > 1L) stats::var(x) / length(x) else 0, numeric(1))
        center <- mean(centers)
        spread <- sqrt(sum(varianceTerms))
    } else {
        # Optional statistically clearer alternative: a robust global
        # center with the ordinary SD of all fitted marker values. This is still
        # a heuristic threshold, not a calibrated false-positive probability.
        allFitted <- unlist(valid, use.names = FALSE)
        if (length(allFitted) < 2L) stop("Insufficient finite LOESS values for global-SD thresholding")
        center <- stats::median(allFitted)
        spread <- stats::sd(allFitted)
    }
    list(cutoff = center + k * spread, center = center, spread = spread, method = method)
}


prePeak <- function(mmapprData) {
    mmapprData@peaks <- list()
    cutoffInfo <- .calculatePeakCutoff(
        mmapprData@snpDistance,
        method = peakCutoffMethod(mmapprData@param),
        k = peakCutoffSd(mmapprData@param)
    )

    for (i in seq_along(mmapprData@snpDistance)) {
        chr <- mmapprData@snpDistance[[i]]
        if (!is.list(chr) || !inherits(chr$loess, "loess")) next
        loessForChr <- chr$loess
        if (length(loessForChr$x) < 50L) next
        # NA fitted values no longer make any() return NA and break if().
        containsPeak <- any(loessForChr$fitted > cutoffInfo$cutoff, na.rm = TRUE)
        chrName <- names(mmapprData@snpDistance)[[i]]
        if (containsPeak) {
            mmapprData@peaks[[chrName]] <- list(
                seqname = chrName,
                cutoff = cutoffInfo$cutoff,
                cutoffCenter = cutoffInfo$center,
                cutoffSpread = cutoffInfo$spread,
                cutoffMethod = cutoffInfo$method
            )
        }
    }

    .messageAndLog(sprintf("Peak cutoff (%s): center %.6g + %.3g*spread %.6g = %.6g",
                           cutoffInfo$method, cutoffInfo$center,
                           peakCutoffSd(mmapprData@param), cutoffInfo$spread,
                           cutoffInfo$cutoff),
                   outputFolder(mmapprData@param))
    mmapprData
}
