#' 
#' @title Characterize Euclidean distance peaks using resampling simulation
#'
#' @name peakRefinement
#'
#' @usage Follows the \code{\link{prePeak}} step and precedes
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
    mmapprData@peaks <-
        bplapply(mmapprData@peaks,
                 .peakRefinementChr,
                 mmapprData = mmapprData)
    return(mmapprData)
}

.getSubsampleLoessMax <- function(rawData, loessSpan) {
    tempData <- rawData[sample(seq_len(nrow(rawData)), 
                               size = nrow(rawData)*.5),]
    tempData <- tempData[order(tempData$pos),]
    loessData <- suppressWarnings(
        loess(euclideanDistance~pos, data = tempData,
              span = loessSpan, family = c("symmetric")))
    return(loessData$x[which.max(loessData$fitted)])
}

.getPeakFromTopP <- function(data, topP) {
    # Takes a dataframe of x and y and identfies the top
    # values totaling `topP` or greater. It returns
    # a peak defined by this region
    stopifnot(ncol(data) == 2)
    names(data) <- c('x', 'y')
    # this arrange part is slow.
    data <- dplyr::arrange(data, dplyr::desc(data$y))

    rollingSum <- cumsum(data$y)
    cutoffValue <- data$y[rollingSum >= topP][1]

    data <- dplyr::filter(data, data$y >= cutoffValue)

    minPos = min(data$x)
    maxPos = max(data$x)
    peakPos <- data$x[which.max(data$y)]

    return(list(minPos = minPos, maxPos = maxPos, peakPos = peakPos))
}

.peakRefinementChr <- function(inputList, mmapprData) {
    stopifnot('seqname' %in% names(inputList))
    seqname <- inputList$seqname


    loessSpan <- mmapprData@snpDistance[[seqname]]$loess$pars$span
    pos <- mmapprData@snpDistance[[seqname]]$loess$x
    euclideanDistance <- mmapprData@snpDistance[[seqname]]$loess$y
    rawData <- data.frame(pos, euclideanDistance)

    # get peak values of loess fits of 1000 subsamples
    maxValues <- replicate(1000, .getSubsampleLoessMax(rawData = rawData,
                                                       loessSpan = loessSpan))

    densityData <- density.default(maxValues)
    densityFunction <- approxfun(x = densityData$x, y = densityData$y)

    xMin <- min(densityData$x)
    xMax <- max(densityData$x)

    densityRank <- data.frame(seq(xMin, xMax))
    names(densityRank) <- "pos"
    densityRank <- dplyr::mutate(densityRank,
                      'densityValue' = densityFunction(seq(xMin,xMax)))

    peak <- .getPeakFromTopP(densityRank, mmapprData@param@peakIntervalWidth)

    outputList <- list()
    outputList$seqname <- seqname
    outputList$start <- peak$minPos
    outputList$end <- peak$maxPos
    outputList$densityFunction <- densityFunction
    outputList$peakPosition <- peak$peakPos
    outputList$densityData <- densityData

    return(outputList)
}


.stDevForChr <- function(chr) {
    tryCatch({
        return(var(chr$loess$fitted)/length(chr$loess$fitted))
    }, error = function(e) {
        return(0)
    })
}


.meanForChr <- function(chr) {
    tryCatch({
        return(mean(chr$loess$fitted))
    }, error = function(e) {
        return(NA)
    })
}


#' @title Identify chromosomes containing peaks
#'
#' @name prePeak
#'
#' @usage Follows the \code{\link{loessFit}} step and precedes
#' \code{\link{peakRefinement}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \linkS4class{MmapprData} object with the \code{peaks}
#'   slot initalized.
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
#'
#' postPrePeakMD <- prePeak(postLoessMD)
#' }
#' 
NULL


prePeak <- function(mmapprData) {
    mmapprData@peaks <- list()

    #need to calculate standard dev of all chromosomes for cutoff
    combinedStDev <- vapply(mmapprData@snpDistance, .stDevForChr, numeric(1))
    combinedStDev <- sum(combinedStDev)^(1/2)

    distancemean <- vapply(mmapprData@snpDistance, .meanForChr, numeric(1))
    # mean of list of chr means
    distancemean <- mean(distancemean, na.rm = TRUE)
    cutoff <- 3*combinedStDev + distancemean

    #get which peaks have values above cutoff, initialize them in mmapprData
    for(i in seq_along(mmapprData@snpDistance)){
        if(!is(mmapprData@snpDistance[[i]], 'list')) next
        if(!is(mmapprData@snpDistance[[i]]$loess, 'loess')) next

        loessForChr <- mmapprData@snpDistance[[i]]$loess
        if (length(loessForChr$x) < 50) next
        containsPeak <- any(loessForChr$fitted > cutoff)
        chrName <- names(mmapprData@snpDistance)[[i]]
        if (containsPeak) {
            mmapprData@peaks[[chrName]] <- list(seqname = chrName)
        }
    }
    return(mmapprData)
}
