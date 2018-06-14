peakRefinement <- function(mmapprData){
    mmapprData@peaks <- 
        .runFunctionInParallel(mmapprData@peaks,
                               .peakRefinementChr,
                               mmapprData=mmapprData)
    return(mmapprData)
}

.peakRefinementChr <- function(inputList, mmapprData) {
    stopifnot('seqname' %in% names(inputList))
    seqname <- inputList$seqname
    
    options(warn=-1)
    
    loessSpan <- mmapprData@distance[[seqname]]$loess$pars$span
    pos <- mmapprData@distance[[seqname]]$loess$x
    euclideanDistance <- mmapprData@distance[[seqname]]$loess$y
    rawData <- data.frame(pos, euclideanDistance)
    
    maxValues <- rep(NA, 1000)
    for (i in 1:1000){
        tempData <- rawData[sample(1:nrow(rawData), size=nrow(rawData)/2),]
        tempData <- tempData[order(tempData$pos),]
        loessData <- suppressWarnings(
            loess(euclideanDistance~pos, data=tempData,
                  span=loessSpan, family=c("symmetric")))
        maxValues[i] <- loessData$x[which.max(loessData$fitted)]
    }
    
    densityData <- density.default(maxValues)
    
    densityFunction <- approxfun(x=densityData$x, y=densityData$y)
    xMin <- min(densityData$x)
    xMax <- max(densityData$x)
    
    densityRank <- data.frame(seq(xMin, xMax))
    names(densityRank) <- "pos"
    densityRank <- dplyr::mutate(densityRank, 
                      'densityValue'=densityFunction(seq(xMin,xMax)))
    densityRank <- dplyr::arrange(densityRank, dplyr::desc('densityValue'))
    
    rollingSum <- cumsum(densityRank$densityValue)
    cutoffPosition <- which(rollingSum > mmapprData@param@peakIntervalWidth)[1]
    cutoffValue <- densityRank$densityValue[cutoffPosition]
    
    densityRank <- dplyr::filter(densityRank, 'densityValue' >= cutoffValue)
    
    minPos = min(densityRank$pos)
    maxPos = max(densityRank$pos)
    
    densityPlot <- plot(densityData)
    abline(v = c(minPos, maxPos))
    
    peakPosition <- densityRank$pos[which.max(densityRank$densityValue)]
    
    outputList <- list()
    outputList$seqname <- seqname
    outputList$start <- minPos
    outputList$end <- maxPos
    outputList$densityFunction <- densityFunction
    outputList$peakPosition <- peakPosition
    
    
    options(warn=0)
    return(outputList)
}


prePeak <- function(mmapprData) {
    mmapprData@peaks <- list()
    
    #need to calculate standard dev of all chromosomes for cutoff
    combinedStDev <- sapply(mmapprData@distance, FUN = function(chr){
        if(class(chr) == "list"){
            var(chr$loess$fitted)/length(chr$loess$fitted)
        } else {
            return(0)
        }
    })
    combinedStDev <- sum(combinedStDev)^(1/2)
    distanceMedian <- sapply(mmapprData@distance, function(chr) {
        if(class(chr)=="list"){
            median(chr$loess$fitted)}else{
                return(NA)
            }})
    distanceMedian <- median(distanceMedian,na.rm=TRUE) #median of list of chr medians
    cutoff <- 3*combinedStDev + distanceMedian
    
    cat("Using", round(cutoff, digits=4), "as cutoff.\n")
    
    #get which peaks have values above cutoff, initialize them in mmapprData
    for(i in seq_along(mmapprData@distance)){
        if(class(mmapprData@distance[[i]]) != "list"){
            next
        }
        loessForChr <- mmapprData@distance[[i]]$loess
        if (length(loessForChr$x) < 50) next
        containsPeak <- any(loessForChr$fitted > cutoff)
        chrName <- names(mmapprData@distance)[[i]]
        if (containsPeak) {
            mmapprData@peaks[[chrName]] <- list(seqname = chrName)
            cat("Sequence", chrName, "contains peak.\n")
        }
    }
    
    return(mmapprData)
}
