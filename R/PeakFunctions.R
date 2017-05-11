library(dplyr)

PeakRefinement <- function(mmapprData){
  options(warn=-1)
  
  for(chr in names(mmapprData@peaks)){
    loessSpan <- mmapprData@distance[[chr]]$loess$pars$span
    pos <- mmapprData@distance[[chr]]$loess$x
    euclideanDistance <- mmapprData@distance[[chr]]$loess$fitted
    rawData <- data.frame(pos, euclideanDistance)
    
    maxValues <- rep(NA, 1000)
    for (i in 1:1000){
      tempData <- rawData[sample(1:nrow(rawData), size = nrow(rawData)/2),]
      tempData <- tempData[order(tempData$pos),]
      loessData <- loess(euclideanDistance~pos,data=tempData,span = loessSpan, family = c("symmetric"))
      maxValues[i] <- loessData$x[which.max(loessData$fitted)]
    }
    
    densityData <- density.default(maxValues)
    
    densityFunction <- approxfun(x=densityData$x,y=densityData$y)
    xMin <- min(densityData$x)
    xMax <- max(densityData$x)
    
    densityRank <- data.frame(seq(xMin,xMax))
    names(densityRank) <- "pos"
    densityRank <- mutate(densityRank, densityValue = densityFunction(seq(xMin,xMax))) %>% arrange(desc(densityValue))
    
    rollingSum <- cumsum(densityRank$densityValue)
    cutoffPosition <- which(rollingSum > 0.95)[1]
    cutoffValue <- densityRank$densityValue[cutoffPosition]
    
    densityRank <- filter(densityRank,densityValue >= cutoffValue)
    
    min = min(densityRank$pos)
    max = max(densityRank$pos)
    
    plot(densityData)
    abline(v = c(min,max))
    
    peakPosition <- densityRank$pos[which.max(densityRank$densityValue)]

    mmapprData@peaks[[chr]]$start <- min
    mmapprData@peaks[[chr]]$end <- max
    mmapprData@peaks[[chr]]$densityFunction <- densityFunction
    mmapprData@peaks[[chr]]$peakPosition <- peakPosition
    
    options(warn = 0)
  } 
  return (mmapprData)
}


PrePeak <- function(mmapprData) {
  mmapprData@peaks <- list()
  
  #need to calculate standard dev of all chromosomes for cutoff
  combinedStDev <- sapply(mmapprData@distance, FUN = function(chr){
    var(chr$loess$fitted)/length(chr$loess$fitted)})
  combinedStDev <- sum(combinedStDev)^(1/2)
  distanceMedian <- sapply(mmapprData@distance, function(chr) {
    median(chr$loess$fitted)})
  distanceMedian <- median(distanceMedian) #median of list of chr medians
  cutoff <- 3*combinedStDev + distanceMedian
  
  cat("Using", round(cutoff, digits=4), "as cutoff.\n")
  
  #get which peaks have values above cutoff, initialize them in mmapprData
  for(i in seq_along(mmapprData@distance)){
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
