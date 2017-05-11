library(dplyr)

PeakRefinement <- function(mmapprData){
  options(warn=-1)
  
  for(chr in names(mmapprData@peaks)){
    loessSpan <- mmapprData@distance[[chr]]$loess$pars$span
    pos <- mmapprData@distance[[chr]]$loess$x
    euclideanDistance <- mmapprData@distance[[chr]]$loess$fitted
    rawData <- data.frame(pos, euc)
  
  maxValues <- rep(NA, 1000)
  for (i in 1:1000){
    tempData <- rawData[sample(1:nrow(rawData), size = nrow(rawData)/2),]
    tempData <- tempData[order(tempData$pos),]
    loessData <- loess(euclideanDistance~pos,data=tempData,span = loessSpan, family = c("symmetric"))
    maxValues[i] <- loessData$x[which.max(loessData$fitted)]
  }
  
  densityData <- density.default(maxValues)

  plot(densityData)

  densityFunction <- approxfun(x=density_data$x,y=density_data$y)
  x_min <- min(densityData$x)
  x_max <- max(densityData$x)
  
  densityRank <- data.frame(seq(x_min,x_max))
  names(densityRank) <- "pos"
  densityRank <- mutate(densityRank, densityValue = density_function(seq(x_min,x_max))) %>% arrange(desc(densityValue))

  rollingSum <- cumsum(data$densityValue)
  cutoffPosition <- which(rollingSum > 0.95)[1]
  cutoffValue <- data$densityValue[cutoffPosition]
  
  densityRank <- filter(densityRank,densityValue >= cutoffValue)
  
  min = min(densityRank$pos)
  max = max(densityRank$pos)
  abline(v = c(min,max))
  #set these equal to mapprrdata
  
  mmapprData@peaks[[chr]]$start <- min
  mmapprData@peaks[[chr]]$end <- max
  mmapprData@peaks[[chr]]$densityFunction <- density_function
  
  options(warn = 0)
  } #end of for loop
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
