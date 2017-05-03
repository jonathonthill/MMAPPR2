library(dplyr)

#actual program will not need start, stop, max, or actual. 
#Raw_Data: Data frame with raw loess_data in it
#loess_span: optimized loess span value found previously

peak_function <- function(raw_data, start, stop, max, actual, loess_span){
    options(warn=-1)
  #simulation running loess 1000 times on all data keeping track of the peak
  max_values <- rep(NA, 1000) #container vector variable to hold all max values
  
  for (i in 1:1000){
    temp_data <- raw_data[sample(1:nrow(raw_data), size = nrow(raw_data)/2),]
    temp_data <- temp_data[order(temp_data$pos),]
    loess_dat <- loess(euc^4~pos,data=temp_data,span = loess_span, family = c("symmetric")) #span will change for each data set
    max_values[i] <- loess_dat$x[which.max(loess_dat$fitted)]
  }
  density_data <- density(max_values)
  max_of_density<-density_data$x[which.max(density_data$y)]
  
  #plot:REMOVE
  plot(density(max_values),main = "Density of Max Values Generated from Sampling Simulation on All Data", xlab = "Chromosome 14 Position")
  abline(v = actual, col = "red") #actual mutation position
  abline(v =max_of_density, col = "blue",lty=3 )#max of density plot
  abline(v=start, col = "blue" ) #STARTS
  abline(v=stop, col = "blue") #STOPS
  abline( v= max, col ="green") #MAX
  
  legend("topright", 
         c(paste("Actual mutation position(",actual,")",sep = ""),
           paste("New Max of Density Plot(",trunc(max_of_density),")",sep = ""),
           paste("Original Max(",max,")",sep = ""),
           paste("Original Start(",start,")",sep = ""),
           paste("Original stop(",stop,")",sep = "")),
         lty=c(1,3,1,1,1),
         col= c("red", "blue", "green","blue", "blue" ),cex = 1)
  
  options(warn = 0)
  return (density_data)
}

PrePeak <- function(mmapprData) {
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
    containsPeak <- any(loessForChr$fitted > cutoff)
    chrName <- names(mmapprData@distance)[[i]]
    if (containsPeak) {
      mmapprData@peaks[[chrName]] <- list(seqname = chrName)
      cat("Sequence", chrName, "contains peak.\n")
    }	
  }
  
  return(mmapprData)
}