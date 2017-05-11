library(dplyr)

refinement_function <- function(mmapprData){
  options(warn=-1)
  
  for(chr in names(mmapprData@peaks)){
    loess_span_value <- mmapprData@distance[[chr]]$loess$pars$span
    pos <- mmapprData@distance[[chr]]$loess$x
    euc <- mmapprData@distance[[chr]]$loess$fitted
    raw_data <- data.frame(pos, euc)
  
  max_values <- rep(NA, 1000)
  for (i in 1:1000){
    temp_data <- raw_data[sample(1:nrow(raw_data), size = nrow(raw_data)/2),]
    temp_data <- temp_data[order(temp_data$pos),]
    loess_dat <- loess(euc~pos,data=temp_data,span = loess_span_value, family = c("symmetric"))
    max_values[i] <- loess_dat$x[which.max(loess_dat$fitted)]
  }
  
  density_data <- density.default(max_values)

  plot(density_data)

  density_function <- approxfun(x=density_data$x,y=density_data$y)
  x_min <- min(density_data$x)
  x_max <- max(density_data$x)
  
  data <- data.frame(seq(x_min,x_max))
  names(data) <- "pos"
  data <- mutate(data, density_value = density_function(seq(x_min,x_max))) %>% arrange(desc(density_value))

  rolling_sum <- cumsum(data$density_value)
  cutoff_pos <- which(rolling_sum > 0.95)[1]
  cutoff_val <- data$density_value[cutoff_pos]
  
  data <- filter(data,density_value >= cutoff_val)
  
  min = min(data$pos)
  max = max(data$pos)
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
