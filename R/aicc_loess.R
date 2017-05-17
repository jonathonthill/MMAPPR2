LoessFit <- function(mmapprData, silent = FALSE) {
  loessOptResolution <- mmapprData@param@loessOptResolution
  loessOptCutFactor <- mmapprData@param@loessOptCutFactor
  
  #each item (chr) of distance list has mutCounts, wtCounts, distanceDf going in
  mmapprData@distance <- RunFunctionInParallel(mmapprData@distance, functionToRun = .loessFitForChr, silent = silent,
                                               secondInput = loessOptResolution, thirdInput = loessOptCutFactor,
                                               numCores = mmapprData@param@numCores, packages = c("MMAPPR2"))
  #LoessFitForChr returns list with mutCounts, wtCounts, loess, aicc
  
  return(mmapprData)
}


.GetLoess <- function(s, pos, euc_dist, ...){
  x <- try(loess(euc_dist ~ pos, span=s, degree=1, 
                 family="symmetric", ...), silent=T)
  return(x)
}

#gets greater of the two differences to either side of a span in a vector of spans
.localResolution <- function(spans, span) {
  stopifnot(is.numeric(spans))
  stopifnot(all(spans >= 0 & spans <= 1))
  
  spans <- unique(sort(spans))
  index <- which(spans %in% span)
  stopifnot(length(index) == 1)
  
  result <- max(c(diff(spans)[index-1], diff(spans)[index]), na.rm = T)
  stopifnot(is.numeric(result))
  return(result)
}

#returns a list of aicc and span values as well as the time it took
#define needed functions
.AiccOpt <- function(distanceDf, spans, resolution, cutFactor, showDebug = FALSE){
  aiccValues <- sapply(spans, .Aicc, euc_dist = distanceDf$distance, pos = distanceDf$pos)
  if (length(spans) != length(aiccValues)) stop("AICc values and spans don't match")
  aiccDf <- data.frame(spans, aiccValues)
  
  #finds two lowest local minima for finer attempt
  minVals <- .MinTwo(.LocalMin(aiccDf$aiccValues))
  #gets first column of dataframe after selecting for rows that match the two lowest localMins
  minSpans <- aiccDf[aiccDf$aiccValues %in% minVals, 1]
  
  #goes through each minimum and goes deeper if resolution isn't fine enough
  for (minSpan in minSpans) {
    #get resolution at that point, which is greater of differences to either side
    localResolution <- round(.localResolution(spans, minSpan), digits = .NumDecimals(resolution))
    stopifnot(is.numeric(localResolution))
    if (abs(localResolution) > resolution){
      addVector <- ((localResolution * cutFactor) * 1:9)
      newSpans <- c(minSpan - addVector, minSpan + addVector)
      newSpans <- newSpans[newSpans > 0 & newSpans <= 1]
      newSpans <- round(newSpans, digits = .NumDecimals(resolution))
      newSpans <- unique(newSpans)
      if (showDebug) message("# new spans = ", length(newSpans))
      #recursive call
      aiccDf <- rbind(aiccDf, .AiccOpt(distanceDf, spans = newSpans, 
                                       resolution = resolution, cutFactor = cutFactor))
    }
  }
  return(aiccDf)
}

#returns minimum two elements (useful for aicc_opt)
.MinTwo <- function(x){
  len <- length(x)
  if(len<2){
    warning('len < 2; returning x')
    return(x)
  }
  sort(x,partial=c(1, 2))[c(1, 2)]
}

#returns all local minima (problem if repeated local maxima on end)
.LocalMin <- function(x){
  indices <- which(diff(c(FALSE,diff(x)>0,TRUE))>0)
  return(x[indices])
}

.NumDecimals <- function(x) {
  stopifnot(class(x)=="numeric")
  if (!grepl("[.]", x)) return(0)
  x <- sub("0+$","",x)
  x <- sub("^.+[.]","",x)
  nchar(x)
}

.Aicc <- function (s, euc_dist, pos) {
  # extract values from loess object
  x <- .GetLoess(s, pos, euc_dist)
  if(class(x)=="try-error") return(NA)
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum( x$residuals^2 ) / (n-1)
  delta1 <- x$one.delta
  delta2 <- x$two.delta
  enp <- x$enp
  #return aicc value
  return(log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2))
  #what I understand:
  #return(n * log(sigma2) + 2*traceL + 2*traceL*(traceL + 1) / (n - traceL - 1))
}


#the function that gets run for each chromosome
#takes element of mmapprData@distance (with mutCounts, wtCounts, distanceDf)
#outputs complete element of mmapprData@distance (mutcounts, wtCounts, loess, aicc)
.loessFitForChr <- function(resultList, loessOptResolution, loessOptCutFactor,
                           showDebug = FALSE){
  
  startTime <- proc.time()
  tryCatch({
    if(class(resultList) == 'character') stop('--Loess fit failed')
    
    #returns dataframe with spans and aicc values for each loess
    startSpans <- .01 * c(1:16, 1 + 10*2:9)
    resultList$aicc <- MMAPPR2:::.AiccOpt(distanceDf = resultList$distanceDf, 
                                          spans = startSpans,
                                          resolution = loessOptResolution,
                                          cutFactor = loessOptCutFactor,
                                          showDebug = showDebug)
    
    #now get loess for best aicc
    bestSpan <- round(resultList$aicc[resultList$aicc$aiccValues == min(resultList$aicc$aiccValues, na.rm = T), 'spans'],
                      digits = .NumDecimals(loessOptResolution))
    bestSpan <- mean(bestSpan, na.rm = TRUE)
    resultList$loess <- MMAPPR2:::.GetLoess(bestSpan, resultList$distanceDf$pos, 
                                 resultList$distanceDf$distance, surface = "direct")
    
    precision <- .NumDecimals(loessOptResolution)
    message(sprintf(paste0("%s: LoessFit complete with optimized span of %.", precision, "f"), resultList$seqname, bestSpan))
    
    #no longer needed
    resultList$distanceDf <- NULL
    resultList$seqname <- NULL
    
    resultList$loessTime <- proc.time() - startTime
    return(resultList)
  },
  error = function(e) {
    if (class(resultList) == "character")
      return(paste0(resultList, e$message))
    else return(e$message)
  }
  )
}

