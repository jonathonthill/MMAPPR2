LoessFit <- function(mmapprData) {
  loessOptResolution <- mmapprData@input$loessOptResolution
  loessOptCutFactor <- mmapprData@input$loessOptCutFactor
  
  #each item (chr) of distance list has mutCounts, wtCounts, distanceDf going in
  mmapprData@distance <- RunFunctionInParallel(mmapprData@distance, LoessFitForChr)
  #LoessFitForChr returns list with mutCounts, wtCounts, loess, aicc
  
  return(mmapprData)
}


#the function that gets run for each chromosome
#takes element of mmapprData@distance (with mutCounts, wtCounts, distanceDf)
#outputs complete element of mmapprData@distance (mutcounts, wtCounts, loess, aicc)
LoessFitForChr <- function(resultList){
  #returns a list of aicc and span values as well as the time it took
  #define needed functions
  AiccOpt <- function(distanceDf, spans = (0.01 + 0.1*(0:9)), 
                      resolution = loessOptResolution, cutFactor = loessOptCutFactor){
    aiccValues <- sapply(spans, Aicc, euc_dist = distanceDf$distance, pos = distanceDf$pos)
    aiccDf <- data.frame(spans, aiccValues)
    currentResolution <- round(min(abs(diff(spans))), digits = NumDecimals(resolution))
    if (abs(currentResolution) > resolution){
      #finds two lowest local minima for finer attempt
      minVals <- MinTwo(LocalMin(aiccDf$aiccValues))
      #gets first column of dataframe after selecting for rows that match the two lowest localMins
      minSpans <- aiccDf[aiccDf$aiccValues %in% minVals, 1]
      addVector <- ((currentResolution * cutFactor) * 1:9)
      newSpans <- c()
      for (x in minSpans){
        newSpans <- append(newSpans, x - addVector)
        newSpans <- append(newSpans, x + addVector)
      }
      newSpans <- newSpans[newSpans > 0 & newSpans <= 1]
      newSpans <- round(newSpans, digits = NumDecimals(resolution))
      newSpans <- unique(newSpans)
      message("# new spans = ", length(newSpans))
      return(rbind(aiccDf,AiccOpt(distanceDf, spans = newSpans, resolution = resolution, cutFactor = cutFactor)))
    }
    else return(aiccDf)
  }
  
  #returns minimum two elements (useful for aicc_opt)
  MinTwo <- function(x){
    len <- length(x)
    if(len<2){
      warning('len < 2; returning x')
      return(x)
    }
    sort(x,partial=c(1, 2))[c(1, 2)]
  }
  
  #returns all local minima (problem if repeated local maxima on end)
  LocalMin <- function(x){
    indices <- which(diff(c(FALSE,diff(x)>0,TRUE))>0)
    return(x[indices])
  }
  
  NumDecimals <- function(x) {
    stopifnot(class(x)=="numeric")
    x <- sub("0+$","",x)
    x <- sub("^.+[.]","",x)
    nchar(x)
  }
  
  GetLoess <- function(s, euc_dist, pos){
    x <- try(loess(euc_dist ~ pos, span=s, degree=1, family="symmetric", surface='direct'), silent=T)
    return(x)
  }
  
  Aicc <- function (s, euc_dist, pos) {
    # extract values from loess object
    x <- GetLoess(s, euc_dist, pos)
    if(class(x)=="try-error"){return(NA)}
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
  
  startTime <- proc.time()
  tryCatch({
    if(class(resultList) != 'list') stop('Unable to perform loess fit for sequence--distance data missing.')
    
    resultList$aicc <- AiccOpt(resultList$distanceDf) #returns dataframe with spans and aicc values for each loess
    
    #now get loess for best aicc
    bestSpan <- resultList$aicc[resultList$aicc$aiccValues == min(resultList$aicc$aiccValues, na.rm = T), 'spans']
    resultList$loess <- GetLoess(bestSpan, resultList$distanceDf$distance, 
                                 resultList$distanceDf$pos)
    
    resultList$distanceDf <- NULL #no longer needed
    
    message("LoessFit for chr ", names(chrRange), " took: ", proc.time() - startTime)
    return(resultList)
  },
  error = function(e) {
    return(e$message)
  }
  )
}

