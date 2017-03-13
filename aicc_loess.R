
#the function that gets run for each chromosome
#takes chromosome range (in class GRanges)
#returns list of 1: Loess fit object (with optimum span) and 2: dataframe of AICc calculations
LoessFit <- function(chrRange){
  lf.start <- proc.time()
  dist_df = GetDistanceDf(chrRange)
  if(class(dist_df) == 'try-error') stop('Error: empty dataframe')
  
  aicc_start <- proc.time()
  aicc_df <- aicc_opt(dist_df) #returns dataframe with spans and aicc values for each loess
  message("aicc for chr ", names(chrRange), " took: ", proc.time() - aicc_start)
  
  #now get loess for best aicc
  best_span <- aicc_df[aicc_df$aicc_v == min(aicc_df$aicc_v, na.rm = T),'spans']
  chr.loess.fit <- get.loess(best_span, dist_df$distance, dist_df$pos)

  
  #loess fit(df, aiccParam) #keep all and pick best of list? or takes too much memory
  #  recursive: (min, max, interval parameters)
  # each chromosome has list: $fit(loess object WITH ACGTcvg FOR MUTANT), $aicc(vector of aicc data)
  message("whole LoessFit for chr ", names(chrRange), " took: ", proc.time() - aicc_start)
  return(list(fit = chr.loess.fit, aicc = aicc_df)) #actually return list of fit and aicc
}

#returns a list of aicc and span values as well as the time it took
aicc_opt <- function(in_df, spans = (0.01 + 0.1*(0:9)), resolution = loess.opt.resolution, cut.factor = loess.opt.cut.factor){
  aicc_v <- sapply(spans, aicc, euc_dist = in_df$distance, pos = in_df$pos)
  df <- data.frame(spans, aicc_v)
  current_resolution <- round(min(abs(diff(spans))), digits = num.decimals(resolution))
  if (abs(current_resolution) > resolution){
    #finds two lowest local minima for finer attempt
    min_vals <- minTwo(localMin(df$aicc_v))
    #gets first column of dataframe after selecting for rows that match the two lowest localMins
    min_spans <- df[df$aicc_v %in% min_vals, 1]
    add_vector <- ((current_resolution * cut.factor) * 1:9)
    new_spans <- c()
    for (x in min_spans){
      new_spans <- append(new_spans, x - add_vector)
      new_spans <- append(new_spans, x + add_vector)
    }
    new_spans <- new_spans[new_spans > 0]
    new_spans <- new_spans[new_spans <= 1]
    new_spans <- round(new_spans, digits = num.decimals(resolution))
    new_spans <- unique(new_spans)
    message("# new spans = ", length(new_spans))
    return(rbind(df,aicc_opt(in_df, spans = new_spans, resolution = resolution, cut.factor = cut.factor)))
  }
  else return(df)
}

#returns minimum two elements (useful for aicc_opt)
minTwo <- function(x){
  len <- length(x)
  if(len<2){
    warning('len < 2; returning x')
    return(x)
  }
  sort(x,partial=c(1, 2))[c(1, 2)]
}

#returns all local minima (problem if repeated local maxima on end)
localMin <- function(x){
  indices <- which(diff(c(FALSE,diff(x)>0,TRUE))>0)
  return(x[indices])
}

num.decimals <- function(x) {
  stopifnot(class(x)=="numeric")
  x <- sub("0+$","",x)
  x <- sub("^.+[.]","",x)
  nchar(x)
}

get.loess <- function(s, euc_dist, pos){
  x <- try(loess(euc_dist ~ pos, span=s, degree=1, family="symmetric", surface='direct'), silent=T)
  return(x)
}

aicc <- function (s, euc_dist, pos) {
  # extract values from loess object
  x <- get.loess(s, euc_dist, pos)
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