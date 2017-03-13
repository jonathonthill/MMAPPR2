###contains peak caller function from original MMAPPR
###used as stand-in to develop variant caller

GetPeak <- function(loess.fit.list) {
  #need to calculate standard dev of all chromosomes for cutoff
  comb.sd <- sapply(loess.fit.list, FUN = function(chr){
    return(var(chr$fit$fitted)/length(chr$fit$fitted))})
  comb.sd <- sum(comb.sd)^(1/2)
  distance.median <- sapply(loess.fit.list, function(x) {return(median(x$fit$fitted))})
  distance.median <- median(distance.median) #median of list of chr medians
  cutoff <- 3*comb.sd + distance.median
  
  #get peak intervals
  cat("Using", round(cutoff, digits=4), "as cutoff.\n")
  starts <- NULL
  stops <- NULL
  peakchr <- NULL
  for(i in seq_along(loess.fit.list)){
    chr <- loess.fit.list[[i]]
    # chr.active <- plot.df[plot.df$chr==chr & !is.na(plot.df$fitted),]
    # for(j in 1:nrow(chr.active)){
    #   if (j==0){next}
    #   if (is.na(chr.active$fitted[j])){next}
    #   if (j==nrow(chr.active) & chr.active$fitted[j] > cutoff) {stops <- append(stops, chr.active$pos[j]); next}
    #   if (j==1 & chr.active$fitted[j] > cutoff) {peakchr <- append(peakchr, chr); starts <- append(starts, 0); next}
    #   if (j!=nrow(chr.active) & chr.active$fitted[j] < cutoff & chr.active$fitted[j+1] > cutoff) {peakchr <- append(peakchr, chr); starts <- append(starts, chr.active$pos[j+1]); next}
    #   if (j!=nrow(chr.active) & chr.active$fitted[j] > cutoff & chr.active$fitted[j+1] < cutoff) {stops <- append(stops, chr.active$pos[j]); next}		
    peak.domain <- chr$fit$x[chr$fit$fitted > cutoff]
    if (length(peak.domain) > 0) {
      peakchr <- append(peakchr, names(loess.fit.list)[i])
      starts <- append(starts, peak.domain[1])
      stops <- append(stops, peak.domain[length(peak.domain)])
    }	
  }
    str(peakchr)
  
  peak.df <- data.frame(cbind(peakchr,starts,stops), stringsAsFactors = F)
  colnames(peak.df) <- c('chr','starts','stops')
  peak.df <- transform(peak.df, starts = as.numeric(starts), stops = as.numeric(stops))
  peakmax <- NULL
  if (nrow(peak.df) < 1) {
    cat("No peaks found")
    dev.off()
    quit()
  }
  
#   #puts position and height of each peak into dataframe
#   for(i in 1:nrow(peak.df)){
#     peak.df$peakmax[i] <- comb$pos[comb$fitted==max(comb$fitted[comb$pos > peak.df$starts[i] & comb$pos < peak.df$stops[i]], na.rm=T) & comb$pos > peak.df$starts[i] & comb$pos < peak.df$stops[i]]
#     peak.df$height[i] <- max(comb$fitted[comb$pos > peak.df$starts[i] & comb$pos < peak.df$stops[i]], na.rm=T) - median(comb$fitted)
#   }
#   cat("Peaks:\n")
#   write.table(peak.df, quote=FALSE, row.names=FALSE, sep="\t")
#   for(chr in unique(peak.df$chr)){
#     plot(plot.df$pos[plot.df$chr==chr]/1000000, plot.df$fitted[plot.df$chr==chr], type='l', ylim=c(min(plot.df$fitted, na.rm=T),1.1*max(plot.df$fitted, na.rm=T)), ylab=substitute("ED"^p~ ~"(Loess fit)", list(p=power)), xlab=paste(chr,"Base Position (MB)"), xaxs='i' )
#     chr.a <- transform(peak.df[peak.df$chr==chr,], starts = as.integer(starts), stops = as.integer(stops))
#     for(i in 1:nrow(chr.a)){
#       shade.x <- c(chr.a$starts[i],plot.df$pos[plot.df$pos >= chr.a$starts[i] & plot.df$pos <= chr.a$stops[i] & plot.df$chr==chr & !is.na(plot.df$fitte)],chr.a$stops[i])
#       shade.y <- c(-5, plot.df$fitted[plot.df$pos >= chr.a$starts[i] & plot.df$pos <= chr.a$stops[i] & plot.df$chr==chr & !is.na(plot.df$fitted)], -5)
#       polygon(shade.x/1000000,shade.y, col=rgb(0,1,1,.75), border=NA)
#     }
#     plot(plot.df$pos[plot.df$chr==chr]/1000000, plot.df$unfitted[plot.df$chr==chr], pch=16, cex=.6, ylim=c(min(plot.df$unfitted, na.rm=T),1.1*max(plot.df$unfitted, na.rm=T)), ylab=substitute("ED"^p, list(p=power)), xlab=paste(chr,"Base Position (MB)"), xaxs='i' )
#   }
#   noreport <- dev.off()
# }
  
  cat("Making allele file.", "\n")
  
  genotype <- function(a, c, g, t){
    maxbase <- max(c( a, c, g, t))
    if(maxbase == a){return("A")}
    if(maxbase == c){return("C")}
    if(maxbase == g){return("G")}
    if(maxbase == t){return("T")}
  }
  
  #starts empty dataframe with right columns
  alleler <- data.frame(t(rep(NA, 5)))
  names(alleler) <- c("chr", "pos", "pos.1", "call", "notes")
  
  for (i in 1:nrow(peak.df)) {
    
    #we need replacement for comb2, chr unnecessary since region will be just one chromosome
    #we need pos, euc.distance, and ACGTcvg for mutant pool
    #rescan file
    current.peak.range <- GRanges(seqnames = peak.df[i, 'chr'], ranges = IRanges(peak.df$starts[i], peak.df$stops[i]))
    current.peak.snps <- GetDistanceDf(current.peak.range, peak.version = T)
    
    #what's going on here?
    #x is control, y is mutant
    #gets in peak range, tests for homozygosity in mutant pool, tests coverage > 1, tests distance > 0.5
    tempalleler <- filter(current.peak.snps, (A.mut > 0.75 | C.mut > 0.75 | G.mut > 0.75 | T.mut > 0.75), distance > 0.5)
    tempalleler <- select(tempalleler, -matches("distance"))
    #tempalleler <- comb2[comb2$chr==peak.df$chr[i] & comb2$pos > peak.df$starts[i] & comb2$pos < peak.df$stops[i] & (comb2$A.y>.75 | comb2$G.y>.75 | comb2$C.y>.75 | comb2$T.y>.75) & comb2$cov.y > 1 & comb2$euc > .5, c('chr','pos','pos', 'ref.y', 'cov.y', 'A.y', 'C.y', 'G.y', 'T.y')]
    tempalleler <- cbind(peak.df[i, 'chr'], tempalleler)
    names(tempalleler) <- c("chr", "pos", "cov", "A", "C", "G", "T")
    #why are there two pos columns?
    tempalleler$pos.1 <- tempalleler$pos
    tempalleler$pos <- tempalleler$pos-1
    str(tempalleler)
    
    tempalleler$call <- mapply(genotype, tempalleler$A, tempalleler$C, tempalleler$G, tempalleler$T)
    #shouldn't need because of variant effect predictor, unless it's slow
    #tempalleler <- tempalleler[tempalleler$call != tempalleler$ref, ]
    
    tempalleler$notes <- mapply(function(x,a,c,g,t) paste("Coverage:A:C:G:T = ", x, ":", a, ":", c, ":", g, ":", t, sep=""), tempalleler$cov, as.integer(tempalleler$A*tempalleler$cov), as.integer(tempalleler$C*tempalleler$cov), as.integer(tempalleler$G*tempalleler$cov), as.integer(tempalleler$T*tempalleler$cov))
    alleler <- rbind(alleler, tempalleler[c('chr', 'pos', 'pos.1', 'call','notes')])
  }

  return(alleler)
}
alleler <- alleler[!is.na(alleler$chr), ]

