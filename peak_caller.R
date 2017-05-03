###contains peak caller function from original MMAPPR
###used as stand-in to develop variant caller

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
    chrName <- names(mmapprData@distance)
    if (containsPeak) {
      mmapprData@peaks[[chrName]] <- list(seqname = chrName)
    }	
  }

  return(mmapprData)
}

#   
# #   #puts position and height of each peak into dataframe
# #   for(i in 1:nrow(peak.df)){
# #     peak.df$peakmax[i] <- comb$pos[comb$fitted==max(comb$fitted[comb$pos > peak.df$starts[i] & comb$pos < peak.df$stops[i]], na.rm=T) & comb$pos > peak.df$starts[i] & comb$pos < peak.df$stops[i]]
# #     peak.df$height[i] <- max(comb$fitted[comb$pos > peak.df$starts[i] & comb$pos < peak.df$stops[i]], na.rm=T) - median(comb$fitted)
# #   }
# #   cat("Peaks:\n")
# #   write.table(peak.df, quote=FALSE, row.names=FALSE, sep="\t")
# #   for(chr in unique(peak.df$chr)){
# #     plot(plot.df$pos[plot.df$chr==chr]/1000000, plot.df$fitted[plot.df$chr==chr], type='l', ylim=c(min(plot.df$fitted, na.rm=T),1.1*max(plot.df$fitted, na.rm=T)), ylab=substitute("ED"^p~ ~"(Loess fit)", list(p=power)), xlab=paste(chr,"Base Position (MB)"), xaxs='i' )
# #     chr.a <- transform(peak.df[peak.df$chr==chr,], starts = as.integer(starts), stops = as.integer(stops))
# #     for(i in 1:nrow(chr.a)){
# #       shade.x <- c(chr.a$starts[i],plot.df$pos[plot.df$pos >= chr.a$starts[i] & plot.df$pos <= chr.a$stops[i] & plot.df$chr==chr & !is.na(plot.df$fitte)],chr.a$stops[i])
# #       shade.y <- c(-5, plot.df$fitted[plot.df$pos >= chr.a$starts[i] & plot.df$pos <= chr.a$stops[i] & plot.df$chr==chr & !is.na(plot.df$fitted)], -5)
# #       polygon(shade.x/1000000,shade.y, col=rgb(0,1,1,.75), border=NA)
# #     }
# #     plot(plot.df$pos[plot.df$chr==chr]/1000000, plot.df$unfitted[plot.df$chr==chr], pch=16, cex=.6, ylim=c(min(plot.df$unfitted, na.rm=T),1.1*max(plot.df$unfitted, na.rm=T)), ylab=substitute("ED"^p, list(p=power)), xlab=paste(chr,"Base Position (MB)"), xaxs='i' )
# #   }
# #   noreport <- dev.off()
# # }
#   
#   cat("Making allele file.", "\n")
#   
#   genotype <- function(a, c, g, t){
#     maxbase <- max(c( a, c, g, t))
#     if(maxbase == a){return("A")}
#     if(maxbase == c){return("C")}
#     if(maxbase == g){return("G")}
#     if(maxbase == t){return("T")}
#   }
#   
#   #starts empty dataframe with right columns
#   alleler <- data.frame(t(rep(NA, 5)))
#   names(alleler) <- c("chr", "pos", "pos.1", "call", "notes")
#   
#   for (i in 1:nrow(peak.df)) {
#     
#     #we need replacement for comb2, chr unnecessary since region will be just one chromosome
#     #we need pos, euc.distance, and ACGTcvg for mutant pool
#     #rescan file
#     current.peak.range <- GRanges(seqnames = peak.df[i, 'chr'], ranges = IRanges(peak.df$starts[i], peak.df$stops[i]))
#     current.peak.snps <- GetDistanceDf(current.peak.range, peak.version = T)
#     
#     #what's going on here?
#     #x is control, y is mutant
#     #gets in peak range, tests for homozygosity in mutant pool, tests coverage > 1, tests distance > 0.5
#     tempalleler <- filter(current.peak.snps, (A.mut > 0.75 | C.mut > 0.75 | G.mut > 0.75 | T.mut > 0.75), distance > 0.5)
#     tempalleler <- select(tempalleler, -matches("distance"))
#     #tempalleler <- comb2[comb2$chr==peak.df$chr[i] & comb2$pos > peak.df$starts[i] & comb2$pos < peak.df$stops[i] & (comb2$A.y>.75 | comb2$G.y>.75 | comb2$C.y>.75 | comb2$T.y>.75) & comb2$cov.y > 1 & comb2$euc > .5, c('chr','pos','pos', 'ref.y', 'cov.y', 'A.y', 'C.y', 'G.y', 'T.y')]
#     tempalleler <- cbind(peak.df[i, 'chr'], tempalleler)
#     names(tempalleler) <- c("chr", "pos", "cov", "A", "C", "G", "T")
#     #why are there two pos columns?
#     tempalleler$pos.1 <- tempalleler$pos
#     tempalleler$pos <- tempalleler$pos-1
#     str(tempalleler)
#     
#     tempalleler$call <- mapply(genotype, tempalleler$A, tempalleler$C, tempalleler$G, tempalleler$T)
#     #shouldn't need because of variant effect predictor, unless it's slow
#     #tempalleler <- tempalleler[tempalleler$call != tempalleler$ref, ]
#     
#     tempalleler$notes <- mapply(function(x,a,c,g,t) paste("Coverage:A:C:G:T = ", x, ":", a, ":", c, ":", g, ":", t, sep=""), tempalleler$cov, as.integer(tempalleler$A*tempalleler$cov), as.integer(tempalleler$C*tempalleler$cov), as.integer(tempalleler$G*tempalleler$cov), as.integer(tempalleler$T*tempalleler$cov))
#     alleler <- rbind(alleler, tempalleler[c('chr', 'pos', 'pos.1', 'call','notes')])
#   }
# 
#   # should get rid of 1 NA row created because of rbind with empty dataframe
#   alleler <- alleler[!is.na(alleler$chr), ]
#   return(alleler)
# }

