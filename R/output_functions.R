OutputMmapprData <- function(mmapprData, plotAicc = FALSE) {
  if (class(mmapprData) != "MmapprData"){
    stop("Input object not of 'MmapprData' type")
  }
  
  if (!dir.exists('mmappr_results')) dir.create("mmappr_results")
  
  PlotGenomeDistance(mmapprData@distance)
  PlotPeaks(mmapprData)
  if (plotAicc) PlotAicc(mmapprData@distance)
  
  candidateVCF <- asVCF(mmapprData)
}

PlotGenomeDistance <- function(mmapprData) {
  #generate one big dataframe for plots, along with break and label points
  tailPos <- 0
  plotDf <- NULL
  breaks <- tailPos
  labelpos <- NULL
  for (i in orderSeqlevels(names(mmapprData@distance))) {
    chrLoess <- mmapprData@distance[[i]]$loess
    chrDf <- data.frame(pos = chrLoess$x + tailPos, seqname = names(mmapprData@distance)[i], 
                        fitted = chrLoess$fitted, unfitted = chrLoess$y)
    plotDf <- rbind(plotDf, chrDf)
    tailPos <- chrDf$pos[nrow(chrDf)]
    breaks <- c(breaks, tailPos)
    labelpos <- c(labelpos, (breaks[length(breaks)] + breaks[length(breaks) - 1]) / 2)
  }
  
  pdf("mmappr_results/genome_plots.pdf", width=11, height=8.5)
  par(mfrow=c(2,1))
  plot(x = plotDf$pos, y = plotDf$fitted, type='l', 
       ylim=c(min(plotDf$fitted, na.rm=T),
              1.1*max(plotDf$fitted, na.rm=T)), 
       ylab=substitute("ED"^p~ ~"(Loess fit)", list(p=mmapprData@input$distancePower)), 
       xaxt='n', xaxs='i', xlab="Chromosome" )
  abline(v=(breaks[1:length(breaks)-1]+2), col="grey")
  mtext(unique(gsub("chr", "", plotDf$seqname[!is.na(plotDf$seqname)])), 
        at = labelpos, side=1, cex=.6)
  
  plot(x = plotDf$pos, y = plotDf$unfitted, pch=16, cex=.8, col="#999999AA", 
       ylim=c(min(plotDf$unfitted, na.rm=T), 
              1.1*max(plotDf$unfitted, na.rm=T)), 
       ylab=substitute("ED"^p, list(p=mmapprData@input$distancePower)), 
       xaxt='n', xaxs='i', xlab="Chromosome" )
  abline(v=(breaks), col="grey")
  mtext(unique(gsub("chr", "", plotDf$seqname[!is.na(plotDf$seqname)])), 
        at = labelpos, side=1, cex=.6)
  
  dev.off()
}

PlotPeaks <- function(mmapprData) {
  pdf("mmappr_results/peak_plots.pdf", width=11, height=8.5)
  par(mfrow=c(2,1))
  
  for (seqname in names(mmapprData@peaks)){
    chrLoess <- mmapprData@distance[[seqname]]$loess
    start <- mmapprData@peaks[[seqname]]$start
    end <- mmapprData@peaks[[seqname]]$end
    
    plot(chrLoess$x/1000000, chrLoess$fitted, type='l', 
         ylim=c(min(chrLoess$fitted, na.rm=T),
                1.1*max(chrLoess$fitted, na.rm=T)), 
         ylab=substitute("ED"^p~ ~"(Loess fit)", list(p=mmapprData@input$distancePower)), 
         xlab=paste(seqname,"Base Position (MB)"), xaxs='i' )
    shadeX <- c(start, 
                 chrLoess$x[chrLoess$x >= start & chrLoess$x <= end & 
                               !is.na(chrLoess$fitted)],
                 end)
    shadeY <- c(-5, 
                 chrLoess$fitted[chrLoess$x >= start & chrLoess$x <= end & 
                                  !is.na(chrLoess$fitted)],
                 -5)
    polygon(shadeX/1000000, shadeY, col=rgb(0,1,1,.75), border=NA)

    plot(chrLoess$x/1000000, chrLoess$y, pch=16, cex=.6,
         ylim=c(min(chrLoess$y, na.rm=T),
                1.1*max(chrLoess$y, na.rm=T)),
         ylab=substitute("ED"^p, list(p=mmapprData@input$distancePower)),
         xlab=paste(seqname,"Base Position (MB)"), xaxs='i' )
  }
  
    dev.off()
}

WriteCandidateVCFs <- function(mmapprData) {
  for (seqname in names(mmapprData@candidates)) {
    candidateVCF <- asVCF(mmapprData@candidates[[seqname]])
    writeVcf(candidateVCF, paste0("mmappr_results", seqname, ".vcf"))
  }
}