#' @title Generate plots and tables from MMAPPR2 data
#' 
#' @name outputMmapprData
#' 
#' @param mmapprData The \linkS4class{MmapprData} object to be output
#'
#' @return A \linkS4class{MmapprData} object after writing output files
#'   to the folder specified in the \code{outputFolder} slot of the
#'   \code{link{MmapprParam}} used.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#' 
#' \dontrun{
#' md <- mmapprData(mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#'
#' outputMmapprData(postCandidatesMD)
#' }
#' 
NULL


outputMmapprData <- function(mmapprData) {
  stopifnot(is(mmapprData, "MmapprData"))
  
  if (!dir.exists(outputFolder(param(mmapprData)))) {
    mmapprData <- .prepareOutputFolder(mmapprData)
  }
  current_devs = dev.list() # get open graphics devices to help cleanup
  if (length(snpDistance(mmapprData)) > 0) {
    tryCatch({
      .plotGenomeDistance(mmapprData)
      .plotPeaks(mmapprData)
    }, finally = {
      opened_devs <- dev.list()[!(dev.list() %in% current_devs)]
      if (length(opened_devs) > 0) dev.off(opened_devs)
    })
  }
  
  if (length(candidates(mmapprData)) > 0) {
    .writeCandidateTables(mmapprData@candidates,
                          outputFolder(param(mmapprData)))
  }
  
  invisible(1)
}


.defaultOutputFolder <- function()
  paste0("mmappr2_", format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))


#' @title Generate temporary output folder
#' @name tempOutputFolder
#' @usage Conveniently creates a timestamp-named temporary directory
#'
#' @return The path to the temporary directory
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- MmapprParam(refFasta = MMAPPR2data::goldenFasta(),
#'                                 wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 species = "danio_rerio",
#'                                 outputFolder = tempOutputFolder())
#' }
NULL


tempOutputFolder <- function() {
  file.path(tempdir(), .defaultOutputFolder())
}


.prepareOutputFolder <- function(outputFolder) {
  if(dir.exists(outputFolder)){
    message(sprintf("Output folder %s has already been created",
                    outputFolder))
    answer <- readline(
        "Would you like to overwrite previously created results folder (Y/n)?\n")
    if (answer == "n") {
      newOutputFolder <-
        readline("Please enter name for new output folder or press enter for default: ")
      if (is.na(newOutputFolder))
        outputFolder <- .defaultOutputFolder()
      else outputFolder <- newOutputFolder
      dir.create(outputFolder)
    } else {
      unlink(file.path(outputFolder, '*'))
    }
  } else {
    dir.create(outputFolder, recursive = TRUE)
  }
  
  file.create(file.path(outputFolder, 'mmappr2.log'))
  
  return(outputFolder)
}


.plotGenomeDistance <- function(mmapprData, savePdf = TRUE) {
  #generate one big dataframe for plots, along with break and label points
  tailPos <- 0
  plotDf <- NULL
  breaks <- tailPos
  labelpos <- NULL
  for (i in orderSeqlevels(names(mmapprData@snpDistance))) {
    if (!is(mmapprData@snpDistance[[i]], 'list'))
      next
    else if (!("loess" %in% names(mmapprData@snpDistance[[i]])))
      stop("Distance list for sequence %s missing loess fit data",
           names(mmapprData@snpDistance)[i])
    
    chrLoess <- mmapprData@snpDistance[[i]]$loess
    chrDf <- data.frame(pos = chrLoess$x + tailPos,
                        seqname = names(mmapprData@snpDistance)[i],
                        fitted = chrLoess$fitted, unfitted = chrLoess$y)
    plotDf <- rbind(plotDf, chrDf)
    tailPos <- chrDf$pos[nrow(chrDf)]
    breaks <- c(breaks, tailPos)
    labelpos <-
      c(labelpos, (breaks[length(breaks)]+breaks[length(breaks)-1])/2)
  }
  
  if (savePdf)
    pdf(file.path(mmapprData@param@outputFolder, "genome_plots.pdf"),
        width = 11, height = 8.5)
  par(mfrow = c(2,1))
  
  plot(x = plotDf$pos, y = plotDf$fitted, type = 'l',
       ylim = c(min(plotDf$fitted, na.rm = TRUE),
                1.1*max(plotDf$fitted, na.rm = TRUE)),
       ylab = NA,
       xaxt = 'n', xaxs = 'i', xlab = "Chromosome")
  abline(v = (breaks[seq_len(length(breaks))-1]+2), col = "grey")
  mtext(unique(gsub("chr", "", plotDf$seqname[!is.na(plotDf$seqname)])),
        at = labelpos, side = 1, cex = .6)
  mtext(substitute("ED"^p~ ~"(Loess fit)",
                   list(p = mmapprData@param@distancePower)), 
        side = 2, line = 2)
  
  plot(x = plotDf$pos, y = plotDf$unfitted, 
       pch = 16, cex = .8, col = "#999999AA",
       ylim = c(min(plotDf$unfitted, na.rm = TRUE),
                1.1*max(plotDf$unfitted, na.rm = TRUE)),
       ylab = NA,
       xaxt = 'n', xaxs = 'i', xlab = "Chromosome" )
  abline(v = (breaks), col = "grey")
  mtext(unique(gsub("chr", "", plotDf$seqname[!is.na(plotDf$seqname)])),
        at = labelpos, side = 1, cex = .6)
  mtext(substitute("ED"^p~ ~"(Loess fit)", 
                   list(p = mmapprData@param@distancePower)),
        side = 2, line = 2)
  
  if (savePdf) dev.off()
}


.plotPeaks <- function(mmapprData) {
  pdf(file.path(mmapprData@param@outputFolder, "peak_plots.pdf"),
      width = 11, height = 8.5)
  par(mfrow = c(2,1), mar = c(5, 4, 4, 4))
  
  for (seqname in names(mmapprData@peaks)){
    chrLoess <- mmapprData@snpDistance[[seqname]]$loess
    start <- mmapprData@peaks[[seqname]]$start
    end <- mmapprData@peaks[[seqname]]$end
    densityData <- mmapprData@peaks[[seqname]]$densityData
    
    plot(chrLoess$x/1000000, chrLoess$fitted, type = 'l',
         ylim = c(min(chrLoess$fitted, na.rm = TRUE),
                  1.1*max(chrLoess$fitted, na.rm = TRUE)),
         xlim = c(chrLoess$x[1], tail(chrLoess$x, n = 1)) * 1E-6,
         xlab = paste(seqname,"Base Position (MB)"),
         ylab = NA, xaxs = 'i')
    mtext(substitute("ED"^p~ ~"(Loess fit)",
                     list(p = mmapprData@param@distancePower)),
          side = 2, line = 2)
    shadeX <- c(start,
                chrLoess$x[chrLoess$x >= start & chrLoess$x <= end &
                             !is.na(chrLoess$fitted)],
                end)
    shadeY <- c(-5,
                chrLoess$fitted[chrLoess$x >= start & chrLoess$x <= end &
                                  !is.na(chrLoess$fitted)],
                -5)
    polygon(shadeX/1000000, shadeY, col = '#2ecc71', border = NA)
    
    # overlay density curve
    if (length(densityData) > 0) {
      posMatch <- (densityData$x >= chrLoess$x[1]) &
        densityData$x <= chrLoess$x[length(chrLoess$x)]
      densityData$x <- densityData$x[posMatch]
      densityData$y <- densityData$y[posMatch]
      par(new = TRUE)
      plot(densityData, type = 'l',
           ylim = c(min(densityData$y, na.rm = TRUE),
                    1.1*max(densityData$y, na.rm = TRUE)),
           xlim = c(chrLoess$x[1], tail(chrLoess$x, n = 1)),
           ann = FALSE, xaxs = 'i',
           xaxt = 'n', yaxt = 'n', col = '#502ecc')
      axis(side = 4, col = '#502ecc')
      mtext(side = 4, line = 2, 'Probability', col = '#502ecc')
      legend('topright',
             legend = c("Fitted Distance Curve",
                        "Peak Resampling Distribution"),
             col = c("black", "#502ecc"), lty = c(1, 1), cex = 0.8, lwd = 3)
    }
    
    # SNPs plot
    if (length(chrLoess) > 0) {
      plot(chrLoess$x/1000000, chrLoess$y, pch = 16, cex = .6,
           ylim = c(min(chrLoess$y, na.rm = TRUE),
                    1.1*max(chrLoess$y, na.rm = TRUE)),
           ylab = NA,
           xlab = paste(seqname,"Base Position (MB)"), xaxs = 'i' )
      mtext(substitute("ED"^p, list(p = mmapprData@param@distancePower)), 
            side = 2, line = 2)
    }
  }
  dev.off()
}


.writeCandidateTables <- function(candList, outputFolder){
  for (seqname in names(candList$effects)) {
    listData <- as.data.frame(candList$effects[[seqname]])
    filename <- paste0("DetectedMutationsFor", seqname, '.tsv')
    listData$PROTEINLOC <- sapply(listData$PROTEINLOC, function(x) {x[1]})
    outdata <- listData[, c("seqnames", "start", "end", "width", "strand", "ref", "alt", "refDepth", "altDepth", "peakDensity", "GENEID", "TXID", "PROTEINLOC", "CONSEQUENCE", "REFAA", "VARAA")]
    write.table(x = outdata, 
                file = file.path(outputFolder, filename), 
                sep = '\t', quote = FALSE, row.names = FALSE)
    filename <- paste0("DifferentiallyExpressedGenesFor", seqname, '.tsv')
    write.table(x = candList$diff[[seqname]], 
                file = file.path(outputFolder, filename), 
                sep = '\t', quote = FALSE, row.names = FALSE)
  }
}