#' Generate plots and tables from MMAPPR2 data
#'
#' @param mmapprData The \linkS4class{MmapprData} object to be output
#'
#' @return A \linkS4class{MmapprData} object after writing output files
#'   to the folder specified in the \code{outputFolder} slot of the
#'   \code{link{MmapprParam}} used.
#' @export
#'
#' @examples 
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)
#'         & Sys.which('vep') != '') {
#'     mmapprParam <- MmapprParam(refFasta = MMAPPR2data::goldenFasta(),
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio",
#'                                outputFolder = tempOutputFolder())
#' }
#' \dontrun{
#' md <- new('MmapprData', param = mmapprParam)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' 
#' outputMmapprData(postCandidatesMD)
#' }
outputMmapprData <- function(mmapprData) {
    stopifnot(is(mmapprData, "MmapprData"))
    
    if (!dir.exists(outputFolder(param(mmapprData)))) {
        mmapprData <- .prepareOutputFolder(mmapprData)
    }
    
    tryCatch({
        .plotGenomeDistance(mmapprData)
        .plotPeaks(mmapprData)
    }, finally={
        graphics.off()
    })
    
    .writeCandidateTables(mmapprData@candidates,
                          outputFolder(param(mmapprData)))  
}


.defaultOutputFolder <- function()
    paste0("mmappr2_", format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))


#' Generate temporary output folder
#' 
#' Conveniently creates a timestamp-named temporary directory 
#'
#' @return The path to the temporary directory
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)
#'         & Sys.which('vep') != '') {
#'     mmapprParam <- MmapprParam(refFasta = MMAPPR2data::goldenFasta(),
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio",
#'                                outputFolder = tempOutputFolder())
#' }
tempOutputFolder <- function() {
    file.path(tempdir(), .defaultOutputFolder())
}

            
.prepareOutputFolder <- function(mmapprData) {
    if (outputFolder(mmapprData@param) == 'DEFAULT')
        outputFolder(mmapprData@param) <- .defaultOutputFolder()
    
    if(dir.exists(outputFolder(mmapprData@param))){
        message(sprintf("Output folder %s has already been created",
                        outputFolder(param(mmapprData))))
        answer <- " "
        while(answer != "y" & answer != "n"){
            answer <- readline(
                "Would you like to overwrite previously created results folder (y/n)?\n")
        }
        if (answer == "n") {
            newOutputFolder <-
                readline("Please enter name for new output folder or press enter for default: ")
            if (is.na(newOutputFolder))
                outputFolder(mmapprData@param) <- .defaultOutputFolder()
            else outputFolder(mmapprData@param) <- newOutputFolder
            dir.create(outputFolder(mmapprData@param))
        } else {
            unlink(file.path(outputFolder(param(mmapprData)), '*'))
        }
    } else {
        dir.create(outputFolder(mmapprData@param), recursive=TRUE)
    }
    
    file.create(file.path(outputFolder(mmapprData@param), 'mmappr2.log'))
    
    return(mmapprData)
}


.plotGenomeDistance <- function(mmapprData, savePdf=TRUE) {
    #generate one big dataframe for plots, along with break and label points
    tailPos <- 0
    plotDf <- NULL
    breaks <- tailPos
    labelpos <- NULL
    for (i in GenomeInfoDb::orderSeqlevels(names(mmapprData@distance))) {
        if (!is(mmapprData@distance[[i]], 'list'))
            next
        else if (!("loess" %in% names(mmapprData@distance[[i]])))
            stop("Distance list for sequence %s missing loess fit data",
                 names(mmapprData@distance)[i])
        
        chrLoess <- mmapprData@distance[[i]]$loess
        chrDf <- data.frame(pos = chrLoess$x + tailPos,
                            seqname = names(mmapprData@distance)[i], 
                            fitted = chrLoess$fitted, unfitted = chrLoess$y)
        plotDf <- rbind(plotDf, chrDf)
        tailPos <- chrDf$pos[nrow(chrDf)]
        breaks <- c(breaks, tailPos)
        labelpos <-
            c(labelpos, (breaks[length(breaks)]+breaks[length(breaks)-1])/2)
    }
    
    if (savePdf)
        pdf(file.path(mmapprData@param@outputFolder, "genome_plots.pdf"),
            width=11, height=8.5)
    par(mfrow=c(2,1))
    plot(x = plotDf$pos, y = plotDf$fitted, type='l', 
         ylim=c(min(plotDf$fitted, na.rm=TRUE),
                1.1*max(plotDf$fitted, na.rm=TRUE)), 
         ylab=substitute("ED"^p~ ~"(Loess fit)",
                         list(p=mmapprData@param@distancePower)), 
         xaxt='n', xaxs='i', xlab="Chromosome" )
    abline(v=(breaks[seq_len(length(breaks))-1]+2), col="grey")
    mtext(unique(gsub("chr", "", plotDf$seqname[!is.na(plotDf$seqname)])), 
          at = labelpos, side=1, cex=.6)
    
    plot(x = plotDf$pos, y = plotDf$unfitted, pch=16, cex=.8, col="#999999AA", 
         ylim=c(min(plotDf$unfitted, na.rm=TRUE), 
                1.1*max(plotDf$unfitted, na.rm=TRUE)), 
         ylab=substitute("ED"^p, list(p=mmapprData@param@distancePower)), 
         xaxt='n', xaxs='i', xlab="Chromosome" )
    abline(v=(breaks), col="grey")
    mtext(unique(gsub("chr", "", plotDf$seqname[!is.na(plotDf$seqname)])), 
          at = labelpos, side=1, cex=.6)
    
    if (savePdf) dev.off()
}


.plotPeaks <- function(mmapprData) {
    pdf(file.path(mmapprData@param@outputFolder, "peak_plots.pdf"),
        width=11, height=8.5)
    par(mfrow=c(2,1), mar=c(5, 4, 4, 4))
    
    for (seqname in names(mmapprData@peaks)){
        chrLoess <- mmapprData@distance[[seqname]]$loess
        start <- mmapprData@peaks[[seqname]]$start
        end <- mmapprData@peaks[[seqname]]$end
        densityData <- mmapprData@peaks[[seqname]]$densityData
        
        plot(chrLoess$x/1000000, chrLoess$fitted, type='l', 
             ylim=c(min(chrLoess$fitted, na.rm=TRUE),
                    1.1*max(chrLoess$fitted, na.rm=TRUE)),
             xlim=c(chrLoess$x[1], tail(chrLoess$x, n=1)) * 1E-6,
             xlab=paste(seqname,"Base Position (MB)"),
             ylab=NA, xaxs='i')
        mtext(substitute("ED"^p~ ~"(Loess fit)",
                         list(p=mmapprData@param@distancePower)),
              side=2, line=2)
        shadeX <- c(start, 
                    chrLoess$x[chrLoess$x >= start & chrLoess$x <= end & 
                                   !is.na(chrLoess$fitted)],
                    end)
        shadeY <- c(-5, 
                    chrLoess$fitted[chrLoess$x >= start & chrLoess$x <= end & 
                                        !is.na(chrLoess$fitted)],
                    -5)
        polygon(shadeX/1000000, shadeY, col='#2ecc71', border=NA)
        
        # overlay density curve
        posMatch <- (densityData$x >= chrLoess$x[1]) &
            densityData$x <= chrLoess$x[length(chrLoess$x)]
        densityData$x <- densityData$x[posMatch]
        densityData$y <- densityData$y[posMatch]
        par(new=TRUE)
        plot(densityData, type='l',
             ylim=c(min(densityData$y, na.rm=TRUE),
                    1.1*max(densityData$y, na.rm=TRUE)),
             xlim=c(chrLoess$x[1], tail(chrLoess$x, n=1)),
             ann=FALSE, xaxs='i',
             xaxt='n', yaxt='n', col='#502ecc')
        axis(side=4, col='#502ecc')
        mtext(side=4, line=2, 'Probability', col='#502ecc')
        legend('topright',
               legend=c("Fitted Distance Curve",
                        "Peak Resampling Distribution"),
               col=c("black", "#502ecc"), lty=c(1, 1), cex=0.8, lwd=3)
        
        # SNPs plot
        plot(chrLoess$x/1000000, chrLoess$y, pch=16, cex=.6,
             ylim=c(min(chrLoess$y, na.rm=TRUE),
                    1.1*max(chrLoess$y, na.rm=TRUE)),
             ylab=substitute("ED"^p, list(p=mmapprData@param@distancePower)),
             xlab=paste(seqname,"Base Position (MB)"), xaxs='i' )
    }
    
    dev.off()
}


.writeCandidateTables <- function(candList, outputFolder){
    for (seqname in names(candList)) {
        listData <- candList[[seqname]]@elementMetadata@listData
        output <- 
            data.frame(
                "Position"=candList[[seqname]]@ranges@start, 
                "Symbol"=listData$SYMBOL, 
                "Impact"=listData$IMPACT,
                "Consequence"=listData$Consequence,
                "DensityScore"=listData$peakDensity,
                "Allele"=listData$Allele,
                "AminoAcid"=listData$Amino_acids,
                "Feature"=listData$Feature
            )
        filename <- paste0(seqname, '.tsv')
        write.table(output, file=file.path(outputFolder, filename),
                    sep='\t', row.names=FALSE, quote=FALSE)
    }
}
