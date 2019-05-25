#' Mutation Mapping Analysis Pipeline for Pooled RNA-Seq
#' 
#' MMAPPR2 is designed to map the causative mutation in a forward genetics
#' screen. It analyzes aligned sequence files, calculates the per-base
#' Euclidean distance between the mutant and wild-type pools, performs
#' a Loess regression on that distance, and generates candidate variants
#' in regions of peak distance.
#'
#' @param mmapprParam A \code{\linkS4class{MmapprParam}} object containing
#'   desired parameters.
#'
#' @return A \code{\linkS4class{MmapprData}} object containing results
#'   and/or intermediate data.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE) &
#'         Sys.which('vep') != '') {
#'     slc24a5genome <- gmapR::GmapGenome(MMAPPR2data::goldenFasta(),
#'                                        name = 'slc24a5',
#'                                        create = TRUE)
#'
#'     # Specify parameters:
#'     mmapprParam <- MmapprParam(refGenome = slc24a5genome,
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio")
#'                                
#'     # Run pipeline:
#'     mmapprData <- mmappr(mmapprParam)
#' 
#' }
#' \dontrun{
#' ### Alternately, you can navigate the pipeline step by step.
#' ### This may be helpful for debugging.
#' md <- new('MmapprData', param = mmapprParam)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' outputMmapprData(postCandidatesMD)
#' } 
#' 
#' @seealso \code{\link{calculateDistance}}, \code{\link{loessFit}},
#'   \code{\link{prePeak}}, \code{\link{peakRefinement}},
#'   \code{\link{generateCandidates}}, \code{\link{outputMmapprData}}
mmappr <- function(mmapprParam) {
    startTime <- Sys.time()
    message("------------------------------------")
    message("-------- Welcome to MMAPPR2 --------")
    message("------------------------------------\n")
    
    .checkDep('samtools')
    
    md <- new("MmapprData", param=mmapprParam)
    md <- .prepareOutputFolder(md)
    oF <- outputFolder(md@param)
    .messageAndLog(paste('Start time:', Sys.time()), oF)
    .messageAndLog(paste('Output folder:', file.path(getwd(), oF), '\n'), oF)
    .log('Parameters:', oF)
    .log(mmapprParam, oF)
    .log('', oF)
    
    md <- tryCatch({
        .messageAndLog("Reading BAM files and generating Euclidean distance data...", oF)
        md <- calculateDistance(md)
        .messageAndLog("Generating optimal Loess curves for each chromosome...", oF)
        md <- loessFit(md)
        .messageAndLog("Identifying chromosome(s) harboring linkage region...", oF)
        md <- prePeak(md)
        if (length(peaks(md)) > 0)
            .messageAndLog("Peak regions succesfully identified", oF)
        else stop("No peak regions identified")
        .messageAndLog("Refining peak characterization using SNP resampling...", oF)
        md <- peakRefinement(md)
        .messageAndLog('Generating, analyzing, and ranking candidate variants...', oF)        
        md <- generateCandidates(md)
        .messageAndLog("Writing output plots and tables...", oF)
        outputMmapprData(md)
        
        return(md)  # return for use after block
    }, 
    error = function(e) {
        .messageAndLog(paste('ERROR:', e$message), oF)
        .messageAndLog("MmapprData object up to the failing step is returned.", oF)
        .messageAndLog(paste0("You can also recover this object ",
            "from 'mmappr_data.RDS' in the '", outputFolder(param(md)),
            "' output folder"), oF)
        return(md)
    })
    
    endTime <- Sys.time()
    .messageAndLog(paste('\nEnd time:', endTime), oF)
    .messageAndLog(paste("MMAPPR2 runtime:", format(endTime - startTime)), oF)
    saveRDS(md, file.path(md@param@outputFolder, "mmappr_data.RDS"))
    .log('\nsessionInfo()', oF)
    .log(sessionInfo(), oF)
    
    return(md)
}


.messageAndLog <- function(msg, outputFolder) {
    logFile <- file.path(outputFolder, 'mmappr2.log')
    if (!is.character(msg)) msg <- capture.output(msg)
    cat(msg, file=logFile, sep='\n', append=TRUE)
    msg <- paste(msg, collapse='\n')
    message(msg)
}


.log <- function(msg, outputFolder) {
    logFile <- file.path(outputFolder, 'mmappr2.log')
    if (!is.character(msg)) msg <- capture.output(msg)
    cat(msg, file=logFile, sep='\n', append=TRUE)
}


.checkDep <- function(program) {
    if (Sys.which(program) == '' || is.null(Sys.which(program))) {
        stop(program + ' dependency is not installed.')
    }
}