.runFunctionInParallel <- function(inputList, functionToRun, ..., BPPARAM=NULL) {
    if (is.null(BPPARAM)) bpParam <- BiocParallel::bpparam()
    BiocParallel::bpprogressbar(bpParam) <- TRUE
    resultList <- BiocParallel::bplapply(inputList, functionToRun, ...,
                                    BPPARAM=bpParam)
    BiocParallel::bpstop(bpParam)
    return(resultList)
}


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
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     ## Ignore these lines:
#'     MMAPPR2:::.insertFakeVEPintoPath()
#'     genDir <- gmapR::GmapGenomeDirectory('example', create=TRUE)
#' 
#'     # Specify parameters:
#'     mmapprParam <- MmapprParam(refGenome = gmapR::GmapGenome("GRCz11", genDir),
#'                                wtFiles = MMAPPR2data::zy13wtBam(),
#'                                mutFiles = MMAPPR2data::zy13mutBam(),
#'                                species = "danio_rerio")
#'                                
#'     # Run pipeline:
#'     mmapprData <- mmappr(mmapprParam)
#' 
#'     ### Alternately, you can navigate the pipeline step by step.
#'     ### This may be helpful for debugging.
#'     md <- new('MmapprData', param = mmapprParam)
#'     postCalcDistMD <- calculateDistance(md)
#'     postLoessMD <- loessFit(postCalcDistMD)
#'     postPrePeakMD <- prePeak(postLoessMD)
#'     postPeakRefMD <- peakRefinement(postPrePeakMD)
#'     \dontrun{postCandidatesMD <- generateCandidates(postPeakRefMD)}
#'     outputMmapprData(postCandidatesMD)
#' }
#' 
#' @seealso \code{\link{calculateDistance}}, \code{\link{loessFit}},
#'   \code{\link{prePeak}}, \code{\link{peakRefinement}},
#'   \code{\link{generateCandidates}}, \code{\link{outputMmapprData}}
mmappr <- function(mmapprParam) {
    startTime <- Sys.time()
    message("------------------------------------")
    message("-------- Welcome to MMAPPR2 --------")
    message("------------------------------------")
    message('')
    
    mmapprData <- new("MmapprData", param=mmapprParam)
    mmapprData <- .prepareOutputFolder(mmapprData)
    oF <- outputFolder(mmapprData@param)
    .messageAndLog(paste('Start time:', Sys.time()), oF)
    .messageAndLog(paste('Output folder:', file.path(getwd(), oF), '\n'), oF)
    
    .log('Parameters:', oF)
    .log(mmapprParam, oF)
    .log('', oF)
    
    mmapprData <- tryCatch({
        .messageAndLog("Reading BAM files and generating Euclidean distance data...", oF)
        mmapprData <- calculateDistance(mmapprData)
        .messageAndLog("Done", oF)
        
        .messageAndLog("Generating optimal Loess curves for each chromosome...", oF)
        mmapprData <- loessFit(mmapprData)
        .messageAndLog("Done", oF)
        
        .messageAndLog("Identifying chromosome(s) harboring linkage region...", oF)
        mmapprData <- prePeak(mmapprData)
        if (length(mmapprData@peaks) > 0)
            .messageAndLog("Peak regions succesfully identified", oF)
        else {
            stop("No peak regions identified")
        }
        .messageAndLog("Refining peak characterization using SNP resampling...", oF)
        mmapprData <- peakRefinement(mmapprData)
        .messageAndLog("Done", oF)

        .messageAndLog('Generating, analyzing, and ranking candidate variants...', oF)        
        mmapprData <- generateCandidates(mmapprData)
        .messageAndLog("Done", oF)
        
        .messageAndLog("Writing output plots and tables...", oF)
        outputMmapprData(mmapprData)
        .messageAndLog("Done", oF)
        
        mmapprData  # return for use after block
    }, 
    error = function(e) {
        .messageAndLog(paste('ERROR:', e$message), oF)
        .messageAndLog("MmapprData object up to the failing step is returned.", oF)
        .messageAndLog(paste0("You can also recover this object ",
                              "from 'mmappr_data.RDS' in the '", 
                              outputFolder(param(mmapprData)),
                              "' output folder"), oF)
        return(mmapprData)
    })
    
    endTime <- Sys.time()
    .messageAndLog(paste('\nEnd time:', endTime), oF)
    runtime <- format(endTime - startTime)
    .messageAndLog(paste("MMAPPR2 runtime:", runtime), oF)
    saveRDS(mmapprData, file.path(mmapprData@param@outputFolder, "mmappr_data.RDS"))
    
    .log('\nsessionInfo()', oF)
    .log(sessionInfo(), oF)
    
    return(mmapprData)
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


.insertFakeVEPintoPath <- function() {
    unlink('/tmp/bin', recursive=TRUE)
    dir.create('/tmp/bin')
    file.create('/tmp/bin/vep')
    system('chmod 777 /tmp/bin/vep')
    Sys.setenv('PATH'=paste(Sys.getenv('PATH'), '/tmp/bin', sep=':'))
}