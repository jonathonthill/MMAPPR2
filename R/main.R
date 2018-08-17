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
#' This describes what we're doing
#'
#' @param mmapprParam A \linkS4class{MmapprParam} object containing desired parameters
#'
#' @return A \linkS4class{MmapprData} object containing results
#' @export
#'
#' @examples
#' \dontrun{
#' vepFlags <- ensemblVEP::VEPFlags(
#'                  flags=list(species = "danio_rerio_merged", 
#'                             format = "vcf",
#'                             database = TRUE))
#' mmapprParam <- MmapprParam(GmapGenome("danio_rerio_10"),
#'                            "wild_type.sorted.bam",
#'                            "mutant.sorted.bam",
#'                            vepFlags)
#' mmapprData <- mmappr(mmapprParam)
#' }
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
        mmapprData <- readInFiles(mmapprData)
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
        
        .messageAndLog("Output files generated", oF)
        outputMmapprData(mmapprData)
        
        return(mmapprData)
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


.addBamFileIndex <- function(bf) {
    path <- bf$path
    index <- paste0(path, ".bai")
    if (!file.exists(index)){
        Rsamtools::indexBam(bf)
    }
    return(Rsamtools::BamFile(path, index = index))
}

