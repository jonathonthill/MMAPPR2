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
#' @param showDebug If set to \code{TRUE}, intermediate dataframe sizes are
#'     displayed during file reading step
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
    startTime <- proc.time()
    message("------------------------------------")
    message("-------- Welcome to MMAPPR2 --------")
    message("------------------------------------")
    
    mmapprData <- new("MmapprData", param=mmapprParam)
    mmapprData <- .prepareOutputFolder(mmapprData)
    
    mmapprData <- tryCatch({
        mmapprData <- readInFiles(mmapprData)
        message("File reading successfully completed and Euclidian distance data generated")
        
        mmapprData <- loessFit(mmapprData)
        message("Loess regression successfully completed")
        
        mmapprData <- prePeak(mmapprData)
        mmapprData <- peakRefinement(mmapprData)
        if (length(mmapprData@peaks) > 0)
            message("Peak regions succesfully identified")
        else {
            stop("No peak regions identified")
        }
        
        mmapprData <- generateCandidates(mmapprData)
        if (length(mmapprData@candidates) > 0)
            message("Candidate variants generated, analyzed, and ranked")
        
        outputMmapprData(mmapprData)
        message("Output PDF files generated")
        
        return(mmapprData)
    }, 
    error = function(e) {
        traceback()
        message(e$message)
        message("MmapprData object is returned up until the failing step")
        message(paste0("You can also recover the object after the latest successful step from 'mmappr_data.RDS' in the '", 
                outputFolder(param(mmapprData)),
                "' output folder"))
        return(mmapprData)
    })
    
    
    
    message("Mmappr runtime:")
    print(proc.time() - startTime)
    saveRDS(mmapprData, file.path(mmapprData@param@outputFolder, "mmappr_data.RDS"))
    return(mmapprData)
}

.addBamFileIndex <- function(bf) {
    path <- bf$path
    index <- paste0(path, ".bai")
    if (!file.exists(index)){
        Rsamtools::indexBam(bf)
    }
    return(Rsamtools::BamFile(path, index = index))
}

