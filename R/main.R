#' @title Mutation Mapping Analysis Pipeline for Pooled RNA-Seq
#'
#' @name mmappr
#'
#' @usage MMAPPR2 is designed to map the causative mutation in a forward genetics
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
#'
#'     # Specify parameters:
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#'
#'     # Run pipeline:
#'     mmapprData <- mmappr(mmappr_param)
#' }
#' 
#' \dontrun{
#' ### Alternately, you can navigate the pipeline step by step.
#' ### This may be helpful for debugging.
#' md <- MmapprData(mmapprParam)
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
#'   
NULL

mmappr <- function(mmapprParam) {
    startTime <- Sys.time()
    message("\n")
    message("[/][/][/][/][/][/][/][/][/][/][/][/][/][/][/]")
    message("[\\][\\][\\][\\] Welcome to MMAPPR2! [\\][\\][\\][\\]")
    message("[/][/][/][/][/][/][/][/][/][/][/][/][/][/][/]\n")

    .checkDep('samtools')

    md <- mmapprData(mmapprParam)
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

        md  # return for use after block
    },
    error = function(e) {
        .messageAndLog(paste('ERROR:', e$message), oF)
        .messageAndLog("MmapprData object up to the failing step is returned.", oF)
        .messageAndLog(paste0("You can also recover this object ",
            "from 'mmappr_data.RDS' in the '", outputFolder(param(md)),
            "' output folder"), oF)
        md
    })

    endTime <- Sys.time()
        .messageAndLog(paste('\nEnd time:', endTime), oF)
        .messageAndLog(paste("MMAPPR2 runtime:", format(endTime - startTime)), oF)
        saveRDS(md, file.path(md@param@outputFolder, "mmappr_data.RDS"))
        .log('\nsessionInfo()', oF)
        .log(sessionInfo(), oF)
        
        # Automatically reorder the log
        .reorderLogFile(oF)
        
        return(md)
    }

.messageAndLog <- function(msg, outputFolder) {
    logFile <- file.path(outputFolder, 'mmappr2.log')
    if (!is.character(msg)) msg <- capture.output(msg)
    cat(msg, file = logFile, sep = '\n', append = TRUE)
    msg <- paste(msg, collapse = '\n')
    message(msg)
}


.log <- function(msg, outputFolder) {
    logFile <- file.path(outputFolder, 'mmappr2.log')
    if (!is.character(msg)) msg <- capture.output(msg)
    cat(msg, file = logFile, sep = '\n', append = TRUE)
}

.reorderLogFile <- function(outputFolder) {

    logFile <- file.path(outputFolder, "mmappr2.log")
    if (!file.exists(logFile)) {
        return(invisible(NULL))
    }

    lines <- readLines(logFile)

    # Peak Summary
    peak_start <- grep("^Peak Summary:", lines)
    peak_block <- character()
    peak_range <- integer(0)

    if (length(peak_start) > 0L) {
        dens_line <- grep("^Value of maximum density:", lines)
        dens_line <- dens_line[dens_line >= peak_start[1L]]

        if (length(dens_line) > 0L) {
            peak_end   <- dens_line[1L]
            peak_range <- peak_start[1L]:peak_end
            peak_block <- lines[peak_range]
        }
    }

    # End Time
    end_start  <- grep("^End time:", lines)
    tail_block <- character()
    tail_range <- integer(0)

    if (length(end_start) > 0L) {
        tail_range <- end_start[1L]:length(lines)
        tail_block <- lines[tail_range]
    }

    # Chronological Log
    remove_idx <- unique(c(peak_range, tail_range))
    remaining  <- lines[setdiff(seq_along(lines), remove_idx)]

    # Build the reordered file
    new_lines <- c(
        "===== MMAPPR2 SUMMARY =====",
        "",
        if (length(peak_block)) c(peak_block, "") else character(0L),
        if (length(tail_block)) c(tail_block, "") else character(0L),
        "===== FULL LOG BELOW (chronological) =====",
        remaining
    )

    # Overwrite the file so there aren't two logs created
    writeLines(new_lines, logFile)
    invisible(NULL)
}



.checkDep <- function(program) {
    if (Sys.which(program) == '' || is.null(Sys.which(program))) {
        stop(paste(program, 'dependency is not installed (or at least not in path).'))
    } else {
        return(TRUE)
    }
}
