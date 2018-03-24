## load bam file using the which argument to ScanBamParam

bamfile1 <- Rsamtools::BamFile("../shared/mmappr_tests/zy14/hisat2-zv9/8187X3_110524_SN141_0355_AD0D99ABXX_2.hisat2.sorted.fixed.bam")
bamfilem1 <- Rsamtools::BamFile("../shared/mmappr_tests/zy14/hisat2-zv9/8187X4_110524_SN141_0355_AD0D99ABXX_2.hisat2.sorted.fixed.bam")

wt_list <- Rsamtools::BamFileList(bamfile1)#, bamfile2, bamfile3))
mut_list <- Rsamtools::BamFileList(bamfilem1)#, bamfilem2))#, bamfilem3))


#automatic core number generation (only on Linux)
#mem_req is memory required per base for distance list item, as returned in ReadInFiles
.coreCalc <- function(mem_req = 33.3){
    tryCatch({
        # #this gets free RAM (on Linux)
        # free_ram <- system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE) #gets mem info in MB
        # #get numeric, in bytes (not MB)
        # free_ram <- 1000 * as.numeric(free_ram)
        # #get longest chromosome, calculate needed memory, and divide into free RAM
        # num_cores <- ceiling(free_ram / (max(seqlengths(myrange)) * mem_req))
        # message("num_cores = ", num_cores)
        # return(num_cores)
        
        return(if (parallel::detectCores() > 2) parallel::detectCores() - 2 else 1)
    },
    warning = function(w) {
        message('Warning: arbitrary number of clusters used (are you not on Linux?)')
        return(ceiling(detectCores() / 3))
    }
    )
}

.repList <- function(object, times) {
    stopifnot(times >= 0)
    result <- list()
    if (times == 0) return(result)
    for (i in 1:times) {
        result <- append(result, object)
    }
    return(result)
}

.runFunctionInParallel <- function(inputList, functionToRun, numCores, silent=FALSE,
                                  packages=c(), secondInput=NULL, thirdInput=NULL) {
    require(doParallel, quietly=TRUE)
    if (length(inputList) < numCores) numCores <- length(inputList)
    
    
    #cluster generation
    if (numCores > 1) {
        outfile <- if (silent) "/dev/null" else ""
        cl <- parallel::makeCluster(numCores, type='SOCK', outfile=outfile)
        #clusterEvalQ(cl, devtools::dev_mode()) #TODO: for development only
        
        
        # register the cluster
        doParallel::registerDoParallel(cl)
    }
    
    
    
    tryCatch({
        # find number of parameters and perform calculation
        message(sprintf("Beginning parallel calculation with %i core(s)", numCores))
        if (!is.null(thirdInput)) {
            numParams <- 3
            thirdInput <- .repList(thirdInput, length(inputList))
            secondInput <- .repList(secondInput, length(inputList))
            if (numCores > 1)
                resultList <- foreach(a = inputList, b = secondInput, c = thirdInput,
                                  .packages = packages) %dopar% functionToRun(a, b, c)
            else
                resultList <- 
                    lapply(inputList, 
                           function(x) functionToRun(x, secondInput, thirdInput))
        }
        else if (!is.null(secondInput)) {
            numParams <- 2
            secondInput <- .repList(secondInput, length(inputList))
            if (numCores > 1)
                resultList <- foreach(a = inputList, b = secondInput,
                                  .packages = packages) %dopar% functionToRun(a, b)
            else
                resultList <- 
                    lapply(inputList, function(x) functionToRun(x, secondInput))
        }
        else {
            numParams <- 1
            if (numCores > 1)
                resultList <- foreach(a = inputList,
                                  .packages = packages) %dopar% functionToRun(a)
            else
                resultList <- lapply(inputList, functionToRun)
        }
        
        names(resultList) <- names(inputList)
        
    }, finally = {
        if (numCores > 1) {
            stopCluster(cl)
            foreach::registerDoSEQ()
            invisible(gc)
        }
    })
    
    return(resultList)
}



#' Title
#'
#' @param mmapprParam 
#' @param showDebug 
#'
#' @return
#' @export
#'
#' @examples
mmappr <- function(mmapprParam, showDebug=FALSE) {
    startTime <- proc.time()
    
    mmapprData <- new("MmapprData", param=mmapprParam)
    mmapprData <- .prepareOutputFolder(mmapprData)
    
    mmapprData <- tryCatch({
        mmapprData <- readInFiles(mmapprData, showDebug=showDebug)
        message("File reading successfully completed and Euclidian distance data generated")
        
        mmapprData <- loessFit(mmapprData)
        message("Loess regression successfully completed")
        
        mmapprData <- prePeak(mmapprData)
        mmapprData <- peakRefinement(mmapprData)
        if (length(mmapprData@peaks) > 0)
            message("Peak regions succesfully identified")
        else {
            error("No peak regions identified")
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
        indexBam(bf)
    }
    return(Rsamtools::BamFile(path, index = index))
}

