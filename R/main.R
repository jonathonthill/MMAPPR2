require(doParallel)
require(foreach)
require(Rsamtools)

## load bam file using the which argument to ScanBamParam
# bamfile3 <- "STAR/11647X5_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
# bamfile2 <- "STAR/11647X3_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
# bamfile1 <- "STAR/11647X1_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
# bamfilem1 <- "STAR/11647X2_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
# bamfilem2 <- "STAR/11647X4_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
# bamfilem3 <- "STAR/11647X6_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"

bamfile1 <- BamFile("../shared/MMAPPRtests/zy13/8187X1_110524_SN141_0355_AD0D99ABXX_1Aligned_chr_fix.bam")
bamfilem1 <- BamFile("../shared/MMAPPRtests/zy13/8187X2_110524_SN141_0355_AD0D99ABXX_1Aligned_chr_fix.bam")

wt_list <- BamFileList(bamfile1)#, bamfile2, bamfile3))
mut_list <- BamFileList(bamfilem1)#, bamfilem2))#, bamfilem3))


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
        
        return(detectCores() - 2)
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

.runFunctionInParallel <- function(inputList, functionToRun, numCores, silent = FALSE,
                                  packages = c(), secondInput = NULL, thirdInput = NULL) {
    if (length(inputList) < numCores) numCores <- length(inputList)
    
    #cluster generation
    outfile <- if (silent) "/dev/null" else ""
    cl <- makeCluster(numCores, type = "SOCK", outfile = outfile)
    
    # register the cluster
    registerDoParallel(cl)
    
    tryCatch({
        # find number of parameters and perform calculation
        message(sprintf("Beginning parallel calculation with %i core(s)", numCores))
        if (!is.null(thirdInput)) {
            numParams <- 3
            thirdInput <- .repList(thirdInput, length(inputList))
            secondInput <- .repList(secondInput, length(inputList))
            resultList <- foreach(a = inputList, b = secondInput, c = thirdInput,
                                  .packages = packages) %dopar% functionToRun(a, b, c)
        }
        else if (!is.null(secondInput)) {
            numParams <- 2
            secondInput <- .repList(secondInput, length(inputList))
            resultList <- foreach(a = inputList, b = secondInput,
                                  .packages = packages) %dopar% functionToRun(a, b)
        }
        else {
            numParams <- 1
            resultList <- foreach(a = inputList,
                                  .packages = packages) %dopar% functionToRun(a)
        }
        
        names(resultList) <- names(inputList)
        
    }, finally = {
        stopCluster(cl)
        registerDoSEQ()
        invisible(gc)
    })
    
    return(resultList)
}

mmappr <- function(mmapprParam) {
    startTime <- proc.time()
    
    mmapprData = new("MmapprData", param = mmapprParam)
    
    if(dir.exists(mmapprParam@outputFolder)){
        cat(sprintf("Output folder %s has already been created", mmapprParam@outputFolder))
        answer <- character()
        while(answer != "y" & answer != "n"){
            answer <- readline("Would you like to overwrite previously created results folder (y\n)?  ")
        }
        if (answer == "n") {
            newOutputFolder <- readline("Please enter name for new output folder or press enter for default: ")
            if (is.na(newOutputFolder)) mmapprParam@outputFolder <- "DEFAULT"
            else mmapprParam@outputFolder <- newOutputFolder
        }
    }
    
    mmapprData <- tryCatch({
        mmapprData <- readInFiles(mmapprData)
        message("File reading successfully completed and Euclidian distance data generated")
        
        mmapprData <- loessFit(mmapprData)
        message("Loess regression successfully completed")
        
        mmapprData <- prePeak(mmapprData)
        mmapprData <- peakRefinement(mmapprData)
        if (length(mmapprData@peaks) > 0)
            message("Peak regions succesfully identified")
        else message("No peak regions identified")
        
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
        message("You can also recover the object after the latest successful step from 'mmappr_data.RDS' in the output folder")
        return(mmapprData)
    })
    
    
    
    message("Mmappr runtime:")
    print(proc.time() - startTime)
    .prepareOutputFolder(mmapprData)
    saveRDS(mmapprData, file.path(mmapprData@param@outputFolder, "mmappr_data.RDS"))
    return(mmapprData)
}

.addBamFileIndex <- function(bf) {
    path <- bf$path
    index <- paste0(path, ".bai")
    if (!file.exists(index)){
        indexBam(bf)
    }
    return(BamFile(path, index = index))
}

