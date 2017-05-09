#source function definitions or have just one document
require(doParallel)
require(foreach)

## load bam file using the which argument to ScanBamParam
bamfile3 <- "STAR/11647X5_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfile2 <- "STAR/11647X3_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfile1 <- "STAR/11647X1_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfilem1 <- "STAR/11647X2_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfilem2 <- "STAR/11647X4_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfilem3 <- "STAR/11647X6_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"

bamfile1 <- BamFile("../shared/MMAPPRtests/zy13/8187X1_110524_SN141_0355_AD0D99ABXX_1Aligned_chr_fix.bam")
bamfilem1 <- BamFile("../shared/MMAPPRtests/zy13/8187X2_110524_SN141_0355_AD0D99ABXX_1Aligned_chr_fix.bam")

wt_list <- BamFileList(bamfile1)#, bamfile2, bamfile3))
mut_list <- BamFileList(bamfilem1)#, bamfilem2))#, bamfilem3))


#automatic core number generation (only on Linux)
#mem_req is memory required per base for distance list item, as returned in ReadInFiles
CoreCalc <- function(mem_req = 33.3){
  tryCatch({
    #this gets free RAM (on Linux)
    free_ram <- system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE) #gets mem info in MB
    #get numeric, in bytes (not MB)
    free_ram <- 1000 * as.numeric(free_ram)
    #get longest chromosome, calculate needed memory, and divide into free RAM
    num_cores <- ceiling(free_ram / (max(seqlengths(myrange)) * mem_req))
    message("num_cores = ", num_cores)
    return(num_cores)
    },
  warning = function(w) {
    message('Warning: arbitrary number of clusters used (are you not on Linux?)')
    return(ceiling(detectCores() / 3))
  }
  )
}

RunFunctionInParallel <- function(inputList, functionToRun, packages = c(), secondInput = NULL, thirdInput = NULL) {
  #cluster generation
  cl <- makeCluster(CoreCalc(), type = "SOCK", outfile = "")
  
  # register the cluster
  registerDoParallel(cl)
  
  # find number of parameters and perform calculation
  message("Beginning parallel calculation")
  if (!is.null(thirdInput)) {
    numParams <- 3
    thirdInput <- rep(thirdInput, length(inputList))
    secondInput <- rep(secondInput, length(inputList))
    resultList <- foreach(a = inputList, b = secondInput, c = thirdInput,
                          .packages = packages) %dopar% functionToRun(a, b, c)
  }
  else if (!is.null(secondInput)) {
    numParams <- 2
    secondInput <- rep(secondInput, length(inputList))
    resultList <- foreach(a = inputList, b = secondInput,
                          .packages = packages) %dopar% functionToRun(a, b)
  }
  else {
    numParams <- 1
    resultList <- foreach(a = inputList,
                          .packages = packages) %dopar% functionToRun(a)
  }
  
  
  names(resultList) <- names(inputList)
  
  stopCluster(cl)
  registerDoSEQ()
  
  invisible(gc)
  
  return(resultList)
}

Mmappr <- function(mmapprData = NULL, wtFiles = NULL, mutFiles = NULL) {
  startTime <- proc.time()
  
  if (is.null(mmapprData)) mmapprData = new("MmapprData")
  if (!is.null(wtFiles)) mmapprData@input$wtFiles <- wtFiles
  if (!is.null(mutFiles)) mmapprData@input$mutFiles <- mutFiles
  
  mmapprData <- ReadInFiles(mmapprData)
  mmapprData <- LoessFit(mmapprData)
  # mmapprData <- PrePeak(mmapprData)
  # # mmapprData <- PeakDensity(mmapprData)
  # # mmapprData <- PeakRefinement(mmapprData)
  # mmapprData <- GenerateCandidates(mmapprData)
  # 
  # OutputMmapprData(mmapprData)
  print(proc.time() - startTime)
  return(mmapprData)
}



