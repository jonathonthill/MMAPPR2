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

source('~/MMAPPR2/read_bam.r')

start <- proc.time()


#automatic core number generation (only on Linux)
#mem_req is memory required per base per replicate (per two files) in MB
#I calculated it myself in a test using only the memory taken up by 
#the apply_pileup and dividing it by reps and sequence length
#in other words, it doesn't take into account all the memory used by 
#each core as the LoessFit (hence the getChrDf) function is run
core_calc <- function(mem_req = 33.3){
  tryCatch({
    #this gets free RAM (on Linux)
    free_ram <- system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE) #gets mem info in MB
    #get numeric, in bytes (not MB)
    free_ram <- 1000 * as.numeric(free_ram)
    #get longest chromosome, calculate needed memory, and divide into free RAM
    num_cores <- ceiling(free_ram / (max(seqlengths(myrange)) * mem_req * length(wt_list)))
    message("num_cores = ", num_cores)
    return(num_cores)
    },
  warning = function(w) {
    message('Warning: arbitrary number of clusters used (are you not on Linux?)')
    return(ceiling(detectCores() / 3))
  }
  )
}

RunFunctionInParallel <- function(inputList, functionToRun, packages = c()) {
  #cluster generation
  cl <- makeCluster(core_calc(), type = "SOCK", outfile = "")
  
  # register the cluster
  registerDoParallel(cl)
  
  #perform calculation
  message("Beginning parallel calculation")
  resultList <- foreach(i = inputList,
                        .packages = packages) %dopar% functionToRun(i)
  names(resultList) <- names(inputList)
  
  stopCluster(cl)
  # insert serial backend, otherwise error in repetetive tasks
  registerDoSEQ()
  
  # clean up a bit.
  invisible(gc)
  
  return(resultList)
}





