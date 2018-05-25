context("BAM file reading")

DebugSkip <- function() { 
    if (F) skip("Skipping file reading tests to save time") 
}


numCores <- parallel::detectCores() / 2
param <- MmapprParam(new("GmapGenome"), "./test_data/bam_files/zy14_wt_cut_filt.bam", 
                     "./test_data/bam_files/zy14_mut_cut_filt.bam", numCores=numCores,
                     vepFlags = ensemblVEP::VEPFlags(flags=list(format='vcf')))
mmapprData <- new("MmapprData", param = param)


test_that("correct ranges are being read", {
    DebugSkip()
    chrList <- .getFileReadChrList(mmapprData)
    expect_equal(length(chrList), 25)
    expect_type(chrList, "list")
    expect_type(chrList[[1]], "list")
    
    expect_named(chrList$chr5, c("range", "param"), ignore.order=TRUE)
    expect_s4_class(chrList$chr5$range, "GRanges")
    expect_s4_class(chrList$chr5$param, "MmapprParam")
})


test_that("single chromosome is read correctly", {
    DebugSkip()
    inputList <- .getFileReadChrList(mmapprData)[['chr5']]
    expect_type(inputList, "list")
    
    result <- .readFilesForChr(inputList)
    expect_named(result, c("wtCounts", "mutCounts", "distanceDf", "seqname"),
                 ignore.order=TRUE)
    expect_gt(nrow(result$wtCounts), 0)
    expect_gt(nrow(result$mutCounts), 0)
    expect_gt(nrow(result$distanceDf), 0)
    expect_named(result$distanceDf, c("pos", "distance"))
    expect_known_value(result$distanceDf, "test_data/objects/chr5_distance.RDS",
                              update=FALSE)
})


test_that("whole genome is read correctly", {
    DebugSkip()
    # with dummy files
    wtFiles(mmapprData@param) <- 'test_data/bam_files/zy14_dummy.bam'
    mutFiles(mmapprData@param) <- 'test_data/bam_files/zy14_dummy.bam'
    cat("before\n")
    mmapprData <- readInFiles(mmapprData)
    cat("after\n")
    expect_known_value(mmapprData@distance, "test_data/objects/post_file_read_dummy_distance.RDS",
                       update=FALSE)
})
