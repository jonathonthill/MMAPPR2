context("BAM file reading")

DebugSkip <- function() { 
    if (F) skip("Skipping file reading tests to save time") 
}


param <- MmapprParam(new("GmapGenome"), "./test_data/bam_files/zy14_wt_cut_filt.bam", 
                     "./test_data/bam_files/zy14_mut_cut_filt.bam",
                     vepFlags = ensemblVEP::VEPFlags())
mmapprData <- new("MmapprData", param = param)


test_that("correct ranges are being read", {
    DebugSkip()
    chrList <- .getFileReadChrList(mmapprData)
    expect_equal(length(chrList), 26)
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
    expect_equal_to_reference(result$distanceDf, "test_data/reference/read_bam_distance.RDS")
})


test_that("whole genome is read correctly", {
    DebugSkip()
    cat("before\n")
    mmapprData <- readInFiles(mmapprData, showDebug=T, silent = F)
    cat("after\n")
    expect_equal(length(mmapprData@distance), 26)
    
    classes <- lapply(mmapprData@distance, class)
    expect_equal(sum(classes == "character"), 25)
    expect_type(mmapprData@distance$chr5, "list")
    
    # saveRDS(mmapprData, "test_data/intermediate_MDs/post_file_read.RDS")
})
