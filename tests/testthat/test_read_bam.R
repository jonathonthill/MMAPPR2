context("BAM file reading")
Sys.unsetenv("R_TESTS")

DebugSkip <- function() {
    if (F)
        skip("Skipping file reading tests to save time")
}


vepFlags <- readRDS('test_data/objects/vep_flags.RDS')
numCores <- ceiling(parallel::detectCores() / 2)
param <-
    MmapprParam(
        new("GmapGenome"),
        "./test_data/bam_files/zy14_wt_cut_filt.bam",
        "./test_data/bam_files/zy14_mut_cut_filt.bam",
        numCores = numCores,
        vepFlags = vepFlags
    )
mmapprData <- new("MmapprData", param = param)

# with dummy files
test_that("whole genome is read correctly", {
    wtFiles(mmapprData@param) <- 'test_data/bam_files/zy14_dummy.bam'
    mutFiles(mmapprData@param) <- 'test_data/bam_files/zy14_dummy.bam'
    DebugSkip()
    mmapprData <- readInFiles(mmapprData)
    expect_known_value(
        mmapprData@distance,
        "test_data/objects/post_file_read_dummy_distance.RDS",
        update = FALSE
    )
})


test_that("correct ranges are being read", {
    DebugSkip()
    chrList <- .getFileReadChrList(mmapprData)
    expect_equal(length(chrList), 25)
    expect_type(chrList, "list")
    expect_s4_class(chrList$chr5, "GRanges")

    rm(chrList)
})


test_that("single chromosome is read correctly", {
    DebugSkip()
    gc()
    inputRange <- .getFileReadChrList(mmapprData)[['chr5']]
    expect_s4_class(inputRange, "GRanges")

    result <- .readFilesForChr(inputRange, param=param)
    expect_true(all(
        c("wtCounts", "mutCounts", "distanceDf", "seqname") %in%
            names(result)
    ))
    expect_type(result, 'list')
    expect_gt(nrow(result$wtCounts), 0)
    expect_gt(nrow(result$mutCounts), 0)
    expect_gt(nrow(result$distanceDf), 0)
    expect_named(result$distanceDf, c("pos", "distance"))
    expect_known_value(result$distanceDf,
                       "test_data/objects/chr5_distance.RDS",
                       update = FALSE)
    rm(result)
})


