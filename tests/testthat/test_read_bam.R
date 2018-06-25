context("BAM file reading")
Sys.unsetenv("R_TESTS")


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
    skip_if_not_travis_or_bioc()
    wtFiles(mmapprData@param) <- 'test_data/bam_files/zy14_dummy.bam'
    mutFiles(mmapprData@param) <- 'test_data/bam_files/zy14_dummy.bam'
    mmapprData <- readInFiles(mmapprData)
    expect_known_value(
        mmapprData@distance,
        "test_data/objects/post_file_read_dummy_distance.RDS",
        update = FALSE
    )
})


test_that("correct ranges are being read", {
    skip_if_not_travis_or_bioc()
    chrList <- .getFileReadChrList(mmapprData)
    expect_equal(length(chrList), 25)
    expect_type(chrList, "list")
    expect_s4_class(chrList$chr5, "GRanges")

    rm(chrList)
})


test_that("single chromosome is read correctly", {
    skip_if_not_installed('mockery')

    inputRange <-
        GenomicRanges::GRanges('chr5',
                               IRanges::IRanges(start=1, width=75682077))
    expect_s4_class(inputRange, "GRanges")

    infoVec <- rep(c(1, 1, 1, 1, 4)*10, 10) # coverage=40
    mockAP <- mockery::mock(
        list(list(pos=1:10, # WT
             info=matrix(infoVec, ncol=1))),
        list(list(pos=2:11, # MUT
             info=matrix(infoVec, ncol=1)))
    )
    mockery::stub(.readFilesForChr, 'Rsamtools::applyPileups', mockAP)
    result <- .readFilesForChr(inputRange, param=param)
    mockery::expect_called(mockAP, 2)
    expect_true(all(
        c("wtCounts", "mutCounts", "distanceDf", "seqname") %in%
            names(result)
    ))
    expect_type(result, 'list')
    expect_gt(nrow(result$wtCounts), 0)
    expect_gt(nrow(result$mutCounts), 0)
    expect_gt(nrow(result$distanceDf), 0)
    expect_named(result$distanceDf, c("pos", "distance"))
    expect_is(result$distanceDf$pos, 'integer')
    expect_is(result$distanceDf$distance, 'numeric')
})
