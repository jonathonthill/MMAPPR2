context("BAM file reading")
Sys.unsetenv("R_TESTS")


vepFlags <- readRDS('test_data/objects/vep_flags.RDS')
param <-
    MmapprParam(
        new("GmapGenome"),
        "./test_data/bam_files/zy14_wt_cut_filt.bam",
        "./test_data/bam_files/zy14_mut_cut_filt.bam",
        species='danio_rerio',
        vepFlags=vepFlags
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

test_that('chrM and MT are dropped and sequences are reordered', {
    skip_if_not_installed('mockery')
    ucscNames <- c('chr1', 'chr2', 'chrM')
    ensemblNames <- c('1', '2', 'MT')
    mockGR <- mockery::mock(GenomicRanges::GRanges(ucscNames,
                                                   IRanges::IRanges(0, 1)
                            ),
                            GenomicRanges::GRanges(ensemblNames,
                                                   IRanges::IRanges(0, 1)
                            )
    )
    mockery::stub(.getFileReadChrList, 'as', mockGR)
    chrList <- .getFileReadChrList(mmapprData)
    expect_equal(names(chrList), c('chr1', 'chr2'))
    
    chrList <- .getFileReadChrList(mmapprData)
    expect_equal(names(chrList), c('1', '2'))
    
    mockery::expect_called(mockGR, 2)
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


test_that('.avgFiles works as expected with sum fileAggregation', {
    inDf <- data.frame(pos=rep(c(1, 2, 3), each=5),
                       nucleotide=rep(c('A', 'C', 'G', 'T', 'cvg'), 3),
                       file1=rep(c(4, 0, 0, 0, 4), 3),
                       file2=c(0, 6, 0, 0, 6,
                               0, 0, 23, 23, 46,
                               48, 0, 0, 48, 96)
            )
    result <- .avgFiles(inDf, fileAggregation='sum')
    expected <- data.frame(pos=rep(c(1, 2, 3), each=5),
                           nucleotide=rep(c('A', 'C', 'G', 'T', 'cvg'), 3),
                           avgCount=c(.4, .6, 0, 0, 10,
                                      .08, 0, .46, .46, 50,
                                      .52, 0, 0, .48, 100))
    expect_identical(result, expected)
})


test_that('.avgFiles works as expected with mean fileAggregation', {
    inDf <- data.frame(pos=rep(c(1, 2, 3), each=5),
                       nucleotide=rep(c('A', 'C', 'G', 'T', 'cvg'), 3),
                       file1=rep(c(4, 0, 0, 0, 4), 3),
                       file2=c(0, 6, 0, 0, 6,
                               0, 0, 23, 23, 46,
                               48, 0, 0, 48, 96)
            )
    result <- .avgFiles(inDf, fileAggregation='mean')
    expected <- data.frame(pos=rep(c(1, 2, 3), each=5),
                           nucleotide=rep(c('A', 'C', 'G', 'T', 'cvg'), 3),
                           avgCount=c(.5, .5, 0, 0, 5,
                                      .5, 0, .25, .25, 25,
                                      .75, 0, 0, .25, 50))
    expect_identical(result, expected)
})


test_that("single chromosome is read correctly with replicates", {
    skip_if_not_installed('mockery')

    inputRange <-
        GenomicRanges::GRanges('chr5',
                               IRanges::IRanges(start=1, width=75682077))
    expect_s4_class(inputRange, "GRanges")

    infoVec <- rep(c(1, 1, 1, 1, 4)*10, 10) # coverage=40
    mockAP <- mockery::mock(
        list(list(pos=1:10, # WT
             info=matrix(rep(infoVec, 2), ncol=2))),
        list(list(pos=2:11, # MUT
             info=matrix(rep(infoVec, 3), ncol=3)))
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
