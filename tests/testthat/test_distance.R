context("Calculate distance data from BAM files")
Sys.unsetenv("R_TESTS")


vepFlags <- readRDS('test_data/objects/vep_flags.RDS')
param <-
    MmapprParam(
        refFasta='test_data/dummy.fasta',
        wtFiles='test_data/bam_files/zy14_dummy.bam',
        mutFiles='test_data/bam_files/zy14_dummy.bam',
        species='danio_rerio',
        vepFlags=vepFlags,
        refGenome=new("GmapGenome")
    )
mmapprData <- new("MmapprData", param = param)

# with dummy files
test_that("whole genome is read correctly", {
    #skip_if_not_travis_or_bioc()
    skip_if_not_in_path('samtools')
    mmapprData <- calculateDistance(mmapprData)
    expect_identical(
        mmapprData@distance,
        readRDS('test_data/objects/post_calcdist_dummy_md.RDS')@distance
    )
})


test_that("files get indexed automatically if needed", {
    skip_if_not_installed('mockery')
    mockery::stub(calculateDistance, '.getFileReadChrList',
                  function(...) stop('fail'))
    mockery::stub(calculateDistance, 'Rsamtools::indexBam',
                  function(...) write.table(matrix(), '/tmp/test.bam.bai'))

    file.copy('test_data/bam_files/zy14_dummy.bam', '/tmp/test.bam')

    expect_false(file.exists('/tmp/test.bam.bai'))
    wtFiles(mmapprData@param) <- '/tmp/test.bam'
    mutFiles(mmapprData@param) <- '/tmp/test.bam'
    expect_error(calculateDistance(mmapprData), 'fail')
    expect_true(file.exists('/tmp/test.bam.bai'))

    unlink('/tmp/test.bam*')
})


test_that('files in BamFileList also get indexed automatically if needed', {
    skip_if_not_installed('mockery')
    mockery::stub(calculateDistance, '.getFileReadChrList', function(...)
        stop('fail'))
    mockery::stub(calculateDistance, 'Rsamtools::indexBam',
                  function(...) write.table(matrix(), '/tmp/test.bam.bai'))

    file.copy('test_data/bam_files/zy14_dummy.bam', '/tmp/test.bam')
    file.copy('test_data/bam_files/zy14_dummy.bam', '/tmp/test2.bam')

    expect_false(file.exists('/tmp/test.bam.bai'))
    wtFiles(mmapprData@param) <- c('/tmp/test.bam', '/tmp/test2.bam')
    mutFiles(mmapprData@param) <- '/tmp/test.bam'
    expect_error(calculateDistance(mmapprData), 'fail')
    expect_true(file.exists('/tmp/test.bam.bai'))
    expect_true(file.exists('/tmp/test2.bam.bai'))

    unlink('/tmp/test*.bam*')

})


test_that("correct ranges are being read", {
    chrList <- .getFileReadChrList(mmapprData)
    expect_equal(length(chrList), 25)
    expect_type(chrList, "list")
    expect_s4_class(chrList$chr5, "GRanges")
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
    skip_if_not_in_path('samtools')

    inputRange <-
        GenomicRanges::GRanges('7',
                               IRanges::IRanges(start=1, width=999))

    infoVec <- rep(c(1, 1, 1, 1, 4)*10, 10) # coverage=40
    mockPileup <- mockery::mock(
        data.frame(pos=rep(1:10, each=5), nucleotide=c('A', 'C', 'G', 'T', 'cvg'),
                   file1=infoVec),  # WT
        data.frame(pos=rep(2:11, each=5), nucleotide=c('A', 'C', 'G', 'T', 'cvg'),
                   file1=infoVec),  # MUT
    )
    mockery::stub(.calcDistForChr, '.samtoolsPileup', mockPileup)
    result <- .calcDistForChr(inputRange, param=param)
    mockery::expect_called(mockPileup, 2)
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
    skip_if_not_in_path('samtools')

    inputRange <-
        GenomicRanges::GRanges('7',
                               IRanges::IRanges(start=1, width=999))

    infoVec <- rep(c(1, 1, 1, 1, 4)*10, 10) # coverage=40
    mockPileup <- mockery::mock(
        data.frame(pos=rep(1:10, each=5), nucleotide=c('A', 'C', 'G', 'T', 'cvg'),
                   file1=infoVec),  # WT
        data.frame(pos=rep(2:11, each=5), nucleotide=c('A', 'C', 'G', 'T', 'cvg'),
                   file1=infoVec),  # MUT
    )
    mockery::stub(.calcDistForChr, '.samtoolsPileup', mockPileup)
    result <- .calcDistForChr(inputRange, param=param)
    mockery::expect_called(mockPileup, 2)
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
