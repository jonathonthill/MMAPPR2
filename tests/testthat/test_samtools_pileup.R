context('Samtools Pileup')


test_that("Simple pileup works", {
    skip_if_not_installed('mockery')
    output <- paste(
        '1    1    A    5    ....C    AAAAA',
        '1    2    A    5    ...C.    AAAAA',
        '1    3    A    5    ..C..    AAAAA',
        sep='\n'
    )
    mockOutput <- mockery::stub(.samtoolsPileup, 'system2', output)

    files <- Rsamtools::BamFileList(c(Rsamtools::BamFile("f")))
    pileup <- .samtoolsPileup(files, new('MmapprParam', minDepth=0), new("GRanges"))
    expect_true(all(names(pileup) %in% c('nucleotide', 'pos', 'f')))
    expect_equal(nrow(pileup), 15)
    expect_equal(as.character(pileup$nucleotide[1:5]), c('A', 'C', 'G', 'T', 'cvg'))
})


test_that("Pileup filters by depth", {
    skip_if_not_installed('mockery')
    output <- paste(
        '1    1    A    5    ....C    AAAAA',
        '1    2    A    5    ...C.    AAAAA',
        '1    3    A    6    ..C...    AAAAAA',
        sep='\n'
    )
    mockOutput <- mockery::stub(.samtoolsPileup, 'system2', output)

    files <- Rsamtools::BamFileList(c(Rsamtools::BamFile("f")))
    pileup <- .samtoolsPileup(files, new('MmapprParam', minDepth=6), new("GRanges"))
    expect_true(all(names(pileup) %in% c('nucleotide', 'pos', 'f')))
    expect_equal(nrow(pileup), 5)
    expect_equal(as.character(pileup$nucleotide[1:5]), c('A', 'C', 'G', 'T', 'cvg'))
})


test_that("Pileup for multiple files is read correctly", {
    skip_if_not_installed('mockery')

    output <- paste(
        '1    1    A    5    ....C    AAAAA',
        '1    2    A    5    ...C.    AAAAA',
        '1    3    A    5    ..C..    AAAAA',
        sep='\n'
    )
    mockOutput <- mockery::mock(output, output)
    mockery::stub(.samtoolsPileup, 'system2', mockOutput)

    files <- Rsamtools::BamFileList(c(Rsamtools::BamFile("f1"), Rsamtools::BamFile("f2")))
    pileup <- .samtoolsPileup(files, new('MmapprParam', minDepth=5), new("GRanges"))

    expect_true(all(names(pileup) %in% c('nucleotide', 'pos', 'f1', 'f2')))
    expect_equal(nrow(pileup), 15)
    expect_equal(as.character(pileup$nucleotide[1:5]), c('A', 'C', 'G', 'T', 'cvg'))

    mockery::expect_called(mockOutput, 2)
})


test_that("Depth filter leaves NAs with outer join", {
    skip_if_not_installed('mockery')
    
    output1 <- paste(
        '1    1    A    4    ....     AAAA',
        '1    2    A    5    ...C.    AAAAA',
        '1    3    A    4    ..C.     AAAA',
        sep='\n'
    )
    output2 <- paste(
        '1    1    A    5    ....C    AAAAA',
        '1    2    A    4    ...C     AAAA',
        '1    3    A    4    ..C.     AAAA',
        sep='\n'
    )
    mockOutput <- mockery::mock(output1, output2)
    mockery::stub(.samtoolsPileup, 'system2', mockOutput)
    
    files <- Rsamtools::BamFileList(c(Rsamtools::BamFile("f1"), Rsamtools::BamFile("f2")))
    pileup <- .samtoolsPileup(files, new('MmapprParam', minDepth=5), new("GRanges"))
    
    expect_true(all(names(pileup) %in% c('nucleotide', 'pos', 'f1', 'f2')))
    expect_equal(nrow(pileup), 10)
    expect_equal(as.character(pileup$nucleotide[1:5]), c('A', 'C', 'G', 'T', 'cvg'))
    expect_equal(sum(pileup$pos == 1), 5)
    expect_equal(sum(pileup$pos == 2), 5)
    expect_equal(as.numeric(pileup$f1[pileup$pos == 1]), rep(as.numeric(NA), times=5))
    expect_equal(as.numeric(pileup$f2[pileup$pos == 2]), rep(as.numeric(NA), times=5))
    
    mockery::expect_called(mockOutput, 2)
    
    mockery::stub(.calcDistForChr, '.samtoolsPileup', mockery::mock(pileup, pileup))
})


test_that("Bad data throws error", {
    skip_if_not_installed('mockery')
    output <- paste(
        '1    1    A    5   $$%^#Q ....C    AAAAA',
        '1    2  3sdfsd  5    ...C.    AAAAA',
        'Z    3    GARBAGE..C..    AAAAA',
        sep='\n'
    )
    mockOutput <- mockery::stub(.samtoolsPileup, 'system2', output)
    
    expect_error(samtoolsPileup(
        Rsamtools::BamFile("f"), new('MmapprParam', minDepth=0), new("GRanges"))
    )
})