context('Output functions')

mmapprData <- readRDS('test_data/objects/post_peakref_dummy_md.RDS')


test_that('default output folder is prepared correctly', {
    skip_if_not_installed('mockery')
    unlink('tmp', recursive=TRUE)
    outputFolder(mmapprData@param) <- 'DEFAULT'
    mockery::stub(.prepareOutputFolder, '.defaultOutputFolder', 'tmp')
    mmapprData <- .prepareOutputFolder(mmapprData)
    expect_true(dir.exists('tmp'))
    unlink('tmp', recursive=TRUE)
})

test_that('user can overwrite previous output folder when prompted', {
    skip_if_not_installed('mockery')
    unlink('tmp', recursive=TRUE)
    dir.create('tmp')
    write.table(matrix(), file='tmp/tmp.txt')
    expect_true(file.exists('tmp/tmp.txt'))
    outputFolder(mmapprData@param) <- 'tmp'
    mockery::stub(.prepareOutputFolder, 'readline', 'y')
    mmapprData <- .prepareOutputFolder(mmapprData)
    expect_equal(mmapprData@param@outputFolder, 'tmp')
    expect_true(dir.exists('tmp'))
    expect_false(file.exists('tmp/tmp.txt'))
    unlink('tmp', recursive=TRUE)
})

test_that('user can change output folder name to avoid overwriting', {
    skip_if_not_installed('mockery')
    unlink(c('tmp', 'tmp2'), recursive=TRUE)
    dir.create('tmp')
    outputFolder(mmapprData@param) <- 'tmp'
    mock_readline <- mockery::mock('n', 'tmp2')
    mockery::stub(.prepareOutputFolder, 'readline', mock_readline)
    mmapprData <- .prepareOutputFolder(mmapprData)
    expect_equal(mmapprData@param@outputFolder, 'tmp2')
    expect_true(dir.exists('tmp'))
    expect_true(dir.exists('tmp2'))
    unlink(c('tmp', 'tmp2'), recursive=TRUE)
})

test_that('user can choose default output folder to avoid overwriting', {
    skip_if_not_installed('mockery')
    unlink(c('tmp', 'default'), recursive=TRUE)
    dir.create('tmp')
    outputFolder(mmapprData@param) <- 'tmp'
    mock_readline <- mockery::mock('n', NA)
    mockery::stub(.prepareOutputFolder, 'readline', mock_readline)
    mockery::stub(.prepareOutputFolder, '.defaultOutputFolder', 'default')
    mmapprData <- .prepareOutputFolder(mmapprData)
    expect_equal(mmapprData@param@outputFolder, 'default')
    expect_true(dir.exists('tmp'))
    expect_true(dir.exists('default'))
    unlink(c('tmp', 'default'), recursive=TRUE)
})

test_that('multiple candidate tables get written', {
    unlink('tmp', recursive=TRUE)
    skip_if_not_installed('mockery')
    gr <- GenomicRanges::GRanges()
    candList <- list(chr1=gr, chr2=gr)
    candDF <- data.frame("Position" = 1,
                         "Feature" = 'del',
                         "Symbol" = 'GENE',
                         "Allele" = 'C',
                         "Consequence" = 'missense',
                         "AminoAcid" = 'Y',
                         "Impact" = 'HIGH',
                         "DensityScore" = 1234
    )
    mockery::stub(.writeCandidateTables, 'data.frame', candDF)
    .writeCandidateTables(candList, tempdir())
    expect_true(file.exists(file.path(tempdir(), 'chr1.tsv')))
    expect_true(file.exists(file.path(tempdir(), 'chr2.tsv')))
})
