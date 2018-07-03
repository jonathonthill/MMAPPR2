context('Output functions')

mmapprData <- readRDS('test_data/objects/post_candidates_chr5_md.RDS')

test_that('default output folder is prepared correctly', {
    outputFolder(mmapprData@param) <- 'DEFAULT'
    mockery::stub(.prepareOutputFolder, '.defaultOutputFolder', 'tmp')
    expect_
})