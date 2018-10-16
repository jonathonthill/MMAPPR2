context("Main and helper functions")

vepFlags <- readRDS('test_data/objects/vep_flags.RDS')
param <- MmapprParam(new("GmapGenome"),
                     'test_data/bam_files/zy14_dummy.bam',
                     'test_data/bam_files/zy14_dummy.bam',
                     species='danio_rerio',
                     vepFlags=vepFlags)

test_that(".runFunctionInParallel works on list with single item", {
  input <- list(a='test')
  output <- .runFunctionInParallel(input, function(x) paste(x, x))
  expect_identical(output, list(a='test test'))
})

test_that(".runFunctionInParallel works on list with multiple items", {
    input <- list(a='test', b='test')
    expected <- list(a='test test', b='test test')
    output <- .runFunctionInParallel(input, function(x) paste(x, x))
    expect_identical(output, expected)
})

test_that(".runFunctionInParallel work with package functions", {
    testFunction <- function(x) {
        return(rep(x, 2))
    }
    x <- list(a=1, b=1, c=1)
    y <- .runFunctionInParallel(x, testFunction)
    expResult <- list(a=c(1, 1), b=c(1, 1), c=c(1, 1))
    expect_identical(y, expResult)
})

test_that(".runFunctionInParellel works with extra inputs", {
    testFun <- function(a, b=0, c=0) {
        return(c(a, b, c))
    }
    x <- list(1)
    expect_identical(.runFunctionInParallel(x, testFun),
                     list(c(1, 0, 0)))
    expect_identical(.runFunctionInParallel(x, testFun, b=1),
                     list(c(1, 1, 0)))
    expect_identical(.runFunctionInParallel(x, testFun, b=1, c=1),
                     list(c(1, 1, 1)))
})

test_that("MmapprParam has default VEPFlags when one isn't provided", {
    # put fake vep in path since VEPFlags requires it
    dir.create('/tmp/ensembl-vep')
    system('touch /tmp/ensembl-vep/vep; chmod 777 /tmp/ensembl-vep/vep')
    originalPath <- Sys.getenv('PATH')
    Sys.setenv('PATH'=paste0(originalPath, ':', '/tmp/ensembl-vep'))
    
    mp <- MmapprParam(new("GmapGenome"),
                      'test_data/bam_files/zy14_dummy.bam',
                      'test_data/bam_files/zy14_dummy.bam',
                      species='danio_rerio')
    expect_false(is.null(mp@vepFlags))
    expect_true(ensemblVEP::flags(vepFlags(mp))$database == FALSE)
    expect_true(ensemblVEP::flags(vepFlags(mp))$species == 'danio_rerio')
    
    # cleanup
    Sys.setenv('PATH'=originalPath)
    unlink('/tmp/ensembl-vep', recursive=TRUE)
})


test_that("MmapprParam takes character, BamFile, or BamFileList of real files", {
    fn_wt <- 'test_data/bam_files/zy14_dummy.bam'
    fn_mut <- 'test_data/bam_files/zy14_dummy.bam'

    # non-existing filename shouldn't work
    expect_error(MmapprParam(new("GmapGenome"), "wt", 'mut',
                             'danio_rerio', vepFlags))

    # existing filename, just character, should work
    param <- MmapprParam(new("GmapGenome"), fn_wt, fn_mut,
                         'danio_rerio', vepFlags)
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")
    
    # character with length>1
    param <- MmapprParam(new("GmapGenome"), c(fn_wt, fn_wt), fn_mut,
                         'danio_rerio', vepFlags)
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    bf <- Rsamtools::BamFile(fn_mut)
    param <- MmapprParam(new("GmapGenome"), bf, bf, 'danio_rerio', vepFlags)
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    bfl <- Rsamtools::BamFileList(fn_mut)
    param <- MmapprParam(new("GmapGenome"), bfl, bfl, 'danio_rerio', vepFlags)
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

})

test_that('file setters for param can change format automatically', {
    fn <- 'test_data/bam_files/zy14_dummy.bam'
    param <- MmapprParam(new("GmapGenome"), fn, fn, 'danio_rerio', vepFlags)
    
    wtFiles(param) <- c(fn, fn)
    mutFiles(param) <- c(fn, fn)
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")

    bf <- Rsamtools::BamFile(fn)
    wtFiles(param) <- c(bf, bf)
    mutFiles(param) <- c(bf, bf)
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")
})

test_that('log writes to default folder', {
    skip_if_not_installed('mockery')
    newParam <- param
    unlink('/tmp/m2', recursive=TRUE)
    dir.create('/tmp/m2')
    newParam@outputFolder <- '/tmp/m2'
    mockery::stub(mmappr, '.prepareOutputFolder', new('MmapprData', param=newParam))
    mockery::stub(mmappr, 'tryCatch', function(...) stop('expected failure'))
    expect_error(mmappr(param), 'expected failure')
    expect_true(file.exists('/tmp/m2/mmappr2.log'))
    unlink('/tmp/m2', recursive=TRUE)
})
