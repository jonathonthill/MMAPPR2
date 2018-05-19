context("Main and helper functions")

vepFlags <- ensemblVEP::VEPFlags(flags=list(format="vcf"))
param <- MmapprParam(new("GmapGenome"), "./test_data/bam_files/zy14_wt_cut_filt.bam",
                     "test_data/bam_files/zy14_mut_cut_filt.bam",
                     vepFlags = vepFlags)

test_that(".runFunctionInParallel works on single core", {
  input <- list(a='test')
  output <- .runFunctionInParallel(input, function(x) paste(x, x), 1)
  expect_identical(output, list(a='test test'))
})

test_that(".runFunctionInParllel works on multiple cores", {
    input <- list(a='test', b='test')
    expected <- list(a='test test', b='test test')
    output <- .runFunctionInParallel(input, function(x) paste(x, x), 2)
    expect_identical(output, expected)
})

test_that(".runFunctionInParallel work with packages", {
    testFunction <- function(x) {
        return(.repList(x, 2))
    }
    x <- list(a=1, b=1, c=1)
    y <- .runFunctionInParallel(x, testFunction, silent=TRUE,
                                         numCores=MMAPPR2:::.coreCalc(),
                                         packages=c("MMAPPR2"))
    expResult <- list(a=list(1, 1), b=list(1, 1), c=list(1, 1))
    expect_identical(y, expResult)
})

test_that(".runFunctionInParellel works with extra inputs", {
    testFun <- function(a, b=0, c=0) {
        return(c(a, b, c))
    }
    x <- list(1)
    expect_identical(.runFunctionInParallel(x, testFun, 1, silent=TRUE),
                     list(c(1, 0, 0)))
    expect_identical(.runFunctionInParallel(x, testFun, 1, silent=TRUE,
                                                     secondInput=1),
                     list(c(1, 1, 0)))
    expect_identical(.runFunctionInParallel(x, testFun, 1, silent=TRUE,
                                                     secondInput=1, thirdInput=1),
                     list(c(1, 1, 1)))
})

test_that(".repList function works", {
    expect_equal(length(.repList(param, 5)), 5)
    expect_equal(length(.repList(param, 0)), 0)
    expect_equal(length(.repList(param, 1)), 1)
    expect_equal(length(.repList(param, 2)), 2)
    expect_error(length(.repList(param, -1)))
    x <- .repList(param, 2)
    expect_identical(x[[1]], x[[2]])
})

test_that("MmapprParam takes character, BamFile, or BamFileList of real files", {
    fn_wt <- "test_data/bam_files/zy14_wt_cut_filt.bam"
    fn_mut <- "test_data/bam_files/zy14_mut_cut_filt.bam"

    # non-existing filename shouldn't work
    expect_error(MmapprParam(new("GmapGenome"), "wt", 'mut', vepFlags))

    # existing filename, just character, should work
    param1 <- MmapprParam(new("GmapGenome"), fn_wt, fn_mut, vepFlags)
    expect_s4_class(param1, "MmapprParam")
    expect_s4_class(param1@wtFiles, "BamFileList")
    expect_s4_class(param1@mutFiles, "BamFileList")

    bf <- Rsamtools::BamFile("test_data/bam_files/zy14_mut_cut_filt.bam")
    param2 <- MmapprParam(new("GmapGenome"), bf, bf, vepFlags=vepFlags)
    expect_s4_class(param2, "MmapprParam")
    expect_s4_class(param2@wtFiles, "BamFileList")
    expect_s4_class(param2@mutFiles, "BamFileList")

    bfl <- Rsamtools::BamFileList("test_data/bam_files/zy14_mut_cut_filt.bam")
    param3 <- MmapprParam(new("GmapGenome"), bfl, bfl, vepFlags=vepFlags)
    expect_s4_class(param3, "MmapprParam")
    expect_s4_class(param3@wtFiles, "BamFileList")
    expect_s4_class(param3@mutFiles, "BamFileList")

})