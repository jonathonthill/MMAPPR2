context("Loess fit/AICc optimization")


md <- readRDS("test_data/intermediate_MDs/post_file_read.RDS")
md@param@loessOptResolution <- 0.001
md@param@loessOptCutFactor <- 0.01
SkipDebug <- function() if(T) skip("skipping to save time")

test_that(".minTwo works right", {
    x <- c(1, 2, 3)
    expect_equal(.minTwo(x), c(1, 2))
    x <- c(1, 3, 2)
    expect_equal(.minTwo(x), c(1, 2))
    x <- c(3, 2, 1)
    expect_equal(.minTwo(x), c(1, 2))
    x <- c(1.00000001, 1.000001, 1.0001)
    expect_identical(.minTwo(x), c(1.00000001, 1.000001))
    x <- c(-1, -2, -3)
    expect_identical(.minTwo(x), c(-3, -2))
})

test_that(".numDecimals works right", {
    expect_equal(.numDecimals(.005), 3)
    expect_equal(.numDecimals(.05), 2)
    expect_equal(.numDecimals(.5), 1)
    expect_equal(.numDecimals(5), 0)
    expect_equal(.numDecimals(50), 0)
})

test_that(".localMin works right", {
    expect_equal(.localMin(c(1, 3, 1, 3)), c(1, 1))
    expect_equal(.localMin(c(1, 3, 2, 3)), c(1, 2))
    expect_equal(.localMin(c(1, 3, 2.999999, 3)), c(1, 2.999999))
    expect_equal(.localMin(c(1, 3, 2, 3, 1)), c(1, 2, 1))
    expect_equal(.localMin(c(-1, -3, -2, -3, -1)), c(-3, -3))
})

test_that(".localResolution works right", {
    spans <- c(.31, .32, .34, .36, .36, .37, .40)
    expect_equal(.localResolution(spans, .31), .01)
    expect_equal(.localResolution(spans, .32), .02)
    expect_error(.localResolution(spans, .33))
    expect_equal(.localResolution(spans, .34), .02)
    expect_equal(.localResolution(spans, .37), .03)
    expect_equal(.localResolution(spans, .40), .03)

    spans <- c(-.01, 0, .01)
    expect_error(.localResolution(spans, .01))

    spans <- c(.32, .31, .34, .40, .36, .36, .37)
    expect_equal(.localResolution(spans, .31), .01)
    expect_equal(.localResolution(spans, .32), .02)
    expect_error(.localResolution(spans, .33))

})

test_that("normal chromosome performs fit correctly", {
    SkipDebug()
    chr5List <- .loessFitForChr(md@distance$chr5,
                                loessOptResolution = md@param@loessOptResolution,
                                loessOptCutFactor = md@param@loessOptCutFactor)
    expect_type(chr5List, "list")
    expect_true(all(c("wtCounts", "mutCounts", "loess", "aicc") %in% names(chr5List)))
    expect_gt(length(chr5List$loess$x), 5000)
    expect_equal(chr5List$loess$pars$span, 0.03)
})

test_that("Empty chromosome is skipped", {
    chr4error <- .loessFitForChr(md@distance$chr4,
                                 loessOptResolution = md@param@loessOptResolution,
                                 loessOptCutFactor = md@param@loessOptCutFactor)
    expect_type(chr4error, "character")
})

test_that("loessFit runs properly for whole mmapprData", {
    SkipDebug()
    md2 <- loessFit(md, silent=T)
    successes <- sapply(md2@distance, function(seq) class(seq) == "list")
    expect_equal(length(successes), 26)
    expect_equal(sum(successes), 1)


    expect_equal_to_reference(md2, "test_data/intermediate_MDs/post_loess.RDS")
})

