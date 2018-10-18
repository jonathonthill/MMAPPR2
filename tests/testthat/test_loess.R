context("Loess fit/AICc optimization")


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

test_that("chromosome is fit correctly", {
    set.seed(1)
    chr7list <- list(distanceDf=data.frame(
        pos=seq_len(1000),
        distance=stats::rnorm(1000, 0.01, sd = 0.0005)
    ), wtCounts=1, mutCounts=1)
    chr7list <- expect_warning(.loessFitForChr(
        chr7list,
        loessOptResolution = .001,
        loessOptCutFactor = 0.01
    ))
    expect_type(chr7list, "list")
    expect_true(all(
        c("wtCounts", "mutCounts", "loess", "aicc") %in%
            names(chr7list)
    ))
    expect_equal(length(chr7list$loess$x), 1000)
    expect_equal(chr7list$loess$pars$span, 0.919)
})


test_that("loessFit runs properly for whole mmapprData", {
    md <- readRDS('test_data/objects/post_calcdist_dummy_md.RDS')
    md <- loessFit(md)
    successes <-
        sapply(md@distance, function(seq)
            is(seq, 'list'))
    expect_equal(length(successes), 25)
    expect_equal(sum(successes), 0)
    
    expected <- readRDS('test_data/objects/post_loess_dummy_md.RDS')
    
    expect_identical(md@distance, expected@distance)
})
