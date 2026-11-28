test_that("LOESS helper returns a loess fit for ordinary data", {
  pos <- 1:50
  y <- sin(pos / 8)^2 + (pos %% 5) / 100
  fit <- MMAPPR2:::.getLoess(.4, pos, y)
  expect_s3_class(fit, "loess")
})

test_that("local span resolution handles boundaries, missing spans, and singleton grids", {
  f <- MMAPPR2:::.localResolution
  expect_equal(f(c(.1, .2, .4), .1), .1)
  expect_equal(f(c(.1, .2, .4), .2), .2)
  expect_equal(f(c(.1, .2, .4), .4), .2)
  expect_equal(f(c(.1, .2, .4), .3), Inf)
  expect_equal(f(.2, .2), Inf)
  expect_equal(f(c(NA, .1, .2), .1), .1)
})

test_that("AICc helper returns NA for failed, overfit, and zero-variance fits", {
  testthat::local_mocked_bindings(.getLoess = function(...) structure("fail", class = "try-error"), .package = "MMAPPR2")
  expect_true(is.na(MMAPPR2:::.aicc(.2, 1:5, 1:5)))

  testthat::local_mocked_bindings(.getLoess = function(...) list(n = 3, trace.hat = 2, residuals = c(1, 1, 1)), .package = "MMAPPR2")
  expect_true(is.na(MMAPPR2:::.aicc(.2, 1:3, 1:3)))

  testthat::local_mocked_bindings(.getLoess = function(...) list(n = 10, trace.hat = 2, residuals = rep(0, 10)), .package = "MMAPPR2")
  expect_true(is.na(MMAPPR2:::.aicc(.2, 1:10, 1:10)))
})

test_that("AICc optimizer refines around promising spans and returns only valid sorted spans", {
  testthat::local_mocked_bindings(
    .aicc = function(s, eucDist, pos) (s - .37)^2 + 1,
    .package = "MMAPPR2"
  )
  d <- data.frame(POS = 1:20, DISTANCE = rep(1, 20))
  out <- MMAPPR2:::.aiccOpt(d, spans = c(-.1, .1, .5, 1.1), resolution = .01, cutFactor = .5)
  expect_true(all(diff(out$spans) > 0))
  expect_true(all(out$spans > 0 & out$spans <= 1))
  expect_gt(nrow(out), 2)
  expect_lt(min(abs(out$spans - .37)), .1)
})

test_that("AICc optimizer errors when every candidate fit fails", {
  testthat::local_mocked_bindings(.aicc = function(...) NA_real_, .package = "MMAPPR2")
  d <- data.frame(POS = 1:20, DISTANCE = rep(1, 20))
  expect_error(MMAPPR2:::.aiccOpt(d, spans = c(.1, .2), resolution = .01, cutFactor = .5),
               "All candidate LOESS/AICc fits failed")
})

test_that("local minimum detection handles boundaries, ties, and nonfinite values", {
  f <- MMAPPR2:::.localMinIndices
  expect_identical(f(numeric()), integer())
  expect_identical(f(c(1, 2, 3)), 1L)
  expect_identical(f(c(2, 1, 1, 2)), c(2L, 3L))
  expect_identical(f(c(NA, 2, Inf, 1)), c(2L, 4L))
})

test_that("decimal precision rejects invalid resolutions and handles small values", {
  expect_identical(MMAPPR2:::.numDecimals(.00001), 5L)
  expect_error(MMAPPR2:::.numDecimals(0))
  expect_error(MMAPPR2:::.numDecimals(-.1))
  expect_error(MMAPPR2:::.numDecimals(Inf))
})

test_that("best AICc span errors without finite values and selects the unique optimum", {
  expect_error(MMAPPR2:::.chooseBestAiccSpan(data.frame(spans = c(.1, .2), aiccValues = c(NA, Inf)), .01),
               "No finite AICc")
  expect_equal(MMAPPR2:::.chooseBestAiccSpan(data.frame(spans = c(.1, .2, .3), aiccValues = c(2, 1, 3)), .01), .2)
})

test_that("per-chromosome LOESS fitting has success and isolated failure contracts", {
  pos <- 1:80
  d <- data.frame(POS = pos, DISTANCE = sin(pos / 9)^2 + (pos %% 7) / 100)
  input <- list(wtCounts = data.table::data.table(), mutCounts = data.table::data.table(), distanceDf = d)
  out <- MMAPPR2:::.loessFitForChr(input, .05, .5)
  expect_true(is.list(out))
  expect_s3_class(out$loess, "loess")
  expect_true(is.data.frame(out$aicc))
  expect_true(is.finite(out$bestSpan))
  expect_null(out$distanceDf)
  expect_true("loessTime" %in% names(out))

  expect_match(MMAPPR2:::.loessFitForChr("chr1: upstream failure", .05, .5), "upstream failure")
  few <- list(distanceDf = data.frame(POS = 1:4, DISTANCE = 1:4))
  expect_match(MMAPPR2:::.loessFitForChr(few, .05, .5), "Too few informative")
})

test_that("loessFit directly transforms successful chromosomes while isolating failed ones", {
  p <- .test_param(loessOptResolution = .05, loessOptCutFactor = .5)
  md <- mmapprData(p)
  pos <- 1:80
  md@snpDistance <- list(
    chr1 = list(wtCounts = data.table::data.table(), mutCounts = data.table::data.table(),
                distanceDf = data.frame(POS = pos, DISTANCE = sin(pos / 9)^2 + (pos %% 7) / 100)),
    chr2 = "chr2: no data"
  )
  out <- loessFit(md)
  expect_named(snpDistance(out), c("chr1", "chr2"))
  expect_s3_class(snpDistance(out)$chr1$loess, "loess")
  expect_true(is.data.frame(snpDistance(out)$chr1$aicc))
  expect_type(snpDistance(out)$chr2, "character")
})
