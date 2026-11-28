test_that("decimal precision helper behaves on common resolutions", {
  expect_identical(MMAPPR2:::.numDecimals(1), 0L)
  expect_identical(MMAPPR2:::.numDecimals(.1), 1L)
  expect_identical(MMAPPR2:::.numDecimals(.001), 3L)
})

test_that("local minimum detection finds separated minima", {
  expect_identical(MMAPPR2:::.localMinIndices(c(3, 2, 3, 1, 2)), c(2L, 4L))
})

test_that("tied AICc optimum is one of the evaluated spans", {
  tab <- data.frame(spans = c(.1, .2, .3), aiccValues = c(1, 2, 1))
  best <- MMAPPR2:::.chooseBestAiccSpan(tab, .001)
  expect_true(best %in% tab$spans)
  expect_equal(best, .1)
})
