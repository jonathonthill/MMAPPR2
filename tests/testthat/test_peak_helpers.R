test_that("density point masses are normalized on nonuniform grids", {
  x <- c(0, 1, 3, 6)
  y <- c(1, 2, 2, 1)
  mass <- MMAPPR2:::.densityPointMass(x, y)
  expect_length(mass, length(x))
  expect_true(all(mass >= 0))
  expect_equal(sum(mass), 1, tolerance = 1e-12)
})

test_that("hpd_span conservatively spans multiple high-density modes", {
  d <- data.frame(x = 1:12,
                  y = c(5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3))
  p <- MMAPPR2:::.getPeakFromTopP(d, 0.80, method = "hpd_span")
  expect_equal(p$minPos, 1)
  expect_equal(p$maxPos, 12)
})

test_that("shortest_contiguous can isolate a dominant mode", {
  d <- data.frame(x = 1:12,
                  y = c(5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3))
  p <- MMAPPR2:::.getPeakFromTopP(d, 0.50,
                                  method = "shortest_contiguous")
  expect_lte(p$maxPos, 2)
})

test_that("peak apex uses the center of a flat top", {
  expect_equal(MMAPPR2:::.peakApexPosition(1:5, c(0, 2, 2, 2, 0)), 3)
  expect_true(is.na(MMAPPR2:::.peakApexPosition(1:2, c(NA, NA))))
})

test_that("legacy-current and global-SD peak cutoffs are explicit", {
  fake <- list(
    chr1 = list(loess = structure(list(fitted = c(1, 2, 3)), class = "loess")),
    chr2 = list(loess = structure(list(fitted = c(2, 4, 6)), class = "loess"))
  )
  legacy <- MMAPPR2:::.calculatePeakCutoff(fake, "legacy_current", k = 3)
  expect_equal(legacy$center, 3)
  expect_equal(legacy$spread, sqrt(5 / 3))
  expect_equal(legacy$cutoff, 3 + 3 * sqrt(5 / 3))

  global <- MMAPPR2:::.calculatePeakCutoff(fake, "global_sd", k = 3)
  all_values <- c(1, 2, 3, 2, 4, 6)
  expect_equal(global$center, stats::median(all_values))
  expect_equal(global$spread, stats::sd(all_values))
  expect_equal(global$cutoff,
               stats::median(all_values) + 3 * stats::sd(all_values))
})

test_that("derived chromosome seed is deterministic and chromosome-specific", {
  a1 <- MMAPPR2:::.derivedSeed(1, "chr1")
  a2 <- MMAPPR2:::.derivedSeed(1, "chr1")
  b <- MMAPPR2:::.derivedSeed(1, "chr2")
  expect_identical(a1, a2)
  expect_false(identical(a1, b))
})
