test_that("subsample LOESS peak helper handles insufficient and ordinary data", {
  small <- data.frame(pos = 1:4, euclideanDistance = 1:4)
  expect_true(is.na(MMAPPR2:::.getSubsampleLoessMax(small, .5)))
  x <- 1:80
  d <- data.frame(pos = x, euclideanDistance = exp(-((x - 40) / 12)^2) + (x %% 3) / 100)
  set.seed(1)
  got <- MMAPPR2:::.getSubsampleLoessMax(d, .4)
  expect_true(is.finite(got))
  expect_true(got >= min(x) && got <= max(x))
})

test_that("density point mass validates ordering and nonempty mass", {
  expect_identical(MMAPPR2:::.densityPointMass(1, 5), 1)
  expect_error(MMAPPR2:::.densityPointMass(c(1, 1), c(1, 1)), "strictly increasing")
  expect_error(MMAPPR2:::.densityPointMass(c(2, 1), c(1, 1)), "strictly increasing")
  expect_error(MMAPPR2:::.densityPointMass(c(1, 2), c(0, 0)), "mass is empty")
  expect_equal(sum(MMAPPR2:::.densityPointMass(c(1, 2, 3), c(-1, 2, 1))), 1)
})

test_that("peak interval extraction sorts input, discards nonfinite rows, and supports full mass", {
  d <- data.frame(x = c(3, 1, 2, NA), y = c(1, 1, 2, 4))
  out <- MMAPPR2:::.getPeakFromTopP(d, 1, "hpd_span")
  expect_equal(out$minPos, 1)
  expect_equal(out$maxPos, 3)
  expect_equal(out$peakPos, 2)
  expect_error(MMAPPR2:::.getPeakFromTopP(data.frame(x = 1:2, y = c(0, 0)), .8), "empty")
})

test_that("global-SD cutoff requires at least two finite fitted values", {
  fit <- stats::loess(c(1, 1.1, 1.2, 1.3, 1.4) ~ c(1, 2, 3, 4, 5), span = 1)
  fit$fitted[] <- NA_real_
  fit$fitted[1] <- 1
  expect_error(MMAPPR2:::.calculatePeakCutoff(list(chr1 = list(loess = fit)), "global_sd", 3),
               "Insufficient finite")
})

test_that("peak cutoff ignores failed chromosomes and nonfinite fitted values", {
  a <- .make_loess_result()$loess
  a$fitted[1:2] <- NA_real_
  out <- MMAPPR2:::.calculatePeakCutoff(list(chr1 = list(loess = a), chr2 = "failed"),
                                        "legacy_current", 2)
  expect_true(all(is.finite(unlist(out[c("cutoff", "center", "spread")]))))
  expect_identical(out$method, "legacy_current")
  expect_error(MMAPPR2:::.calculatePeakCutoff(list(chr1 = "failed"), "legacy_current", 2),
               "No finite LOESS")
})

test_that("prePeak directly identifies qualifying chromosomes and records cutoff metadata", {
  p <- .signal_param(peakCutoffSd = 0)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", peak = TRUE),
                         chr2 = .make_loess_result("chr2", peak = FALSE))
  out <- prePeak(md)
  expect_true("chr1" %in% names(peaks(out)))
  expect_false("chr2" %in% names(peaks(out)))
  pk <- peaks(out)$chr1
  expect_true(all(c("seqname", "cutoff", "cutoffCenter", "cutoffSpread", "cutoffMethod") %in% names(pk)))
})

test_that("prePeak uses strict above-cutoff comparison and skips short chromosomes", {
  p <- .signal_param()
  md <- mmapprData(p)
  equalFit <- .make_loess_result("chr1")$loess
  equalFit$fitted[] <- 1
  shortFit <- .make_loess_result("chr2", n = 20)$loess
  shortFit$fitted[] <- 2
  md@snpDistance <- list(chr1 = list(loess = equalFit), chr2 = list(loess = shortFit))
  testthat::local_mocked_bindings(
    .calculatePeakCutoff = function(...) list(cutoff = 1, center = 1, spread = 0, method = "legacy_current"),
    .package = "MMAPPR2"
  )
  out <- prePeak(md)
  expect_length(peaks(out), 0)
})

test_that("peak refinement is deterministic and restores the caller RNG state", {
  p <- .signal_param(peakResampleIterations = 20, randomSeed = 123)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  input <- list(seqname = "chr1", cutoff = .2, cutoffCenter = .1,
                cutoffSpread = .03, cutoffMethod = "legacy_current")

  set.seed(999)
  before <- .Random.seed
  a <- MMAPPR2:::.peakRefinementChr(input, md)
  expect_identical(.Random.seed, before)
  b <- MMAPPR2:::.peakRefinementChr(input, md)
  expect_equal(a$start, b$start)
  expect_equal(a$end, b$end)
  expect_equal(a$densityPeakPosition, b$densityPeakPosition)
  expect_identical(a$resampleSeed, b$resampleSeed)
  expect_true(a$resampleSuccessRate >= .5)
})

test_that("peak refinement restores absence of a global RNG seed", {
  p <- .signal_param(peakResampleIterations = 10)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  input <- list(seqname = "chr1", cutoff = .2, cutoffCenter = .1,
                cutoffSpread = .03, cutoffMethod = "legacy_current")
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) saved <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had) assign(".Random.seed", saved, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  MMAPPR2:::.peakRefinementChr(input, md)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("peak refinement rejects too few markers and too many failed resamples", {
  p <- .signal_param(peakResampleIterations = 10)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 9))
  input <- list(seqname = "chr1", cutoff = .2, cutoffCenter = .1,
                cutoffSpread = .03, cutoffMethod = "legacy_current")
  expect_error(MMAPPR2:::.peakRefinementChr(input, md), "Too few finite markers")

  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  testthat::local_mocked_bindings(.getSubsampleLoessMax = function(...) NA_real_, .package = "MMAPPR2")
  expect_error(MMAPPR2:::.peakRefinementChr(input, md), "Too many failed")
})

test_that("degenerate peak resampling produces a one-point physical interval", {
  p <- .signal_param(peakResampleIterations = 10)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  input <- list(seqname = "chr1", cutoff = .2, cutoffCenter = .1,
                cutoffSpread = .03, cutoffMethod = "legacy_current")
  testthat::local_mocked_bindings(.getSubsampleLoessMax = function(...) 500, .package = "MMAPPR2")
  out <- MMAPPR2:::.peakRefinementChr(input, md)
  expect_equal(out$start, 500L)
  expect_equal(out$end, 500L)
  expect_equal(out$densityFunction(c(499, 500, 501)), c(0, 1, 0))
})

test_that("peakRefinement directly populates every initialized peak", {
  p <- .signal_param(peakResampleIterations = 10)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  md@peaks <- list(chr1 = list(seqname = "chr1", cutoff = .2, cutoffCenter = .1,
                               cutoffSpread = .03, cutoffMethod = "legacy_current"))
  out <- peakRefinement(md)
  expect_true(all(c("start", "end", "densityFunction", "loessPeakPosition",
                    "densityPeakPosition", "resampleSuccessRate", "resampleSeed") %in%
                  names(peaks(out)$chr1)))
  expect_true(peaks(out)$chr1$start >= 1)
  expect_true(peaks(out)$chr1$end <= 1200)
})

test_that("public peak refinement is reproducible across serial and SOCK parallel backends", {
  testthat::skip_if(Sys.getenv("MMAPPR2_RUN_PARALLEL_TESTS") != "true",
                    "parallel backend validation is run as a separate installed-package CI job")
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)

  p <- .signal_param(peakResampleIterations = 20, randomSeed = 314)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  md@peaks <- list(chr1 = list(seqname = "chr1", cutoff = .2, cutoffCenter = .1,
                               cutoffSpread = .03, cutoffMethod = "legacy_current"))

  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  serial <- peakRefinement(md)
  BiocParallel::register(BiocParallel::SnowParam(workers = 2, type = "SOCK"), default = TRUE)
  parallel <- peakRefinement(md)

  a <- peaks(serial)$chr1
  b <- peaks(parallel)$chr1
  for (nm in c("start", "end", "peakPosition", "densityPeakPosition",
               "loessPeakPosition", "resampleSuccessRate", "resampleSeed",
               "intervalMethod")) {
    expect_equal(a[[nm]], b[[nm]], info = nm)
  }
  expect_equal(a$densityData$x, b$densityData$x)
  expect_equal(a$densityData$y, b$densityData$y)
})
