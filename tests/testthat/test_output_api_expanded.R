test_that("default output folder names are unique and do not consume the RNG stream", {
  set.seed(42)
  before <- .Random.seed
  a <- MMAPPR2:::.defaultOutputFolder()
  b <- MMAPPR2:::.defaultOutputFolder()
  expect_identical(.Random.seed, before)
  expect_false(identical(a, b))
  expect_match(a, "^mmappr2_[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2}_")
})

test_that("output folder preparation creates logs and requires explicit overwrite for prior results", {
  td <- tempfile("prepare_output_")
  out <- MMAPPR2:::.prepareOutputFolder(td)
  expect_true(dir.exists(out))
  expect_true(file.exists(file.path(out, "mmappr2.log")))
  writeLines("keep", file.path(out, "result.txt"))
  expect_error(MMAPPR2:::.prepareOutputFolder(out, overwrite = FALSE), "not empty")
  expect_silent(out2 <- MMAPPR2:::.prepareOutputFolder(out, overwrite = TRUE))
  expect_false(file.exists(file.path(out2, "result.txt")))
  expect_true(file.exists(file.path(out2, "mmappr2.log")))
})

test_that("output folder preparation refuses filesystem roots and the current working directory", {
  expect_error(MMAPPR2:::.prepareOutputFolder(getwd(), overwrite = TRUE), "protected")
  root <- normalizePath("/", mustWork = TRUE)
  expect_error(MMAPPR2:::.prepareOutputFolder(root, overwrite = TRUE), "protected")
})

test_that("genome plot rejects absent or incomplete fitted chromosome data", {
  p <- .test_param()
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = "failed")
  expect_error(MMAPPR2:::.plotGenomeDistance(md), "No successfully fitted")
  md@snpDistance <- list(chr1 = list(dummy = TRUE))
  expect_error(MMAPPR2:::.plotGenomeDistance(md), "missing loess")
})

test_that("genome, peak, and AICc plotters write nonempty PDFs for a populated result", {
  p <- .signal_param(exportAiccPlots = TRUE)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  md@peaks <- list(chr1 = .make_peak_entry("chr1", 300, 700, 500))
  expect_silent(MMAPPR2:::.plotGenomeDistance(md))
  expect_silent(MMAPPR2:::.plotPeaks(md))
  expect_silent(MMAPPR2:::.plotAicc(md))
  for (nm in c("genome_plots.pdf", "peak_plots.pdf", "aicc_plots.pdf")) {
    path <- file.path(outputFolder(p), nm)
    expect_true(file.exists(path), info = nm)
    expect_true(file.info(path)$size > 0, info = nm)
  }
})

test_that("AICc export warns cleanly when no finite search results exist", {
  p <- .test_param(exportAiccPlots = TRUE)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = list(aicc = data.frame(spans = .1, aiccValues = NA_real_)))
  expect_warning(expect_null(MMAPPR2:::.plotAicc(md)), "no finite AICc")
  expect_false(file.exists(file.path(outputFolder(p), "aicc_plots.pdf")))
})

test_that("peak plotting is a no-op when no refined peaks exist", {
  md <- mmapprData(.test_param())
  expect_invisible(MMAPPR2:::.plotPeaks(md))
})

test_that("outputMmapprData writes configured plots and returns the same object invisibly", {
  p <- .signal_param(exportAiccPlots = TRUE)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  md@peaks <- list(chr1 = .make_peak_entry("chr1", 300, 700, 500))
  out <- outputMmapprData(md)
  expect_identical(out, md)
  expect_true(file.exists(file.path(outputFolder(p), "genome_plots.pdf")))
  expect_true(file.exists(file.path(outputFolder(p), "peak_plots.pdf")))
  expect_true(file.exists(file.path(outputFolder(p), "aicc_plots.pdf")))
})

test_that("outputMmapprData closes only graphics devices it opened after a plotting error", {
  p <- .test_param()
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1"))
  userPdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(userPdf)
  userDev <- grDevices::dev.cur()
  on.exit({
    if (!is.null(grDevices::dev.list()) && userDev %in% grDevices::dev.list()) grDevices::dev.off(userDev)
  }, add = TRUE)

  testthat::local_mocked_bindings(
    .plotGenomeDistance = function(...) {
      grDevices::pdf(tempfile(fileext = ".pdf"))
      stop("synthetic plotting failure")
    },
    .package = "MMAPPR2"
  )
  expect_error(outputMmapprData(md), "synthetic plotting failure")
  expect_true(userDev %in% grDevices::dev.list())
  expect_identical(grDevices::dev.cur(), userDev)
})

test_that("TSV sanitization handles all-NA list cells and carriage returns", {
  x <- data.frame(id = 1:2)
  x$vals <- I(list(c(NA_character_, NA_character_), c("a\rb", "c\td")))
  safe <- MMAPPR2:::.tsvSafeDataFrame(x)
  expect_true(is.na(safe$vals[1]))
  expect_identical(safe$vals[2], "a b,c d")
})

test_that("atomic RDS save rejects a missing destination directory", {
  path <- file.path(tempfile("missing_parent_"), "x.rds")
  expect_error(MMAPPR2:::.atomicSaveRDS(1, path), "destination directory does not exist")
})

test_that("outputMmapprData recreates a removed configured output directory", {
  p <- .test_param()
  md <- mmapprData(p)
  outDir <- outputFolder(p)
  unlink(outDir, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(outDir))
  out <- outputMmapprData(md)
  expect_identical(out, md)
  expect_true(dir.exists(outDir))
})
