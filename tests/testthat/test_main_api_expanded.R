test_that("logging helpers append messages and capture non-character objects", {
  td <- tempfile("logging_"); dir.create(td)
  file.create(file.path(td, "mmappr2.log"))
  expect_message(MMAPPR2:::.messageAndLog("first", td), "first")
  MMAPPR2:::.log(list(a = 1), td)
  MMAPPR2:::.log("last", td)
  txt <- readLines(file.path(td, "mmappr2.log"))
  expect_true(any(txt == "first"))
  expect_true(any(grepl("a", txt)))
  expect_true(any(txt == "last"))
})

test_that("mmappr orchestrates public stages in the documented order and saves final state", {
  p <- .test_param()
  calls <- character()
  add <- function(x) calls <<- c(calls, x)
  testthat::local_mocked_bindings(
    calculateDistance = function(md) { add("calculateDistance"); md@snpDistance <- list(chr1 = list(dummy = TRUE)); md },
    loessFit = function(md) { add("loessFit"); md },
    prePeak = function(md) { add("prePeak"); md@peaks <- list(chr1 = list(seqname = "chr1")); md },
    peakRefinement = function(md) { add("peakRefinement"); md@peaks$chr1 <- .make_peak_entry(); md },
    generateCandidates = function(md) { add("generateCandidates"); md@candidates <- list(snps = list(chr1 = GenomicRanges::GRanges())); md },
    outputMmapprData = function(md) { add("outputMmapprData"); invisible(md) },
    .package = "MMAPPR2"
  )
  expect_message(out <- mmappr(p), "Welcome to MMAPPR2")
  expect_identical(calls, c("calculateDistance", "loessFit", "prePeak", "peakRefinement", "generateCandidates", "outputMmapprData"))
  expect_s4_class(out, "MmapprData")
  saved <- file.path(outputFolder(p), "mmappr_data.RDS")
  expect_true(file.exists(saved))
  expect_s4_class(readRDS(saved), "MmapprData")
  log <- readLines(file.path(outputFolder(p), "mmappr2.log"))
  expect_true(any(grepl("sessionInfo", log, fixed = TRUE)))
  expect_true(any(grepl("MMAPPR2 runtime", log, fixed = TRUE)))
})

test_that("mmappr returns and serializes the most recently completed partial state after failure", {
  p <- .test_param()
  testthat::local_mocked_bindings(
    calculateDistance = function(md) { md@snpDistance <- list(chr1 = list(marker = "distance complete")); md },
    loessFit = function(md) stop("synthetic loess failure"),
    .package = "MMAPPR2"
  )
  expect_message(out <- mmappr(p), "synthetic loess failure")
  expect_identical(snpDistance(out)$chr1$marker, "distance complete")
  saved <- readRDS(file.path(outputFolder(p), "mmappr_data.RDS"))
  expect_identical(snpDistance(saved)$chr1$marker, "distance complete")
  log <- paste(readLines(file.path(outputFolder(p), "mmappr2.log")), collapse = "\n")
  expect_match(log, "ERROR: synthetic loess failure")
  expect_match(log, "failing step is returned")
})

test_that("mmappr treats absence of peaks as a recoverable pipeline failure", {
  p <- .test_param()
  testthat::local_mocked_bindings(
    calculateDistance = function(md) { md@snpDistance <- list(chr1 = list(dummy = TRUE)); md },
    loessFit = function(md) md,
    prePeak = function(md) { md@peaks <- list(); md },
    .package = "MMAPPR2"
  )
  expect_message(out <- mmappr(p), "No peak regions identified")
  expect_length(peaks(out), 0)
  expect_true(file.exists(file.path(outputFolder(p), "mmappr_data.RDS")))
})
