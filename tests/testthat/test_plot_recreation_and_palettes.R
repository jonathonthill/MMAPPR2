# Replot API: no analysis reruns, no table rewrites, and independent figure styling.

test_that("named figure palettes preserve default roles and validate overrides", {
  defaults <- MMAPPR2:::.plotPaletteColors()
  expect_identical(defaults[["peakFill"]], "#2ecc71")
  expect_identical(defaults[["genomeRaw"]], "#999999AA")
  expect_identical(defaults[["aiccSelected"]], "black")
  expect_identical(names(defaults), names(MMAPPR2:::.plotPaletteColors("grayscale")))
  expect_identical(names(defaults), names(MMAPPR2:::.plotPaletteColors("colorblind")))

  custom <- MMAPPR2:::.plotPaletteColors("colorblind",
                                        c(peakFill = "#11223399", density = "steelblue"))
  expect_identical(custom[["peakFill"]], "#11223399")
  expect_identical(custom[["density"]], "steelblue")
  expect_identical(custom[["fitted"]], "#222222")
  expect_error(MMAPPR2:::.plotPaletteColors("sunset"), "palette must be one of")
  expect_error(MMAPPR2:::.plotPaletteColors(NA_character_), "palette must be one of")
  expect_error(MMAPPR2:::.plotPaletteColors(colors = c(bogus = "red")), "Unknown plot color")
  expect_error(MMAPPR2:::.plotPaletteColors(colors = c(density = "not a color")), "Invalid R color")
  expect_error(MMAPPR2:::.plotPaletteColors(colors = c("red")), "named character")
  expect_error(MMAPPR2:::.plotPaletteColors(colors = c(fitted = "red", fitted = "blue")),
               "unique")
  expect_error(MMAPPR2:::.plotPaletteColors(colors = c(density = NA_character_)),
               "named character")
})

.test_replot_md <- function() {
  p <- .signal_param(exportAiccPlots = FALSE)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = .make_loess_result("chr1", n = 80))
  md@peaks <- list(chr1 = .make_peak_entry("chr1", 300, 700, 500))
  md@candidates <- list(snps = list(chr1 = data.frame(
    seqnames = "chr1", start = 500L, ref = "A", alt = "C"
  )))
  md
}

test_that("saved RDS replots all figures despite disabled original AICc export", {
  md <- .test_replot_md()
  rdsDir <- tempfile("replot_rds_")
  dir.create(rdsDir)
  on.exit(unlink(rdsDir, recursive = TRUE, force = TRUE), add = TRUE)
  rdsFile <- file.path(rdsDir, "mmappr_data.RDS")
  # Simulate original input resources being unavailable without touching
  # the shared synthetic-test files used by other test cases.
  md@param@refFasta <- file.path(rdsDir, "now_missing_reference.fa")
  md@param@gtf <- file.path(rdsDir, "now_missing_annotation.gtf")
  saveRDS(md, rdsFile)
  original <- readBin(rdsFile, "raw", n = file.size(rdsFile))

  # These operations are never allowed while replotting saved data.
  testthat::local_mocked_bindings(
    calculateDistance = function(...) stop("recomputation attempted"),
    loessFit = function(...) stop("recomputation attempted"),
    peakRefinement = function(...) stop("recomputation attempted"),
    generateCandidates = function(...) stop("recomputation attempted"),
    .writeCandidateTables = function(...) stop("TSV rewrite attempted"),
    .package = "MMAPPR2"
  )
  outputs <- plotMmapprData(rdsFile)
  expect_identical(names(outputs), c("genome", "peaks", "aicc"))
  expect_true(all(file.exists(outputs)))
  expect_true(all(file.info(outputs)$size > 0L))
  expect_identical(dirname(outputs[["genome"]]),
                   normalizePath(file.path(rdsDir, "replots"), mustWork = TRUE))
  expect_identical(readBin(rdsFile, "raw", n = file.size(rdsFile)), original)
  expect_length(list.files(rdsDir, pattern = "\\.tsv$", recursive = TRUE), 0)
  expect_false(file.exists(file.path(outputFolder(param(md)), "aicc_plots.pdf")))
})

test_that("selective replotting from MmapprData does not rewrite TSV or RDS", {
  md <- .test_replot_md()
  out <- tempfile("selected_figures_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("must not change", file.path(out, "AllCandidateVariantsForchr1.tsv"))
  saveRDS(list(tag = "untouched"), file.path(out, "mmappr_data.RDS"))

  outputs <- plotMmapprData(md, outputDir = out, plots = c("genome", "aicc"),
                            palette = "grayscale",
                            colors = c(aiccSelected = "#E69F00", fitted = "#224466"))
  expect_identical(names(outputs), c("genome", "aicc"))
  expect_true(all(file.exists(outputs)))
  expect_false(file.exists(file.path(out, "peak_plots.pdf")))
  expect_identical(readLines(file.path(out, "AllCandidateVariantsForchr1.tsv")),
                   "must not change")
  expect_identical(readRDS(file.path(out, "mmappr_data.RDS")), list(tag = "untouched"))
})

test_that("existing PDFs are protected by default and overwrite replaces only selection", {
  md <- .test_replot_md()
  out <- tempfile("replot_overwrite_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("old PDF placeholder", file.path(out, "genome_plots.pdf"))
  writeLines("other figure", file.path(out, "peak_plots.pdf"))
  writeLines("other results", file.path(out, "candidate.tsv"))
  expect_error(plotMmapprData(md, outputDir = out, plots = "genome"),
               "Refusing to overwrite")
  expect_identical(readLines(file.path(out, "genome_plots.pdf")), "old PDF placeholder")
  paths <- plotMmapprData(md, outputDir = out, plots = "genome", overwrite = TRUE)
  expect_true(file.info(paths[["genome"]])$size > 1000L)
  expect_identical(readLines(file.path(out, "peak_plots.pdf")), "other figure")
  expect_identical(readLines(file.path(out, "candidate.tsv")), "other results")
  expect_length(list.files(out, all.files = TRUE,
                           pattern = "\\.backup\\.|^\\.mmappr2_"), 0L)
})

test_that("a replot graphics failure preserves all existing PDFs and devices", {
  md <- .test_replot_md()
  out <- tempfile("replot_failure_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("old genome", file.path(out, "genome_plots.pdf"))
  writeLines("old peaks", file.path(out, "peak_plots.pdf"))

  userPdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(userPdf)
  userDev <- grDevices::dev.cur()
  on.exit({
    if (!is.null(grDevices::dev.list()) && userDev %in% grDevices::dev.list())
      grDevices::dev.off(userDev)
  }, add = TRUE)
  testthat::local_mocked_bindings(
    .plotPeaks = function(...) {
      grDevices::pdf(tempfile(fileext = ".pdf"))
      stop("simulated PDF failure")
    },
    .package = "MMAPPR2"
  )
  expect_error(plotMmapprData(md, outputDir = out, plots = c("genome", "peaks"),
                              overwrite = TRUE), "simulated PDF failure")
  expect_identical(readLines(file.path(out, "genome_plots.pdf")), "old genome")
  expect_identical(readLines(file.path(out, "peak_plots.pdf")), "old peaks")
  expect_length(list.files(out, all.files = TRUE, pattern = "^\\.mmappr2_"), 0L)
  expect_true(userDev %in% grDevices::dev.list())
  expect_identical(grDevices::dev.cur(), userDev)
})

test_that("a partial saved result plots available types and skips unavailable ones", {
  md <- .test_replot_md()
  md@snpDistance$failed <- "error: missing chromosome data"
  md@snpDistance$unfitted <- list(dummy = TRUE)
  md@peaks <- list()
  out <- tempfile("partial_replot_")
  expect_warning(paths <- plotMmapprData(md, outputDir = out), "Skipping unavailable figure")
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  expect_identical(names(paths), c("genome", "aicc"))
  expect_true(all(file.exists(paths)))
  expect_false(file.exists(file.path(out, "peak_plots.pdf")))

  md@snpDistance <- list(failed = "error")
  expect_error(suppressWarnings(plotMmapprData(md, outputDir = tempfile("no_figures_"))),
               "No requested figures")
})

test_that("invalid replot parameters fail without creating files", {
  md <- .test_replot_md()
  out <- tempfile("invalid_replot_")
  expect_error(plotMmapprData(md, outputDir = out, plots = "unknown"), "plots must contain")
  expect_error(plotMmapprData(md, outputDir = out, plots = character()), "plots must contain")
  expect_error(plotMmapprData(md, outputDir = out, palette = "holographic"),
               "palette must be")
  expect_error(plotMmapprData(md, outputDir = out, colors = c(peakFill = "bogus-color")),
               "Invalid R color")
  expect_error(plotMmapprData(md, outputDir = out, overwrite = NA), "overwrite must be")
  expect_error(plotMmapprData(tempfile("nonexistent_rds_")), "existing MMAPPR2 RDS")
  wrong <- tempfile(fileext = ".RDS")
  saveRDS(list(foo = 1), wrong)
  on.exit(unlink(wrong), add = TRUE)
  expect_error(plotMmapprData(wrong), "does not contain a MmapprData")
  expect_false(dir.exists(out))
})

test_that("normal output with styled figures still writes candidate tables", {
  md <- .test_replot_md()
  out <- outputFolder(param(md))
  expect_identical(outputMmapprData(md, palette = "colorblind",
                                    colors = c(peakFill = "#E69F00AA")), md)
  expect_true(file.exists(file.path(out, "genome_plots.pdf")))
  expect_true(file.exists(file.path(out, "peak_plots.pdf")))
  expect_false(file.exists(file.path(out, "aicc_plots.pdf")))
  expect_true(file.exists(file.path(out, "AllCandidateVariantsForchr1.tsv")))
})

test_that("normal output still writes candidate TSVs if figure rendering fails", {
  md <- .test_replot_md()
  out <- outputFolder(param(md))
  testthat::local_mocked_bindings(
    .plotGenomeDistance = function(...) stop("broken PDF device"),
    .package = "MMAPPR2"
  )
  expect_error(outputMmapprData(md), "broken PDF device")
  expect_true(file.exists(file.path(out, "AllCandidateVariantsForchr1.tsv")))
})

test_that("mmappr accepts styling without changing downstream processing", {
  p <- .test_param()
  args <- NULL
  testthat::local_mocked_bindings(
    calculateDistance = function(md) md,
    loessFit = function(md) md,
    prePeak = function(md) { md@peaks <- list(chr1 = list(seqname = "chr1")); md },
    peakRefinement = function(md) md,
    generateCandidates = function(md) md,
    outputMmapprData = function(md, palette = "default", colors = NULL) {
      args <<- list(palette = palette, colors = colors)
      invisible(md)
    },
    .package = "MMAPPR2"
  )
  md <- suppressMessages(mmappr(p, palette = "grayscale",
                               colors = c(fitted = "#010203")))
  expect_s4_class(md, "MmapprData")
  expect_identical(args, list(palette = "grayscale", colors = c(fitted = "#010203")))
  expect_true(file.exists(file.path(outputFolder(p), "mmappr_data.RDS")))
})

test_that("saved AICc evaluations are sufficient even without a fitted LOESS", {
  md <- mmapprData(.test_param(exportAiccPlots = FALSE))
  md@snpDistance <- list(chr1 = list(
    bestSpan = 0.2,
    aicc = data.frame(spans = c(.1, .2, .3), aiccValues = c(9, 2, 5))
  ))
  out <- tempfile("aicc_only_")
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  paths <- plotMmapprData(md, outputDir = out, plots = "aicc")
  expect_identical(names(paths), "aicc")
  expect_true(file.info(paths[["aicc"]])$size > 0L)
})

test_that("normal styled output respects existing AICc export parameter", {
  md <- .test_replot_md()
  md@param@exportAiccPlots <- TRUE
  expect_identical(outputMmapprData(md, palette = "grayscale",
                                    colors = c(aiccSelected = "#D55E00")), md)
  out <- outputFolder(param(md))
  expect_true(file.exists(file.path(out, "aicc_plots.pdf")))
  expect_true(file.exists(file.path(out, "AllCandidateVariantsForchr1.tsv")))
})

test_that("replotting never follows a symlink at a selected PDF destination", {
  md <- .test_replot_md()
  out <- tempfile("replot_symlink_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  sentinel <- file.path(out, "unrelated_data.txt")
  writeLines("unrelated and protected", sentinel)
  fig <- file.path(out, "genome_plots.pdf")
  if (!isTRUE(file.symlink(sentinel, fig))) skip("Symlinks not available on this filesystem")
  expect_error(plotMmapprData(md, outputDir = out, plots = "genome",
                              overwrite = TRUE), "symlink")
  expect_identical(readLines(sentinel), "unrelated and protected")
})

test_that("malformed AICc entries do not prevent rendering a valid chromosome", {
  md <- mmapprData(.test_param(exportAiccPlots = FALSE))
  md@snpDistance <- list(
    bad_spans = list(aicc = data.frame(spans = c("bad", "data"),
                                       aiccValues = c(2, 3))),
    good = list(aicc = data.frame(spans = c(0.1, 0.2, 0.3),
                                  aiccValues = c(4, 1, 5)),
                bestSpan = "invalid metadata"),
    bad_values = list(aicc = data.frame(spans = c(0.1, 0.2),
                                        aiccValues = c("1", "2")))
  )
  expect_identical(MMAPPR2:::.usableAiccSequences(md@snpDistance), "good")
  out <- tempfile("partial_aicc_")
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  pdfs <- plotMmapprData(md, outputDir = out, plots = "aicc")
  expect_identical(names(pdfs), "aicc")
  expect_true(file.info(pdfs[["aicc"]])$size > 0L)
})

test_that("replot preflight omits incomplete LOESS and peak entries", {
  md <- .test_replot_md()
  broken <- .make_loess_result("broken", n = 80)
  broken$loess$x[1] <- NA_real_  # An unusable genomic coordinate.
  md@snpDistance$broken <- broken
  md@peaks$broken <- .make_peak_entry("broken", 300, 700, 500)
  expect_identical(MMAPPR2:::.usableLoessSequences(md@snpDistance), "chr1")
  expect_identical(MMAPPR2:::.usablePeakSequences(md@peaks, "chr1"), "chr1")

  out <- tempfile("incomplete_loess_")
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  pdfs <- plotMmapprData(md, outputDir = out, plots = c("genome", "peaks"))
  expect_identical(names(pdfs), c("genome", "peaks"))
  expect_true(all(file.info(pdfs)$size > 0L))
})

test_that("PDF installation restores the original file if replacement fails", {
  out <- tempfile("plot_install_rollback_")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
  old <- file.path(out, "genome_plots.pdf")
  staged <- file.path(out, ".new_genome.pdf")
  writeLines("original PDF placeholder", old)
  writeLines("replacement PDF placeholder", staged)

  calls <- 0L
  failInstall <- function(from, to) {
    calls <<- calls + 1L
    if (calls == 2L) return(FALSE)
    file.rename(from, to)
  }
  expect_error(MMAPPR2:::.installPlotPdf(staged, old, .renameFile = failInstall),
               "previous file was restored")
  expect_identical(calls, 3L)
  expect_identical(readLines(old), "original PDF placeholder")
  expect_true(file.exists(staged))
  expect_length(list.files(out, all.files = TRUE, pattern = "\\.backup\\."), 0L)
})
