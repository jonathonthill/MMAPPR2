.peak_log_fields <- function(seqname) {
  c(paste0("seqname: ", seqname),
    "refined interval: 468-541",
    "full-data LOESS apex: 500",
    "resampling-density apex: 516",
    "resampling success rate: 1",
    "resampling seed: 863")
}

.reorder_test_log <- function(lines) {
  td <- tempfile("reorder_log_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))
  path <- file.path(td, "mmappr2.log")
  writeLines(lines, path)
  MMAPPR2:::.reorderLogFile(td)
  readLines(path)
}

.log_summary_header <- "===== MMAPPR2 SUMMARY ====="
.log_full_header <- "===== FULL LOG BELOW (chronological) ====="
.log_footer <- c("End time: 2026-10-23 10:00:00",
                 "MMAPPR2 runtime: 1 secs", "", "sessionInfo()", "R version 4.5.0")

test_that("log reordering moves the entire current peak summary for one or more chromosomes", {
  before <- c("Start time: 2026-10-23 09:59:59",
              "Getting variants in refined peak interval(s)")
  after <- c("Predicting coding variant effects",
             "Summarizing expression changes in peak interval(s)",
             "Scoring candidates by peak position", "Writing output plots and tables...", "")
  for (n in c(1L, 3L)) {
    peak <- c("Peak Summary:", .peak_log_fields("chr1"))
    if (n > 1L) {
      for (i in 2:n)
        peak <- c(peak, "", .peak_log_fields(paste0("chr", i)))
    }
    actual <- .reorder_test_log(c(before, peak, "", after, .log_footer))
    expect_identical(actual, c(.log_summary_header, "", peak, "", .log_footer, "",
                               .log_full_header, before, "", after))
    expect_equal(sum(actual == "Peak Summary:"), 1L)
    expect_equal(sum(grepl("^seqname:", actual)), n)
  }
})

test_that("log reordering keeps errors and unrelated messages in chronological order", {
  before <- c("Start time: 2026-10-23 09:59:59", "Getting variants in refined peak interval(s)")
  peak <- c("Peak Summary:", .peak_log_fields("chr1"), "", "seqname: chr2",
            "refined interval: 40-70")
  after <- c("ERROR: synthetic candidate failure", "MmapprData object up to the failing step is returned.", "")
  actual <- .reorder_test_log(c(before, peak, after, .log_footer))
  expect_identical(actual, c(.log_summary_header, "", peak, "", .log_footer, "",
                             .log_full_header, before, after))
})

test_that("log reordering handles failure before any peak summary is written", {
  before <- c("Start time: 2026-10-23 09:59:59", "ERROR: No peak regions identified", "")
  actual <- .reorder_test_log(c(before, .log_footer))
  expect_identical(actual, c(.log_summary_header, "", .log_footer, "", .log_full_header, before))
})

test_that("log reordering handles an empty peak heading and a missing runtime footer", {
  before <- c("Getting variants in refined peak interval(s)", "Peak Summary:", "ERROR: synthetic failure")
  expect_identical(.reorder_test_log(before), c(.log_summary_header, "", .log_full_header, before))
  peak <- c("Peak Summary:", .peak_log_fields("chr1"))
  actual <- .reorder_test_log(c("Start time: 2026-10-23 09:59:59", peak))
  expect_identical(actual, c(.log_summary_header, "", peak, "", .log_full_header,
                             "Start time: 2026-10-23 09:59:59"))
})
