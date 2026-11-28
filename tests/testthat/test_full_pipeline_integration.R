test_that("a real synthetic analysis traverses every public pipeline stage and retains the planted candidate", {
  fx <- .make_signal_fixture()
  p <- .signal_param(exportAiccPlots = TRUE)
  md <- mmapprData(p)

  md <- calculateDistance(md)
  expect_true(is.list(snpDistance(md)$chr1))
  expect_gte(nrow(snpDistance(md)$chr1$distanceDf), 50)

  md <- loessFit(md)
  expect_s3_class(snpDistance(md)$chr1$loess, "loess")
  expect_true(is.data.frame(snpDistance(md)$chr1$aicc))

  md <- prePeak(md)
  expect_true("chr1" %in% names(peaks(md)))

  md <- peakRefinement(md)
  expect_lte(peaks(md)$chr1$start, fx$candidate_pos)
  expect_gte(peaks(md)$chr1$end, fx$candidate_pos)

  md <- generateCandidates(md)
  expect_true(any(BiocGenerics::start(candidates(md)$snps$chr1) == fx$candidate_pos))
  planted <- candidates(md)$snps$chr1[BiocGenerics::start(candidates(md)$snps$chr1) == fx$candidate_pos]
  expect_true(as.numeric(S4Vectors::mcols(planted)$mutAltFreq) > .8)
  expect_equal(as.numeric(S4Vectors::mcols(planted)$wtAltFreq), .4, tolerance = .001)

  expect_identical(outputMmapprData(md), md)
  expected <- c("genome_plots.pdf", "peak_plots.pdf", "aicc_plots.pdf",
                "AllCandidateVariantsForchr1.tsv", "DetectedMutationsForchr1.tsv",
                "DifferentiallyExpressedGenesForchr1.tsv")
  expect_true(all(file.exists(file.path(outputFolder(p), expected))))
  expect_true(all(file.info(file.path(outputFolder(p), expected))$size > 0))
})

test_that("the public mmappr wrapper completes a real synthetic run and serializes an equivalent result", {
  fx <- .make_signal_fixture()
  p <- .signal_param(exportAiccPlots = FALSE)
  expect_message(md <- mmappr(p), "Welcome to MMAPPR2")
  expect_s4_class(md, "MmapprData")
  expect_true("chr1" %in% names(peaks(md)))
  expect_true(any(BiocGenerics::start(candidates(md)$snps$chr1) == fx$candidate_pos))
  saved <- readRDS(file.path(outputFolder(p), "mmappr_data.RDS"))
  expect_s4_class(saved, "MmapprData")
  expect_equal(BiocGenerics::start(candidates(saved)$snps$chr1),
               BiocGenerics::start(candidates(md)$snps$chr1))
  expect_true(file.exists(file.path(outputFolder(p), "genome_plots.pdf")))
  expect_true(file.exists(file.path(outputFolder(p), "peak_plots.pdf")))
  expect_false(file.exists(file.path(outputFolder(p), "aicc_plots.pdf")))

  log <- readLines(file.path(outputFolder(p), "mmappr2.log"))
  boundary <- which(log == "===== FULL LOG BELOW (chronological) =====")
  expect_length(boundary, 1L)
  expect_identical(log[1L], "===== MMAPPR2 SUMMARY =====")
  expect_length(which(log == "Peak Summary:"), 1L)
  expect_lt(which(log == "Peak Summary:"), boundary)
  expect_lt(which(log == "seqname: chr1"), boundary)
  expect_lt(which(grepl("^resampling seed:", log)), boundary)
  expect_lt(which(grepl("^End time:", log)), boundary)
  expect_gt(which(log == "Writing output plots and tables..."), boundary)
})
