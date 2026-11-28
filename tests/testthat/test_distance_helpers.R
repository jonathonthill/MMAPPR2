test_that("replicate aggregation counts missing-coverage files in average depth", {
  x <- data.table::data.table(
    CHROM = "chr1", POS = 10L,
    A.FREQ = 1, C.FREQ = 0, G.FREQ = 0, T.FREQ = 0,
    CVG = 20, FILE_ID = 1L
  )
  out <- MMAPPR2:::.avgFiles(x, fileAggregation = "simple", nFiles = 2L)
  expect_equal(out$AVE.CVG, 10)
  expect_equal(out$AVE.A.FREQ, 1)
})

test_that("weighted aggregation pools reads while simple weights covered replicates equally", {
  x <- data.table::data.table(
    CHROM = c("chr1", "chr1"), POS = c(10L, 10L),
    A.FREQ = c(1, 0), C.FREQ = c(0, 0), G.FREQ = c(0, 1), T.FREQ = c(0, 0),
    CVG = c(90, 10), FILE_ID = c(1L, 2L)
  )
  simple <- MMAPPR2:::.avgFiles(data.table::copy(x), "simple", 2L)
  weighted <- MMAPPR2:::.avgFiles(data.table::copy(x), "weighted", 2L)
  expect_equal(simple$AVE.A.FREQ, 0.5)
  expect_equal(weighted$AVE.A.FREQ, 0.9)
  expect_equal(simple$AVE.CVG, 50)
  expect_equal(weighted$AVE.CVG, 50)
})

test_that("WT homozygosity filter uses the four frequency columns, not colon arithmetic", {
  x <- data.table::data.table(
    AVE.A.FREQ.WT = c(.96, .50, NA),
    AVE.C.FREQ.WT = c(.04, .50, NA),
    AVE.G.FREQ.WT = c(0, 0, NA),
    AVE.T.FREQ.WT = c(0, 0, NA)
  )
  cols <- names(x)
  out <- MMAPPR2:::.wtHomozygousRows(x, cols, cutoff = .95)
  expect_identical(out, c(TRUE, FALSE, TRUE))
})

test_that("empty mapping pileup has a stable schema", {
  x <- MMAPPR2:::.emptyPileupTable()
  expect_named(x, c("CHROM", "POS", "A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ", "CVG"))
  expect_equal(nrow(x), 0)
})
