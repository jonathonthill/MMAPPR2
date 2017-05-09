library(MMAPPR2)
library(gmapR)
library(Rsamtools)
library(rtracklayer)

context("BAM file reading")


param <- MmapprParam(new("GmapGenome"), "test_data/bam_files/zy14_wt_cut.bam", "test_data/bam_files/zy14_mut_cut.bam",
                    vepParam = VEPParam())
mmapprData <- new("MmapprData", param = param)
mmapprData <- ReadInFiles(mmapprData)

test_that("correct data was produced", {
  expect_equal(nrow(mmapprData@distance$wtCounts) > 0, TRUE)
  expect_equal(nrow(mmapprData@distance$mutCounts) > 0, TRUE)
  expect_equal(nrow(mmapprData@distance$distanceDf) > 0, TRUE)
})