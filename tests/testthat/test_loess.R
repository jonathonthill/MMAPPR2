context("Loess fit/AICc optimization")

md <- readRDS("test_data/intermediate_MDs/post_file_read.RDS")

test_that("normal chromosome performs fit correctly", {
  chr5List <- LoessFitForChr(md@distance$chr5, 
                             loessOptResolution = md@param@loessOptResolution,
                             loessOptCutFactor = md@param@loessOptCutFactor)
  expect_true(all(c("wtCounts", "mutCounts", "loess", "aicc") %in% names(chr5List)))
})

test_that("error chromosome is skipped", {
  
})

test_that("whole mmapprData is properly processed", {
  md2 <- LoessFit(md)
  saveRDS(md2, "test_data/intermediate_MDs/post_loess.RDS")
})

