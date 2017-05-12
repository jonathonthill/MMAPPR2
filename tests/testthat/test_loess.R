context("Loess fit/AICc optimization")

md <- readRDS("test_data/intermediate_MDs/post_file_read.RDS")

test_that("normal chromosome performs fit correctly", {
  chr5List <- LoessFitForChr(md@distance$chr5, 
                             loessOptResolution = md@param@loessOptResolution,
                             loessOptCutFactor = md@param@loessOptCutFactor)
  expect_true(all(c("wtCounts", "mutCounts", "loess", "aicc") %in% names(chr5List)))
})

test_that("error chromosome is skipped", {
  chr4error <- LoessFitForChr(md@distance$chr4,
                              loessOptResolution = md@param@loessOptResolution,
                              loessOptCutFactor = md@param@loessOptCutFactor))

})

test_that("LoessFit runs properly for whole mmapprData", {
  md2 <- LoessFit(md)
  successes <- sapply(md2@distance, function(seq) class(seq) == "list")
  expect_equal(length(successes), 26)
  expect_equal(sum(successes), 1)
  
  
  saveRDS(md2, "test_data/intermediate_MDs/post_loess.RDS")
})

