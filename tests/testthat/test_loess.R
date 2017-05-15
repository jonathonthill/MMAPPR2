context("Loess fit/AICc optimization")

md <- readRDS("test_data/intermediate_MDs/post_file_read.RDS")
SkipDebug <- function() if(T) skip("skipping to save time")

test_that("normal chromosome performs fit correctly", {
  chr5List <- LoessFitForChr(md@distance$chr5, 
                             loessOptResolution = md@param@loessOptResolution,
                             loessOptCutFactor = md@param@loessOptCutFactor)
  expect_type(chr5List, "list")
  expect_true(all(c("wtCounts", "mutCounts", "loess", "aicc") %in% names(chr5List)))
  expect_gt(length(chr5List$loess$x), 5000)
})

test_that("Empty chromosome is skipped", {
  chr4error <- LoessFitForChr(md@distance$chr4,
                              loessOptResolution = md@param@loessOptResolution,
                              loessOptCutFactor = md@param@loessOptCutFactor)
  expect_type(chr4error, "character")
})

test_that("LoessFit runs properly for whole mmapprData", {
  SkipDebug()
  md2 <- LoessFit(md)
  successes <- sapply(md2@distance, function(seq) class(seq) == "list")
  expect_equal(length(successes), 26)
  expect_equal(sum(successes), 1)
  
  
  saveRDS(md2, "test_data/intermediate_MDs/post_loess.RDS")
})

