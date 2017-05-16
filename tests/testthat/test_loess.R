context("Loess fit/AICc optimization")

md <- readRDS("test_data/intermediate_MDs/post_file_read.RDS")
md@param@loessOptResolution <- 0.001
md@param@loessOptCutFactor <- 0.01
SkipDebug <- function() if(F) skip("skipping to save time")

test_that(".MinTwo works right", {
  x <- c(1, 2, 3)
  expect_equal(.MinTwo(x), c(1, 2))
  x <- c(1, 3, 2)
  expect_equal(.MinTwo(x), c(1, 2))
  x <- c(3, 2, 1)
  expect_equal(.MinTwo(x), c(1, 2))
  x <- c(1.00000001, 1.000001, 1.0001)
  expect_identical(.MinTwo(x), c(1.00000001, 1.000001))
  x <- c(-1, -2, -3)
  expect_identical(.MinTwo(x), c(-3, -2))
})

test_that(".NumDecimals works right", {
  expect_equal(.NumDecimals(.005), 3)
  expect_equal(.NumDecimals(.05), 2)
  expect_equal(.NumDecimals(.5), 1)
  expect_equal(.NumDecimals(5), 0)
  expect_equal(.NumDecimals(50), 0)
})

test_that(".LocalMin works right", {
  expect_equal(.LocalMin(c(1, 3, 1, 3)), c(1, 1))
  expect_equal(.LocalMin(c(1, 3, 2, 3)), c(1, 2))
  expect_equal(.LocalMin(c(1, 3, 2.999999, 3)), c(1, 2.999999))
  expect_equal(.LocalMin(c(1, 3, 2, 3, 1)), c(1, 2, 1))
  expect_equal(.LocalMin(c(-1, -3, -2, -3, -1)), c(-3, -3))
})

test_that("normal chromosome performs fit correctly", {
  SkipDebug()
  chr5List <- .LoessFitForChr(md@distance$chr5, 
                             loessOptResolution = md@param@loessOptResolution,
                             loessOptCutFactor = md@param@loessOptCutFactor)
  expect_type(chr5List, "list")
  expect_true(all(c("wtCounts", "mutCounts", "loess", "aicc") %in% names(chr5List)))
  expect_gt(length(chr5List$loess$x), 5000)
  expect_equal(chr5List$loess$span, 0.03)
})

test_that("Empty chromosome is skipped", {
  chr4error <- .LoessFitForChr(md@distance$chr4,
                              loessOptResolution = md@param@loessOptResolution,
                              loessOptCutFactor = md@param@loessOptCutFactor)
  expect_type(chr4error, "character")
})

test_that("LoessFit runs properly for whole mmapprData", {
  SkipDebug()
  md2 <- LoessFit(md, silent=T)
  successes <- sapply(md2@distance, function(seq) class(seq) == "list")
  str(md2@distance$chr5)
  expect_equal(length(successes), 26)
  expect_equal(sum(successes), 1)
  
  
  saveRDS(md2, "test_data/intermediate_MDs/post_loess.RDS")
})

