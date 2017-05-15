library("testthat")


test_that("MapperData Read Properly", {
  md <- readRDS("tests/testthat/test_data/intermediate_MDs/post_loess.RDS")
  expect_s4_class(md, "MmapprData")
})

md <- PrePeak(md)
md <- PeakRefinement(md)

test_that("PrePeakTest", {
  expect_equal(length(names(md@peaks)), 1)
})

test_that("PeakRefinementTest", {
  expect_true(md@peaks$chr5$peakPosition>35200000)
  
})




