library("testthat")

if (F){
  md <- readRDS("test_data/intermediate_MDs/post_loess.RDS")
  
  test_that("MapperData Read Properly", {
    expect_s4_class(md, "MmapprData")
  })
  
  md <- PrePeak(md)
  md <- PeakRefinement(md)
  
  test_that("PrePeakTest", {
    expect_equal(length(names(md@peaks)), 1)
  })
  
  test_that("PeakRefinementTest", {
    expect_true(md@peaks$chr5$peakPosition>35200000)
    expect_true(md@peaks$chr5$peakPosition<36000000)
  })
  
}
