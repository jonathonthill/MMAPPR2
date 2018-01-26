library("testthat")

if (F){
    md <- readRDS("test_data/intermediate_MDs/post_loess.RDS")
    
    test_that("MapperData Read Properly", {
        expect_s4_class(md, "MmapprData")
    })
    
    md <- prePeak(md)
    md <- peakRefinement(md)
    
    test_that("prePeakTest", {
        expect_equal(length(names(md@peaks)), 1)
    })
    
    test_that("peakRefinementTest", {
        expect_true(md@peaks$chr5$peakPosition>19000000)
        expect_true(md@peaks$chr5$peakPosition<25000000)
    })
    
}
