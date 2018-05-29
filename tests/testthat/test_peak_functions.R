context("Peak functions")

if (F){
    md <- readRDS("test_data/objects/post_loess_dummy_md.RDS")
    md <- prePeak(md)
    md <- peakRefinement(md)
    
    test_that("no peaks are found for dummy files", {
        expect_equal(length(names(md@peaks)), 0)
    })
    
    md <- readRDS("test_data/objects/post_loess_dummy_md.RDS")
    md@distance <- readRDS('test_data/objects/post_loess_chr5_list.RDS')
    md <- prePeak(md)
    
    test_that('peak is identified for chromosome 5', {
        expect_equal(length(names(md@peaks)), 1)
        expect_equal(names(md@peaks)[0], 'chr5')
    })
    
    md <- peakRefinement(md)
    
    
    test_that("peakRefinementTest", {
        expect_true(md@peaks$chr5$peakPosition>19000000)
        expect_true(md@peaks$chr5$peakPosition<25000000)
    })
    
}
