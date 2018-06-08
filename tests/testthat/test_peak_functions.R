context("Peak functions")

md <- readRDS("test_data/objects/post_loess_dummy_md.RDS")
md <- prePeak(md)
md <- peakRefinement(md)

test_that("no peaks are found for dummy files", {
    expect_equal(length(names(md@peaks)), 0)
})

md <- readRDS("test_data/objects/post_loess_chr5_md.RDS")
md <- prePeak(md)

test_that('peak is identified for chromosome 5', {
    expect_equal(length(names(md@peaks)), 1)
    expect_equal(names(md@peaks)[1], 'chr5')
})


test_that(".peakRefinementChr works right on chr 5", {
    inputList <- list(seqname='chr5')
    chr5list <- .peakRefinementChr(inputList, md)
    expect_true(all(c('start', 'end', 'densityFunction',
                      'peakPosition', 'seqname') %in%
                        names(chr5list)))
    expect_true(chr5list$peakPosition > 19000000)
    expect_true(chr5list$peakPosition < 25000000)
})

