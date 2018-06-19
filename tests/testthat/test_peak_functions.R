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

mockSubsampleLoessMax <- function(rawData, loessSpan) {
    # return max pos without performing loess fit
    tempData <- rawData[sample(1:nrow(rawData), size=nrow(rawData)/2),]
    return(tempData$pos[which.max(tempData$euclideanDistance)])
}


test_that(".peakRefinementChr works right on chr 5", {
    skip_if_not_travis_or_bioc()
    inputList <- list(seqname='chr5')
    chr5list <- with_mock(
        .getSubsampleLoessMax=mockSubsampleLoessMax,
        .peakRefinementChr(inputList, md)
    )
    expect_true(all(c('start', 'end', 'densityFunction',
                      'peakPosition', 'seqname') %in%
                        names(chr5list)))
    expect_true(chr5list$peakPosition > 19000000)
    expect_true(chr5list$peakPosition < 25000000)
})

