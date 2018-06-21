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

test_that('.getPeakFromTopP works on double peak', {
    data <- data.frame(pos=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                       density=c(0.25, 0.25, 1.5, 2, 0, 1, 3, 2, 0, 0)/10)
    stopifnot(sum(data$density) == 1)
    result <- .getPeakFromTopP(data, 0.95)
    expect_identical(result, list(minPos=3, maxPos=8, peakPos=7))
})

test_that('.getPeakFromTopP works on edge peak', {
    data <- data.frame(pos=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                       density=c(0.25, 0.25, 1.5, 2, 0, 1, 0, 2, 0, 3)/10)
    stopifnot(sum(data$density) == 1)
    result <- .getPeakFromTopP(data, 0.95)
    expect_identical(result, list(minPos=3, maxPos=10, peakPos=10))
})

.mockSubsampleLoessMax <- function(rawData, loessSpan) {
    # return max pos without performing loess fit
    tempData <- rawData[sample(1:nrow(rawData), size=nrow(rawData)/2),]
    return(tempData$pos[which.max(tempData$euclideanDistance)])
}

# saves time since the dplyr::arrange call is super slow
.mockPeakFromTopP <- function(data, topP) {
    names(data) <- c('x', 'y')
    nrows <- nrow(data)
    result <- list()
    result$peakPos = data$x[which.max(data$y)]
    result$minPos = result$peakPos - (nrows/15)
    result$maxPos = result$peakPos + (nrows/15)
    return(result)
}

test_that(".peakRefinementChr works right on chr 5", {
    inputList <- list(seqname='chr5')
    chr5list <- mockr::with_mock(
        .getSubsampleLoessMax=.mockSubsampleLoessMax,
        .getPeakFromTopP=.mockPeakFromTopP,
        .env=as.environment('package:MMAPPR2'),
        .peakRefinementChr(inputList, md)
    )
    expect_true(all(c('start', 'end', 'densityFunction',
                      'peakPosition', 'seqname') %in%
                        names(chr5list)))
    expect_true(chr5list$peakPosition > 19000000)
    expect_true(chr5list$peakPosition < 25000000)
})

