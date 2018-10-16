context("Peak functions")

md <- readRDS("test_data/objects/post_loess_dummy_md.RDS")
md <- prePeak(md)
md <- peakRefinement(md)

test_that("no peaks are found for dummy files", {
    expect_equal(length(names(md@peaks)), 0)
})

test_that('prePeak handles failed Loess', {
    md <- readRDS('test_data/objects/post_loess_dummy_md.RDS')
    md@distance$chr1 = list(wtCounts=data.frame(), mutCounts=data.frame(),
                            loess=structure('loess failed', class='try-error'))
    expect_equal(length(names(md@peaks)), 0)
})


test_that('peak is identified with prePeak', {
    md <- readRDS('test_data/objects/post_loess_dummy_md.RDS')
    md@distance$chr5 <- list(
        loess=structure(
            list(fitted=c(rep(0, 99), 100), x=seq_len(100)), class='loess'
        )
    )
    md <- prePeak(md)
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
    # return max pos of subsample without performing loess fit
    tempData <- rawData[sample(1:nrow(rawData), size=nrow(rawData)/2),]
    return(tempData$pos[which.max(tempData$euclideanDistance)])
}

### There was an integration test here, but it used big data which we
### wanted to take out of the package. Look at
### v0.98.10 if you want to revive it in the future.
