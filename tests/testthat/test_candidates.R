context("Variant calling and effect prediction")

mmapprData <- readRDS("test_data/objects/post_peakref_chr5_md.RDS")


test_that("ranges for peaks are prepared", {
    candList <- .getPeakRange(mmapprData@peaks$chr5)
    expect_s4_class(candList, 'GRanges')
    chr5PeakList <- mmapprData@peaks$chr5
    expect_equal(BiocGenerics::width(candList),
                 chr5PeakList$end - chr5PeakList$start + 1)
})



test_that('.getVariantsForRange returns VRanges with sample names', {
    mockVariantCall <- mockery::mock(
        VariantAnnotation::VRanges(seqnames='chr5',
                                   ranges=IRanges::IRanges(start=c(1,2,3), width=1),
                                   sampleNames='zy14_mut_cut_filt.bam',
                                   ref=c('A', 'A', 'A')
        )
    )
    
    inputRange <-
        GenomicRanges::GRanges('chr5',
                               IRanges::IRanges(start = 1, width = 75682077))
    
    mockery::stub(
        .getVariantsForRange,
        'callSampleSpecificVariants',
        mockVariantCall
    )
    mockery::stub(.getVariantsForRange, 'TallyVariantsParam', 42)
    result <- .getVariantsForRange(inputRange, mmapprData@param)
    
    mockery::expect_called(mockVariantCall, 1)
    expect_s4_class(result, 'VRanges')
    expect_equal(as.character(GenomicRanges::seqnames(result))[1], 'chr5')
    expect_equal(as.character(Biobase::sampleNames(result))[1],
                 'zy14_mut_cut_filt.bam')
})


# TODO: get cached variant GRanges saved
# TODO: test filter function
# TODO: test order function