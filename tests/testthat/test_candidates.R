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


test_that('.getVariantsForRange handles replicates', {
    skip_if_not_installed('mockery')
    
    mockVariantCall <- mockery::mock(
        VariantAnnotation::VRanges(seqnames='chr5',
                                   ranges=IRanges::IRanges(start=c(1), width=1),
                                   ref=c('A')
        )
    )
    param <- mmapprData@param
    bf <- Rsamtools::BamFile('test_data/bam_files/zy14_dummy.bam')
    mockMergeBam <- mockery::mock(Rsamtools::BamFile('tmp_wt.bam'),
                                  Rsamtools::BamFile('tmp_mut.bam'))
    wtFiles(param) <- c(bf, bf)
    mutFiles(param) <- c(bf, bf)
    inputRange <-
        GenomicRanges::GRanges('chr5',
                               IRanges::IRanges(start = 1, width = 75682077))
    mockery::stub(
        .getVariantsForRange,
        'callSampleSpecificVariants',
        mockVariantCall
    )
    mockery::stub(.getVariantsForRange, 'mergeBam', mockMergeBam)
    mockery::stub(.getVariantsForRange, 'TallyVariantsParam', 42)
    result <- .getVariantsForRange(inputRange, param)
    mockery::expect_called(mockVariantCall, 1)
    mockery::expect_called(mockMergeBam, 2)
})


test_that('.filterVariants removes variants with "LOW" impact', {
    gr <- GenomicRanges::GRanges('chr5', IRanges::IRanges(start=c(100, 200, 300), width=1),
                                 IMPACT=c('HIGH', 'LOW', NA))
    result <- .filterVariants(gr)
    expected <- GenomicRanges::GRanges('chr5', IRanges::IRanges(start=c(100, 300), width=1),
                                       IMPACT=c('HIGH', NA))
    expect_identical(result, expected)
})


test_that('.densityScoreAndOrderVariants works', {
    gr <- GenomicRanges::GRanges('chr5', IRanges::IRanges(start=c(100, 200, 300, 400, 500), width=1))
    densFunc <- function(pos) -(pos-310)^2
    result <- .densityScoreAndOrderVariants(gr, densFunc)
    expPositions <- c(300, 400, 200, 500, 100)
    expected <- GenomicRanges::GRanges('chr5', IRanges::IRanges(start=expPositions, width=1),
                                       peakDensity=vapply(expPositions, densFunc, 0))
    expect_identical(result, expected)
})
