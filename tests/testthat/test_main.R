library(ensemblVEP)

context("Main and helper functions")

vepParam <- VEPParam(input = c(format="vcf"))
param <- MmapprParam(new("GmapGenome"), "./test_data/bam_files/zy14_wt_cut_filt.bam", 
                     "test_data/bam_files/zy14_mut_cut_filt.bam",
                     vepParam = vepParam)

### TODO: test .runFunctionInParallel

test_that(".repList function works", {
    expect_equal(length(.repList(param, 5)), 5)
    expect_equal(length(.repList(param, 0)), 0)
    expect_equal(length(.repList(param, 1)), 1)
    expect_equal(length(.repList(param, 2)), 2)
    expect_error(length(.repList(param, -1)))
    x <- .repList(param, 2)
    expect_identical(x[[1]], x[[2]])
})

test_that("MmapprParam takes character, BamFile, or BamFileList of real files", {
    fn_wt <- "test_data/bam_files/zy14_wt_cut_filt.bam"
    fn_mut <- "test_data/bam_files/zy14_mut_cut_filt.bam"
    
    # non-existing filename shouldn't work
    expect_error(MmapprParam(new("GmapGenome"), "wt", 'mut', vepParam))
    
    # existing filename, just character, should work
    param1 <- MmapprParam(new("GmapGenome"), fn_wt, fn_mut, vepParam)
    expect_s4_class(param1, "MmapprParam")
    expect_s4_class(param1@wtFiles, "BamFileList")
    expect_s4_class(param1@mutFiles, "BamFileList")
    
    bf <- BamFile("test_data/bam_files/zy14_mut_cut_filt.bam")
    expect_warning(param2 = MmapprParam(new("GmapGenome"), bf, bf, vepParam=vepParam))
    expect_s4_class(param2@wtFiles, "BamFileList")
    expect_s4_class(param2@mutFiles, "BamFileList")
    
    bfl <- BamFileList("test_data/bam_files/zy14_mut_cut_filt.bam")
    expect_warning(param3 = MmapprParam(new("GmapGenome"), bfl, bfl, vepParam=vepParam))
    expect_s4_class(param3@wtFiles, "BamFileList")
    expect_s4_class(param3@mutFiles, "BamFileList")
    
})