library(ensemblVEP)

context("Main and helper functions")

vepParam <- VEPParam(input = c(format="vcf"))

test_that("RepList function works", {
    param = expect_warning(MmapprParam(refGenome = new("GmapGenome"), "yo", "hi", VEPParam()))
    expect_equal(length(RepList(param, 5)), 5)
    expect_equal(length(RepList(param, 0)), 0)
    expect_equal(length(RepList(param, 1)), 1)
    expect_equal(length(RepList(param, 2)), 2)
    expect_error(length(RepList(param, -1)))
    x <- RepList(param, 2)
    expect_identical(x[[1]], x[[2]])
})

test_that("MmapprParam takes character, BamFile, or BamFileList", {
    param = expect_warning(MmapprParam(new("GmapGenome"), "wt", 'mut', vepParam))
    expect_s4_class(param, "MmapprParam")
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")
    
    bf <- BamFile("tests/testthat/test_data/bam_files/zy14_mut_chr5.bam")
    param = expect_warning(MmapprParam(new("GmapGenome"), bf, bf, vepParam))
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")
    
    bfl <- BamFileList("tests/testthat/test_data/bam_files/zy14_mut_chr5.bam")
    param = expect_warning(MmapprParam(new("GmapGenome"), bfl, bfl, vepParam))
    expect_s4_class(param@wtFiles, "BamFileList")
    expect_s4_class(param@mutFiles, "BamFileList")
    
})