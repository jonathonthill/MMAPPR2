context("Variant calling")
library(gmapR)
library(VariantTools)
library(rtracklayer)

refGenFasta <- FastaFile(normalizePath("./test_data/ref_genome_fastas/chr5.fa"))
#TODO get right ref genome for tests
testGenome <- GmapGenome(refGenFasta, name = "dan_rerio_zv9_chr5", create = TRUE)

mmapprData <- readRDS("test_data/intermediate_MDs/post_peak.RDS")
mmapprData@param@vepParam <- VEPParam(
    scriptPath = "./test_data/variant_effect_predictor_86.pl",
    input = c(species = "danio_rerio_merged", format = "vcf"),
    database = c(database = TRUE)
)
mmapprData@param@refGenome <- testGenome


test_that("GmapGenome was created properly", {
    expect_s4_class(testGenome, "GmapGenome")
})

mmapprData@candidates$chr5 <- GetPeakRange(mmapprData@peaks$chr5)

test_that("ranges for peaks are prepared", {
    chr5PeakList <- mmapprData@peaks$chr5
    expect_equal(width(mmapprData@candidates$chr5), 
                 chr5PeakList$end - chr5PeakList$start + 1)
})

mmapprData@candidates$chr5 <- GetVariantsForRange(mmapprData@candidates$chr5,
                                                  mmapprData@param)

test_that("variants are called for peak", {
    expect_s4_class(mmapprData@candidates$chr5, "VRanges")
    expect_gt(length(mmapprData@candidates$chr5), 0)
})



