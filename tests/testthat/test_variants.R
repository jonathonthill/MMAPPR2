context("Variant calling")

refGenFasta <- FastaFile(normalizePath("./test_data/ref_genome_fastas/chr5.fa"))
testGenome <- GmapGenome(refGenFasta, name = "dan_rerio_chr5", create = TRUE)

vepParam <- VEPParam(scriptPath = "./test_data/variant_effect_predictor_86.pl",
                     input = c(species = "danio_rerio_merged", format = "vcf"),
                     database = c(database = TRUE)
                     )
mmapprData <- readRDS("test_data/intermediate_MDs/post_peak.RDS")


test_that("GmapGenome was created properly", {
  expect_s4_class(testGenome, "GmapGenome")
})

mmapprData@peaks$chr5 <- AddPeakRange(mmapprData@peaks$chr5)
  
test_that("ranges for peaks are prepared", {
  chr5PeakList <- mmapprData@peaks$chr5
  expect_true("range" %in% names(chr5PeakList))
  expect_equal(width(chr5PeakList$range), chr5PeakList$end - chr5PeakList$start + 1)
})

mmapprData@peaks$chr5 <- AddVariantsForPeak(mmapprData@peaks$chr5)

test_that("variants are called for peak", {
  chr5PeakList <- mmapprData@peaks$chr5
  expect_true("variants" %in% names(chr5PeakList))
  expect_s4_class(chr5PeakList, "VRanges")
})


