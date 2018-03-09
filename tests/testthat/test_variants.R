# context("Variant calling and effect prediction")
# 
# refGenFasta <- rtracklayer::FastaFile(normalizePath("./test_data/ref_genome_fastas/chr5.fa"))
# invisible(capture.output(testGenome <- 
#                              gmapR::GmapGenome(refGenFasta,
#                                                name = "dan_rerio_zv9_chr5", 
#                                                create = TRUE)))
# 
# mmapprData <- readRDS("test_data/intermediate_MDs/post_peak.RDS")
# mmapprData@param@vepParam <- ensemblVEP::VEPFlags(
#     flags = list(species = "danio_rerio_merged", format = "vcf",
#                  database = TRUE)
# )
# mmapprData@param@refGenome <- testGenome
# 
# 
# test_that("GmapGenome was created properly", {
#     expect_s4_class(testGenome, "GmapGenome")
# })
# 
# mmapprData@candidates$chr5 <- .getPeakRange(mmapprData@peaks$chr5)
# 
# test_that("ranges for peaks are prepared", {
#     chr5PeakList <- mmapprData@peaks$chr5
#     expect_equal(width(mmapprData@candidates$chr5), 
#                  chr5PeakList$end - chr5PeakList$start + 1)
# })
# 
# mmapprData@candidates$chr5 <- .getVariantsForRange(mmapprData@candidates$chr5,
#                                                   mmapprData@param)
# 
# test_that("variants are called for peak", {
#     expect_s4_class(mmapprData@candidates$chr5, "VRanges")
#     expect_gt(length(mmapprData@candidates$chr5), 0)
# })
# 
# 
# 
