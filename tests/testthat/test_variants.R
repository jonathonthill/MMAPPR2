context("Variant calling")

testGenome <- GmapGenome(FastaFile("test_data/ref_genome_fastas/chr5.fa"), name = "dan_rerio_chr5", create = TRUE)

test_that("GmapGenome was created properly", {
  expect_equal(class(testGenome), "GmapGenome")
})
