test_that("spliced RNA-seq reads survive mapping pileup", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = file.path(fx$dir, "out"),
                   includeScaffolds = TRUE,
                   minDepth = 1, candidateMinDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 200))
  x <- MMAPPR2:::.getPileup(wtFiles(p)[[1]], p, region)
  expect_true(any(x$POS == 10L))
  expect_true(any(x$POS == 120L))
  expect_true(any(x$POS == 50L)) # duplicate status remains unspecified
  # Alternative alignments must not be counted in addition to the
  # primary alignment. These positions exist only in the synthetic secondary/
  # supplementary records above.
  expect_false(any(x$POS == 160L))
  expect_false(any(x$POS == 180L))
})

test_that("secondary and supplementary alignments are excluded from candidate pileup", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 200))
  fa <- Rsamtools::FaFile(fx$fasta)
  x <- MMAPPR2:::.pooledBasePileup(
    bams = fx$bam, genome = fa, which = region,
    minBaseQuality = 0, minMapQuality = 0, maxDepth = 1000
  )
  expect_true(any(x$pos == 10L))
  expect_true(any(x$pos == 120L))
  expect_true(any(x$pos == 50L)) # duplicate status remains unspecified
  expect_false(any(x$pos == 160L))
  expect_false(any(x$pos == 180L))
})


test_that("duplicate status remains unspecified in the shared read filter", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  sb <- Rsamtools::scanBam(
    fx$bam,
    param = Rsamtools::ScanBamParam(
      flag = MMAPPR2:::.primaryMappedScanFlag(),
      what = "qname"
    )
  )[[1L]]$qname
  expect_true("rdup" %in% sb)
  expect_false("r2" %in% sb)
  expect_false("r3" %in% sb)
})

test_that("overlapping annotated genes are reduced before BAM queries", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = file.path(fx$dir, "out"),
                   includeScaffolds = TRUE,
                   minDepth = 1, candidateMinDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)
  chr <- MMAPPR2:::.getFileReadChrList(p)$chr1
  expect_true(GenomicRanges::isDisjoint(chr, ignore.strand = TRUE))
  expect_equal(length(chr), 1L)
  expect_equal(BiocGenerics::start(chr), 1L)
  expect_equal(BiocGenerics::end(chr), 200L)
})
