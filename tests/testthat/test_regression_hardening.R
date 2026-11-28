test_that("annotation gene helper preserves explicit gene_id and gene_name", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  # Add names to a dedicated GTF so the assertion is independent of the base fixture.
  gtf <- file.path(fx$dir, "named.gtf")
  writeLines(c(
    'chr1\tMMAPPR2\tgene\t1\t150\t.\t+\t.\tgene_id "g1"; gene_name "GeneOne";',
    'chr1\tMMAPPR2\ttranscript\t1\t150\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";'
  ), gtf)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = gtf,
                   outputFolder = file.path(fx$dir, "out_named"),
                   includeScaffolds = TRUE, minDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)
  tx <- MMAPPR2:::.buildTxDb(p)
  on.exit(MMAPPR2:::.disconnectTxDb(tx), add = TRUE)
  genes <- MMAPPR2:::.annotationGeneRanges(p, tx)
  expect_equal(as.character(S4Vectors::mcols(genes)$gene_id), "g1")
  expect_equal(as.character(S4Vectors::mcols(genes)$gene_name), "GeneOne")
  expect_equal(BiocGenerics::start(genes), 1L)
  expect_equal(BiocGenerics::end(genes), 150L)
})

test_that("transient TxDb connections can be closed deterministically", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = file.path(fx$dir, "out_txdb"),
                   includeScaffolds = TRUE, minDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)
  tx <- MMAPPR2:::.buildTxDb(p)
  con <- BiocGenerics::dbconn(tx)
  expect_true(DBI::dbIsValid(con))
  MMAPPR2:::.disconnectTxDb(tx)
  expect_false(DBI::dbIsValid(con))
})


test_that("annotation gene helper normalizes unknown strand to unstranded", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  gtf <- file.path(fx$dir, "unstranded.gtf")
  writeLines(
    'chr1\tMMAPPR2\tgene\t10\t90\t.\t.\t.\tgene_id "gdot"; gene_name "DotGene";',
    gtf
  )
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = gtf,
                   outputFolder = file.path(fx$dir, "out_unstranded"),
                   includeScaffolds = TRUE, minDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)
  genes <- MMAPPR2:::.annotationGeneRanges(p)
  expect_equal(as.character(GenomicRanges::strand(genes)), "*")
  expect_equal(as.character(S4Vectors::mcols(genes)$gene_id), "gdot")
})


test_that("annotation gene helper preserves GFF3 ID and Name attributes", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  gff <- file.path(fx$dir, "named.gff3")
  writeLines(c(
    "##gff-version 3",
    "chr1\tMMAPPR2\tgene\t5\t120\t.\t+\t.\tID=gff_gene_1;Name=GffGeneOne"
  ), gff)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = gff,
                   outputFolder = file.path(fx$dir, "out_gff_named"),
                   includeScaffolds = TRUE, minDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)
  genes <- MMAPPR2:::.annotationGeneRanges(p)
  expect_equal(as.character(S4Vectors::mcols(genes)$gene_id), "gff_gene_1")
  expect_equal(as.character(S4Vectors::mcols(genes)$gene_name), "GffGeneOne")
  expect_equal(BiocGenerics::start(genes), 5L)
  expect_equal(BiocGenerics::end(genes), 120L)
})


test_that("annotation attribute parser tolerates harmless surrounding whitespace", {
  id <- MMAPPR2:::.annotationGeneIdentity(c(
    '   gene_id "g1";   gene_name "Gene One";',
    ' ID = g2 ; Name = GeneTwo '
  ))
  expect_identical(id$gene_id, c("g1", "g2"))
  expect_identical(id$gene_name, c("Gene One", "GeneTwo"))
})


test_that("lightweight annotation preflight ignores whitespace-prefixed comments", {
  path <- tempfile(fileext = ".gff3")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "   # deliberately indented comment",
    "chr1\tMMAPPR2\tgene\t10\t25\t.\t+\t.\tID=g1"
  ), path)
  observed <- MMAPPR2:::.annotationObservedRanges(path)
  expect_identical(observed$seqname, "chr1")
  expect_equal(observed$max_end, 25)
})


test_that("BAM preflight rejects style-only sequence-name matches", {
  refSi <- GenomeInfoDb::Seqinfo(seqnames = "chr1", seqlengths = 1000L)
  bamSi <- GenomeInfoDb::Seqinfo(seqnames = "1", seqlengths = 1000L)
  expect_error(
    MMAPPR2:::.compareBamSeqinfoFrozen("Mutant BAM 1", bamSi, refSi),
    "sequence naming.*must match exactly"
  )

  exactSi <- GenomeInfoDb::Seqinfo(seqnames = "chr1", seqlengths = 1000L)
  expect_identical(
    MMAPPR2:::.compareBamSeqinfoFrozen("Mutant BAM 1", exactSi, refSi),
    "chr1"
  )
})
