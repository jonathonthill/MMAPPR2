test_that("sequence-style helpers harmonize compatible Seqinfo objects and tolerate no style hint", {
  ref <- GenomeInfoDb::Seqinfo("chr1", seqlengths = 100)
  ann <- GenomeInfoDb::Seqinfo("1", seqlengths = 100)
  expect_identical(MMAPPR2:::.choose_target_style(ref), "UCSC")
  got <- MMAPPR2:::.harmonizeSeqinfoFrozen(ann, ref)
  expect_true("chr1" %in% GenomeInfoDb::seqlevels(got))
  gr <- GenomicRanges::GRanges("1", IRanges::IRanges(1, 10))
  styled <- MMAPPR2:::.setSeqlevelsStyleFrozen(gr, "UCSC")
  expect_identical(as.character(GenomeInfoDb::seqnames(styled)), "chr1")
  expect_identical(MMAPPR2:::.setSeqlevelsStyleFrozen(gr, NA_character_), gr)
})

test_that("FASTA seqinfo helper reads indexed lengths", {
  fx <- .test_fixture()
  si <- MMAPPR2:::.faSeqinfo(Rsamtools::FaFile(fx$fasta))
  expect_identical(GenomeInfoDb::seqlevels(si), "chr1")
  expect_equal(as.numeric(GenomeInfoDb::seqlengths(si)), 250)
})

test_that("peak range conversion preserves refined coordinates", {
  gr <- MMAPPR2:::.getPeakRange(list(seqname = "chr2", start = 10.2, end = 20.8))
  expect_s4_class(gr, "GRanges")
  expect_identical(as.character(GenomeInfoDb::seqnames(gr)), "chr2")
  expect_identical(BiocGenerics::start(gr), 10L)
  expect_identical(BiocGenerics::end(gr), 20L)
})

test_that("cgroup directory and available-memory helpers always return usable planning values", {
  for (v in c("v1", "v2")) {
    dirs <- MMAPPR2:::.processCgroupMemoryDirs(v)
    expect_type(dirs, "character")
    expect_gte(length(dirs), 1)
    expect_identical(anyDuplicated(dirs), 0L)
  }
  mem <- MMAPPR2:::.availableMemoryBytes()
  expect_length(mem, 1)
  expect_true(is.na(mem) || (is.finite(mem) && mem > 0))
})

test_that("candidate one-range and range dispatcher agree on a small real interval", {
  p <- .signal_param(candidatePoolMode = "memory")
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(490, 510))
  one <- MMAPPR2:::.candidateSnvsOneRange(gr, p, "memory")
  dispatched <- MMAPPR2:::.candidateSnvsForRange(gr, p)
  expect_equal(as.data.frame(one), as.data.frame(dispatched))
  expect_equal(length(dispatched), 1)
})

test_that("WT chunk support returns indexed evidence only for overlapping variants", {
  p <- .signal_param(candidatePoolMode = "chunked", candidateChunkSize = 20)
  vars <- MMAPPR2:::.candidateSnvsOneRange(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(490, 510)), p, "incremental")
  support <- MMAPPR2:::.wtSupportChunkSafe(vars, p,
                                           GenomicRanges::GRanges("chr1", IRanges::IRanges(490, 510)))
  expect_equal(nrow(support), 1)
  expect_identical(support$index, 1L)
  expect_equal(support$total, 100L)
  none <- MMAPPR2:::.wtSupportChunkSafe(vars, p,
                                        GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10)))
  expect_equal(nrow(none), 0)
})

test_that("parameter validity helpers report corrupted scalar state", {
  p <- .test_param()
  expect_true(isTRUE(MMAPPR2:::.validMmapprParam(p)))
  bad <- p
  methods::slot(bad, "minDepth", check = FALSE) <- 0
  v <- MMAPPR2:::.validMmapprParam(bad)
  expect_false(isTRUE(v))
  expect_match(paste(v, collapse = " "), "minDepth")
  expect_error(MMAPPR2:::.validateAfterSet(bad), "minDepth")
})

test_that("BAM index helper builds a missing index and returns a BamFileList", {
  fx <- .test_fixture()
  td <- tempfile("reindex_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  bam <- file.path(td, "copy.bam")
  file.copy(fx$bam, bam)
  expect_message(bfl <- MMAPPR2:::.indexBamFileList(bam), "Indexing BAM")
  expect_s4_class(bfl, "BamFileList")
  idx <- c(paste0(bam, ".bai"), sub("\\.bam$", ".bai", bam))
  expect_true(any(file.exists(idx)))
})

test_that("TxDb construction helper builds a usable database from the annotation", {
  fx <- .make_signal_fixture()
  tx <- MMAPPR2:::.makeTxDbFromAnnotation(fx$gtf)
  on.exit(MMAPPR2:::.disconnectTxDb(tx), add = TRUE)
  expect_s4_class(tx, "TxDb")
  expect_gt(length(GenomicFeatures::genes(tx)), 0)
})

test_that("resource preflight and generic scalar setter helpers return validated parameter objects", {
  p <- .test_param()
  p2 <- MMAPPR2:::.preflightAfterResourceSet(p)
  expect_s4_class(p2, "MmapprParam")
  p3 <- MMAPPR2:::.setScalarSlot(p, "minDepth", 2)
  expect_identical(minDepth(p3), 2)
})
