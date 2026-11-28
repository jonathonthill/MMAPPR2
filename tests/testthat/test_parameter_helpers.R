.base_scalar_args <- function() {
  list(
    includeScaffolds = FALSE, minDepth = 20, homozygoteCutoff = .95,
    minBaseQuality = 20, minMapQuality = 20, fileAggregation = "simple",
    distancePower = 4, peakIntervalWidth = .80, loessOptResolution = .001,
    loessOptCutFactor = .1, maxPileupDepth = 1000, candidateMinDepth = 1,
    candidateMinAltDepth = 2, candidateMinAltFreq = .80,
    candidateMaxWtAltFreq = 1, candidateMinDeltaAF = 0,
    candidatePoolMode = "auto", candidateChunkSize = 0,
    peakCutoffSd = 3, peakCutoffMethod = "legacy_current",
    peakIntervalMethod = "hpd_span", peakResampleIterations = 1000,
    randomSeed = 1, pairedEnd = FALSE, ignoreStrand = FALSE,
    expressionPseudocount = 0.01, exportAiccPlots = FALSE, nMutFiles = 1L
  )
}


test_that("FASTA validity helper rejects empty and malformed paths", {
  expect_false(isTRUE(MMAPPR2:::.validFastaFile(character())))
  expect_false(isTRUE(MMAPPR2:::.validFastaFile(NA_character_)))
  expect_false(isTRUE(MMAPPR2:::.validFastaFile("")))

  missing <- tempfile(fileext = ".fa")
  expect_false(isTRUE(MMAPPR2:::.validFastaFile(missing)))
})

test_that("valid scalar parameter set passes", {
  expect_length(do.call(MMAPPR2:::.validateScalarValues, .base_scalar_args()), 0)
})

test_that("impossible depth combinations are rejected before I/O", {
  a <- .base_scalar_args(); a$minDepth <- 1001
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, a), collapse = " "),
               "minDepth cannot exceed maxPileupDepth")

  a <- .base_scalar_args(); a$candidateMinAltDepth <- 251
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, a), collapse = " "),
               "candidateMinAltDepth exceeds")
})

test_that("strict candidate ALT threshold rejects an impossible cutoff of one", {
  a <- .base_scalar_args(); a$candidateMinAltFreq <- 1
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, a), collapse = " "),
               "candidateMinAltFreq")
  a$candidateMinAltFreq <- 0.999
  expect_length(do.call(MMAPPR2:::.validateScalarValues, a), 0)
})

test_that("invalid scientific-method names are rejected", {
  a <- .base_scalar_args(); a$peakCutoffMethod <- "mystery"
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, a), collapse = " "),
               "peakCutoffMethod")
  a <- .base_scalar_args(); a$peakIntervalMethod <- "mystery"
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, a), collapse = " "),
               "peakIntervalMethod")
})

test_that("candidate memory controls are validated", {
  a <- .base_scalar_args(); a$candidatePoolMode <- "mystery"
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, a), collapse = " "),
               "candidatePoolMode")
  a <- .base_scalar_args(); a$candidateChunkSize <- -1
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, a), collapse = " "),
               "candidateChunkSize")
})


test_that("invalid scalar settings fail before creating the output directory", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  out <- file.path(fx$dir, "must_not_be_created")
  expect_error(
    mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                refFasta = fx$fasta, gtf = fx$gtf,
                outputFolder = out,
                includeScaffolds = TRUE,
                minDepth = 1001, maxPileupDepth = 1000),
    "minDepth cannot exceed maxPileupDepth"
  )
  expect_false(dir.exists(out))
})

test_that("refFasta replacement keeps the cached FaFile synchronized", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = file.path(fx$dir, "out"),
                   includeScaffolds = TRUE,
                   minDepth = 1, candidateMinDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)

  fa2 <- file.path(fx$dir, "ref2.fa")
  writeLines(c(">chr1", paste(rep("C", 250), collapse = "")), fa2)
  Rsamtools::indexFa(fa2)
  refFasta(p) <- fa2

  expect_identical(refFasta(p), normalizePath(fa2))
  expect_identical(normalizePath(BiocGenerics::path(p@refGenome)), normalizePath(fa2))
})

test_that("replacing MmapprData parameters invalidates derived results", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = file.path(fx$dir, "out"),
                   includeScaffolds = TRUE,
                   minDepth = 1, candidateMinDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = list(dummy = TRUE))
  md@peaks <- list(chr1 = list(seqname = "chr1"))
  md@candidates <- list(snps = list(chr1 = GenomicRanges::GRanges()))

  expect_warning(param(md) <- p, "invalidates existing derived results")
  expect_length(snpDistance(md), 0)
  expect_length(peaks(md), 0)
  expect_length(candidates(md), 0)
})


test_that("BAM validation rejects files that merely have BAM/BAI-looking names", {
  td <- tempfile("bad_bam_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  bam <- file.path(td, "broken.bam")
  bai <- paste0(bam, ".bai")
  writeLines("not a BAM", bam)
  writeLines("not an index", bai)
  bfl <- Rsamtools::BamFileList(Rsamtools::BamFile(bam, index = bai))
  result <- MMAPPR2:::.validBamFiles(bfl, deep = TRUE)
  expect_false(isTRUE(result))
  expect_match(paste(result, collapse = " "), "unreadable BAM")
})


test_that("invalid BAM input cannot clear an existing output directory", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  bad <- file.path(fx$dir, "broken_input.bam")
  badIndex <- paste0(bad, ".bai")
  writeLines("not a BAM", bad)
  writeLines("not an index", badIndex)

  out <- file.path(fx$dir, "existing_results")
  dir.create(out)
  sentinel <- file.path(out, "keep_me.txt")
  writeLines("important prior result", sentinel)

  expect_error(
    mmapprParam(wtFiles = bad, mutFiles = bad,
                refFasta = fx$fasta, gtf = fx$gtf,
                outputFolder = out, overwrite = TRUE,
                includeScaffolds = TRUE,
                minDepth = 1, candidateMinDepth = 1,
                minBaseQuality = 0, minMapQuality = 0),
    "unreadable BAM"
  )
  expect_true(file.exists(sentinel))
  expect_identical(readLines(sentinel), "important prior result")
})


test_that("public defaults include safe candidate pooling, AICc export, and the expression pseudocount", {
  f <- formals(MMAPPR2::mmapprParam)
  expect_equal(eval(f$candidateMinDepth), 1)
  expect_identical(eval(f$candidatePoolMode)[1], "auto")
  expect_equal(eval(f$candidateChunkSize), 0)
  expect_false(eval(f$ignoreStrand))
  expect_equal(eval(f$expressionPseudocount), 0.01)
  expect_false(eval(f$exportAiccPlots))
})

test_that("BAM/reference mismatch cannot clear prior output", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)

  wrongFa <- file.path(fx$dir, "wrong.fa")
  writeLines(c(">chr1", paste(rep("A", 249), collapse = "")), wrongFa)
  Rsamtools::indexFa(wrongFa)

  out <- file.path(fx$dir, "existing_reference_results")
  dir.create(out)
  sentinel <- file.path(out, "keep_me.txt")
  writeLines("important prior result", sentinel)

  expect_error(
    mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                refFasta = wrongFa, gtf = fx$gtf,
                outputFolder = out, overwrite = TRUE,
                includeScaffolds = TRUE, minDepth = 1,
                minBaseQuality = 0, minMapQuality = 0),
    "sequence-length conflict"
  )
  expect_true(file.exists(sentinel))
})

test_that("annotation bounds mismatch cannot clear prior output", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  badGtf <- file.path(fx$dir, "bad.gtf")
  writeLines('chr1\tMMAPPR2\tgene\t1\t300\t.\t+\t.\tgene_id "too_long";', badGtf)

  out <- file.path(fx$dir, "existing_annotation_results")
  dir.create(out)
  sentinel <- file.path(out, "keep_me.txt")
  writeLines("important prior result", sentinel)

  expect_error(
    mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                refFasta = fx$fasta, gtf = badGtf,
                outputFolder = out, overwrite = TRUE,
                includeScaffolds = TRUE, minDepth = 1,
                minBaseQuality = 0, minMapQuality = 0),
    "Annotation contains coordinates outside"
  )
  expect_true(file.exists(sentinel))
})


test_that("resource setters preserve BAM/FASTA/annotation concordance", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = file.path(fx$dir, "out_setter_preflight"),
                   includeScaffolds = TRUE, minDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0)

  wrongFa <- file.path(fx$dir, "setter_wrong.fa")
  writeLines(c(">chr1", paste(rep("A", 249), collapse = "")), wrongFa)
  Rsamtools::indexFa(wrongFa)
  expect_error(refFasta(p) <- wrongFa, "sequence-length conflict")
  expect_identical(refFasta(p), normalizePath(fx$fasta))

  badGtf <- file.path(fx$dir, "setter_bad.gtf")
  writeLines('chr1\tMMAPPR2\tgene\t1\t300\t.\t+\t.\tgene_id "too_long";', badGtf)
  expect_error(gtf(p) <- badGtf, "outside the indexed reference FASTA")
  expect_identical(gtf(p), normalizePath(fx$gtf))
})


test_that("sidecar index freshness rejects missing, empty, and stale indexes", {
  td <- tempfile("index_freshness_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  data <- file.path(td, "input.dat")
  idx <- file.path(td, "input.dat.idx")
  writeBin(as.raw(1:10), data)

  expect_false(MMAPPR2:::.indexIsCurrent(data, idx))
  file.create(idx)
  expect_false(MMAPPR2:::.indexIsCurrent(data, idx))

  writeBin(as.raw(1:5), idx)
  now <- Sys.time()
  Sys.setFileTime(data, now)
  Sys.setFileTime(idx, now - 10)
  expect_false(MMAPPR2:::.indexIsCurrent(data, idx))

  Sys.setFileTime(idx, now + 10)
  expect_true(MMAPPR2:::.indexIsCurrent(data, idx))
})
