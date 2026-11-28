test_that("mutant candidate threshold preserves strict original >80 percent boundary", {
  keep <- MMAPPR2:::.candidateMutantKeep(
    totalDepth = c(20, 20, 19, 2, 20, 0),
    altDepth = c(16, 17, 16, 2, 1, 0),
    minDepth = 1, minAltDepth = 2, minAltFreq = .80
  )
  expect_identical(keep, c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("WT candidate filters can be disabled or applied independently", {
  wt <- c(.2, .9, NA)
  delta <- c(.7, .05, NA)
  expect_true(all(MMAPPR2:::.candidateWtKeep(wt, delta, 1, 0)))
  expect_identical(
    MMAPPR2:::.candidateWtKeep(wt, delta, .5, 0),
    c(TRUE, FALSE, FALSE)
  )
  expect_identical(
    MMAPPR2:::.candidateWtKeep(wt, delta, 1, .2),
    c(TRUE, FALSE, FALSE)
  )
})

test_that("base pileup conversion does not mutate input and preserves depths", {
  pile <- data.table::data.table(
    seqnames = c("chr1", "chr1"), pos = c(10L, 20L),
    A = c(2L, 0L), C = c(18L, 0L), G = c(0L, 20L), T = c(0L, 0L),
    totalDepth = c(20L, 20L), ref = c("A", "G")
  )
  before <- data.table::copy(pile)
  vr <- MMAPPR2:::.basePileupToVRanges(pile)
  expect_identical(pile, before)
  expect_false(".row_id" %in% names(pile))
  expect_length(vr, 1)
  expect_equal(as.character(VariantAnnotation::ref(vr)), "A")
  expect_equal(as.character(VariantAnnotation::alt(vr)), "C")
  expect_equal(as.integer(VariantAnnotation::refDepth(vr)), 2L)
  expect_equal(as.integer(VariantAnnotation::altDepth(vr)), 18L)
  expect_equal(as.integer(VariantAnnotation::totalDepth(vr)), 20L)
})

test_that("variant-to-pileup matching uses exact genomic coordinates", {
  vars <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 30), width = 1))
  pile <- data.table::data.table(seqnames = c("chr1", "chr1"), pos = c(30L, 10L))
  idx <- MMAPPR2:::.matchVariantToPileup(vars, pile)
  expect_identical(idx, c(2L, 1L))
})

test_that("expression group means cannot leak samples across groups", {
  mat <- matrix(c(10, 20, 100, 200,
                  30, 50, 300, 500), nrow = 2, byrow = TRUE)
  means <- MMAPPR2:::.poolMeanCounts(mat, num_wt = 2L, num_mut = 2L)
  expect_equal(means$wt, c(15, 40))
  expect_equal(means$mut, c(150, 400))
})

test_that("WT background metadata does not change default candidate ordering", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 20), width = 1))
  S4Vectors::mcols(gr)$peakDensity <- c(0.9, 0.5)
  S4Vectors::mcols(gr)$sharedBackgroundLike <- c(TRUE, FALSE)
  ordered <- MMAPPR2:::.orderVariants(gr)
  expect_equal(BiocGenerics::start(ordered)[1], 10)
})

test_that("reference range validation rejects coordinates beyond FASTA", {
  si <- GenomeInfoDb::Seqinfo("chr1", seqlengths = 100)
  ok <- GenomicRanges::GRanges("chr1", IRanges::IRanges(90, 100))
  bad <- GenomicRanges::GRanges("chr1", IRanges::IRanges(90, 101))
  expect_invisible(MMAPPR2:::.validateRangesWithinReference(ok, si, "test"))
  expect_error(MMAPPR2:::.validateRangesWithinReference(bad, si, "test"),
               "outside the indexed reference FASTA bounds")
})

test_that("expression fold-change helper uses a finite pseudocount by default", {
  default <- MMAPPR2:::.expressionLog2FC(c(20, 10, 0, 0), c(10, 0, 10, 0))
  expect_true(all(is.finite(default)))
  expect_equal(default[4], 0)

  raw <- MMAPPR2:::.expressionLog2FC(c(20, 10, 0, 0), c(10, 0, 10, 0), 0)
  expect_equal(raw[1], 1)
  expect_true(is.infinite(raw[2]) && raw[2] > 0)
  expect_true(is.infinite(raw[3]) && raw[3] < 0)
  expect_true(is.nan(raw[4]))

  opt <- MMAPPR2:::.expressionLog2FC(c(20, 10, 0, 0), c(10, 0, 10, 0), 0.5)
  expect_true(all(is.finite(opt)))
  expect_equal(opt[4], 0)
})

test_that("candidate pooling plan preserves fast path on roomy systems and chunks constrained ones", {
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 20000000))

  roomy <- MMAPPR2:::.candidatePoolPlan(
    region, nBams = 2, mode = "auto", chunkSize = 0,
    availableBytes = 64 * 1024^3
  )
  expect_false(roomy$chunked)

  constrained <- MMAPPR2:::.candidatePoolPlan(
    region, nBams = 2, mode = "auto", chunkSize = 0,
    availableBytes = 512 * 1024^2
  )
  expect_true(constrained$chunked)
  expect_lt(constrained$chunkSize, roomy$chunkSize)
  # Low-memory planning reserves extra space for dense multi-nucleotide pileup
  # rows, while >=8 GiB retains the higher-throughput RNA-seq-sparse estimate.
  expect_lte(constrained$chunkSize, 250000L)
  expect_gte(roomy$chunkSize, 80000000L)

  forced <- MMAPPR2:::.candidatePoolPlan(
    region, nBams = 2, mode = "chunked", chunkSize = 1000000,
    availableBytes = 64 * 1024^3
  )
  expect_true(forced$chunked)
  expect_equal(forced$chunkSize, 1000000L)
})



test_that("candidate auto chunk size stays bounded when memory is nearly exhausted", {
  expect_equal(
    MMAPPR2:::.candidateAutoChunkSize(nBams = 2, availableBytes = 1),
    10000L
  )
})


test_that("cgroup memory helper uses remaining bytes and ignores unlimited sentinels", {
  td <- tempfile("mmappr2_cgroup_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("1000", file.path(td, "memory.max"))
  writeLines("250", file.path(td, "memory.current"))
  expect_equal(
    MMAPPR2:::.cgroupMemoryRemaining(td, "memory.max", "memory.current", "max"),
    750
  )

  writeLines("max", file.path(td, "memory.max"))
  expect_length(
    MMAPPR2:::.cgroupMemoryRemaining(td, "memory.max", "memory.current", "max"),
    0L
  )

  writeLines(format(2^60, scientific = FALSE, trim = TRUE),
             file.path(td, "memory.limit_in_bytes"))
  writeLines("1", file.path(td, "memory.usage_in_bytes"))
  expect_length(
    MMAPPR2:::.cgroupMemoryRemaining(
      td, "memory.limit_in_bytes", "memory.usage_in_bytes"
    ),
    0L
  )
})

test_that("candidate range splitting is non-overlapping and exhaustive", {
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10, 34))
  chunks <- MMAPPR2:::.splitCandidateRanges(region, 10)
  expect_length(chunks, 3)
  got <- do.call(c, chunks)
  expect_equal(sum(BiocGenerics::width(got)), 25)
  expect_true(GenomicRanges::isDisjoint(got))
  expect_equal(BiocGenerics::start(got), c(10L, 20L, 30L))
  expect_equal(BiocGenerics::end(got), c(19L, 29L, 34L))
})


test_that("pooled base counts are unchanged by genomic chunking", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE, force = TRUE), add = TRUE)

  fa <- Rsamtools::FaFile(fx$fasta)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 220))
  whole <- MMAPPR2:::.pooledBasePileup(
    bams = c(fx$bam, fx$bam), genome = fa, which = region,
    minBaseQuality = 0L, minMapQuality = 0L, maxDepth = 250L
  )

  chunks <- MMAPPR2:::.splitCandidateRanges(region, 37L)
  chunked <- data.table::rbindlist(lapply(chunks, function(chunk) {
    MMAPPR2:::.pooledBasePileup(
      bams = c(fx$bam, fx$bam), genome = fa, which = chunk,
      minBaseQuality = 0L, minMapQuality = 0L, maxDepth = 250L
    )
  }))

  data.table::setorder(whole, seqnames, pos)
  data.table::setorder(chunked, seqnames, pos)
  # dcast() keys the one-shot table while rbindlist() intentionally drops that
  # optimization attribute. Candidate evidence must be identical regardless of
  # that data.table implementation detail.
  data.table::setkey(whole, NULL)
  data.table::setkey(chunked, NULL)
  expect_equal(chunked, whole)
})

.make_mock_candidate_vrange <- function(inputRange) {
  VariantAnnotation::VRanges(
    seqnames = S4Vectors::Rle(as.character(GenomicRanges::seqnames(inputRange))[1]),
    ranges = IRanges::IRanges(BiocGenerics::start(inputRange)[1], width = 1L),
    ref = "A", alt = "C",
    refDepth = S4Vectors::Rle(0L),
    altDepth = S4Vectors::Rle(2L),
    totalDepth = S4Vectors::Rle(2L),
    sampleNames = S4Vectors::Rle("pooled_mutant")
  )
}


test_that("candidate chunks recursively shrink after an allocation failure", {
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
  testthat::local_mocked_bindings(
    .candidateSnvsOneRange = function(inputRange, param, poolStrategy = NULL) {
      if (sum(BiocGenerics::width(inputRange)) > 25L)
        stop("cannot allocate vector of size 1.0 Gb")
      .make_mock_candidate_vrange(inputRange)
    },
    .package = "MMAPPR2"
  )

  got <- MMAPPR2:::.candidateSnvsChunkSafe(region, param = NULL)
  expect_length(got, 4L)
  expect_equal(BiocGenerics::start(got), c(1L, 26L, 51L, 76L))
})


test_that("emergency candidate disk staging is cleaned up", {
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
  pattern <- "^mmappr2_candidate_chunk_.*\\.rds$"
  before <- list.files(tempdir(), pattern = pattern, full.names = TRUE)

  testthat::local_mocked_bindings(
    .candidateSnvsChunkSafe = function(inputRange, param) {
      .make_mock_candidate_vrange(inputRange)
    },
    .package = "MMAPPR2"
  )

  got <- MMAPPR2:::.candidateSnvsChunked(
    region, param = NULL, chunkSize = 25L, spillToDisk = TRUE
  )
  after <- list.files(tempdir(), pattern = pattern, full.names = TRUE)

  expect_length(got, 4L)
  expect_equal(BiocGenerics::start(got), c(1L, 26L, 51L, 76L))
  expect_setequal(after, before)
})

test_that("candidate range splitting groups sparse loci without one chunk per site", {
  region <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(seq(1, 901, by = 100), width = 1L)
  )
  chunks <- MMAPPR2:::.splitCandidateRanges(region, 4L)
  expect_length(chunks, 3L)
  expect_equal(vapply(chunks, function(x) sum(BiocGenerics::width(x)), numeric(1)),
               c(4, 4, 2))
  expect_equal(sum(vapply(chunks, length, integer(1))), length(region))
})


test_that("incremental BAM pooling matches the fast in-memory strategy", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE, force = TRUE), add = TRUE)

  fa <- Rsamtools::FaFile(fx$fasta)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 220))
  fast <- MMAPPR2:::.pooledBasePileup(
    bams = c(fx$bam, fx$bam), genome = fa, which = region,
    minBaseQuality = 0L, minMapQuality = 0L, maxDepth = 250L,
    poolStrategy = "memory"
  )
  bounded <- MMAPPR2:::.pooledBasePileup(
    bams = c(fx$bam, fx$bam), genome = fa, which = region,
    minBaseQuality = 0L, minMapQuality = 0L, maxDepth = 250L,
    poolStrategy = "incremental"
  )

  data.table::setkey(fast, NULL)
  data.table::setkey(bounded, NULL)
  expect_equal(bounded, fast)
})


test_that("WT support is unchanged by forced bounded-memory chunking", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE, force = TRUE), add = TRUE)
  p <- suppressWarnings(mmapprParam(wtFiles = c(fx$bam, fx$bam), mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = file.path(fx$dir, "out_wt_chunks"),
                   includeScaffolds = TRUE, minDepth = 1,
                   minBaseQuality = 0, minMapQuality = 0,
                   candidatePoolMode = "memory"))

  vars <- VariantAnnotation::VRanges(
    seqnames = S4Vectors::Rle(rep("chr1", 3)),
    ranges = IRanges::IRanges(c(10L, 50L, 120L), width = 1L),
    ref = rep("A", 3), alt = rep("C", 3),
    refDepth = S4Vectors::Rle(rep(1L, 3)),
    altDepth = S4Vectors::Rle(rep(9L, 3)),
    totalDepth = S4Vectors::Rle(rep(10L, 3)),
    sampleNames = S4Vectors::Rle("pooled_mutant", 3)
  )

  fast <- MMAPPR2:::.addWtAlleleSupport(vars, p)
  candidatePoolMode(p) <- "chunked"
  candidateChunkSize(p) <- 2L
  bounded <- MMAPPR2:::.addWtAlleleSupport(vars, p)

  fields <- c("wtRefDepth", "wtAltDepth", "wtTotalDepth", "wtAltFreq")
  for (field in fields) {
    expect_equal(S4Vectors::mcols(bounded)[[field]],
                 S4Vectors::mcols(fast)[[field]])
  }
})
