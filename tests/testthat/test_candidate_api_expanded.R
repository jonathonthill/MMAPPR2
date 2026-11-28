test_that("mutant candidate filtering rejects malformed depths and honors inclusive depth boundaries", {
  keep <- MMAPPR2:::.candidateMutantKeep(
    totalDepth = c(10, 10, 10, -1, 5, NA, Inf),
    altDepth = c(9, 8, 10, 1, 6, 9, 9),
    minDepth = 10, minAltDepth = 8, minAltFreq = .8
  )
  expect_identical(keep, c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("candidate memory error classifier is specific to allocation failures", {
  f <- MMAPPR2:::.isMemoryAllocationError
  expect_true(f(simpleError("cannot allocate vector of size 2 Gb")))
  expect_true(f(simpleError("std::bad_alloc")))
  expect_true(f(simpleError("vector memory exhausted")))
  expect_false(f(simpleError("BAM index is missing")))
})

test_that("candidate auto chunk planning has documented memory tiers and BAM scaling", {
  gib <- 1024^3
  f <- MMAPPR2:::.candidateAutoChunkSize
  expect_identical(f(1, NA_real_), 10000000L)
  expect_gte(f(1, 1), 10000L)
  expect_lte(f(1, 100 * gib), 250000000L)
  expect_gt(f(1, 4 * gib), f(2, 4 * gib))
  expect_gt(f(1, 8 * gib), f(1, 2 * gib))
})

test_that("candidate pool plan honors forced modes, explicit chunk size, and empty ranges", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000))
  expect_false(MMAPPR2:::.candidatePoolPlan(gr, 1, "memory", 10, 1)$chunked)
  expect_true(MMAPPR2:::.candidatePoolPlan(gr, 1, "chunked", 10000, 1e12)$chunked)
  expect_true(MMAPPR2:::.candidatePoolPlan(gr, 1, "auto", 100, 1e12)$chunked)
  expect_false(MMAPPR2:::.candidatePoolPlan(gr, 1, "auto", 1000, 1)$chunked)
  empty <- GenomicRanges::GRanges()
  plan <- MMAPPR2:::.candidatePoolPlan(empty, 1, "auto", 0, 1e9)
  expect_false(plan$chunked)
  expect_equal(plan$totalWidth, 0)
})

test_that("candidate range splitting handles empty, multiple chromosomes, and preserves queried bases", {
  expect_identical(MMAPPR2:::.splitCandidateRanges(GenomicRanges::GRanges(), 10), list())
  gr <- GenomicRanges::GRanges(c("chr1", "chr1", "chr2"),
                               IRanges::IRanges(c(1, 100, 1), c(25, 104, 7)))
  chunks <- MMAPPR2:::.splitCandidateRanges(gr, 10)
  expect_true(all(vapply(chunks, function(x) sum(BiocGenerics::width(x)) <= 10, logical(1))))
  expect_equal(sum(vapply(chunks, function(x) sum(BiocGenerics::width(x)), numeric(1))),
               sum(BiocGenerics::width(gr)))
})

test_that("base pileup conversion handles reference-only, multiallelic, and non-ACGT reference rows", {
  pile <- data.table::data.table(
    seqnames = c("chr1", "chr1", "chr1"), pos = 1:3,
    A = c(10L, 5L, 1L), C = c(0L, 3L, 2L), G = c(0L, 2L, 3L), T = 0L,
    totalDepth = c(10L, 10L, 6L), ref = c("A", "A", "N")
  )
  vr <- MMAPPR2:::.basePileupToVRanges(pile, sampleName = "x")
  expect_equal(length(vr), 2)
  expect_true(all(BiocGenerics::start(vr) == 2L))
  expect_setequal(as.character(VariantAnnotation::alt(vr)), c("C", "G"))
  expect_true(all(as.character(VariantAnnotation::sampleNames(vr)) == "x"))
})

test_that("candidate SNV conversion returns only variants passing mutant thresholds", {
  p <- .signal_param(candidateMinDepth = 10, candidateMinAltDepth = 10, candidateMinAltFreq = .8)
  pile <- data.table::data.table(
    seqnames = c("chr1", "chr1"), pos = c(1L, 2L),
    A = c(1L, 4L), C = c(19L, 16L), G = 0L, T = 0L,
    totalDepth = c(20L, 20L), ref = c("A", "A")
  )
  vr <- MMAPPR2:::.candidateSnvsFromPileup(pile, p)
  expect_equal(length(vr), 1)
  expect_equal(BiocGenerics::start(vr), 1L)
})

test_that("candidate chunk recursion propagates non-memory errors and unsplittable memory failures", {
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  testthat::local_mocked_bindings(
    .candidateSnvsOneRange = function(...) stop("bad BAM"),
    .package = "MMAPPR2"
  )
  expect_error(MMAPPR2:::.candidateSnvsChunkSafe(region, NULL), "bad BAM")

  one <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, width = 1))
  testthat::local_mocked_bindings(
    .candidateSnvsOneRange = function(...) stop("cannot allocate vector"),
    .package = "MMAPPR2"
  )
  expect_error(MMAPPR2:::.candidateSnvsChunkSafe(one, NULL), "cannot allocate")
})

test_that("candidate VRanges combiner handles empty, singleton, and multiple pieces", {
  empty <- MMAPPR2:::.combineCandidateVranges(list(VariantAnnotation::VRanges()))
  expect_equal(length(empty), 0)
  a <- VariantAnnotation::VRanges(seqnames = "chr1", ranges = IRanges::IRanges(1, width = 1),
                                  ref = "A", alt = "C", refDepth = 1L, altDepth = 9L,
                                  totalDepth = 10L, sampleNames = "x")
  b <- a; BiocGenerics::start(b) <- 2L
  expect_equal(length(MMAPPR2:::.combineCandidateVranges(list(a))), 1)
  expect_equal(length(MMAPPR2:::.combineCandidateVranges(list(a, b))), 2)
})

test_that("WT depth extraction handles empty pileups, absent coordinates, and alternative alleles", {
  vars <- VariantAnnotation::VRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr1")), ranges = IRanges::IRanges(c(10, 20), width = 1),
    ref = c("A", "A"), alt = c("C", "G"), refDepth = c(1L, 1L),
    altDepth = c(9L, 9L), totalDepth = c(10L, 10L), sampleNames = S4Vectors::Rle("m", 2)
  )
  empty <- MMAPPR2:::.wtDepthsFromPile(vars, MMAPPR2:::.emptyPooledBasePileup())
  expect_identical(empty$total, c(0L, 0L))
  pile <- data.table::data.table(seqnames = "chr1", pos = 10L, A = 7L, C = 3L, G = 0L, T = 0L,
                                 totalDepth = 10L, ref = "A")
  got <- MMAPPR2:::.wtDepthsFromPile(vars, pile)
  expect_identical(got$ref, c(7L, 0L))
  expect_identical(got$alt, c(3L, 0L))
  expect_identical(got$total, c(10L, 0L))
})

test_that("default candidate discovery retains a strong mutant SNV with zero WT coverage", {
  fx <- .make_signal_fixture(wt_zero_at_candidate = TRUE)
  p <- .signal_param(wt_zero_at_candidate = TRUE, candidatePoolMode = "memory")
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(450, 550))
  vr <- MMAPPR2:::.getVariantsForRange(region, p)
  expect_equal(length(vr), 1)
  expect_equal(BiocGenerics::start(vr), fx$candidate_pos)
  expect_equal(S4Vectors::mcols(vr)$wtTotalDepth, 0)
  expect_true(is.na(S4Vectors::mcols(vr)$wtAltFreq))
  expect_true(as.numeric(S4Vectors::mcols(vr)$mutAltFreq) > .8)
})

test_that("explicit WT filtering removes candidates lacking WT evidence while defaults do not", {
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(450, 550))
  permissive <- .signal_param(wt_zero_at_candidate = TRUE, candidateMaxWtAltFreq = 1)
  strict <- .signal_param(wt_zero_at_candidate = TRUE, candidateMaxWtAltFreq = .5)
  expect_equal(length(MMAPPR2:::.getVariantsForRange(region, permissive)), 1)
  expect_equal(length(MMAPPR2:::.getVariantsForRange(region, strict)), 0)
})

test_that("WT support metadata flags high-frequency shared background without default exclusion", {
  p <- .signal_param()
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(450, 550))
  vr <- MMAPPR2:::.getVariantsForRange(region, p)
  expect_equal(length(vr), 1)
  expect_true(S4Vectors::mcols(vr)$wtEvidenceSufficient)
  expect_false(S4Vectors::mcols(vr)$sharedBackgroundLike)
  expect_equal(as.numeric(S4Vectors::mcols(vr)$wtAltFreq), .4, tolerance = .001)
  expect_equal(as.numeric(S4Vectors::mcols(vr)$deltaAltFreq), .55, tolerance = .001)
})

test_that("known internal range warning matcher is call-sensitive", {
  w <- simpleWarning("GRanges object contains 1 out-of-bound ranges located on sequence chr1",
                     call = quote(valid.GenomicRanges.seqinfo(x)))
  expect_true(MMAPPR2:::.isKnownInternalRangeWarning(w))
  w2 <- simpleWarning("GRanges object contains 1 out-of-bound ranges located on sequence chr1",
                      call = quote(other(x)))
  expect_false(MMAPPR2:::.isKnownInternalRangeWarning(w2))
  expect_false(MMAPPR2:::.isKnownInternalRangeWarning(simpleWarning("other")))
})

test_that("annotation attribute extractors handle present, absent, quoted, and equals forms", {
  x <- c('gene_id "g1"; gene_name "Name 1";', "ID=g2;Name=Name2", "foo=bar")
  expect_identical(MMAPPR2:::.extractQuotedAnnotationAttribute(x, "gene_id"), c("g1", NA, NA))
  expect_identical(MMAPPR2:::.extractEqualsAnnotationAttribute(x, "ID"), c(NA, "g2", NA))
  id <- MMAPPR2:::.annotationGeneIdentity(x)
  expect_identical(id$gene_id, c("g1", "g2", NA))
  expect_identical(id$gene_name, c("Name 1", "Name2", NA))
})

test_that("annotation gene ranges fall back from missing gene names to IDs", {
  fx <- .test_fixture()
  p <- .test_param(fx)
  genes <- MMAPPR2:::.annotationGeneRanges(p)
  expect_true(length(genes) >= 1)
  expect_true(all(c("gene_id", "gene_name") %in% colnames(S4Vectors::mcols(genes))))
  expect_true(all(!is.na(S4Vectors::mcols(genes)$gene_name)))
})

test_that("coding prediction propagates candidate evidence by QUERYID", {
  fx <- .make_signal_fixture()
  p <- .signal_param()
  vr <- MMAPPR2:::.getVariantsForRange(GenomicRanges::GRanges("chr1", IRanges::IRanges(490, 510)), p)
  tx <- MMAPPR2:::.buildTxDb(p)
  on.exit(MMAPPR2:::.disconnectTxDb(tx), add = TRUE)
  eff <- MMAPPR2:::.predictEffects(vr, p, tx)
  expect_true(length(eff) >= 1)
  for (nm in c("mutAltFreq", "wtAltFreq", "deltaAltFreq", "wtEvidenceSufficient")) {
    expect_true(nm %in% colnames(S4Vectors::mcols(eff)))
  }
  expect_true(all(is.finite(S4Vectors::mcols(eff)$mutAltFreq)))
})

test_that("coding prediction rejects candidate coordinates outside the reference", {
  p <- .signal_param()
  tx <- MMAPPR2:::.buildTxDb(p)
  on.exit(MMAPPR2:::.disconnectTxDb(tx), add = TRUE)
  vr <- VariantAnnotation::VRanges(seqnames = "chr1", ranges = IRanges::IRanges(1300, width = 1),
                                  ref = "A", alt = "C", refDepth = 1L, altDepth = 9L,
                                  totalDepth = 10L, sampleNames = "x")
  expect_error(MMAPPR2:::.predictEffects(vr, p, tx), "outside the indexed reference")
})

test_that("descriptive expression returns empty results when no peak gene passes fold-change criteria", {
  p <- .signal_param()
  tx <- MMAPPR2:::.buildTxDb(p)
  on.exit(MMAPPR2:::.disconnectTxDb(tx), add = TRUE)
  peak <- GenomicRanges::GRanges("chr1", IRanges::IRanges(300, 700))
  genes <- MMAPPR2:::.annotationGeneRanges(p, tx)
  out <- MMAPPR2:::.addDiff(peak, p, tx, annotationGenes = genes)
  expect_s4_class(out, "GRanges")
  expect_equal(length(out), 0)
})

test_that("candidate scoring uses feature midpoint density and coding severity before density", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 20, 30, 40), width = 3))
  S4Vectors::mcols(gr)$CONSEQUENCE <- c("synonymous", "nonsynonymous", "frameshift", "nonsense")
  cand <- list(effects = list(chr1 = gr))
  peaks <- list(chr1 = list(densityFunction = function(x) x / 100))
  out <- MMAPPR2:::.scoreVariants(cand, peaks)
  expect_equal(S4Vectors::mcols(out$effects$chr1)$CONSEQUENCE[1], "nonsense")
  expect_equal(S4Vectors::mcols(out$effects$chr1)$peakDensity[1], 41 / 100)
})

test_that("candidate ordering warns on invalid classes and places missing density last", {
  expect_warning(out <- MMAPPR2:::.orderVariants(data.frame(x = 1)), "Invalid data type")
  expect_s3_class(out, "data.frame")
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:3, width = 1))
  S4Vectors::mcols(gr)$peakDensity <- c(.2, NA, .5)
  ordered <- MMAPPR2:::.orderVariants(gr)
  expect_equal(BiocGenerics::start(ordered), c(3L, 1L, 2L))
})

test_that("generateCandidates directly populates SNV, effect, and expression groups", {
  fx <- .make_signal_fixture()
  p <- .signal_param()
  md <- mmapprData(p)
  md@peaks <- list(chr1 = .make_peak_entry("chr1", 300, 700, fx$candidate_pos))
  out <- generateCandidates(md)
  expect_setequal(names(candidates(out)), c("snps", "effects", "diff"))
  expect_equal(length(candidates(out)$snps$chr1), 1)
  expect_equal(BiocGenerics::start(candidates(out)$snps$chr1), fx$candidate_pos)
  expect_true("peakDensity" %in% colnames(S4Vectors::mcols(candidates(out)$snps$chr1)))
  expect_true(length(candidates(out)$effects$chr1) >= 1)
})

test_that("generateCandidates handles a refined interval with no passing SNVs", {
  p <- .signal_param()
  md <- mmapprData(p)
  md@peaks <- list(chr1 = .make_peak_entry("chr1", 100, 150, 125))
  out <- generateCandidates(md)
  expect_equal(length(candidates(out)$snps$chr1), 0)
  expect_equal(length(candidates(out)$effects$chr1), 0)
  expect_true("diff" %in% names(candidates(out)))
})

test_that("generateCandidates refuses to run without refined peaks", {
  expect_error(generateCandidates(mmapprData(.signal_param())), "No refined peaks")
})

test_that("pooled base pileup keeps all nucleotide count columns integer and converts without type warnings", {
  p <- .signal_param(candidatePoolMode = "memory")
  pile <- MMAPPR2:::.pooledBasePileup(
    mutFiles(p), p@refGenome,
    which = GenomicRanges::GRanges("chr1", IRanges::IRanges(490, 510)),
    minBaseQuality = 0, minMapQuality = 0, maxDepth = 250
  )
  for (base in c("A", "C", "G", "T")) expect_type(pile[[base]], "integer")
  expect_silent(MMAPPR2:::.basePileupToVRanges(pile))
})


test_that("whole candidate filtering keeps background-like WT evidence by default and honors opt-in WT/delta filters", {
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(490, 510))
  permissive <- .signal_param(candidateMaxWtAltFreq = 1, candidateMinDeltaAF = 0)
  wtStrict <- .signal_param(candidateMaxWtAltFreq = .90, candidateMinDeltaAF = 0)
  deltaStrict <- .signal_param(candidateMaxWtAltFreq = 1, candidateMinDeltaAF = .10)

  highWtSupport <- function(variants, param, which = NULL) {
    S4Vectors::mcols(variants)$wtRefDepth <- rep(5L, length(variants))
    S4Vectors::mcols(variants)$wtAltDepth <- rep(95L, length(variants))
    S4Vectors::mcols(variants)$wtTotalDepth <- rep(100L, length(variants))
    S4Vectors::mcols(variants)$wtAltFreq <- rep(.95, length(variants))
    variants
  }
  testthat::local_mocked_bindings(
    .addWtAlleleSupport = highWtSupport,
    .package = "MMAPPR2"
  )

  kept <- MMAPPR2:::.getVariantsForRange(region, permissive)
  expect_equal(length(kept), 1)
  expect_true(S4Vectors::mcols(kept)$sharedBackgroundLike)
  expect_equal(as.numeric(S4Vectors::mcols(kept)$deltaAltFreq), 0, tolerance = .001)

  expect_equal(length(MMAPPR2:::.getVariantsForRange(region, wtStrict)), 0)
  expect_equal(length(MMAPPR2:::.getVariantsForRange(region, deltaStrict)), 0)
})

test_that("pooled base pileup accepts supported BAM container forms and optional whole-file queries", {
  fx <- .make_quality_fixture()
  genome <- Rsamtools::FaFile(fx$fasta)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 30))
  bf <- Rsamtools::BamFile(fx$bam)
  bfl <- Rsamtools::BamFileList(bf)

  fromChar <- MMAPPR2:::.pooledBasePileup(fx$bam, genome, which = region,
                                          minBaseQuality = 0, minMapQuality = 0,
                                          maxDepth = 100)
  fromBam <- MMAPPR2:::.pooledBasePileup(bf, genome, which = region,
                                         minBaseQuality = 0, minMapQuality = 0,
                                         maxDepth = 100)
  fromList <- MMAPPR2:::.pooledBasePileup(bfl, genome, which = region,
                                          minBaseQuality = 0, minMapQuality = 0,
                                          maxDepth = 100)
  expect_equal(fromBam, fromChar)
  expect_equal(fromList, fromChar)

  whole <- MMAPPR2:::.pooledBasePileup(fx$bam, genome,
                                       minBaseQuality = 0, minMapQuality = 0,
                                       maxDepth = 100)
  expect_setequal(whole$pos, c(10L, 20L))
  expect_error(MMAPPR2:::.pooledBasePileup(character(), genome),
               "at least one BAM")
  expect_error(MMAPPR2:::.pooledBasePileup(list(fx$bam), genome),
               "at least one BAM")
})

test_that("annotation gene ranges fall back to TxDb genes when explicit gene rows are absent", {
  fx <- .make_signal_fixture()
  noGene <- file.path(fx$dir, "no-explicit-gene.gtf")
  writeLines(c(
    'chr1\tMMAPPR2\ttranscript\t1\t1100\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";',
    'chr1\tMMAPPR2\texon\t1\t1100\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; exon_id "e1";',
    'chr1\tMMAPPR2\tCDS\t2\t1099\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";'
  ), noGene)
  p <- .signal_param(gtf = noGene)
  genes <- MMAPPR2:::.annotationGeneRanges(p)
  expect_equal(length(genes), 1)
  expect_identical(as.character(S4Vectors::mcols(genes)$gene_id), "g1")
  expect_identical(as.character(S4Vectors::mcols(genes)$gene_name), "g1")
  expect_equal(BiocGenerics::start(genes), 1L)
  expect_equal(BiocGenerics::end(genes), 1100L)
})
