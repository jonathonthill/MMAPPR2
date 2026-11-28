test_that("mapping pileup frequencies sum to one and no-read regions return the stable empty schema", {
  fx <- .test_fixture()
  p <- .test_param(fx, includeScaffolds = TRUE)
  covered <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10, 29))
  x <- MMAPPR2:::.getPileup(wtFiles(p)[[1]], p, covered)
  expect_gt(nrow(x), 0)
  expect_equal(x$A.FREQ + x$C.FREQ + x$G.FREQ + x$T.FREQ, rep(1, nrow(x)))
  expect_true(all(x$CVG > 0))

  empty <- MMAPPR2:::.getPileup(wtFiles(p)[[1]], p,
                                GenomicRanges::GRanges("chr1", IRanges::IRanges(220, 240)))
  expect_identical(names(empty), names(MMAPPR2:::.emptyPileupTable()))
  expect_equal(nrow(empty), 0)
})

test_that("WT homozygous filtering keeps the exact cutoff boundary and rejects unusable rows", {
  x <- data.table::data.table(
    A = c(.95, .951, NA), C = c(.05, .049, NA), G = c(0, 0, NA), T = c(0, 0, NA)
  )
  cols <- names(x)
  got <- MMAPPR2:::.wtHomozygousRows(x, cols, .95)
  expect_identical(got, c(FALSE, TRUE, TRUE))
})

test_that("replicate aggregation validates mode and file count", {
  x <- data.table::data.table(CHROM = "chr1", POS = 1L, FILE_ID = 1L,
                              A.FREQ = 1, C.FREQ = 0, G.FREQ = 0, T.FREQ = 0, CVG = 10)
  expect_error(MMAPPR2:::.avgFiles(data.table::copy(x), "bad", 1), "fileAggregation")
  expect_error(MMAPPR2:::.avgFiles(data.table::copy(x), "simple", 0), "nFiles")
})

test_that("per-chromosome distance uses the Euclidean frequency formula and distance power", {
  p <- .signal_param(minDepth = 1, distancePower = 2)
  wtPath <- normalizePath(.make_signal_fixture()$wt)
  mutPath <- normalizePath(.make_signal_fixture()$mut)
  wt <- data.table::data.table(CHROM = "chr1", POS = 100L,
                              A.FREQ = .6, C.FREQ = .4, G.FREQ = 0, T.FREQ = 0, CVG = 20)
  mt <- data.table::data.table(CHROM = "chr1", POS = 100L,
                              A.FREQ = .2, C.FREQ = .8, G.FREQ = 0, T.FREQ = 0, CVG = 20)
  testthat::local_mocked_bindings(
    .getPileup = function(file, param, chrRange) {
      path <- normalizePath(BiocGenerics::path(file))
      if (identical(path, wtPath)) data.table::copy(wt) else if (identical(path, mutPath)) data.table::copy(mt) else stop("unexpected BAM")
    },
    .package = "MMAPPR2"
  )
  ans <- MMAPPR2:::.calcDistForChr(GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000)), p)
  expect_true(is.list(ans))
  expected <- (sqrt((.6 - .2)^2 + (.4 - .8)^2))^2
  expect_equal(ans$distanceDf$DISTANCE, expected)
})

test_that("per-chromosome distance isolates common failure modes as informative strings", {
  p <- .signal_param(minDepth = 10)
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000))
  empty <- MMAPPR2:::.emptyPileupTable()
  testthat::local_mocked_bindings(
    .getPileup = function(file, param, chrRange) empty,
    .package = "MMAPPR2"
  )
  ans <- MMAPPR2:::.calcDistForChr(region, p)
  expect_type(ans, "character")
  expect_match(ans, "Insufficient data in wild-type")
})

test_that("per-chromosome distance reports an empty WT-mutant coordinate join", {
  p <- .signal_param(minDepth = 1)
  wtPath <- normalizePath(.make_signal_fixture()$wt)
  wt <- data.table::data.table(CHROM = "chr1", POS = 100L,
                              A.FREQ = .6, C.FREQ = .4, G.FREQ = 0, T.FREQ = 0, CVG = 20)
  mt <- data.table::data.table(CHROM = "chr1", POS = 200L,
                              A.FREQ = .6, C.FREQ = .4, G.FREQ = 0, T.FREQ = 0, CVG = 20)
  testthat::local_mocked_bindings(
    .getPileup = function(file, param, chrRange) {
      if (identical(normalizePath(BiocGenerics::path(file)), wtPath)) data.table::copy(wt) else data.table::copy(mt)
    }, .package = "MMAPPR2"
  )
  expect_match(MMAPPR2:::.calcDistForChr(GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000)), p),
               "Empty dataframe after joining")
})

test_that("calculateDistance directly populates the public MmapprData distance slot", {
  p <- .signal_param()
  md <- mmapprData(p)
  out <- calculateDistance(md)
  expect_s4_class(out, "MmapprData")
  expect_named(snpDistance(out), "chr1")
  expect_true(is.list(snpDistance(out)$chr1))
  expect_true(all(c("wtCounts", "mutCounts", "distanceDf") %in% names(snpDistance(out)$chr1)))
  expect_gte(nrow(snpDistance(out)$chr1$distanceDf), 50)
  expect_true(all(is.finite(snpDistance(out)$chr1$distanceDf$DISTANCE)))
})

test_that("calculateDistance rejects an annotation stage with no query ranges", {
  p <- .test_param()
  md <- mmapprData(p)
  testthat::local_mocked_bindings(.getFileReadChrList = function(param) list(), .package = "MMAPPR2")
  expect_error(calculateDistance(md), "No annotated genomic ranges")
})

test_that("mapping query ranges apply scaffold and mitochondrial policies explicitly", {
  # Construct real parameter objects before mocking reference metadata; otherwise
  # constructor preflight would (correctly) compare the fixture BAM against the mock.
  pStandard <- .test_param(includeScaffolds = FALSE)
  pScaffolds <- .test_param(includeScaffolds = TRUE)
  genes <- GenomicRanges::GRanges(
    c("chr1", "chrM", "chrUn_gl0001"),
    IRanges::IRanges(c(1, 1, 1), c(50, 50, 50))
  )
  refSi <- GenomeInfoDb::Seqinfo(c("chr1", "chrM", "chrUn_gl0001"),
                                 seqlengths = c(100, 100, 100))
  testthat::local_mocked_bindings(
    .annotationGeneRanges = function(param) genes,
    .faSeqinfo = function(fa) refSi,
    .choose_target_style = function(si) NA_character_,
    .package = "MMAPPR2"
  )

  standard <- MMAPPR2:::.getFileReadChrList(pStandard)
  expect_identical(names(standard), "chr1")

  withScaffolds <- MMAPPR2:::.getFileReadChrList(pScaffolds)
  expect_setequal(names(withScaffolds), c("chr1", "chrUn_gl0001"))
  expect_false("chrM" %in% names(withScaffolds))
})

test_that("mapping query range construction reports empty, incompatible, custom-only, and mitochondrial-only annotations", {
  # Build all real parameter objects before installing test-local mocks because
  # mocked FASTA metadata would otherwise be used by constructor preflight.
  p <- .test_param(includeScaffolds = TRUE)
  pNoScaffolds <- .test_param(includeScaffolds = FALSE)
  pMito <- .test_param(includeScaffolds = TRUE)

  testthat::local_mocked_bindings(
    .annotationGeneRanges = function(param) GenomicRanges::GRanges(),
    .package = "MMAPPR2"
  )
  expect_error(MMAPPR2:::.getFileReadChrList(p), "No gene ranges")

  testthat::local_mocked_bindings(
    .annotationGeneRanges = function(param) GenomicRanges::GRanges("other", IRanges::IRanges(1, 10)),
    .faSeqinfo = function(fa) GenomeInfoDb::Seqinfo("chr1", seqlengths = 100),
    .choose_target_style = function(si) NA_character_,
    .package = "MMAPPR2"
  )
  expect_error(MMAPPR2:::.getFileReadChrList(p), "do not share any sequence names")

  testthat::local_mocked_bindings(
    .annotationGeneRanges = function(param) GenomicRanges::GRanges("chrUn_gl0001", IRanges::IRanges(1, 10)),
    .faSeqinfo = function(fa) GenomeInfoDb::Seqinfo("chrUn_gl0001", seqlengths = 100),
    .choose_target_style = function(si) NA_character_,
    .package = "MMAPPR2"
  )
  expect_error(MMAPPR2:::.getFileReadChrList(pNoScaffolds),
               "No standard chromosomes remained")

  testthat::local_mocked_bindings(
    .annotationGeneRanges = function(param) GenomicRanges::GRanges("chrM", IRanges::IRanges(1, 10)),
    .faSeqinfo = function(fa) GenomeInfoDb::Seqinfo("chrM", seqlengths = 100),
    .choose_target_style = function(si) NA_character_,
    .package = "MMAPPR2"
  )
  expect_error(MMAPPR2:::.getFileReadChrList(pMito),
               "No non-mitochondrial")
})

test_that("mapping pileup enforces mapping quality, base quality, and maximum depth", {
  fx <- .make_quality_fixture()
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 30))

  strict <- .test_param(fx, includeScaffolds = TRUE,
                        minBaseQuality = 20, minMapQuality = 20,
                        maxPileupDepth = 100)
  x <- MMAPPR2:::.getPileup(wtFiles(strict)[[1]], strict, region)
  p10 <- x[x$POS == 10L]
  expect_equal(p10$CVG, 1)
  expect_equal(p10$A.FREQ, 1)
  expect_equal(p10$C.FREQ, 0)
  expect_equal(p10$G.FREQ, 0)

  permissive <- .test_param(fx, includeScaffolds = TRUE,
                            minBaseQuality = 0, minMapQuality = 0,
                            maxPileupDepth = 100)
  y <- MMAPPR2:::.getPileup(wtFiles(permissive)[[1]], permissive, region)
  q10 <- y[y$POS == 10L]
  expect_equal(q10$CVG, 3)
  expect_equal(c(q10$A.FREQ, q10$C.FREQ, q10$G.FREQ), rep(1 / 3, 3))

  capped <- .test_param(fx, includeScaffolds = TRUE,
                        minBaseQuality = 0, minMapQuality = 0,
                        maxPileupDepth = 1)
  z <- MMAPPR2:::.getPileup(wtFiles(capped)[[1]], capped, region)
  expect_equal(z[POS == 20L]$CVG, 1)
})
