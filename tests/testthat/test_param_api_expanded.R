test_that("scalar numeric validation respects open, closed, finite, and integer bounds", {
  f <- MMAPPR2:::.scalarNumeric
  expect_true(f(0, lower = 0))
  expect_false(f(0, lower = 0, lowerOpen = TRUE))
  expect_true(f(1, upper = 1))
  expect_false(f(1, upper = 1, upperOpen = TRUE))
  expect_true(f(2, integerish = TRUE))
  expect_false(f(2.5, integerish = TRUE))
  for (x in list(numeric(), c(1, 2), NA_real_, NaN, Inf, "1")) {
    expect_false(f(x))
  }
})

test_that("all scalar parameter setters round-trip valid values and reject invalid values", {
  p <- .test_param()
  cases <- list(
    includeScaffolds = list(TRUE, NA),
    minDepth = list(2, 0),
    homozygoteCutoff = list(.8, 1.1),
    minBaseQuality = list(1, -1),
    minMapQuality = list(1, -1),
    fileAggregation = list("weighted", "bad"),
    distancePower = list(2, 0),
    peakIntervalWidth = list(.5, 0),
    loessOptResolution = list(.01, 0),
    loessOptCutFactor = list(.2, 1),
    maxPileupDepth = list(200, 0),
    candidateMinDepth = list(2, 0),
    candidateMinAltDepth = list(3, 0),
    candidateMinAltFreq = list(.7, 1),
    candidateMaxWtAltFreq = list(.8, 1.1),
    candidateMinDeltaAF = list(.1, -1),
    candidatePoolMode = list("chunked", "bad"),
    candidateChunkSize = list(100, -1),
    peakCutoffSd = list(2, -1),
    peakCutoffMethod = list("global_sd", "bad"),
    peakIntervalMethod = list("shortest_contiguous", "bad"),
    peakResampleIterations = list(20, 9),
    randomSeed = list(123, -1),
    pairedEnd = list(TRUE, NA),
    ignoreStrand = list(TRUE, NA),
    expressionPseudocount = list(.5, -1),
    exportAiccPlots = list(TRUE, NA)
  )

  for (nm in names(cases)) {
    getter <- get(nm, mode = "function")
    setter <- get(paste0(nm, "<-"), mode = "function")
    p2 <- setter(p, cases[[nm]][[1]])
    expect_identical(getter(p2), cases[[nm]][[1]], info = nm)
    expect_error(setter(p, cases[[nm]][[2]]), info = nm)
  }
})

test_that("depth capacity validation scales with the number of mutant BAMs", {
  base <- list(
    includeScaffolds = FALSE, minDepth = 1, homozygoteCutoff = .95,
    minBaseQuality = 0, minMapQuality = 0, fileAggregation = "simple",
    distancePower = 4, peakIntervalWidth = .8, loessOptResolution = .001,
    loessOptCutFactor = .1, maxPileupDepth = 1000,
    candidateMinDepth = 500, candidateMinAltDepth = 500,
    candidateMinAltFreq = .8, candidateMaxWtAltFreq = 1,
    candidateMinDeltaAF = 0, candidatePoolMode = "auto",
    candidateChunkSize = 0, peakCutoffSd = 3,
    peakCutoffMethod = "legacy_current", peakIntervalMethod = "hpd_span",
    peakResampleIterations = 10, randomSeed = 1, pairedEnd = FALSE,
    ignoreStrand = FALSE, expressionPseudocount = .01,
    exportAiccPlots = FALSE, nMutFiles = 2L
  )
  expect_length(do.call(MMAPPR2:::.validateScalarValues, base), 0)
  base$candidateMinDepth <- 501
  expect_match(paste(do.call(MMAPPR2:::.validateScalarValues, base), collapse = " "),
               "candidateMinDepth exceeds")
})

test_that("BAM input coercion accepts supported forms and rejects empty input", {
  fx <- .test_fixture()
  bf <- Rsamtools::BamFile(fx$bam)
  bfl <- Rsamtools::BamFileList(bf)
  expect_identical(MMAPPR2:::.asBamPaths(fx$bam), fx$bam)
  expect_identical(MMAPPR2:::.asBamPaths(bf), fx$bam)
  expect_identical(MMAPPR2:::.asBamPaths(bfl), fx$bam)
  expect_error(MMAPPR2:::.asBamPaths(character()), "BAM inputs must be character paths")
  expect_error(MMAPPR2:::.asBamPaths(1), "BAM")
})

test_that("mmapprData constructor and accessors expose empty initial state", {
  p <- .test_param()
  md <- mmapprData(p)
  expect_s4_class(md, "MmapprData")
  expect_identical(param(md), p)
  expect_identical(snpDistance(md), list())
  expect_identical(peaks(md), list())
  expect_identical(candidates(md), list())
  expect_error(mmapprData("not a parameter"), "MmapprParam")
  expect_error(`param<-`(md, "bad"), "MmapprParam")
})

test_that("replacing parameters before analysis does not warn or fabricate derived state", {
  p1 <- .test_param()
  p2 <- .test_param()
  md <- mmapprData(p1)
  expect_silent(md <- `param<-`(md, p2))
  expect_identical(param(md), p2)
  expect_length(snpDistance(md), 0)
  expect_length(peaks(md), 0)
  expect_length(candidates(md), 0)
})

test_that("MmapprParam and MmapprData show methods handle empty and populated state", {
  p <- .test_param()
  expect_output(show(p), "MmapprParam")
  md <- mmapprData(p)
  expect_output(show(md), "Euclidean distance data for 0")
  md@snpDistance <- list(chr1 = .make_loess_result())
  md@peaks <- list(chr1 = .make_peak_entry())
  md@candidates <- list(snps = list(chr1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(500, width = 1))))
  txt <- capture.output(show(md))
  expect_true(any(grepl("Loess regression data for 1", txt)))
  expect_true(any(grepl("start = 300, end = 700", txt)))
})

test_that("custom printer truncates long components without failing", {
  txt <- capture.output(MMAPPR2:::.customPrint(list(x = 1:20), lineMax = 3))
  expect_true(length(txt) > 0)
})

test_that("BAM seqinfo extraction returns the indexed sequence dictionary", {
  fx <- .test_fixture()
  si <- MMAPPR2:::.bamSeqinfoFrozen(Rsamtools::BamFile(fx$bam))
  expect_identical(GenomeInfoDb::seqlevels(si), "chr1")
  expect_equal(as.numeric(GenomeInfoDb::seqlengths(si)["chr1"]), 250)
})

test_that("BAM seqinfo comparison catches missing names and conflicting lengths", {
  ref <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), seqlengths = c(100, 200))
  expect_error(MMAPPR2:::.compareBamSeqinfoFrozen(
    "bam", GenomeInfoDb::Seqinfo("chr3", seqlengths = 100), ref),
    "shares no exact sequence names")
  expect_error(MMAPPR2:::.compareBamSeqinfoFrozen(
    "bam", GenomeInfoDb::Seqinfo("chr1", seqlengths = 101), ref),
    "sequence-length conflict")
  expect_identical(MMAPPR2:::.compareBamSeqinfoFrozen(
    "bam", GenomeInfoDb::Seqinfo("chr1", seqlengths = 100), ref), "chr1")
})

test_that("annotation preflight rejects malformed coordinates and accepts the reference boundary", {
  td <- tempfile("ann_bounds_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  nonnumeric <- file.path(td, "nonnumeric.gtf")
  writeLines('chr1\tx\tgene\tX\t10\t.\t+\t.\tgene_id "g";', nonnumeric)
  expect_error(MMAPPR2:::.annotationObservedRanges(nonnumeric), "non-numeric")
  invalid <- file.path(td, "invalid.gtf")
  writeLines('chr1\tx\tgene\t10\t9\t.\t+\t.\tgene_id "g";', invalid)
  expect_error(MMAPPR2:::.annotationObservedRanges(invalid), "invalid genomic coordinate")

  ok <- file.path(td, "ok.gtf")
  writeLines('chr1\tx\tgene\t1\t100\t.\t+\t.\tgene_id "g";', ok)
  ref <- GenomeInfoDb::Seqinfo("chr1", seqlengths = 100)
  expect_identical(MMAPPR2:::.checkAnnotationAgainstReferenceFrozen(ok, ref), "chr1")
  bad <- file.path(td, "bad.gtf")
  writeLines('chr1\tx\tgene\t1\t101\t.\t+\t.\tgene_id "g";', bad)
  expect_error(MMAPPR2:::.checkAnnotationAgainstReferenceFrozen(bad, ref), "outside")
})

test_that("annotation sequence style can be harmonized to the FASTA", {
  td <- tempfile("ann_style_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  gtf <- file.path(td, "genes.gtf")
  writeLines('1\tx\tgene\t1\t90\t.\t+\t.\tgene_id "g";', gtf)
  ref <- GenomeInfoDb::Seqinfo("chr1", seqlengths = 100)
  expect_identical(MMAPPR2:::.checkAnnotationAgainstReferenceFrozen(gtf, ref), "chr1")
})

test_that("valid resources pass the complete frozen preflight", {
  fx <- .test_fixture()
  wt <- Rsamtools::BamFileList(Rsamtools::BamFile(fx$bam))
  mut <- Rsamtools::BamFileList(Rsamtools::BamFile(fx$bam))
  ref <- Rsamtools::FaFile(fx$fasta)
  expect_invisible(MMAPPR2:::.preflightInputResourcesFrozen(wt, mut, ref, fx$gtf))
})

test_that("resource setters accept compatible BAM replacements and output folder creation", {
  p <- .test_param()
  old <- wtFiles(p)
  expect_message(wtFiles(p) <- old, "Current index found")
  expect_message(mutFiles(p) <- mutFiles(p), "Current index found")
  newOut <- tempfile("setter_output_")
  expect_false(dir.exists(newOut))
  outputFolder(p) <- newOut
  expect_true(dir.exists(newOut))
  expect_identical(outputFolder(p), normalizePath(newOut))
})
