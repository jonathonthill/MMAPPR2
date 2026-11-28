test_that("safe plotting limits remain finite for constant or missing data", {
  for (fun in list(MMAPPR2:::.safeXLim, MMAPPR2:::.safeYLim)) {
    lim <- fun(rep(0, 5))
    expect_length(lim, 2)
    expect_true(all(is.finite(lim)))
    expect_lt(lim[1], lim[2])
    expect_equal(fun(c(NA_real_, Inf)), c(0, 1))
  }
})

test_that("temporary output paths do not collide", {
  a <- tempOutputFolder()
  b <- tempOutputFolder()
  expect_false(identical(a, b))
  expect_false(dir.exists(a))
  expect_false(dir.exists(b))
})

test_that("overwrite refuses R session temp root", {
  expect_error(MMAPPR2:::.prepareOutputFolder(tempdir(), overwrite = TRUE),
               "protected directory")
})


test_that("TSV serializer flattens multi-valued/list columns without changing numeric data", {
  x <- data.frame(score = c(1.25, 2.5), label = c("a\tb", "c\nd"),
                  stringsAsFactors = FALSE)
  x$annotations <- I(list(c("missense", "splice_region"), character()))
  x$positions <- I(list(c(10L, 11L), 20L))

  safe <- MMAPPR2:::.tsvSafeDataFrame(x)
  expect_false(any(vapply(safe, is.list, logical(1))))
  expect_equal(safe$score, c(1.25, 2.5))
  expect_identical(safe$label, c("a b", "c d"))
  expect_identical(safe$annotations, c("missense,splice_region", ""))
  expect_identical(safe$positions, c("10,11", "20"))

  out <- tempfile(fileext = ".tsv")
  on.exit(unlink(out), add = TRUE)
  expect_silent(MMAPPR2:::.writeTsv(x, out))
  lines <- readLines(out, warn = FALSE)
  expect_match(lines[1], "annotations")
  expect_match(lines[2], "missense,splice_region", fixed = TRUE)
})


test_that("atomic RDS save replaces output only after successful serialization", {
  td <- tempfile("atomic_rds_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  path <- file.path(td, "state.rds")
  saveRDS(list(version = "old"), path)
  expect_invisible(MMAPPR2:::.atomicSaveRDS(list(version = "new", x = 1:3), path))
  expect_identical(readRDS(path), list(version = "new", x = 1:3))
  expect_length(list.files(td, pattern = "^\\.state\\.rds\\."), 0L)
})


test_that("legacy detected-mutation table keeps original column schema", {
  x <- data.frame(
    seqnames = "chr1", start = 10L, end = 10L, width = 1L, strand = "+",
    ref = "A", alt = "G", mutRefDepth = 2L, mutAltDepth = 18L,
    peakDensity = 0.5, GENEID = "g1", TXID = "tx1",
    CONSEQUENCE = "nonsynonymous", REFAA = "K", VARAA = "R",
    extraModernMetadata = "retained elsewhere", stringsAsFactors = FALSE
  )
  x$PROTEINLOC <- I(list(c(10L, 11L)))
  got <- MMAPPR2:::.originalDetectedMutationTable(x)
  expect_identical(names(got), c(
    "seqnames", "start", "end", "width", "strand", "ref", "alt",
    "refDepth", "altDepth", "peakDensity", "GENEID", "TXID",
    "PROTEINLOC", "CONSEQUENCE", "REFAA", "VARAA"
  ))
  expect_identical(got$refDepth, 2L)
  expect_identical(got$altDepth, 18L)
  expect_identical(got$PROTEINLOC, "10")
})


test_that("candidate writer preserves legacy filenames and additive full-candidate output", {
  td <- tempfile("candidate_output_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  snps <- data.frame(seqnames = "chr1", start = 10L, end = 10L,
                     ref = "A", alt = "G", wtAltFreq = 0.1,
                     stringsAsFactors = FALSE)
  effects <- data.frame(
    seqnames = "chr1", start = 10L, end = 10L, width = 1L, strand = "+",
    ref = "A", alt = "G", mutRefDepth = 2L, mutAltDepth = 18L,
    peakDensity = 0.5, GENEID = "g1", TXID = "tx1", PROTEINLOC = 10L,
    CONSEQUENCE = "nonsynonymous", REFAA = "K", VARAA = "R",
    wtAltFreq = 0.1, stringsAsFactors = FALSE
  )
  diff <- data.frame(seqnames = "chr1", start = 1L, end = 100L,
                     gene_id = "g1", gene_name = "GeneOne",
                     ave_wt = 20, ave_mt = 50, log2FC = 1.322,
                     stringsAsFactors = FALSE)
  cand <- list(snps = list(chr1 = snps),
               effects = list(chr1 = effects),
               diff = list(chr1 = diff))

  expect_invisible(MMAPPR2:::.writeCandidateTables(cand, td))
  expect_true(file.exists(file.path(td, "DetectedMutationsForchr1.tsv")))
  expect_true(file.exists(file.path(td, "DifferentiallyExpressedGenesForchr1.tsv")))
  expect_true(file.exists(file.path(td, "AllCandidateVariantsForchr1.tsv")))

  legacy <- utils::read.delim(file.path(td, "DetectedMutationsForchr1.tsv"),
                              check.names = FALSE)
  expect_identical(names(legacy), c(
    "seqnames", "start", "end", "width", "strand", "ref", "alt",
    "refDepth", "altDepth", "peakDensity", "GENEID", "TXID",
    "PROTEINLOC", "CONSEQUENCE", "REFAA", "VARAA"
  ))
  full <- utils::read.delim(file.path(td, "AllCandidateVariantsForchr1.tsv"),
                            check.names = FALSE)
  expect_true("wtAltFreq" %in% names(full))
})


test_that("atomic RDS fallback replaces an existing file without deleting it first", {
  td <- tempfile("atomic_rds_fallback_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  path <- file.path(td, "state.rds")
  saveRDS(list(version = "old"), path)

  calls <- 0L
  renameFallback <- function(from, to) {
    calls <<- calls + 1L
    if (calls == 1L) return(FALSE)  # emulate Windows refusing overwrite-by-rename
    file.rename(from, to)
  }

  expect_invisible(MMAPPR2:::.atomicSaveRDS(
    list(version = "new"), path, .renameFile = renameFallback
  ))
  expect_identical(readRDS(path), list(version = "new"))
  expect_identical(calls, 3L)
  expect_length(list.files(td, all.files = TRUE,
                           pattern = "^\\.state\\.rds\\."), 0L)
})


test_that("atomic RDS fallback restores the previous file if installation fails", {
  td <- tempfile("atomic_rds_rollback_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  path <- file.path(td, "state.rds")
  old <- list(version = "old", sentinel = 42L)
  saveRDS(old, path)

  calls <- 0L
  renameRollback <- function(from, to) {
    calls <<- calls + 1L
    if (calls %in% c(1L, 3L)) return(FALSE)
    file.rename(from, to)
  }

  expect_error(
    MMAPPR2:::.atomicSaveRDS(list(version = "new"), path,
                             .renameFile = renameRollback),
    "previous file was restored"
  )
  expect_identical(readRDS(path), old)
  expect_identical(calls, 4L)
  expect_length(list.files(td, all.files = TRUE,
                           pattern = "^\\.state\\.rds\\."), 0L)
})


test_that("legacy mutation adapter handles Bioconductor List protein locations", {
  x <- S4Vectors::DataFrame(
    seqnames = "chr1", start = 10L, end = 10L, width = 1L, strand = "+",
    ref = "A", alt = "G", mutRefDepth = 2L, mutAltDepth = 18L,
    peakDensity = 0.5, GENEID = "g1", TXID = "tx1",
    PROTEINLOC = IRanges::IntegerList(c(10L, 11L)),
    CONSEQUENCE = "nonsynonymous", REFAA = "K", VARAA = "R"
  )
  got <- MMAPPR2:::.originalDetectedMutationTable(x)
  expect_identical(got$PROTEINLOC, "10")
})


test_that("legacy mutation adapter is safe for empty coding-effect results", {
  got <- MMAPPR2:::.originalDetectedMutationTable(data.frame())
  expect_identical(nrow(got), 0L)
  expect_identical(names(got), c(
    "seqnames", "start", "end", "width", "strand", "ref", "alt",
    "refDepth", "altDepth", "peakDensity", "GENEID", "TXID",
    "PROTEINLOC", "CONSEQUENCE", "REFAA", "VARAA"
  ))

  out <- tempfile(fileext = ".tsv")
  on.exit(unlink(out), add = TRUE)
  expect_invisible(MMAPPR2:::.writeTsv(got, out))
  header <- readLines(out, n = 1L, warn = FALSE)
  expect_match(header, "PROTEINLOC", fixed = TRUE)
})


test_that("TSV serializer flattens list columns even for zero-row tables", {
  x <- data.frame(id = integer())
  x$annotations <- I(vector("list", 0L))
  safe <- MMAPPR2:::.tsvSafeDataFrame(x)
  expect_identical(nrow(safe), 0L)
  expect_false(is.list(safe$annotations))

  out <- tempfile(fileext = ".tsv")
  on.exit(unlink(out), add = TRUE)
  expect_invisible(MMAPPR2:::.writeTsv(x, out))
  expect_match(readLines(out, n = 1L, warn = FALSE), "annotations", fixed = TRUE)
})


test_that("AICc plot export uses retained search results without refitting", {
  fx <- .make_synthetic_fixture()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  out <- file.path(fx$dir, "out_aicc")
  p <- mmapprParam(wtFiles = fx$bam, mutFiles = fx$bam,
                   refFasta = fx$fasta, gtf = fx$gtf,
                   outputFolder = out, includeScaffolds = TRUE,
                   minDepth = 1, minBaseQuality = 0, minMapQuality = 0,
                   exportAiccPlots = TRUE)
  md <- mmapprData(p)
  md@snpDistance <- list(chr1 = list(
    aicc = data.frame(spans = c(.1, .2, .3), aiccValues = c(4, 2, 3)),
    bestSpan = .2
  ))

  expect_invisible(MMAPPR2:::.plotAicc(md))
  path <- file.path(out, "aicc_plots.pdf")
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 0)
})
