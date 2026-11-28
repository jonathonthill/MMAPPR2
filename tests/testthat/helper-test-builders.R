.mmappr2_test_cache <- new.env(parent = emptyenv())

# Keep ordinary tests deterministic and avoid inheriting a host-specific
# parallel backend. Dedicated tests explicitly switch to a parallel backend when
# validating serial/parallel equivalence.
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

.test_fixture <- function() {
  if (!exists("basic", envir = .mmappr2_test_cache, inherits = FALSE)) {
    assign("basic", .make_synthetic_fixture(), envir = .mmappr2_test_cache)
  }
  get("basic", envir = .mmappr2_test_cache, inherits = FALSE)
}

.test_param <- function(fx = .test_fixture(), ...) {
  args <- list(
    wtFiles = fx$bam,
    mutFiles = fx$bam,
    refFasta = fx$fasta,
    gtf = fx$gtf,
    outputFolder = tempfile("mmappr2_test_output_"),
    minDepth = 1,
    minBaseQuality = 0,
    minMapQuality = 0,
    peakResampleIterations = 10
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(MMAPPR2::mmapprParam, args)
}

.write_one_base_sam <- function(path, seqname, seqlength, positions, alt_counts,
                                depth = 20L, ref = "A", alt = "C",
                                skip_positions = integer()) {
  stopifnot(length(positions) == length(alt_counts))
  header <- c("@HD\tVN:1.6\tSO:coordinate",
              sprintf("@SQ\tSN:%s\tLN:%d", seqname, seqlength))
  rows <- character()
  q <- 0L
  for (i in seq_along(positions)) {
    pos <- positions[[i]]
    if (pos %in% skip_positions) next
    n_alt <- as.integer(alt_counts[[i]])
    n_ref <- as.integer(depth - n_alt)
    if (n_ref < 0L || n_alt < 0L) stop("invalid synthetic allele counts")
    bases <- c(rep(ref, n_ref), rep(alt, n_alt))
    for (base in bases) {
      q <- q + 1L
      rows <- c(rows, paste0("r", q, "\t0\t", seqname, "\t", pos,
                            "\t60\t1M\t*\t0\t0\t", base, "\tI"))
    }
  }
  writeLines(c(header, rows), path)
}

.make_signal_fixture <- function(wt_zero_at_candidate = FALSE) {
  key <- if (wt_zero_at_candidate) "signal_zero_wt" else "signal"
  if (exists(key, envir = .mmappr2_test_cache, inherits = FALSE))
    return(get(key, envir = .mmappr2_test_cache, inherits = FALSE))

  td <- tempfile(paste0("mmappr2_", key, "_"))
  dir.create(td, recursive = TRUE)
  seqname <- "chr1"
  seqlength <- 1200L
  candidate_pos <- 500L

  fa <- file.path(td, "ref.fa")
  writeLines(c(">chr1", paste(rep("A", seqlength), collapse = "")), fa)
  Rsamtools::indexFa(fa)

  gtf <- file.path(td, "genes.gtf")
  writeLines(c(
    'chr1\tMMAPPR2\tgene\t1\t1100\t.\t+\t.\tgene_id "g1"; gene_name "Gene1";',
    'chr1\tMMAPPR2\ttranscript\t1\t1100\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";',
    'chr1\tMMAPPR2\texon\t1\t1100\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; exon_id "e1";',
    'chr1\tMMAPPR2\tCDS\t2\t1099\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";'
  ), gtf)

  positions <- seq.int(100L, 890L, by = 10L) # 80 mapping markers
  depth <- 100L
  wt_alt <- rep(40L, length(positions))       # 40% ALT, safely polymorphic
  d <- abs(positions - candidate_pos)
  # A smooth, broad linkage signal keeps half-sample robust LOESS fits stable
  # while reserving >80% mutant AF for the planted candidate itself.
  mut_alt <- as.integer(round(40 + 39 * exp(-(d / 160)^2)))
  mut_alt[positions == candidate_pos] <- 95L  # unique > 0.80 candidate

  wt_sam <- file.path(td, "wt.sam")
  mut_sam <- file.path(td, "mut.sam")
  .write_one_base_sam(wt_sam, seqname, seqlength, positions, wt_alt, depth = depth,
                      skip_positions = if (wt_zero_at_candidate) candidate_pos else integer())
  .write_one_base_sam(mut_sam, seqname, seqlength, positions, mut_alt, depth = depth)

  wt_bam <- as.character(Rsamtools::asBam(wt_sam, destination = file.path(td, "wt"), overwrite = TRUE))[1]
  mut_bam <- as.character(Rsamtools::asBam(mut_sam, destination = file.path(td, "mut"), overwrite = TRUE))[1]
  if (!file.exists(paste0(wt_bam, ".bai")) && !file.exists(sub("\\.bam$", ".bai", wt_bam))) Rsamtools::indexBam(wt_bam)
  if (!file.exists(paste0(mut_bam, ".bai")) && !file.exists(sub("\\.bam$", ".bai", mut_bam))) Rsamtools::indexBam(mut_bam)

  out <- list(dir = td, fasta = fa, gtf = gtf, wt = wt_bam, mut = mut_bam,
              candidate_pos = candidate_pos, positions = positions)
  assign(key, out, envir = .mmappr2_test_cache)
  out
}

.signal_param <- function(wt_zero_at_candidate = FALSE, ...) {
  fx <- .make_signal_fixture(wt_zero_at_candidate)
  args <- list(
    wtFiles = fx$wt,
    mutFiles = fx$mut,
    refFasta = fx$fasta,
    gtf = fx$gtf,
    outputFolder = tempfile("mmappr2_signal_output_"),
    minDepth = 10,
    homozygoteCutoff = 0.95,
    minBaseQuality = 0,
    minMapQuality = 0,
    distancePower = 1,
    peakCutoffSd = 0,
    peakResampleIterations = 20,
    randomSeed = 17,
    candidateMinDepth = 10,
    candidateMinAltDepth = 10,
    candidateMinAltFreq = 0.80,
    expressionPseudocount = 0.01,
    exportAiccPlots = TRUE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(MMAPPR2::mmapprParam, args)
}

.make_loess_result <- function(seqname = "chr1", n = 80L,
                               peak = TRUE, offset = 0) {
  pos <- seq_len(n) * 10 + offset
  y <- if (peak) exp(-((pos - stats::median(pos)) / 90)^2) else rep(0.05, n)
  fit <- stats::loess(y ~ pos, span = 0.75, degree = 1, family = "gaussian")
  list(seqname = seqname, loess = fit,
       bestSpan = fit$pars$span,
       aicc = data.frame(spans = c(.2, .35, .5), aiccValues = c(3, 1, 2)),
       wtCounts = data.table::data.table(), mutCounts = data.table::data.table())
}

.make_peak_entry <- function(seqname = "chr1", start = 300L, end = 700L,
                             center = 500) {
  x <- seq(start, end, length.out = 51)
  y <- stats::dnorm(x, mean = center, sd = max(1, (end - start) / 8))
  list(seqname = seqname, start = start, end = end,
       densityFunction = stats::approxfun(x, y, yleft = 0, yright = 0),
       peakPosition = center, densityPeakPosition = center,
       loessPeakPosition = center, densityData = list(x = x, y = y),
       resampleSuccessRate = 1, resampleSeed = 1L,
       cutoff = 0.2, cutoffCenter = 0.1, cutoffSpread = 0.05,
       cutoffMethod = "legacy_current", intervalMethod = "hpd_span")
}

.make_quality_fixture <- function() {
  if (exists("quality", envir = .mmappr2_test_cache, inherits = FALSE))
    return(get("quality", envir = .mmappr2_test_cache, inherits = FALSE))

  td <- tempfile("mmappr2_quality_")
  dir.create(td, recursive = TRUE)
  fa <- file.path(td, "ref.fa")
  writeLines(c(">chr1", paste(rep("A", 100), collapse = "")), fa)
  Rsamtools::indexFa(fa)
  gtf <- file.path(td, "genes.gtf")
  writeLines('chr1\tMMAPPR2\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; gene_name "G1";', gtf)

  sam <- file.path(td, "reads.sam")
  rows <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:100",
    # Three observations at position 10 isolate mapping-quality and base-quality filters.
    "highA\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\tI",
    "lowMapC\t0\tchr1\t10\t5\t1M\t*\t0\t0\tC\tI",
    "lowBaseG\t0\tchr1\t10\t60\t1M\t*\t0\t0\tG\t!",
    # Five reads at position 20 exercise Rsamtools max_depth behavior.
    paste0("depth", 1:5, "\t0\tchr1\t20\t60\t1M\t*\t0\t0\tT\tI")
  )
  writeLines(rows, sam)
  bam <- as.character(Rsamtools::asBam(sam, destination = file.path(td, "reads"), overwrite = TRUE))[1]
  if (!file.exists(paste0(bam, ".bai")) && !file.exists(sub("\\.bam$", ".bai", bam)))
    Rsamtools::indexBam(bam)

  out <- list(dir = td, fasta = fa, gtf = gtf, bam = bam)
  assign("quality", out, envir = .mmappr2_test_cache)
  out
}
