#' @title Generate candidate mutations and consequences in peak regions
#'
#' @name generateCandidates
#' @description
#' Follows the \code{\link{peakRefinement}} step and produces a
#' \code{\linkS4class{MmapprData}} object ready for
#' \code{\link{outputMmapprData}}.
#'
#' @param md The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{candidates}
#'   slot filled with a \code{\link[GenomicRanges]{GRanges}} object for each
#'   peak chromosome containing variants and predicted consequences.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#'
#' \dontrun{
#' md <- mmapprData(mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#'
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' }
#'
NULL

# ---- main API ------------------------------------------------------------------

generateCandidates <- function(md) {
  if (length(md@peaks) == 0L) stop("No refined peaks are available for candidate generation")
  .messageAndLog("Getting variants in refined peak interval(s)", outputFolder(param(md)))
  peakGRanges <- lapply(md@peaks, .getPeakRange)

  .messageAndLog("Peak Summary:", outputFolder(param(md)))
  for (seqname in names(md@peaks)) {
    peak <- md@peaks[[seqname]]
    log_text <- paste0(
      "seqname: ", seqname, "\n",
      "refined interval: ", as.integer(peak$start), "-", as.integer(peak$end), "\n",
      "full-data LOESS apex: ", as.integer(round(peak$loessPeakPosition)), "\n",
      "resampling-density apex: ", as.integer(round(peak$densityPeakPosition)), "\n",
      "resampling success rate: ", round(peak$resampleSuccessRate, 4), "\n",
      "resampling seed: ", peak$resampleSeed, "\n"
    )
    .messageAndLog(log_text, outputFolder(param(md)))
  }

  md@candidates$snps <- lapply(peakGRanges, .getVariantsForRange, param = md@param)
  haveVariants <- any(vapply(md@candidates$snps, length, integer(1)) > 0L)

  txdb <- NULL
  if (haveVariants) {
    txdb <- .buildTxDb(md@param)
    on.exit(.disconnectTxDb(txdb), add = TRUE)
  }
  annotationGenes <- .annotationGeneRanges(md@param, txdb)

  .messageAndLog("Predicting coding variant effects", outputFolder(param(md)))
  if (haveVariants) {
    md@candidates$effects <- lapply(md@candidates$snps, .predictEffects,
                                    param = md@param, txdb = txdb)
  } else {
    md@candidates$effects <- lapply(md@candidates$snps,
                                    function(x) GenomicRanges::GRanges())
  }

  .messageAndLog("Summarizing expression changes in peak interval(s)", outputFolder(param(md)))
  md@candidates$diff <- lapply(peakGRanges, .addDiff,
                               param = md@param, txdb = txdb,
                               annotationGenes = annotationGenes)

  .messageAndLog("Scoring candidates by peak position", outputFolder(param(md)))
  md@candidates <- .scoreVariants(md@candidates, md@peaks)
  md
}


# ---- internals -----------------------------------------------------------------

.getPeakRange <- function(peakList) {
  GenomicRanges::GRanges(seqnames = peakList$seqname,
                         ranges = IRanges::IRanges(start = as.integer(peakList$start),
                                                  end = as.integer(peakList$end)))
}


.getVariantsForRange <- function(inputRange, param) {
  fa <- param@refGenome
  fa_si <- .faSeqinfo(fa)
  tstyle <- .choose_target_style(fa_si)
  if (!is.na(tstyle))
    inputRange <- tryCatch(.setSeqlevelsStyleFrozen(inputRange, tstyle),
                           error = function(e) inputRange)

  # Candidate pileup uses an adaptive bounded-memory caller. In the common case
  # (a modest peak on a machine with ample RAM) this is still one Rsamtools
  # pileup over the full interval, so supporting small machines does not impose
  # chunk-loop or temporary-file overhead on ordinary/high-memory runs. Large
  # peaks are processed in genomic chunks and reduced to passing candidate SNVs
  # immediately, so full-interval pileup tables never need to coexist in RAM.
  resultVr <- .candidateSnvsForRange(inputRange, param)
  if (length(resultVr) == 0L) return(resultVr)
  mutAF <- VariantAnnotation::altDepth(resultVr) / VariantAnnotation::totalDepth(resultVr)

  # Candidate discovery used to compare mutant reads only with the
  # reference. Add the phenotypically WT F2 pool as evidence at every candidate.
  # We report WT depth/AF for all candidates and expose optional WT/delta-AF
  # filters. Defaults are deliberately permissive because a recessive F2 WT pool
  # legitimately contains heterozygotes (~1/3 mutant-allele frequency at the locus).
  resultVr <- .addWtAlleleSupport(resultVr, param)
  S4Vectors::mcols(resultVr)$mutAltFreq <- mutAF
  S4Vectors::mcols(resultVr)$mutRefDepth <- VariantAnnotation::refDepth(resultVr)
  S4Vectors::mcols(resultVr)$mutAltDepth <- VariantAnnotation::altDepth(resultVr)
  S4Vectors::mcols(resultVr)$mutTotalDepth <- VariantAnnotation::totalDepth(resultVr)
  S4Vectors::mcols(resultVr)$deltaAltFreq <- mutAF - S4Vectors::mcols(resultVr)$wtAltFreq
  # Preserve the WT evidence even when the optional hard filters are
  # disabled. A high-frequency ALT allele in both pools is a strong hint that the
  # site is simply a background difference from the reference genome.
  wtEnough <- S4Vectors::mcols(resultVr)$wtTotalDepth >= minDepth(param)
  S4Vectors::mcols(resultVr)$wtEvidenceSufficient <- wtEnough
  S4Vectors::mcols(resultVr)$sharedBackgroundLike <-
    wtEnough & is.finite(S4Vectors::mcols(resultVr)$wtAltFreq) &
    S4Vectors::mcols(resultVr)$wtAltFreq > candidateMinAltFreq(param)

  wtAF <- S4Vectors::mcols(resultVr)$wtAltFreq
  deltaAF <- S4Vectors::mcols(resultVr)$deltaAltFreq
  wtFilterRequested <- candidateMaxWtAltFreq(param) < 1
  deltaFilterRequested <- candidateMinDeltaAF(param) > 0

  keep <- .candidateWtKeep(
    wtAF = wtAF, deltaAF = deltaAF,
    maxWtAltFreq = candidateMaxWtAltFreq(param),
    minDeltaAF = candidateMinDeltaAF(param),
    wtFilterRequested = wtFilterRequested,
    deltaFilterRequested = deltaFilterRequested
  )
  resultVr[keep]
}


.candidateMutantKeep <- function(totalDepth, altDepth, minDepth, minAltDepth, minAltFreq) {
  af <- ifelse(totalDepth > 0, altDepth / totalDepth, NA_real_)
  is.finite(totalDepth) & is.finite(altDepth) & is.finite(af) &
    totalDepth >= 0 & altDepth >= 0 & altDepth <= totalDepth &
    totalDepth >= minDepth & altDepth >= minAltDepth & af > minAltFreq
}

# Optional hard filters on WT ALT frequency and mutant-minus-WT ALT frequency
# are disabled by default (WT maximum 1; delta minimum 0), so candidate inclusion
# is not narrowed unless the user opts in.
.candidateWtKeep <- function(wtAF, deltaAF, maxWtAltFreq, minDeltaAF,
                             wtFilterRequested = maxWtAltFreq < 1,
                             deltaFilterRequested = minDeltaAF > 0) {
  keep <- rep(TRUE, length(wtAF))
  if (wtFilterRequested) keep <- keep & is.finite(wtAF) & wtAF <= maxWtAltFreq
  if (deltaFilterRequested) keep <- keep & is.finite(deltaAF) & deltaAF >= minDeltaAF
  keep
}


# Return candidate cgroup directories for this process. Container runtimes often
# mount the process cgroup at /sys/fs/cgroup itself, while systemd/CI hosts can
# expose a nested path from /proc/self/cgroup. Check both forms rather than
# assuming that memory.max lives at the mount root.
.processCgroupMemoryDirs <- function(version = c("v2", "v1")) {
  version <- match.arg(version)
  dirs <- if (identical(version, "v2")) "/sys/fs/cgroup" else "/sys/fs/cgroup/memory"
  lines <- tryCatch(readLines("/proc/self/cgroup", warn = FALSE),
                    error = function(e) character())
  if (!length(lines)) return(dirs)

  if (identical(version, "v2")) {
    hit <- grep("^0::", lines, value = TRUE)
    if (length(hit)) {
      rel <- sub("^0::/?", "", hit[[1L]])
      if (nzchar(rel)) dirs <- c(file.path("/sys/fs/cgroup", rel), dirs)
    }
  } else {
    parts <- strsplit(lines, ":", fixed = TRUE)
    hit <- Filter(function(x) length(x) >= 3L &&
                    "memory" %in% strsplit(x[[2L]], ",", fixed = TRUE)[[1L]],
                  parts)
    if (length(hit)) {
      rel <- sub("^/", "", hit[[1L]][[3L]])
      if (nzchar(rel)) {
        dirs <- c(file.path("/sys/fs/cgroup/memory", rel),
                  file.path("/sys/fs/cgroup", rel), dirs)
      }
    }
  }
  unique(dirs)
}


.cgroupMemoryRemaining <- function(dirs, limitName, currentName,
                                   unlimited = character()) {
  remaining <- numeric()
  for (dir in dirs) {
    limitFile <- file.path(dir, limitName)
    currentFile <- file.path(dir, currentName)
    if (!file.exists(limitFile) || !file.exists(currentFile)) next
    limTxt <- tryCatch(readLines(limitFile, n = 1L, warn = FALSE),
                       error = function(e) character())
    curTxt <- tryCatch(readLines(currentFile, n = 1L, warn = FALSE),
                       error = function(e) character())
    if (!length(limTxt) || !length(curTxt) || limTxt[[1L]] %in% unlimited) next
    limit <- suppressWarnings(as.numeric(limTxt[[1L]]))
    current <- suppressWarnings(as.numeric(curTxt[[1L]]))
    if (!is.finite(limit) || !is.finite(current) || limit <= 0) next
    # cgroup v1 commonly uses an enormous sentinel in place of "unlimited".
    if (identical(limitName, "memory.limit_in_bytes") && limit >= 2^60) next
    remaining <- c(remaining, max(1, limit - current))
  }
  remaining
}


# Return a conservative estimate of memory currently available to this process.
# Linux containers frequently expose the host's MemAvailable value even when a
# cgroup imposes a much smaller limit, so use the minimum of host and cgroup
# availability when both are present. NA means that the platform does not expose
# a reliable value and the caller should use a portable fallback.
.availableMemoryBytes <- function() {
  candidates <- numeric()

  if (file.exists("/proc/meminfo")) {
    mem <- tryCatch(readLines("/proc/meminfo", warn = FALSE),
                    error = function(e) character())
    line <- grep("^MemAvailable:", mem, value = TRUE)
    if (length(line)) {
      kb <- suppressWarnings(as.numeric(sub("^MemAvailable:\\s+([0-9]+).*", "\\1", line[[1L]])))
      if (is.finite(kb) && kb > 0) candidates <- c(candidates, kb * 1024)
    }
  }

  candidates <- c(candidates, .cgroupMemoryRemaining(
    .processCgroupMemoryDirs("v2"), "memory.max", "memory.current",
    unlimited = "max"
  ))
  candidates <- c(candidates, .cgroupMemoryRemaining(
    .processCgroupMemoryDirs("v1"), "memory.limit_in_bytes", "memory.usage_in_bytes"
  ))

  candidates <- candidates[is.finite(candidates) & candidates > 0]
  if (length(candidates)) min(candidates) else NA_real_
}


# Choose a genomic chunk width from available memory. Candidate calling is
# RNA-seq based and pileup rows are usually sparse relative to genomic width, so
# this is deliberately a planning heuristic rather than a claim about exact
# allocation. A 30% memory budget leaves headroom for Rsamtools, data.table
# aggregation/casting, the R heap, and the rest of the MmapprData object. The
# per-base allowance is intentionally stricter when less than 8 GiB is currently
# available because a densely covered locus can yield several nucleotide rows
# per base. At >=8 GiB the RNA-seq-sparse allowance keeps typical tens-of-Mb
# peaks on the one-shot path, avoiding a chunk-loop penalty on workstations and
# servers. Any unexpected allocation failure is still caught and retried in
# progressively smaller chunks.
.candidateAutoChunkSize <- function(nBams, availableBytes = .availableMemoryBytes()) {
  nBams <- max(1, as.integer(nBams))
  if (!is.finite(availableBytes) || availableBytes <= 0)
    return(10000000L)

  gib <- 1024^3
  bytesPerBasePerBam <- if (availableBytes < 2 * gib) {
    384
  } else if (availableBytes < 8 * gib) {
    256
  } else {
    128
  }
  budget <- availableBytes * 0.30
  estimate <- floor(budget / (bytesPerBasePerBam * nBams))
  estimate <- max(10000, min(250000000, estimate))
  as.integer(estimate)
}


.candidatePoolPlan <- function(which, nBams, mode = "auto", chunkSize = 0,
                               availableBytes = .availableMemoryBytes()) {
  mode <- match.arg(mode, c("auto", "memory", "chunked"))
  totalWidth <- sum(as.numeric(BiocGenerics::width(which)))
  if (!is.finite(totalWidth) || totalWidth < 1)
    return(list(chunked = FALSE, chunkSize = 1L, totalWidth = 0))

  adaptive <- if (is.numeric(chunkSize) && length(chunkSize) == 1L &&
                  is.finite(chunkSize) && chunkSize > 0) {
    as.integer(min(chunkSize, .Machine$integer.max))
  } else {
    .candidateAutoChunkSize(nBams, availableBytes = availableBytes)
  }

  list(
    chunked = identical(mode, "chunked") ||
      (identical(mode, "auto") && totalWidth > adaptive),
    chunkSize = max(1L, adaptive),
    totalWidth = totalWidth
  )
}


.splitCandidateRanges <- function(which, maxWidth) {
  if (length(which) == 0L) return(list())
  maxWidth <- max(1, as.integer(maxWidth))

  # First split only ranges that individually exceed the memory budget. Most
  # calls contain one refined peak, while WT-support calls may contain many
  # disjoint one-base candidate loci. Keeping short ranges in one GRanges object
  # avoids creating millions of one-element list entries for sparse candidates.
  widths <- as.numeric(BiocGenerics::width(which))
  longIdx <- which(widths > maxWidth)
  longRanges <- which[longIdx]
  segments <- if (length(longIdx)) which[-longIdx] else which
  if (length(longRanges)) {
    splitLong <- lapply(seq_along(longRanges), function(i) {
      left <- as.numeric(BiocGenerics::start(longRanges)[i])
      right <- as.numeric(BiocGenerics::end(longRanges)[i])
      starts <- seq(from = left, to = right, by = maxWidth)
      GenomicRanges::GRanges(
        seqnames = rep(as.character(GenomicRanges::seqnames(longRanges)[i]), length(starts)),
        ranges = IRanges::IRanges(start = starts,
                                  end = pmin(right, starts + maxWidth - 1L))
      )
    })
    splitLong <- do.call(c, splitLong)
    segments <- c(segments, splitLong)
  }
  if (!length(segments)) return(list())
  segments <- GenomicRanges::sort(segments, ignore.strand = TRUE)

  # Greedily group disjoint segments while keeping the sum of queried bases in
  # each chunk at or below maxWidth. This is important for WT support, where a
  # chunk can efficiently contain many sparse candidate loci.
  segWidths <- as.numeric(BiocGenerics::width(segments))
  group <- integer(length(segments))
  g <- 1L
  used <- 0
  for (i in seq_along(segWidths)) {
    if (used > 0 && used + segWidths[i] > maxWidth) {
      g <- g + 1L
      used <- 0
    }
    group[i] <- g
    used <- used + segWidths[i]
  }
  idx <- split(seq_along(segments), group)
  unname(lapply(idx, function(i) segments[i]))
}


.isMemoryAllocationError <- function(e) {
  msg <- conditionMessage(e)
  grepl(paste(c("cannot allocate", "vector memory exhausted", "memory exhausted",
                "std::bad_alloc", "cannot allocate memory"), collapse = "|"),
        msg, ignore.case = TRUE)
}


.candidateSnvsFromPileup <- function(pile, param) {
  if (nrow(pile) == 0L) return(VariantAnnotation::VRanges())
  vr <- .basePileupToVRanges(pile, sampleName = "pooled_mutant")
  if (length(vr) == 0L) return(vr)

  keep <- .candidateMutantKeep(
    totalDepth = VariantAnnotation::totalDepth(vr),
    altDepth = VariantAnnotation::altDepth(vr),
    minDepth = candidateMinDepth(param),
    minAltDepth = candidateMinAltDepth(param),
    minAltFreq = candidateMinAltFreq(param)
  )
  vr[keep]
}


.candidateSnvsOneRange <- function(inputRange, param,
                                   poolStrategy = c("memory", "incremental")) {
  poolStrategy <- match.arg(poolStrategy)
  pile <- .pooledBasePileup(
    bams = mutFiles(param), genome = param@refGenome, which = inputRange,
    minBaseQuality = minBaseQuality(param),
    minMapQuality = minMapQuality(param),
    maxDepth = .CANDIDATE_MAX_DEPTH,
    poolStrategy = poolStrategy
  )
  .candidateSnvsFromPileup(pile, param)
}


.combineCandidateVranges <- function(x) {
  x <- x[vapply(x, length, integer(1)) > 0L]
  if (!length(x)) return(VariantAnnotation::VRanges())
  out <- x[[1L]]
  if (length(x) > 1L) {
    for (i in 2:length(x)) out <- c(out, x[[i]])
  }
  out
}


# Process one chunk and recursively subdivide only if the allocation itself
# proves too large. This catches unusual high-coverage intervals even when the
# platform-level memory estimate looked generous.
.candidateSnvsChunkSafe <- function(inputRange, param) {
  tryCatch(
    .candidateSnvsOneRange(inputRange, param, poolStrategy = "incremental"),
    error = function(e) {
      if (!.isMemoryAllocationError(e) ||
          sum(as.numeric(BiocGenerics::width(inputRange))) <= 1)
        stop(e)

      half <- max(1L, as.integer(floor(sum(as.numeric(BiocGenerics::width(inputRange))) / 2)))
      subranges <- .splitCandidateRanges(inputRange, half)
      if (length(subranges) <= 1L) stop(e)
      .combineCandidateVranges(lapply(subranges, .candidateSnvsChunkSafe, param = param))
    }
  )
}


.candidateSnvsChunked <- function(inputRange, param, chunkSize,
                                  spillToDisk = FALSE) {
  chunks <- .splitCandidateRanges(inputRange, chunkSize)
  if (!length(chunks)) return(VariantAnnotation::VRanges())

  if (!isTRUE(spillToDisk)) {
    return(.combineCandidateVranges(
      lapply(chunks, .candidateSnvsChunkSafe, param = param)
    ))
  }

  # Disk staging is reserved for the emergency retry path. Only already-filtered
  # candidate VRanges are written, so temporary storage is much smaller than the
  # pileup tables that triggered the allocation problem. Unique tempfile() names
  # make concurrent MMAPPR2 runs safe.
  tempFiles <- character()
  on.exit(if (length(tempFiles)) unlink(tempFiles, force = TRUE), add = TRUE)
  for (chunk in chunks) {
    vr <- .candidateSnvsChunkSafe(chunk, param)
    if (!length(vr)) next
    path <- tempfile(pattern = "mmappr2_candidate_chunk_", fileext = ".rds")
    # Register the path before serialization so a partial file is still removed
    # if the temporary filesystem fills or saveRDS() is interrupted by an error.
    tempFiles <- c(tempFiles, path)
    tryCatch(
      saveRDS(vr, path),
      error = function(e) stop(
        "Unable to stage candidate results in the temporary directory '",
        tempdir(), "': ", conditionMessage(e), call. = FALSE
      )
    )
    rm(vr)
  }
  if (!length(tempFiles)) return(VariantAnnotation::VRanges())

  out <- readRDS(tempFiles[[1L]])
  if (length(tempFiles) > 1L) {
    for (path in tempFiles[-1L]) out <- c(out, readRDS(path))
  }
  out
}


.candidateSnvsForRange <- function(inputRange, param) {
  plan <- .candidatePoolPlan(
    inputRange, nBams = length(mutFiles(param)),
    mode = candidatePoolMode(param),
    chunkSize = candidateChunkSize(param)
  )

  if (isTRUE(plan$chunked)) {
    .messageAndLog(
      sprintf("Candidate pileup using bounded-memory chunks (<= %s bp)",
              format(plan$chunkSize, big.mark = ",", scientific = FALSE)),
      outputFolder(param)
    )
    return(.candidateSnvsChunked(inputRange, param, plan$chunkSize))
  }

  # Fast path: exactly one pooled pileup for the interval. If allocation still
  # fails, recover instead of aborting the run: retry in smaller chunks and stage
  # filtered chunk results in temp files so both the full pileup and all chunk
  # results are never resident simultaneously.
  tryCatch(
    .candidateSnvsOneRange(inputRange, param, poolStrategy = "memory"),
    error = function(e) {
      if (!.isMemoryAllocationError(e)) stop(e)
      fallbackSize <- max(1L, as.integer(min(
        plan$chunkSize,
        max(1, ceiling(plan$totalWidth / 4))
      )))
      .messageAndLog(
        paste0("Candidate pileup exceeded available memory; retrying with <= ",
               format(fallbackSize, big.mark = ",", scientific = FALSE),
               " bp chunks and temporary disk staging"),
        outputFolder(param)
      )
      .candidateSnvsChunked(inputRange, param, fallbackSize, spillToDisk = TRUE)
    }
  )
}


# Pooled A/C/G/T counts across one or more BAMs. Counts are summed before
# frequencies are calculated, which is the natural representation of a pooled
# sequencing library and avoids losing zero-ALT samples during replicate merge.
# Each call piles BAMs separately and sums A/C/G/T counts without creating a
# merged BAM. The caller keeps this fast in-memory operation for appropriately
# sized intervals and bounds its memory footprint by genomic chunking when the
# interval/system requires it. Counts are pooled before ALT selection, so reads
# with the reference base still contribute to total depth; for multiple BAMs,
# the explicit 250-read candidate cap applies per BAM.
.pooledBasePileup <- function(bams, genome, which = NULL,
                              minBaseQuality = 0L,
                              minMapQuality = 0L,
                              maxDepth = 1000L,
                              poolStrategy = c("memory", "incremental")) {
  poolStrategy <- match.arg(poolStrategy)
  # BamFileList() is dots-style in Rsamtools. Construct it explicitly for
  # character vectors so multiple BAM paths cannot be interpreted as one list-like
  # argument on versions with stricter S4 coercion.
  if (is.character(bams)) {
    bams <- do.call(Rsamtools::BamFileList, lapply(bams, Rsamtools::BamFile))
  }
  if (is(bams, "BamFile")) bams <- Rsamtools::BamFileList(bams)
  if (!is(bams, "BamFileList") || length(bams) == 0L)
    stop("bams must contain at least one BAM file")

  pupar <- Rsamtools::PileupParam(
    max_depth = as.integer(maxDepth),
    distinguish_nucleotides = TRUE,
    distinguish_strands = FALSE,
    min_base_quality = as.integer(minBaseQuality),
    min_mapq = as.integer(minMapQuality),
    ignore_query_Ns = TRUE,
    include_insertions = FALSE,
    include_deletions = FALSE
  )
  # Apply the same primary-alignment policy used by linkage pileup so a
  # chimeric/split read cannot contribute both its primary and supplementary
  # alignments to candidate allele depth. Secondary alignments are excluded too.
  sbpar <- if (!is.null(which))
    Rsamtools::ScanBamParam(flag = .primaryMappedScanFlag(),
                            which = which, simpleCigar = FALSE,
                            mapqFilter = as.integer(minMapQuality))
  else Rsamtools::ScanBamParam(flag = .primaryMappedScanFlag(),
                               simpleCigar = FALSE,
                               mapqFilter = as.integer(minMapQuality))

  pileOne <- function(i, aggregate = FALSE) {
    x <- Rsamtools::pileup(bams[[i]], scanBamParam = sbpar, pileupParam = pupar)
    if (NROW(x) == 0L) return(NULL)
    x <- data.table::as.data.table(x)
    x <- x[nucleotide %in% c("A", "C", "G", "T"),
           .(seqnames, pos, nucleotide, count)]
    if (!nrow(x)) return(NULL)
    if (isTRUE(aggregate))
      x <- x[, .(count = sum(count)), .(seqnames, pos, nucleotide)]
    x
  }

  if (identical(poolStrategy, "memory")) {
    # Fast path for machines with enough headroom: collect per-BAM tables and
    # aggregate them once. This preserves the lowest-overhead behavior for
    # ordinary analyses on workstations and servers.
    pieces <- lapply(seq_along(bams), pileOne, aggregate = FALSE)
    pieces <- pieces[!vapply(pieces, is.null, logical(1))]
    if (length(pieces) == 0L) return(.emptyPooledBasePileup())
    long <- data.table::rbindlist(pieces)
    long <- long[, .(count = sum(count)), .(seqnames, pos, nucleotide)]
  } else {
    # Low-memory path: keep at most the accumulated counts plus one BAM's
    # pileup table resident at a time. Re-aggregation after each BAM is slower
    # than the fast path, so this strategy is used only for bounded-memory work.
    long <- NULL
    for (i in seq_along(bams)) {
      x <- pileOne(i, aggregate = TRUE)
      if (is.null(x)) next
      if (is.null(long)) {
        long <- x
      } else {
        long <- data.table::rbindlist(list(long, x), use.names = TRUE)
        long <- long[, .(count = sum(count)), .(seqnames, pos, nucleotide)]
      }
      rm(x)
    }
    if (is.null(long) || !nrow(long)) return(.emptyPooledBasePileup())
  }
  wide <- data.table::dcast(long, seqnames + pos ~ nucleotide,
                            value.var = "count", fun.aggregate = sum, fill = 0)
  for (base in c("A", "C", "G", "T")) if (!base %in% names(wide)) wide[, (base) := 0L]
  wide <- wide[, .(seqnames, pos, A, C, G, T)]
  wide[, totalDepth := A + C + G + T]
  wide <- wide[totalDepth > 0]
  if (nrow(wide) == 0L) return(.emptyPooledBasePileup())

  gr <- GenomicRanges::GRanges(seqnames = wide$seqnames,
                               ranges = IRanges::IRanges(wide$pos, width = 1L))
  wide[, ref := as.character(Rsamtools::getSeq(genome, gr))]
  wide
}

.emptyPooledBasePileup <- function() {
  data.table::data.table(seqnames = character(), pos = integer(),
                         A = integer(), C = integer(), G = integer(), T = integer(),
                         totalDepth = integer(), ref = character())
}


.basePileupToVRanges <- function(pile, sampleName = "pooled") {
  if (nrow(pile) == 0L) return(VariantAnnotation::VRanges())
  bases <- c("A", "C", "G", "T")
  countMatrix <- as.matrix(pile[, ..bases])
  refIndex <- match(pile$ref, bases)
  refDepthVec <- rep(0L, nrow(pile))
  validRef <- !is.na(refIndex)
  if (any(validRef))
    refDepthVec[validRef] <- countMatrix[cbind(which(validRef), refIndex[validRef])]

  # data.table modifies by reference. Work on a copy so converting a pileup
  # to VRanges does not unexpectedly add an internal .row_id column to the caller's
  # object (important for tests and for future reuse of the pileup table).
  pile <- data.table::copy(pile)
  pile[, .row_id := .I]
  long <- data.table::melt(pile,
                           id.vars = c(".row_id", "seqnames", "pos", "totalDepth", "ref"),
                           measure.vars = bases,
                           variable.name = "alt", value.name = "altDepth")
  long[, alt := as.character(alt)]
  long <- long[alt != ref & altDepth > 0 & ref %in% bases]
  if (nrow(long) == 0L) return(VariantAnnotation::VRanges())
  # Keep the source vector name distinct from the data.table column name.
  # This avoids data.table scope ambiguity during := evaluation.
  long[, refDepth := refDepthVec[.row_id]]

  VariantAnnotation::VRanges(
    seqnames = S4Vectors::Rle(long$seqnames),
    ranges = IRanges::IRanges(start = long$pos, width = 1L),
    ref = long$ref,
    alt = as.character(long$alt),
    refDepth = S4Vectors::Rle(as.integer(long$refDepth)),
    altDepth = S4Vectors::Rle(as.integer(long$altDepth)),
    totalDepth = S4Vectors::Rle(as.integer(long$totalDepth)),
    sampleNames = S4Vectors::Rle(sampleName, nrow(long))
  )
}


.matchVariantToPileup <- function(variants, pile) {
  if (length(variants) == 0L) return(integer())
  if (nrow(pile) == 0L) return(rep(NA_integer_, length(variants)))
  varGr <- GenomicRanges::GRanges(
    seqnames = as.character(GenomicRanges::seqnames(variants)),
    ranges = IRanges::IRanges(BiocGenerics::start(variants), width = 1L)
  )
  pileGr <- GenomicRanges::GRanges(
    seqnames = as.character(pile$seqnames),
    ranges = IRanges::IRanges(as.integer(pile$pos), width = 1L)
  )
  as.integer(GenomicRanges::findOverlaps(varGr, pileGr, type = "equal",
                                         select = "first", ignore.strand = TRUE))
}


.wtDepthsFromPile <- function(variants, pile) {
  n <- length(variants)
  wtAltDepth <- integer(n)
  wtRefDepth <- integer(n)
  wtTotalDepth <- integer(n)
  if (!n || nrow(pile) == 0L)
    return(list(ref = wtRefDepth, alt = wtAltDepth, total = wtTotalDepth))

  idx <- .matchVariantToPileup(variants, pile)
  found <- !is.na(idx)
  if (any(found)) {
    bases <- c("A", "C", "G", "T")
    counts <- as.matrix(pile[idx[found], ..bases])
    altIdx <- match(as.character(VariantAnnotation::alt(variants))[found], bases)
    refIdx <- match(as.character(VariantAnnotation::ref(variants))[found], bases)
    rr <- seq_len(nrow(counts))
    okAlt <- !is.na(altIdx)
    okRef <- !is.na(refIdx)
    tempAlt <- integer(nrow(counts)); tempRef <- integer(nrow(counts))
    if (any(okAlt)) tempAlt[okAlt] <- counts[cbind(rr[okAlt], altIdx[okAlt])]
    if (any(okRef)) tempRef[okRef] <- counts[cbind(rr[okRef], refIdx[okRef])]
    wtAltDepth[found] <- tempAlt
    wtRefDepth[found] <- tempRef
    wtTotalDepth[found] <- pile$totalDepth[idx[found]]
  }
  list(ref = wtRefDepth, alt = wtAltDepth, total = wtTotalDepth)
}


.wtSupportChunkSafe <- function(variants, param, inputRange) {
  hits <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(
    variants, inputRange, ignore.strand = TRUE
  )))
  if (!length(hits))
    return(data.frame(index = integer(), ref = integer(), alt = integer(),
                      total = integer()))

  tryCatch({
    pile <- .pooledBasePileup(
      bams = wtFiles(param), genome = param@refGenome, which = inputRange,
      minBaseQuality = minBaseQuality(param),
      minMapQuality = minMapQuality(param),
      maxDepth = .CANDIDATE_MAX_DEPTH,
      poolStrategy = "incremental"
    )
    depths <- .wtDepthsFromPile(variants[hits], pile)
    data.frame(index = hits, ref = depths$ref, alt = depths$alt,
               total = depths$total)
  }, error = function(e) {
    if (!.isMemoryAllocationError(e) ||
        sum(as.numeric(BiocGenerics::width(inputRange))) <= 1)
      stop(e)

    half <- max(1L, as.integer(floor(
      sum(as.numeric(BiocGenerics::width(inputRange))) / 2
    )))
    subranges <- .splitCandidateRanges(inputRange, half)
    if (length(subranges) <= 1L) stop(e)
    do.call(rbind, lapply(subranges, function(x)
      .wtSupportChunkSafe(variants, param, x)))
  })
}


.addWtAlleleSupport <- function(variants, param, which = NULL) {
  if (length(variants) == 0L) return(variants)
  if (is.null(which)) {
    which <- unique(GenomicRanges::GRanges(
      seqnames = as.character(GenomicRanges::seqnames(variants)),
      ranges = IRanges::IRanges(BiocGenerics::start(variants), width = 1L)
    ))
    which <- GenomicRanges::reduce(which, ignore.strand = TRUE)
  }

  n <- length(variants)
  wtAltDepth <- integer(n)
  wtRefDepth <- integer(n)
  wtTotalDepth <- integer(n)

  plan <- .candidatePoolPlan(
    which, nBams = length(wtFiles(param)),
    mode = candidatePoolMode(param), chunkSize = candidateChunkSize(param)
  )

  applyChunkSupport <- function(chunks) {
    for (chunk in chunks) {
      support <- .wtSupportChunkSafe(variants, param, chunk)
      if (!nrow(support)) next
      wtRefDepth[support$index] <<- support$ref
      wtAltDepth[support$index] <<- support$alt
      wtTotalDepth[support$index] <<- support$total
    }
    invisible(NULL)
  }

  if (isTRUE(plan$chunked)) {
    .messageAndLog(
      sprintf("WT candidate support using bounded-memory chunks (<= %s queried bp)",
              format(plan$chunkSize, big.mark = ",", scientific = FALSE)),
      outputFolder(param)
    )
    applyChunkSupport(.splitCandidateRanges(which, plan$chunkSize))
  } else {
    fast <- tryCatch({
      wtPile <- .pooledBasePileup(
        bams = wtFiles(param), genome = param@refGenome, which = which,
        minBaseQuality = minBaseQuality(param),
        minMapQuality = minMapQuality(param),
        maxDepth = .CANDIDATE_MAX_DEPTH,
        poolStrategy = "memory"
      )
      .wtDepthsFromPile(variants, wtPile)
    }, error = function(e) {
      if (!.isMemoryAllocationError(e)) stop(e)
      NULL
    })

    if (!is.null(fast)) {
      wtRefDepth <- fast$ref
      wtAltDepth <- fast$alt
      wtTotalDepth <- fast$total
    } else {
      fallbackSize <- max(1L, as.integer(min(
        plan$chunkSize,
        max(1, ceiling(plan$totalWidth / 4))
      )))
      .messageAndLog(
        paste0("WT candidate support exceeded available memory; retrying with <= ",
               format(fallbackSize, big.mark = ",", scientific = FALSE),
               " queried-bp chunks"),
        outputFolder(param)
      )
      applyChunkSupport(.splitCandidateRanges(which, fallbackSize))
    }
  }

  wtAltFreq <- ifelse(wtTotalDepth > 0, wtAltDepth / wtTotalDepth, NA_real_)
  S4Vectors::mcols(variants)$wtRefDepth <- wtRefDepth
  S4Vectors::mcols(variants)$wtAltDepth <- wtAltDepth
  S4Vectors::mcols(variants)$wtTotalDepth <- wtTotalDepth
  S4Vectors::mcols(variants)$wtAltFreq <- wtAltFreq
  variants
}


.isKnownInternalRangeWarning <- function(w) {
  msg <- conditionMessage(w)
  call <- conditionCall(w)
  grepl("^GRanges object contains [0-9]+ out-of-bound ranges located on sequence ", msg) &&
    !is.null(call) && identical(as.character(call[[1L]]), "valid.GenomicRanges.seqinfo")
}

.validateRangesWithinReference <- function(x, refSeqinfo, label) {
  if (!length(x)) return(invisible(TRUE))
  seqn <- as.character(GenomeInfoDb::seqnames(x))
  refNames <- GenomeInfoDb::seqlevels(refSeqinfo)
  lens <- GenomeInfoDb::seqlengths(refSeqinfo)
  idx <- match(seqn, refNames)
  lim <- as.numeric(lens[idx])
  bad <- is.na(idx) | !is.finite(lim) |
    BiocGenerics::start(x) < 1L | BiocGenerics::end(x) > lim
  if (any(bad)) {
    example <- unique(paste0(seqn[bad], ":", BiocGenerics::start(x)[bad],
                             "-", BiocGenerics::end(x)[bad]))
    stop(label, " contains range(s) outside the indexed reference FASTA bounds: ",
         paste(utils::head(example, 5L), collapse = ", "))
  }
  invisible(TRUE)
}

.predictEffects <- function(inputVariants, param, txdb) {
  if (length(inputVariants) == 0L) return(GenomicRanges::GRanges())
  fa <- param@refGenome
  fa_si <- .faSeqinfo(fa)
  tstyle <- .choose_target_style(fa_si)

  # Harmonize sequence naming *before* intersecting seqlevels. Intersecting
  # first would incorrectly discard everything when, for example, the GTF uses
  # "1" while the FASTA/candidates use "chr1".
  txdbUse <- txdb
  vars <- inputVariants
  if (!is.na(tstyle)) {
    txdbUse <- tryCatch(.setSeqlevelsStyleFrozen(txdbUse, tstyle),
                        error = function(e) txdbUse)
    vars <- tryCatch(.setSeqlevelsStyleFrozen(vars, tstyle),
                     error = function(e) vars)
  }
  # Restrict the common set to all THREE objects. A FASTA/TxDb-wide list
  # can contain chromosomes that are not seqlevels of this peak's VRanges, and
  # passing those extra names to keepSeqlevels(vars, ...) can itself error.
  keep <- Reduce(intersect, list(GenomeInfoDb::seqlevels(txdbUse),
                                 GenomeInfoDb::seqlevels(vars),
                                 GenomeInfoDb::seqlevels(fa_si)))
  if (length(keep) == 0L)
    stop("No common sequence names between the annotation, reference FASTA, and candidate variants")
  txdbUse <- GenomeInfoDb::keepSeqlevels(txdbUse, keep, pruning.mode = "coarse")
  vars <- GenomeInfoDb::keepSeqlevels(vars, keep, pruning.mode = "coarse")
  if (length(vars) == 0L) return(GenomicRanges::GRanges())

  # Validate external genomic coordinates before predictCoding(). This turns
  # wrong-reference annotation problems into clear errors rather than opaque
  # internal transcript-mapping warnings.
  .validateRangesWithinReference(vars, fa_si, "Candidate variants")
  cds <- GenomicFeatures::cds(txdbUse)
  .validateRangesWithinReference(cds, fa_si, "Annotation CDS")

  # predictCoding() only describes coding consequences. The full SNV
  # candidate table is therefore still written separately for noncoding variants.
  effects <- withCallingHandlers(
    VariantAnnotation::predictCoding(
      query = vars,
      subject = txdbUse,
      seqSource = fa,
      varAllele = Biostrings::DNAStringSet(VariantAnnotation::alt(vars))
    ),
    warning = function(w) {
      if (.isKnownInternalRangeWarning(w)) invokeRestart("muffleWarning")
    }
  )

  # Carry mutant/WT evidence into the coding-effect table. predictCoding()
  # identifies each source query with QUERYID but does not guarantee that custom
  # metadata added during candidate calling survives in the result.
  if (length(effects) > 0L && "QUERYID" %in% colnames(S4Vectors::mcols(effects))) {
    qid <- as.integer(S4Vectors::mcols(effects)$QUERYID)
    sourceMeta <- S4Vectors::mcols(vars)
    evidenceCols <- intersect(
      c("mutAltFreq", "mutRefDepth", "mutAltDepth", "mutTotalDepth",
        "wtRefDepth", "wtAltDepth", "wtTotalDepth", "wtAltFreq",
        "deltaAltFreq", "wtEvidenceSufficient", "sharedBackgroundLike"),
      colnames(sourceMeta)
    )
    validQid <- !is.na(qid) & qid >= 1L & qid <= length(vars)
    for (nm in evidenceCols) {
      out <- rep(NA, length(effects))
      out[validQid] <- sourceMeta[[nm]][qid[validQid]]
      S4Vectors::mcols(effects)[[nm]] <- out
    }
  }
  effects
}


# Calculate descriptive log2 fold change from group mean counts. A small default
# pseudocount avoids Inf/NaN values when one or both groups have zero counts; callers
# can explicitly pass 0 when the raw count ratio is desired.
.expressionLog2FC <- function(ave_mt, ave_wt, pseudocount = 0.01) {
  if (pseudocount > 0)
    log2((ave_mt + pseudocount) / (ave_wt + pseudocount))
  else
    log2(ave_mt / ave_wt)
}


.poolMeanCounts <- function(countMat, num_wt, num_mut) {
  if (!is.matrix(countMat)) countMat <- as.matrix(countMat)
  if (num_wt < 1L || num_mut < 1L || ncol(countMat) != num_wt + num_mut)
    stop("countMat columns must equal num_wt + num_mut, with both groups non-empty")
  wt_idx <- seq_len(num_wt)
  mut_idx <- num_wt + seq_len(num_mut)
  list(wt = rowMeans(countMat[, wt_idx, drop = FALSE]),
       mut = rowMeans(countMat[, mut_idx, drop = FALSE]))
}


.extractQuotedAnnotationAttribute <- function(x, key) {
  pattern <- paste0("(?:^|;)[[:space:]]*", key, "[[:space:]]+\"([^\"]+)\"")
  hit <- regexec(pattern, x, perl = TRUE)
  parts <- regmatches(x, hit)
  vapply(parts, function(z) if (length(z) >= 2L) z[[2L]] else NA_character_,
         FUN.VALUE = character(1))
}

.extractEqualsAnnotationAttribute <- function(x, key) {
  pattern <- paste0("(?:^|;)[[:space:]]*", key, "[[:space:]]*=[[:space:]]*([^;]+)")
  hit <- regexec(pattern, x, perl = TRUE)
  parts <- regmatches(x, hit)
  out <- vapply(parts, function(z) if (length(z) >= 2L) z[[2L]] else NA_character_,
                FUN.VALUE = character(1))
  trimws(out)
}

.annotationGeneIdentity <- function(attributes) {
  # GTF convention first; fill missing identifiers from standard GFF3 ID/Name.
  gene_id <- .extractQuotedAnnotationAttribute(attributes, "gene_id")
  gff_id <- .extractEqualsAnnotationAttribute(attributes, "ID")
  gene_id[is.na(gene_id) | !nzchar(gene_id)] <- gff_id[is.na(gene_id) | !nzchar(gene_id)]

  gene_name <- .extractQuotedAnnotationAttribute(attributes, "gene_name")
  gff_name <- .extractEqualsAnnotationAttribute(attributes, "Name")
  missing_name <- is.na(gene_name) | !nzchar(gene_name)
  gene_name[missing_name] <- gff_name[missing_name]
  missing_name <- is.na(gene_name) | !nzchar(gene_name)
  gene_name[missing_name] <- gene_id[missing_name]
  list(gene_id = gene_id, gene_name = gene_name)
}

.annotationGeneRanges <- function(param, txdb = NULL) {
  tab <- tryCatch(
    data.table::fread(gtf(param), sep = "\t", header = FALSE,
                      select = c(1L, 3L, 4L, 5L, 7L, 9L), fill = TRUE,
                      quote = "", showProgress = FALSE, data.table = TRUE),
    error = function(e) NULL
  )
  if (!is.null(tab) && ncol(tab) == 6L) {
    data.table::setnames(tab,
                         c("seqnames", "type", "start", "end", "strand", "attributes"))
    tab[, seqnames := as.character(seqnames)]
    tab[, type := as.character(type)]
    tab[, start := suppressWarnings(as.numeric(start))]
    tab[, end := suppressWarnings(as.numeric(end))]
    tab[, strand := as.character(strand)]
    # GTF/GFF commonly uses "." for an unknown/unstranded feature. GRanges
    # represents that state as "*"; normalize it explicitly rather than allowing
    # an atypical annotation to fail during GRanges construction.
    tab[is.na(strand) | !(strand %in% c("+", "-", "*")), strand := "*"]
    tab[, attributes := as.character(attributes)]
    tab <- tab[!is.na(seqnames) & nzchar(trimws(seqnames)) &
               !grepl("^[[:space:]]*#", seqnames) & !is.na(type) & type == "gene" &
               is.finite(start) & is.finite(end) & start >= 1 & end >= start]
    if (nrow(tab)) {
      identity <- .annotationGeneIdentity(tab$attributes)
      df <- data.frame(seqnames = tab$seqnames,
                       start = as.integer(tab$start),
                       end = as.integer(tab$end),
                       strand = tab$strand,
                       gene_id = identity$gene_id,
                       gene_name = identity$gene_name,
                       stringsAsFactors = FALSE)
      return(GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE))
    }
  }

  # New fallback for GFF3, compressed/atypical annotations, or files that fread
  # cannot parse directly. Explicit gene rows remain preferred here as well.
  ann <- tryCatch(rtracklayer::import(gtf(param)), error = function(e) NULL)
  if (!is.null(ann) && length(ann)) {
    type <- S4Vectors::mcols(ann)$type
    if (!is.null(type)) {
      type <- as.character(type)
      genes <- ann[!is.na(type) & type == "gene"]
      if (length(genes)) {
        mc <- S4Vectors::mcols(genes)
        ids <- if ("gene_id" %in% colnames(mc)) as.character(mc$gene_id) else rep(NA_character_, length(genes))
        fallbackIds <- if ("ID" %in% colnames(mc)) as.character(mc$ID) else names(genes)
        if (is.null(fallbackIds)) fallbackIds <- rep(NA_character_, length(genes))
        missingId <- is.na(ids) | !nzchar(ids)
        ids[missingId] <- fallbackIds[missingId]

        geneNames <- if ("gene_name" %in% colnames(mc)) as.character(mc$gene_name) else rep(NA_character_, length(genes))
        fallbackNames <- if ("Name" %in% colnames(mc)) as.character(mc$Name) else rep(NA_character_, length(genes))
        missingName <- is.na(geneNames) | !nzchar(geneNames)
        geneNames[missingName] <- fallbackNames[missingName]
        missingName <- is.na(geneNames) | !nzchar(geneNames)
        geneNames[missingName] <- ids[missingName]

        S4Vectors::mcols(genes)$gene_id <- ids
        S4Vectors::mcols(genes)$gene_name <- geneNames
        return(genes)
      }
    }
  }

  # Last-resort hardening for annotations without explicit gene rows. Build a
  # TxDb only if the caller does not already have one; ordinary GTFs therefore
  # avoid this comparatively expensive SQLite construction in the linkage stage.
  ownTxdb <- is.null(txdb)
  if (ownTxdb) {
    txdb <- .buildTxDb(param)
    on.exit(.disconnectTxDb(txdb), add = TRUE)
  }
  genes <- GenomicFeatures::genes(txdb)
  if (!length(genes)) return(genes)
  S4Vectors::mcols(genes)$gene_id <- names(genes)
  S4Vectors::mcols(genes)$gene_name <- names(genes)
  genes
}

# The descriptive expression summary uses annotation gene spans within the refined peak,
# harmonizes seqnames, applies the primary-read filter, and honors paired-end and strand
# controls. Group means are kept separate, and log2 fold change uses the configured
# pseudocount (0.01 by default). This remains a descriptive raw-count summary, not
# formal differential-expression inference.
.addDiff <- function(peakGRange, param, txdb, annotationGenes = NULL) {
  genes <- if (is.null(annotationGenes)) .annotationGeneRanges(param, txdb) else annotationGenes
  if (!length(genes)) return(GenomicRanges::GRanges())

  # Match annotation/peak naming before overlap tests; old gene-span features
  # are retained instead of replacing them with exon-group features.
  peakUse <- peakGRange
  geneUse <- genes
  peakStyle <- .choose_target_style(.faSeqinfo(param@refGenome))
  if (!is.na(peakStyle)) {
    geneUse <- tryCatch(.setSeqlevelsStyleFrozen(geneUse, peakStyle), error = function(e) geneUse)
    peakUse <- tryCatch(.setSeqlevelsStyleFrozen(peakUse, peakStyle), error = function(e) peakUse)
  }
  common <- intersect(GenomeInfoDb::seqlevels(geneUse), GenomeInfoDb::seqlevels(peakUse))
  if (!length(common)) return(GenomicRanges::GRanges())
  geneUse <- GenomeInfoDb::keepSeqlevels(geneUse, common, pruning.mode = "coarse")
  peakUse <- GenomeInfoDb::keepSeqlevels(peakUse, common, pruning.mode = "coarse")
  geneUse <- IRanges::subsetByOverlaps(geneUse, peakUse, ignore.strand = TRUE)
  if (!length(geneUse)) return(geneUse)

  bamObjects <- c(
    lapply(seq_along(param@wtFiles), function(i) param@wtFiles[[i]]),
    lapply(seq_along(param@mutFiles), function(i) param@mutFiles[[i]])
  )
  readfiles <- do.call(Rsamtools::BamFileList, bamObjects)

  # Preserve old gene-span/peak-window counting semantics. The only counting
  # hardening here is excluding secondary/supplementary/unmapped records so the
  # same physical alignment cannot inflate expression evidence.
  counts <- GenomicAlignments::summarizeOverlaps(
    features = geneUse,
    reads = readfiles,
    mode = "Union",
    ignore.strand = ignoreStrand(param),
    singleEnd = !pairedEnd(param),
    fragments = pairedEnd(param),
    param = Rsamtools::ScanBamParam(which = peakUse,
                                    flag = .primaryMappedScanFlag())
  )
  countMat <- SummarizedExperiment::assays(counts)$counts
  num_wt <- length(param@wtFiles)
  num_mut <- length(param@mutFiles)
  if (ncol(countMat) != num_wt + num_mut)
    stop("Unexpected number of columns returned by summarizeOverlaps")

  means <- .poolMeanCounts(countMat, num_wt = num_wt, num_mut = num_mut)
  ave_wt <- means$wt
  ave_mt <- means$mut
  # Add the configured pseudocount to both group means before taking the ratio.
  # The default is 0.01 to keep zero-count genes finite; setting it to 0 restores
  # the raw ratio when that behavior is explicitly desired.
  pc <- expressionPseudocount(param)
  log2FC <- .expressionLog2FC(ave_mt, ave_wt, pseudocount = pc)

  S4Vectors::mcols(geneUse)$ave_wt <- ave_wt
  S4Vectors::mcols(geneUse)$ave_mt <- ave_mt
  S4Vectors::mcols(geneUse)$log2FC <- round(log2FC, 3)

  geneUse[(abs(S4Vectors::mcols(geneUse)$log2FC) > 1 |
             is.na(S4Vectors::mcols(geneUse)$log2FC)) &
            (S4Vectors::mcols(geneUse)$ave_wt > 10 |
             S4Vectors::mcols(geneUse)$ave_mt > 10)]
}


.scoreVariants <- function(candList, peaks) {
  for (groupName in names(candList)) {
    group <- candList[[groupName]]
    if (!is.list(group)) next
    for (seqname in names(group)) {
      if (is.null(peaks[[seqname]]) || length(group[[seqname]]) == 0L) next
      densityFunc <- peaks[[seqname]]$densityFunction
      positions <- BiocGenerics::start(group[[seqname]]) +
        ((BiocGenerics::width(group[[seqname]]) - 1) / 2)
      densityCol <- vapply(positions, densityFunc, FUN.VALUE = numeric(1))
      S4Vectors::mcols(candList[[groupName]][[seqname]])$peakDensity <- densityCol
      candList[[groupName]][[seqname]] <- .orderVariants(candList[[groupName]][[seqname]])
    }
  }
  candList
}


.orderVariants <- function(candidateGRanges) {
  if (!inherits(candidateGRanges, c("VRanges", "GRanges"))) {
    warning("Invalid data type for candidate sorting: ", class(candidateGRanges)[1])
    return(candidateGRanges)
  }
  if (length(candidateGRanges) == 0L) return(candidateGRanges)

  consequence <- S4Vectors::mcols(candidateGRanges)$CONSEQUENCE
  density <- S4Vectors::mcols(candidateGRanges)$peakDensity
  if (!is.null(consequence)) {
    # This remains a deliberately simple coding-consequence ordering; a
    # full ontology/functional-impact model is a larger scientific redesign.
    impactLevels <- c("synonymous", "nonsynonymous", "frameshift", "nonsense")
    severity <- match(as.character(consequence), impactLevels)
    severity[is.na(severity)] <- 0L
    orderVec <- order(severity, density, na.last = TRUE, decreasing = TRUE)
  } else {
    orderVec <- order(density, na.last = TRUE, decreasing = TRUE)
  }
  candidateGRanges[orderVec]
}


