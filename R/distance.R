#' @title Read BAM files and generate Euclidean distance data
#'
#' @name calculateDistance
#' @description
#' First step in the MMAPPR2 pipeline. Precedes the \code{\link{loessFit}}
#' step.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{snpDistance}
#'   slot filled.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#'
#'     md <- mmapprData(mmappr_param)
#'     postCalcDistMD <- calculateDistance(md)
#' }
#' 
NULL

calculateDistance <- function(mmapprData) {
  chrList <- suppressWarnings(.getFileReadChrList(param(mmapprData)))

  if (length(chrList) == 0L) stop("No annotated genomic ranges were available for pileup")

  mmapprData@snpDistance <-
    BiocParallel::bplapply(chrList, .calcDistForChr, param = mmapprData@param)

  return(mmapprData)
}


.getFileReadChrList <- function(param) {
  genes <- .annotationGeneRanges(param)

  if (length(genes) == 0L) stop("No gene ranges could be derived from the annotation")

  # Harmonize annotation naming to the FASTA where possible.
  fa_si <- .faSeqinfo(param@refGenome)
  targetStyle <- .choose_target_style(fa_si)
  if (!is.na(targetStyle)) {
    genes <- tryCatch(.setSeqlevelsStyleFrozen(genes, targetStyle),
                      error = function(e) genes)
  }

  common <- intersect(GenomeInfoDb::seqlevels(genes), GenomeInfoDb::seqlevels(fa_si))
  if (length(common) == 0L)
    stop("Annotation and FASTA do not share any sequence names after style harmonization")
  genes <- GenomeInfoDb::keepSeqlevels(genes, common, pruning.mode = "coarse")

  if (!isTRUE(param@includeScaffolds)) {
    # Standard-chromosome heuristics are useful for common model organisms
    # but can remove valid contigs from custom assemblies. Users of such assemblies
    # should set includeScaffolds=TRUE.
    genes <- GenomeInfoDb::keepStandardChromosomes(genes, pruning.mode = "coarse")
    if (length(genes) == 0L)
      stop("No standard chromosomes remained after filtering. For a custom/non-model assembly, try includeScaffolds=TRUE.")
  }

  mitoNames <- intersect(GenomeInfoDb::seqlevels(genes),
                         c("chrM", "MT", "M", "MtDNA", "mitochondrion_genome"))
  if (length(mitoNames) > 0L)
    genes <- GenomeInfoDb::dropSeqlevels(genes, mitoNames, pruning.mode = "coarse")
  if (length(genes) == 0L)
    stop("No non-mitochondrial annotated gene ranges remain for mapping")

  genes <- GenomicRanges::reduce(genes, ignore.strand = TRUE)
  split(genes, as.character(GenomicRanges::seqnames(genes)))
}


.calcDistForChr <- function(chrRange, param){
  tryCatch({
    stopifnot(length(unique(GenomicRanges::seqnames(chrRange))) == 1L)
    stopifnot(is(param, "MmapprParam"))
    stopifnot(length(chrRange) > 0L)

    # Keep sample identity long enough to make replicate aggregation semantics
    # explicit. In particular, coverage is averaged across *all* input files so
    # a replicate with zero coverage cannot silently inflate AVE.CVG.
    pileupWT <- lapply(seq_along(wtFiles(param)), function(i) {
      x <- .getPileup(wtFiles(param)[[i]], param, chrRange)
      if (nrow(x)) x[, FILE_ID := i]
      x
    })
    pileupWT <- data.table::rbindlist(pileupWT, use.names = TRUE, fill = TRUE)
    if (nrow(pileupWT) == 0L) stop("Insufficient data in wild-type file(s)")
    pileupWT <- .avgFiles(pileupWT,
                          fileAggregation = fileAggregation(param),
                          nFiles = length(wtFiles(param)))
    # "minimum depth" now means >= threshold rather than strictly >.
    pileupWT <- pileupWT[AVE.CVG >= minDepth(param)]
    if (nrow(pileupWT) == 0L) stop("Insufficient wild-type depth after filtering")

    pileupMut <- lapply(seq_along(mutFiles(param)), function(i) {
      x <- .getPileup(mutFiles(param)[[i]], param, chrRange)
      if (nrow(x)) x[, FILE_ID := i]
      x
    })
    pileupMut <- data.table::rbindlist(pileupMut, use.names = TRUE, fill = TRUE)
    if (nrow(pileupMut) == 0L) stop("Insufficient data in mutant file(s)")
    pileupMut <- .avgFiles(pileupMut,
                           fileAggregation = fileAggregation(param),
                           nFiles = length(mutFiles(param)))
    pileupMut <- pileupMut[AVE.CVG >= minDepth(param)]
    if (nrow(pileupMut) == 0L) stop("Insufficient mutant depth after filtering")

    data.table::setkey(pileupWT, CHROM, POS)
    data.table::setkey(pileupMut, CHROM, POS)
    distanceDf <- merge(pileupWT, pileupMut, by = c("CHROM", "POS"),
                        suffixes = c(".WT", ".MT"))

    # length(data.frame) counts columns, not rows.
    if (nrow(distanceDf) == 0L)
      stop("Empty dataframe after joining WT and mutant count tables")

    distanceDf[, DISTANCE := sqrt((AVE.A.FREQ.WT - AVE.A.FREQ.MT)^2 +
                                  (AVE.C.FREQ.WT - AVE.C.FREQ.MT)^2 +
                                  (AVE.G.FREQ.WT - AVE.G.FREQ.MT)^2 +
                                  (AVE.T.FREQ.WT - AVE.T.FREQ.MT)^2) ^
                                  distancePower(param)]

    # Retain markers that are polymorphic in the phenotypically WT F2 pool.
    # Select the four WT frequency columns explicitly. In data.table,
    # `AVE.A.FREQ.WT:AVE.T.FREQ.WT` evaluates `:` on the column vectors; it is
    # not a safe column-range selector and can produce nonsensical input to apply().
    wtFreqCols <- c("AVE.A.FREQ.WT", "AVE.C.FREQ.WT",
                    "AVE.G.FREQ.WT", "AVE.T.FREQ.WT")
    homozygousWT <- .wtHomozygousRows(distanceDf, wtFreqCols,
                                      cutoff = homozygoteCutoff(param))
    distanceDf <- distanceDf[!homozygousWT & is.finite(DISTANCE)]
    if (nrow(distanceDf) == 0L) stop("No informative markers remain after filtering")

    resultList <- list(wtCounts = pileupWT,
                       mutCounts = pileupMut,
                       distanceDf = distanceDf)

    .messageAndLog(paste("Finished", as.character(GenomicRanges::seqnames(chrRange)[1])),
                   outputFolder = outputFolder(param))
    resultList
  }, error = function(e) {
    paste(toString(unique(GenomicRanges::seqnames(chrRange))), e$message, sep = ": ")
  })
}


.wtHomozygousRows <- function(distanceDf, wtFreqCols, cutoff) {
  stopifnot(length(wtFreqCols) == 4L, all(wtFreqCols %in% names(distanceDf)))
  mat <- as.matrix(distanceDf[, ..wtFreqCols])
  # Frequencies should already be finite, but make malformed rows conservative:
  # a row with no finite frequency evidence is not useful as a mapping marker.
  major <- apply(mat, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) Inf else max(x)
  })
  major > cutoff
}


.primaryMappedScanFlag <- function() {
  Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE
  )
}


.getPileup <- function(file, param, chrRange) {
  stopifnot(length(file) == 1L)

  # simpleCigar=TRUE discarded ordinary spliced RNA-seq reads containing
  # CIGAR N operations, as well as clipped/indel-containing reads. pileup itself
  # understands reference skips, so do not throw away those reads wholesale.
  scanParam <- Rsamtools::ScanBamParam(
    flag = .primaryMappedScanFlag(),
    simpleCigar = FALSE,
    which = chrRange,
    mapqFilter = as.integer(minMapQuality(param))
  )

  pParam <- Rsamtools::PileupParam(
    max_depth = as.integer(maxPileupDepth(param)),
    min_mapq = as.integer(minMapQuality(param)),
    min_base_quality = as.integer(minBaseQuality(param)),
    distinguish_strands = FALSE,
    distinguish_nucleotides = TRUE,
    ignore_query_Ns = TRUE,
    # The linkage statistic is intentionally A/C/G/T based. Indel support
    # belongs in candidate calling and requires exact allele reconstruction; merely
    # counting '+'/'-' pileup symbols would not provide a valid normalized indel.
    include_deletions = FALSE,
    include_insertions = FALSE
  )

  raw <- Rsamtools::pileup(file, scanBamParam = scanParam, pileupParam = pParam)
  if (NROW(raw) == 0L) return(.emptyPileupTable())

  pData <- data.table::as.data.table(raw)
  pData <- pData[nucleotide %in% c("A", "C", "G", "T")]
  if (nrow(pData) == 0L) return(.emptyPileupTable())

  pData <- data.table::dcast(pData, seqnames + pos ~ nucleotide,
                             value.var = "count", fun.aggregate = sum, fill = 0)

  data.table::setnames(pData, c("seqnames", "pos"), c("CHROM", "POS"))
  # Small regions do not necessarily contain all four nucleotide categories.
  # Create absent columns explicitly rather than renaming by column position.
  for (base in c("A", "C", "G", "T")) {
    if (!base %in% names(pData)) pData[, (base) := 0]
  }
  data.table::setnames(pData, c("A", "C", "G", "T"),
                       c("A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ"))
  pData <- pData[, .(CHROM, POS, A.FREQ, C.FREQ, G.FREQ, T.FREQ)]

  pData[, CVG := A.FREQ + C.FREQ + G.FREQ + T.FREQ]
  pData <- pData[CVG > 0]
  pData[, c("A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ") :=
          list(A.FREQ/CVG, C.FREQ/CVG, G.FREQ/CVG, T.FREQ/CVG)]
  pData
}

.emptyPileupTable <- function() {
  data.table::data.table(CHROM = character(), POS = integer(),
                         A.FREQ = numeric(), C.FREQ = numeric(),
                         G.FREQ = numeric(), T.FREQ = numeric(),
                         CVG = numeric())
}


.avgFiles <- function(chrDf, fileAggregation, nFiles) {
  stopifnot(fileAggregation %in% c("simple", "weighted"))
  stopifnot(nFiles >= 1L)
  data.table::setkey(chrDf, CHROM, POS)

  if (fileAggregation == "weighted") {
    # Count-weighted allele frequencies are equivalent to pooling the reads.
    chrDf <- chrDf[, .(
      AVE.CVG = sum(CVG) / nFiles,
      AVE.A.FREQ = sum(A.FREQ * CVG) / sum(CVG),
      AVE.C.FREQ = sum(C.FREQ * CVG) / sum(CVG),
      AVE.G.FREQ = sum(G.FREQ * CVG) / sum(CVG),
      AVE.T.FREQ = sum(T.FREQ * CVG) / sum(CVG)
    ), .(CHROM, POS)]
  } else {
    # Equal-replicate weighting: frequency is averaged over replicates with data,
    # while coverage is averaged over all supplied files (missing coverage = 0).
    chrDf <- chrDf[, .(
      AVE.CVG = sum(CVG) / nFiles,
      AVE.A.FREQ = mean(A.FREQ),
      AVE.C.FREQ = mean(C.FREQ),
      AVE.G.FREQ = mean(G.FREQ),
      AVE.T.FREQ = mean(T.FREQ)
    ), .(CHROM, POS)]
  }
  chrDf
}
