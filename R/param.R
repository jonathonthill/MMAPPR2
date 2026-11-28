#' @title MmapprParam Class
#' 
#' @name MmapprParam
#' @description
#' \code{MmapprParam} stores parameters for running \code{\link{mmappr}}.
#'
#' @slot wtFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the wild-type pool to be analyzed.
#' @slot mutFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the mutant pool to be analyzed.
#' @slot refFasta The path to the reference FASTA file.
#' @slot gtf The path to a gtf-formatted annotation file. 
#' @slot outputFolder Length-one character vector specifying where to save
#'   output, including a \code{\linkS4class{MmapprData}} stored as
#'   \code{mmappr_data.RDS}, \code{mmappr2.log}, a \code{.tsv} file
#'   for each peak chromosome containing candidate mutations, and PDF plots
#'   of both the entire genome and peak chromosomes. Defaults to an
#'   automatically generated \code{mmappr2_<timestamp>}.
#' @slot includeScaffolds Logical indicating whether non-standard chromosomes
#'   should be included. FALSE (default) trims the seqnames of the reference to 
#'   include only standard chromosome names based on UCSC or Ensembl 
#'   conventions. In all cases, the mitochondrial chromosome is removed.
#' @slot minDepth Length-one integer vector determining minimum depth
#'   required for a position to
#'   be considered in the analysis. Defaults to 20.
#' @slot homozygoteCutoff Length-one numeric vector between \code{0} and
#'   \code{1} specifying threshold for throwing out base pairs on account
#'   of homozygosity. Positions with high major allele frequency in the
#'   wild-type pool are unlikely to exhibit polymorphism and are thus thrown
#'   out when they exceed this cutoff. Defaults to \code{0.95}.
#' @slot minBaseQuality Length-one numeric vector indicating minimum base
#'   call quality to consider in analysis. Read positions with qualities
#'   below this score will be thrown out. Defaults to 20.
#' @slot minMapQuality Length-one numeric vector indicating minimum read
#'   mapping quality to consider in analysis. Reads with qualities below
#'   this score will be thrown out. Defaults to 20.
#' @slot fileAggregation A length-one character vector determining strategy
#'   for aggregating base calls when multiple wild-type or multiple mutant
#'   files are provided.
#'   When 'weighted', allele frequencies are weighted by read depth (pooled-read
#'   behavior). When 'simple', replicate allele frequencies are weighted equally
#'   while coverage is averaged across all supplied files.
#' @slot distancePower Length-one numeric vector determining to what power
#'   Euclidean distance values are raised before fitting. Higher powers tend
#'   to increase high values and decrease low values, exaggerating the
#'   variation in the data. Default of 4.
#' @slot peakIntervalWidth Length-one numeric vector between \code{0} and
#'   \code{1} specifying desired width of linkage region(s). The default value
#'   of \code{0.80} targets 80\% of the resampling-density mass. The exact
#'   interval construction is controlled by \code{peakIntervalMethod}; the default
#'   \code{hpd_span} conservatively spans all selected high-density modes.
#' @slot loessOptResolution Length-one numeric vector between \code{0} and
#'   \code{1} specifying
#'   desired resolution for Loess fit optimization. The default of \code{0.001},
#'   for example, indicates that the span ultimately chosen will perform better
#'   than its neighbor values at \code{+-0.001}.
#' @slot loessOptCutFactor Length-one numeric vector between \code{0} and
#'   \code{1} specifying how aggressively the Loess
#'   optimization algorithm proceeds. For example, with a default of \code{0.1}
#'   different spans at intervals of \code{0.001} would be evaluated after
#'   intervals of \code{0.01}.
#' @slot maxPileupDepth Maximum reads considered per BAM at one pileup position.
#' @slot candidateMinDepth Optional pooled mutant depth floor for final candidates. The default (1) is effectively inert given the two-ALT-read rule.
#' @slot candidateMinAltDepth Minimum pooled mutant reads supporting the alternative allele.
#' @slot candidateMinAltFreq Strict lower bound on mutant alternative-allele frequency for a candidate; the observed frequency must be greater than this value.
#' @slot candidateMaxWtAltFreq Optional maximum WT alternative-allele frequency; 1 disables the filter.
#' @slot candidateMinDeltaAF Optional minimum mutant-minus-WT allele-frequency difference; 0 disables the filter.
#' @slot candidatePoolMode Candidate-pileup memory strategy. \code{auto} uses
#'   available-memory-aware chunking only when useful, \code{memory} prefers a
#'   one-shot in-memory pileup, and \code{chunked} always processes bounded
#'   genomic chunks. Allocation failures are retried with smaller chunks.
#' @slot candidateChunkSize Maximum candidate-pileup chunk width in bases.
#'   A value of 0 lets the package choose an adaptive size from available memory
#'   and the number of BAM files.
#' @slot peakCutoffSd Multiplier applied to the variability term used for initial peak detection.
#' @slot peakCutoffMethod Peak-threshold calculation: 'legacy_current' uses the old mean-of-chromosome-means plus sqrt(sum(var/n)) formula; 'global_sd' uses the median plus a genome-wide SD multiple.
#' @slot peakIntervalMethod Refined-interval strategy: 'hpd_span' conservatively spans all high-density modes; 'shortest_contiguous' selects the shortest single interval with the requested mass.
#' @slot peakResampleIterations Number of half-marker LOESS resamples used for peak refinement.
#' @slot randomSeed Base seed used to make peak resampling reproducible.
#' @slot pairedEnd Whether RNA-seq reads should be counted as paired-end fragments in the descriptive expression summary.
#' @slot ignoreStrand Whether strand should be ignored in the descriptive expression summary. FALSE uses strand-aware counting.
#' @slot expressionPseudocount Non-negative pseudocount added to both WT and mutant mean counts before calculating descriptive log2 fold change. Defaults to 0.01 to avoid infinite values at zero counts; set to 0 to request the raw count ratio.
#' @slot exportAiccPlots Whether AICc span-search diagnostics should be written
#'   to \code{aicc_plots.pdf} in the output folder.
#' @slot refGenome An indexed \code{\link[Rsamtools]{FaFile}} generated
#'   internally from refFasta.
#' @rdname MmapprParam
#' 
NULL

# Candidate calling historically used Rsamtools' max_depth=250. Keep that stage-specific
# value separate from the configurable linkage maxPileupDepth setting. Because mutant BAMs
# are piled up separately and pooled afterward, this cap applies per mutant BAM; the
# linkage-stage depth setting remains independent.
.CANDIDATE_MAX_DEPTH <- 250L

setClass("MmapprParam",
         representation(
           wtFiles = "BamFileList",
           mutFiles = "BamFileList",
           refFasta = "character",
           refGenome = "FaFile",
           gtf = "character",
           outputFolder = "character",
           includeScaffolds = "logical",
           minDepth = "numeric",
           homozygoteCutoff = "numeric",
           minBaseQuality = "numeric",
           minMapQuality = "numeric",
           fileAggregation = "character",
           distancePower = "numeric",
           peakIntervalWidth = "numeric",
           loessOptResolution = "numeric",
           loessOptCutFactor = "numeric",
           maxPileupDepth = "numeric",
           candidateMinDepth = "numeric",
           candidateMinAltDepth = "numeric",
           candidateMinAltFreq = "numeric",
           candidateMaxWtAltFreq = "numeric",
           candidateMinDeltaAF = "numeric",
           candidatePoolMode = "character",
           candidateChunkSize = "numeric",
           peakCutoffSd = "numeric",
           peakCutoffMethod = "character",
           peakIntervalMethod = "character",
           peakResampleIterations = "numeric",
           randomSeed = "numeric",
           pairedEnd = "logical",
           ignoreStrand = "logical",
           expressionPseudocount = "numeric",
           exportAiccPlots = "logical"
         )
)



#' @title MmapprParam Constructor
#' 
#' @name mmapprParam
#' @description
#' Creates a new instance of a \code{\linkS4class{MmapprParam}} class object.
#'
#' @param wtFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the wild-type pool to be analyzed.
#' @param mutFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the mutant pool to be analyzed.
#' @param refFasta The path to the reference FASTA file.
#' @param gtf The path to a gtf-formatted annotation file. 
#' @param outputFolder Length-one character vector specifying where to save
#'   output, including a \code{\linkS4class{MmapprData}} stored as
#'   \code{mmappr_data.RDS}, \code{mmappr2.log}, a \code{.tsv} file
#'   for each peak chromosome containing candidate mutations, and PDF plots
#'   of both the entire genome and peak chromosomes. Defaults to an
#'   automatically generated \code{mmappr2_<timestamp>}.
#' @param includeScaffolds Logical indicating whether non-standard chromosomes
#'   should be included. FALSE (default) trims the seqnames of the reference to 
#'   include only standard chromosome names based on UCSC or Ensembl 
#'   conventions. In all cases, the mitochondrial chromosome is removed.
#' @param minDepth Length-one integer vector determining minimum depth
#'   required for a position to
#'   be considered in the analysis. Defaults to 20.
#' @param homozygoteCutoff Length-one numeric vector between \code{0} and
#'   \code{1} specifying threshold for throwing out base pairs on account
#'   of homozygosity. Positions with high major allele frequency in the
#'   wild-type pool are unlikely to exhibit polymorphism and are thus thrown
#'   out when they exceed this cutoff. Defaults to \code{0.95}.
#' @param minBaseQuality Length-one numeric vector indicating minimum base
#'   call quality to consider in analysis. Read positions with qualities
#'   below this score will be thrown out. Defaults to 20.
#' @param minMapQuality Length-one numeric vector indicating minimum read
#'   mapping quality to consider in analysis. Reads with qualities below
#'   this score will be thrown out. Defaults to 20.
#' @param fileAggregation A length-one character vector determining strategy
#'   for aggregating base calls when multiple wild-type or multiple mutant
#'   files are provided.
#'   When 'weighted', allele frequencies are weighted by read depth (pooled-read
#'   behavior). When 'simple', replicate allele frequencies are weighted equally
#'   while coverage is averaged across all supplied files.
#' @param distancePower Length-one numeric vector determining to what power
#'   Euclidean distance values are raised before fitting. Higher powers tend
#'   to increase high values and decrease low values, exaggerating the
#'   variation in the data. Default of 4.
#' @param peakIntervalWidth Length-one numeric vector between \code{0} and
#'   \code{1} specifying desired width of linkage region(s). The default value
#'   of \code{0.80} targets 80\% of the resampling-density mass. The exact
#'   interval construction is controlled by \code{peakIntervalMethod}; the default
#'   \code{hpd_span} conservatively spans all selected high-density modes.
#' @param loessOptResolution Length-one numeric vector between \code{0} and
#'   \code{1} specifying
#'   desired resolution for Loess fit optimization. The default of \code{0.001},
#'   for example, indicates that the span ultimately chosen will perform better
#'   than its neighbor values at \code{+-0.001}.
#' @param loessOptCutFactor Length-one numeric vector between \code{0} and
#'   \code{1} specifying how aggressively the Loess
#'   optimization algorithm proceeds. For example, with a default of \code{0.1}
#'   different spans at intervals of \code{0.001} would be evaluated after
#'   intervals of \code{0.01}.
#' @param maxPileupDepth Maximum reads considered at one linkage-stage pileup position. Defaults to 1000; candidate-stage pileup uses a separate explicit cap of 250 reads per mutant BAM.
#' @param candidateMinDepth Optional total mutant depth floor for a final candidate. Defaults to 1, which is effectively inert given the default two-ALT-read minimum and >80% ALT-frequency threshold.
#' @param candidateMinAltDepth Minimum mutant reads supporting the alternative allele; defaults to 2.
#' @param candidateMinAltFreq Strict lower bound on mutant alternative-allele frequency; candidates must be greater than this value. Defaults to 0.80.
#' @param candidateMaxWtAltFreq Optional maximum WT alternative-allele frequency. Default 1 disables this filter.
#' @param candidateMinDeltaAF Optional minimum mutant-minus-WT alternative-allele-frequency difference. Default 0 disables this filter.
#' @param candidatePoolMode Candidate-pileup memory strategy. \code{"auto"}
#'   (default) keeps the fast one-shot in-memory path when the requested peak is
#'   comfortably sized for available RAM and otherwise processes bounded genomic
#'   chunks. \code{"memory"} prefers the one-shot path and \code{"chunked"}
#'   always uses bounded chunks. Allocation failures are retried with smaller
#'   chunks and temporary on-disk staging.
#' @param candidateChunkSize Maximum candidate-pileup chunk width in bases.
#'   The default 0 chooses an adaptive size from available memory and BAM count.
#'   Positive values provide an explicit chunk width when chunking is used.
#' @param peakCutoffSd Multiplier applied to the selected peak-cutoff spread term. Defaults to 3.
#' @param peakCutoffMethod Initial peak-threshold method. \code{legacy_current}
#'   uses the compatibility mean-of-chromosome-means plus \code{sqrt(sum(var/n))} formula;
#'   \code{global_sd} uses a genome-wide median plus an ordinary SD multiple.
#' @param peakIntervalMethod Refined interval method. \code{hpd_span} (default)
#'   spans all high-density KDE grid points needed to reach the requested mass;
#'   \code{shortest_contiguous} chooses the narrowest single interval containing it.
#' @param peakResampleIterations Number of half-marker LOESS resamples used in peak refinement. Defaults to 1000.
#' @param randomSeed Base seed used to make peak refinement reproducible.
#' @param pairedEnd Logical indicating paired-end RNA-seq for the descriptive expression summary.
#' @param ignoreStrand Logical indicating whether strand is ignored in the descriptive expression summary. Defaults to FALSE for strand-aware \code{summarizeOverlaps} counting.
#' @param expressionPseudocount Non-negative pseudocount added to both WT and mutant mean counts before calculating descriptive log2 fold change. Defaults to 0.01 to avoid Inf/NaN values at zero counts; set to 0 to request the raw count ratio.
#' @param exportAiccPlots Logical. If TRUE, \code{outputMmapprData()} writes a
#'   multi-page \code{aicc_plots.pdf} showing every evaluated LOESS span and
#'   the selected optimum for each successfully fitted chromosome. Defaults to FALSE.
#' @param overwrite Logical. If TRUE, an existing non-empty output folder may be cleared explicitly.
#' @return A \code{MmapprParam} object.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#' 
NULL

mmapprParam <- function(wtFiles,
                        mutFiles,
                        refFasta,
                        gtf,
                        outputFolder = 'DEFAULT',
                        includeScaffolds = FALSE,
                        minDepth = 20,
                        homozygoteCutoff = 0.95,
                        minBaseQuality = 20,
                        minMapQuality = 20,
                        fileAggregation = c('simple', 'weighted'),
                        distancePower = 4,
                        peakIntervalWidth = 0.80,
                        loessOptResolution = 0.001,
                        loessOptCutFactor = 0.1,
                        maxPileupDepth = 1000,
                        candidateMinDepth = 1,
                        candidateMinAltDepth = 2,
                        candidateMinAltFreq = 0.80,
                        candidateMaxWtAltFreq = 1.00,
                        candidateMinDeltaAF = 0.00,
                        candidatePoolMode = c("auto", "memory", "chunked"),
                        candidateChunkSize = 0,
                        peakCutoffSd = 3,
                        peakCutoffMethod = c("legacy_current", "global_sd"),
                        peakIntervalMethod = c("hpd_span", "shortest_contiguous"),
                        peakResampleIterations = 1000,
                        randomSeed = 1,
                        pairedEnd = FALSE,
                        ignoreStrand = FALSE,
                        expressionPseudocount = 0.01,
                        exportAiccPlots = FALSE,
                        overwrite = FALSE) {

    wtPaths <- .asBamPaths(wtFiles)
    mutPaths <- .asBamPaths(mutFiles)

    # Validate control-flow parameters before using them in scalar `if`
    # statements. Character vectors of length >1 or NA used to fail with obscure
    # base-R condition errors instead of a useful MMAPPR2 message.
    if (!is.character(outputFolder) || length(outputFolder) != 1L || is.na(outputFolder) || !nzchar(outputFolder))
      stop("outputFolder must be one non-empty character path or 'DEFAULT'")
    if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
      stop("overwrite must be TRUE or FALSE")

    if (length(refFasta) != 1L || !file.exists(refFasta))
      stop("refFasta must name one existing FASTA file")
    if (length(gtf) != 1L || !file.exists(gtf))
      stop("gtf must name one existing GTF/GFF annotation file")
    if (!all(file.exists(wtPaths))) stop("One or more wild-type BAM files do not exist")
    if (!all(file.exists(mutPaths))) stop("One or more mutant BAM files do not exist")

    refFasta <- normalizePath(refFasta, mustWork = TRUE)
    gtf <- normalizePath(gtf, mustWork = TRUE)
    wtPaths <- normalizePath(wtPaths, mustWork = TRUE)
    mutPaths <- normalizePath(mutPaths, mustWork = TRUE)
    fileAggregation <- match.arg(fileAggregation)
    candidatePoolMode <- match.arg(candidatePoolMode)
    peakCutoffMethod <- match.arg(peakCutoffMethod)
    peakIntervalMethod <- match.arg(peakIntervalMethod)

    # Validate scalar analysis settings BEFORE creating/clearing the
    # output directory or indexing inputs. Invalid thresholds should fail without
    # leaving filesystem side effects behind.
    scalarErrors <- .validateScalarValues(
      includeScaffolds = includeScaffolds, minDepth = minDepth,
      homozygoteCutoff = homozygoteCutoff, minBaseQuality = minBaseQuality,
      minMapQuality = minMapQuality, fileAggregation = fileAggregation,
      distancePower = distancePower, peakIntervalWidth = peakIntervalWidth,
      loessOptResolution = loessOptResolution, loessOptCutFactor = loessOptCutFactor,
      maxPileupDepth = maxPileupDepth, candidateMinDepth = candidateMinDepth,
      candidateMinAltDepth = candidateMinAltDepth, candidateMinAltFreq = candidateMinAltFreq,
      candidateMaxWtAltFreq = candidateMaxWtAltFreq, candidateMinDeltaAF = candidateMinDeltaAF,
      candidatePoolMode = candidatePoolMode, candidateChunkSize = candidateChunkSize,
      peakCutoffSd = peakCutoffSd, peakCutoffMethod = peakCutoffMethod,
      peakIntervalMethod = peakIntervalMethod, peakResampleIterations = peakResampleIterations,
      randomSeed = randomSeed, pairedEnd = pairedEnd, ignoreStrand = ignoreStrand,
      expressionPseudocount = expressionPseudocount,
      exportAiccPlots = exportAiccPlots, nMutFiles = length(mutPaths)
    )
    if (length(scalarErrors)) stop(paste(scalarErrors, collapse = "\n  "))

    if (anyDuplicated(wtPaths)) warning("Duplicate wild-type BAM path(s) supplied; they will be counted repeatedly")
    if (anyDuplicated(mutPaths)) warning("Duplicate mutant BAM path(s) supplied; they will be counted repeatedly")

    wtFiles <- .indexBamFileList(wtPaths, NULL)
    mutFiles <- .indexBamFileList(mutPaths, NULL)
    bamChecks <- list(.validBamFiles(wtFiles, deep = TRUE),
                      .validBamFiles(mutFiles, deep = TRUE))
    failedBamChecks <- bamChecks[!vapply(bamChecks, isTRUE, logical(1))]
    if (length(failedBamChecks))
      stop(paste(unlist(failedBamChecks, use.names = FALSE), collapse = "\n  "))

    # Reuse a current FASTA index, but rebuild a missing, empty, or
    # stale sidecar before any output-folder side effects.
    fastaIndex <- paste0(refFasta, ".fai")
    if (!.indexIsCurrent(refFasta, fastaIndex)) Rsamtools::indexFa(refFasta)
    refGenome <- Rsamtools::FaFile(refFasta)
    fastaCheck <- tryCatch({
      si <- GenomeInfoDb::seqinfo(refGenome)
      if (length(si) == 0L) stop("reference FASTA has no indexed sequences")
      TRUE
    }, error = function(e) e$message)
    if (!isTRUE(fastaCheck)) stop("Reference FASTA/index is unreadable: ", fastaCheck)

    .preflightInputResourcesFrozen(wtFiles, mutFiles, refGenome, gtf)

    if (outputFolder == 'DEFAULT') outputFolder <- .defaultOutputFolder()
    # Never delete an existing result directory because of an interactive
    # prompt. Batch jobs cannot answer safely, and accidental deletion is worse
    # than requiring an explicit overwrite=TRUE.
    outputFolder <- .prepareOutputFolder(outputFolder, overwrite = overwrite)


    param <- new("MmapprParam",
                 wtFiles = wtFiles,
                 mutFiles = mutFiles,
                 refFasta = refFasta,
                 refGenome = refGenome,
                 gtf = gtf,
                 outputFolder = outputFolder,
                 includeScaffolds = includeScaffolds,
                 minDepth = minDepth,
                 homozygoteCutoff = homozygoteCutoff,
                 minBaseQuality = minBaseQuality,
                 minMapQuality = minMapQuality,
                 fileAggregation = fileAggregation,
                 distancePower = distancePower,
                 peakIntervalWidth = peakIntervalWidth,
                 loessOptResolution = loessOptResolution,
                 loessOptCutFactor = loessOptCutFactor,
                 maxPileupDepth = maxPileupDepth,
                 candidateMinDepth = candidateMinDepth,
                 candidateMinAltDepth = candidateMinAltDepth,
                 candidateMinAltFreq = candidateMinAltFreq,
                 candidateMaxWtAltFreq = candidateMaxWtAltFreq,
                 candidateMinDeltaAF = candidateMinDeltaAF,
                 candidatePoolMode = candidatePoolMode,
                 candidateChunkSize = candidateChunkSize,
                 peakCutoffSd = peakCutoffSd,
                 peakCutoffMethod = peakCutoffMethod,
                 peakIntervalMethod = peakIntervalMethod,
                 peakResampleIterations = peakResampleIterations,
                 randomSeed = randomSeed,
                 pairedEnd = pairedEnd,
                 ignoreStrand = ignoreStrand,
                 expressionPseudocount = expressionPseudocount,
                 exportAiccPlots = exportAiccPlots)

    # Deep BAM/index checks already succeeded before the output directory was
    # touched. Normal object validity here stays lightweight.
    validity <- .validMmapprParam(param, deepBam = FALSE)
    if (isTRUE(validity)) param else stop(paste(validity, collapse = '\n  '))
}



### VALIDITY FUNCTIONS

.scalarNumeric <- function(x, lower = -Inf, upper = Inf, lowerOpen = FALSE,
                           upperOpen = FALSE, integerish = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) return(FALSE)
    if (lowerOpen) { if (x <= lower) return(FALSE) } else if (x < lower) return(FALSE)
    if (upperOpen) { if (x >= upper) return(FALSE) } else if (x > upper) return(FALSE)
    if (integerish && x != floor(x)) return(FALSE)
    TRUE
}

.validateScalarValues <- function(includeScaffolds, minDepth, homozygoteCutoff,
                                  minBaseQuality, minMapQuality, fileAggregation,
                                  distancePower, peakIntervalWidth, loessOptResolution,
                                  loessOptCutFactor, maxPileupDepth, candidateMinDepth,
                                  candidateMinAltDepth, candidateMinAltFreq,
                                  candidateMaxWtAltFreq, candidateMinDeltaAF,
                                  candidatePoolMode, candidateChunkSize,
                                  peakCutoffSd, peakCutoffMethod, peakIntervalMethod,
                                  peakResampleIterations, randomSeed, pairedEnd,
                                  ignoreStrand, expressionPseudocount,
                                  exportAiccPlots, nMutFiles = 1L) {
    errors <- character()
    add <- function(ok, msg) if (!isTRUE(ok)) errors <<- c(errors, msg)
    add(is.logical(includeScaffolds) && length(includeScaffolds) == 1L && !is.na(includeScaffolds), "includeScaffolds must be TRUE or FALSE")
    add(.scalarNumeric(minDepth, 1, integerish = TRUE), "minDepth must be a positive integer")
    add(.scalarNumeric(homozygoteCutoff, 0, 1), "homozygoteCutoff must be between 0 and 1")
    add(.scalarNumeric(minBaseQuality, 0, integerish = TRUE), "minBaseQuality must be a non-negative integer")
    add(.scalarNumeric(minMapQuality, 0, integerish = TRUE), "minMapQuality must be a non-negative integer")
    add(is.character(fileAggregation) && length(fileAggregation) == 1L && fileAggregation %in% c("simple", "weighted"), "fileAggregation must be 'simple' or 'weighted'")
    add(.scalarNumeric(distancePower, 0, lowerOpen = TRUE), "distancePower must be > 0")
    add(.scalarNumeric(peakIntervalWidth, 0, 1, lowerOpen = TRUE), "peakIntervalWidth must be in (0, 1]")
    add(.scalarNumeric(loessOptResolution, 0, 1, lowerOpen = TRUE), "loessOptResolution must be in (0, 1]")
    add(.scalarNumeric(loessOptCutFactor, 0, 1, lowerOpen = TRUE, upperOpen = TRUE), "loessOptCutFactor must be in (0, 1)")
    add(.scalarNumeric(maxPileupDepth, 1, integerish = TRUE), "maxPileupDepth must be a positive integer")
    add(.scalarNumeric(candidateMinDepth, 1, integerish = TRUE), "candidateMinDepth must be a positive integer")
    add(.scalarNumeric(candidateMinAltDepth, 1, integerish = TRUE), "candidateMinAltDepth must be a positive integer")
    add(.scalarNumeric(candidateMinAltFreq, 0, 1, upperOpen = TRUE), "candidateMinAltFreq must be in [0, 1)")
    add(.scalarNumeric(candidateMaxWtAltFreq, 0, 1), "candidateMaxWtAltFreq must be between 0 and 1")
    add(.scalarNumeric(candidateMinDeltaAF, 0, 1), "candidateMinDeltaAF must be between 0 and 1")
    add(is.character(candidatePoolMode) && length(candidatePoolMode) == 1L &&
          !is.na(candidatePoolMode) &&
          candidatePoolMode %in% c("auto", "memory", "chunked"),
        "candidatePoolMode must be 'auto', 'memory', or 'chunked'")
    add(.scalarNumeric(candidateChunkSize, 0, integerish = TRUE),
        "candidateChunkSize must be 0 or a positive integer number of bases")
    add(.scalarNumeric(peakCutoffSd, 0), "peakCutoffSd must be non-negative")
    add(is.character(peakCutoffMethod) && length(peakCutoffMethod) == 1L && peakCutoffMethod %in% c("legacy_current", "global_sd"), "peakCutoffMethod must be 'legacy_current' or 'global_sd'")
    add(is.character(peakIntervalMethod) && length(peakIntervalMethod) == 1L && peakIntervalMethod %in% c("hpd_span", "shortest_contiguous"), "peakIntervalMethod must be 'hpd_span' or 'shortest_contiguous'")
    add(.scalarNumeric(peakResampleIterations, 10, integerish = TRUE), "peakResampleIterations must be an integer >= 10")
    add(.scalarNumeric(randomSeed, 0, .Machine$integer.max, integerish = TRUE), "randomSeed must be a non-negative integer")
    add(is.logical(pairedEnd) && length(pairedEnd) == 1L && !is.na(pairedEnd), "pairedEnd must be TRUE or FALSE")
    add(is.logical(ignoreStrand) && length(ignoreStrand) == 1L && !is.na(ignoreStrand), "ignoreStrand must be TRUE or FALSE")
    add(.scalarNumeric(expressionPseudocount, 0, lowerOpen = FALSE), "expressionPseudocount must be >= 0")
    add(is.logical(exportAiccPlots) && length(exportAiccPlots) == 1L && !is.na(exportAiccPlots),
        "exportAiccPlots must be TRUE or FALSE")

    # Catch parameter combinations that guarantee an empty result.
    if (.scalarNumeric(minDepth, 1, integerish = TRUE) && .scalarNumeric(maxPileupDepth, 1, integerish = TRUE))
        add(minDepth <= maxPileupDepth, "minDepth cannot exceed maxPileupDepth (no mapping position could pass)")
    # Candidate pileup uses an explicit Rsamtools cap of 250 reads per mutant BAM.
    # maxPileupDepth controls the linkage stage only.
    if (.scalarNumeric(candidateMinDepth, 1, integerish = TRUE) &&
        .scalarNumeric(nMutFiles, 1, integerish = TRUE))
        add(candidateMinDepth <= .CANDIDATE_MAX_DEPTH * nMutFiles, "candidateMinDepth exceeds the candidate-pileup depth capacity (250 reads per mutant BAM)")
    if (.scalarNumeric(candidateMinAltDepth, 1, integerish = TRUE) &&
        .scalarNumeric(nMutFiles, 1, integerish = TRUE))
        add(candidateMinAltDepth <= .CANDIDATE_MAX_DEPTH * nMutFiles, "candidateMinAltDepth exceeds the candidate-pileup depth capacity (250 reads per mutant BAM)")
    errors
}

.validMmapprParam <- function(param, deepBam = FALSE) {
    errors <- character()
    add <- function(ok, msg) if (!isTRUE(ok)) errors <<- c(errors, msg)

    fasta_ok <- .validFastaFile(refFasta(param))
    if (!isTRUE(fasta_ok)) errors <- c(errors, fasta_ok)
    wt_ok <- .validBamFiles(wtFiles(param), deep = deepBam); if (!isTRUE(wt_ok)) errors <- c(errors, wt_ok)
    mut_ok <- .validBamFiles(mutFiles(param), deep = deepBam); if (!isTRUE(mut_ok)) errors <- c(errors, mut_ok)
    add(length(gtf(param)) == 1L && file.exists(gtf(param)), "gtf must name one existing annotation file")
    add(length(outputFolder(param)) == 1L && nzchar(outputFolder(param)) && dir.exists(outputFolder(param)),
        "outputFolder must be one existing directory")
    # refFasta and refGenome represent the same resource and must remain in
    # sync. This catches hand-constructed/corrupted S4 objects in addition to the
    # corrected refFasta<- setter.
    refGenomePath <- tryCatch(BiocGenerics::path(param@refGenome), error = function(e) NA_character_)
    add(length(refGenomePath) == 1L && !is.na(refGenomePath) &&
          normalizePath(refGenomePath, mustWork = FALSE) == normalizePath(refFasta(param), mustWork = FALSE),
        "refGenome does not point to refFasta")
    scalarErrors <- .validateScalarValues(
      includeScaffolds = includeScaffolds(param), minDepth = minDepth(param),
      homozygoteCutoff = homozygoteCutoff(param), minBaseQuality = minBaseQuality(param),
      minMapQuality = minMapQuality(param), fileAggregation = fileAggregation(param),
      distancePower = distancePower(param), peakIntervalWidth = peakIntervalWidth(param),
      loessOptResolution = loessOptResolution(param), loessOptCutFactor = loessOptCutFactor(param),
      maxPileupDepth = maxPileupDepth(param), candidateMinDepth = candidateMinDepth(param),
      candidateMinAltDepth = candidateMinAltDepth(param), candidateMinAltFreq = candidateMinAltFreq(param),
      candidateMaxWtAltFreq = candidateMaxWtAltFreq(param), candidateMinDeltaAF = candidateMinDeltaAF(param),
      candidatePoolMode = candidatePoolMode(param), candidateChunkSize = candidateChunkSize(param),
      peakCutoffSd = peakCutoffSd(param), peakCutoffMethod = peakCutoffMethod(param),
      peakIntervalMethod = peakIntervalMethod(param), peakResampleIterations = peakResampleIterations(param),
      randomSeed = randomSeed(param), pairedEnd = pairedEnd(param), ignoreStrand = ignoreStrand(param),
      expressionPseudocount = expressionPseudocount(param),
      exportAiccPlots = exportAiccPlots(param), nMutFiles = max(1L, length(mutFiles(param)))
    )
    if (length(scalarErrors)) errors <- c(errors, scalarErrors)

    if (length(errors) == 0L) TRUE else errors
}

# Register lightweight structural/scalar validity for direct S4 object
# construction. Deep BAM/index and cross-resource build concordance remain public
# constructor/resource-setter checks so ordinary scalar validation stays cheap.
setValidity("MmapprParam", function(object) .validMmapprParam(object, deepBam = FALSE))

.validFastaFile <- function(filepath) {
    # Validate path shape/content before touching the filesystem so malformed
    # inputs produce a stable diagnostic instead of an incidental `if()` error.
    if (!is.character(filepath) || length(filepath) != 1L ||
        is.na(filepath) || !nzchar(filepath))
        return("refFasta must name one existing FASTA file")
    if (!file.exists(filepath))
        return(paste(filepath, "does not exist"))
    TRUE
}

.validBamFiles <- function(files, deep = FALSE) {
    errors <- character()
    if (!is(files, 'BamFileList')) return("Input is not a BamFileList object")
    if (length(files) == 0L) return("At least one BAM file is required")
    for (i in seq_along(files)) {
        bam <- files[[i]]
        p <- BiocGenerics::path(bam)
        if (!file.exists(p)) {
            errors <- c(errors, paste0(p, " does not exist"))
            next
        }
        if (isTRUE(deep)) {
            # File existence alone does not prove that a .bam is a
            # readable BAM or that its .bai is usable. A stray/corrupt BAM+BAI pair
            # otherwise survives construction and fails much later during regional
            # pileup. Deep checks are reserved for construction/BAM replacement so
            # changing an unrelated scalar parameter does not repeatedly touch disk.
            headerOk <- tryCatch({ Rsamtools::scanBamHeader(bam); TRUE },
                                 error = function(e) e$message)
            if (!isTRUE(headerOk))
                errors <- c(errors, paste0(p, " has an unreadable BAM header: ", headerOk))
            indexOk <- tryCatch({ Rsamtools::idxstatsBam(bam); TRUE },
                                error = function(e) e$message)
            if (!isTRUE(indexOk))
                errors <- c(errors, paste0(p, " has an unreadable BAM index: ", indexOk))
        }
    }
    if (length(errors) == 0L) TRUE else errors
}

.asBamPaths <- function(x) {
    if (is(x, "BamFileList")) {
        if (length(x) == 0L) stop("At least one BAM file is required")
        return(vapply(x, BiocGenerics::path, character(1)))
    }
    if (is(x, "BamFile")) return(BiocGenerics::path(x))
    if (is.character(x) && length(x) > 0L) return(x)
    stop("BAM inputs must be character paths, a BamFile, or a BamFileList")
}

.indexIsCurrent <- function(dataPath, indexPath) {
    if (!file.exists(dataPath) || !file.exists(indexPath)) return(FALSE)
    di <- file.info(dataPath)
    ii <- file.info(indexPath)
    if (is.na(ii$size) || ii$size <= 0L || is.na(di$mtime) || is.na(ii$mtime))
        return(FALSE)
    ii$mtime >= di$mtime
}

.indexBamFileList <- function(bfl, oF = NULL) {
    emit <- function(msg) {
        if (is.null(oF)) message(msg) else .messageAndLog(msg, oF)
    }
    indexed_bfl <- lapply(bfl, function(bam_file) {
        candidates <- paste0(bam_file, ".bai")
        if (grepl("\\.bam$", bam_file, ignore.case = TRUE)) {
            candidates <- c(candidates, sub("\\.bam$", ".bai", bam_file, ignore.case = TRUE))
        }
        candidates <- unique(candidates)
        current <- candidates[vapply(candidates, function(idx)
            .indexIsCurrent(bam_file, idx), logical(1))]
        if (length(current) == 0L) {
            existing <- candidates[file.exists(candidates)]
            reason <- if (length(existing)) "Index missing/stale. Re-indexing BAM now." else "No index found. Indexing BAM now."
            emit(paste(bam_file, "--", reason))
            bam_index <- Rsamtools::indexBam(bam_file)
        } else {
            emit(paste(bam_file, "-- Current index found. Skipping index step."))
            bam_index <- current[[1]]
        }
        Rsamtools::BamFile(bam_file, index = bam_index)
    })
    # BamFileList() is a dots-style constructor. Splice the BamFile
    # objects into `...` explicitly rather than passing the ordinary list as one
    # element, which is version-sensitive in Rsamtools/S4Vectors.
    do.call(Rsamtools::BamFileList, indexed_bfl)
}

.choose_target_style <- function(si) {
    # Unknown/custom assemblies can make seqlevelsStyle() unable to
    # infer a naming convention. Treat that as "no style hint" and fall back to
    # exact sequence-name matching instead of failing before the real comparison.
    st <- tryCatch(unique(unlist(GenomeInfoDb::seqlevelsStyle(si))),
                   error = function(e) character())
    st <- st[!is.na(st) & nzchar(st)]
    if (length(st)) st[[1]] else NA_character_
}

.faSeqinfo <- function(fa) {
    # FaFile has a seqinfo() method backed directly by the FASTA index;
    # use it rather than reconstructing Seqinfo indirectly from scanFaIndex().
    GenomeInfoDb::seqinfo(fa)
}

.disconnectTxDb <- function(txdb) {
    con <- tryCatch(BiocGenerics::dbconn(txdb), error = function(e) NULL)
    if (!is.null(con) && tryCatch(DBI::dbIsValid(con), error = function(e) FALSE))
        try(DBI::dbDisconnect(con), silent = TRUE)
    invisible(NULL)
}

.makeTxDbFromAnnotation <- function(path) {
    withCallingHandlers(
        txdbmaker::makeTxDbFromGFF(file = path, format = "auto"),
        warning = function(w) {
            if (identical(conditionMessage(w),
                          "genome version information is not available for this TxDb object"))
                invokeRestart("muffleWarning")
        }
    )
}

.buildTxDb <- function(param) .makeTxDbFromAnnotation(gtf(param))


setMethod("show", "MmapprParam", function(object) {
    margin <- "   "
    cat("MmapprParam object with following values:\n")
    cat("Reference fasta file:\n", sep = "")
    cat(paste0(margin, object@refFasta, '\n'))
    cat("GTF file:\n", sep = "")
    cat(paste0(margin, object@gtf, '\n'))
    cat("wtFiles:\n", sep = "")
    .customPrint(object@wtFiles, margin)
    cat("mutFiles:\n", sep = "")
    .customPrint(object@mutFiles, margin)

    cat("Other parameters:\n")
    # Print every tunable scalar parameter, including newly exposed
    # candidate/resampling settings, instead of relying on brittle numeric slots.
    slotNames <- setdiff(methods::slotNames("MmapprParam"),
                         c("wtFiles", "mutFiles", "refFasta", "refGenome", "gtf"))
    slotValues <- vapply(slotNames,
                         function(name, object) {
                           as.character(slot(object, name)[1])
                         },
                         FUN.VALUE = character(1), object)

    names(slotValues) <- slotNames
    print(slotValues, quote = FALSE)
})


.customPrint <- function(obj, margin = "  ", lineMax = getOption("max.print")) {
  lines <- capture.output(obj)
  lines <- strsplit(lines, split = "\n")
  lines <- vapply(lines, function(x) paste0(margin, x), character(1))
  if (lineMax > length(lines)) lineMax = length(lines)
  cat(lines[seq_len(lineMax)], sep = "\n")
}


#' MmapprParam Getters and Setters
#'
#' Access and assign slots of \code{\link{MmapprParam}} object.
#'
#' @name MmapprParam-functions
#' @aliases
#'   wtFiles wtFiles<-
#'   mutFiles mutFiles<-
#'   refFasta refFasta<-
#'   gtf gtf<-
#'   outputFolder outputFolder<-
#'   includeScaffolds includeScaffolds<-
#'   minDepth minDepth<-
#'   homozygoteCutoff homozygoteCutoff<-
#'   minBaseQuality minBaseQuality<-
#'   minMapQuality minMapQuality<-
#'   fileAggregation fileAggregation<-
#'   distancePower distancePower<-
#'   peakIntervalWidth peakIntervalWidth<-
#'   loessOptResolution loessOptResolution<-
#'   loessOptCutFactor loessOptCutFactor<-
#'   maxPileupDepth maxPileupDepth<-
#'   candidateMinDepth candidateMinDepth<-
#'   candidateMinAltDepth candidateMinAltDepth<-
#'   candidateMinAltFreq candidateMinAltFreq<-
#'   candidateMaxWtAltFreq candidateMaxWtAltFreq<-
#'   candidateMinDeltaAF candidateMinDeltaAF<-
#'   candidatePoolMode candidatePoolMode<-
#'   candidateChunkSize candidateChunkSize<-
#'   peakCutoffSd peakCutoffSd<-
#'   peakCutoffMethod peakCutoffMethod<-
#'   peakIntervalMethod peakIntervalMethod<-
#'   peakResampleIterations peakResampleIterations<-
#'   randomSeed randomSeed<-
#'   pairedEnd pairedEnd<-
#'   ignoreStrand ignoreStrand<-
#'   expressionPseudocount expressionPseudocount<-
#'   exportAiccPlots exportAiccPlots<-
#'
#' @param obj Desired \code{\link{MmapprParam}} object.
#' @param value Value to replace desired attribute.
#'
#' @return The desired \code{\link{MmapprParam}} attribute.
#'
#' @seealso \code{\link{MmapprParam}}
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                          mutFiles = MMAPPR2data::exampleMutBam(),
#'                          refFasta = MMAPPR2data::goldenFasta(),
#'                          gtf = MMAPPR2data::gtf())
#'
#'     outputFolder(mmappr_param) <- 'mmappr2_test_1'
#'     minBaseQuality(mmappr_param) <- 25
#'     candidateMinAltFreq(mmappr_param)
#' }
NULL


### GETTERS
#' @rdname MmapprParam-functions
#' @export
setMethod("wtFiles", "MmapprParam", function(obj) obj@wtFiles)
#' @rdname MmapprParam-functions
#' @export
setMethod("mutFiles", "MmapprParam", function(obj) obj@mutFiles)
#' @rdname MmapprParam-functions
#' @export
setMethod("refFasta", "MmapprParam", function(obj) obj@refFasta)
#' @rdname MmapprParam-functions
#' @export
setMethod("gtf", "MmapprParam", function(obj) obj@gtf)
#' @rdname MmapprParam-functions
#' @export
setMethod("outputFolder", "MmapprParam", function(obj) obj@outputFolder)
#' @rdname MmapprParam-functions
#' @export
setMethod("includeScaffolds", "MmapprParam", function(obj) obj@includeScaffolds)
#' @rdname MmapprParam-functions
#' @export
setMethod("minDepth", "MmapprParam", function(obj) obj@minDepth)
#' @rdname MmapprParam-functions
#' @export
setMethod("homozygoteCutoff", "MmapprParam", function(obj) obj@homozygoteCutoff)
#' @rdname MmapprParam-functions
#' @export
setMethod("minBaseQuality", "MmapprParam", function(obj) obj@minBaseQuality)
#' @rdname MmapprParam-functions
#' @export
setMethod("minMapQuality", "MmapprParam", function(obj) obj@minMapQuality)
#' @rdname MmapprParam-functions
#' @export
setMethod("fileAggregation", "MmapprParam", function(obj) obj@fileAggregation)
#' @rdname MmapprParam-functions
#' @export
setMethod("distancePower", "MmapprParam", function(obj) obj@distancePower)
#' @rdname MmapprParam-functions
#' @export
setMethod("peakIntervalWidth", "MmapprParam", function(obj) obj@peakIntervalWidth)
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptResolution", "MmapprParam", function(obj) obj@loessOptResolution)
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptCutFactor", "MmapprParam", function(obj) obj@loessOptCutFactor)
#' @rdname MmapprParam-functions
#' @export
setMethod("maxPileupDepth", "MmapprParam", function(obj) obj@maxPileupDepth)
#' @rdname MmapprParam-functions
#' @export
setMethod("candidateMinDepth", "MmapprParam", function(obj) obj@candidateMinDepth)
#' @rdname MmapprParam-functions
#' @export
setMethod("candidateMinAltDepth", "MmapprParam", function(obj) obj@candidateMinAltDepth)
#' @rdname MmapprParam-functions
#' @export
setMethod("candidateMinAltFreq", "MmapprParam", function(obj) obj@candidateMinAltFreq)
#' @rdname MmapprParam-functions
#' @export
setMethod("candidateMaxWtAltFreq", "MmapprParam", function(obj) obj@candidateMaxWtAltFreq)
#' @rdname MmapprParam-functions
#' @export
setMethod("candidateMinDeltaAF", "MmapprParam", function(obj) obj@candidateMinDeltaAF)
#' @rdname MmapprParam-functions
#' @export
setMethod("candidatePoolMode", "MmapprParam", function(obj) obj@candidatePoolMode)
#' @rdname MmapprParam-functions
#' @export
setMethod("candidateChunkSize", "MmapprParam", function(obj) obj@candidateChunkSize)
#' @rdname MmapprParam-functions
#' @export
setMethod("peakCutoffSd", "MmapprParam", function(obj) obj@peakCutoffSd)
#' @rdname MmapprParam-functions
#' @export
setMethod("peakCutoffMethod", "MmapprParam", function(obj) obj@peakCutoffMethod)
#' @rdname MmapprParam-functions
#' @export
setMethod("peakIntervalMethod", "MmapprParam", function(obj) obj@peakIntervalMethod)
#' @rdname MmapprParam-functions
#' @export
setMethod("peakResampleIterations", "MmapprParam", function(obj) obj@peakResampleIterations)
#' @rdname MmapprParam-functions
#' @export
setMethod("randomSeed", "MmapprParam", function(obj) obj@randomSeed)
#' @rdname MmapprParam-functions
#' @export
setMethod("pairedEnd", "MmapprParam", function(obj) obj@pairedEnd)
#' @rdname MmapprParam-functions
#' @export
setMethod("ignoreStrand", "MmapprParam", function(obj) obj@ignoreStrand)
#' @rdname MmapprParam-functions
#' @export
setMethod("expressionPseudocount", "MmapprParam", function(obj) obj@expressionPseudocount)
#' @rdname MmapprParam-functions
#' @export
setMethod("exportAiccPlots", "MmapprParam", function(obj) obj@exportAiccPlots)

### SETTERS

.validateAfterSet <- function(obj, deepBam = FALSE) {
    v <- .validMmapprParam(obj, deepBam = deepBam)
    if (isTRUE(v)) obj else stop(paste(v, collapse = "\n  "))
}

.preflightAfterResourceSet <- function(obj) {
    .preflightInputResourcesFrozen(obj@wtFiles, obj@mutFiles, obj@refGenome, obj@gtf)
    obj
}

setMethod("wtFiles<-", "MmapprParam", function(obj, value) {
    obj@wtFiles <- .indexBamFileList(normalizePath(.asBamPaths(value), mustWork = TRUE),
                                     outputFolder(obj))
    obj <- .validateAfterSet(obj, deepBam = TRUE)
    .preflightAfterResourceSet(obj)
})

setMethod("mutFiles<-", "MmapprParam", function(obj, value) {
    obj@mutFiles <- .indexBamFileList(normalizePath(.asBamPaths(value), mustWork = TRUE),
                                      outputFolder(obj))
    obj <- .validateAfterSet(obj, deepBam = TRUE)
    .preflightAfterResourceSet(obj)
})

setMethod("refFasta<-", "MmapprParam", function(obj, value) {
    if (length(value) != 1L || !file.exists(value)) stop("refFasta must exist")
    value <- normalizePath(value, mustWork = TRUE)
    if (!.indexIsCurrent(value, paste0(value, ".fai"))) Rsamtools::indexFa(value)
    obj@refFasta <- value
    obj@refGenome <- Rsamtools::FaFile(value)
    obj <- .validateAfterSet(obj)
    .preflightAfterResourceSet(obj)
})

setMethod("gtf<-", "MmapprParam", function(obj, value) {
    if (length(value) != 1L || !file.exists(value)) stop("gtf must exist")
    obj@gtf <- normalizePath(value, mustWork = TRUE)
    obj <- .validateAfterSet(obj)
    .preflightAfterResourceSet(obj)
})

setMethod("outputFolder<-", "MmapprParam", function(obj, value) {
    if (!is.character(value) || length(value) != 1L || !nzchar(value))
        stop("outputFolder must be one non-empty path")
    if (!dir.exists(value)) dir.create(value, recursive = TRUE, showWarnings = FALSE)
    obj@outputFolder <- normalizePath(value, mustWork = TRUE)
    .validateAfterSet(obj)
})

.setScalarSlot <- function(obj, slotName, value) {
    slot(obj, slotName) <- value
    .validateAfterSet(obj)
}

setMethod("includeScaffolds<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "includeScaffolds", value))
setMethod("minDepth<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "minDepth", value))
setMethod("homozygoteCutoff<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "homozygoteCutoff", value))
setMethod("minBaseQuality<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "minBaseQuality", value))
setMethod("minMapQuality<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "minMapQuality", value))
setMethod("fileAggregation<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "fileAggregation", match.arg(value, c("simple", "weighted"))))
setMethod("distancePower<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "distancePower", value))
setMethod("peakIntervalWidth<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "peakIntervalWidth", value))
setMethod("loessOptResolution<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "loessOptResolution", value))
setMethod("loessOptCutFactor<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "loessOptCutFactor", value))
setMethod("maxPileupDepth<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "maxPileupDepth", value))
setMethod("candidateMinDepth<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "candidateMinDepth", value))
setMethod("candidateMinAltDepth<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "candidateMinAltDepth", value))
setMethod("candidateMinAltFreq<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "candidateMinAltFreq", value))
setMethod("candidateMaxWtAltFreq<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "candidateMaxWtAltFreq", value))
setMethod("candidateMinDeltaAF<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "candidateMinDeltaAF", value))
setMethod("candidatePoolMode<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "candidatePoolMode", match.arg(value, c("auto", "memory", "chunked"))))
setMethod("candidateChunkSize<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "candidateChunkSize", value))
setMethod("peakCutoffSd<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "peakCutoffSd", value))
setMethod("peakCutoffMethod<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "peakCutoffMethod", match.arg(value, c("legacy_current", "global_sd"))))
setMethod("peakIntervalMethod<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "peakIntervalMethod", match.arg(value, c("hpd_span", "shortest_contiguous"))))
setMethod("peakResampleIterations<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "peakResampleIterations", value))
setMethod("randomSeed<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "randomSeed", value))
setMethod("pairedEnd<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "pairedEnd", value))
setMethod("ignoreStrand<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "ignoreStrand", value))
setMethod("expressionPseudocount<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "expressionPseudocount", value))
setMethod("exportAiccPlots<-", "MmapprParam", function(obj, value) .setScalarSlot(obj, "exportAiccPlots", value))
