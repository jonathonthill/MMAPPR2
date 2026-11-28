% Manually synchronized with R/param.R.
\name{mmapprParam}
\alias{mmapprParam}
\title{Construct MMAPPR2 analysis parameters}
\usage{
mmapprParam(
  wtFiles,
  mutFiles,
  refFasta,
  gtf,
  outputFolder = "DEFAULT",
  includeScaffolds = FALSE,
  minDepth = 20,
  homozygoteCutoff = 0.95,
  minBaseQuality = 20,
  minMapQuality = 20,
  fileAggregation = c("simple", "weighted"),
  distancePower = 4,
  peakIntervalWidth = 0.8,
  loessOptResolution = 0.001,
  loessOptCutFactor = 0.1,
  maxPileupDepth = 1000,
  candidateMinDepth = 1,
  candidateMinAltDepth = 2,
  candidateMinAltFreq = 0.8,
  candidateMaxWtAltFreq = 1,
  candidateMinDeltaAF = 0,
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
  overwrite = FALSE
)
}
\arguments{
\item{wtFiles}{Character paths, \code{BamFile}, or \code{BamFileList} for the phenotypically wild-type F2 pool. At least one BAM is required.}
\item{mutFiles}{Character paths, \code{BamFile}, or \code{BamFileList} for the mutant F2 pool. At least one BAM is required.}
\item{refFasta}{Path to the reference-genome FASTA used for alignment.}
\item{gtf}{Path to a GTF or GFF3 annotation readable by Bioconductor.}
\item{outputFolder}{Output directory. \code{"DEFAULT"} creates a collision-resistant timestamped directory.}
\item{includeScaffolds}{Whether non-standard chromosomes/contigs should be retained. Mitochondrial sequence is excluded.}
\item{minDepth}{Minimum average pool depth used for linkage markers. Default 20.}
\item{homozygoteCutoff}{Discard mapping positions whose largest WT A/C/G/T frequency exceeds this value. Default 0.95.}
\item{minBaseQuality}{Minimum base quality used in pileups. Default 20.}
\item{minMapQuality}{Minimum read mapping quality used in pileups. Default 20.}
\item{fileAggregation}{\code{"simple"} gives equal weight to replicate allele-frequency estimates; \code{"weighted"} pools allele evidence by read depth.}
\item{distancePower}{Power applied to Euclidean distance before LOESS. Default 4.}
\item{peakIntervalWidth}{Target fraction of the resampling-density mass retained during peak refinement. Default 0.80. Interpretation depends on \code{peakIntervalMethod}.}
\item{loessOptResolution}{Target resolution of LOESS-span optimization. Default 0.001.}
\item{loessOptCutFactor}{Coarse-to-fine LOESS optimization refinement factor. Default 0.1.}
\item{maxPileupDepth}{Maximum reads considered per BAM at one linkage-stage pileup position. Default 1000; final candidate pileup retains the historical 250-read cap.}
\item{candidateMinDepth}{Optional pooled mutant depth floor for a candidate SNV. Default 1 is effectively inert given the two-ALT-read rule.}
\item{candidateMinAltDepth}{Minimum pooled mutant reads supporting the alternative allele. Default 2.}
\item{candidateMinAltFreq}{Minimum mutant alternative-allele frequency for a candidate. Default 0.80.}
\item{candidateMaxWtAltFreq}{Optional maximum WT alternative-allele frequency. Default 1 disables this hard filter while WT evidence is still reported.}
\item{candidateMinDeltaAF}{Optional minimum mutant-minus-WT alternative-allele-frequency difference. Default 0 disables this hard filter.}
\item{candidatePoolMode}{Candidate-pileup memory strategy. \code{"auto"} (default) keeps the one-shot in-memory path when the interval is comfortably sized for available memory and otherwise uses bounded genomic chunks. \code{"memory"} prefers the one-shot path; \code{"chunked"} always uses bounded chunks. Allocation failures are retried with smaller chunks; the emergency mutant-candidate retry stages already-filtered chunk results in the R temporary directory and removes them automatically.}
\item{candidateChunkSize}{Maximum candidate-pileup chunk width in bases. Default 0 chooses an adaptive size from available memory and the number of mutant BAM files.}
\item{peakCutoffSd}{Multiplier applied to the variability term used by the selected initial peak-cutoff method. Default 3.}
\item{peakCutoffMethod}{\code{"legacy_current"} (default) uses the compatibility cutoff calculation: mean of chromosome fitted-value means plus \code{peakCutoffSd} times the square root of the summed per-chromosome \code{var(fitted)/n}. \code{"global_sd"} uses the genome-wide fitted-value median plus \code{peakCutoffSd} times the ordinary fitted-value SD.}
\item{peakIntervalMethod}{\code{"hpd_span"} (default) conservatively spans all KDE grid positions selected by the high-density mass rule, including separated modes. \code{"shortest_contiguous"} chooses the shortest single contiguous interval containing the requested mass.}
\item{peakResampleIterations}{Number of half-marker LOESS resamples for peak refinement. Default 1000.}
\item{randomSeed}{Base random seed used to make peak refinement reproducible. Default 1.}
\item{pairedEnd}{Whether expression-summary counting should treat the RNA-seq data as paired-end fragments.}
\item{ignoreStrand}{Whether strand is ignored in the descriptive expression summary. Default FALSE uses strand-aware counting.}
\item{expressionPseudocount}{Non-negative pseudocount added to WT and mutant mean counts before descriptive log2 fold change. Default 0.01 avoids Inf/NaN values at zero counts; set to 0 for the raw ratio.}
\item{exportAiccPlots}{Whether \code{outputMmapprData()} should write \code{aicc_plots.pdf} from the AICc span-search evaluations retained by \code{loessFit()}. Default FALSE.}
\item{overwrite}{If TRUE, explicitly permit clearing an existing non-empty output directory. Protected directories are still refused. Default FALSE.}
}
\value{A validated \code{MmapprParam} object.}
\description{Constructs the parameter object used by the MMAPPR2 workflow. Scalar settings are validated before filesystem side effects. BAM and FASTA indexes are created when missing and rebuilt when empty or older than the indexed input; the annotation input is preserved rather than destructively rewritten. BAM sequence names must match the reference FASTA exactly for regional queries.}
\examples{
\dontrun{
param <- mmapprParam(
  wtFiles = "wt.bam",
  mutFiles = "mut.bam",
  refFasta = "reference.fa",
  gtf = "genes.gtf"
)
}
}
