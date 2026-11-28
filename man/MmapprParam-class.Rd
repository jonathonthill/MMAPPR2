% Manually synchronized with the S4 class definition in R/param.R.
\name{MmapprParam-class}
\alias{MmapprParam-class}
\alias{MmapprParam}
\title{MmapprParam class}
\description{Stores validated inputs and analysis settings for an MMAPPR2 run. Prefer \code{\link{mmapprParam}} over calling \code{new("MmapprParam")} directly because the constructor performs indexing and cross-parameter validation.}
\section{Slots}{
\describe{
\item{\code{wtFiles}}{Wild-type-pool \code{BamFileList}.}
\item{\code{mutFiles}}{Mutant-pool \code{BamFileList}.}
\item{\code{refFasta}}{Reference FASTA path.}
\item{\code{refGenome}}{Indexed \code{FaFile} corresponding to \code{refFasta}.}
\item{\code{gtf}}{GTF/GFF annotation path.}
\item{\code{outputFolder}}{Result directory.}
\item{\code{includeScaffolds}}{Whether non-standard contigs are included.}
\item{\code{minDepth}}{Minimum mapping-marker depth.}
\item{\code{homozygoteCutoff}}{WT major-allele-frequency threshold for excluding uninformative mapping positions.}
\item{\code{minBaseQuality}}{Minimum base quality.}
\item{\code{minMapQuality}}{Minimum mapping quality.}
\item{\code{fileAggregation}}{Replicate aggregation strategy.}
\item{\code{distancePower}}{Power applied to Euclidean distance.}
\item{\code{peakIntervalWidth}}{Target resampling-density mass for the refined interval.}
\item{\code{loessOptResolution}}{LOESS span-search resolution.}
\item{\code{loessOptCutFactor}}{LOESS coarse-to-fine refinement factor.}
\item{\code{maxPileupDepth}}{Per-BAM linkage-stage pileup cap; candidate calling uses a separate 250-read-per-BAM cap.}
\item{\code{candidateMinDepth}}{Optional mutant candidate depth floor; default 1 is inert given the two-ALT-read rule.}
\item{\code{candidateMinAltDepth}}{Minimum mutant ALT-supporting reads.}
\item{\code{candidateMinAltFreq}}{Minimum mutant ALT frequency.}
\item{\code{candidateMaxWtAltFreq}}{Optional maximum WT ALT frequency.}
\item{\code{candidateMinDeltaAF}}{Optional minimum mutant-minus-WT ALT-frequency difference.}
\item{\code{candidatePoolMode}}{Candidate-pileup memory strategy: adaptive \code{auto}, one-shot-preferred \code{memory}, or always-bounded \code{chunked}.}
\item{\code{candidateChunkSize}}{Maximum candidate-pileup chunk width in bases; 0 selects an adaptive width.}
\item{\code{peakCutoffSd}}{Peak-cutoff variability multiplier.}
\item{\code{peakCutoffMethod}}{Initial peak-cutoff method.}
\item{\code{peakIntervalMethod}}{Refined-interval method.}
\item{\code{peakResampleIterations}}{Number of peak-resampling iterations.}
\item{\code{randomSeed}}{Base peak-resampling seed.}
\item{\code{pairedEnd}}{Paired-end expression-counting flag.}
\item{\code{ignoreStrand}}{Expression-counting strand flag; default FALSE uses strand-aware counting.}
\item{\code{expressionPseudocount}}{Non-negative pseudocount added to WT and mutant mean counts before descriptive log2 fold change. Default 0.01 avoids infinite values at zero counts; set to 0 for the raw ratio.}
\item{\code{exportAiccPlots}}{Whether AICc span-search diagnostics should be written to \code{aicc_plots.pdf}.}
}}
\seealso{\code{\link{mmapprParam}}, \code{\link{MmapprParam-functions}}}
