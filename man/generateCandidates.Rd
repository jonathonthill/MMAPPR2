% Manually synchronized with R/candidates.R.
\name{generateCandidates}
\alias{generateCandidates}
\title{Generate and annotate candidate variants in refined peak regions}
\usage{generateCandidates(md)}
\arguments{\item{md}{A \code{MmapprData} object after \code{peakRefinement}.}}
\value{The object with \code{candidates$snps}, \code{candidates$effects}, and \code{candidates$diff} populated for peak chromosomes.}
\description{Builds an explicit SNV candidate set inside each refined interval using pooled mutant depth, mutant ALT-supporting read count, and mutant ALT frequency; attaches WT allele support and optional WT/delta-AF filters; predicts coding consequences locally with \code{predictCoding}; and calculates a descriptive RNA raw-count fold-change summary. Exact indel/haplotype calling and formal differential-expression statistics are intentionally not implemented.}
\seealso{\code{\link{peakRefinement}}, \code{\link{mmapprParam}}}
