% Manually synchronized with R/peaks.R.
\name{prePeak}
\alias{prePeak}
\title{Identify chromosomes containing linkage peaks}
\usage{prePeak(mmapprData)}
\arguments{\item{mmapprData}{A \code{MmapprData} object after \code{loessFit}.}}
\value{The object with its \code{peaks} list initialized for chromosomes whose fitted linkage signal crosses the selected genome-wide threshold.}
\description{Calculates an initial genome-wide linkage cutoff using \code{peakCutoffMethod}. The default \code{"legacy_current"} preserves the old cutoff calculation; optional \code{"global_sd"} uses the genome-wide fitted-value median plus \code{peakCutoffSd} times the ordinary SD. Chromosomes whose fitted LOESS signal crosses the resulting cutoff are retained for refinement.}
\seealso{\code{\link{peakRefinement}}, \code{\link{mmapprParam}}}
