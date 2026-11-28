% Manually synchronized with R/peaks.R.
\name{peakRefinement}
\alias{peakRefinement}
\title{Refine linkage peaks using marker resampling}
\usage{peakRefinement(mmapprData)}
\arguments{\item{mmapprData}{A \code{MmapprData} object after \code{prePeak}.}}
\value{The object with per-chromosome refined intervals, resampling density, full-data and resampled peak positions, success rate, interval method, and deterministic seed stored in \code{peaks}.}
\description{Repeatedly refits LOESS to random half-samples of markers and estimates a density for the resulting peak locations. The default \code{peakIntervalMethod="hpd_span"} conservatively spans all grid locations selected by the high-density mass rule; \code{"shortest_contiguous"} is an opt-in single-interval alternative that can exclude secondary modes.}
\seealso{\code{\link{prePeak}}, \code{\link{generateCandidates}}, \code{\link{mmapprParam}}}
