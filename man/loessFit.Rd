% Manually synchronized with R/loess.R.
\name{loessFit}
\alias{loessFit}
\title{Perform optimized LOESS regression for each chromosome}
\usage{loessFit(mmapprData)}
\arguments{\item{mmapprData}{A \code{MmapprData} object after \code{calculateDistance}.}}
\value{The object with an optimized \code{loess} fit stored for each successful chromosome in \code{snpDistance}.}
\description{Uses a coarse-to-fine AICc search to select a chromosome-specific LOESS span and smooth the noisy marker-level Euclidean-distance signal.}
