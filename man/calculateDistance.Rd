% Manually synchronized with R/distance.R.
\name{calculateDistance}
\alias{calculateDistance}
\title{Read BAM files and generate Euclidean-distance linkage data}
\usage{calculateDistance(mmapprData)}
\arguments{\item{mmapprData}{A \code{MmapprData} object.}}
\value{The object with its \code{snpDistance} slot populated.}
\description{Pileups A/C/G/T observations in annotated gene ranges, aggregates WT and mutant pools, filters low-depth/uninformative positions, and computes powered Euclidean distance between pool nucleotide-frequency vectors.}
