% Manually synchronized with R/mmapprData.R.
\name{MmapprData-Class}
\alias{MmapprData-Class}
\alias{MmapprData}
\title{MmapprData Class}
\description{Stores inputs and intermediate/results data from each stage of the MMAPPR2 pipeline.}
\section{Slots}{
\describe{
\item{\code{param}}{A \code{MmapprParam} object.}
\item{\code{snpDistance}}{Per-chromosome linkage data. Before LOESS this can contain WT/mutant counts and a distance table; after LOESS it contains the fitted model and optimization information.}
\item{\code{peaks}}{Per-chromosome refined linkage-peak information, including interval coordinates and resampling density.}
\item{\code{candidates}}{Nested lists containing all candidate SNVs, coding-effect annotations, and descriptive expression summaries for peak regions.}
}}
\seealso{\code{\link{mmappr}}, \code{\link{mmapprData}}}
