% Manually synchronized with the roxygen documentation in R/output.R.
\name{plotMmapprData}
\alias{plotMmapprData}
\title{Regenerate MMAPPR2 figures from saved analysis results}
\description{
Write figure PDFs from a saved \code{MmapprData} RDS file or an in-memory
\code{MmapprData} object. Does not load original BAM/FASTA files, recompute
analysis, modify the RDS, or export candidate TSVs.
}
\usage{
plotMmapprData(x, outputDir = NULL,
               plots = c("genome", "peaks", "aicc"),
               palette = "default", colors = NULL,
               overwrite = FALSE)
}
\arguments{
  \item{x}{A \code{MmapprData} object or path to a saved \code{mmappr_data.RDS}.}
  \item{outputDir}{Figure directory. Defaults to a \code{replots} directory beside
    the RDS, or within the analysis directory when passing an object.}
  \item{plots}{One or more of \code{"genome"}, \code{"peaks"}, and
    \code{"aicc"}. Requests all three by default. Unavailable types are skipped.}
  \item{palette}{One of \code{"default"}, \code{"colorblind"}, or
    \code{"grayscale"}.}
  \item{colors}{Optional named character vector of overrides. The roles are
    \code{fitted}, \code{genomeRaw}, \code{peakRaw}, \code{peakFill},
    \code{density}, \code{cutoff}, \code{boundary}, \code{aiccCurve}, and
    \code{aiccSelected}. Accepts R color names or hexadecimal colors, including
    \code{#RRGGBBAA} alpha.}
  \item{overwrite}{Whether the selected figure PDFs may replace existing files;
    defaults to \code{FALSE}. Never clears the directory or replaces unselected files.}
}
\value{
Invisibly returns a named character vector of written figure paths.
Errors if none of the requested figure types can be generated.
}
\details{
Replotting uses only stored LOESS fits, peak intervals and density results,
and AICc span evaluations. It can create AICc plots even if
\code{exportAiccPlots} was \code{FALSE} in the original run. Chromosomes with
unusable stored data are omitted from the corresponding plots. An unavailable
figure type is skipped with a warning.

Figures are first rendered to temporary PDFs, so a rendering failure leaves
existing PDFs intact. Any individually failed overwrite restores its previous
file when possible. This function does not create or update candidate tables
or the saved RDS file.
}
\examples{
\dontrun{
plotMmapprData("results/mmappr_data.RDS")
plotMmapprData("results/mmappr_data.RDS", outputDir = "paper_figures",
               plots = c("genome", "peaks"), palette = "colorblind",
               colors = c(peakFill = "#E69F00AA", density = "#0072B2"))
}
}
\seealso{\code{\link{mmappr}}, \code{\link{outputMmapprData}}}
