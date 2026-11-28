#' @title Generate plots and tables from MMAPPR2 data
#'
#' @name outputMmapprData
#'
#' @param mmapprData The \linkS4class{MmapprData} object to be output
#' @param palette Named figure palette: \code{"default"}, \code{"colorblind"},
#'   or \code{"grayscale"}. Styling has no effect on analytical results.
#' @param colors Optional named character vector of individual figure-color
#'   overrides (see \code{\link{plotMmapprData}} for supported roles).
#'
#' @return A \linkS4class{MmapprData} object after writing output files
#'   to the folder specified in the \code{outputFolder} slot of the
#'   \code{\link{MmapprParam}} used.
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
#' \dontrun{
#' md <- mmapprData(mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#'
#' outputMmapprData(postCandidatesMD)
#' }
#'
NULL


# Resolve all styling before making any output changes. Keep graphical settings
# separate from MmapprParam so existing RDS objects remain readable unchanged.
.plotPaletteColors <- function(palette = "default", colors = NULL) {
  available <- c("default", "colorblind", "grayscale")
  if (!is.character(palette) || length(palette) != 1L ||
      is.na(palette) || !(palette %in% available))
    stop("palette must be one of: ", paste(available, collapse = ", "))

  base <- switch(palette,
    default = c(fitted = "black", genomeRaw = "#999999AA", peakRaw = "black",
                peakFill = "#2ecc71", density = "#502ecc", cutoff = "grey",
                boundary = "grey", aiccCurve = "black", aiccSelected = "black"),
    colorblind = c(fitted = "#222222", genomeRaw = "#77777799", peakRaw = "#555555",
                   peakFill = "#E69F00AA", density = "#0072B2", cutoff = "#666666",
                   boundary = "#999999", aiccCurve = "#333333", aiccSelected = "#D55E00"),
    grayscale = c(fitted = "#151515", genomeRaw = "#99999999", peakRaw = "#555555",
                  peakFill = "#C5C5C5", density = "#454545", cutoff = "#777777",
                  boundary = "#B0B0B0", aiccCurve = "#333333", aiccSelected = "#111111")
  )

  if (is.null(colors)) return(base)
  if (!is.character(colors) || !length(colors) || is.null(names(colors)) ||
      anyNA(colors) || anyNA(names(colors)) || any(!nzchar(names(colors))) ||
      anyDuplicated(names(colors)))
    stop("colors must be a named character vector with unique, nonempty names and values")
  unknown <- setdiff(names(colors), names(base))
  if (length(unknown))
    stop("Unknown plot color role(s): ", paste(unknown, collapse = ", "),
         ". Supported roles: ", paste(names(base), collapse = ", "))
  for (name in names(colors)) {
    valid <- tryCatch({
      grDevices::col2rgb(colors[[name]], alpha = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (!valid) stop("Invalid R color for '", name, "': ", colors[[name]])
  }
  base[names(colors)] <- colors
  base
}

# Keep the capability checks shared between preflight and the PDF rendering
# functions. Partial pipeline results may contain error strings instead of
# chromosome lists; do not let a malformed entry invalidate good chromosomes.
.usableLoessSequences <- function(snpDistance) {
  names(snpDistance)[vapply(snpDistance, function(entry) {
    if (!is.list(entry) || !is.list(entry$loess)) return(FALSE)
    fit <- entry$loess
    if (!all(c("x", "y", "fitted") %in% names(fit)) ||
        !is.numeric(fit$x) || !is.numeric(fit$y) ||
        !is.numeric(fit$fitted)) return(FALSE)
    n <- length(fit$x)
    n > 0L && length(fit$y) == n && length(fit$fitted) == n &&
      all(is.finite(fit$x)) &&
      any(is.finite(fit$fitted)) && any(is.finite(fit$y))
  }, logical(1))]
}

.usablePeakSequences <- function(peaks, fitted) {
  names(peaks)[vapply(names(peaks), function(name) {
    if (!(name %in% fitted)) return(FALSE)
    entry <- peaks[[name]]
    if (!is.list(entry) || !is.numeric(entry$start) ||
        !is.numeric(entry$end) || length(entry$start) != 1L ||
        length(entry$end) != 1L || !is.finite(entry$start) ||
        !is.finite(entry$end) || entry$start > entry$end) return(FALSE)
    if (!is.null(entry$cutoff) && (!is.numeric(entry$cutoff) ||
                                   length(entry$cutoff) != 1L)) return(FALSE)
    density <- entry$densityData
    is.null(density) || (is.list(density) &&
      all(c("x", "y") %in% names(density)) &&
      is.numeric(density$x) && is.numeric(density$y) &&
      length(density$x) == length(density$y))
  }, logical(1))]
}

.usableAiccSequences <- function(snpDistance) {
  names(snpDistance)[vapply(snpDistance, function(entry) {
    if (!is.list(entry) || !is.data.frame(entry$aicc)) return(FALSE)
    tab <- entry$aicc
    all(c("spans", "aiccValues") %in% names(tab)) &&
      is.numeric(tab$spans) && is.numeric(tab$aiccValues) &&
      any(is.finite(tab$spans) & is.finite(tab$aiccValues))
  }, logical(1))]
}

# Only graphics are optional during normal analysis; candidate tables must still
# be written even when a plot-rendering error is encountered.
outputMmapprData <- function(mmapprData, palette = "default", colors = NULL) {
  stopifnot(is(mmapprData, "MmapprData"))
  plotColors <- .plotPaletteColors(palette, colors)

  if (!dir.exists(outputFolder(param(mmapprData)))) {
    dir.create(outputFolder(param(mmapprData)), recursive = TRUE, showWarnings = FALSE)
  }
  current_devs <- grDevices::dev.list()
  plotError <- NULL
  if (length(snpDistance(mmapprData)) > 0) {
    tryCatch({
      .plotGenomeDistance(mmapprData, plotColors = plotColors)
      .plotPeaks(mmapprData, plotColors = plotColors)
      if (isTRUE(exportAiccPlots(param(mmapprData))))
        .plotAicc(mmapprData, plotColors = plotColors)
    }, error = function(e) {
      plotError <<- e
    }, finally = {
      .closePlotDevices(current_devs)
    })
  }

  if (length(candidates(mmapprData)) > 0) {
    .writeCandidateTables(mmapprData@candidates,
                          outputFolder(param(mmapprData)))
  }

  if (!is.null(plotError)) stop(plotError)
  invisible(mmapprData)
}

.closePlotDevices <- function(existing) {
  now <- grDevices::dev.list()
  opened <- if (is.null(now)) integer() else now[!(now %in% existing)]
  for (d in rev(opened)) try(grDevices::dev.off(which = d), silent = TRUE)
  invisible(NULL)
}

#' @title Recreate MMAPPR2 figures from an RDS file or analysis object
#' @description
#' Writes only PDF figures from completed analysis data, without reading input
#' BAMs, recomputing results, writing candidate TSVs, or updating the RDS file.
#' This is independent of the saved \code{exportAiccPlots} setting. Unavailable
#' figures in partial analysis results are skipped with a warning.
#'
#' @param x A \code{MmapprData} object or path to a saved \code{mmappr_data.RDS}.
#' @param outputDir Destination for figures; by default a \code{replots}
#'   directory adjacent to the RDS, or inside the object's analysis directory.
#' @param plots Figure types to create: any nonempty selection from
#'   \code{c("genome", "peaks", "aicc")}. All three are requested by default.
#' @param palette Named figure palette: \code{default} (original colors),
#'   \code{colorblind}, or \code{grayscale}.
#' @param colors Optional named character vector overriding individual colors.
#'   Supported roles: \code{fitted}, \code{genomeRaw}, \code{peakRaw},
#'   \code{peakFill}, \code{density}, \code{cutoff}, \code{boundary},
#'   \code{aiccCurve}, and \code{aiccSelected}. Use any base R color name or
#'   hex color; alpha-channel hex colors (\code{#RRGGBBAA}) are accepted.
#' @param overwrite Whether existing selected figure PDFs may be replaced.
#'   Defaults to \code{FALSE}. No other files are touched.
#' @details
#' Figures are rendered to temporary PDFs first, so a plotting failure does
#' not replace existing figures. Chromosomes with unusable saved data are
#' omitted from the corresponding figures. File installation preserves an
#' individual previous PDF when its replacement fails; replotting does not
#' modify the saved data or export candidate tables.
#'
#' @return Invisibly, a named character vector of generated PDF paths. If no
#'   requested figure has its required results, the function raises an error.
#' @seealso \code{\link{outputMmapprData}}, \code{\link{mmappr}}
#' @export
#' @examples
#' \dontrun{
#' plotMmapprData("analysis/mmappr_data.RDS", outputDir = "replots")
#' plotMmapprData("analysis/mmappr_data.RDS", plots = "aicc",
#'                palette = "colorblind", colors = c(aiccSelected = "#D55E00"))
#' }
NULL

plotMmapprData <- function(x, outputDir = NULL,
                           plots = c("genome", "peaks", "aicc"),
                           palette = "default", colors = NULL,
                           overwrite = FALSE) {
  plotColors <- .plotPaletteColors(palette, colors)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
    stop("overwrite must be TRUE or FALSE")
  if (!is.character(plots) || !length(plots) || anyNA(plots) ||
      any(!(plots %in% c("genome", "peaks", "aicc"))))
    stop("plots must contain one or more of: genome, peaks, aicc")
  plots <- unique(plots)

  if (is.character(x)) {
    if (length(x) != 1L || is.na(x) || !nzchar(x) || !file.exists(x) || dir.exists(x))
      stop("x must name an existing MMAPPR2 RDS file")
    inputPath <- normalizePath(x, mustWork = TRUE)
    md <- readRDS(inputPath)
    if (is.null(outputDir)) outputDir <- file.path(dirname(inputPath), "replots")
  } else {
    md <- x
    if (!is(md, "MmapprData")) stop("x must be a MmapprData object or an RDS path")
    if (is.null(outputDir)) outputDir <- file.path(outputFolder(param(md)), "replots")
  }
  if (!is(md, "MmapprData")) stop("RDS file does not contain a MmapprData object")
  if (!is.character(outputDir) || length(outputDir) != 1L ||
      is.na(outputDir) || !nzchar(outputDir))
    stop("outputDir must be one nonempty directory path")
  if (file.exists(outputDir) && !dir.exists(outputDir))
    stop("outputDir exists but is not a directory: ", outputDir)

  # Select valid chromosomes without modifying the original analysis object.
  # This allows a partial run to produce every figure it has enough data for.
  fitted <- .usableLoessSequences(md@snpDistance)
  refined <- .usablePeakSequences(md@peaks, fitted)
  validAicc <- length(.usableAiccSequences(md@snpDistance)) > 0L
  available <- c(genome = length(fitted) > 0L, peaks = length(refined) > 0L,
                 aicc = validAicc)
  unavailable <- plots[!available[plots]]
  if (length(unavailable))
    warning("Skipping unavailable figure(s): ", paste(unavailable, collapse = ", "),
            ". Required saved analysis data are missing.")
  selected <- plots[available[plots]]
  if (!length(selected)) stop("No requested figures can be generated from this analysis")

  if (!dir.exists(outputDir) && !dir.create(outputDir, recursive = TRUE,
                                            showWarnings = FALSE) && !dir.exists(outputDir))
    stop("Could not create outputDir: ", outputDir)
  outputDir <- normalizePath(outputDir, mustWork = TRUE)
  filenames <- c(genome = "genome_plots.pdf", peaks = "peak_plots.pdf",
                 aicc = "aicc_plots.pdf")
  destinations <- file.path(outputDir, filenames[selected])
  names(destinations) <- selected
  existing <- destinations[file.exists(destinations)]
  if (!overwrite && length(existing))
    stop("Refusing to overwrite existing figure(s): ",
         paste(basename(existing), collapse = ", "), ". Set overwrite=TRUE to replace them.")
  # Check all destination paths, including dangling symlinks. An overwrite
  # must never follow a link into an unrelated file or directory.
  links <- Sys.readlink(destinations)
  if (any(dir.exists(destinations) | (!is.na(links) & nzchar(links))))
    stop("Refusing to overwrite a directory or symlink at a figure path")

  # Render into sibling temporary files first. On a graphics failure, old PDFs
  # remain intact and we remove all staged files without touching TSVs or RDS.
  staged <- vapply(selected, function(type) {
    tempfile(pattern = paste0(".mmappr2_", type, "_"),
             tmpdir = outputDir, fileext = ".pdf")
  }, character(1))
  on.exit(unlink(staged[file.exists(staged)], force = TRUE), add = TRUE)
  current_devs <- grDevices::dev.list()
  on.exit(.closePlotDevices(current_devs), add = TRUE)

  for (type in selected) {
    if (type == "genome") {
      plotMd <- md
      plotMd@snpDistance <- md@snpDistance[fitted]
      .plotGenomeDistance(plotMd, file = staged[[type]], plotColors = plotColors)
    } else if (type == "peaks") {
      plotMd <- md
      plotMd@peaks <- md@peaks[refined]
      .plotPeaks(plotMd, file = staged[[type]], plotColors = plotColors)
    } else {
      .plotAicc(md, file = staged[[type]], plotColors = plotColors)
    }
  }
  for (type in selected) .installPlotPdf(staged[[type]], destinations[[type]])
  invisible(destinations)
}

.installPlotPdf <- function(staged, destination, .renameFile = file.rename) {
  if (!file.exists(staged) || is.na(file.info(staged)$size) ||
      file.info(staged)$size <= 0)
    stop("Figure rendering did not produce a nonempty PDF: ", destination)
  if (!file.exists(destination)) {
    if (!.renameFile(staged, destination)) stop("Could not install PDF: ", destination)
    return(invisible(destination))
  }
  # Preserve the previous PDF in case installing the new one fails (including
  # Windows, which does not allow rename-over-existing semantics).
  previous <- tempfile(pattern = paste0(".", basename(destination), ".backup."),
                       tmpdir = dirname(destination))
  if (!.renameFile(destination, previous))
    stop("Could not move existing PDF aside: ", destination)
  if (!.renameFile(staged, destination)) {
    if (!.renameFile(previous, destination))
      stop("Could not install new PDF; previous PDF remains at: ", previous)
    stop("Could not install new PDF; previous file was restored: ", destination)
  }
  if (file.exists(previous) && unlink(previous, force = TRUE) != 0L)
    warning("Previous figure backup could not be removed: ", previous)
  invisible(destination)
}


# Export the AICc evaluations retained by loessFit(). One page is written per
# successfully fitted chromosome, with all finite evaluated spans shown and the
# selected span highlighted. This is output-only: it never reruns a LOESS fit or
# changes the optimized result stored in the MmapprData object.
.plotAicc <- function(mmapprData, file = NULL,
                      plotColors = .plotPaletteColors()) {
  usable <- .usableAiccSequences(mmapprData@snpDistance)
  if (!length(usable)) {
    warning("AICc plot export was requested, but no finite AICc search results are available")
    return(invisible(NULL))
  }

  if (is.null(file)) file <- file.path(outputFolder(param(mmapprData)), "aicc_plots.pdf")
  grDevices::pdf(file, width = 8.5, height = 6.5)
  deviceOpen <- TRUE
  on.exit(if (deviceOpen) try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  for (seqname in usable) {
    result <- mmapprData@snpDistance[[seqname]]
    tab <- result$aicc
    ok <- is.finite(tab$spans) & is.finite(tab$aiccValues)
    tab <- tab[ok, , drop = FALSE]
    tab <- tab[order(tab$spans), , drop = FALSE]
    if (!nrow(tab)) next

    graphics::plot(tab$spans, tab$aiccValues, type = "b", pch = 16,
                   xlab = "LOESS span", ylab = "AICc",
                   main = paste(seqname, "LOESS span optimization"),
                   xlim = .safeXLim(tab$spans), ylim = .safeYLim(tab$aiccValues),
                   col = plotColors[["aiccCurve"]])

    best <- result$bestSpan
    if (is.numeric(best) && length(best) == 1L && is.finite(best)) {
      nearest <- which.min(abs(tab$spans - best))
      graphics::abline(v = best, lty = 2, col = plotColors[["aiccSelected"]])
      graphics::points(tab$spans[nearest], tab$aiccValues[nearest],
                       pch = 19, cex = 1.4, col = plotColors[["aiccSelected"]])
      graphics::legend("topright",
                       legend = sprintf("selected span = %s", format(best)),
                       lty = 2, pch = 19, bty = "n", col = plotColors[["aiccSelected"]])
    }
  }

  grDevices::dev.off()
  deviceOpen <- FALSE
  invisible(NULL)
}


.atomicSaveRDS <- function(object, file, .renameFile = file.rename) {
  dir <- dirname(file)
  if (!dir.exists(dir)) stop("RDS destination directory does not exist: ", dir)
  tmp <- tempfile(pattern = paste0(".", basename(file), "."), tmpdir = dir)
  on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)

  saveRDS(object, tmp)
  info <- file.info(tmp)
  if (!file.exists(tmp) || is.na(info$size) || info$size <= 0)
    stop("Temporary RDS serialization failed: ", tmp)

  if (!.renameFile(tmp, file)) {
    if (!file.exists(file))
      stop("Could not move completed RDS into place: ", file)

    backup <- tempfile(pattern = paste0(".", basename(file), ".previous."), tmpdir = dir)
    if (file.exists(backup)) unlink(backup, force = TRUE)
    if (!.renameFile(file, backup))
      stop("Could not move existing RDS aside for replacement: ", file)

    installed <- .renameFile(tmp, file)
    if (!installed) {
      restored <- .renameFile(backup, file)
      if (!restored) {
        stop("Could not install completed RDS and could not restore the previous file. ",
             "Previous state remains at: ", backup)
      }
      stop("Could not install completed RDS; previous file was restored: ", file)
    }

    if (file.exists(backup)) {
      rc <- unlink(backup, force = TRUE)
      if (!identical(as.integer(rc), 0L) || file.exists(backup))
        warning("Replacement succeeded but previous RDS backup could not be removed: ", backup)
    }
  }
  invisible(file)
}


.defaultOutputFolder <- function() {
  stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # tempfile() adds a collision-resistant suffix without consuming R's
  # random-number stream, so two analyses started in the same second do not clash.
  basename(tempfile(pattern = paste0("mmappr2_", stamp, "_"), tmpdir = getwd()))
}


#' @title Generate temporary output folder
#' @name tempOutputFolder
#' @description
#' Returns a unique temporary output-directory path. The directory is created
#' by \code{mmapprParam()} when the parameter object is constructed.
#'
#' @return The path to the temporary directory
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(refFasta = MMAPPR2data::goldenFasta(),
#'                                 wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
NULL


tempOutputFolder <- function() {
  # Return a unique path; mmapprParam() creates it. This avoids same-second name
  # collisions in tests and parallel workflows.
  tempfile(pattern = paste0("mmappr2_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_"),
           tmpdir = tempdir())
}


.prepareOutputFolder <- function(outputFolder, overwrite = FALSE) {
  if (!is.character(outputFolder) || length(outputFolder) != 1L || !nzchar(outputFolder))
    stop("outputFolder must be one non-empty path")

  if (dir.exists(outputFolder)) {
    existingPath <- normalizePath(outputFolder, mustWork = TRUE)

    if (isTRUE(overwrite)) {
      protected <- unique(normalizePath(c(path.expand("~"), getwd(), tempdir()),
                                        mustWork = TRUE))
      isRoot <- identical(dirname(existingPath), existingPath)
      parent <- dirname(existingPath)
      isTopLevel <- !isRoot && identical(dirname(parent), parent)
      if (isRoot || isTopLevel || existingPath %in% protected)
        stop("Refusing to overwrite protected directory: ", existingPath)
    }

    contents <- list.files(existingPath, all.files = TRUE, no.. = TRUE)
    if (length(contents) > 0L) {
      if (!isTRUE(overwrite)) {
        stop("Output folder already exists and is not empty: ", outputFolder,
             ". Choose a new folder or set overwrite=TRUE explicitly.")
      }
      unlink(file.path(existingPath, contents), recursive = TRUE, force = TRUE)
      stillThere <- list.files(existingPath, all.files = TRUE, no.. = TRUE)
      if (length(stillThere) > 0L)
        stop("Could not fully clear output folder: ", existingPath)
    }
  } else {
    ok <- dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(outputFolder)) stop("Could not create output folder: ", outputFolder)
  }

  outputFolder <- normalizePath(outputFolder, mustWork = TRUE)
  logPath <- file.path(outputFolder, "mmappr2.log")
  if (!file.exists(logPath) && !file.create(logPath))
    stop("Could not create log file in output folder: ", outputFolder)
  outputFolder
}



.safeYLim <- function(x, upperPad = 0.10) {
  # Base plot() rejects zero-width or non-finite limits. Constant/degenerate
  # fits can occur on small datasets, so construct a finite visible range.
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(0, 1))
  lo <- min(x); hi <- max(x)
  if (lo == hi) {
    pad <- if (lo == 0) 1 else max(abs(lo) * upperPad, .Machine$double.eps^0.5)
    return(c(lo - pad, hi + pad))
  }
  c(lo, hi + (hi - lo) * upperPad)
}

.safeXLim <- function(x, padFraction = 0.01) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(0, 1))
  lo <- min(x); hi <- max(x)
  if (lo == hi) {
    pad <- if (lo == 0) 1 else max(abs(lo) * padFraction, 1)
    return(c(lo - pad, hi + pad))
  }
  c(lo, hi)
}

# Genome-level plotting uses safe limits, normalized numeric LOESS coordinates, and
# namespaced sequence ordering. Peak-cutoff/refined-region diagnostics are drawn by
# .plotPeaks(); device handling here avoids closing unrelated user devices.
.plotGenomeDistance <- function(mmapprData, savePdf = TRUE, file = NULL,
                                plotColors = .plotPaletteColors()) {
  #generate one big dataframe for plots, along with break and label points
  tailPos <- 0
  plotDf <- NULL
  breaks <- tailPos
  labelpos <- NULL
  for (i in GenomeInfoDb::orderSeqlevels(names(mmapprData@snpDistance))) {
    if (!is(mmapprData@snpDistance[[i]], 'list'))
      next
    else if (!("loess" %in% names(mmapprData@snpDistance[[i]])))
      stop(sprintf("Distance list for sequence %s missing loess fit data",
                   names(mmapprData@snpDistance)[i]))

    chrLoess <- mmapprData@snpDistance[[i]]$loess
    chrX <- as.numeric(chrLoess$x)
    chrDf <- data.frame(pos = chrX + tailPos,
                        seqname = names(mmapprData@snpDistance)[i],
                        fitted = as.numeric(chrLoess$fitted),
                        unfitted = as.numeric(chrLoess$y))
    plotDf <- rbind(plotDf, chrDf)
    tailPos <- chrDf$pos[nrow(chrDf)]
    breaks <- c(breaks, tailPos)
    labelpos <-
      c(labelpos, (breaks[length(breaks)]+breaks[length(breaks)-1])/2)
  }

  if (is.null(plotDf) || nrow(plotDf) == 0L)
    stop("No successfully fitted chromosome data are available for genome plotting")

  openedPdf <- FALSE
  if (savePdf) {
    if (is.null(file)) file <- file.path(outputFolder(param(mmapprData)), "genome_plots.pdf")
    grDevices::pdf(file, width = 11, height = 8.5)
    openedPdf <- TRUE
    on.exit(if (openedPdf) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  }
  par(mfrow = c(2,1))

  plot(x = plotDf$pos, y = plotDf$fitted, type = 'l',
       ylim = .safeYLim(plotDf$fitted),
       ylab = NA, col = plotColors[["fitted"]],
       xaxt = 'n', xaxs = 'i', xlab = "Chromosome")
  abline(v = breaks[seq_len(length(breaks) - 1L)], col = plotColors[["boundary"]])
  mtext(unique(sub("^chr", "", plotDf$seqname[!is.na(plotDf$seqname)])),
        at = labelpos, side = 1, cex = .6)
  mtext(substitute("ED"^p~ ~"(Loess fit)",
                   list(p = mmapprData@param@distancePower)),
        side = 2, line = 2)

  plot(x = plotDf$pos, y = plotDf$unfitted,
       pch = 16, cex = .8, col = plotColors[["genomeRaw"]],
       ylim = .safeYLim(plotDf$unfitted),
       ylab = NA,
       xaxt = 'n', xaxs = 'i', xlab = "Chromosome" )
  abline(v = (breaks), col = plotColors[["boundary"]])
  mtext(unique(sub("^chr", "", plotDf$seqname[!is.na(plotDf$seqname)])),
        at = labelpos, side = 1, cex = .6)
  mtext(substitute("ED"^p,
                   list(p = mmapprData@param@distancePower)),
        side = 2, line = 2)

  if (openedPdf) {
    grDevices::dev.off()
    openedPdf <- FALSE
  }
}


.plotPeaks <- function(mmapprData, file = NULL,
                       plotColors = .plotPaletteColors()) {
  if (length(mmapprData@peaks) == 0L) return(invisible(NULL))
  if (is.null(file)) file <- file.path(outputFolder(param(mmapprData)), "peak_plots.pdf")
  grDevices::pdf(file, width = 11, height = 8.5)
  deviceOpen <- TRUE
  on.exit(if (deviceOpen) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  par(mfrow = c(2,1), mar = c(5, 4, 4, 4))

  for (seqname in names(mmapprData@peaks)) {
    chrLoess <- mmapprData@snpDistance[[seqname]]$loess
    chrX <- as.numeric(chrLoess$x)
    fitted <- as.numeric(chrLoess$fitted)
    rawY <- as.numeric(chrLoess$y)
    start <- mmapprData@peaks[[seqname]]$start
    end <- mmapprData@peaks[[seqname]]$end
    densityData <- mmapprData@peaks[[seqname]]$densityData
    yLim <- .safeYLim(fitted)
    xLim <- .safeXLim(chrX) * 1E-6

    plot(chrX/1000000, fitted, type = 'l',
         ylim = yLim, xlim = xLim,
         xlab = paste(seqname, "Base Position (MB)"),
         ylab = NA, xaxs = 'i', col = plotColors[["fitted"]])
    mtext(substitute("ED"^p~ ~"(Loess fit)",
                     list(p = mmapprData@param@distancePower)),
          side = 2, line = 2)

    inside <- is.finite(chrX) & is.finite(fitted) &
              chrX >= start & chrX <= end
    if (any(inside)) {
      # Shade to the actual plot baseline rather than an arbitrary -5,
      # which could distort clipping when the plotted range is small/positive.
      shadeX <- c(start, chrX[inside], end)
      shadeY <- c(yLim[1], fitted[inside], yLim[1])
      polygon(shadeX/1000000, shadeY, col = plotColors[["peakFill"]], border = NA)
    }
    cutoff <- mmapprData@peaks[[seqname]]$cutoff
    if (length(cutoff) == 1L && is.finite(cutoff))
      abline(h = cutoff, lty = 2, col = plotColors[["cutoff"]])

    # Overlay the resampling density only when finite support remains inside the
    # plotted chromosome range. Degenerate one-point densities use a point rather
    # than an invisible zero-length line.
    if (length(densityData) > 0L && length(densityData$x) > 0L) {
      dx <- as.numeric(densityData$x)
      dy <- as.numeric(densityData$y)
      posMatch <- is.finite(dx) & is.finite(dy) & dx >= min(chrX, na.rm = TRUE) &
                  dx <= max(chrX, na.rm = TRUE)
      if (any(posMatch)) {
        dx <- dx[posMatch]; dy <- dy[posMatch]
        par(new = TRUE)
        plot(dx, dy, type = if (length(dx) > 1L) 'l' else 'p',
             ylim = .safeYLim(dy), xlim = .safeXLim(chrX),
             ann = FALSE, xaxs = 'i', xaxt = 'n', yaxt = 'n', col = plotColors[["density"]])
        axis(side = 4, col = plotColors[["density"]])
        mtext(side = 4, line = 2, 'Density', col = plotColors[["density"]])
        legend('topright',
               legend = c("Fitted Distance Curve", "Peak Resampling Distribution"),
               col = c(plotColors[["fitted"]], plotColors[["density"]]),
               lty = c(1, 1), cex = 0.8, lwd = 3)
      }
    }

    if (length(chrX) > 0L) {
      plot(chrX/1000000, rawY, pch = 16, cex = .6,
           col = plotColors[["peakRaw"]],
           ylim = .safeYLim(rawY), xlim = xLim,
           ylab = NA, xlab = paste(seqname, "Base Position (MB)"), xaxs = 'i')
      mtext(substitute("ED"^p, list(p = mmapprData@param@distancePower)),
            side = 2, line = 2)
    }
  }
  grDevices::dev.off()
  deviceOpen <- FALSE
  invisible(NULL)
}



.tsvSafeDataFrame <- function(x) {
  df <- as.data.frame(x)
  # Even a zero-row table can contain a list/List column whose class makes
  # write.table() fail before it writes the header. Only the truly zero-column
  # case can bypass column sanitization safely.
  if (ncol(df) == 0L) return(df)

  sanitizeText <- function(x) {
    # A literal tab/newline inside an annotation value would
    # corrupt an unquoted TSV row. GTF-derived labels should not contain these,
    # but sanitize defensively so output remains one record per line.
    out <- x
    ok <- !is.na(out)
    out[ok] <- gsub("[\t\r\n]+", " ", out[ok])
    out
  }

  collapseCell <- function(x) {
    if (length(x) == 0L) return("")
    vals <- tryCatch(as.character(x), error = function(e) character())
    if (length(vals) == 0L) return("")
    if (all(is.na(vals))) return(NA_character_)
    paste(sanitizeText(vals[!is.na(vals)]), collapse = ",")
  }

  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.character(col)) {
      df[[nm]] <- sanitizeText(col)
      next
    }
    # Leave ordinary scalar vectors alone; only flatten columns that base
    # write.table() cannot safely serialize. The [[i]] path is important for
    # CharacterList/IntegerList/list columns because one table row can contain
    # multiple annotation values.
    listLike <- is.list(col) || methods::is(col, "List") || !is.atomic(col)
    if (!listLike) next
    df[[nm]] <- vapply(seq_len(nrow(df)), function(i) {
      cell <- tryCatch(col[[i]], error = function(e) col[i])
      collapseCell(cell)
    }, character(1))
  }
  df
}

.writeTsv <- function(x, file) {
  utils::write.table(.tsvSafeDataFrame(x), file = file, sep = "\t",
                     quote = FALSE, row.names = FALSE, na = "NA")
  invisible(file)
}


.originalDetectedMutationTable <- function(x) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  originalCols <- c("seqnames", "start", "end", "width", "strand",
                    "ref", "alt", "refDepth", "altDepth", "peakDensity",
                    "GENEID", "TXID", "PROTEINLOC", "CONSEQUENCE",
                    "REFAA", "VARAA")

  if (!"refDepth" %in% names(df) && "mutRefDepth" %in% names(df))
    df$refDepth <- df$mutRefDepth
  if (!"altDepth" %in% names(df) && "mutAltDepth" %in% names(df))
    df$altDepth <- df$mutAltDepth

  if ("PROTEINLOC" %in% names(df)) {
    proteinLoc <- df$PROTEINLOC
    listLike <- is.list(proteinLoc) || methods::is(proteinLoc, "List") || !is.atomic(proteinLoc)
    if (listLike) {
      df$PROTEINLOC <- vapply(seq_len(nrow(df)), function(i) {
        value <- tryCatch(proteinLoc[[i]], error = function(e) proteinLoc[i])
        value <- tryCatch(as.character(value), error = function(e) character())
        if (length(value)) value[[1L]] else NA_character_
      }, character(1))
    }
  }

  # Assign a zero-length vector when x has zero rows. Assigning scalar NA to a
  # zero-row data.frame can otherwise create a replacement-length error instead
  # of producing the desired header-only old-format table.
  for (nm in setdiff(originalCols, names(df))) df[[nm]] <- rep(NA, nrow(df))
  df[, originalCols, drop = FALSE]
}


.writeCandidateTables <- function(candList, outputFolder){
  seqnames <- unique(c(names(candList$snps), names(candList$effects), names(candList$diff)))
  for (seqname in seqnames) {
    # Always write the full candidate SNV set, including noncoding
    # variants and the new WT-vs-mutant allele-frequency evidence.
    if (!is.null(candList$snps[[seqname]])) {
      .writeTsv(candList$snps[[seqname]],
                file.path(outputFolder, paste0("AllCandidateVariantsFor", seqname, ".tsv")))
    }

    if (!is.null(candList$effects[[seqname]])) {
      .writeTsv(.originalDetectedMutationTable(candList$effects[[seqname]]),
                file.path(outputFolder, paste0("DetectedMutationsFor", seqname, ".tsv")))
    }

    if (!is.null(candList$diff[[seqname]])) {
      .writeTsv(candList$diff[[seqname]],
                file.path(outputFolder, paste0("DifferentiallyExpressedGenesFor", seqname, ".tsv")))
    }
  }
}

