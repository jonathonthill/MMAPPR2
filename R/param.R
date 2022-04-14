#' @title MmapprParam Class
#' 
#' @name MmapprParam
#' 
#' @usage \code{MmapprParam} stores parameters for running \code{\link{mmappr}}.
#'
#' @slot wtFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the wild-type pool to be analyzed.
#' @slot mutFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the mutant pool to be analyzed.
#' @slot refFasta The path to the fasta file genome.
#' @slot gtf The path to a gtf-formatted annotation file. 
#' @slot outputFolder Length-one character vector specifying where to save
#'   output, including a \code{\linkS4class{MmapprData}} stored as
#'   \code{mmappr_data.RDS}, \code{mmappr2.log}, a \code{.tsv} file
#'   for each peak chromosome containing candidate mutations, and PDF plots
#'   of both the entire genome and peak chromosomes. Defaults to an
#'   automatically generated \code{mmappr2_<timestamp>}.
#' @slot includeScaffolds Logical indicating whether non-standard chromosomes
#'   should be included. False (default) trims the seqnames of the reference to 
#'   include only standard chromosome names based on UCSC or Ensembl 
#'   conventions. In all cases, the mitochondrial chromosome is removed.
#' @slot minDepth Length-one integer vector determining minimum depth
#'   required for a position to
#'   be considered in the analysis. Defaults to 20.
#' @slot homozygoteCutoff Length-one numeric vector between \code{0} and
#'   \code{1} specifying threshold for throwing out base pairs on account
#'   of homozygosity. Positions with high major allele frequency in the
#'   wild-type pool are unlikely to exhibit polymorphism and are thus thrown
#'   out when they exceed this cutoff. Defaults to \code{0.95}.
#' @slot minBaseQuality Length-one numeric vector indicating minimum base
#'   call quality to consider in analysis. Read positions with qualities
#'   below this score will be thrown out. Defaults to 20.
#' @slot minMapQuality Length-one numeric vector indicating minimum read
#'   mapping quality to consider in analysis. Reads with qualities below
#'   this score will be thrown out. Defaults to 30.
#' @slot fileAggregation A length-one character vector determining strategy
#'   for aggregating base calls when multiple wild-type or multiple mutant
#'   files are provided.
#'   When 'sum', average base call proportions are computed using
#'   the read counts from each file, effectively weighting files
#'   with higher counts at a given position. When equal to 'mean', the
#'   base call proportions as well as read depths, rather than the absolute count,
#'   are averaged across files, which is useful when you want to weight each
#'   replicate evenly without regards to differing depth.
#' @slot distancePower Length-one numeric vector determing to what power
#'   Euclidean distance values are raised before fitting. Higher powers tend
#'   to increase high values and decrease low values, exaggerating the
#'   variation in the data. Default of 4.
#' @slot peakIntervalWidth Length-one numeric vector between \code{0} and
#'   \code{1} specifying desired width of linkage region(s). The default value
#'   of \code{0.95}, for example, yields peak regions defined as including the
#'   top 95\% of SNPs in the peak region, as determined by the peak
#'   resampling distribution.
#' @slot loessOptResolution Length-one numeric vector between \code{0} and
#'   \code{1} specifying
#'   desired resolution for Loess fit optimization. The default of \code{0.001},
#'   for example, indicates that the span ultimately chosen will perform better
#'   than its neighbor values at \code{+-0.001}.
#' @slot loessOptCutFactor Length-one numeric vector between \code{0} and
#'   \code{1} specifying how aggressively the Loess
#'   optimization algorithm proceeds. For example, with a default of \code{0.1}
#'   different spans at intervals of \code{0.001} would be evaluated after
#'   intervals of \code{0.01}.
#' @slot refGenome A \code{\link[GmapGenome]{GmapGenome}} object. This is 
#'   generated internally from the provided refFasta file to ensure 
#'   consistency.
#' @rdname MmapprParam
#' 
NULL

setClass("MmapprParam",
         representation(
           wtFiles = "BamFileList", 
           mutFiles = "BamFileList",
           refFasta = "character",
           refGenome = "GmapGenome",
           gtf = "character", 
           outputFolder = "character", 
           includeScaffolds = "logical",
           minDepth = "numeric",
           homozygoteCutoff = "numeric", 
           minBaseQuality = "numeric",
           minMapQuality = "numeric",
           fileAggregation = "character", 
           distancePower = "numeric",
           peakIntervalWidth = "numeric", 
           loessOptResolution = "numeric",
           loessOptCutFactor = "numeric"
         )
)


#' @title MmapprParam Constructor
#' 
#' @name mmapprParam
#'
#' @usage Creates a new instance of a \code{\linkS4class{MmapprParam}} class object.
#'
#' @param wtFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the wild-type pool to be analyzed.
#' @param mutFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the mutant pool to be analyzed.
#' @param refFasta The path to the fasta file genome.
#' @param gtf The path to a gtf-formatted annotation file. 
#' @param outputFolder Length-one character vector specifying where to save
#'   output, including a \code{\linkS4class{MmapprData}} stored as
#'   \code{mmappr_data.RDS}, \code{mmappr2.log}, a \code{.tsv} file
#'   for each peak chromosome containing candidate mutations, and PDF plots
#'   of both the entire genome and peak chromosomes. Defaults to an
#'   automatically generated \code{mmappr2_<timestamp>}.
#' @param includeScaffolds Logical indicating whether non-standard chromosomes
#'   should be included. False (default) trims the seqnames of the reference to 
#'   include only standard chromosome names based on UCSC or Ensembl 
#'   conventions. In all cases, the mitochondrial chromosome is removed.
#' @param minDepth Length-one integer vector determining minimum depth
#'   required for a position to
#'   be considered in the analysis. Defaults to 20.
#' @param homozygoteCutoff Length-one numeric vector between \code{0} and
#'   \code{1} specifying threshold for throwing out base pairs on account
#'   of homozygosity. Positions with high major allele frequency in the
#'   wild-type pool are unlikely to exhibit polymorphism and are thus thrown
#'   out when they exceed this cutoff. Defaults to \code{0.95}.
#' @param minBaseQuality Length-one numeric vector indicating minimum base
#'   call quality to consider in analysis. Read positions with qualities
#'   below this score will be thrown out. Defaults to 20.
#' @param minMapQuality Length-one numeric vector indicating minimum read
#'   mapping quality to consider in analysis. Reads with qualities below
#'   this score will be thrown out. Defaults to 30.
#' @param fileAggregation A length-one character vector determining strategy
#'   for aggregating base calls when multiple wild-type or multiple mutant
#'   files are provided.
#'   When 'sum', average base call proportions are computed using
#'   the read counts from each file, effectively weighting files
#'   with higher counts at a given position. When equal to 'mean', the
#'   base call proportions as well as read depths, rather than the absolute count,
#'   are averaged across files, which is useful when you want to weight each
#'   replicate evenly without regards to differing depth.
#' @param distancePower Length-one numeric vector determing to what power
#'   Euclidean distance values are raised before fitting. Higher powers tend
#'   to increase high values and decrease low values, exaggerating the
#'   variation in the data. Default of 4.
#' @param peakIntervalWidth Length-one numeric vector between \code{0} and
#'   \code{1} specifying desired width of linkage region(s). The default value
#'   of \code{0.95}, for example, yields peak regions defined as including the
#'   top 95\% of SNPs in the peak region, as determined by the peak
#'   resampling distribution.
#' @param loessOptResolution Length-one numeric vector between \code{0} and
#'   \code{1} specifying
#'   desired resolution for Loess fit optimization. The default of \code{0.001},
#'   for example, indicates that the span ultimately chosen will perform better
#'   than its neighbor values at \code{+-0.001}.
#' @param loessOptCutFactor Length-one numeric vector between \code{0} and
#'   \code{1} specifying how aggressively the Loess
#'   optimization algorithm proceeds. For example, with a default of \code{0.1}
#'   different spans at intervals of \code{0.001} would be evaluated after
#'   intervals of \code{0.01}.
#' @return A \code{MmapprParam} object.
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
NULL

mmapprParam <- function(wtFiles, 
                        mutFiles,
                        refFasta, 
                        gtf, 
                        outputFolder = 'DEFAULT', 
                        includeScaffolds = FALSE,
                        minDepth = 20,
                        homozygoteCutoff = 0.95, 
                        minBaseQuality = 20,
                        minMapQuality = 30,
                        fileAggregation = c('simple', 'weighted'), 
                        distancePower = 4,
                        peakIntervalWidth = 0.80, 
                        loessOptResolution = 0.001,
                        loessOptCutFactor = 0.1) {
    
    refFasta <- normalizePath(refFasta)
    genomeName <- tools::file_path_sans_ext(basename(refFasta), 
                                            compression = TRUE)
    gtf <- normalizePath(gtf)
    wtFiles <- normalizePath(wtFiles)
    mutFiles <- normalizePath(mutFiles)
    
    stopifnot(file.exists(refFasta), 
              file.exists(gtf), 
              file.exists(wtFiles), 
              file.exists(mutFiles))
    
    if (outputFolder == 'DEFAULT')
      outputFolder <- .defaultOutputFolder()
    outputFolder <- .prepareOutputFolder(outputFolder)
    
    # Need to prep bam files if no index present
    wtFiles <- .indexBamFileList(wtFiles, outputFolder)
    mutFiles <- .indexBamFileList(mutFiles, outputFolder)
    
    # Need to convert Fasta to gmapGenome
    indexFa(refFasta)
    refGenome <- gmapR::GmapGenome(refFasta, 
                                   name = genomeName,
                                   create = TRUE)
    
    # Need to prep gtf file if no index present
    if (!file.exists(paste0(gtf, ".tbi"))) {
      if (grepl("\\.gz$", gtf)) {
        system2("gunzip", gtf)
        gtf <- gsub("\\.gz$", "", gtf)
      } 
      system(paste0("grep -v '#' ", gtf, 
                    " | sort -k1,1 -k4,4n -k5,5n -t '\t'",
                    "| bgzip -c > ", gtf, ".bgz"))
      gtf <- paste0(gtf, ".bgz") 
      index = Rsamtools::indexTabix(gtf, 
                                    seq = 1, start = 4, end = 5, comment = "#")
    } else {
      .messageAndLog("Found GTF index, skipping index step.", outputFolder)
    }

    param <- new("MmapprParam", 
                 wtFiles = wtFiles, 
                 mutFiles = mutFiles,
                 refFasta = refFasta, 
                 refGenome = refGenome,
                 gtf = gtf, 
                 outputFolder = outputFolder, 
                 includeScaffolds = includeScaffolds,
                 minDepth = minDepth,
                 homozygoteCutoff = homozygoteCutoff, 
                 minBaseQuality = minBaseQuality,
                 minMapQuality = minMapQuality,
                 fileAggregation = fileAggregation[1], 
                 distancePower = distancePower,
                 peakIntervalWidth = peakIntervalWidth, 
                 loessOptResolution = loessOptResolution,
                 loessOptCutFactor = loessOptCutFactor)

    validity <- .validMmapprParam(param)
    if (typeof(validity) == "logical") {
      return(param)
    }
    else {
      stop(paste(validity, collapse = '\n  '))
    }
}


### VALIDITY FUNCTIONS

.validMmapprParam <- function(param) {
    errors <- character()

    .validityErrors <- function(fxn, value, errors) {
        result <- fxn(value)
        if (typeof(result) != 'logical')
            return(c(errors, result))
        else
            return(errors)
    }

    errors <- .validityErrors(.validFastaFile, refFasta(param), errors)
    errors <- .validityErrors(.validBamFiles, wtFiles(param), errors)
    errors <- .validityErrors(.validBamFiles, mutFiles(param), errors)

    if (length(errors) == 0) TRUE else errors
}


.validFastaFile <- function(filepath) {
    errors <- c()
    if (!file.exists(filepath))
        errors <- c(errors, paste(filepath, "does not exist"))
    if (length(errors) == 0) TRUE else errors
}


.validBamFiles <- function(files) {
    errors <- c()
    if (!is(files, 'BamFileList'))
        errors <- c(errors, paste0(files, " is not a BamFileList object"))
    for (i in seq_along(files)) {
        file <- files[[i]]
        if (!file.exists(file$path)) {
            errors <- c(errors, paste0(file$path, " does not exist"))
        }
    }
    if (length(errors) == 0) TRUE else errors
}

# Index Bams
.indexBamFileList <- function(bfl, oF) {
  indexed_bfl <- list()
  for (bam_file in bfl) {
    bam_path <- gsub("(.*/).*.bam$", "\\1", bam_file)
    bam_base <- gsub(".*/(.*).bam$", "\\1", bam_file)
    bam_index <- list.files(path = bam_path, 
                            pattern = paste0(bam_base, ".*bai$"), 
                            full.names = TRUE)
    if (length(bam_index) == 0) {
      .messageAndLog(paste(bam_file, "-- No index found. Indexing Bam now."), oF)
      bam_index <- indexBam(bam_file)
      indexed_bfl <- append(indexed_bfl, BamFile(bam_file, bam_index))
    } else {
      .messageAndLog(paste(bam_file, "-- Index found. Skipping index step."), oF)
      indexed_bfl <- append(indexed_bfl, BamFile(bam_file, bam_index[1]))
    }
  }
  BamFileList(indexed_bfl)
}


setMethod("show", "MmapprParam", function(object) {
    margin <- "   "
    cat("MmapprParam object with following values:\n")
    cat("Reference fasta file:\n", sep = "")
    cat(paste0(margin, object@refFasta, '\n'))
    cat("GTF file:\n", sep = "")
    cat(paste0(margin, object@gtf, '\n'))
    cat("wtFiles:\n", sep = "")
    .customPrint(object@wtFiles, margin)
    cat("mutFiles:\n", sep = "")
    .customPrint(object@mutFiles, margin)

    cat("Other parameters:\n")
    slotNames <- slotNames("MmapprParam")[c(6:16)]
    slotValues <- vapply(slotNames,
                         function(name, object) {
                           as.character(slot(object, name)[1])
                         },
                         FUN.VALUE = character(1), object)

    names(slotValues) <- slotNames
    print(slotValues, quote = FALSE)
})


.customPrint <- function(obj, margin = "  ", lineMax = getOption("max.print")) {
  lines <- capture.output(obj)
  lines <- strsplit(lines, split = "\n")
  lines <- vapply(lines, function(x) paste0(margin, x), character(1))
  if (lineMax > length(lines)) lineMax = length(lines)
  cat(lines[seq_len(lineMax)], sep = "\n")
}


#' MmapprParam Getters and Setters
#'
#' Access and assign slots of \code{\link{MmapprParam}} object.
#'
#' @name MmapprParam-functions
#' @aliases
#'   wtFiles wtFiles<-
#'   mutFiles mutFiles<-
#'   refFasta refFasta<-
#'   gtf gtf<-
#'   outputFolder outputFolder<-
#'   includeScaffolds includeScaffolds<-
#'   minDepth minDepth<-
#'   homozygoteCutoff homozygoteCutoff<-
#'   minBaseQuality minBaseQuality<-
#'   minMapQuality minMapQuality<-
#'   fileAggregation fileAggregation<-
#'   distancePower distancePower<-
#'   peakIntervalWidth peakIntervalWidth<-
#'   loessOptResolution loessOptResolution<-
#'   loessOptCutFactor loessOptCutFactor<-
#'
#' @param obj Desired \code{\link{MmapprParam}} object.
#' @param value Value to replace desired attribute.
#'
#' @return The desired \code{\link{MmapprParam}} attribute.
#'
#' @seealso \code{\link{MmapprParam}}
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                          mutFiles = MMAPPR2data::exampleMutBam(),
#'                          refFasta = MMAPPR2data::goldenFasta(),
#'                          gtf = MMAPPR2data::gtf())
#'
#'     outputFolder(mmappr_param) <- 'mmappr2_test_1'
#'     minBaseQuality(mmappr_param) <- 25
#'     vepFlags(mmappr_param)
#' }
NULL


### GETTERS
#' @rdname MmapprParam-functions
#' @export
setMethod("wtFiles", "MmapprParam", function(obj) obj@wtFiles)
#' @rdname MmapprParam-functions
#' @export
setMethod("mutFiles", "MmapprParam", function(obj) obj@mutFiles)
#' @rdname MmapprParam-functions
#' @export
setMethod("refFasta", "MmapprParam", function(obj) obj@refFasta)
#' @rdname MmapprParam-functions
#' @export
setMethod("gtf", "MmapprParam", function(obj) obj@gtf)
#' @rdname MmapprParam-functions
#' @export
setMethod("outputFolder", "MmapprParam", function(obj) obj@outputFolder)
#' @rdname MmapprParam-functions
#' @export
setMethod("includeScaffolds", "MmapprParam", function(obj) obj@includeScaffolds)
#' @rdname MmapprParam-functions
#' @export
setMethod("minDepth", "MmapprParam", function(obj) obj@minDepth)
#' @rdname MmapprParam-functions
#' @export
setMethod("homozygoteCutoff", "MmapprParam", function(obj) obj@homozygoteCutoff)
#' @rdname MmapprParam-functions
#' @export
setMethod("minBaseQuality", "MmapprParam", function(obj) obj@minBaseQuality)
#' @rdname MmapprParam-functions
#' @export
setMethod("minMapQuality", "MmapprParam", function(obj) obj@minMapQuality)
#' @rdname MmapprParam-functions
#' @export
setMethod("fileAggregation", "MmapprParam", function(obj) obj@fileAggregation)
#' @rdname MmapprParam-functions
#' @export
setMethod("distancePower", "MmapprParam", function(obj) obj@distancePower)
#' @rdname MmapprParam-functions
#' @export
setMethod("peakIntervalWidth", "MmapprParam", 
          function(obj) obj@peakIntervalWidth)
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptResolution", "MmapprParam", 
          function(obj) obj@loessOptResolution)
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptCutFactor", "MmapprParam", 
          function(obj) obj@loessOptCutFactor)

### SETTERS

#' @rdname MmapprParam-functions
#' @export
setMethod("wtFiles<-", "MmapprParam",
          function(obj, value) {
            obj@wtFiles <- Rsamtools::BamFileList(value)
            v <- .validBamFiles(obj@wtFiles)
            if (typeof(v) == 'logical') obj else stop(v)
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("mutFiles<-", "MmapprParam",
          function(obj, value) {
            obj@mutFiles <- Rsamtools::BamFileList(value)
            v <- .validBamFiles(obj@wtFiles)
            if (typeof(v) == 'logical') obj else stop(v)
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("refFasta<-", "MmapprParam",
          function(obj, value) {
            obj@refFasta <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("gtf<-", "MmapprParam",
          function(obj, value) {
            obj@gtf <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("outputFolder<-", "MmapprParam",
          function(obj, value) {
            obj@outputFolder <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("includeScaffolds<-", "MmapprParam",
          function(obj, value) {
            obj@includeScaffolds <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("minDepth<-", "MmapprParam",
          function(obj, value) {
            obj@minDepth <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("homozygoteCutoff<-", "MmapprParam",
          function(obj, value) {
            obj@homozygoteCutoff <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("minBaseQuality<-", "MmapprParam",
          function(obj, value) {
            obj@minBaseQuality <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("minMapQuality<-", "MmapprParam",
          function(obj, value) {
            obj@minMapQuality <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("fileAggregation<-", "MmapprParam",
          function(obj, value) {
            obj@fileAggregation <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("distancePower<-", "MmapprParam",
          function(obj, value) {
            obj@distancePower <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("peakIntervalWidth<-", "MmapprParam",
          function(obj, value) {
            obj@peakIntervalWidth <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptResolution<-", "MmapprParam",
          function(obj, value) {
            obj@loessOptResolution <- value
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptCutFactor<-", "MmapprParam",
          function(obj, value) {
            obj@loessOptCutFactor <- value
            obj
          })

