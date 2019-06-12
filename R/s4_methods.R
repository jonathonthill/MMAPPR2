#' @name MmapprParam
#'
#' @param refGenome \code{\link[gmapR]{GmapGenome}}
#'   object storing reference genome to be used in variant calling.
#'   Make sure it is the same genome aligned to and used installed with VEP.
#' @param wtFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the wild-type pool to be analyzed.
#' @param mutFiles Character vector,
#'   \code{\link[Rsamtools]{BamFile}}, or
#'   \code{\link[Rsamtools]{BamFileList}} containing
#'   BAM files for the mutant pool to be analyzed.
#' @param species Length-one character vector of name of species under
#'   analysis. Used only in generating default
#'   \code{\link[ensemblVEP]{VEPFlags}} object.
#' @param vepFlags Optional \code{\link[ensemblVEP]{VEPFlags}}
#'   object containing runtime options for Ensembl's Variant Effect Predictor.
#'   See vignette for details.
#' @param fasta The path to the fasta file genome, which will be referenced
#'   in reading the BAM files.
#' @param outputFolder Length-one character vector specifying where to save
#'   output, including a \code{\linkS4class{MmapprData}} stored as
#'   \code{mmappr_data.RDS}, \code{mmappr2.log}, a \code{.tsv} file
#'   for each peak chromosome containing candidate mutations, and PDF plots
#'   of both the entire genome and peak chromosomes. Defaults to an
#'   automatically generated \code{mmappr2_<timestamp>}.
#' @param distancePower Length-one numeric vector determing to what power
#'   Euclidean distance values are raised before fitting. Higher powers tend
#'   to increase high values and decrease low values, exaggerating the
#'   variation in the data. Default of 4.
#' @param peakIntervalWidth Length-one numeric vector between \code{0} and
#'   \code{1} specifying desired width of linkage region(s). The default value
#'   of \code{0.95}, for example, yields peak regions defined as including the
#'   top 95\% of SNPs in the peak region, as determined by the peak
#'   resampling distribution.
#' @param minDepth Length-one integer vector determining minimum depth
#'   required for a position to
#'   be considered in the analysis. Defaults to 10.
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
#' @param naCutoff Integer specifying the most NAs to accept at a given
#'   position--that is, the number of files without data for that position.
#'   Defaults to 0.
#' @param fileAggregation A length-one character vector determining strategy
#'   for aggregating base calls when multiple wild-type or multiple mutant
#'   files are provided.
#'   When 'sum', average base call proportions are computed using
#'   the read counts from each file, effectively weighting files
#'   with higher counts at a given position. When equal to 'mean', the
#'   base call proportions as well as read depths, rather than the absolute count,
#'   are averaged across files, which is useful when you want to weight each
#'   replicate evenly without regards to differing depth.
#'
#' @return A \code{MmapprParam} object.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE) &
#'         Sys.which('vep') != '') {
#'     genDir <- gmapR::GmapGenomeDirectory(tempdir(), create=TRUE)
#'     
#'     mmapprParam <- MmapprParam(refGenome = gmapR::GmapGenome("GRCz11", genDir),
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio")
#' }
MmapprParam <- function(refGenome, wtFiles, mutFiles, species, vepFlags=NULL,
                        fasta, outputFolder=NULL, distancePower=4,
                        peakIntervalWidth=0.95, minDepth=10,
                        homozygoteCutoff=0.95, minBaseQuality=20,
                        minMapQuality=30, loessOptResolution=0.001,
                        loessOptCutFactor=0.1, naCutoff=0, 
                        fileAggregation=c('sum', 'mean')) {
    
    wtFiles <- Rsamtools::BamFileList(wtFiles)
    mutFiles <- Rsamtools::BamFileList(mutFiles)
    
    # TODO: fasta
    
    if (is.null(vepFlags))
        vepFlags <- ensemblVEP::VEPFlags(flags = list(
            format = 'vcf',  # <-- this is necessary
            vcf = FALSE,  # <-- as well as this
            species = species,
            database = FALSE,
            cache = TRUE,
            filter_common = TRUE,
            coding_only = TRUE  # assuming RNA-seq data
        ))    
    
    if (is.null(outputFolder)) outputFolder <- 'DEFAULT'
    
    param <- new("MmapprParam", refGenome=refGenome, wtFiles=wtFiles, 
                 mutFiles=mutFiles, species=species, vepFlags=vepFlags, 
                 fasta=fasta, distancePower=distancePower,
                 peakIntervalWidth=peakIntervalWidth,
                 minDepth=minDepth,
                 homozygoteCutoff=homozygoteCutoff,
                 minBaseQuality=minBaseQuality, 
                 minMapQuality=minMapQuality,
                 loessOptResolution=loessOptResolution,
                 loessOptCutFactor=loessOptCutFactor, naCutoff=naCutoff, 
                 outputFolder=outputFolder,
                 fileAggregation=match.arg(fileAggregation))
    
    validity <- .validMmapprParam(param)
    if (typeof(validity) == "logical") param
    else stop(paste(validity, collapse='\n  '))
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

    errors <- .validityErrors(.validFastaFile, fasta(param), errors)
    errors <- .validityErrors(.validBamFiles, wtFiles(param), errors)
    errors <- .validityErrors(.validBamFiles, mutFiles(param), errors)
    errors <- .validityErrors(.validVepFlags, vepFlags(param), errors)

    
    if (length(errors) == 0) TRUE else errors
}

.validFastaFile <- function(filepath) {
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

.validVepFlags <- function(vepFlags) {
    vepFormat <- ensemblVEP::flags(vepFlags)$format
    # makes next conditional statement work:
    if (is.null(vepFormat)) vepFormat <- ""
    if (vepFormat != 'vcf'){
        return(paste0("VEPFlags format flag must be 'vcf'\n",
                      "  e.g., flags(vepFlags)$format <- 'vcf'"))
    } 
    return(TRUE)
}


setMethod("show", "MmapprParam", function(object) {
    margin <- "   "
    cat("MmapprParam object with following values:\n")
    cat("refGenome:\n")
    .customPrint(object@refGenome, margin)
    cat("wtFiles:\n", sep="")
    .customPrint(object@wtFiles, margin)
    cat("mutFiles:\n", sep="")
    .customPrint(object@mutFiles, margin)
    cat("vepFlags:\n", sep="")
    .customPrint(object@vepFlags, margin)
    cat("Reference fasta file:\n", sep="")
    .customPrint(as.character(object@fasta), margin)
    
    cat("Other parameters:\n")
    slotNames <- slotNames("MmapprParam")[7:length(slotNames("MmapprParam"))]
    slotNames <- c('species', slotNames)
    slotValues <- vapply(slotNames,
                         function(name) as.character(slot(object, name)),
                         FUN.VALUE=character(1)
    )
    names(slotValues) <- slotNames
    print(slotValues, quote=FALSE)
})

setMethod("show", "MmapprData", function(object) {
    margin <- "  "
    cat("MmapprData object with following slots:\n")
    cat("param:\n")
    .customPrint(object@param, margin)
    
    cat("distance:\n")
    classes <- vapply(object@distance, class, character(1))
    successes <- classes == "list"
    cat(margin, sprintf(
        "Contains Euclidian distance data for %i sequence(s)\n", 
        sum(successes)), sep="")
    loessFits <- 0
    try({loessFits <- sum(vapply(object@distance[successes],
                                 FUN.VALUE=logical(1),
                                 FUN=function(seq) {
                                     if (!is.null(seq$loess)) return(TRUE)
                                     else return(FALSE)
                                 }
                          ))
        }, silent=TRUE
    )
    cat(margin, sprintf(
        "and Loess regression data for %i of those\n", loessFits
    ), sep="")
    distanceSize <- object.size(object@distance)
    cat(margin, sprintf("Memory usage = %.0f MB\n", 
                        distanceSize/1000000), sep="")
    
    cat("peaks:\n")
    for (peak in object@peaks){
        cat(margin, sprintf('%s: ', peak$seqname), sep='')
        # this will fail and skip if start and end aren't calculated
        cat(sprintf('start = %.0f, end = %.0f',
                    peak$start, peak$end), sep="")
        cat('\n')
        if (!is.null(peak$densityFunction))
            cat(margin, margin, "Density function calculated\n", sep="")
    }
    
    cat("candidates:\n")
    for (i in seq_along(object@candidates))
        .customPrint(object@candidates[i], margin, lineMax=5)
})

.customPrint <- function(obj, margin="  ", lineMax=getOption("max.print")) {
  lines <- capture.output(obj)
  lines <- strsplit(lines, split="\n")
  lines <- vapply(lines, function(x) paste0(margin, x), character(1))
  if (lineMax > length(lines)) lineMax=length(lines)
  cat(lines[seq_len(lineMax)], sep="\n")
}


### GETTERS
#' MmapprParam Getters and Setters
#' 
#' Access and assign slots of \code{\link{MmapprParam}} object.
#' 
#' @name MmapprParam-functions
#' @aliases fileAggregation fileAggregation<-
#'   refGenome refGenome<-
#'   wtFiles wtFiles<-
#'   mutFiles mutFiles<-
#'   species species<-
#'   vepFlags vepFlags<-
#'   fasta fasta<-
#'   homozygoteCutoff homozygoteCutoff<-
#'   distancePower distancePower<-
#'   peakIntervalWidth peakIntervalWidth<-
#'   minDepth minDepth<-
#'   minBaseQuality minBaseQuality<-
#'   minMapQuality minMapQuality<-
#'   loessOptResolution loessOptResolution<-
#'   loessOptCutFactor loessOptCutFactor<-
#'   naCutoff naCutoff<-
#'   outputFolder outputFolder<-
#'   fileAggregation fileAggregation<-
#'
#' @param obj Desired \code{\link{MmapprParam}} object.
#' @param value Value to replace desired attribute.
#' 
#' @return The desired \code{\link{MmapprParam}} attribute.
#'   
#' @seealso \code{\link{MmapprParam}}
#' 
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE) &
#'         Sys.which('vep') != '') {
#'     genDir <- gmapR::GmapGenomeDirectory(tempdir(), create=TRUE)
#'     
#'     param <- MmapprParam(refGenome = gmapR::GmapGenome("GRCz11", genDir),
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio")
#'
#'     outputFolder(param) <- 'mmappr2_test_1'
#'     minBaseQuality(param) <- 25
#'     vepFlags(param)
#' }

NULL

#' @rdname MmapprParam-functions
#' @export
setMethod("refGenome", "MmapprParam", function(obj) obj@refGenome)
#' @rdname MmapprParam-functions
#' @export
setMethod("wtFiles", "MmapprParam", function(obj) obj@wtFiles)
#' @rdname MmapprParam-functions
#' @export
setMethod("mutFiles", "MmapprParam", function(obj) obj@mutFiles)
#' @rdname MmapprParam-functions
#' @export
setMethod("species", "MmapprParam", function(obj) obj@species)
#' @rdname MmapprParam-functions
#' @export
setMethod("vepFlags", "MmapprParam", function(obj) obj@vepFlags)
#' @rdname MmapprParam-functions
#' @export
setMethod("fasta", "MmapprParam", function(obj) obj@fasta)
#' @rdname MmapprParam-functions
#' @export
setMethod("homozygoteCutoff", "MmapprParam", function(obj) obj@homozygoteCutoff)
#' @rdname MmapprParam-functions
#' @export
setMethod("distancePower", "MmapprParam", function(obj) obj@distancePower)
#' @rdname MmapprParam-functions
#' @export
setMethod("peakIntervalWidth", "MmapprParam", function(obj) obj@peakIntervalWidth)
#' @rdname MmapprParam-functions
#' @export
setMethod("minDepth", "MmapprParam", function(obj) obj@minDepth)
#' @rdname MmapprParam-functions
#' @export
setMethod("minBaseQuality", "MmapprParam", function(obj) obj@minBaseQuality)
#' @rdname MmapprParam-functions
#' @export
setMethod("minMapQuality", "MmapprParam", function(obj) obj@minMapQuality)
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptResolution", "MmapprParam", function(obj) obj@loessOptResolution)
#' @rdname MmapprParam-functions
#' @export
setMethod("loessOptCutFactor", "MmapprParam", function(obj) obj@loessOptCutFactor)
#' @rdname MmapprParam-functions
#' @export
setMethod("naCutoff", "MmapprParam", function(obj) obj@naCutoff)
#' @rdname MmapprParam-functions
#' @export
setMethod("outputFolder", "MmapprParam", function(obj) obj@outputFolder)
#' @rdname MmapprParam-functions
#' @export
setMethod("fileAggregation", "MmapprParam", function(obj) obj@fileAggregation)

#' MmapprData Getters
#' 
#' Access slots of \code{\linkS4class{MmapprData}} object.
#' 
#' @param obj Desired \code{\linkS4class{MmapprData}} object.
#' 
#' @return Desired attribute.
#' 
#' @name MmapprData-getters
#' @aliases param distance peaks candidates
#' @seealso \code{\linkS4class{MmapprData}}
#' 
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE) &
#'         Sys.which('vep') != '') {
#'     genDir <- gmapR::GmapGenomeDirectory(tempdir(), create=TRUE)
#'     
#'     param <- MmapprParam(refGenome = gmapR::GmapGenome("GRCz11", genDir),
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio")
#'
#'     md <- new('MmapprData', param = param)
#' 
#'     param(md)
#'     distance(md)
#'     peaks(md)
#'     candidates(md)
#' }
NULL

#' @rdname MmapprData-getters
#' @export
setMethod("param", "MmapprData", function(obj) obj@param)
#' @rdname MmapprData-getters
#' @export
setMethod("distance", "MmapprData", function(obj) obj@distance)
#' @rdname MmapprData-getters
#' @export
setMethod("peaks", "MmapprData", function(obj) obj@peaks)
#' @rdname MmapprData-getters
#' @export
setMethod("candidates", "MmapprData", function(obj) obj@candidates)


### SETTERS

#' @rdname MmapprParam-functions
#' @export
setMethod("refGenome<-", "MmapprParam",
          function(obj, value) {
            obj@refGenome <- value 
            obj
          })
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
setMethod("vepFlags<-", "MmapprParam",
          function(obj, value) {
            obj@vepFlags <- value 
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("fasta<-", "MmapprParam",
          function(obj, value) {
            obj@fasta <- value 
            obj
          })
#' @rdname MmapprParam-functions
#' @export
setMethod("species<-", "MmapprParam",
          function(obj, value) {
            obj@species <- value 
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
setMethod("minDepth<-", "MmapprParam",
          function(obj, value) {
            obj@minDepth <- value 
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
#' @rdname MmapprParam-functions
#' @export
setMethod("naCutoff<-", "MmapprParam",
          function(obj, value) {
            obj@naCutoff <- value 
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