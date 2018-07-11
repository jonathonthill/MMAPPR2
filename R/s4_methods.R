
MmapprParam <- function(refGenome, wtFiles, mutFiles, vepFlags,
                        outputFolder="DEFAULT", distancePower=4,
                        peakIntervalWidth=0.95, minDepth=10,
                        homozygoteCutoff=0.95, minBaseQuality=20,
                        minMapQuality=30, loessOptResolution=0.001,
                        loessOptCutFactor=0.1, naCutoff=0, 
                        fileAggregation=c('sum', 'mean')) {
    
    wtFiles <- Rsamtools::BamFileList(wtFiles)
    mutFiles <- Rsamtools::BamFileList(mutFiles)
    
    param <- new("MmapprParam", refGenome=refGenome, wtFiles=wtFiles, 
                 mutFiles=mutFiles, vepFlags=vepFlags, 
                 distancePower=distancePower, peakIntervalWidth=peakIntervalWidth,
                 minDepth=minDepth,
                 homozygoteCutoff=homozygoteCutoff,
                 minBaseQuality=minBaseQuality, 
                 minMapQuality=minMapQuality,
                 loessOptResolution=loessOptResolution,
                 loessOptCutFactor=loessOptCutFactor, naCutoff=naCutoff, 
                 outputFolder=outputFolder,
                 fileAggregation=match.arg(fileAggregation))
    
    validity <- .validMmapprParam(param)
    if (typeof(validity) == "logical") param else stop(paste(validity, collapse='\n  '))
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

    errors <- .validityErrors(.validFiles, param@wtFiles, errors)
    errors <- .validityErrors(.validFiles, param@mutFiles, errors)
    errors <- .validityErrors(.validVepFlags, param@vepFlags, errors)

    
    if (length(errors) == 0) TRUE else errors
}

.validFiles <- function(files) {
    errors <- c()
    if (class(files) != "BamFileList") 
        errors <- c(errors, paste0(files, " is not a BamFileList object"))
    for (i in seq_along(files)) {
        file <- files[[i]]
        if (file.exists(file$path)) {
            if (length(Rsamtools::index(file)) == 0) {
                file <- .addBamFileIndex(file)
                if (length(Rsamtools::index(file)) == 0)
                    warning(paste0(file$path), " in wtFiles has no index file")
            }
        } else {
            errors <- c(errors, paste0(file$path, " does not exist"))
        }
    }
    if (length(errors) == 0) TRUE else errors
}

.validVepFlags <- function(vepFlags) {
    vepFormat <- ensemblVEP::flags(vepFlags)$format
    if (is.null(vepFormat)) vepFormat <- "" # makes next conditional statement work
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
    
    cat("Other parameters:\n")
    slotNames <- slotNames("MmapprParam")[5:length(slotNames("MmapprParam"))]
    slotValues <- sapply(slotNames, function(name) slot(object, name))
    names(slotValues) <- slotNames
    print(slotValues, quote=FALSE)
})

setMethod("show", "MmapprData", function(object) {
    margin <- "  "
    cat("MmapprData object with following slots:\n")
    cat("param:\n")
    .customPrint(object@param, margin)
    
    cat("distance:\n")
    classes <- sapply(object@distance, class)
    successes <- classes == "list"
    cat(margin, sprintf(
        "Contains Euclidian distance data for %i sequence(s)\n", 
        sum(successes)), sep="")
    loessFits <- 0
    try(
        loessFits <- sum(sapply(object@distance[successes], function(seq) {
            if (!is.null(seq$loess)) return(TRUE)
            else return(FALSE)
        })), silent=TRUE
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
    for (i in seq_along(object@candidates)) .customPrint(object@candidates[i], margin, lineMax=5)
})

.customPrint <- function(obj, margin="  ", lineMax=getOption("max.print")) {
  lines <- capture.output(obj)
  lines <- strsplit(lines, split="\n")
  lines <- sapply(lines, function(x) paste0(margin, x))
  if (lineMax > length(lines)) lineMax=length(lines)
  cat(lines[1:lineMax], sep="\n")
}

setMethod("refGenome", "MmapprParam", function(obj) obj@refGenome)
setMethod("wtFiles", "MmapprParam", function(obj) obj@wtFiles)
setMethod("mutFiles", "MmapprParam", function(obj) obj@mutFiles)
setMethod("homozygoteCutoff", "MmapprParam", function(obj) obj@homozygoteCutoff)
setMethod("vepFlags", "MmapprParam", function(obj) obj@vepFlags)
setMethod("distancePower", "MmapprParam", function(obj) obj@distancePower)
setMethod("peakIntervalWidth", "MmapprParam", function(obj) obj@peakIntervalWidth)
setMethod("minDepth", "MmapprParam", function(obj) obj@minDepth)
setMethod("minBaseQuality", "MmapprParam", function(obj) obj@minBaseQuality)
setMethod("minMapQuality", "MmapprParam", function(obj) obj@minMapQuality)
setMethod("loessOptResolution", "MmapprParam", function(obj) obj@loessOptResolution)
setMethod("loessOptCutFactor", "MmapprParam", function(obj) obj@loessOptCutFactor)
setMethod("naCutoff", "MmapprParam", function(obj) obj@naCutoff)
setMethod("outputFolder", "MmapprParam", function(obj) obj@outputFolder)
setMethod("fileAggregation", "MmapprParam", function(obj) obj@fileAggregation)
setMethod("param", "MmapprData", function(obj) obj@param)
setMethod("distance", "MmapprData", function(obj) obj@distance)
setMethod("peaks", "MmapprData", function(obj) obj@peaks)
setMethod("candidates", "MmapprData", function(obj) obj@candidates)


setMethod("refGenome<-", "MmapprParam",
          function(obj, value) {
            obj@refGenome <- value 
            obj
          })
setMethod("wtFiles<-", "MmapprParam",
          function(obj, value) {
              obj@wtFiles <- Rsamtools::BamFileList(value)
              v <- .validFiles(obj@wtFiles)
              if (typeof(v) == 'logical') obj else v
          })
setMethod("mutFiles<-", "MmapprParam",
          function(obj, value) {
              obj@mutFiles <- Rsamtools::BamFileList(value)
              v <- .validFiles(obj@wtFiles)
              if (typeof(v) == 'logical') obj else v
          })
setMethod("vepFlags<-", "MmapprParam",
          function(obj, value) {
            obj@vepFlags <- value 
            obj
          })
setMethod("homozygoteCutoff<-", "MmapprParam",
          function(obj, value) {
            obj@homozygoteCutoff <- value 
            obj
          })
setMethod("distancePower<-", "MmapprParam",
          function(obj, value) {
            obj@distancePower <- value 
            obj
          })
setMethod("peakIntervalWidth<-", "MmapprParam",
          function(obj, value) {
            obj@peakIntervalWidth <- value 
            obj
          })
setMethod("minDepth<-", "MmapprParam",
          function(obj, value) {
            obj@minDepth <- value 
            obj
          })
setMethod("minBaseQuality<-", "MmapprParam",
          function(obj, value) {
            obj@minBaseQuality <- value 
            obj
          })
setMethod("loessOptResolution<-", "MmapprParam",
          function(obj, value) {
            obj@loessOptResolution <- value 
            obj
          })
setMethod("loessOptCutFactor<-", "MmapprParam",
          function(obj, value) {
            obj@loessOptCutFactor <- value 
            obj
          })
setMethod("naCutoff<-", "MmapprParam",
          function(obj, value) {
            obj@naCutoff <- value 
            obj
          })
setMethod("outputFolder<-", "MmapprParam",
          function(obj, value) {
            obj@outputFolder <- value 
            obj
          })
setMethod("minMapQuality<-", "MmapprParam",
          function(obj, value) {
            obj@minMapQuality <- value 
            obj
          })
setMethod("fileAggregation<-", "MmapprParam",
          function(obj, value) {
            obj@fileAggregation <- value 
            obj
          })