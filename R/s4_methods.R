MmapprParam <- function(refGenome, wtFiles, mutFiles, vepParam,
                        distancePower = 4, peakIntervalWidth = 0.95, minDepth = 10,
                        homozygoteCutoff = 0.8, numCores = 4, minBaseQuality = 20,
                        minMapQuality = 30, loessOptResolution = 0.001,
                        loessOptCutFactor = 0.1, naCutoff = 0, outputFolder = "DEFAULT") {
  
  if (class(wtFiles) == "character" | class(wtFiles) == "BamFile") wtFiles <- BamFileList(wtFiles)
  if (class(mutFiles) == "character" | class(mutFiles) == "BamFile") mutFiles <- BamFileList(mutFiles)

  if (outputFolder == "DEFAULT") outputFolder <- paste0("mmappr_results_", unclass(Sys.time()))
  
  param <- new("MmapprParam", refGenome = refGenome, wtFiles = wtFiles, mutFiles = mutFiles,
      vepParam = vepParam, distancePower = distancePower, minDepth = minDepth,
      homozygoteCutoff = homozygoteCutoff, numCores = numCores, minBaseQuality = minBaseQuality,
      minMapQuality = minMapQuality, loessOptResolution = loessOptResolution,
      loessOptCutFactor = loessOptCutFactor, naCutoff = naCutoff, outputFolder = outputFolder)
  
  validity <- .validMmapprParam(param)
  if (typeof(validity) == "logical") param else validity
}



.validMmapprParam <- function(param) {
  errors <- character()
  
  vepInputFormat <- input(param@vepParam)$format
  if (length(vepInputFormat) == 0) {
    msg <- "VEPParam requires 'vcf' input format"
    errors <- c(errors, msg)
  }
  else if (vepInputFormat != "vcf") {
    msg <- "VEPParam input format must be 'vcf'"
    errors <- c(errors, msg)
  }
  
  for (i in 1:length(param@wtFiles)) {
    file <- param@wtFiles[i]
    if (length(index(file)) == 0) {
      file <- .addBamFileIndex(file)
      if (length(index(file)) == 0) warning(paste0(file$path), " in wtFiles has no index file")
    }
  }
  for (i in 1:length(param@mutFiles)) {
    file <- param@wtFiles[i]
    if (length(index(file)) == 0) {
      file <- .addBamFileIndex(file)
      if (length(index(file)) == 0) warning(paste0(file$path), " in mutFiles has no index file")
    }
  }
  
  if (length(errors) == 0) TRUE else errors
}

setMethod("show", "MmapprParam", function(object) {
  margin <- "   "
  cat("MmapprParam object with following values:\n")
  cat("refGenome:\n")
  PrintWithMargin(object@refGenome, margin)
  cat("wtFiles:\n", sep="")
  PrintWithMargin(object@wtFiles, margin)
  cat("mutFiles:\n", sep="")
  PrintWithMargin(object@mutFiles, margin)
  cat("vepParam:\n", sep="")
  PrintWithMargin(object@vepParam, margin)
  
  cat("Other parameters:\n")
  slotNames <- names(
    getSlots("MmapprParam")[5:length(getSlots("MmapprParam"))])
  slotValues <- sapply(slotNames, function(name) slot(object, name))
  names(slotValues) <- slotNames
  print(slotValues, quote = FALSE)
})

setMethod("show", "MmapprData", function(object) {
  margin <- "  "
  cat("MmapprData object with following slots:\n")
  cat("param:\n")
  PrintWithMargin(object@param, margin)
  
  cat("distance:\n")
  classes <- sapply(object@distance, class)
  successes <- classes == "list"
  cat(margin, sprintf(
    "Contains Euclidian distance data for %i sequences\n", 
      sum(successes)), sep="")
  loessFits <- 0
  try(
    loessFits <- sum(sapply(object@distance[successes], function(seq) {
      if (!is.null(seq$loess)) return(TRUE)
      else return(FALSE)
    })), silent = TRUE
  )
  cat(margin, sprintf(
    "and Loess regression data for %i of those\n", loessFits
  ), sep="")
  distanceSize <- object.size(object@distance)
  cat(margin, sprintf("Memory usage = %.0f MB\n", 
                      distanceSize/1000000), sep = "")
  
  cat("peaks:\n")
  for (peak in object@peaks){
    cat(margin, sprintf("%s: start = %.0f, end = %.0f\n", 
                       peak$seqname, peak$start, peak$end),sep="")
    if (!is.null(peak$densityFunction))
      cat(margin, margin, "Density function calculated\n", sep="")
  }
  
  cat("candidates:\n")
  PrintWithMargin(object@candidates, margin)
})

PrintWithMargin<- function(obj, margin = "  ") {
  lines <- capture.output(obj)
  lines <- strsplit(lines, split = "\n")
  lines <- sapply(lines, function(x) paste0(margin, x))
  cat(lines, sep="\n")
}