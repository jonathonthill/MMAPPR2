MmapprParam <- function(refGenome, wtFiles, mutFiles, vepParam,
                        distancePower = 4, peakIntervalWidth = 0.95, minDepth = 10,
                        homozygoteCutoff = 0.8, numCores = 4, minBaseQuality = 20,
                        minMapQuality = 30, loessOptResolution = 0.001,
                        loessOptCutFactor = 0.1, naCutoff = 0, outputFolder = "DEFAULT") {
  
  if (class(wtFiles) == "character") wtFiles <- BamFileList(wtFiles)
  if (class(mutFiles) == "character") mutFiles <- BamFileList(mutFiles)
  if (outputFolder == "DEFAULT") outputFolder <- paste0("mmappr_results_", unclass(Sys.time()))
  
  new("MmapprParam", refGenome = refGenome, wtFiles = wtFiles, mutFiles = mutFiles,
      vepParam = vepParam, distancePower = distancePower, minDepth = minDepth,
      homozygoteCutoff = homozygoteCutoff, numCores = numCores, minBaseQuality = minBaseQuality,
      minMapQuality = minMapQuality, loessOptResolution = loessOptResolution,
      loessOptCutFactor = loessOptCutFactor, naCutoff = naCutoff, outputFolder = outputFolder)
  
}

CheckMmapprParam <- function(mmapprParam) {
  errors <- character()
  
  if (input(mmapprParam@vepParam)$format != "vcf") {
    msg <- "VEP param input format must be 'vcf'"
    errors <- c(errors, msg)
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

setMethod("refGenome", "MmapprParam", function(obj) obj@refGenome)
setMethod("wtFiles", "MmapprParam", function(obj) obj@wtFiles)
setMethod("mutFiles", "MmapprParam", function(obj) obj@mutFiles)
setMethod("homozygoteCutoff", "MmapprParam", function(obj) obj@homozygoteCutoff)
setMethod("vepParam", "MmapprParam", function(obj) obj@vepParam)
setMethod("distancePower", "MmapprParam", function(obj) obj@distancePower)
setMethod("peakIntervalWidth", "MmapprParam", function(obj) obj@peakIntervalWidth)
setMethod("minDepth", "MmapprParam", function(obj) obj@minDepth)
setMethod("numCores", "MmapprParam", function(obj) obj@numCores)
setMethod("minBaseQuality", "MmapprParam", function(obj) obj@minBaseQuality)
setMethod("minMapQuality", "MmapprParam", function(obj) obj@minMapQuality)
setMethod("loessOptResolution", "MmapprParam", function(obj) obj@loessOptResolution)
setMethod("loessOptCutFactor", "MmapprParam", function(obj) obj@loessOptCutFactor)
setMethod("naCutoff", "MmapprParam", function(obj) obj@naCutoff)
setMethod("outputFolder", "MmapprParam", function(obj) obj@outputFolder)

setMethod("refGenome<-", "MmapprParam",
          function(obj, value) {
            obj@refGenome <- value 
            obj
          })
setMethod("wtFiles<-", "MmapprParam",
          function(obj, value) {
            obj@wtFiles <- value 
            obj
          })
setMethod("mutFiles<-", "MmapprParam",
          function(obj, value) {
            obj@mutFiles <- value 
            obj
          })
setMethod("vepParam<-", "MmapprParam",
          function(obj, value) {
            obj@vepParam <- value 
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
setMethod("numCores<-", "MmapprParam",
          function(obj, value) {
            obj@numCores <- value 
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
            obj@outputFolder <- value 
            obj
          })




