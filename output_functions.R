OutputMmapprData <- function(mmapprData, plotAicc = FALSE) {
  if (class(mmapprData) != "MmapprData"){
    stop("Input object not of 'MmapprData' type")
  }
  
  PlotGenomeDistance(mmapprData@distance)
  PlotPeaks(mmapprData)
  if (plotAicc) PlotAicc(mmapprData@distance)
}

PlotGenomeDistance <- function(mmapprData) {
  
}