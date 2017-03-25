# take peak region dataframe and run VEP for each region
# peak region df has cols chr, starts, and stops
ProcessPeakRegionDf <- function(peak_region_df){
  for (i in 1:nrow(peak_region_df)){
    pileup_param <- PileupParam(include_insertions = TRUE)
  }
}
