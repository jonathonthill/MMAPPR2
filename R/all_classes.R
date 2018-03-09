setClass("MmapprParam",
         representation(
             refGenome = "GmapGenome",
             wtFiles = "BamFileList",
             mutFiles = "BamFileList",
             vepFlags = "VEPFlags",
             distancePower = "numeric",
             peakIntervalWidth = "numeric",
             minDepth = "numeric",
             homozygoteCutoff = "numeric", #the maximum WT allele frequency we'll accept for candidates
             numCores = "numeric",
             minBaseQuality = "numeric",
             minMapQuality = "numeric",
             #resolution at which AICc will be calculated to find optimum Loess fit span
             loessOptResolution = "numeric",
             #factor between rounds of Loess fit optimization (e.g., factor of 0.1 results in spans of 0.1 apart, then 0.01 apart, etc.)
             loessOptCutFactor = "numeric",
             naCutoff = "numeric", # the most NAs we'll accept, that is, the number of files without data for that position
             outputFolder = "character"
         )
)



setClass("MmapprData",
         representation(
             param="MmapprParam",
             distance="list",
             peaks="list",
             candidates="list"
         )
)

