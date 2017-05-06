setClass("MmapprData",
         representation(
           input="list",
           distance="list",
           peaks="list",
           candidates="list"
         ),
         prototype(
           input=list(
             wtFiles = wt_list,
             mutFiles = mut_list,
             refGenome = GmapGenome("danRer10"),
             #resolution at which AICc will be calculated to find optimum Loess fit span
             loessOptResolution = 0.01,
             #factor between rounds of Loess fit optimization (e.g., factor of 0.1 results in spans of 0.1 apart, then 0.01 apart, etc.)
             loessOptCutFactor = 0.1,
             minMapQuality = 30,
             minBaseQuality = 20,
             minDepth = 10,
             distancePower = 4,
             naCutoff = 0, # the most NAs we'll accept, that is, the number of files without data for that position
             homozygoteCutoff = 0.8 #the maximum WT allele frequency we'll accept for candidates
           )
         )
)