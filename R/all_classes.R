library(ensemblVEP)
library(gmapR)

setClass("MmapprParam",
         representation(
           refGenome = "GmapGenome",
           wtFiles = "BamFileList",
           mutFiles = "BamFileList",
           vepParam = "VEPParam",
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
         ),
         prototype(
           vepParam = VEPParam(scriptPath =
                                 "ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl",
                               input = c(species = 'danio_rerio_merged', format = "vcf"),
                               cache = c(cache = TRUE, offline = TRUE),
                               database = c(database = FALSE)),
           distancePower = 4,
           peakIntervalWidth = 0.95,
           minDepth = 10,
           homozygoteCutoff = 0.8, 
           numCores = 4,
           minBaseQuality = 20,
           minMapQuality = 30,
           loessOptResolution = 0.01,
           loessOptCutFactor = 0.1,
           naCutoff = 0,
           outputFolder = "mmappr_result"
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

