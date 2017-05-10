library(ensemblVEP)

setClass("MmapprParam",
         representation(
           refGenome = "GmapGenome",
           wtFiles = "BamFileList",
           mutFiles = "BamFileList",
           vepParam = "VEPParam",
           distancePower = "numeric",
           minDepth = "numeric",
           homozygoteCutoff = "numeric", #the maximum WT allele frequency we'll accept for candidates
           numCores = "numeric",
           minBaseQuality = "numeric",
           minMapQuality = "numeric",
           #resolution at which AICc will be calculated to find optimum Loess fit span
           loessOptResolution = "numeric",
           #factor between rounds of Loess fit optimization (e.g., factor of 0.1 results in spans of 0.1 apart, then 0.01 apart, etc.)
           loessOptCutFactor = "numeric",
           naCutoff = "numeric" # the most NAs we'll accept, that is, the number of files without data for that position
         ),
         prototype(
           vepParam = VEPParam(scriptPath =
                                 "ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl",
                               input = c(species = 'danio_rerio_merged', format = "vcf"),
                               cache = c(cache = TRUE, offline = TRUE),
                               database = c(database = FALSE)),
           distancePower = 4,
           minDepth = 10,
           homozygoteCutoff = 0.8, 
           minBaseQuality = 20,
           minMapQuality = 30,
           loessOptResolution = 0.01,
           loessOptCutFactor = 0.1,
           naCutoff = 0
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

MmapprParam <- function(refGenome, wtFiles, mutFiles, vepParam,
                          distancePower = 4, minDepth = 10,
                          homozygoteCutoff = 0.8, numCores = 4, minBaseQuality = 20,
                          minMapQuality = 30, loessOptResolution = 0.01,
                          loessOptCutFactor = 0.1, naCutoff = 0) {
  
  if (class(wtFiles) == "character") wtFiles <- BamFileList(wtFiles)
  if (class(mutFiles) == "character") mutFiles <- BamFileList(mutFiles)
  
  new("MmapprParam", refGenome = refGenome, wtFiles = wtFiles, mutFiles = mutFiles,
      vepParam = vepParam, distancePower = distancePower, minDepth = minDepth,
      homozygoteCutoff = homozygoteCutoff, numCores = numCores, minBaseQuality = minBaseQuality,
      minMapQuality = minMapQuality, loessOptResolution = loessOptResolution,
      loessOptCutFactor = loessOptCutFactor, naCutoff = naCutoff)
  
}

CheckMmapprParam <- function(mmapprParam) {
  errors <- character()
  
  if (input(mmapprParam@vepParam)$format != "vcf") {
    msg <- "VEP param input format must be 'vcf'"
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else 0
}