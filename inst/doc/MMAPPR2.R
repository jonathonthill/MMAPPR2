# param <- mmapprParam(
#   wtFiles = c("wt_pool.bam"),
#   mutFiles = c("mutant_pool.bam"),
#   refFasta = "reference.fa",
#   gtf = "annotation.gtf",
#   outputFolder = "mmappr2_results"
# )

# param <- mmapprParam(
#   wtFiles = "wt_pool.bam",
#   mutFiles = "mutant_pool.bam",
#   refFasta = "reference.fa",
#   gtf = "annotation.gtf",
#   minDepth = 20,
#   minBaseQuality = 20,
#   minMapQuality = 20,
#   maxPileupDepth = 1000,
#   candidateMinDepth = 1,
#   candidateMinAltDepth = 2,
#   candidateMinAltFreq = 0.80,
#   # These WT-specific hard filters are optional by default. WT evidence is
#   # still reported for every candidate but does not change the default candidate ordering.
#   candidateMaxWtAltFreq = 1.00,
#   candidateMinDeltaAF = 0.00,
#   # Keep the fast one-shot pileup when it fits; otherwise chunk adaptively.
#   candidatePoolMode = "auto",
#   candidateChunkSize = 0,
#   peakCutoffSd = 3,
#   # Default: use the compatibility cutoff; "global_sd" selects the genome-wide-SD alternative.
#   peakCutoffMethod = "legacy_current",
#   peakIntervalWidth = 0.80,
#   # Conservative default: span all selected high-density modes.
#   peakIntervalMethod = "hpd_span",
#   peakResampleIterations = 1000,
#   randomSeed = 1,
#   pairedEnd = FALSE,
#   ignoreStrand = FALSE,
#   expressionPseudocount = 0.01,
#   # Export the evaluated AICc span curve for each fitted chromosome.
#   exportAiccPlots = FALSE
# )

# md <- mmappr(param)

# md <- mmapprData(param)
# md <- calculateDistance(md)  # pooled A/C/G/T frequencies and ED^p
# md <- loessFit(md)           # chromosome-wise robust LOESS
# md <- prePeak(md)            # genome-wide peak threshold
# md <- peakRefinement(md)     # reproducible half-marker resampling
# md <- generateCandidates(md) # mutant SNVs + WT evidence + coding effects
# outputMmapprData(md)

# md <- mmappr(param, palette = "colorblind",
#              colors = c(peakFill = "#E69F00AA", density = "#0072B2"))

# plotMmapprData("analysis_results/mmappr_data.RDS")
# plotMmapprData("analysis_results/mmappr_data.RDS",
#                outputDir = "publication_figures",
#                plots = c("genome", "peaks", "aicc"),
#                palette = "grayscale",
#                colors = c(fitted = "#222222", peakFill = "#AAAAAA"))
