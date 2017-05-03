# take peak region dataframe and run VEP for each region
# peak region df has cols chr, starts, and stops

GenerateCandidates <- function(mmapprData) {
  
  #for each peak region
  mmapprData@candidates <- lapply(mmapprData@peaks, function(peak) {
    #get GRanges representation of peak
    i_ranges <- IRanges(start = as.numeric(peak$start),
                        end = as.numeric(peak$end),
                        names = peak$seqname)
    
    g_ranges <- GRanges(seqnames = names(i_ranges), 
                        ranges = i_ranges)
    
    #call variants in peak
    variants <- GetPeakVariants(g_ranges, mmapprData@input$refGenome)
    
    #run VEP
    variants <- RunVEPForVariants(variants)
    
    #filter out low impact variants
    variants <- FilterVariants(variants)
    
    #density score and order variants
    variants <- DensityScoreAndOrderVariants(variants, peak$density_function)
    
    return(variants)
  })
  
  return(mmapprData)
}

GetPeakVariants <- function(peak_granges, refGenome){
  require(tidyr)
  require(dplyr)
  require(Rsamtools)
  require(VariantTools)
  
  result_vranges <- NULL
  
  # merge mutant bam files in desired regions
  if (length(mut_list) < 2) mut_bam <- mut_list[[1]]
  else{
    mut_bam <- mergeBam(mut_list, destination = "tmp_m.bam", region = g_ranges)
  }
  #merge wt files
  if (length(wt_list) < 2) wt_bam <- wt_list[[1]]
  else{
    wt_bam <- mergeBam(wt_list, destination = "tmp_wt.bam", region = g_ranges)
  }
  
  # create param for variant calling 
  tally_param <- TallyVariantsParam(genome = refGenome, 
                                    which = peak_granges,
                                    indels = TRUE,
                                    minimum_mapq = 1L
  )
  
  result_vranges <- (callSampleSpecificVariants(mut_bam, wt_bam, 
                                                tally.param = tally_param))
  
  if (file.exists("tmp_wt.bm")) file.remove("tmp_wt.bam")
  if (file.exists("tmp_m.bm")) file.remove("tmp_m.bam")
  
  if (length(result_vranges) > 0){
    # need sampleNames to convert to VCF; using mutant file names
    sampleNames(result_vranges) <- paste0(sapply(mut_list, path), collapse = " -- ")
    mcols(result_vranges) <- NULL
    return(result_vranges)
  } 
  else return(NULL)
}

RunVEPForVariants <- function(inputVariants){
  require(ensemblVEP)
  
  peakVcfFile <- "peak.vcf"
  
  #write file
  writeVcf(inputVariants, peakVcfFile)
  
  param <- VEPParam(scriptPath =
                      "ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl",
                    input = c(species = 'danio_rerio_merged', format = "vcf"),
                    cache = c(cache = TRUE, offline = TRUE),
                    database = c(database = FALSE))
                    # dataformat = c(vcf = TRUE),
  
  gr <- ensemblVEP(peakVcfFile, param)
  
  if (file.exists(peakVcfFile)) file.remove(peakVcfFile)
  
  #output granges
  return(gr)
}

FilterVariants <- function(candidate_granges) {
  filter <- elementMetadata(candidate_granges)$IMPACT != 'LOW'
  return(candidate_granges[filter])
}

DensityScoreAndOrderVariants <- function(candidate_granges, density_function) {
  #density calculation
  positions <- start(candidate_granges) + ((width(candidate_granges) - 1) / 2)
  density_col <- sapply(positions, density_function)
  mcols(candidate_granges)$density <- density_col
  
  #re-order
  order_vec <- order(density_col, decreasing = TRUE)
  candidate_granges <- candidate_granges[order_vec]
  
  return(candidate_granges)
}
