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
    variants <- GetPeakVariants(g_ranges, mmapprData@param)
    
    #run VEP
    variants <- RunVEPForVariants(variants, mmapprData@param@vepParam)
    
    #filter out low impact variants
    variants <- FilterVariants(variants)
    
    #density score and order variants
    variants <- DensityScoreAndOrderVariants(variants, peak$density_function)
    
    return(variants)
  })
  
  return(mmapprData)
}

GetPeakVariants <- function(peak_granges, params){
  require(tidyr)
  require(dplyr)
  require(Rsamtools)
  require(VariantTools)
  
  result_vranges <- NULL
  
  # merge mutant bam files in desired regions
  if (length(params$mutFiles) < 2) mutBam <- params$mutFiles[[1]]
  else{
    mutBam <- mergeBam(params$mutFiles, destination = "tmp_m.bam", region = g_ranges)
  }
  #merge wt files
  if (length(params$wtFiles) < 2) wtBam <- params$wtFiles[[1]]
  else{
    wtBam <- mergeBam(params$wtFiles, destination = "tmp_wt.bam", region = g_ranges)
  }
  
  # create param for variant calling 
  tally_param <- TallyVariantsParam(genome = refGenome, 
                                    which = peak_granges,
                                    indels = TRUE,
                                    minimum_mapq = 1L
  )
  
  result_vranges <- (callSampleSpecificVariants(mutBam, wtBam, 
                                                tally.param = tally_param))
  
  if (file.exists("tmp_wt.bm")) file.remove("tmp_wt.bam")
  if (file.exists("tmp_m.bm")) file.remove("tmp_m.bam")
  
  if (length(result_vranges) > 0){
    # need sampleNames to convert to VCF; using mutant file names
    sampleNames(result_vranges) <- paste0(sapply(params$mutFiles, path), collapse = " -- ")
    mcols(result_vranges) <- NULL
    return(result_vranges)
  } 
  else return(NULL)
}

RunVEPForVariants <- function(inputVariants, vepParam){
  require(ensemblVEP)
  
  peakVcfFile <- "peak.vcf"
  
  #write file
  writeVcf(inputVariants, peakVcfFile)
  
  param <- vepParam
  
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
