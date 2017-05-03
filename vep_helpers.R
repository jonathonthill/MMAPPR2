# take peak region dataframe and run VEP for each region
# peak region df has cols chr, starts, and stops
WritePeakPileupFile <- function(peak_region_df, output_file = "peak_region_pileup.txt"){
  require(tidyr)
  require(dplyr)
  require(Rsamtools)
  
  result <- as.data.frame(t(rep(NA, 5)))
  names(result) <- c('seqnames', 'pos', 'ref_nuc', 'nucleotide', 'count')
  
  # for each region, merge bam file, get pileup version and append to pileup table
  for (i in 1:nrow(peak_region_df)){
    # get desired region
    i_ranges <- IRanges(start = as.numeric(peak_region_df$starts[i]),
                        end = as.numeric(peak_region_df$stops[i]),
                        names = peak_region_df$chr[i])
    
    g_ranges <- GRanges(seqnames = names(i_ranges), 
                        ranges = i_ranges)
    scan_bam_param <- ScanBamParam(which = g_ranges)
    
    # merge mutant bam files in desired regions
    if (length(mut_list) < 2) input_bam <- mut_list[[1]]
    else{
      input_bam <- mergeBam(mut_list, destination = "tmp.bam", region = g_ranges)
    }
    
    pileup_param <- PileupParam(include_insertions = TRUE, min_base_quality = my_minBaseQuality,
                                min_nucleotide_depth = my_minDepth, min_mapq = my_minMapQuality)
    new_pileup_rows <- pileup(file = input_bam, scanBamParam = scan_bam_param,
                              pileupParam = pileup_param)
    new_pileup_rows <- transmute(new_pileup_rows, seqnames, pos, strand, nucleotide, count)
    new_pileup_rows$strand <- 'A'
    names(new_pileup_rows)[3] <- "ref_nuc"
    result <- rbind(result, new_pileup_rows)
  }
  
  result <- result[!is.na(result$seqnames), ]
  
  # for testing, cut table
  result <- result[1:1000, ]
  
  write.table(result, file = output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  if (file.exists("tmp.bm")) file.remove("tmp.bam")
  str(result)
}

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
