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

GenerateCandidates <- function(mmappr_run) {
  #for each peak region
  for (peak in mmappr_run@peaks){
    #get GRanges representation of peak
    i_ranges <- IRanges(start = as.numeric(peak$start),
                        end = as.numeric(peak$end),
                        names = peak$chr)
    
    g_ranges <- GRanges(seqnames = names(i_ranges), 
                             ranges = i_ranges)
    
    #call variants in peak
    variants <- GetPeakVariants(g_ranges)
    
    #run VEP
    variants <- RunVEPForVCF(variants)
    
    #filter out low impact variants
    variants <- FilterVariants(variants)
    
    #density score and order variants
    variants <- DensityScoreAndOrderVariants(variants, peak$density_function)
    
    #add peak variants to result object
    mmappr_run@candidates <- append(mmappr_run@candidates, variants)
  }
  
  return(mmappr_run)
}

GetReferenceGenome <- function(genome_name) {
  require(gmapR)
  # require(genome_package_name, character.only = TRUE)
  
  #TODO fix calling right genome based on name
  result <- GmapGenome(genome_name, create = TRUE)
  return(result)
}

GetPeakVariants <- function(peak_granges, genome_name = "danRer10"){
  require(tidyr)
  require(dplyr)
  require(Rsamtools)
  require(VariantTools)
  
  reference_genome <- GetReferenceGenome(genome_name)
  
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
  tally_param <- TallyVariantsParam(genome = reference_genome, 
                                    which = peak_granges,
                                    indels = TRUE,
                                    minimum_mapq = 1L
  )
  
  # calling filters
  # calling_filters <- VariantCallingFilters(p.lower = 0.5)
  
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

RunVEPForVCF <- function(inputVCF){
  #TODO put package in filename thing
  peak_vcf_file <- "peak.vcf"
  
  #write file
  writeVcf(inputVCF, peak_vcf_file)
  
  require(ensemblVEP)
  
  param <- VEPParam(scriptPath =
                      "ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl",
                    input = c(species = 'danio_rerio_merged', format = "vcf"),
                    cache = c(cache = TRUE, offline = TRUE),
                    database = c(database = FALSE))
                    # dataformat = c(vcf = TRUE),
  
  gr <- ensemblVEP(peak_vcf_file, param)
  
  if (file.exists(peak_vcf_file)) file.remove(peak_vcf_file)
  
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
