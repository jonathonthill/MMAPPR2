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

GetReferenceGenome <- function(genome_package_name) {
  require(gmapR)
  require(genome_package_name, character.only = TRUE)
  
  #TODO fix calling right genome based on name
  return(GmapGenome(genome_package_name, create = TRUE))
}

WritePeakVCF<- function(peak_region_df, output_file = "peak_region.vcf", 
                        genome_package_name = "BSgenome.Drerio.UCSC.danRer10"){
  require(tidyr)
  require(dplyr)
  require(Rsamtools)
  require(VariantTools)
  
  # reference_genome <- GetReferenceGenome(genome_package_name)
  reference_genome <- test_genome
  
  result_vcf <- NULL
  
  # for each region, merge bam files for wt and control
  for (i in 1:nrow(peak_region_df)){
    # get desired region
    i_ranges <- IRanges(start = as.numeric(peak_region_df$starts[i]),
                        end = as.numeric(peak_region_df$stops[i]),
                        names = peak_region_df$chr[i])
    
    g_ranges <- GRanges(seqnames = names(i_ranges), 
                        ranges = i_ranges)
    
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
                                      which = g_ranges + 1e4,
                                      indels = TRUE)
    
    # calling filters
    calling_filters <- VariantCallingFilters(read.count = my_minDepth)
    
    #TODO test call variants for case and control separately to see where problem is
    vcf_for_region <- callSampleSpecificVariants(mut_bam, wt_bam, 
                                               tally.param = tally_param,
                                               calling.filters = calling_filters)
    
    result_vcf <- rbind(result_vcf, vcf_for_region)
  }
  
  
  # write.table(result_vcf, file = output_file, sep = "\t", quote = FALSE, 
  #             row.names = FALSE, col.names = FALSE)
  if (file.exists("tmp_wt.bm")) file.remove("tmp_wt.bam")
  if (file.exists("tmp_m.bm")) file.remove("tmp_m.bam")
  return(result_vcf)
}

RunVEPForPileupFile <- function(peak_pileup_file = "peak_region.vcf"){
  #TODO put package in filename thing
  peak_pileup_file <- normalizePath(peak_pileup_file, mustWork = TRUE)
  
  require(ensemblVEP)
  
  param <- VEPParam(scriptPath =
                      "ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl",
                    input = c(species = 'danio_rerio_merged', format = "pileup"),
                    cache = c(cache = TRUE, offline = TRUE),
                    database = c(database = FALSE))
  
  gr <- ensemblVEP(peak_pileup_file, param)
  return(gr)
}
