# take peak region dataframe and run VEP for each region
# peak region df has cols chr, starts, and stops

generateCandidates <- function(mmapprData) {
    
    #get GRanges representation of peak
    mmapprData@candidates <- lapply(mmapprData@peaks, .getPeakRange)
    
    #call variants in peak
    mmapprData@candidates <- lapply(mmapprData@candidates, FUN = .getVariantsForRange, 
                                    param = mmapprData@param)
    
    #run VEP
    mmapprData@candidates <- lapply(mmapprData@candidates, FUN = .runVEPForVariants,
                                    vepParam = mmapprData@param@vepParam)
    
    #filter out low impact variants
    mmapprData@candidates <- lapply(mmapprData@candidates, .filterVariants)
    
    #density score and order variants
    mmapprData@candidates <- lapply(names(mmapprData@candidates), function(seqname) {
        densityFunction <- mmapprData@peaks[[seqname]]$densityFunction
        variants <- mmapprData@candidates[[seqname]]
        variants <- DensityScoreAndOrderVariants(variants, densityFunction)
        return(variants)
    })
    
    #transfer names
    names(mmapprData@candidates) <- names(mmapprData@peaks)
    
    
    return(mmapprData)
}

.getPeakRange <- function(peakList) {
    ir <- IRanges(start = as.numeric(peakList$start),
                  end = as.numeric(peakList$end),
                  names = peakList$seqname)
    
    gr <- GRanges(seqnames = names(ir), 
                  ranges = ir)
    return(gr)
}

.getVariantsForRange <- function(inputRange, param){
    require(tidyr)
    require(dplyr)
    require(Rsamtools)
    require(VariantTools)
    
    # merge mutant bam files in desired regions
    if (length(param@mutFiles) < 2) mutBam <- param@mutFiles[[1]]
    else{
        mutBam <- mergeBam(param@mutFiles, destination = "tmp_m.bam", region = g_ranges)
    }
    #merge wt files
    if (length(param@wtFiles) < 2) wtBam <- param@wtFiles[[1]]
    else{
        wtBam <- mergeBam(param@wtFiles, destination = "tmp_wt.bam", region = g_ranges)
    }
    
    # create param for variant calling 
    tally_param <- TallyVariantsParam(genome = param@refGenome, 
                                      which = inputRange,
                                      indels = TRUE
    )
    
    resultVr <- (callSampleSpecificVariants(mutBam, wtBam, 
                                            tally.param = tally_param))
    
    if (file.exists("tmp_wt.bm")) file.remove("tmp_wt.bam")
    if (file.exists("tmp_m.bm")) file.remove("tmp_m.bam")
    
    if (length(resultVr) > 0){
        # need sampleNames to convert to VCF; using mutant file names
        sampleNames(resultVr) <- paste0(sapply(param@mutFiles, path), collapse = " -- ")
        mcols(resultVr) <- NULL
        return(resultVr)
    } 
    else return(NULL)
}

.runVEPForVariants <- function(inputVariants, vepParam){
    require(ensemblVEP)
    stopifnot(is(vepParam, "VEPParam"))
    
    peakVcfFile <- "peak.vcf"
    
    #write file
    writeVcf(inputVariants, peakVcfFile)
    
    param <- vepParam
    
    resultGRanges <- ensemblVEP(peakVcfFile, param)
    
    if (file.exists(peakVcfFile)) file.remove(peakVcfFile)
    
    #output granges
    return(resultGRanges)
}

.filterVariants <- function(candidate_granges) {
    filter <- elementMetadata(candidate_granges)$IMPACT != 'LOW'
    filter[is.na(filter)] <- TRUE
    return(candidate_granges[filter])
}

.densityScoreAndOrderVariants <- function(candidate_granges, density_function) {
    #density calculation
    positions <- start(candidate_granges) + ((width(candidate_granges) - 1) / 2)
    density_col <- sapply(positions, density_function)
    mcols(candidate_granges)$peakDensity <- density_col
    
    #re-order
    order_vec <- order(density_col, decreasing = TRUE)
    candidate_granges <- candidate_granges[order_vec]
    
    return(candidate_granges)
}
