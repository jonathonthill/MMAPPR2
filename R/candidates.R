# take peak region dataframe and run VEP for each region
# peak region df has cols chr, starts, and stops

generateCandidates <- function(mmapprData) {
    
    #get GRanges representation of peak
    mmapprData@candidates <- lapply(mmapprData@peaks, .getPeakRange)
    
    #call variants in peak
    mmapprData@candidates <- lapply(mmapprData@candidates, FUN=.getVariantsForRange, 
                                    param=mmapprData@param)
    
    #run VEP
    mmapprData@candidates <- lapply(mmapprData@candidates, FUN=.runVEPForVariants,
                                    vepParam=mmapprData@param@vepParam)
    
    #filter out low impact variants
    mmapprData@candidates <- lapply(mmapprData@candidates, .filterVariants)
    
    #density score and order variants
    mmapprData@candidates <- lapply(names(mmapprData@candidates), function(seqname) {
        densityFunction <- mmapprData@peaks[[seqname]]$densityFunction
        stopifnot(!is.null(densityFunction))
        variants <- mmapprData@candidates[[seqname]]
        variants <- .densityScoreAndOrderVariants(variants, densityFunction)
        return(variants)
    })
    
    #transfer names
    names(mmapprData@candidates) <- names(mmapprData@peaks)
    
    
    return(mmapprData)
}

.getPeakRange <- function(peakList) {
    ir <- IRanges::IRanges(start=as.numeric(peakList$start),
                  end=as.numeric(peakList$end),
                  names=peakList$seqname)
    
    gr <- GenomicRanges::GRanges(seqnames=names(ir), 
                  ranges=ir)
    return(gr)
}

.getVariantsForRange <- function(inputRange, param) {
    # merge mutant bam files in desired regions
    if (length(param@mutFiles) < 2) mutBam <- param@mutFiles[[1]]
    else{
        mutBam <- Rsamtools::mergeBam(param@mutFiles, destination="tmp_m.bam", region=inputRange)
    }
    #merge wt files
    if (length(param@wtFiles) < 2) wtBam <- param@wtFiles[[1]]
    else{
        wtBam <- Rsamtools::mergeBam(param@wtFiles, destination="tmp_wt.bam", region=inputRange)
    }

    # create param for variant calling
    # lame fix for mockery double colon bug; normally use :: import
    require(VariantTools, quietly=TRUE)
    tallyParam <- TallyVariantsParam(genome=param@refGenome,
                                      which=inputRange,
                                      indels=TRUE
    )

    resultVr <- callSampleSpecificVariants(
        mutBam, wtBam, tally.param=tallyParam)

    if (file.exists("tmp_wt.bm")) file.remove("tmp_wt.bam")
    if (file.exists("tmp_m.bm")) file.remove("tmp_m.bam")

    if (length(resultVr) > 0) {
        # need sampleNames to convert to VCF; using mutant file names
        Biobase::sampleNames(resultVr) <-
            paste0(names(param@mutFiles),
                   collapse = " -- ")
        S4Vectors::mcols(resultVr) <- NULL
        return(resultVr)
    }
    else return(NULL)
}


.runVEPForVariants <- function(inputVariants, vepParam){
    stopifnot(is(vepParam, "VEPParam"))
    
    peakVcfFile <- "peak.vcf"
    
    #write file
    VariantAnnotation::writeVcf(inputVariants, peakVcfFile)
    
    param <- vepParam
    
    resultGRanges <- ensemblVEP::ensemblVEP(peakVcfFile, param)
    
    if (file.exists(peakVcfFile)) file.remove(peakVcfFile)
    
    #output GRanges
    return(resultGRanges)
}

.filterVariants <- function(candidateGRanges) {
    filter <-
        GenomicRanges::mcols(candidateGRanges)$IMPACT != 'LOW'
    filter[is.na(filter)] <- TRUE
    return(candidateGRanges[filter])
}

.densityScoreAndOrderVariants <- function(candidateGRanges, densityFunction) {
    #density calculation
    positions <- BiocGenerics::start(candidateGRanges) + 
        ((BiocGenerics::width(candidateGRanges) - 1) / 2)
    densityCol <- sapply(positions, densityFunction)
    GenomicRanges::mcols(candidateGRanges)$peakDensity <- densityCol
    
    #re-order
    orderVec <- order(densityCol, decreasing=TRUE)
    candidateGRanges <- candidateGRanges[orderVec]
    
    return(candidateGRanges)
}
