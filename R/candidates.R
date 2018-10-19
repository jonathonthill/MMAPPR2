#' Generate potential causative mutations and consequences in peak regions
#' 
#' Follows the \code{\link{peakRefinement}} step and produces a
#' \code{\linkS4class{MmapprData}} object ready for
#' \code{\link{outputMmapprData}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{candidates}
#'   slot filled with a \code{\link[GenomicRanges]{GRanges}} object for each
#'   peak chromosome containing variants and predicted consequences from
#'   Ensembl's Variant Effect Predictor.
#' @export
#'
#' @examples
#' \dontrun{
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' }
generateCandidates <- function(mmapprData) {
    
    #get GRanges representation of peak
    mmapprData@candidates <- lapply(mmapprData@peaks, .getPeakRange)
    
    #call variants in peak
    mmapprData@candidates <- lapply(mmapprData@candidates, FUN=.getVariantsForRange, 
                                    param=mmapprData@param)
    
    #run VEP
    mmapprData@candidates <- lapply(mmapprData@candidates, FUN=.runVEPForVariants,
                                    vepFlags=mmapprData@param@vepFlags)
    
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
    # merge files in desired region if there are multiple
    if (length(param@mutFiles) < 2) mutBam <- param@mutFiles[[1]]
    else{
        mutBam <- mergeBam(param@mutFiles, destination="tmp.bam", region=inputRange)
    }

    # create param for variant calling
    tallyParam <- TallyVariantsParam(genome=param@refGenome,
                                     which=inputRange,
                                     indels=TRUE
    )

    resultVr <- callVariants(mutBam, tally.param=tallyParam)

    if (file.exists("tmp.bam")) file.remove("tmp.bam")

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


.runVEPForVariants <- function(inputVariants, vepFlags){
    stopifnot(is(vepFlags, "VEPFlags"))
    stopifnot(is(inputVariants, 'VRanges'))
    
    peakVcfFile <- "/tmp/peak.vcf"
    tryCatch({
        VariantAnnotation::writeVcf(inputVariants, peakVcfFile)
        resultGRanges <- ensemblVEP::ensemblVEP(peakVcfFile, vepFlags)
    }, error=function(e) {
        stop(e)
    },finally={
        if (file.exists(peakVcfFile)) file.remove(peakVcfFile)
    })
    
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
    densityCol <- vapply(positions, densityFunction, FUN.VALUE=numeric(1))
    GenomicRanges::mcols(candidateGRanges)$peakDensity <- densityCol
    
    #re-order
    orderVec <- order(densityCol, decreasing=TRUE)
    candidateGRanges <- candidateGRanges[orderVec]
    
    return(candidateGRanges)
}
