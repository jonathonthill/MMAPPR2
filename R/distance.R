#' Read BAM files and generate Euclidean distance data
#'
#' First step in the MMAPPR2 pipeline. Precedes the \code{\link{loessFit}}
#' step.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{distance}
#'   slot filled.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)
#'         & all(Sys.which(c("samtools", "vep")) != "") {
#'     mmappr_param <- MmapprParam(refFasta = MMAPPR2data::goldenFasta(),
#'                                wtFiles = MMAPPR2data::exampleWTbam(),
#'                                mutFiles = MMAPPR2data::exampleMutBam(),
#'                                species = "danio_rerio",
#'                                outputFolder = tempOutputFolder())
#'
#'     md <- new('MmapprData', param = mmappr_param)
#'     postCalcDistMD <- calculateDistance(md)
#' }

calculateDistance <- function(mmapprData) {
    .indexBamFileList(wtFiles(param(mmapprData)))
    .indexBamFileList(mutFiles(param(mmapprData)))

    chrList <- suppressWarnings(.getFileReadChrList(mmapprData))

    mmapprData@distance <-
        BiocParallel::bplapply(chrList, .calcDistForChr, param = mmapprData@param)

    return(mmapprData)
}


.indexBamFileList <- function(bfl) {
    for (i in seq_along(bfl)) {
        bamFile <- bfl[[i]]
        if (is.na(Rsamtools::index(bamFile)))
            Rsamtools::indexBam(bamFile)
    }
}


.calcDistForChr <- function(chrRange, param){
    startTime <- proc.time()
    tryCatch({
        #parameter check
        stopifnot(length(GenomeInfoDb::seqnames(chrRange)) == 1)
        stopifnot(is(param, "MmapprParam"))

        stopifnot(!is.null(chrRange))

        #apply functions to wild type pool
        #CAUTION: functions must be applied in this order to work right
        tryCatch({
            pileupWT <- .samtoolsPileup(files = wtFiles(param), param, chrRange)

            wtCounts <- pileupWT %>%
                .naFilter(naCutoff = naCutoff(param)) %>%
                .avgFiles(fileAggregation = fileAggregation(param)) %>%
                tidyr::spread(key = 'nucleotide', value = 'avgCount') %>%
                #homoz filter only on wt pool
                .homozygoteFilter(homozygoteCutoff=homozygoteCutoff(param))

            }, error = function(e) {
                msg <- 'Insufficient data in wild-type file(s)'
                print(e)
                stop(msg)
            }
        )
        colnames(wtCounts)[2:6] <- c('A.wt', 'C.wt', 'cvg.wt', 'G.wt', 'T.wt')
        rm(pileupWT)
        gc()


        tryCatch({
            pileupMut <- .samtoolsPileup(files = mutFiles(param), param, chrRange)
            mutCounts <- pileupMut %>%
                .naFilter(naCutoff = naCutoff(param)) %>%
                .avgFiles(fileAggregation = fileAggregation(param)) %>%
                tidyr::spread(key = 'nucleotide', value = 'avgCount')

            }, error = function(e) {
                msg <- 'Insufficient data in mutant file(s)'
                stop(msg)
            }
        )

        colnames(mutCounts)[2:6] <-
            c('A.mut', 'C.mut', 'cvg.mut', 'G.mut', 'T.mut')
        rm(pileupMut)
        gc()

        #inner_join already removes rows without a match
        distanceDf <- dplyr::inner_join(wtCounts, mutCounts, by=c('pos'))
        if (length(distanceDf) == 0)
            stop('Empty dataframe after joining WT and Mut count tables')

        #calculate Euclidian distance
        distanceDf$A <- (distanceDf$A.wt - distanceDf$A.mut)^2
        distanceDf$C <- (distanceDf$C.wt - distanceDf$C.mut)^2
        distanceDf$G <- (distanceDf$G.wt - distanceDf$G.mut)^2
        distanceDf$T <- (distanceDf$T.wt - distanceDf$T.mut)^2

        distanceDf <- dplyr::transmute(distanceDf,
                                       pos = distanceDf$pos,
                                       distance =
                                           (distanceDf$A +
                                                distanceDf$C +
                                                distanceDf$G +
                                                distanceDf$T)^(1/2))
        distanceDf$distance <- distanceDf$distance ^ param@distancePower

        stopifnot(nrow(distanceDf) > 0)

        resultList <- list(wtCounts = wtCounts, mutCounts = mutCounts,
                           distanceDf = distanceDf)
        resultList$seqname <- as.character(GenomeInfoDb::seqnames(chrRange))

        return(resultList)
    },

    error = function(e) {
        msg <- paste(toString(GenomeInfoDb::seqnames(chrRange)),
                     e$message,
                     sep = ': ')
        return(msg)
    }
    )
}

.getFileReadChrList <- function(mmapprData) {suppressWarnings({
    bams <- Rsamtools::BamFileList(c(mmapprData@param@wtFiles,
                                     mmapprData@param@mutFiles))

    bamInfo <- Rsamtools::seqinfo(bams)
    chrRanges <- as(bamInfo, "GRanges")
    #cut to standard chromosomes
    chrRanges <- GenomeInfoDb::keepStandardChromosomes(chrRanges,
                                                       pruning.mode = 'coarse')
    chrRanges <- GenomeInfoDb::dropSeqlevels(chrRanges, 'chrM',
                                             pruning.mode = 'coarse')
    chrRanges <- GenomeInfoDb::dropSeqlevels(chrRanges, 'MT',
                                             pruning.mode = 'coarse')

    chrList <- list()
    # store range for each chromosome as list item
    for (i in suppressWarnings(
        GenomeInfoDb::orderSeqlevels(
            as.character(GenomeInfoDb::seqnames(chrRanges))))) {

        chrList[[toString(GenomeInfoDb::seqnames(chrRanges[i]))]] <-
            chrRanges[i]
    }


    return(chrList)
})}

###Filter for NAs
#function we can use in lapply on results list, takes dataframe,
# returns it after filtering for NAs
#format for df (input and output) should be pos (5 for each),
# nuc (ACGTcvg), file1, file2...
.naFilter <- function(chrDf, naCutoff){
    #if only one file and na_cutoff isn't 0
    if ( sum( !(names(chrDf) %in% c("pos", "nucleotide")) ) == 1
         & naCutoff != 0){
        warning(paste('naCutoff for a single-file pool must be 0.',
                      'Continuing with cutoff of 0.'))
        naCutoff <- 0
    }
    # vector returns true on rows with cutoff or less NaNs or NAs
    # (is.na accounts for both)
    # need drop=F so it works if only 1 column is passed to rowSums
    filter_vec <-
        (rowSums(is.na(chrDf[, 3:length(chrDf), drop = FALSE])) <= naCutoff)
    chrDf <- chrDf[filter_vec, ]
    return(chrDf)
}

###Get averages between files and divides by coverage
#function takes dataframe with pos, nuc (ACGTcvg--with counts),
#  file1, file2...returns it with files combined into one mean column
#so output should be df with pos, nuc (ACGTcvg--with base proportions), avgCount
.avgFiles <- function(chrDf, fileAggregation="sum") {
    stopifnot(fileAggregation %in% c('sum', 'mean'))
    #throw away non-mean columns and return
    if (fileAggregation == 'sum') {
        chrDf$avgCount <-
            rowSums(chrDf[, 3:length(chrDf), drop=FALSE], na.rm=TRUE)
        #takes non-cvg rows and divides them by cvg
        chrDf[-seq(5, nrow(chrDf), by=5), 'avgCount'] <-
            chrDf[-seq(5, nrow(chrDf), by=5), 'avgCount'] /
            chrDf[rep(seq(5, nrow(chrDf), by=5), each = 4), 'avgCount']
    }
    else if (fileAggregation == 'mean') {
        #takes non-cvg rows and divides them by cvg
        chrDf[-seq(5, nrow(chrDf), by=5), 3:length(chrDf)] <-
            chrDf[-seq(5, nrow(chrDf), by=5), 3:length(chrDf)] /
            chrDf[rep(seq(5, nrow(chrDf), by=5), each = 4), 3:length(chrDf)]
        chrDf$avgCount <-
            rowMeans(chrDf[, 3:length(chrDf), drop=FALSE], na.rm=TRUE)
    }
    return(chrDf[, c("pos", "nucleotide", 'avgCount')])
}

#filter homozygotes wt pool, so when later inner_joined,
#rows will be taken out of both pools
#removes rows where major allele frequency exceeds cutoff
#takes df with cols pos, A, C, cvg, G, T; returns same
.homozygoteFilter <- function(chrDf, homozygoteCutoff){
    rows_to_keep <- apply((chrDf[, -which(names(chrDf) %in% c("pos", "cvg"))]
                           <= homozygoteCutoff), MARGIN = 1, all)
    chrDf <- chrDf[rows_to_keep, ]
    return(chrDf)
}
