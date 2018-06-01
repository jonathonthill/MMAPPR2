.getFileReadChrList <- function(mmapprData) {suppressWarnings({
    wtFiles <- mmapprData@param@wtFiles
    mutFiles <- mmapprData@param@mutFiles
    
    chrRanges <- as(GenomeInfoDb::seqinfo(Rsamtools::BamFileList(c(wtFiles, mutFiles))), "GRanges")
    #cut to standard chromosomes
    chrRanges <- GenomeInfoDb::keepStandardChromosomes(chrRanges, pruning.mode='coarse')
    chrRanges <- GenomeInfoDb::dropSeqlevels(chrRanges, 'chrM', pruning.mode='coarse')
    
    chrList <- list()
    # store range for each chromosome as list item
    for (i in suppressWarnings(GenomeInfoDb::orderSeqlevels(names(chrRanges)))) {
        chrList[[toString(GenomeInfoDb::seqnames(chrRanges[i]))]] <- chrRanges[i]
    }
    
    return(chrList)
})}


readInFiles <- function(mmapprData, showDebug=FALSE) {
    message("Reading BAM files")
    
    chrList <- suppressWarnings(.getFileReadChrList(mmapprData))
    
    mmapprData@distance <-
    .runFunctionInParallel(chrList, .readFilesForChr, param=mmapprData@param,
                           showDebug=showDebug)
    
    return(mmapprData)
}


.readFilesForChr <- function(chrRange, param, showDebug=FALSE){
    require(magrittr, quietly=TRUE)
    startTime <- proc.time()
    tryCatch({
        #parameter check
        stopifnot(length(GenomeInfoDb::seqnames(chrRange)) == 1)
        stopifnot(class(param) == "MmapprParam")
        
        stopifnot(!is.null(chrRange))
        wtFiles <- param@wtFiles
        mutFiles <- param@mutFiles
        
        #DEBUG
        # width(chrRange) <- 10000
        
        pf <- Rsamtools::PileupFiles(wtFiles)
        pf_mut <- Rsamtools::PileupFiles(mutFiles)
        
        #Function for use in applyPileups: gets list for each pileup position
        CalcInfo <-
            function(x)
            {
                ## x == information for each file
                nuc <- c('A', 'C', 'G', 'T')
                #this apply goes over each column of 'seq' array, which is 
                # each file
                #that means 'y' is a matrix of nucleotide by position
                info <- apply(x[["seq"]], 2, function(y) {
                    y <- y[c("A", "C", "G", "T"),,drop=FALSE]
                    cvg <- colSums(y)
                    #y <- y / cvg[col(y)] # normalizes reads: don't do that yet, 
                    # since we'll combine files later
                    #add coverage row that stores read depth
                    y <- rbind(y, cvg)
                    return(y)
                })
                #info is matrix of ACGTcvg by position
                list(pos=x[["pos"]], info=info) 
            }
        
        apParam <- Rsamtools::ApplyPileupsParam(which = chrRange, what="seq",
                                 minBaseQuality = param@minBaseQuality,
                                 minMapQuality = param@minMapQuality,
                                 minDepth = param@minDepth,
                                 maxDepth = 8000L,
                                 flag = Rsamtools::scanBamFlag(
                                     isSecondaryAlignment=FALSE))
        
        applyPileupWT <- Rsamtools::applyPileups(pf, FUN=CalcInfo, 
                                                 param=apParam)
        
        if (showDebug){
            for (ap in applyPileupWT)
                message(paste0("WT pileup with ", length(ap$pos), " rows"))
        }
        
        
        #make_df_for_chromsome: makes a function that turns each info list 
        # (from pileup) into dataframe
        #So it takes a list and returns a dataframe
        make_df_for_chromosome <- function(infoList){
                x <- data.frame('pos' = rep(infoList$pos, each = 5), 
                                'nucleotide' = c("A", "C", "G", "T", "cvg"))
                #gets right number of columns(files) by column bind
                #because $info is info by file, and info, which was originally
                #a matrix, gets simplified into a vector which matches
                x <- cbind(x, infoList$info)
                if (showDebug) message("after make_df_for_chr size is ", nrow(x))
                return(x)
        }
        
        
        ###depth filter
        #takes df with columns pos, nuc (ACGTcvg), file1, file2 ...
        #returns df with columns pos, nuc (ACGTcvg), file1, file2...
        #assumption: for each file, take out positions with reads under cutoff
        #accomplished by replacing position's data with NAs for each column(file)
        #NAs will then be dealt with appropriately by NA filter
        depth_filter <- function(chrDf){
            #gets only cvg rows (every 5th), then only info cols
            #returns true for undercovered rows
            #can only subset 2 dimensions with a matrix, not a df. Don't know why.
            filter_mat <- as.matrix(chrDf[seq(5, to=nrow(chrDf), by=5), 
                                          3:ncol(chrDf)] < param@minDepth)
            
            #expand filter_mat to cover all 5 rows for each position
            new_mat <- c()
            for (i in 1:ncol(filter_mat)){
                new_col <- rep(filter_mat[,i], each=5)
                new_mat <- cbind(new_mat, new_col)
            }
            filter_mat <- new_mat
            
            #apply filter to dataframe
            chrDf[, 3:ncol(chrDf)][filter_mat] <- NA
            
            if (showDebug) message("after depth_filter size is ", nrow(chrDf))
            return(chrDf)      
        }
        
        
        ###Filter for NAs
        #function we can use in lapply on results list, takes dataframe, 
        # returns it after filtering for NAs
        #format for df (input and output) should be pos (5 for each), 
        # nuc (ACGTcvg), file1, file2...
        na_filter <- function(chrDf){
            #if only one file and na_cutoff isn't 0
            if ( sum( !(names(chrDf) %in% c("pos", "nucleotide")) ) == 1 
                 & param@naCutoff != 0){
                warning('naCutoff for a single-file pool must be 0. Continuing with cutoff of 0.')
                param@naCutoff <- 0
            }
            #vector returns true on rows with cutoff or less NaNs or NAs (is.na accounts for both)
            #need drop=F so it works if only 1 column is passed to rowSums
            filter_vec <- (rowSums(is.na(chrDf[, 3:length(chrDf), drop = FALSE])) <= param@naCutoff)
            chrDf <- chrDf[filter_vec, ]
            if (showDebug) message("after na_filter size is ", nrow(chrDf))
            return(chrDf)
        }
        
        ###Get averages between files and divides by coverage
        #function takes dataframe with pos, nuc (ACGTcvg), file1, file2...returns it with files combined into one mean column
        #so output should be df with pos, nuc (ACGT), avg_count -- also removes coverage
        #by adding each file together then dividing by total coverage, we are essentially weighting each file by coverage
        avg_files <- function(chrDf){
            chrDf$avg_count <- rowSums(chrDf[,3:length(chrDf), drop = FALSE], na.rm = TRUE)
            #takes non-cvg rows and divides them by cvg
            chrDf[-seq(5, nrow(chrDf), by=5), "avg_count"] <- 
                chrDf[-seq(5, nrow(chrDf), by=5), "avg_count"] / chrDf[rep(seq(5, nrow(chrDf), by=5), each = 4), "avg_count"]
            
            #throw away non-mean columns and return
            if (showDebug) message("after avg_files size is ", nrow(chrDf))
            return(chrDf[, c("pos", "nucleotide", "avg_count")])
        }
        
        #filter homozygotes wt pool, so when later inner_joined,
        #rows will be taken out of both pools
        #removes rows where major allele frequency exceeds cutoff
        #takes df with cols pos, A, C, cvg, G, T; returns same
        homozygote_filter <- function(chrDf){
            rows_to_keep <- apply((chrDf[, -which(names(chrDf) %in% c("pos", "cvg"))] 
                                   <= param@homozygoteCutoff), MARGIN = 1, all)
            chrDf <- chrDf[rows_to_keep, ]
            if (showDebug) message("after hz_filter size is ", nrow(chrDf))
            return(chrDf)
        }
        
        
        #apply functions to wild type pool
        #CAUTION: functions must be applied in this order to work right
        tryCatch(
            wtCounts <- applyPileupWT[[1]] %>%
                make_df_for_chromosome() %>%
                depth_filter %>%
                na_filter() %>%
                avg_files() %>%
                tidyr::spread(key = nucleotide, value = avg_count) %>%
                #homoz filter only on wt pool
                homozygote_filter(),
            
            error = function(e) {
                msg <- 'Insufficient data in wild-type file(s)'
                stop(msg)
            }
        )
        colnames(wtCounts)[2:6] <- c('A.wt', 'C.wt', 'cvg.wt', 'G.wt', 'T.wt')
        rm(applyPileupWT)
        
        #apply functions to mutant pool
        applyPileupMut <- Rsamtools::applyPileups(pf_mut, FUN = CalcInfo, param = apParam)
        tryCatch(
            mutCounts <- applyPileupMut[[1]] %>%
                make_df_for_chromosome() %>%
                depth_filter %>%
                na_filter() %>%
                avg_files() %>%
                tidyr::spread(key = nucleotide, value = avg_count),
            
            error = function(e) {
                msg <- 'Insufficient data in mutant file(s)'
                stop(msg)
            }
        )
        
        colnames(mutCounts)[2:6] <- 
            c('A.mut', 'C.mut', 'cvg.mut', 'G.mut', 'T.mut')
        rm(applyPileupMut)
        
        #inner_join already removes rows without a match
        distanceDf <- dplyr::inner_join(wtCounts, mutCounts, by=c('pos'))
        if (length(distanceDf) == 0) 
            stop('Empty dataframe after joining WT and Mut count tables')
        
        #calculate Euclidian distance
        distanceDf$A.wt <- (distanceDf$A.wt - distanceDf$A.mut)^2
        distanceDf$C.wt <- (distanceDf$C.wt - distanceDf$C.mut)^2
        distanceDf$G.wt <- (distanceDf$G.wt - distanceDf$G.mut)^2
        distanceDf$T.wt <- (distanceDf$T.wt - distanceDf$T.mut)^2
        
        distanceDf <- dplyr::transmute(distanceDf, 
                                pos = pos,
                                distance = (A.wt + C.wt + G.wt + T.wt)^(1/2))
        distanceDf$distance <- distanceDf$distance ^ param@distancePower
        
        stopifnot(nrow(distanceDf) > 0)
        #print(proc.time() - startTime)
        
        resultList <- list(wtCounts = wtCounts, mutCounts = mutCounts, 
                           distanceDf = distanceDf)
        #debug: option to return pileup (memory / base pairs / reps) 
        # in order to predict memory
        readMemReqPerBP <- object.size(resultList) / nrow(wtCounts)
        resultList$readMemReqPerBP <- readMemReqPerBP  
        
        resultList$seqname <- as.character(GenomeInfoDb::seqnames(chrRange))
        
        return(resultList)
    },
    
    error = function(e) {
        msg <- paste0(toString(GenomeInfoDb::seqnames(chrRange)), ": ", e$message)
        str(msg)
        return(msg)
    }
    )
}
