###this whole script should return a list with two dataframes
# one of the Euclidian distance for a chromosome
# and another of the original proportions and distance
# arguments: GRanges object with desired chromosome, (cutoffs?, )


my_minMapQuality <- 30
my_minBaseQuality <- 20
my_minDepth <- 10
my_power <- 4
na_cutoff <- 0 # the most NAs we'll accept, that is, the number of files without data for that position
#careful, na_cutoff of 1 doesn't work if we only have one replicate
homozygoteCutoff <- 0.8 #the maximum WT allele frequency we'll accept for candidates


pf <- PileupFiles(wt_list)
pf_mut <- PileupFiles(mut_list)

#Function for use in applyPileups: gets list for each pileup position
calcInfo <-
  function(x)
  {
    ## x == information for each file
    nuc <- c('A', 'C', 'G', 'T')
    #this apply goes over each column of 'seq' array, which is each file
    #that means 'y' is a matrix of nucleotide by position
    info <- apply(x[["seq"]], 2, function(y) {
      y <- y[c("A", "C", "G", "T"),,drop=FALSE]
      cvg <- colSums(y)
      #y <- y / cvg[col(y)] # normalizes reads: don't do that yet, since we'll combine files later
      #add coverage row that stores read depth
      y <- rbind(y, cvg)
      return(y)
    })
    #info is matrix of ACGTcvg by position
    list(pos=x[["pos"]], info=info) 
  }

GetDistanceDf <- function(myChrRange, peak.version = F){
  start_gdf <- proc.time()
  library(dplyr)
  library(tidyr)
  library(Rsamtools)
  tryCatch(
    {
    message('chromosome name: ', toString(seqnames(myChrRange)))
      gc()
    
    
    param <- ApplyPileupsParam(which = myChrRange, what="seq", 
                               minBaseQuality = my_minBaseQuality,
                               minMapQuality = my_minMapQuality,
                               minDepth = my_minDepth)
    
    gc()
    apply_pileup <- applyPileups(pf, FUN = calcInfo, param = param)

    #debug: option to return pileup (memory / base pairs / reps) in order to predict memory
    # return(object.size(x = apply_pileup)/length(apply_pileup[[1]]$pos)/length(pf))
    #apply_pileup_mut is calculated later, to be more conservative with memory
    
    
    #make_df_for_chromsome: makes a function that turns each info list (from pileup) into dataframe
    #So it taks a list and returns a dataframe
    make_df_for_chromosome <-
      function(my_info_list){
        x <- data.frame('pos' = rep(my_info_list$pos, each = 5), 
                   'nucleotide' = c("A", "C", "G", "T", "cvg"))
        #gets right number of columns(files) by column bind
        #because $info is info by file, and info, which was originally
        #a matrix, gets simplified into a vector which matches
        x <- cbind(x, my_info_list$info)
        message("after make_df_for_chr size is ", nrow(x))
        return(x)
    }
    
    
    ###depth filter
    #takes df with columns pos, nuc (ACGTcvg), file1, file2 ...
    #returns df with columns pos, nuc (without cvg!!), file1, file2...
    #assumption: for each file, take out positions with reads under cutoff
    #accomplished by replacing position's data with NAs for each column(file)
    #NAs will then be dealt with appropriately by NA filter
    depth_filter <- function(chrDf){
      #DEBUG: for bypassing depth_filter
      # chrDf <- chrDf[-seq(5, to=nrow(chrDf), by=5),]
      # chrDf <- droplevels(chrDf)
      # return(chrDf)

      #gets only cvg rows (every 5th), then only info cols
      #returns true for undercovered rows
      #can only subset 2 dimensions with a matrix, not a df. Don't know why.
      filter_mat <- as.matrix(chrDf[seq(5, to=nrow(chrDf), by=5), 3:ncol(chrDf)] < my_minDepth)
      
      #expand filter_mat to cover all 5 rows for each position
      new_mat <- c()
      for (i in 1:ncol(filter_mat)){
        new_col <- rep(filter_mat[,i], each=5)
        new_mat <- cbind(new_mat, new_col)
      }
      filter_mat <- new_mat
      
      #apply filter to dataframe
      chrDf[,3:ncol(chrDf)][filter_mat] <- NA
      
     
      message("after depth_filter size is ", nrow(chrDf))
      return(chrDf)      
    }
    
    
    ###Filter for NAs
    #function we can use in lapply on results list, takes dataframe, returns it after filtering for nas
    #format for df (input and output) should be pos (5 for each), nuc (ACGTcvg), file1, file2...
    na_filter <- function(chrDf){
      #if only one file and na_cutoff isn't 0
      if ( sum( !(names(chrDf) %in% c("pos", "nucleotide")) ) == 1 & na_cutoff != 0){
        message('Warning: na_cutoff for a single-file pool must be 0. Continuing with cutoff of 0.')
        na_cutoff = 0
      }
      #vector returns true on rows with cutoff or less NaNs or NAs (is.na accounts for both)
      #need drop=F so it works if only 1 column is passed to rowSums
      filter_vec <- (rowSums(is.na(chrDf[,3:length(chrDf), drop = FALSE])) <= na_cutoff)
      chrDf <- chrDf[filter_vec,]
      message("after na_filter size is ", nrow(chrDf))
      return(chrDf)
    }
    
    ###Get averages between files and divides by coverage
    #function takes dataframe with pos, nuc (ACGTcvg), file1, file2...returns it with files combined into one mean column
    #so output should be df with pos, nuc (ACGT), avg_count -- also removes coverage
    #by adding each file together then dividing by total coverage, we are essentially weighting each file by coverage
    avg_files <- function(chrDf){
      # if (length(chrDf) == 3){
      #   names(chrDf)[3] <- "avg_count"
      #   return(chrDf)
      # }
      chrDf$avg_count <- rowSums(chrDf[,3:length(chrDf), drop = FALSE], na.rm = TRUE)
      #takes non-cvg rows and divides them by cvg
      chrDf[-seq(5, nrow(chrDf), by=5), "avg_count"] <- 
        chrDf[-seq(5, nrow(chrDf), by=5), "avg_count"] / chrDf[rep(seq(5, nrow(chrDf), by=5), each = 4), "avg_count"]
      
      #remove coverage rows to make compatible with future functions
      #unless used for PeakCaller version
      if (!peak.version){
        chrDf <- chrDf[-seq(5, to=nrow(chrDf), by=5),]
        #remove cvg factor
        chrDf <- droplevels(chrDf)
      }
      
      #throw away non-mean columns and return
      message("after avg_files size is ", nrow(chrDf))
      return(chrDf[,c("pos", "nucleotide", "avg_count")])
    }
    
    #filter homozygotes wt pool, so when later inner_joined,
    #rows will be taken out of both pools
    #removes rows where major allele frequency exceeds cutoff
    #takes df with cols pos, A, C, G, T; returns same
    homozygote_filter <- function(chrDf){
      rows_to_keep <- apply((chrDf[, -which(names(chrDf) %in% c("pos", "cvg"))] <= homozygoteCutoff), MARGIN = 1, all)
      chrDf <- chrDf[rows_to_keep,]
      message("after hz_filter size is ", nrow(chrDf))
      return(chrDf)
    }
    
    
    #apply functions to wild type pool
    #CAUTION: functions must be applied in this order to work right
    tryCatch(
      wt_df <- apply_pileup[[1]] %>%
        make_df_for_chromosome() %>%
        depth_filter %>%
        na_filter() %>%
        avg_files() %>%
        spread(key = nucleotide, value = avg_count) %>%
        #homoz filter only on wt pool
        homozygote_filter(),
        
      error = function(e) {stop('empty dataframe')}
    )
    colnames(wt_df)[2:6] <- c('A.wt', 'C.wt', 'cvg.wt', 'G.wt', 'T.wt')
    rm(apply_pileup)
    
    #apply functions to mutant pool
    apply_pileup_mut <- applyPileups(pf_mut, FUN = calcInfo, param = param)
    tryCatch(
      mut_df <- apply_pileup_mut[[1]] %>%
        make_df_for_chromosome() %>%
        depth_filter %>%
        na_filter() %>%
        avg_files() %>%
        spread(key = nucleotide, value = avg_count),
      
      error = function(e) {stop('empty dataframe')}
    )
    colnames(mut_df)[2:6] <- c('A.mut', 'C.mut', 'cvg.mut', 'G.mut', 'T.mut')
    rm(apply_pileup_mut)
    
    #inner_join already removes rows without a match
    distance_df <- inner_join(wt_df, mut_df, by=c('pos'))
    if (length(distance_df) == 0) stop('empty dataframe')
    rm(wt_df)
    rm(mut_df)
    
    #calculate Euclidian distance
    distance_df$A.wt <- (distance_df$A.wt - distance_df$A.mut)^2
    distance_df$C.wt <- (distance_df$C.wt - distance_df$C.mut)^2
    distance_df$G.wt <- (distance_df$G.wt - distance_df$G.mut)^2
    distance_df$T.wt <- (distance_df$T.wt - distance_df$T.mut)^2
    #TODO change depending on peak.version
    distance_df <- transmute(distance_df, 
                             pos = pos,
                             distance = (A.wt + C.wt + G.wt + T.wt)^(1/2),
                             cvg.mut = cvg.mut,
                             A.mut = A.mut,
                             C.mut = C.mut,
                             G.mut = G.mut,
                             T.mut = T.mut
    )
    distance_df$distance <- distance_df$distance ^ my_power
    
    stopifnot(nrow(distance_df) > 0)
    print(proc.time() - start_gdf)
    return(distance_df)
  },
    
    error = function(e) {
      msg <- paste0(toString(seqnames(myChrRange))," had an error: ",e)
      return(msg)
    }
  )
}
