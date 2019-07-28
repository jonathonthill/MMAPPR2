.samtoolsPileup <- function(files, param, chrRange) {
    # Build command
    args <- paste("mpileup -ERI",   #Redo Baq, ignore readgroups, and skip indels
                 "-f", refFasta(param),
                 "-C 50",
                 "--min-MQ", minMapQuality(param),
                 "--min-BQ", minBaseQuality(param),
                 "--region", as.character(chrRange, ignore.strand=TRUE))
    # Run command and read results for each file
    for (file in BiocGenerics::path(files)) {
        pileupDT <- data.table::fread(text=system2(command="samtools",
                                        args=paste(args, file),
                                        stdout=TRUE,
                                        stderr=FALSE))
        
        # Clean extra symbols for beginning and end of reads and sub in letters
        pileupDT$V5 <- toupper(pileupDT$V5)
        pileupDT$V5 <- stringr::str_remove_all(pileupDT$V5, "\\^.|\\$")
        pileupDT$V5 <- stringr::str_replace_all(pileupDT$V5, "\\.|,", pileupDT$V3)
        
        # Filter for min depth
        pileupDT <- pileupDT[pileupDT$V4 >= minDepth(param), ]
        
        # Calculate allele freqs
        pos <- pileupDT$V2
        counts <- rbind(stringr::str_count(pileupDT$V5, "A"),
                       stringr::str_count(pileupDT$V5, "C"),
                       stringr::str_count(pileupDT$V5, "G"),
                       stringr::str_count(pileupDT$V5, "T"),
                       cvg = pileupDT$V4)
        counts <- as.vector(counts)
        stopifnot(length(counts) > 0)
        countsDf <- data.frame(nucleotide=c("A", "C", "G", "T", "cvg"), pos=rep(pos, each=5), counts=counts)
        
        # Assemble
        if (exists("outdat")) {
            outdat <- merge(outdat, countsDf, by=c("pos", "nucleotide"), all=TRUE, sort=FALSE)
        } else {
            outdat <- countsDf
        }
    }
    stopifnot(is.data.frame(outdat))
    colnames(outdat) <- c(colnames(outdat)[1:2], names(files))
    return(outdat)
}
