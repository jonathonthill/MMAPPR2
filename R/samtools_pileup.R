.samtoolsPileup <- function(files, param, chrRange) {
    # Build command
    args <- paste("mpileup -ERI",   #Redo Baq, ignore readgroups, and skip indels
                 "-f", param@fasta,
                 "-C 50",
                 "--min-MQ", param@minMapQuality,
                 "--min-BQ", param@minBaseQuality,
                 "--region", as.character(chrRange, ignore.strand = TRUE))
    # Run command and read results for each file
    for (file in BiocGenerics::path(files)) {
        pile_dat <- data.table::fread(text = system2(command = "samtools",
                                        args = paste(args, file),
                                        stdout = TRUE))
        
        # Clean extra symbols for beginning and end of reads and sub in letters
        pile_dat$V5 <- toupper(pile_dat$V5)
        pile_dat$V5 <- stringr::str_remove_all(pile_dat$V5, "\\^.|\\$")
        pile_dat$V5 <- stringr::str_replace_all(pile_dat$V5, "\\.|,", pile_dat$V3)
        
        # Filter for min depth
        pile_dat <- pile_dat[pile_dat$V4 >= minDepth(param), ]
        
        # Calculate allele freqs
        pos <- pile_dat$V2
        counts <- rbind(stringr::str_count(pile_dat$V5, "A"),
                       stringr::str_count(pile_dat$V5, "C"),
                       stringr::str_count(pile_dat$V5, "G"),
                       stringr::str_count(pile_dat$V5, "T"),
                       cvg = pile_dat$V4)
        counts <- as.vector(counts)
        countsdf <- data.frame(nucleotide=c("A", "C", "G", "T", "cvg"), pos = rep(pos, each = 5), counts = counts)
        
        # Assemble
        if (exists("outdat")) {
            outdat <- merge(outdat, countsdf, by = c("pos", "nucleotide"), all = TRUE, sort = FALSE)
        } else {
            outdat <- countsdf
        }
    }
    return(outdat)
}
