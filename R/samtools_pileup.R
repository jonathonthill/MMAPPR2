.samtoolsPileup <- function(files, param, chrRange) {
  # Build command
  args <- paste("mpileup -EI6",   #Redo Baq, ignore readgroups, and skip indels
                "-f", refFasta(param),
                "-a AD", 
                "--min-MQ", minMapQuality(param),
                "--min-BQ", minBaseQuality(param),
                "--region", as.character(chrRange, ignore.strand=TRUE))
  
  # Run command and read results for each file
  for (file in BiocGenerics::path(files)) {
    pileupDT <- fread(text=system2(command="bcftools",
                                   args=paste(args, file, "| bcftools call -Am"),
                                   stdout=TRUE,
                                   stderr=FALSE),
                      skip = "#CHROM",
                      select = c(1,2,4,5,8))
    
    setnames(pileupDT, "#CHROM", "CHROM")
    ref1m <- "\\d+(?=(,\\d+){3};)"                                               # 4th from last number
    ref2m <- "\\d+(?=(,\\d+){2};)"                                               # 3rd from last number
    alt1m <- "\\d+(?=(,\\d+){1};)"                                               # 2nd from last number
    alt2m <- "\\d+(?=;MQ=)"                                                      # last number
    pileupDT <- pileupDT[, REF1 := as.numeric(regmatches(INFO, regexpr(ref1m, INFO, perl = TRUE)))
                         ][, REF2 := as.numeric(regmatches(INFO, regexpr(ref2m, INFO, perl = TRUE)))
                         ][, ALT1 := as.numeric(regmatches(INFO, regexpr(alt1m, INFO, perl = TRUE)))
                         ][, ALT2 := as.numeric(regmatches(INFO, regexpr(alt2m, INFO, perl = TRUE)))
                         ][, CVG := REF1 + REF2 + ALT1 + ALT2                    # Calculate coverage depth
                         ][CVG > minDepth(param)                                 # Filter for minDepth
                         ][ALT == ".", ALT := REF                                # Replace placeholder <*> with ALT allele
                         ][ALT %in% c("A", "C", "G", "T")                        # Select rows with unambiguous calls
                         ][, REF_FREQ := (REF1 + REF2)/CVG                       # Calculate relative freqs
                         ][, ALT_FREQ := (ALT1 + ALT2)/CVG
                         ][, .(CHROM, POS, REF, ALT, CVG, REF_FREQ, ALT_FREQ)]   # Select columns 
    
    # Assemble Multiple Files
    if (exists("outdat")) {
      outdat <- rbindlist(list(outdat, pileupDT))
    } else {
      outdat <- pileupDT
    }
  }
  stopifnot(is.data.frame(outdat))
  
  return(outdat)
}
