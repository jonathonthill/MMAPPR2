.make_synthetic_fixture <- function() {
  td <- tempfile("mmappr2_fixture_")
  dir.create(td, recursive = TRUE)

  fa <- file.path(td, "ref.fa")
  writeLines(c(">chr1", paste(rep("A", 250), collapse = "")), fa)
  Rsamtools::indexFa(fa)

  gtf <- file.path(td, "genes.gtf")
  writeLines(c(
    'chr1\tMMAPPR2\tgene\t1\t150\t.\t+\t.\tgene_id "g1";',
    'chr1\tMMAPPR2\ttranscript\t1\t150\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";',
    'chr1\tMMAPPR2\texon\t1\t40\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; exon_id "e1";',
    'chr1\tMMAPPR2\texon\t120\t150\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; exon_id "e2";',
    'chr1\tMMAPPR2\tgene\t100\t200\t.\t+\t.\tgene_id "g2";',
    'chr1\tMMAPPR2\ttranscript\t100\t200\t.\t+\t.\tgene_id "g2"; transcript_id "tx2";',
    'chr1\tMMAPPR2\texon\t100\t200\t.\t+\t.\tgene_id "g2"; transcript_id "tx2"; exon_id "e3";'
  ), gtf)

  sam <- file.path(td, "reads.sam")
  seq20 <- paste(rep("A", 20), collapse = "")
  qual20 <- paste(rep("I", 20), collapse = "")
  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:250",
    # Primary spliced RNA-seq alignment: this MUST survive pileup.
    paste("r1", 0, "chr1", 10, 60, "10M100N10M", "*", 0, 0,
          seq20, qual20, sep = "\t"),
    paste("rdup", 1024, "chr1", 50, 60, "20M", "*", 0, 0,
          seq20, qual20, sep = "\t"),
    # Secondary (0x100) and supplementary (0x800) alignments: these represent
    # alternative placements/segments and MUST NOT add a second vote to allele depth.
    paste("r2", 256, "chr1", 160, 60, "20M", "*", 0, 0,
          seq20, qual20, sep = "\t"),
    paste("r3", 2048, "chr1", 180, 60, "20M", "*", 0, 0,
          seq20, qual20, sep = "\t")
  ), sam)

  bam_prefix <- file.path(td, "reads")
  bam <- tryCatch(
    Rsamtools::asBam(sam, destination = bam_prefix, overwrite = TRUE),
    error = function(e) e
  )
  if (inherits(bam, "error")) testthat::skip(paste("SAM-to-BAM conversion unavailable:", bam$message))
  bam <- as.character(bam)[1]
  if (!file.exists(paste0(bam, ".bai")) &&
      !file.exists(sub("\\.bam$", ".bai", bam))) Rsamtools::indexBam(bam)

  list(dir = td, fasta = fa, gtf = gtf, bam = bam)
}
