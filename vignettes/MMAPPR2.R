## ----installVEP, eval=FALSE------------------------------------------------
#  git clone https://github.com/Ensembl/ensembl-vep.git
#  cd ensembl-vep
#  perl INSTALL.pl -a ac -s {my_species}

## ----param-----------------------------------------------------------------
BiocParallel::register(BiocParallel::MulticoreParam())  ## see below for explanation of BiocParallel
library(MMAPPR2, quietly = TRUE)
library(MMAPPR2data, quietly = TRUE)
library(Rsamtools, quietly = TRUE)

# This is normally configured automatically:
vepFlags <- ensemblVEP::VEPFlags(flags = list(
    format = 'vcf',  # <-- this is necessary
    vcf = FALSE,  # <-- as well as this
    species = 'danio_rerio',
    database = FALSE,  # <-- these three arguments allow us to run VEP offline,
    fasta = goldenFasta(),  # <-╯|  which you probably won't need
    gff = goldenGFF(),  # <------╯
    filter_common = TRUE,
    coding_only = TRUE  # assuming RNA-seq data
))

param <- MmapprParam(refFasta = goldenFasta(),
                     wtFiles = exampleWTbam(),
                     mutFiles = exampleMutBam(),
                     species = 'danio_rerio',
                     vepFlags = vepFlags,  ## optional
                     outputFolder = tempOutputFolder())  ## optional

## ----mmappr----------------------------------------------------------------
mmapprData <- mmappr(param)

