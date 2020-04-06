## ---- echo = FALSE------------------------------------------------------------
dataDir <- system.file('extdata', package = 'MMAPPR2data')
WTpileupFile <- file.path(dataDir, 'exwt.plp')
MTpileupFile <- file.path(dataDir, 'exmut.plp')
samtoolsScript <- file('/tmp/samtools', "a")
writeLines(c(
  'if [[ ${@:$#} == *"wt.bam"* ]];',
  'then',
  paste('cat', WTpileupFile),
  'else',
  paste('cat', MTpileupFile),
  'fi'
  ), samtoolsScript)
close(samtoolsScript)

origPath <- Sys.getenv('PATH')
Sys.setenv(PATH = paste(origPath, '/tmp', sep = ':')) 

## ----installVEP, eval=FALSE---------------------------------------------------
#  git clone https://github.com/Ensembl/ensembl-vep.git
#  cd ensembl-vep
#  perl INSTALL.pl -a ac -s {my_species}

## ---- eval = FALSE------------------------------------------------------------
#  Sys.setenv(PATH=paste("/Path/to/Perlbrew", Sys.getenv("PATH"), sep=":"))

## ----param--------------------------------------------------------------------
BiocParallel::register(BiocParallel::MulticoreParam())  ## see below for explanation of BiocParallel
library(MMAPPR2, quietly = TRUE)
library(MMAPPR2data, quietly = TRUE)

param <- MmapprParam(refFasta = goldenFasta(),
                     wtFiles = exampleWTbam(),
                     mutFiles = exampleMutBam(),
                     species = 'danio_rerio',
                     outputFolder = tempOutputFolder(),
                     minDepth = 10)  ## optional

## ----mmappr-------------------------------------------------------------------
mmapprData <- mmappr(param)

