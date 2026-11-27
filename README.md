# MMAPPR2
[![Build Status](https://travis-ci.org/kjohnsen/MMAPPR2.svg?branch=master)](https://travis-ci.org/kjohnsen/MMAPPR2)

## Mutation Mapping Analysis Pipeline for Pooled RNA-Seq
### Kyle Johnsen, Nathaniel Jenkins, Jonathon Hill

### Introduction
MMAPPR2 maps mutations resulting from pooled RNA-seq data from the F2
cross of forward genetic screens. Its predecessor is described in a paper published
in Genome Research (Hill et al. 2013). MMAPPR2 accepts aligned BAM files as well as
a reference genome as input, identifies loci of high sequence disparity between the
control and mutant RNA sequences, predicts variant effects, 
and outputs a ranked list of candidate mutations.

[See vignette for instructions](vignettes/MMAPPR2.Rmd)

Publication for the [original MMAPPR](http://genome.cshlp.org/content/23/4/687.full.pdf)

## Installation Notes
MMAPPR2 depends on Samtools to function. It must be installed and in the PATH to be found by the appropriate functions.

### Installing Samtools
Instructions to install samtools can be found at https://github.com/samtools/samtools and installation instructions are in the INSTALL file included with samtools.

# Installation guide via UBUNTU OS
### type the following code directly into ubuntu operating system after installation to begin using MMAPPR2

## installing miniforge for environment
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

## Make it executable
chmod +x Miniforge3-Linux-x86_64.sh

## Run it
./Miniforge3-Linux-x86_64.sh

## Create environment
mamba create -n mmappr2 -c conda-forge -c bioconda \
  r-base=4.3 \
  r-biocmanager \
  r-devtools \
  r-curl r-xml2 r-openssl r-v8 \
  r-systemfonts r-textshaping r-ragg r-svglite r-magick \
  r-gh r-gert r-usethis r-pkgdown r-rcmdcheck r-roxygen2 \
  bioconductor-biocparallel \
  bioconductor-rsamtools \
  bioconductor-variantannotation \
  bioconductor-biostrings \
  bioconductor-xvector \
  samtools \
  git \
  -y


## Activate environment
mamba activate
(base) username@____:~$ mamba activate mmappr2

## Update apt dependencies
sudo apt-get update && sudo apt-get install -y \
build-essential \
libcurl4-openssl-dev \
libssl-dev \
libxml2-dev \
libv8-dev \
libfreetype6-dev \
libpng-dev \
libtiff-dev \
libjpeg-dev \
libharfbuzz-dev \
libfribidi-dev \
libfontconfig1-dev \
libcairo2-dev \
libmagick++-dev \
libwebp-dev \
samtools

## Activate R
R --vanilla
## OR
R --no-restore --no-save

### inside of R:
install.packages(c(
  "urlchecker",
  "rversions"
))

### install Rsamtools
BiocManager::install("Rsamtools", force=TRUE)

******** Select Update NONE ********

## Pull MMAPPR2 and install from Github
devtools::install_github("jonathonthill/MMAPPR2")

******** UPDATE NONE (SELECT 3) ********

## Run in R again:
library(MMAPPR2)

## Setup your parameters based on your machine path, see example ->
mmappr_param <- mmapprParam(wtFiles="/mnt/c/Users/Dell/Desktop/MMAPPR2_Files/Tbx1_aligned_wt.sorted.bam", 
mutFiles="/mnt/c/Users/Dell/Desktop/MMAPPR2_Files/Tbx1_aligned_mut.sorted.bam", refFasta="/mnt/c/Users/Dell/Desktop/MMAPPR2_Files/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa", gtf="/mnt/c/Users/Dell/Desktop/MMAPPR2_Files/Danio_rerio.GRCz11.100.gtf",
outputFolder="/mnt/c/Users/Dell/Desktop/MMAPPR2_Files/MMAPPR2_out")

mmapprdata <- mmappr(mmappr_param)

## To deactivate/remove/check installed system environments:
mamba deactivate
mamba env remove -n MMAPPR2
mamba env list

