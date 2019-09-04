# MMAPPR2
[![Build Status](https://travis-ci.org/kjohnsen/MMAPPR2.svg?branch=master)](https://travis-ci.org/kjohnsen/MMAPPR2)

## Mutation Mapping Analysis Pipeline for Pooled RNA-Seq
### Kyle Johnsen, Nathaniel Jenkins, Jonathon Hill

### Introduction
MMAPPR2 maps mutations resulting from pooled RNA-seq data from the F2
cross of forward genetic screens. Its predecessor is described in a paper published
in Genome Research (Hill et al. 2013). MMAPPR2 accepts aligned BAM files as well as
a reference genome as input, identifies loci of high sequence disparity between the
control and mutant RNA sequences, predicts variant effects using Ensembl's Variant
Effect Predictor, and outputs a ranked list of candidate mutations.

[See vignette for instructions](vignettes/MMAPPR2.Rmd)

Publication for the [original MMAPPR](http://genome.cshlp.org/content/23/4/687.full.pdf)

## Installation Notes
MMAPPR2 depends on two system tools to function: Samtools and VEP. Both must be installed and in the PATH to be found by the appropriate functions.

### Installing Samtools
Instructions to install samtools can be found at https://github.com/samtools/samtools and installation instructions are in the INSTALL file included with samtools.

### Installing VEP
You'll need Ensembl VEP, which you can install like this, replacing `my_species` with your species (e.g., `danio_rerio`):

    git clone https://github.com/Ensembl/ensembl-vep.git
    cd ensembl-vep
    perl INSTALL.pl -a ac -s {my_species}

This installs the most recent VEP and allows you to create a cache for your desired species, which is what MMAPPR2 expects by default.
If you depart from the installation shown here, or if things don't go smoothly, see [Ensembl's instructions](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer)
and make sure any differences are accounted for in the
[`VEPFlags`](#configure-vepflags-object) object.

*Note:* If you have any trouble installing VEP, using
[their Docker image](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#docker)
may save you a lot of hassle.

*Note:* We have found that R sometimes has issues finding VEP, especially when perlbrew is used. If you encounter errors at the path to your perl installation to the .Rprofile file. For example:

    Sys.setenv(PATH=paste("/Path/to/Perlbrew", Sys.getenv("PATH"), sep=":"))
