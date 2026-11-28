# MMAPPR2

*Mutation Mapping Analysis Pipeline for Pooled RNA-Seq*

Kyle Johnsen, Nathaniel Jenkins, and Jonathon Hill

## Introduction

MMAPPR2 identifies candidate mutations using pooled RNA-seq data from F2 crosses
in forward genetic screens. The pipeline accepts aligned BAM files, a reference
genome, and matching gene annotation as input. It identifies loci of high
sequence disparity between wild-type (WT) and mutant pools, predicts variant
effects, and produces a ranked list of candidate single-nucleotide variants
(SNVs).

The original MMAPPR method was described in
[Hill et al. (2013), *Genome Research*](https://genome.cshlp.org/content/23/4/687.full.pdf).
See the [full workflow guide](vignettes/MMAPPR2.Rmd) for additional details and
analysis options.

## Installation

MMAPPR2 requires **R 4.4.0 or newer**. The installation instructions below
handle the required CRAN and Bioconductor packages. MMAPPR2 no longer requires
the standalone Samtools executable.

MMAPPR2 was [removed from Bioconductor](https://bioconductor.org/about/removed-packages/)
in release 3.20, so `BiocManager::install("MMAPPR2")` will not install the current
version. Bioconductor is still used to obtain its required dependencies.

### Standard R installation

#### Install from GitHub

In R, run:

```r
install.packages(c("remotes", "BiocManager"))

remotes::install_github(
  "jonathonthill/MMAPPR2",
  dependencies = NA,
  upgrade = "never",
  repos = BiocManager::repositories()
)
```

This installs MMAPPR2 from GitHub and resolves its required R dependencies
using the CRAN and Bioconductor repositories appropriate for your R version.
Optional development and documentation packages are not installed, and existing
packages are not automatically upgraded. A fresh R library or isolated environment
is recommended if you encounter conflicting package versions.

#### Install from downloaded source

Download and extract the desired version of the MMAPPR2 source package. Set R's
working directory to the extracted directory containing `DESCRIPTION`, then run:

```r
install.packages(c("remotes", "BiocManager"))

remotes::install_local(
  ".",
  dependencies = NA,
  upgrade = "never",
  repos = BiocManager::repositories()
)
```

This method is useful for installing a specific source version or working with
locally modified code.

#### Verify the installation

Confirm that the package is installed and can be loaded:

```r
library(MMAPPR2)
packageVersion("MMAPPR2")
```

### Installation with Miniforge and Mamba

An isolated Conda environment can simplify installation on Linux by providing
precompiled R packages and their system libraries. Ubuntu and Miniforge are
not required if the standard R installation works on your system.

The following example uses R 4.5. The complete Mamba installation route has not
yet been tested. Run the shell commands in your Linux terminal, not in the R
console.

#### 1. Install Miniforge

If you do not already have Conda or Mamba, download and run the Miniforge
installer. For a 64-bit x86 Linux system:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

Follow the installer prompts, allow shell initialization, and reopen your terminal.
For other architectures, see the [Miniforge project](https://github.com/conda-forge/miniforge).

#### 2. Create and activate the environment

```bash
mamba create -n mmappr2 \
  -c conda-forge -c bioconda \
  --strict-channel-priority \
  r-base=4.5 \
  r-remotes \
  r-biocmanager \
  r-data.table \
  r-dbi \
  bioconductor-rsamtools \
  bioconductor-variantannotation \
  bioconductor-biocparallel \
  bioconductor-biocgenerics \
  bioconductor-genomeinfodb \
  bioconductor-genomicranges \
  bioconductor-genomicfeatures \
  bioconductor-txdbmaker \
  bioconductor-genomicalignments \
  bioconductor-iranges \
  bioconductor-s4vectors \
  bioconductor-summarizedexperiment \
  bioconductor-biostrings \
  bioconductor-rtracklayer \
  -y

mamba activate mmappr2
R
```

If using Conda instead of Mamba, replace `mamba` with `conda` in these commands.

#### 3. Install MMAPPR2 inside R

```r
remotes::install_github(
  "jonathonthill/MMAPPR2",
  dependencies = FALSE,
  upgrade = "never"
)

library(MMAPPR2)
packageVersion("MMAPPR2")
```

All required R dependencies are specified in the environment creation command,
so `dependencies = FALSE` avoids replacing Conda-managed packages with packages
built separately from source. To install a downloaded version instead, set the
working directory to its extracted source directory and use
`remotes::install_local(".", dependencies = FALSE, upgrade = "never")`.

#### Manage the environment (optional)

After exiting R with `q()`, run the following commands in your terminal as needed:

```bash
# Leave the environment
mamba deactivate

# List existing environments
mamba env list

# Return to the environment for another analysis
mamba activate mmappr2
R
```

If the environment is no longer needed, exit R, deactivate it, and run
`mamba env remove -n mmappr2`. This removes the environment and its installed
packages; keep your input data and analysis results outside the environment.

### Troubleshooting installation

MMAPPR2 itself does not require compilation, but some of its R dependencies do.
If you are installing directly into a native Ubuntu R installation and packages
fail to compile from source, the following system libraries may help:

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  zlib1g-dev libbz2-dev liblzma-dev
```

These are a starting point for source-build errors, not an exhaustive list of
system dependencies. Consult the failing package's installation log for any
additional requirements. APT packages are not normally needed for the Mamba
installation above.

For a standard R installation, `BiocManager::valid()` can help diagnose
incompatible Bioconductor package versions. Prefer a consistent R/Bioconductor
environment instead of forcing individual packages to reinstall.

## Prepare input data

MMAPPR2 requires the following files:

| Input | Requirements |
| --- | --- |
| Wild-type BAM files | Coordinate-sorted RNA-seq alignments from the wild-type pool. |
| Mutant BAM files | Coordinate-sorted RNA-seq alignments from the mutant pool. |
| Reference FASTA | The same reference genome assembly used to align the reads. |
| Gene annotation | A matching GTF or GFF3 file for coding-consequence prediction and gene-expression summaries. |

One or more BAM files can be provided for each pool. All BAMs must already be
coordinate-sorted; MMAPPR2 does not sort or align reads. Missing, empty, or stale
BAM and FASTA indexes are created or rebuilt automatically, so the corresponding
input directories must be writable when indexing is needed.

BAM and FASTA sequence names and lengths must agree for the BAM sequences (for
example, `chr1` and `1` are not interchangeable). Use annotations corresponding
to the same reference assembly and chromosome naming scheme.

## Run MMAPPR2

Replace the example filenames with the paths to your own files. Relative paths
are resolved from the current R working directory; absolute paths can also be
used.

```r
library(MMAPPR2)

param <- mmapprParam(
  wtFiles = "wt_pool.bam",
  mutFiles = "mutant_pool.bam",
  refFasta = "reference.fa",
  gtf = "annotation.gtf",
  outputFolder = "mmappr2_results"
)

results <- mmappr(param)
```

Use a new or empty output directory to avoid accidentally replacing earlier
results. To supply multiple BAM files per pool, use a character vector such as
`wtFiles = c("wt1.bam", "wt2.bam")` and
`mutFiles = c("mut1.bam", "mut2.bam")`.

See `?mmapprParam` and the [full workflow guide](vignettes/MMAPPR2.Rmd) for
additional analysis options.

## Output files and checking results

Results are written to `mmappr2_results` (or your chosen output directory).
Some files are produced only when sufficient data or linkage evidence are
available.

| File | Description |
| --- | --- |
| `AllCandidateVariantsFor<seq>.tsv` | Candidate SNVs, with evidence from both the WT and mutant pools. |
| `DetectedMutationsFor<seq>.tsv` | Predicted effects of candidate variants in protein-coding genes. |
| `DifferentiallyExpressedGenesFor<seq>.tsv` | Gene-expression differences; **not** a statistical differential-expression test. |
| `genome_plots.pdf` and `peak_plots.pdf` | Figures showing regions with the strongest evidence of linkage. |
| `aicc_plots.pdf` | Optional figures showing how the smoothing setting was chosen. |
| `mmappr_data.RDS` | Serialized analysis data for later inspection and regeneration of figures. |
| `mmappr2.log` | Analysis parameters, progress, messages, errors, and R session information. |

`<seq>` is the chromosome or other reference sequence name. TSV files are
tab-separated tables that can be opened in a spreadsheet or imported into R.
Begin by reviewing the genome and peak plots, then examine the candidate tables.
To save `aicc_plots.pdf`, add `exportAiccPlots = TRUE` to `mmapprParam()`.

The pipeline can return and save partial analysis results if a stage fails.
Check `mmappr2.log` for errors before interpreting the results; the presence of
`mmappr_data.RDS` alone does not establish that the analysis completed.

## Figure colors and replotting

### Change colors for a new analysis

To select a color palette for a new analysis, replace `mmappr(param)` with:

```r
results <- mmappr(param, palette = "colorblind")
```

The available palettes are `default`, `colorblind`, and `grayscale`.

### Replot saved results

You can regenerate PDF figures from saved analysis results **without rerunning
MMAPPR2**:

```r
plotMmapprData(
  "mmappr2_results/mmappr_data.RDS",
  plots = c("genome", "peaks"),
  palette = "grayscale"
)
```

New figures are written to a `replots` subdirectory by default. The saved RDS
and candidate tables remain unchanged. See `?plotMmapprData` for individual
color overrides and other options.

**Be careful with `overwrite = TRUE`:** in `plotMmapprData()` it permits replacing
only the selected PDF figures; in `mmapprParam()` it permits clearing an existing
results directory.

## Limitations

MMAPPR2 detects candidate single-nucleotide variants, not insertions, deletions,
or larger changes. Uneven RNA-seq coverage can cause variants to be missed;
gene expression and read alignment can also affect results. The analysis reports
one main linked region per chromosome. The gene-expression summary is
not a statistical differential-expression test.

## Documentation, citation, and support

Read the [full workflow guide](vignettes/MMAPPR2.Rmd), or use `?mmapprParam`,
`?mmappr`, and `?plotMmapprData` in R for detailed help.

The predecessor of MMAPPR2 was described in the
[original MMAPPR paper (Hill et al., 2013)](https://genome.cshlp.org/content/23/4/687.full.pdf).
To cite this package, run `citation("MMAPPR2")` in R.

Authors: Kyle Johnsen, Nathaniel Jenkins, and Jonathon Hill.

License: [GPL-3](https://www.gnu.org/licenses/gpl-3.0.html).

Report problems through [GitHub issues](https://github.com/jonathonthill/MMAPPR2/issues).
