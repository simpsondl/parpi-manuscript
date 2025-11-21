# Mapping the Genetic Interaction Network of PARPi Response

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.12.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io) [![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

This repository contains the reproducible analysis code and data for our PARP inhibitor genetic interaction manuscript. It is a frozen snapshot of an early GI Nexus build (v1.0-beta) specifically configured to regenerate all analyses and key figures presented [in our publication](https://www.biorxiv.org/content/10.1101/2023.08.19.553986v1).

**📊 Interactive Data Portal**: https://parpi.princeton.edu/map

## Purpose

This repository is designed for:
- **Reproducing manuscript figures and analyses** using the exact configurations and data from [our publication](https://www.biorxiv.org/content/10.1101/2023.08.19.553986v1)
- **Understanding the methods** used to process our screens
- **Accessing the raw screen data** and intermediate analysis outputs

**🔧 General-Purpose Pipeline**: For analyzing your own genetic interaction screens, please consider [GI Nexus](https://github.com/simpsondl/gi-nexus), a more-developed and feature-rich version of the analysis pipeline implemented here.

## Repository Content Highlights

- `manuscript_data/` - Raw count tables and annotations for the two PARP inhibitor screens (2022 and 2023)
- `tutorial/` - R markdown document walking through the processing steps in the pipeline, function by function
- `config/config.yaml` - Exact configuration parameters used for the manuscript analyses
- `workflow/` - Snakemake workflow, conda environments, and R scripts oh my!
- `workflow/scripts/manuscript_figures/` - Figure generation scripts, for those wondering about the magic

## Quick Start

### 1. Clone this repository
```bash
git clone https://github.com/simpsondl/parpi-manuscript.git
cd parpi-manuscript/workflow
```

### 2. Set up the environment
```bash
# Create and activate a minimal conda environment
# Establishes snakemake/python installation only
conda env create -f workflow/envs/min-env.yaml
conda activate minmanuscript
```

### 3. Run the complete analysis
```bash
# Dry run to preview execution plan
snakemake --use-conda -n process_screens

# Execute the full pipeline to generate all outputs and figures
snakemake --use-conda --cores 6 process_screens
```

The pipeline will process both screens and generate all data required to make the manuscript figures. Expected runtime: ~1 hour on a desktop (16GB RAM, 6 cores). This part of the pipeline does not use any bioconductor packages, and should work on all operating systems including Windows.

### 4. Generate the manuscript figures
```bash
# Dry run to preview execution plan
snakemake --use-conda -n generate_manuscript_figures

# Execute the rule to generate key figures
snakemake --use-conda --cores 6 generate_manuscript_figures
```

This part of the pipeline uses bioconductor packages, which can not be accessed through conda on Windows. It is still possible to use the pipeline on Windows by establishing your own environment first and installing needed packages manually, then adjusting `manuscript_figures.smk` appropriately. All used packages are available on Windows.

Please note that the figure scripts output raw versions of figures which were then paneled in Adobe Illustrator. Label and text sizes, text placement, and legend placement were frequently changed in post-processing to accommodate figure structure.

## Requirements

- **Conda** (recommended for reproducible environments)

That's it. Conda/mamba handles all environments and packages. Dependencies are specified in environment files in `workflow/envs/`.

Other useful things:
- **Disk space** you will need ~20GB free space to generate all files, but ultimately only ~11Gb are used after temporary files are deleted.
- **RAM and cores** a reasonable amount (16Gb, 3-6 cores) is sufficient for relatively quick processing
- **snakemake** you can skip setting up the `minmanuscript` environment if this is available (version $\geq$ 9.12.0)
- **R** useful if you are on Windows and are unable to use conda for figure generation
- **WSL** Windows Subsystem for Linux, also useful if you are on Windows

## Outputs

Results and logs will be generated in the `outputs/` directory:
- `logs/` - Detailed logs for each step
- `misc_results/` - Quality control metrics and supplementary data
- `phenotypes/` - Calculated phenotypes and filtered interaction scores
- `gi_scores/` - Genetic interaction scores at construct and gene levels
- `manuscript_figures/` - Manuscript figures

## Citation

If you use these datasets or workflow in your work, please cite the following reference:

> Simpson D, Ling J, Jing Y, Adamson B. Mapping the Genetic Interaction Network of PARP inhibitor Response. bioRxiv [Preprint]. 2023 Aug 20:2023.08.19.553986. doi: 10.1101/2023.08.19.553986. PMID: 37645833; PMCID: PMC10462155.