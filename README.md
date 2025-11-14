# PARP Inhibitor Genetic Interaction Manuscript

This repository contains the stable, reproducible analysis code and data for the PARP inhibitor genetic interaction manuscript. It is a specifically configured to regenerate all analyses and key figures presented in the publication.

**📊 Interactive Data Portal**: https://parpi.princeton.edu/map

**🔧 General-Purpose Pipeline**: For analyzing your own genetic interaction screens, please consider [GI Nexus](https://github.com/simpsondl/gi-nexus), an actively maintained version of the analysis pipeline.

## Purpose

This repository is designed for:
- **Reproducing manuscript figures and analyses** using the exact configurations and data from the publication
- **Understanding the methods** used to process the PARP inhibitor screens
- **Accessing the raw screen data** and intermediate analysis outputs

For general genetic interaction analysis of new datasets, check out [GI Nexus](https://github.com/simpsondl/gi-nexus).

## Repository Contents

- `manuscript_data/` - Raw count tables and annotations for the two PARP inhibitor screens (2022 and 2023)
- `config/config.yaml` - Exact configuration parameters used for the manuscript analyses
- `workflow/` - Snakemake workflow and R scripts (stable version from publication)
- `workflow/scripts/manuscript_figures/` - Figure generation scripts

## Quick Start

### 1. Clone this repository
```bash
git clone https://github.com/simpsondl/parpi-manuscript.git
cd parpi-manuscript
```

### 2. Set up the environment
```bash
# Create and activate conda environment
conda env create -f workflow/envs/smk-env.yaml
conda activate differential_gi_smk
```

### 3. Run the complete analysis
```bash
# Dry run to preview execution plan
snakemake --use-conda -n --snakefile workflow/Snakefile --configfile config/config.yaml

# Execute the full pipeline to generate all outputs and figures
snakemake --use-conda --cores 4 --snakefile workflow/Snakefile --configfile config/config.yaml
```

The pipeline will process both screens and generate all manuscript figures. Expected runtime: ~1 hour on a desktop (16GB RAM, 6 cores).

## Requirements

- **Snakemake** (≥9.12.0)
- **R** (≥4.0) with tidyverse packages
- **Conda** (recommended for reproducible environments)

All dependencies are specified in `workflow/envs/`.

## Outputs

Results will be generated in the `outputs/` directory:
- `phenotypes/` - Calculated phenotypes and filtered interaction scores
- `gi_scores/` - Genetic interaction scores at construct and gene levels
- `figures/` - Manuscript figures (PDF format)
- `misc_results/` - Quality control metrics and supplementary data

## Citation

If you use these datasets or workflow in your work, please cite the following reference:

> Simpson D, Ling J, Jing Y, Adamson B. Mapping the Genetic Interaction Network of PARP inhibitor Response. bioRxiv [Preprint]. 2023 Aug 20:2023.08.19.553986. doi: 10.1101/2023.08.19.553986. PMID: 37645833; PMCID: PMC10462155.