# Mapping the genetic interaction network of PARPi response

This repository contains the reproducible analysis code and data for our PARP inhibitor genetic interaction manuscript. It is a frozen snapshot of an early GI Nexus build specifically configured to regenerate all analyses and key figures presented in the publication.

**📊 Interactive Data Portal**: https://parpi.princeton.edu/map

## Purpose

This repository is designed for:
- **Reproducing manuscript figures and analyses** using the exact configurations and data from the publication
- **Understanding the methods** used to process the PARP inhibitor screens
- **Accessing the raw screen data** and intermediate analysis outputs

**🔧 General-Purpose Pipeline**: For analyzing your own genetic interaction screens, please consider [GI Nexus](https://github.com/simpsondl/gi-nexus), a more-developed and feature-rich version of the analysis pipeline implemented here and in our manuscript.

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
# Create and activate a minimal conda environment
# Solely establishes snakemake installation
conda env create -f workflow/envs/min-env.yaml
conda activate minmanuscript
```

### 3. Run the complete analysis
```bash
# Dry run to preview execution plan
snakemake -n --snakefile workflow/Snakefile process_screens
# If you installed the minimal environment, add --use-conda flag  
# snakemake --use-conda -n --snakefile workflow/Snakefile process_screens

# Execute the full pipeline to generate all outputs and figures
snakemake --cores 6 --snakefile workflow/Snakefile process_screens
# If you installed the minimal environment, add --use-conda flag  
# snakemake --use-conda --cores 6 --snakefile workflow/Snakefile process_screens
```

The pipeline will process both screens and generate all data required to make the manuscript figures. Expected runtime: ~1 hour on a desktop (16GB RAM, 6 cores).

### 4. Generate the manuscript figures
```bash
# Dry run to preview execution plan
snakemake -n --snakefile workflow/Snakefile generate_manuscript_figures
# If you installed the minimal environment, add --use-conda flag  
# snakemake --use-conda -n --snakefile workflow/Snakefile process_screens

# Execute the full pipeline to generate all outputs and figures
snakemake --cores 6 --snakefile workflow/Snakefile process_screens
# If you installed the minimal environment, add --use-conda flag  
# snakemake --use-conda --cores 6 --snakefile workflow/Snakefile process_screens
```

## Requirements

- **Snakemake** (≥9.12.0)
- **R** (≥4.0) with tidyverse packages
- **Conda** (recommended for reproducible environments)

All dependencies are specified in `workflow/envs/`.

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