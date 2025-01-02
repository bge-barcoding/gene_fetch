#!/bin/bash
#SBATCH --job-name=GoFetch-mito
#SBATCH --partition=day
#SBATCH --output=%x.out
#SBATCH --mem=75G
#SBATCH --cpus-per-task=12
#SBATCH --error=%x.err



source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate fetch



snakemake --snakefile ./workflow/Snakefile2 --configfile ./config/config.yaml --unlock

snakemake \
   --snakefile ./workflow/Snakefile2 \
   --configfile ./config/config.yaml \
   --cores 12 \
   --rerun-incomplete \
   --keep-going

