#!/bin/bash
#SBATCH --job-name=GFU_1
#SBATCH --partition=day
#SBATCH --output=%x.out
#SBATCH --mem=75G
#SBATCH --cpus-per-task=12
#SBATCH --error=%x.err


# Path to conda and activation of environment
source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate fetch


# Initial directory unlock
snakemake --snakefile ./workflow/Snakefile2 --configfile ./config/config.yaml --unlock

# Run snakemake workflow
snakemake \
   --snakefile ./workflow/Snakefile \
   --configfile ./config/config.yaml \
   --cores 12

