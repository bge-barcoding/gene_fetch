#!/bin/bash
#SBATCH --job-name=GF
#SBATCH --partition=day
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1



# Activate conda env
conda init bash
conda activate gene-fetch

# NCBI API parameters
EMAIL=user_email@domain.ac.uk
API_KEY= XYZ

# Gene search variables
GENE=cox1
TYPE=both
#MIN_PROT_SIZE=500    # Default
#MIN_NUC_Size=1000    # Default
#MAX_SEQS=1000

# Input and output directory paths
SAMPLES_CSV=path/to/input/samples.csv
#TAXONOMY_CSV=path/ti/input/samples_taxonomy.csv

OUTPUT_DIR=path/to/out/dir


# Run gene_fetch.py
gene-fetch \
    $GENE \
    $OUTPUT_DIR \
    $SAMPLES_CSV \
    --gene $GENE \
    --type $TYPE
