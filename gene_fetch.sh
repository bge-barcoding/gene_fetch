#!/bin/bash
#SBATCH --job-name=GF
#SBATCH --partition=medium
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1



# Activate conda env
#source path/to/profile.d/conda.sh
#conda activate gene-fetch



# NCBI API parameters
EMAIL=example@example.ac.uk

API_KEY="XYZ"


# Gene search variables
GENE=cox1
TYPE=both
#MIN_PROT_SIZE=500    # Default
#MIN_NUC_Size=1000    # Default


# Input and output directory paths
SAMPLES_CSV=./data/samples.csv
#TAXONOMY_CSV=./data/samples_taxonomy.csv

OUTPUT_DIR=./test_samples



# Run gene_fetch.py
gene-fetch \
    --gene $GENE \
    --type $TYPE \
    --out $OUTPUT_DIR \
    --in $SAMPLES_CSV \
    --email $EMAIL \
    --api-key $API_KEY \
    -gb
