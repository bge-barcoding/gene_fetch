#!/bin/bash
#SBATCH --job-name=GF
#SBATCH --partition=day
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1



# Activate conda env
source path/to/profile.d/conda.sh
conda activate gene-fetch



# NCBI API parameters
EMAIL=user_email@domain.ac.uk
API_KEY= XYZ


# Gene search variables
GENE=cox1
TYPE=both
#MIN_PROT_SIZE=500    # Default
#MIN_NUC_Size=1000    # Default


# Input and output directory paths
SAMPLES_CSV=path/to/input/samples.csv
OUTPUT_DIR=path/to/out/dir



# Run gene_fetch.py
gene-fetch \
    --gene $GENE \
    --type $TYPE
    --out $OUTPUT_DIR \
    --in $SAMPLES_CSV \
	--email $EMAIL \
	--api-key $API_KEY
