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
EMAIL="d.parsons@nhm.ac.uk"
API="1866f9734a06f26bc5895a84387542ac9308"


# Gene search variables
GENE=cbx1
TYPE=nucleotide
#MIN_PROT_SIZE=500    # Default
#MIN_NUC_Size=1000    # Default


# Input and output directory paths
SAMPLES_CSV="test2.csv"
#./data/samples.csv
#TAXONOMY_CSV=./data/samples_taxonomy.csv

OUTPUT_DIR=./test_samples



# Run gene_fetch.py
#gene-fetch \
#    --gene $GENE \
#    --type $TYPE \
#    --out $OUTPUT_DIR \
#    --email $EMAIL \
#    --api-key $API \
#    --in $SAMPLES_CSV
	
gene-fetch -e $EMAIL -k $API \
    -g rbcL -o ./output/out3 -s 3702 \
    --type protein --protein-size 400 --max-sequences 1000

