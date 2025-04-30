#!/bin/bash
#SBATCH --job-name=GF
#SBATCH --partition=day
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8



# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate fetch

# NCBI API parameters
EMAIL=user_email@domain.ac.uk
API_KEY= XYZ

# Gene search variables
GENE=cox1
TYPE=both

# Input and output directory paths
SAMPLES_CSV=path/to/input/samples.csv
#TAXONOMY_CSV=path/ti/input/samples_taxonomy.csv

OUTPUT_DIR=path/to/out/dir

# Optional filtering thresholds
#MIN_PROT_SIZE=
#MIN_NUC_Size=
#MAX_SEQS=


# Run gene_fetch.py
python gene_fetch.py \
	$GENE \
	$OUTPUT_DIR \
	$SAMPLES_CSV \
	--type $TYPE
