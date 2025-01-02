#!/bin/bash
#SBATCH --job-name=GF
#SBATCH --partition=day
#SBATCH --output=%x.out
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --error=%x.err


#activate conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate fetch


GENE=cox1

TYPE=both

OUTPUT_DIR=path/to/out/dir

SAMPLES_CSV=path/to/input/samples.csv

python 1_gene_fetch.py \
	$GENE \
	$OUTPUT_DIR \
	$SAMPLES_CSV \
	--type $TYPE
