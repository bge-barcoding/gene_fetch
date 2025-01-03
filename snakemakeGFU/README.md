# Gene Fetch Ultra
A robust Snakemake workflow for automated retrieval and processing of genetic sequences from NCBI. Designed for both small and large-scale sequence retrieval tasks, with built-in safeguards for API rate limits and comprehensive data validation.

## Overview
This workflow supports two modes of operation:
1. **Gene-specific sequence retrieval**
   - Ideal for barcoding genes and similar targeted sequence retrieval
   - Configurable sequence criteria and filtering
2. **Organelle/ribosomal sequence retrieval**
   - Supports mitochondrial, chloroplast, and ribosomal sequences
   - Automated sequence validation and organisation


## Key Features:
- 🛠️ Configurable via YAML file
- 📊 Batch processing through samples CSV
- 🔄 Robust NCBI API error handling and retry logic
- 📁 Organised output directory structure (see below)
- 🧹 Temporary file management
- 📝 Comprehensive per-sample logging


## Prerequisites
- Activate conda environment (created with fetch.yaml)
- NCBI account (for API key)
- Conda or Mamba package manager


## Installation
Clone the repository:
   ```bash
   git clone https://github.com/yourusername/gene-fetch-ultra.git
   cd gene-fetch-ultra
  conda env create -f envs/fetch.yaml
  conda activate fetch
  ```


## Workflow Structure:
- Input Requirements:
  - Configuration file (config/config.yaml) containing workflow parameters, API credentials, and target sequence criteria
  - Samples CSV file with ID and taxid columns
  - Email and API key for NCBI access
  - Conda environment with required dependencies
- Run using run_snakefile.sh (if using cluster), or ```bash
    snakemake --use-conda --cores <n>
    ```

## Main Rules:
a) fetch_gene_sequences (using gene_fetch.py):
   - Handles gene-specific sequence retrieval
   - Outputs sequence references and logs per gene
b) fetch_organelle_sequences:
   - Processes organelle/ribosomal sequences
   - Creates temporary working directories for reorganised (relative to original go_fetch.py) output directory structure
   - Implements 3-retry logic with random backoff for NCBI API usage
   - Validates outputs and moves to final location
   - Maintains global and sample-specific logs


### Output Organisation (for fetch_organelle_sequences/go_fetch.py):
```
results/
├── {target_type}/
│   └── {run_name}/
│       └── {ID}/
│           └── {taxid}/
│               ├── seed.fasta (or fasta directory)
│               ├── gene.fasta (or genbank directory)
│               ├── annotated_regions/
│               └── go_fetch.log
└── logs/
    └── go_fetch-{run_name}-{target_type}-{ID}-{taxid}.log
```

### Output Organisation (for fetch_gene_sequences/gene_fetch.py):
```
results/
├── {gene}/
│   └── {run_name}/
│       ├── sequence_references.csv
│       └── gene_fetch.log
└── logs/
    └── gene_fetch-{run_name}-{gene}.log
```
