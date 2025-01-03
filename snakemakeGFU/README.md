This Snakemake workflow automates the fetching and processing of genetic sequences from NCBI, suitable for both small and large-scale sequence retrieval tasks, with built-in safeguards for API rate limits and data validation. The workflow supporting two main modes:
- Gene-specific sequence retrieval (e.g., barcoding genes).
- Organelle/ribosomal sequence retrieval (e.g., mitochondrial, chloroplast, or ribosomal).

## Key Features:
- Configurable via YAML file
- Supports batch processing via samples CSV
- Implements error handling and retry logic necessary for robust NCBI API usage
- Maintains organised output directory structure (see below)
- Handles temporary files with cleanup
- Provides comprehensive per-sample logging

## Workflow Structure:
- Input Requirements:
  - Configuration file (config.yaml)
  - Samples CSV file with ID and taxid columns
  - Email and API key for NCBI access
  - Conda environment with required dependencies

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

### Output Organisation (for fetch_gene_sequences/ gene_fetch.py):
```
results/
├── {gene}/
│   └── {run_name}/
│       ├── sequence_references.csv
│       └── gene_fetch.log
└── logs/
    └── gene_fetch-{run_name}-{gene}.log
```
