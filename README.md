<div align="center">
    <img src="./Gene_fetch_logo.svg" width="300" alt="Gene Fetch Logo">
</div>

# Gene_fetch 
This script fetches gene sequences from NCBI databases based on taxonomy IDs (taxids). It can retrieve both protein and nucleotide sequences, with support for various genes including protein-coding genes (e.g., cox1, cox2, cytb, rbcl, matk) and rRNA genes (e.g., 16S, 18S).

## Features
- Fetch protein and/or nucleotide sequences from NCBI databases using taxonomic ID (taxid). Handles both direct nucleotide sequences and protein-linked nucleotide references.
- Support for both protein-coding and rDNA genes.
- Customisable length filtering thresholds
- Automatic taxonomy traversal using NCBI lineage for given taxid when sequences are not found at input taxonomic level. I.e. If sequences are not found at the input taxonomic level (e.g. species), the script searches up higher taxonomic ranks (genus, family, etc.) until one is found.
- By traversing up the fetched NCBI lineage and validating higher taxonomy, potential homonyms are avoided.
- Robust error handling, logging, and NCBI API rate limiting to comply with guidelines (10 requests/second. Requires valid NCBI API key and email for optimal performance).
- Handles complex sequence features (e.g., complement strands, joined sequences, WGS entries) in addition to 'simple' cds extaction (if --type nucleotide/both).
- Single-taxid mode (-s/--single) for retrieving all available sequences for a specific taxon (-i not required)
- 'Checkpointing' available: If a sample fails during a run, the script can be rerun using the same arguments, and it will skip IDs with entries already in the sequence_references.csv and with .fasta files already present in the output directory.


## Contents
 - [Installation](#installation)
 - [Usage](#usage)
 - [Input](#example-input-data)
 - [Output](#output)
 - []()
 - []()
 - []()
 - []()
 - [Supported targets](#support-targets)
 - [Benchmarking](#benchmarking)


## Installation
First, clone the Gene Fetch GitHub repository in your current path
```bash
git clone https://github.com/bge-barcoding/gene_fetch
```

```bash
conda env create -f fetch.yaml
conda activate fetch
```


## Usage
```bash
python gene_fetch.py -g/--gene <gene_name> -o/--out <output_directory> -i/--in <samples.csv> --type <sequence_type>
                        [--protein_size <min_size>] [--nucleotide_size <min_size>] [-s/--single <taxid>] [-i/--in2 <samples_taxonomy.csv>]

Required:
  -e/--email <email_address>: Email address used for NCBI account
  -k/--api-key <key>: API key to use for NCBI API requests
  -g/--gene <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g., cox1/16s/rbcl).
  -o/--out <output_directory>: Path to directory to save output files. The directory will be created if it does not exist.
  -i/--in <samples.csv>: Path to input CSV file containing sample IDs (ID column) and TaxIDs (taxid column).
  -i2/--in2 <samples_taxonomy.csv>: Path to alternative input CSV file containing sample IDs (ID column) and taxonomic heirarchies (phylum, class, order, family, genus, and species columns) for each sample.
  --type: Sequence type to fetch ('protein', 'nucleotide', or 'both')

Optional
--protein_size: Minimum protein sequence length (default: 500).
--nucleotide_size: Minimum nucleotide sequence length (default: 1500).
-s/--single <taxid>: Taxonomic ID for sequence search (-i/--input ignored when -s mode is run).
```

## Input
### Example 'samples.csv' input file
| ID | taxid |
| --- | --- |
| sample-1  | 177658 |
| sample-2 | 177627 |
| sample-3 | 3084599 |

### Example 'samples_taxonomy.csv' input file
| ID | phylum | class | order | family | genus | species |
| --- | --- | --- | --- | --- | --- | --- |
| sample-1  | Arthropoda | Insecta | Diptera | Acroceridae | Astomella | Astomella hispaniae |
| sample-2 | Arthropoda | Insecta | Hemiptera | Cicadellidae | Psammotettix | Psammotettix sabulicola |
| sample-3 | Arthropoda | Insecta | Trichoptera | Limnephilidae | Dicosmoecus | palatus |


## Output
### Output directory structure
#### 'Normal' mode
```
output_dir/
├── nucleotide/             # Nucleotide sequences. Only populated if '--type nucleotide/both' utilised.
│   ├── SAMPLE1_dna.fasta   
│   ├── SAMPLE2_dna.fasta
│   └── ...
├── SAMPLE1.fasta           # Protein sequences.
├── SAMPLE2.fasta
├── sequence_references.csv # Sequence metadata.
├── failed_searches.csv     # Failed search attempts.
└── gene_fetch.log          # Operation log.
```

**sequence_references.csv output example**
| ID | taxid | protein_accession | protein_length | nucleotide_accession | nucleotide_length | matched_rank | ncbi_taxonomy | reference_name | protein_reference_path | nucleotide_reference_path |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BSNHM002-24 | 177658 | AHF21732.1 | 510 | KF756944.1 | 1530 | genus: Apatania | Eukaryota, ..., Apataniinae, Apatania | BSNHM002-24 | abs/path/to/protein_references/BSNHM002-24.fasta | abs/path/to/protein_references/BSNHM002-24_dna.fasta |


#### 'Single-taxid' mode
```
output_dir/
├── nucleotide/                      # Nucleotide sequences. Only populated if '--type nucleotide/both' utilised.
│   ├── ACCESSION1_dna.fasta   
│   ├── ACCESSION2_dna.fasta
│   └── ...
├── ACCESSION1.fasta                 # Protein sequences.
├── ACCESSION2.fasta
├── fetched_nucleotide_sequences.csv # Only populated if '--type nucleotide/both' utilised. Sequence metadata.
├── fetched_protein_sequences.csv    # Only populated if '--type protein/both' utilised. Sequence metadata.
├── failed_searches.csv              # Failed search attempts (ideally empty)
└── gene_fetch.log                   # Operation log
```

**fetched_protein|nucleotide_sequences.csv output example**
| ID | taxid | Description |
| --- | --- | --- |
| BSNHM002-24 | 3086 | cytochrome c oxidase subunit I (mitochondrion) [Pectinodesmus pectinatus] |


### Running gene_fetch.py on a cluster
- See 1_gene_fetch.sh (for running via Slurm).
- 
## Supported targets
- Script functions with other gene/protein targets than those listed below, but has hard-coded synonymns to catch name variations (of the below targets). More targets can be added into script (see 'class config').
- cox1/COI
- cox2/COII
- cox3/COIII
- cytb
- nd1
- rbcl
- matk
- 16s
- 18s
- 28s
- 12s

## Benchmarking
| Sample number | Run mode | target | order | family | genus | species |
| --- | --- | --- | --- | --- | --- | --- |
| BSNHM002-24  | 177658 | --- | --- | --- | --- | --- |
| BSNHM038-24 | 177627 | --- | --- | --- | --- | --- |
| BSNHM046-24 | 3084599 | --- | --- | --- | --- | --- |
