<div align="center">
    <img src="./Gene_fetch_logo.svg" width="300" alt="Gene Fetch Logo">
</div>

# Gene_fetch 
This tool fetches gene sequences from NCBI databases based on taxonomy IDs (taxids) or taxonomic information. It can retrieve both protein and nucleotide sequences for various genes, including protein-coding genes (e.g., cox1, cytb, rbcl, matk) and rRNA genes (e.g., 16S, 18S).

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
 - [Input](#input)
 - [Examples](#Examples)
 - [Output](#output)
 - [Cluster](#running-gene_fetch-on-a-cluster)
 - [Supported targets](#supported-targets)
 - [Benchmarking](#benchmarking)
 - [Contributions and citation](#contributions-and-citations)


## Installation
First, clone the Gene Fetch GitHub repository to your current path, and enter the Gene Fetch installation directory 
```bash
git clone https://github.com/bge-barcoding/gene_fetch

cd gene_fetch
```
Run the commands below to install the necessary dependencies and activate the conda environment. [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) must be installed.
```bash
conda env create -n fetch -f fetch.yaml

conda activate fetch
```

## Usage
```bash
python gene_fetch.py -g/--gene <gene_name> --type <sequence_type> -i/--in <samples.csv> -o/--out <output_directory> 
```
### Options
* `--h/--help`: Show help and exit.
#### Required arguments
* `-g/--gene`: Name of gene to search for in NCBI GenBank database (e.g., cox1/16s/rbcl).
* `--type`: Sequence type to fetch; 'protein', 'nucleotide', or 'both' ('both' will initially search and fetch a protein sequence, and then fetches the corresponding nucleotide CDS for that protein sequence).
* `-i/--in`: Path to input CSV file containing sample IDs and TaxIDs (see [Input](#input) section below).
* `i2/--in2`: Path to alternative input CSV file containing sample IDs and taxonomic information for each sample (see [Input](#input) section below).
* `o/--out`: Path to output directory. The directory will be created if it does not exist.
* `e/--email` and `-k/--api-key`: Email address and associated API key for NCBI account. An NCBI account is required to run this tool (due to otherwise strict API limitations) - information on how to create an NCBI account and find your API key can be found [here](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317).

#### Optional arguments
* `s/--single`: Taxonomic ID for 'single-taxid' sequence search mode (`-i` and `-i2` ignored when run with `-s` mode). 'Single-taxid' mode will fetch all target gene or protein sequences on GenBank for a specific taxonomic ID.
* `--protein_size`: Minimum protein sequence length filter. Applicable to mode 'normal' and 'single-taxid' search modes (default: 500).
* `--nucleotide_size`: Minimum nucleotide sequence length filter. Applicable to mode 'normal' and 'single-taxid' search modes (default: 1500).

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


### Running gene_fetch on a cluster
- See 1_gene_fetch.sh (for running via Slurm).


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
| Sample number | Run mode | target | resources allocated | run time |
| 570 Arthropods | 'normal' | COX1 | 10G, 18 threads | --- |
| ---  | --- | --- | --- | --- |
| --- | --- | --- | --- | --- |
| --- | --- | --- | --- | --- |


## Contributions and citations
GeneFetch was written by Dan Parsons and Ben Price @ NHMUK (2024)

If you use GeneFetch, please cite our publication: XYZ()
