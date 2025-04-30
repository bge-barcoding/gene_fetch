<p align="center">
  <img src="gene_fetch_logo.svg" width="400" alt="gene_fetch_logo">
</p>

# GeneFetch 
This tool enables high-throughput retreival of sequence data from NCBI databases based on taxonomy IDs (taxids) or taxonomic heirarchies. It can retrieve both protein and/or nucleotide sequences for various genes, including protein-coding genes (e.g., cox1, cytb, rbcl, matk) and rRNA genes (e.g., 16S, 18S).


## Highlight features
- Fetch protein and/or nucleotide sequences from NCBI GenBank database.
- Handles both direct nucleotide sequences and protein-linked nucleotide searches (CDS extraction includes fallback mechanisms for atypical annotation formats).
- Support for both protein-coding and rDNA genes.
- Default "batch" mode processes multiple input taxa based on a user specified CSV file.
- Configurable "single" mode (-s/--single) for retrieving a specified number of target sequences for a particular taxon, with default length thresholds reduced (protein: 50aa, nucleotide: 100bp).
- Customisable length filtering thresholds for protein and nucleotide sequences.
- Automatic taxonomy traversal: Uses fetched NCBI taxonomic lineage for a given taxid when sequences are not found at the input taxonomic level. i.e., Search at given taxid level (e.g., species), if no sequences are found, escalate species->phylum until a suitable sequence is found.
- Taxonomic validation: validates fetched sequence taxonomy against input taxonomic heirarchy, avoiding potential taxonomic homonyms (i.e. when the same taxon name is used for different taxa across the tree of life).
- Robust error handling, progress tracking, and logging, with compliance to NCBI API rate limits (10 requests/second). Caches taxonomy lookups for reduced API calls.
- Handles complex sequence features (e.g., complement strands, joined sequences, WGS entries) in addition to 'simple' cds extaction (if --type nucleotide/both). The tool avoids "unverified" sequences and WGS entries not containing sequence data (i.e. master records).
- 'Checkpointing': if a run fails/crashes, the script can be rerun using the same arguments and it will resume from where it stopped.
- When more than 50 matching GenBank records are found for a sample, the tool fetches summary information for all matches (using NCBI esummary API), orders the records by sequence length, and processes the longest sequences first.
- Can output corresponding genbank (.gb) files for each fetched nucleotide and/or protein sequences

## Contents
 - [Installation](#installation)
 - [Usage](#usage)
 - [Examples](#Examples)
 - [Input](#input)
 - [Output](#output)
 - [Cluster](#running-gene_fetch-on-a-cluster)
 - [Supported targets](#supported-targets)
 - [Notes](#notes)
 - [Benchmarking](#benchmarking)
 - [Future developments](#future-developments)
 - [Contributions and citation](#contributions-and-citations)


## Installation
First, clone the GeneFetch GitHub repository to your current path, and enter the Gene Fetch installation directory 
```bash
git clone https://github.com/bge-barcoding/gene_fetch

cd gene_fetch
```
Run the commands below to install the necessary dependencies and activate the Conda environment. [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) must be installed.
```bash
conda env create -n fetch -f fetch.yaml

conda activate fetch
```
Alternatively, you can install the dependencies below directly or in your own Conda environment
```
conda install python>=3.9 pip
pip install ratelimit>=2.2.1
pip install biopython>=1.80
```

## Usage
```bash
python gene_fetch.py -g/--gene <gene_name> --type <sequence_type> -i/--in <samples.csv> -o/--out <output_directory> 
```
* `--h/--help`: Show help and exit.
### Required arguments
* `-g/--gene`: Name of gene to search for in NCBI GenBank database (e.g., cox1/16s/rbcl).
* `--type`: Sequence type to fetch; 'protein', 'nucleotide', or 'both' ('both' will initially search and fetch a protein sequence, and then fetches the corresponding nucleotide CDS for that protein sequence).
* `-i/--in`: Path to input CSV file containing sample IDs and TaxIDs (see [Input](#input) section below).
* `i2/--in2`: Path to alternative input CSV file containing sample IDs and taxonomic information for each sample (see [Input](#input) section below).
* `o/--out`: Path to output directory. The directory will be created if it does not exist.
* `e/--email` and `-k/--api-key`: Email address and associated API key for NCBI account. An NCBI account is required to run this tool (due to otherwise strict API limitations) - information on how to create an NCBI account and find your API key can be found [here](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317).
### Optional arguments
* `--protein-size`: Minimum protein sequence length filter. Applicable to mode 'batch' and 'single' search modes (default: 500).
* `--nucleotide-size`: Minimum nucleotide sequence length filter. Applicable to mode 'batch' and 'single' search modes (default: 1500).
* `s/--single`: Taxonomic ID for 'single' sequence search mode (`-i` and `-i2` are ignored when run with `-s` mode). 'single' mode will fetch all (or N if specifying `--max-sequences`) target gene or protein sequences on GenBank for a specific taxonomic ID.
* `--max-sequences`: Maximum number of sequences to fetch for a specific taxonomic ID (only applies when run in 'single' mode).
* `-b/--genbank`: Saves genbank (.gb) files for fetched nucleotide and/or protein sequences to `genbank/` (applies when run in 'batch' or 'single' mode).


## Examples
Fetch both protein and nucleotide sequences for COI with default sequence length thresholds, and store the corresponding genbank records.
```
python gene_fetch.py -e your.email@domain.com -k your_api_key \
                    -g cox1 -o ./output_dir -i ./samples.csv \
                    --type both --genbank
```

Fetch rbcL nucleotide sequences using sample taxonomic information, applying a minimum nucleotide sequence length of 1000bp
```
python gene_fetch.py -e your.email@domain.com -k your_api_key \
                    -g rbcl -o ./output_dir -i2 ./taxonomy.csv \
                    --type nucleotide --nucleotide_size 1000
```

Retrieve 1000 available matK protein sequences >400aa for _Arabidopsis thaliana_ (taxid: 3702).
```
python gene_fetch.py -e your.email@domain.com -k your_api_key \
                    -g matk -o ./output_dir -s 3702 \
                    --type protein --protein_size 400 --max-sequences 1000
```


## Input
**Example 'samples.csv' input file (-i/--in)**
| ID | taxid |
| --- | --- |
| sample-1  | 177658 |
| sample-2 | 177627 |
| sample-3 | 3084599 |

**Example 'samples_taxonomy.csv' input file (-i2/--in2)**
| ID | phylum | class | order | family | genus | species |
| --- | --- | --- | --- | --- | --- | --- |
| sample-1  | Arthropoda | Insecta | Diptera | Acroceridae | Astomella | n/a |
| sample-2 | Arthropoda | Insecta | Hemiptera | Cicadellidae | Psammotettix | Psammotettix sabulicola |
| sample-3 | Arthropoda | Insecta | Trichoptera | Limnephilidae | Dicosmoecus | Dicosmoecus palatus |


## Output
### 'Batch' mode
```
output_dir/
├── genbank/                    # Genbank (.gb) files for each fetched nucleotide and/or protein sequence.
├── nucleotide/                 # Nucleotide sequences. Only populated if '--type nucleotide/both' utilised.
│   ├── sample-1_dna.fasta   
│   ├── sample-2_dna.fasta
│   └── ...
├── sample-1.fasta              # Protein sequences.
├── sample-2.fasta
├── sequence_references.csv     # Sequence metadata.
├── failed_searches.csv         # Failed search attempts (if any).
└── gene_fetch.log              # Log.
```

**sequence_references.csv output example**
| ID | taxid | protein_accession | protein_length | nucleotide_accession | nucleotide_length | matched_rank | ncbi_taxonomy | reference_name | protein_reference_path | nucleotide_reference_path |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sample-1 | 177658 | AHF21732.1 | 510 | KF756944.1 | 1530 | genus:Apatania | Eukaryota; ...; Apataniinae; Apatania | sample-1 | abs/path/to/protein_references/sample-1.fasta | abs/path/to/protein_references/sample-1_dna.fasta |
| sample-2 | 2719103 | QNE85983.1 | 518 | MT410852.1 | 1557 | species:Isoptena serricornis | Eukaryota; ...; Chloroperlinae; Isoptena | sample-2 | abs/path/to/protein_references/sample-2.fasta | abs/path/to/protein_references/sample-2_dna.fasta |
| sample-3 | 1876143 | YP_009526503.1 | 512 | NC_039659.1 | 1539 | genus:Triaenodes | Eukaryota; ...; Triaenodini; Triaenodes | sample-3 | abs/path/to/protein_references/sample-3.fasta | abs/path/to/protein_references/sample-3_dna.fasta |


### 'Single' mode
```
output_dir/
├── genbank/                         # Genbank (.gb) files for each fetched nucleotide and/or protein sequence.
├── nucleotide/                      # Nucleotide sequences. Only populated if '--type nucleotide/both' utilised.
│   ├── ACCESSION1_dna.fasta   
│   ├── ACCESSION2_dna.fasta
│   └── ...
├── ACCESSION1.fasta                 # Protein sequences.
├── ACCESSION2.fasta
├── fetched_nucleotide_sequences.csv # Only populated if '--type nucleotide/both' utilised. Sequence metadata.
├── fetched_protein_sequences.csv    # Only populated if '--type protein/both' utilised. Sequence metadata.
├── failed_searches.csv              # Failed search attempts (if any).
└── gene_fetch.log                   # Log.
```

**fetched_protein|nucleotide_sequences.csv output example**
| ID | taxid | Description |
| --- | --- | --- |
| PQ645072.1 | 1501 | Ochlerotatus nigripes isolate Pool11 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial |
| PQ645071.1 | 1537 | Ochlerotatus nigripes isolate Pool10 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial |
| PQ645070.1 | 1501 | Ochlerotatus impiger isolate Pool2 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial |
| PQ645069.1 | 1518	| Ochlerotatus impiger isolate Pool1 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial |
| PP355486.1 | 581 | Aedes scutellaris isolate NC.033 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial |


## Running GeneFetch on a cluster
- See '1_gene_fetch.sh' for running gene_fetch.py on a HPC cluster (SLURM job schedular). 
- Edit 'mem' and/or 'cpus-per-task' to set memory and CPU/threads allocation.
- Change paths and variables as needed.
- Run '1_gene_fetch.sh' with:
```
sbatch 1_gene_fetch.sh
```

## Supported targets
GeneFetch will function with other targets than those listed below, but it has hard-coded name variations and 'smarter' searching for the listed targets. More targets can be added into the script if necessary (see 'class config').
- cox1/COI/cytochrome c oxidase subunit I
- cox2/COII/cytochrome c oxidase subunit II
- cox3/COIII/cytochrome c oxidase subunit III
- cytb/cob/cytochrome b
- nd1/NAD1/NADH dehydrogenase subunit 1
- nd2/NAD2/NADH dehydrogenase subunit 2
- rbcL/RuBisCO/ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit
- matK/maturase K/maturase type II intron splicing factor
- 16S ribosomal RNA/16s
- SSU/18s
- LSU/28s
- 12S ribosomal RNA/12s
- ITS (ITS1-5.8S-ITS2)
- ITS1/internal transcribed spacer 1
- ITS2/internal transcribed spacer 2
- tRNA-Leucine/trnL


## Benchmarking
| Sample Description | Run Mode | Target | Input File | Data Type | Memory | CPUs | Run Time |
|--------------------|----------|--------|------------|-----------|--------|------|----------|
| 570 Arthropod samples | Batch | COX1 | taxonomy.csv | Both | 10G | 16 | 1:39:33 |
| 570 Arthropod samples | Batch | COX1 | samples.csv | Nucleotide | 5G | 4 | 02:04:01 |
| 570 Arthropod samples | Batch | COX1 | samples.csv | Protein | 5G | 4 | 01:50:31 |
| 570 Arthropod samples | Batch | 18S | samples.csv | Nucleotide | 10G | 8 | 01:38:16 |
| 570 Arthropod samples | Batch | ND1 | samples.csv | Nucleotide | 10G | 4 | 01:58:35 |
| All (159) _A. thaliana_ sequences >300aa | Single | rbcL | N/A | Protein | 5G | 1 | 00:02:39 |
| 1000 Culicidae sequences >500bp | Single | COX1 | N/A | nucleotide | 20G | 16 | 00:30:36 |
| 1000 _M. tubercolisis_ sequences | Single | 16S | N/A | nucleotide | 20G | 16 | 00:10:33 |

## Future Development
- Add optional alignment of retrieved sequences
- Further improve efficiency of record searching and selecting the longest sequence
- Add support for additional genetic markers beyond the currently supported set


## Contributions and citations
GeneFetch was written by Dan Parsons & Ben Price @ NHMUK (2025).

If you use GeneFetch, please cite our publication: **XYZ()**

If you have any questions or suggested improvements, please do get in touch in the issues!
