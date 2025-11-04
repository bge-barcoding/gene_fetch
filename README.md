<p align="center">
  <img src="gene_fetch_logo.svg" width="400" alt="gene_fetch_logo">
</p>

[![PyPI version](https://img.shields.io/pypi/v/gene-fetch.svg)](https://pypi.org/project/gene-fetch/)
[![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gene-fetch/README.html)
[![JOSS DOI](https://joss.theoj.org/papers/10.21105/joss.08456/status.svg)](https://doi.org/10.21105/joss.08456)
[![Github Action test](https://github.com/bge-barcoding/gene_fetch/workflows/Test%20gene-fetch/badge.svg)](https://github.com/bge-barcoding/gene_fetch/actions)
[![Zenodo archive DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16759414.svg)](https://doi.org/10.5281/zenodo.16759414)


# GeneFetch 
Gene Fetch enables high-throughput retreival of sequence data from NCBI's GenBank sequence database based on taxonomy IDs (taxids) or taxonomic heirarchies. It can retrieve both protein and/or nucleotide sequences for various genes, including protein-coding genes (e.g., cox1, cytb, rbcl, matk) and rRNA genes (e.g., 16S, 18S).


## Highlight features
- Fetch protein and/or nucleotide sequences from NCBI's GenBank database without constructing NCBI search terms.
- Handles both direct nucleotide sequences and protein-linked nucleotide searches (CDS extraction includes fallback mechanisms for atypical annotation formats).
- Support for both protein-coding and rDNA genes.
- Seqeunce matches are made by searching for the target gene and/ protein in the GenBank annotation (feature table). 
- Customisable length filtering thresholds for protein and nucleotide sequences (default: protein=500aa. nucleotide=1000bp).
- Default "batch" mode processes multiple input taxa based on a user-specified CSV file.
- Configurable "single" mode (-s/--single) for retrieving a specified number of target sequences for a particular taxon (default length thresholds can be bypassed by setting the value to zero or a negative number).
- Automatic taxonomy traversal: Uses fetched NCBI taxonomic lineage for a given taxid when sequences are not found at the input taxonomic level. i.e., Search at given taxid level (e.g., species), if no sequences are found, escalate species->phylum until a suitable sequence is found.
- Taxonomic validation: validates fetched sequence taxonomy against input taxonomic heirarchy, avoiding potential taxonomic homonyms (i.e. when the same taxa name is used for different taxa across the tree of life).
- Robust error handling, progress tracking, and logging, with compliance to NCBI API rate limits (10 requests/second). Caches taxonomy lookups for reduced API calls.
- Handles complex sequence features (e.g., complement strands, joined sequences, WGS entries) in addition to 'simple' cds extaction (if --type nucleotide/both). The tool avoids "unverified" sequences and WGS entries not containing sequence data (i.e. master records).
- 'Checkpointing': if a run fails/crashes, gene-fetch can be rerun using the same arguments and parameters, and it will resume from where it stopped (unless `--clean` is specified).
- When more than 50 matching GenBank records are found for a sample, the tool fetches summary information for all matches (using NCBI esummary API), orders the records by sequence length, and processes the longest sequences first.
- Can output corresponding genbank (.gb) files for each fetched nucleotide and/or protein sequences
- Optional detail in FASTA sequence headers of retrieved sequences.

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
- Due to the risk of dependency conflicts, it's recommended to install Gene Fetch in a Conda environment.
- First Conda needs to be installed, which can be done from [here](https://www.anaconda.com/docs/getting-started/miniconda/install).
- Once installed:
```bash
# Create new environment
conda create -n gene-fetch

# Activate environment
conda activate gene-fetch
```

- Gene Fetch and all necessary dependencies can then be installed via [Bioconda](https://anaconda.org/bioconda/gene-fetch), [PyPI](https://pypi.org/project/gene-fetch/#description), or by specifying `environment.yaml`:
```bash
# Install via bioconda
conda install bioconda::gene-fetch

# Or, install via pip
pip install gene-fetch

# Or, via environment specification
conda env update --name gene-fetch -f environment.yaml --prune

# Verify installation
gene-fetch --help
```

- If you would rather clone this repository and run a standalone version of Gene Fetch for some reason, you can do that as follows:
```bash
# Clone the repository
git clone https://github.com/bge-barcoding/gene_fetch.git
cd gene_fetch

# Activate conda environment (once created), and install gene-fetch (+ dependencies) via your preferred method. See `environment.yaml` for list of dependencies.

# Run standalone Gene Fetch:
python /path/to/gene_fetch.py [options]

```
  
## Recommended: Testing
- The Gene Fetch package includes some basic tests for each module that we recommend are run after installation.
```bash
# Clone the repository
git clone https://github.com/bge-barcoding/gene_fetch.git
cd gene_fetch

# Install pytest
pip install pytest

# [Optional] Locally install Gene Fetch in editable mode from source (when inside `gene_fetch`) - enables testing of source code in development
pip install -e .

# Run tests
pytest
```
* This will take a few minutes to run the tests. You will get 1 warning regarding API credentials as these are not provided in the basic tests.

## Usage
```bash
gene-fetch --gene <gene_name> --type <sequence_type> --in <samples.csv> --out <output_directory> --email example@example.co.uk --api-key 1234567890
```
* `--help`: Show usage help and exit.

### Required arguments
* `-g/--gene`: Name of gene to search for in NCBI GenBank database (e.g., cox1/16s/rbcl).
* `-t/--type`: Sequence type to fetch; 'protein', 'nucleotide', or 'both' ('both' will initially search and fetch a protein sequence, and then fetches the corresponding nucleotide CDS for that protein sequence).
* `-i/--in`: Path to input CSV file containing sample IDs and TaxIDs (see [Input](#input) section below).
* `-i2/--in2`: Path to alternative input CSV file containing sample IDs and taxonomic information for each sample (see [Input](#input) section below).
* `o/--out`: Path to output directory. The directory will be created if it does not exist.
* `e/--email` and `-k/--api-key`: Email address and associated API key for NCBI account. An NCBI account is required to run this tool (due to otherwise strict API limitations) - information on how to create an NCBI account and find your API key can be found [here](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317).
### Optional arguments
* `-ps/--protein-size`: Minimum protein sequence length filter. Applicable to mode 'batch' and 'single' search modes (default: 500).
* `-ns/--nucleotide-size`: Minimum nucleotide sequence length filter. Applicable to mode 'batch' and 'single' search modes (default: 1000).
* `s/--single`: Taxonomic ID for 'single' sequence search mode (`-i` and `-i2` are ignored when run with `-s` mode). 'single' mode will fetch all (or N if specifying `--max-sequences`) target gene or protein sequences on GenBank for a specific taxonomic ID.
* `-ms/--max-sequences`: Maximum number of sequences to fetch for a specific taxonomic ID (only applies when run in 'single' mode).
* `-b/--genbank`: Saves genbank (.gb) files for fetched nucleotide and/or protein sequences to `genbank/` (applies when run in 'batch' or 'single' mode).
* `-c/--clear`: Forces clean (re)start by clearing output directory regardless of previous run parameters. If ommiting `--clear` and rerunning gene-fetch with the same arguments and parameters, checkpointing will be enabled.
* `--header`: Dictates the format of sequence headers in output FASTA files. 'basic' = '>ID' (default). 'detailed' = '>ID|taxid|accession_number|genbank_description|length'.


## Examples
Fetch both protein and nucleotide sequences for COI with default sequence length thresholds, and store the corresponding genbank records.
```
gene-fetch -e your.email@domain.com -k your_api_key \
            -g cox1 -o ./output_dir -i ./data/samples.csv \
            --type both --genbank
```

Fetch COI nucleotide sequences using sample taxonomic information, applying a minimum nucleotide sequence length of 1000bp
```
gene-fetch -e your.email@domain.com -k your_api_key \
            -g cox1 -o ./output_dir -i2 ./data/samples_taxonomy.csv \
            --type nucleotide --nucleotide-size 1000
```

Retrieve 100 available rbcL protein sequences >400aa for _Arabidopsis thaliana_ (taxid: 3702).
```
gene-fetch -e your.email@domain.com -k your_api_key \
            -g rbcL -o ./output_dir -s 3702 \
            --type protein --protein-size 400 --max-sequences 100
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
| sample-1  | Arthropoda | Insecta | Diptera | Acroceridae | Astomella | |
| sample-2 | Arthropoda | Insecta | Hemiptera | Cicadellidae | Psammotettix | Psammotettix sabulicola |
| sample-3 | Arthropoda | Insecta | Trichoptera | Limnephilidae | Dicosmoecus | Dicosmoecus palatus |
* Leave blank if taxonomic information not known/needed. At least one rank must be supplied for each sample.

## Output
### 'Batch' mode
```
output_dir/
├── genbank/                    # Genbank (.gb) files for each fetched nucleotide and/or protein sequence.
│   ├── nucleotide/  
│   ├── protein/  
├── nucleotide/                 # Nucleotide sequences. Only populated if '--type nucleotide/both' utilised.
│   ├── sample-1.fasta   
│   ├── sample-2.fasta
│   └── ...
├── protein/                    # Protein sequences. Only populated if '--type protein/both' utilised.
│   ├── sample-1.fasta   
│   ├── sample-2.fasta
│   └── ...
├── sequence_references.csv     # Sequence metadata.
├── failed_searches.csv         # Failed search attempts (if any).
└── gene_fetch.log              # Log.
```

**sequence_references.csv output example**
| ID | input_taxa | first_matched_taxid | first_matched_taxid_rank | protein_accession | protein_length | nucleotide_accession | nucleotide_length | matched_rank | ncbi_taxonomy | reference_name | protein_reference_path | nucleotide_reference_path |
| --- | --- | --- | --- | --- | --- | ---| --- | --- | --- | --- | --- | --- |
| sample-1 | Apatania | 177658 | genus:Apatania | AHF21732.1 | 510 | KF756944.1 | 1530 | genus:Apatania | Eukaryota; ...; Apataniinae; Apatania | sample-1 | abs/path/to/protein_references/sample-1.fasta | abs/path/to/protein_references/sample-1_dna.fasta |
| sample-2 | Isoptena serricornis | 2719103 | species:Isoptena serricornis | QNE85983.1 | 518 | MT410852.1 | 1557 | species:Isoptena serricornis | Eukaryota; ...; Chloroperlinae; Isoptena | sample-2 | abs/path/to/protein_references/sample-2.fasta | abs/path/to/protein_references/sample-2_dna.fasta |
| sample-3 | Triaenodes conspersus | 1876143 | species:Triaenodes conspersus | YP_009526503.1 | 512 | NC_039659.1 | 1539 | genus:Triaenodes | Eukaryota; ...; Triaenodini; Triaenodes | sample-3 | abs/path/to/protein_references/sample-3.fasta | abs/path/to/protein_references/sample-3_dna.fasta |
```
* process_id - The unique identifier (ID) for each sample (from the input CSV)
* input_taxa - The taxon name searched for (e.g., "Apatania" == taxid 177658), or the taxon name for the closest valid taxid found.
* first_matched_taxid - The NCBI taxonomic ID that was searched (same as the taxid from the --in CSV, or the closest valid taxid if using --in2 as input)
* first_matched_taxid_rank - The taxonomic rank and name of the first_matched_taxid (e.g., "genus:Astomella")
* protein_accession - The NCBI accession number of the protein sequence retrieved (if applicable)
* protein_length - Length of the protein sequence in amino acids (if applicable)
* nucleotide_accession - The NCBI accession number of the nucleotide sequence retrieved (if applicable)
* nucleotide_length - Length of the nucleotide sequence in base pairs (if applicable)
* matched_rank - The taxonomic rank where sequences were actually found (e.g., "family:Acroceridae" if no sequences existed at the the proceeding rank, and the search traversed up the taxonomy tree)
* ncbi_taxonomy - The complete NCBI taxonomic lineage for the retrieved sequence (semicolon-separated)
* reference_name - Copy of the process_id (for reference purposes)
* protein_reference_path - Full file path to the saved protein FASTA file (if applicable)
* nucleotide_reference_path - Full file path to the saved nucleotide FASTA file (if applicable)
```

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
├── fetched_nucleotide_sequences.csv # Sequence metadata. Only populated if '--type nucleotide/both' utilised.
├── fetched_protein_sequences.csv    # Sequence metadata. Only populated if '--type protein/both' utilised.
├── failed_searches.csv              # Failed search attempts (if any).
└── gene_fetch.log                   # Log.
```

**fetched_protein|nucleotide_sequences.csv output example**
| ID | length | Description | searched_taxid
| --- | --- | --- | --- |
| PQ645072.1 | 1501 | Ochlerotatus nigripes isolate Pool11 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial | 508662 |
| PQ645071.1 | 1537 | Ochlerotatus nigripes isolate Pool10 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial | 508662 |
| PQ645070.1 | 1501 | Ochlerotatus impiger isolate Pool2 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial | 508662 |
| PQ645069.1 | 1518	| Ochlerotatus impiger isolate Pool1 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial | 508662 |

## Running GeneFetch on a cluster
- See 'gene_fetch.sh' for running gene_fetch.py on a HPC cluster (SLURM job schedular). 
- Edit 'mem' and/or 'cpus-per-task' to set memory and CPU/threads - allocating lots of CPUs is unecessary as Gene Fetch is not paralellised (yet). The tool should run well with 4-10G memory and 1-2 CPUs.
- Change paths and variables as needed.
- Run 'gene_fetch.sh' with:
```
sbatch gene_fetch.sh
```

## Supported targets
GeneFetch includes the following 'hard-coded' search terms with common name variations for 'smarter' searching of the targets listed below. 
- cox1/COI/cytochrome c oxidase subunit I
- cox2/COII/cytochrome c oxidase subunit II
- cox3/COIII/cytochrome c oxidase subunit III
- cytb/cob/cytochrome b
- nd1/NAD1/NADH dehydrogenase subunit 1
- nd2/NAD2/NADH dehydrogenase subunit 2
- rbcL/RuBisCO/ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit
- matK/maturase K
- psbA/photosystem II protein D1
- 16S ribosomal RNA/16s
- SSU/18s
- LSU/28s
- 23s
- 12S ribosomal RNA/12s
- ITS (ITS1-5.8S-ITS2)
- ITS1/internal transcribed spacer 1
- ITS2/internal transcribed spacer 2
- tRNA-Leucine/trnL


Gene/protein targets not listed can also be searched, however, Gene Fetch will implement a more generic search term/strategy with `{target}[Title] OR {target}[Gene] OR {target}[Protein Name]`.
Additional targets can be added if required - see `self._rRNA_genes` and '`self._protein_coding_genes` dictionaries within 'class config' (in `src/gene_fetch/core.py`) for example search terms to construct your own. You are welcome to open an [Issue](https://github.com/bge-barcoding/gene_fetch/issues/new) or create a pull request with your search term for inclusion into the main Gene Fetch release (see [Contributions and guidelines](https://github.com/bge-barcoding/gene_fetch?tab=readme-ov-file#contributions-and-guidelines) section below.

## Benchmarking
| Sample Description | Run Mode | Target | Input File | Data Type | Memory | CPUs | Run Time (hh:mm:ss) |
|--------------------|----------|--------|------------|-----------|--------|------|----------|
| 570 Arthropod samples | Batch | COI | taxonomy.csv | Both | 4G | 1 | 01:34:47 |
| 570 Arthropod samples | Batch | COI | samples.csv | Both (+ genbank) | 4G | 1 | 01:42:37 |
| 570 Arthropod samples | Batch | COI | samples.csv | Nucleotide | 4G | 1 | 1:07:53  |
| 570 Arthropod samples | Batch | ND1 | samples.csv | Nucleotide (>500bp) | 4G | 1 | 1:23:26 |
| All available (30) _A. thaliana_ sequences | Single | rbcL | N/A | Protein (>300aa) | 4G | 1 | 00:00:25 |
| 1000 Culicidae sequences | Single | COI | N/A | nucleotide (>500bp) | 4G | 1 | 0031:05 |
| 1000 _M. tubercolisis_ sequences | Single | 16S | N/A | nucleotide | 4G | 1 | 01:23:54 |
* All benchmarking runs were performed on a SLURM-managed HPC cluster running Debian 12 ("Bookworm), with each job allocated a modest 1 CPU and 4 GB RAM.

## Future Development
- Add optional alignment of retrieved sequences [Ben].
- Further improve efficiency of record searching and selecting the longest sequence [Dan].
- Add support for additional genetic markers beyond the currently supported set [Dan].
- Add BOLD query falback if no 'quality' sequence is found in GenBank [Ben].
- Add optional HMM profile alignment that will attempt to extract the barcode region from certain support target genes (e.g. 658bp COI-5P barcode) [Ben].


## Contributions and guidelines
First off, thanks for taking the time to contribute! ❤️

- If you hav any questions, we assume that you have read the available [Documentation](https://github.com/bge-barcoding/gene_fetch/blob/main/README.md). It may also be worth searching for existing [Issues](https://github.com/bge-barcoding/gene_fetch/issues) that might awnser your question(s).
- If you feel you still need clarification or want to report a possible bug/unexpected behaviour, we recommend opening an [Issue](https://github.com/bge-barcoding/gene_fetch/issues/) and provide as much context as you can about what behaviour you were expecting and the behaviour you're running into.
- If you want to suggest a novel feature or minor improvements to existing functionality, please make your case for the feature/enchanment by opening an [Issue](https://github.com/bge-barcoding/gene_fetch/issues/new) or create a pull request with your contribution (at which point it will be evaluated as a possible addition). We aim to address any issues as soon as possible.

## Authorship & citation
GeneFetch was written by Dan Parsons & Ben Price @ NHMUK (2025).

If you use GeneFetch, please cite our publication: Parsons and Price (2025). Gene Fetch: A Python tool for sequence retrieval from GenBank across the tree of life. Journal of Open Source Software, 10(112), 8456, https://doi.org/10.21105/joss.08456