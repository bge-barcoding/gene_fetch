<div align="center">
    <img src="./Gene_fetch_logo.svg" width="300" alt="Gene Fetch Logo">
</div>

# Gene_fetch 
This script fetches gene sequences from NCBI databases based on taxonomy IDs (taxids). It can retrieve both protein and nucleotide sequences, with support for various genes including protein-coding genes (e.g., cox1, cox2, cytb, rbcl, matk) and rRNA genes (e.g., 16S, 18S).

- Fetch protein and/or nucleotide sequences from NCBI databases using taxonomic ID (taxid). Handles both direct nucleotide sequences and protein-linked nucleotide references.
- Support for both protein-coding and rRNA genes.
- Customisable length filtering thresholds
- Taxonomic traversal: If sequences are not found at the input taxonomic level (e.g. species), searches up higher taxonomic ranks (genus, family, etc.)
- Automatic taxonomy traversal using NCBI lineage for given taxid when sequences are not found at input taxonomic level. By traversing up the fetched NCBI lineage and validating higher taxonomy, potential homonyms are avoided.
- Robust error handling, logging, and NCBI API rate limiting to comply with guidelines (10 requests/second. Requires valid NCBI API key and email for optimal performance).
- Handles complex sequence features (e.g., complement strands, joined sequences, WGS entries) in addition to 'simple' cds extaction (if --type nucleotide/both).
- Single-taxid mode (-s/--single) for retrieving all available sequences for a specific taxon (-i not required)


## Usage
```bash
python 1_gene_fetch.py -g/--gene <gene_name> -o/--out <output_directory> -i/--in <samples.csv> --type <sequence_type>
                        [--protein_size <min_size>] [--nucleotide_size <min_size>] [-s/--single <taxid>]
  -g/--gene <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g., cox1/16s/rbcl).
  -o/--out <output_directory>: Path to directory to save output files. The directory will be created if it does not exist.
  -i/--in <samples.csv>: Path to input CSV file containing Process IDs (ID column) and TaxIDs (taxid column).
  --type: Sequence type to fetch ('protein', 'nucleotide', or 'both')
  Optional: --protein_size: Minimum protein sequence length (default: 500).
  Optional: --nucleotide_size: Minimum nucleotide sequence length (default: 1500).
  Optional: -s/--single <taxid>: Taxonomic ID for sequence search (-i/--input ignored when -s mode is run(.
```
- 'Checkpointing' available: If a sample fails during a run, the script can be rerun using the same arguments, and it will skip IDs with entries already in the sequence_references.csv and with .fasta files already present in the output directory.
- ./snakemakeSFU/workflow/envs/fetch.yaml contains all necessary dependencies to run the script. Can create conda env using commmand below (once conda is installed on your system).
```
conda env create -f fetch.yaml
```


### samples.csv input file example
| ID | taxid |
| --- | --- |
| BSNHM002-24  | 177658 |
| BSNHM038-24 | 177627 |
| BSNHM046-24 | 3084599 |

### Running gene_fetch.py on a cluster
- See 1_gene_fetch.sh (for running via Slurm).

## Output Structure
```
output_dir/
├── nucleotide/             # Only populated if '--type nucleotide/both' utilised.
│   ├── SAMPLE1_dna.fasta   
│   ├── SAMPLE2_dna.fasta
│   └── ...
├── SAMPLE1.fasta           # Protein sequences
├── SAMPLE2.fasta
├── sequence_references.csv # Sequence metadata (for input into MGE snakemake pipeline)
├── failed_searches.csv     # Failed search attempts (ideally empty)
└── gene_fetch.log          # Operation log
```
### sequence_references.csv output example
| ID | taxid | protein_accession | protein_length | nucleotide_accession | nucleotide_length | matched_rank | ncbi_taxonomy | reference_name | protein_reference_path | nucleotide_reference_path |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BSNHM002-24 | 177658 | AHF21732.1 | 510 | KF756944.1 | 1530 | genus: Apatania | Eukaryota, ..., Apataniinae, Apatania | BSNHM002-24 | abs/path/to/protein_references/BSNHM002-24.fasta | abs/path/to/protein_references/BSNHM002-24_dna.fasta |



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
