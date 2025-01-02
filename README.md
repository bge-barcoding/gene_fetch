# gene_fetch
A Python tool for retrieving protein and/or gene sequences from NCBI databases. The script can fetch both protein and nucleotide sequences for a given gene across multiple taxa, with support for traversing taxonomic hierarchies when sequences aren't available at the given taxonomic level (dictated by input taxid).

- Fetch protein and/or nucleotide sequences from NCBI databases using taxonomic ID (taxid). Handles both direct nucleotide sequences and protein-linked nucleotide references.
- Support for both protein-coding and rRNA genes.
- Automatic taxonomy traversal using NCBI lineage for given taxid when sequences aren't found at input taxonomic level. BY traversing up the fetched NCBI lineage, potential homonyms are avoided.
- Configurable sequence length thresholds, with minimum length requirements (can be altered).
- Robust error handling, logging, and NCBI API rate limiting to comply with guidelines (10 requests/second. Requires valid NCBI API key and email for optimal performance).
- Handles complex sequence features (e.g., complement strands, joined sequences) in addition to 'simple' cds extaction (if --type nucleotide/both).


## Usage
```bash
python 1_gene_fetch.py <gene_name> <output_directory> <samples.csv> --type <sequence_type>
  <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g., cox1/16s/rbcl).
  <output_directory>: Path to directory to save output files. The directory will be created if it does not exist.
  <samples.csv>: Path to input CSV file containing Process IDs (ID column) and TaxIDs (taxid column).
  --type: Sequence type to fetch ('protein', 'nucleotide', or 'both')
  --protein_size: Minimum protein sequence length (default: 500). Optional.
  --nucleotide_size: Minimum nucleotide sequence length (default: 1500). Optional.
```
- 'Checkpointing' available: If the script fails during a run, it can be rerun using the same commands/inputs and it will skip IDs with entries already in the sequence_references.csv and with .fasta files already present in the output directory.

### samples.csv input file example
- Only 'ID' and 'taxid' columns are essential.
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

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
