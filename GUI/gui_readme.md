# Gene Fetch GUI

A graphical user interface for [Gene Fetch](https://github.com/bge-barcoding/gene_fetch), a tool that fetches gene sequences from NCBI databases based on taxonomy IDs (taxids) or taxonomic information.

## Features

- User-friendly interface for all Gene Fetch functions
- Two operating modes:
  - **Normal mode:** Fetch sequences for multiple samples from a CSV file
  - **Single taxid mode:** Fetch all available sequences for a specific taxid
- Supports both protein and nucleotide sequence retrieval
- Customizable sequence length thresholds
- Real-time progress tracking
- Detailed logging
- Completion notifications

## Installation

See [INSTALLATION.md](INSTALLATION.md) for detailed installation instructions.

Quick start:
```bash
# Install dependencies
pip install -r requirements.txt

# Run the application
# On Windows:
run_gene_fetch_gui.bat

# On Linux/Mac:
bash run_gene_fetch_gui.sh
```

## Usage

See [USER_GUIDE.md](USER_GUIDE.md) for detailed usage instructions.

## Supported Genes

The GUI supports retrieval of sequences for the following genes:

- **Protein-coding genes:**
  - cox1/COI (cytochrome c oxidase subunit I)
  - cox2/COII (cytochrome c oxidase subunit II)
  - cox3/COIII (cytochrome c oxidase subunit III)
  - cytb/cob (cytochrome b)
  - nd1/NAD1 (NADH dehydrogenase subunit 1)
  - nd2/NAD2 (NADH dehydrogenase subunit 2)
  - rbcL (ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit)
  - matK (maturase K)
  - And more...

- **rRNA genes:**
  - 16S ribosomal RNA
  - 18S ribosomal RNA (SSU)
  - 28S ribosomal RNA (LSU)
  - 12S ribosomal RNA
  - And more...

- **Other markers:**
  - ITS (internal transcribed spacer)
  - ITS1/ITS2
  - trnL (tRNA-Leucine)
  - trnh-psba (intergenic spacer)

## Requirements

- Python 3.9 or higher
- Gooey 1.0.8 or higher
- Biopython 1.80 or higher
- An NCBI account with an API key

## Files Included

- **gene_fetch_gui.py** - The main GUI application
- **run_gene_fetch_gui.bat** - Windows batch file for easy execution
- **run_gene_fetch_gui.sh** - Shell script for Linux/Mac users
- **requirements.txt** - List of dependencies
- **USER_GUIDE.md** - Detailed usage guide
- **INSTALLATION.md** - Installation instructions

## Troubleshooting

If you encounter issues:

1. **Missing dependencies**: Run `pip install -r requirements.txt` to install all required packages
2. **Progress bar not updating**: Make sure you're using the latest version of Gooey
3. **Error launching**: Try running `python gene_fetch_gui.py` directly to see error messages

## Acknowledgements

This GUI is built on top of the Gene Fetch tool created by Dan Parsons & Ben Price @ NHMUK (2024).

## License

MIT

## Screenshots

[Screenshots will be added here] of Gooey supports the progress update patterns
4. **Error launching**: Try running `python gene_fetch_gui.py` directly to see error messages

## Acknowledgements

This GUI is built on top of the Gene Fetch tool created by Dan Parsons & Ben Price @ NHMUK (2024).

## License

MIT

## Screenshots

[Screenshots will be added here]