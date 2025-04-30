# Gene Fetch GUI Installation Guide

This guide will help you install and set up the Gene Fetch GUI application.

## Prerequisites

- Python 3.9 or higher
- pip (Python package installer)
- An NCBI account with an API key (see [NCBI API Key setup](#ncbi-api-key-setup))

## Installation Options

### Option 1: Using pip (Recommended)

1. Install the required dependencies:

```bash
pip install -r requirements.txt
```

2. Run the GUI:

```bash
# On Windows:
run_gene_fetch_gui.bat

# On Linux/Mac:
bash run_gene_fetch_gui.sh
```

### Option 2: Using Conda

1. Create a new conda environment:

```bash
conda create -n genefetch-gui python=3.9
conda activate genefetch-gui
```

2. Install the required dependencies:

```bash
pip install -r requirements.txt
```

3. Run the GUI:

```bash
# On Windows:
run_gene_fetch_gui.bat

# On Linux/Mac:
bash run_gene_fetch_gui.sh
```

## NCBI API Key Setup

Gene Fetch requires an NCBI API key to function properly:

1. Create an NCBI account at [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/) if you don't already have one
2. Log in to your NCBI account
3. Go to your account settings by clicking on your username in the top right corner
4. Click on "Account Settings"
5. Scroll down to the "API Key Management" section
6. Click "Create an API Key"
7. Your new API key will be displayed - copy this key for use in Gene Fetch GUI

## Using the GUI

1. Launch the application by running the appropriate script for your platform
2. Choose the operation mode (Normal or Single Taxid)
3. Fill in the required fields:
   - Gene name
   - Output directory
   - Sequence type
   - Input CSV file (Normal mode) or Taxid (Single Taxid mode)
   - NCBI email and API key
4. Adjust any optional parameters as needed
5. Click the "Start" button to begin fetching sequences

## Troubleshooting

- If you encounter any wxPython installation issues, refer to the official wxPython installation guide: [https://wxpython.org/pages/downloads/](https://wxpython.org/pages/downloads/)
- For Gooey installation issues, check the Gooey GitHub repository: [https://github.com/chriskiehl/Gooey](https://github.com/chriskiehl/Gooey)
- If Gene Fetch fails to run, check the log file for detailed error messages

## Support

If you encounter any issues with the Gene Fetch GUI, please submit an issue on the GitHub repository.