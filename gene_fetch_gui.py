#!/usr/bin/env python
"""
Gene Fetch GUI - A graphical user interface for the Gene Fetch tool

This script provides a user-friendly GUI for the Gene Fetch tool, which fetches gene sequences
from NCBI databases based on taxonomy IDs (taxids) or taxonomic information.

Dependencies:
- Python>=3.9
- Gooey>=1.0.8
- All dependencies of gene_fetch.py (Biopython>=1.80, ratelimit>=2.2.1)

Usage:
    Simply run this script to launch the GUI:
    python gene_fetch_gui.py

Authors: Ben Price and Claude Sonnet 3.7 (AI Assistant)
Version: 1.0
License: MIT
"""

import os
import sys
import subprocess
import argparse
import re
from gooey import Gooey, GooeyParser
import logging
import json

# Set up logging
logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("gene_fetch_gui")

# Define the path to the original gene_fetch.py script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GENE_FETCH_PATH = os.path.join(SCRIPT_DIR, 'gene_fetch.py')

# Define supported genes list
SUPPORTED_GENES = [
    'cox1', 'cox2', 'cox3', 'cytb', 'nd1', 'nd2', 'nd4', 'nd5', 'atp6', 'atp8',
    'rbcl', 'matk', 'psba',
    '16s', '18s', '28s', '12s',
    'its', 'its1', 'its2',
    'trnh-psba', 'trnl'
]

# Function to run the gene_fetch.py script with the given arguments
def run_gene_fetch(args):
    """
    Run the gene_fetch.py script with the given arguments.
    """
    cmd = [sys.executable, GENE_FETCH_PATH]
    
    # Add all the arguments
    for arg_name, arg_value in args.items():
        if arg_value is None or (isinstance(arg_value, bool) and not arg_value):
            continue
            
        if arg_name == 'input_csv':
            cmd.extend(['-i', arg_value])
        elif arg_name == 'input_taxonomy_csv':
            cmd.extend(['-i2', arg_value])
        elif arg_name == 'api_key':
            cmd.extend(['-k', arg_value])
        elif arg_name in ['help', 'h', 'mode']:  # Skip these as they're GUI only
            continue
        else:
            # Convert underscore to hyphen for command line args
            arg_name_cli = arg_name.replace('_', '-')
            if isinstance(arg_value, bool):
                if arg_value:
                    cmd.append(f'--{arg_name_cli}')
            else:
                cmd.extend([f'--{arg_name_cli}' if len(arg_name) > 1 else f'-{arg_name_cli}', str(arg_value)])
    
    # Log the command
    logger.info(f"Running command: {' '.join(cmd)}")
    
    # Run the command with both stdout and stderr redirected to a pipe
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,  # Keep stderr separate from stdout
        universal_newlines=True,
        bufsize=1
    )
    
    # Process output in real-time
    success_message = None
    completion_percentage = 0
    
    # Monitor stdout for progress
    for line in iter(process.stdout.readline, ''):
        # Print the line to the console
        print(line, end='')
        sys.stdout.flush()
        
        # Check for progress patterns - support multiple formats
        if "Progress" in line or "Processing sample" in line:
            # Various progress patterns that might appear in the output
            progress_match = None
            current = None
            total = None
            
            # Pattern 1: "Progress: X/Y"
            match = re.search(r'Progress:\s+(\d+)/(\d+)', line)
            if match:
                current, total = int(match.group(1)), int(match.group(2))
                progress_match = True
            
            # Pattern 2: "===== Progress: X/Y samples processed"
            if not progress_match:
                match = re.search(r'=+\s+Progress:\s+(\d+)/(\d+)\s+samples', line)
                if match:
                    current, total = int(match.group(1)), int(match.group(2))
                    progress_match = True
            
            # Pattern 3: "====== Processing sample X/Y"
            if not progress_match:
                match = re.search(r'=+\s+Processing sample\s+(\d+)/(\d+)', line)
                if match:
                    current, total = int(match.group(1)), int(match.group(2))
                    progress_match = True
            # Calculate progress percentage
            if progress_match and current is not None and total is not None and total > 0:
                completion_percentage = int(100 * current / total)
                print(f"progress: {completion_percentage}%")  # Direct print for Gooey to capture
                logger.info(f"progress: {current}/{total} samples ({completion_percentage}%)")
                
        # Check for completion message
        if "Completed processing:" in line or "completed" in line.lower():
            success_message = line.strip()
    
    # Read any error output
    error_output = process.stderr.read()
    if error_output:
        logger.error(f"Error output from process: {error_output}")
        print(f"Error: {error_output}", file=sys.stderr)
    
    # Wait for the process to finish
    process.wait()
    
    # Check if the process was successful
    if process.returncode != 0:
        logger.error(f"Gene Fetch failed with return code {process.returncode}")
        print(f"Error: Gene Fetch failed with return code {process.returncode}")
        return process.returncode
    else:
        # Log success
        if success_message:
            logger.info(f"Gene Fetch completed successfully: {success_message}")
        else:
            logger.info("Gene Fetch completed successfully!")
            
        # Display notification about completion
        print("\n" + "="*50)
        print("Gene Fetch completed successfully!")
        print("Output files are available in the specified output directory.")
        print("="*50)
        
    return process.returncode

@Gooey(
    program_name="Gene Fetch GUI",
    program_description="GUI for fetching gene sequences from NCBI databases using taxonomy or NCBI taxid",
    image_dir=os.path.join(SCRIPT_DIR),
    default_size=(800, 600),
    required_cols=1,
    optional_cols=2,
    tabbed_groups=True,
    header_bg_color="#3498db",
    body_bg_color="#f9f9f9",
    footer_bg_color="#dcdcdc",
    sidebar_title_color="#ffffff",
    menu=[{
        'name': 'File',
        'items': [{
            'type': 'AboutDialog',
            'menuTitle': 'About',
            'name': 'Gene Fetch GUI',
            'description': 'A graphical user interface for the Gene Fetch tool',
            'version': '1.0',
            'website': 'https://github.com/bge-barcoding/gene_fetch',
            'developer': 'Based on Gene Fetch by Dan Parsons & Ben Price (2025)'
        }]
    }, {
        'name': 'Help',
        'items': [{
            'type': 'Link',
            'menuTitle': 'Documentation',
            'url': 'https://github.com/bge-barcoding/gene_fetch'
        }]
    }],
    # Simplified progress regex that safely handles all output formats
    # progress_regex=r'Progress update: (\d+)%',
    progress_regex=r"^progress: (\d+)%$",
    # progress_expr='int(groups[0])',
    timing=True,
    show_success_modal=True,
    show_failure_modal=True,
    return_to_config=True,
    success_message="Gene Fetch completed successfully!",
    failure_message="Gene Fetch encountered an error. Please check the logs."
)
def main():
    # Create the parser
    parser = GooeyParser(description="Gene Fetch GUI - NCBI Sequence Retrieval Tool")
    
    # Create subparsers for different modes
    subparsers = parser.add_subparsers(dest='mode', help='Select operating mode')
    
    # Batch mode parser
    batch_parser = subparsers.add_parser('batch', help='Batch mode - fetch sequences for multiple samples')
    
    # Single taxid mode parser
    single_parser = subparsers.add_parser('single', help='Single taxid mode - fetch all sequences for a specific taxid')
    
    # Set up the normal mode parser
    setup_batch_mode_parser(batch_parser)
    
    # Set up the single taxid mode parser
    setup_single_taxid_mode_parser(single_parser)
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Convert the args namespace to a dict
    args_dict = vars(args)
    
    # Remove the 'mode' key as it's not needed for gene_fetch.py
    mode = args_dict.pop('mode', 'batch')
    
    # Run the gene_fetch.py script
    return run_gene_fetch(args_dict)

def setup_batch_mode_parser(parser):
    """Set up the parser for batch mode."""
    # Create argument groups for organization
    required_group = parser.add_argument_group('Required Arguments')
    input_group = parser.add_argument_group('Input Options')
    sequence_group = parser.add_argument_group('Sequence Options')
    api_group = parser.add_argument_group('NCBI API Settings')
    
    # Required arguments
    required_group.add_argument(
        '-g', '--gene',
        metavar='Gene Name',
        help='Gene to search NCBI (e.g., cox1, 16s, rbcl)',
        required=True,
        choices=SUPPORTED_GENES,
        default='cox1'
    )
    
    required_group.add_argument(
        '-o', '--out',
        metavar='Output Directory',
        help='Path to directory to save output files',
        required=True,
        widget='DirChooser'
    )
    
    required_group.add_argument(
        '--type',
        metavar='Sequence Type',
        help='Specify sequence type to fetch',
        required=True,
        choices=['protein', 'nucleotide', 'both'],
        default='both'
    )
    
    # Create mutually exclusive group for input files
    input_mutually_exclusive = input_group.add_mutually_exclusive_group(required=True)
    
    input_mutually_exclusive.add_argument(
        '-i', '--input-csv',
        dest='input_csv',
        metavar='Input CSV with TaxIDs',
        help='Input CSV file must have columns: "ID" and "taxid"',
        widget='FileChooser',
        gooey_options={
            'wildcard': "CSV files (*.csv)|*.csv"
        }
    )
    
    input_mutually_exclusive.add_argument(
        '-i2', '--input-taxonomy-csv',
        dest='input_taxonomy_csv',
        metavar='Input CSV with Taxonomy',
        help='Input CSV file must have columns: "ID", "phylum", "class", "order", "family", "genus", "species"',
        widget='FileChooser',
        gooey_options={
            'wildcard': "CSV files (*.csv)|*.csv"
        }
    )
    
    # API settings
    api_group.add_argument(
        '-e', '--email',
        metavar='Email',
        help='Email to use for NCBI API requests (required)',
        required=True
    )
    
    api_group.add_argument(
        '-k', '--api-key',
        metavar='API Key',
        help='API key to use for NCBI API requests (required)',
        required=True
    )
    
    # Sequence options
    sequence_group.add_argument(
        '--protein-size',
        metavar='Min Protein Length',
        help='Default: 500 amino acids',
        type=int,
        default=500
    )
    
    sequence_group.add_argument(
        '--nucleotide-size',
        metavar='Min Nucleotide Length',
        help='Default: 1000 base pairs',
        type=int,
        default=1000
    )

def setup_single_taxid_mode_parser(parser):
    """Set up the parser for single taxid mode."""
    # Create argument groups for organization
    required_group = parser.add_argument_group('Required Arguments')
    sequence_group = parser.add_argument_group('Sequence Options')
    api_group = parser.add_argument_group('NCBI API Settings')
    
    # Required arguments
    required_group.add_argument(
        '-g', '--gene',
        metavar='Gene Name',
        help='Gene to search NCBI (e.g., cox1, 16s, rbcl)',
        required=True,
        choices=SUPPORTED_GENES,
        default='cox1'
    )
    
    required_group.add_argument(
        '-o', '--out',
        metavar='Output Directory',
        help='Path to directory to save output files',
        required=True,
        widget='DirChooser'
    )
    
    required_group.add_argument(
        '--type',
        metavar='Sequence Type',
        help='Sequence type to fetch',
        required=True,
        choices=['protein', 'nucleotide', 'both'],
        default='both'
    )
    
    required_group.add_argument(
        '-s', '--single',
        metavar='Single TaxID',
        help='Will fetch all available sequences',
        required=True
    )
    
    # API settings
    api_group.add_argument(
        '-e', '--email',
        metavar='Email',
        help='Email to use for NCBI API requests (required)',
        required=True
    )
    
    api_group.add_argument(
        '-k', '--api-key',
        metavar='API Key',
        help='API key to use for NCBI API requests (required)',
        required=True
    )
    
    # Sequence options
    sequence_group.add_argument(
        '--protein-size',
        metavar='Min Protein Length',
        help='Default: 100 amino acids in single taxid mode',
        type=int,
        default=100
    )
    
    sequence_group.add_argument(
        '--nucleotide-size',
        metavar='Min Nucleotide Length',
        help='Default: 200 base pairs in single taxid mode',
        type=int,
        default=200
    )
    
    sequence_group.add_argument(
        '--max-sequences',
        metavar='Max Sequences',
        help='Maximum number of sequences to fetch (default: all available)',
        type=int
    )
    
if __name__ == "__main__":
    main()