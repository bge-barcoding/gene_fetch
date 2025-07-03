# src/gene_fetch/main.py
"""
Command-line interface for Gene Fetch.
Provides the entry point and argument parser for the Gene Fetch tool.
"""

import argparse
import csv
import sys
import time
from pathlib import Path
from random import uniform
from time import sleep
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .core import Config, setup_logging, make_out_dir, log_progress, get_process_id_column, logger
from .entrez_handler import EntrezHandler
from .sequence_processor import SequenceProcessor
from .output_manager import OutputManager, save_genbank_file
from .processors import process_sample, process_single_taxid, process_taxid_csv, process_taxonomy_csv


def setup_argument_parser():
    parser = argparse.ArgumentParser(
        description="Fetch gene sequences from NCBI databases."
    )

    parser.add_argument(
        "-g",
        "--gene",
        required=True,
        help="Name of gene to search for in NCBI RefSeq database (e.g., cox1)",
    )

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Path to directory to save output files (will create new directories)",
    )

    # Create mutually exclusive group for input files
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument(
        "-i",
        "--in",
        dest="input_csv",
        help="Path to input CSV file containing taxIDs (must have columns "
        '"taxID" and "ID")',
    )
    input_group.add_argument(
        "-i2",
        "--in2",
        dest="input_taxonomy_csv",
        help="Path to input CSV file containing taxonomic information "
        '(must have columns "ID", "phylum", "class", "order", '
        ' "family", "genus", "species")',
    )

    parser.add_argument(
        "-s",
        "--single",
        type=str,
        help="Single taxID to search and fetch (e.g., 7227)",
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["protein", "nucleotide", "both"],
        help="Specify sequence type to fetch (protein / nucleotide coding sequence / both)",
    )

    parser.add_argument(
        "--protein-size",
        type=int,
        default=500,
        help="Minimum protein sequence length "
        '(default: 500. Can be bypassed by setting to zero/a negative number)',
    )

    parser.add_argument(
        "--nucleotide-size",
        type=int,
        default=1000,
        help="Minimum nucleotide sequence length"
        '(default: 1000. Can be bypassed by setting to zero/a negative number)',
    )

    parser.add_argument(
        "-e",
        "--email",
        type=str,
        required=True,
        help="Email to use for NCBI API requests (required)",
    )

    parser.add_argument(
        "-k",
        "--api-key",
        type=str,
        required=True,
        help="API key to use for NCBI API requests (required)",
    )

    parser.add_argument(
        "--max-sequences",
        type=int,
        default=None,
        help="Maximum number of sequences to fetch (only works "
        "with -s/--single)",
    )

    parser.add_argument(
        "-b",
        "--genbank",
        action="store_true",
        help="Download GenBank (.gb) files corresponding to fetched sequences",
    )

    return parser

def main():
    print("Starting gene_fetch.py")
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # Validate NCBI credentials first
    print(f"DEBUG: Validating NCBI credentials: email='{args.email}', api_key='{args.api_key}'")
    try:
        Config.validate_credentials(args.email, args.api_key)
        print("DEBUG: Credential validation passed")
    except ValueError as e:
        print(f"ERROR: Credential validation failed: {e}")
        print("DEBUG: About to call sys.exit(1)")
        sys.exit(1)
    print("DEBUG: Continuing after credential validation")

    gene_name = args.gene.lower()
    output_dir = Path(args.out)
    sequence_type = args.type.lower()
    save_genbank = args.genbank  # Get genbank flag

    # Setup output directory and logging
    make_out_dir(output_dir)
    setup_logging(output_dir)

    # Log if GenBank download is enabled
    if save_genbank:
        logger.info(
            "GenBank download mode enabled - will save .gb files in genbank/ subdirectory"
        )

    # Initialize components with required email/api_key
    # No try/catch needed here since we already validated credentials
    config = Config(email=args.email, api_key=args.api_key)

    # Always update thresholds based on user input, regardless of mode
    config.update_thresholds(args.protein_size, args.nucleotide_size)

    # In single-taxid mode, log use of user-specified thresholds
    if args.single:
        logger.info(
            f"Single-taxid mode activated: using protein size threshold {args.protein_size} and nucleotide size threshold {args.nucleotide_size}"
        )

    search_type = config.set_gene_search_term(gene_name)

    if sequence_type not in config.valid_sequence_types:
        print(
            f"Invalid sequence type. Choose from: {', '.join(config.valid_sequence_types)}"
        )
        sys.exit(1)

    logger.info(f"Using {search_type} search terms for {gene_name}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Sequence type: {sequence_type}")

    # Initialize remaining components
    entrez = EntrezHandler(config)
    processor = SequenceProcessor(config, entrez)

    # Check if in single-taxid mode
    if args.single:
        logger.info(f"Single-taxid mode activated for taxid: {args.single}")

        if args.max_sequences:
            logger.info(
                f"Maximum number of sequences to fetch: {args.max_sequences}"
            )
            if sequence_type == "both":
                logger.info(
                    "Note: The max_sequences limit will be applied separately to protein and nucleotide sequences"
                )

        process_single_taxid(
            taxid=args.single,
            gene_name=gene_name,
            sequence_type=sequence_type,
            processor=processor,
            output_dir=output_dir,
            max_sequences=args.max_sequences,
            save_genbank=save_genbank,  # Pass genbank flag
        )
        logger.info("Single taxid processing completed")
        sys.exit(0)
    elif args.max_sequences is not None:
        logger.warning(
            "--max-sequences parameter is ignored when not in single taxid mode"
        )

    # Create output manager
    output_manager = OutputManager(output_dir)

    # Check input file requirements
    if args.input_csv is None and args.input_taxonomy_csv is None:
        logger.error(
            "Error: Either input CSV file (-i/--in) or input taxonomy CSV file (-i2/--in2) must be provided"
        )
        sys.exit(1)

    # Process input samples.csv
    if args.input_csv:
        logger.info(
            f"Starting gene fetch for {gene_name} using taxids from {args.input_csv}"
        )
        process_taxid_csv(
            args.input_csv,
            gene_name,
            sequence_type,
            processor,
            output_manager,
            save_genbank,  # Pass genbank flag
        )

    # Process input samples_taxonomy.csv
    elif args.input_taxonomy_csv:
        logger.info(
            f"Starting gene fetch for {gene_name} using taxonomy from {args.input_taxonomy_csv}"
        )
        process_taxonomy_csv(
            args.input_taxonomy_csv,
            gene_name,
            sequence_type,
            processor,
            output_manager,
            entrez,
            save_genbank,
        )

    logger.info("***********************************************************")
    logger.info("              ? ? ? Gene fetch complete ? ? ?              ")
    logger.info("***********************************************************")
