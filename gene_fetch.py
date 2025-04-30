"""
Gene Fetch - NCBI Sequence Retrieval Tool

This script fetches sequence data from NCBI databases
using sample taxonomic information.

Input:
- NCBI account email address and API key (see:
  https://support.nlm.nih.gov/kbArticle/?pn=KA-05317)
- CSV file containing sample IDs, and taxonomy IDs or taxonomic heirarchy
- Gene name (e.g., 'cox1', '16s', 'rbcl', 'matk')
- Output directory path (will create new directories)
- Sequence type ('protein', 'nucleotide', or 'both')
- Optional: Minimum nucleotide/protein sequence length thresholds
- Optional: 'Single' taxid mode (or default 'batch' taxid mode)
- Optional: Maximum number of sequences to fetch ('single' taxid mode only)
- Optional: Retrieve corresponding genbank entries for each sequence

Output:
- FASTA files containing retrieved sequences (named by ID)
- Log file
- CSV file containing fetched sequence metadata
- CSV file documenting failed retrievals with error types
- 'Single' mode: CSV file(s) containing accession number, length, and sequence

Dependencies:
- Python>=3.9
- Biopython>=1.80
- ratelimit>=2.2.1

Usage:
    python gene_fetch.py -g/--gene <gene_name> -o/--out <output_directory>
                    --type <sequence_type> -i/--in <samples.csv>
                    [-i1/--in2 <sameples_taxonomy.csv] [-s/--single <taxid>]
                    [--protein-size <min_size>] [--nucleotide-size <min_size>]
                    [-e/--email <email>] [-k/--api-key <api_key>]
                    [--max-sequences <max_num>]

    # Required arguments:
    -e/--email           Email to use for NCBI API requests
    -k/--api-key         API key to use for NCBI API requests
    -g/--gene            Gene name to search for (e.g., cox1, 16s, rbcl)
    -o/--out             Directory to save output files
    --type               Sequence type to fetch (protein, nucleotide, or both)
    -i/--in              Input CSV file with taxonomy IDs (required unless
                         using --single)
    -i2/--in2            Input CSV file with taxonomic heirarchy (required
                         unless using --single)

    # Optional arguments:
    -s/--single            Single taxID to fetch all available sequences for
    --protein-size         Minimum protein sequence length (default: 500)
    --nucleotide-size      Minimum nucleotide sequence length (default: 1000)
    --max-sequences        Maximum number of sequences to retrieve (in
                           single mode) (default: all)


Future Development:
- Add optional alignment of retrieved sequences
- Add support for direct GenBank submission format output
- Enhance LRU caching for taxonomy lookups to reduce API calls
- Improve efficiency of record searching and selecting the longest sequence
- Add support for additional genetic markers beyond the currently supported set

Author: D. Parsons
Version: 1.1
License: MIT
"""

# =============================================================================
# Imports and setup
# =============================================================================
import csv
import sys
import time
import argparse
import logging
from Bio import Entrez, SeqIO
from functools import wraps
from ratelimit import limits, sleep_and_retry
from time import sleep
from random import uniform
from urllib.error import HTTPError
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from Bio.Seq import Seq
from http.client import IncompleteRead
import re


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
        help="Path to directory to save output files",
    )

    # Create mutually exclusive group for input files
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument(
        "-i",
        "--in",
        dest="input_csv",
        help="Path to input CSV file containing taxIDs (must have columns "
        '"taxid" and "ID")',
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
        help="Single taxID to fetch all available sequences for",
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["protein", "nucleotide", "both"],
        help="Specify sequence type to fetch",
    )

    parser.add_argument(
        "--protein-size",
        type=int,
        default=500,
        help="Minimum protein sequence length (default: 500)",
    )

    parser.add_argument(
        "--nucleotide-size",
        type=int,
        default=1000,
        help="Minimum nucleotide sequence length (default: 1000)",
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
        help="Maximum number of sequences to fetch in single mode (only works "
        "with -s/--single)",
    )
    
    # Add new argument for GenBank file downloads
    parser.add_argument(
        "-b",
        "--genbank",
        action="store_true",
        help="Download GenBank (.gb) files for fetched sequences",
    )

    return parser


# Ensure output directory exists and create if it doesn't
def make_out_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


# Setup logging with both file and console handlers
def setup_logging(output_dir: Path) -> logging.Logger:
    make_out_dir(output_dir)

    # Clear existing handlers
    logger.handlers.clear()
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    log_path = output_dir / "gene_fetch.log"
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info(f"Logging initialised. Log file: {log_path}")
    return logger


# Initialise logger at module level
logger = logging.getLogger("gene_fetch")


# Log progress at specified intervals
def log_progress(current: int, total: int, interval: int = 10) -> None:
    if current == 0:
        logger.info("")
        logger.info(f"Starting processing: 0/{total} samples processed (0%)")
    elif current == total:
        logger.info(f"Processed: {total}/{total} samples processed (100%)")
    elif current % interval == 0:
        percentage = (current / total) * 100
        logger.info(
            f"=====   Progress: {current}/{total} samples processed "
            f"({percentage:.2f}%)"
        )


# Check NCBI service status before searches using einfo endpoint
def check_ncbi_status():
    try:
        handle = Entrez.einfo()
        # Using underscore to indicate intentionally unused variable
        _ = Entrez.read(handle)
        handle.close()
        # If einfo can be successfully queried then services are up
        return True
    except Exception as e:
        logger.warning(f"NCBI service check failed: {str(e)}")
        return False


# Retry decorator with NCBI status checking and expoential retry delay
def enhanced_retry(
    exceptions: tuple,
    tries: int = 4,
    initial_delay: int = 10,
    backoff: int = 2,
    max_delay: int = 240,
):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            mdelay = initial_delay
            for i in range(tries):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    if i == tries - 1:
                        logger.error(f"Final attempt failed: {str(e)}")
                        return None

                    # Check if it's a server error (500)
                    if isinstance(e, HTTPError) and e.code == 500:
                        # Perform NCBI service check
                        service_status = check_ncbi_status()
                        if not service_status:
                            logger.warning("NCBI services may be experiencing issues")
                            # Use a longer delay for service issues
                            mdelay = min(mdelay * 4, max_delay)

                    # Adjust delay based on error type
                    if isinstance(e, HTTPError):
                        if e.code == 429:  # Too Many Requests
                            mdelay = min(mdelay * 3, max_delay)
                        elif e.code >= 500:  # Server errors
                            mdelay = min(mdelay * 4, max_delay)

                    # Add jitter to avoid thundering herd
                    delay_with_jitter = mdelay + uniform(-0.1 * mdelay, 0.1 * mdelay)
                    logger.warning(
                        f"{str(e)}, Retrying in {delay_with_jitter:.2f} " f"seconds..."
                    )
                    sleep(delay_with_jitter)

                    # Progressive backoff
                    mdelay = min(mdelay * backoff, max_delay)

                    # Additional delay for IncompleteRead errors
                    if isinstance(e, IncompleteRead):
                        logger.warning(
                            "Incomplete read detected, adding additional delay"
                        )
                        sleep(uniform(5, 10))  # Add 5-10 seconds extra delay

            return None

        return wrapper

    return decorator


# Identify the ID column from possible variations
def get_process_id_column(header):
    # Print what's in the header
    logger.info("")
    logger.info(f"CSV header detected: {header}")

    valid_names = [
        "ID",
        "process_id",
        "Process ID",
        "process id",
        "Process id",
        "PROCESS ID",
        "sample",
        "SAMPLE",
        "Sample",
    ]

    # Print repr of each header item to see invisible characters
    for i, col in enumerate(header):
        logger.info(f"Header column {i}: {repr(col)}")

    # Try direct comparison first
    for col in header:
        if col in valid_names:
            logger.info(f"Found matching column: {col}")
            return col

    # Try trimming whitespace (in case of spaces)
    for col in header:
        trimmed = col.strip()
        if trimmed in valid_names:
            logger.info(f"Found matching column after trimming: {trimmed}")
            return col

    # Last resort: try case-insensitive comparison
    for col in header:
        if col.upper() in [name.upper() for name in valid_names]:
            logger.info(f"Found matching column case-insensitive: {col}")
            return col

    logger.error(f"No matching column found in {header}")
    return None


# =============================================================================
# Configuration
# =============================================================================
@dataclass
class Config:
    def __init__(self, email, api_key):
        # Email and API key required
        if not email:
            raise ValueError(
                "Email address required for NCBI API requests. Use -e/--email "
                "to provide your email."
            )
        if not api_key:
            raise ValueError(
                "API key required for NCBI API requests. Use -k/--api-key "
                "to provide your API key."
            )

        self.email = email
        self.api_key = api_key

        # With an API key, 10 requests per second can be made
        self.max_calls_per_second = 10

        # Default batch size for fetching sequences
        self.fetch_batch_size = 200

        # Delay between batches (seconds) - uniform random delay between 1-2 seconds
        self.batch_delay = (1, 2)

        # Set search 'type'
        self.valid_sequence_types = frozenset({"protein", "nucleotide", "both"})

        # Minimum nucleotide and protein lengths for 'batch' mode
        self.protein_length_threshold = 500
        self.nucleotide_length_threshold = 1000

        # Minimum nucleotide and protein lengths for 'single' mode
        self.min_nucleotide_size_single_mode = 200
        self.min_protein_size_single_mode = 100

        self.gene_search_term = ""

        # Define gene type categories
        self._rRNA_genes = {
            "16s": [
                "16S ribosomal RNA[Gene]",
                "16S rRNA[Gene]",
                "rrs[Gene]",
                "rrn16[Gene]",
                "16S ribosomal RNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "16S rRNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "16S[Title] NOT methylase[Title] NOT methyltransferase[Title] "
                "NOT pseudouridylate[Title] NOT synthase[Title]",
                "16S ribosomal RNA[All Fields]",
                "rrs[All Fields]",
                "rrn16[All Fields]",
            ],
            "18s": [
                "18S ribosomal RNA[Gene]",
                "18S rRNA[Gene]",
                "rrn18[Gene]",
                "SSU rRNA[Gene]",
                "18S ribosomal RNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "18S rRNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "18S[Title] NOT methylase[Title] NOT methyltransferase[Title] "
                "NOT pseudouridylate[Title] NOT synthase[Title]",
                "18S ribosomal RNA[rRNA]",
                "SSU rRNA[rRNA]",
                "rrn18[rRNA]",
            ],
            "23s": [
                "23S ribosomal RNA[Gene]",
                "23S rRNA[Gene]",
                "rrl[Gene]",
                "rrn23[Gene]",
                "23S ribosomal RNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "23S rRNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "23S[Title] NOT methylase[Title] NOT methyltransferase[Title] "
                "NOT pseudouridylate[Title] NOT synthase[Title]",
                "23S ribosomal RNA[rRNA]",
                "rrl[rRNA]",
                "rrn23[rRNA]",
            ],
            "28s": [
                "28S ribosomal RNA[Gene]",
                "28S rRNA[Gene]",
                "rrn28[Gene]",
                "LSU rRNA[Gene]",
                "28S ribosomal RNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "28S rRNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "28S[Title] NOT methylase[Title] NOT methyltransferase[Title] "
                "NOT pseudouridylate[Title] NOT synthase[Title]",
                "28S ribosomal RNA[rRNA]",
                "LSU rRNA[rRNA]",
                "rrn28[rRNA]",
            ],
            "12s": [
                "12S ribosomal RNA[Gene]",
                "12S rRNA[Gene]",
                "mt-rrn1[Gene]",
                "mt 12S rRNA[Gene]",
                "12S ribosomal RNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "12S rRNA[Title] NOT methylase[Title] NOT "
                "methyltransferase[Title] NOT pseudouridylate[Title] NOT "
                "synthase[Title]",
                "12S[Title] NOT methylase[Title] NOT methyltransferase[Title] "
                "NOT pseudouridylate[Title] NOT synthase[Title]",
                "12S ribosomal RNA[rRNA]",
                "mt-rrn1[rRNA]",
                "mt 12S rRNA[rRNA]",
            ],
            "its1": [
                "ITS1[Title]",
                "internal transcribed spacer 1[Title]",
                "ITS-1[Title]",
                "ITS 1[Title]",
                "ITS1 ribosomal DNA[Title]",
                "ITS1 rDNA[Title]",
            ],
            "its2": [
                "ITS2[Title]",
                "internal transcribed spacer 2[Title]",
                "ITS-2[Title]",
                "ITS 2[Title]",
                "ITS2 ribosomal DNA[Title]",
                "ITS2 rDNA[Title]",
            ],
            "its": [
                "ITS[Title]",
                "internal transcribed spacer[Title]",
                "ITS region[Title]",
                "ITS1-5.8S-ITS2[Title]",
                "ribosomal ITS[Title]",
                "rDNA ITS[Title]",
            ],
            "trnl": [
                "trnL[Title]",
                "trnL gene[Title]",
                "tRNA-Leu[Title]",
                "tRNA-Leucine[Title]",
                "trnL-trnF[Title]",
                "trnL-F[Title]",
                "chloroplast trnL[Title]",
                "cp trnL[Title]",
                "trnL intron[Title]",
                "trnL-UAA[Title]",
                '"trnL gene"[Gene]',
            ],
        }

        self._protein_coding_genes = {
            "cox1": [
                "cox1[Gene]",
                "COI[Gene]",
                '"cytochrome c oxidase subunit 1"[Protein Name]',
                '"cytochrome oxidase subunit 1"[Protein Name]',
                '"cytochrome c oxidase subunit I"[Protein Name]',
                '"COX1"[Protein Name]',
                '"COXI"[Protein Name]',
            ],
            "cox2": [
                "cox2[Gene]",
                "COII[Gene]",
                '"cytochrome c oxidase subunit 2"[Protein Name]',
                '"cytochrome oxidase subunit 2"[Protein Name]',
                '"cytochrome c oxidase subunit II"[Protein Name]',
                '"COX2"[Protein Name]',
                '"COXII"[Protein Name]',
            ],
            "cox3": [
                "cox3[Gene]",
                "COIII[Gene]",
                '"cytochrome c oxidase subunit 3"[Protein Name]',
                '"cytochrome oxidase subunit 3"[Protein Name]',
                '"cytochrome c oxidase subunit III"[Protein Name]',
                '"COX3"[Protein Name]',
                '"COXIII"[Protein Name]',
            ],
            "cytb": [
                "cytb[Gene]",
                "cob[Gene]",
                '"cytochrome b"[Protein Name]',
                '"cytochrome b"[Gene]',
                '"CYTB"[Protein Name]',
            ],
            "nd1": [
                "nd1[Gene]",
                "NAD1[Gene]",
                '"NADH dehydrogenase subunit 1"[Protein Name]',
                '"ND1"[Protein Name]',
            ],
            "nd2": [
                "nd2[Gene]",
                "NAD2[Gene]",
                '"NADH dehydrogenase subunit 2"[Protein Name]',
                '"ND2"[Protein Name]',
            ],
            "rbcl": [
                "rbcL[Gene]",
                "RBCL[Gene]",
                '"ribulose-1,5-bisphosphate carboxylase/oxygenase large '
                'subunit"[Protein Name]',
                '"ribulose 1,5-bisphosphate carboxylase/oxygenase large '
                'subunit"[Protein Name]',
                '"ribulose bisphosphate carboxylase large chain"[Protein Name]',
                '"RuBisCO large subunit"[Protein Name]',
                '"ribulose-1,5-bisphosphate carboxylase/oxygenase small '
                'subunit"[Protein Name]',
                '"ribulose 1,5-bisphosphate carboxylase/oxygenase small '
                'subunit"[Protein Name]',
                '"ribulose-1,5-bisphosphate carboxylase/oxygenase small '
                'chain"[Protein Name]',
                '"RuBisCO small subunit"[Protein Name]',
                '"Ribulose bisphosphate carboxylase/oxygenase ' 'activase"[Gene]',
                '"rbcL gene"[Gene]',
                '"RBCL gene"[Gene]',
                "rbcL[Title]",
                "RBCL[Title]",
            ],
            "matk": [
                "matK[Gene]",
                "MATK[Gene]",
                '"maturase K"[Protein Name]',
                '"maturase K"[Gene]',
                '"maturase-K"[Protein Name]',
                '"maturase-K"[Gene]',
                '"Maturase K"[Protein Name]',
                '"Maturase K"[Gene]',
                '"matK gene"[Gene]',
                '"MATK gene"[Gene]',
                '"trnK-matK"[Gene]',
                '"maturase type II intron splicing factor"[Protein Name]',
                '"chloroplast group II intron splicing factor maturase '
                'K"[Protein Name]',
                '"type II intron maturase K"[Protein Name]',
                '"tRNA-lysine maturase K"[Protein Name]',
            ],
            "psba": [
                "psbA[Gene]",
                "PSBA[Gene]",
                '"photosystem II protein D1"[Protein Name]',
                '"photosystem II protein D1"[Gene]',
                '"PSII D1 protein"[Protein Name]',
                '"PSII D1 protein"[Gene]',
                '"photosystem II reaction center protein D1"[Protein Name]',
                '"photosystem II reaction center protein D1"[Gene]',
                '"photosystem Q(B) protein"[Protein Name]',
                '"32 kDa thylakoid membrane protein"[Protein Name]',
                '"psbA gene"[Gene]',
                '"PSBA gene"[Gene]',
                "psbA[Title]",
                "PSBA[Title]",
                '"trnH-psbA"[Title]',
                '"psbA-trnH"[Title]',
            ],
        }

    # Update sequence length thresholds
    def update_thresholds(self, protein_size: int, nucleotide_size: int):
        self.protein_length_threshold = protein_size
        self.nucleotide_length_threshold = nucleotide_size

    # Set search term based on gene name and type
    def set_gene_search_term(self, gene_name: str):
        gene_name = gene_name.lower()

        # Check if rRNA
        if gene_name in self._rRNA_genes:
            self.gene_search_term = "(" + " OR ".join(self._rRNA_genes[gene_name]) + ")"
            search_type = "rRNA"

        # Check if protein-coding
        elif gene_name in self._protein_coding_genes:
            self.gene_search_term = (
                "(" + " OR ".join(self._protein_coding_genes[gene_name]) + ")"
            )
            search_type = "protein-coding"

        else:
            # Generic search term for non-listed genes
            self.gene_search_term = (
                f"({gene_name}[Title] OR {gene_name}[Gene] OR "
                f'"{gene_name}"[Text Word])'
            )
            search_type = "generic"

        return search_type  # Return the type for logging purposes


# =============================================================================
# Entrez API interactions and error handling
# =============================================================================
class EntrezHandler:
    def __init__(self, config: Config):
        self.config = config
        Entrez.email = config.email
        Entrez.api_key = config.api_key
        self.consecutive_errors = 0
        self.last_service_check = 0
        self.service_check_interval = 60  # Seconds between service checks

        # Creating cache for fetched NCBI taxonomy
        self.taxonomy_cache = {}  # taxid -> (complete_lineage, rank_info, current_rank, taxid_info)

    # Determine if service status check is needed
    def should_check_service_status(self) -> bool:
        current_time = time.time()
        if current_time - self.last_service_check > self.service_check_interval:
            return True
        return False

    # Handle request errors and track consecutive failures
    def handle_request_error(self, error: Exception) -> None:
        self.consecutive_errors += 1
        if self.consecutive_errors >= 3 and self.should_check_service_status():
            service_status = check_ncbi_status()
            self.last_service_check = time.time()
            if not service_status:
                logger.warning(
                    "Multiple consecutive errors - NCBI service appears it might "
                    "be down?"
                )
                sleep(uniform(30, 60))  # Longer delay when service is down

    def handle_request_success(self) -> None:
        """Reset error counter on successful request."""
        self.consecutive_errors = 0

    @enhanced_retry((HTTPError, RuntimeError, IOError, IncompleteRead))
    @sleep_and_retry
    @limits(calls=10, period=1.1)
    # Entrez efetch
    def fetch(self, **kwargs) -> Optional[Any]:
        try:
            result = Entrez.efetch(**kwargs)
            self.handle_request_success()
            return result
        except Exception as e:
            self.handle_request_error(e)
            raise

    # Entrez esearch in batches
    def search(self, **kwargs) -> Optional[Dict]:
        @enhanced_retry((HTTPError, RuntimeError, IOError, IncompleteRead))
        @sleep_and_retry
        @limits(calls=10, period=1.1)
        def _do_search(**kwargs):
            try:
                handle = Entrez.esearch(**kwargs)
                result = Entrez.read(handle)
                handle.close()
                self.handle_request_success()
                return result
            except Exception as e:
                self.handle_request_error(e)
                raise

        # First search to get total count
        initial_result = _do_search(**kwargs)
        if not initial_result:
            return None

        total_count = int(initial_result["Count"])
        all_ids = initial_result["IdList"]

        # If there are more results, fetch them
        batch_size = 100  # Increased from default 20
        if total_count > len(all_ids):
            for start in range(len(all_ids), total_count, batch_size):
                kwargs["retstart"] = start
                kwargs["retmax"] = batch_size
                result = _do_search(**kwargs)
                if result and result.get("IdList"):
                    all_ids.extend(result["IdList"])

                # Add delay between batches
                sleep(uniform(1, 2))

        # Return modified result with all IDs
        initial_result["IdList"] = all_ids
        return initial_result
        
        # Fetch taxonomy information for a given taxid
    @enhanced_retry(
        (HTTPError, RuntimeError, IOError, IncompleteRead),
        tries=5,
        initial_delay=15,
    )
    def fetch_taxonomy(
            self, taxid: str
        ) -> Tuple[List[str], Dict[str, str], str, Dict[str, str]]:
            # Check cache first
            if taxid in self.taxonomy_cache:
                logger.info(f"Using cached taxonomy for taxID: {taxid}")  
                return self.taxonomy_cache[taxid]
            
            logger.info(f"Retrieving taxonomy details from NCBI for taxID: {taxid}")
            
            # First verify taxid format
            taxid = taxid.strip()
            if not taxid.isdigit():
                logger.error(f"Invalid taxID format: {taxid} (must be numerical)")
                return [], {}, "", {}

            # Add initial delay to help avoid rate limiting
            sleep(uniform(0.5, 1.0))

            try:
                # Set up parameters for the request
                params = {
                    "db": "taxonomy",
                    "id": taxid,
                    "email": self.config.email,
                    "api_key": self.config.api_key,
                    "tool": "gene_fetch",
                }

                max_retries = 3
                for attempt in range(max_retries):
                    try:
                        handle = self.fetch(**params)
                        records = Entrez.read(handle)
                        handle.close()
                        break  # If successful, exit retry loop
                    # Deal with spurious HTTP 400 errors
                    except HTTPError as e:
                        if e.code == 400:
                            if attempt < max_retries - 1:  # If not the last attempt
                                delay = (attempt + 1) * 2  # Progressive delay
                                logger.warning(
                                    f"HTTP 400 error for taxID {taxid}, attempt "
                                    f"{attempt + 1}/{max_retries}. Retrying in "
                                    f"{delay} seconds..."
                                )
                                sleep(delay)
                                continue
                            else:
                                logger.error(
                                    f"Failed to fetch taxonomy after {max_retries} "
                                    f"attempts for taxID {taxid}"
                                )
                                return [], {}, "", {}
                        else:
                            raise  # Re-raise other HTTP errors
                else:  # If all retries exhausted
                    return [], {}, "", {}

                if not records:
                    logger.error(f"No taxonomy records found for taxID {taxid}")
                    return [], {}, "", {}

                # Get the first record
                record = records[0]

                # Initialise rank information dictionaries
                rank_info = {}
                taxid_info = {}

                lineage_nodes = record.get("LineageEx", [])
                lineage = []

                # Process lineage nodes
                for node in lineage_nodes:
                    name = node.get("ScientificName", "")
                    rank = node.get("Rank", "no rank")
                    node_taxid = str(node.get("TaxId", ""))

                    # Add to lineage
                    lineage.append(name)

                    # Add rank and taxid info if valid
                    if name and rank != "no rank":
                        rank_info[name] = rank
                        taxid_info[name] = node_taxid

                # Get current taxon information
                current_name = record.get("ScientificName", "")
                current_rank = record.get("Rank", "no rank")

                # Add current taxon to complete lineage
                complete_lineage = lineage + [current_name]

                if current_rank != "no rank":
                    rank_info[current_name] = current_rank
                    taxid_info[current_name] = taxid

                logger.info(f"Successfully retrieved taxonomy information from NCBI for taxID: {taxid}")
                logger.debug(f"Lineage: {complete_lineage}")
                logger.debug(f"Rank info: {rank_info}")

                # Store taxonomy in cache
                self.taxonomy_cache[taxid] = (complete_lineage, rank_info, current_rank, taxid_info)

                return complete_lineage, rank_info, current_rank, taxid_info

            except IncompleteRead as e:
                logger.warning(f"IncompleteRead error for taxID {taxid}: {e}")
                # Re-raise to allow the retry decorator to handle it
                raise
            except Exception as e:
                if isinstance(e, HTTPError) and e.code == 400:
                    logger.error(f"HTTP 400 error for taxID {taxid}, skipping for now")
                else:
                    logger.error(f"Error fetching taxonomy for taxID {taxid}: {e}")
                    logger.error("Full error details:", exc_info=True)
                return [], {}, "", {}

    # Fetch NCBI taxID from input hierarchical taxonomic information
    # Conducts progressive search from most specific level -> phylum
    # Returns taxid for validation
    def fetch_taxid_from_taxonomy(
        self, phylum, class_name, order, family, genus, species
    ):
        # Log fetched sample taxonomy
        logger.info(
            f"Attempting to fetch taxid for Species: {species}, Genus: {genus},"
            f"Family: {family}, Order: {order}, Class: {class_name}, Phylum: {phylum})"
        )
        
        # Store the provided taxonomy for validation
        provided_taxonomy = {
            "phylum": phylum.strip() if phylum else "",
            "class": class_name.strip() if class_name else "",
            "order": order.strip() if order else "",
            "family": family.strip() if family else "",
            "genus": genus.strip() if genus else "",
            "species": species.strip() if species else ""
        }
        
        # Try species level first, and validate
        if genus and species:
            # Check if species already contains genus name to avoid duplication
            if species.startswith(genus):
                # Use species as-is since it already contains the full name
                full_species = species
            else:
                # Combine genus and species for the full species name
                full_species = f"{genus} {species}"
                
            search_term = f"{full_species}[Scientific Name]"
            logger.info(f"Searching for species: {search_term}")
            try:
                result = self.search(db="taxonomy", term=search_term)
                if result and result.get("IdList") and len(result["IdList"]) > 0:
                    taxids = result["IdList"]
                    
                    # Check for multiple matches
                    if len(taxids) > 1:
                        logger.warning(f"Multiple taxids found for species {full_species}: {taxids}")
                        
                        # Validate fetched taxids against the input taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                            if is_valid:
                                valid_taxids.append((taxid, lineage))
                        
                        if not valid_taxids:
                            logger.warning(f"None of the taxids for {full_species} match the provided higher taxonomy")
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(f"Found single valid taxid ({taxid}) for species {full_species}")
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(valid_taxids, key=lambda x: x[1]['match_score'])[0]
                            logger.warning(f"Multiple valid taxids for {full_species}, using best match: {best_taxid}")
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                        if is_valid:
                            logger.info(f"Found taxid ({taxid}) for species {full_species}")
                            return taxid
                        else:
                            logger.warning(f"Taxid {taxid} for {full_species} does not match the provided higher taxonomy")
                            # Continue to next level
            except Exception as e:
                logger.error(f"Error searching for species taxid: {e}")

        # Try genus level next, with validation
        if genus:
            search_term = f"{genus}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching for genus: {search_term}")
            try:
                result = self.search(db="taxonomy", term=search_term)
                if result and result.get("IdList") and len(result["IdList"]) > 0:
                    taxids = result["IdList"]
                    
                    # Check for multiple matches
                    if len(taxids) > 1:
                        logger.warning(f"Multiple taxids found for genus {genus}: {taxids}")
                        
                        # Validate each taxid against the provided taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                            if is_valid:
                                valid_taxids.append((taxid, lineage))
                        
                        if not valid_taxids:
                            logger.warning(f"None of the taxids for genus {genus} match the provided higher taxonomy")
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(f"Found single valid taxid ({taxid}) for genus {genus}")
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(valid_taxids, key=lambda x: x[1]['match_score'])[0]
                            logger.warning(f"Multiple valid taxids for genus {genus}, using best match: {best_taxid}")
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                        if is_valid:
                            logger.info(f"Confirmed taxid {taxid} for {genus} is valid given input taxonomy")
                            return taxid
                        else:
                            logger.warning(f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy")
                            # Continue to next level
            except Exception as e:
                logger.error(f"Error searching for genus taxid: {e}")

        # Try family level with same validation approach
        if family:
            search_term = f"{family}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching for family: {search_term}")
            try:
                result = self.search(db="taxonomy", term=search_term)
                if result and result.get("IdList") and len(result["IdList"]) > 0:
                    taxids = result["IdList"]
                    
                    # Check for multiple matches
                    if len(taxids) > 1:
                        logger.warning(f"Multiple taxids found for family {family}: {taxids}")
                        
                        # Validate each taxid against the provided taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                            if is_valid:
                                valid_taxids.append((taxid, lineage))
                        
                        if not valid_taxids:
                            logger.warning(f"None of the taxids for family {family} match the provided higher taxonomy")
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(f"Found single valid taxid ({taxid}) for family {family}")
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(valid_taxids, key=lambda x: x[1]['match_score'])[0]
                            logger.warning(f"Multiple valid taxids for family {family}, using best match: {best_taxid}")
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                        if is_valid:
                            logger.info(f"Confirmed taxid {taxid} for {family} is valid given input taxonomy")
                            return taxid
                        else:
                            logger.warning(f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy")
                            # Continue to next level
            except Exception as e:
                logger.error(f"Error searching for family taxid: {e}")

        # Try order level with same validation approach
        if order:
            search_term = f"{order}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching for order: {search_term}")
            try:
                result = self.search(db="taxonomy", term=search_term)
                if result and result.get("IdList") and len(result["IdList"]) > 0:
                    taxids = result["IdList"]
                    
                    # Check for multiple matches
                    if len(taxids) > 1:
                        logger.warning(f"Multiple taxids found for order {order}: {taxids}")
                        
                        # Validate each taxid against the provided taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                            if is_valid:
                                valid_taxids.append((taxid, lineage))
                        
                        if not valid_taxids:
                            logger.warning(f"None of the taxids for order {order} match the provided higher taxonomy")
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(f"Found single valid taxid ({taxid}) for order {order}")
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(valid_taxids, key=lambda x: x[1]['match_score'])[0]
                            logger.warning(f"Multiple valid taxids for order {order}, using best match: {best_taxid}")
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                        if is_valid:
                            logger.info(f"Confirmed taxid {taxid} for {order} is valid given input taxonomy")
                            return taxid
                        else:
                            logger.warning(f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy")
                            # Continue to next level
            except Exception as e:
                logger.error(f"Error searching for order taxid: {e}")

        # Try class level with same validation approach
        if class_name:
            search_term = f"{class_name}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching for class: {search_term}")
            try:
                result = self.search(db="taxonomy", term=search_term)
                if result and result.get("IdList") and len(result["IdList"]) > 0:
                    taxids = result["IdList"]
                    
                    # Check for multiple matches
                    if len(taxids) > 1:
                        logger.warning(f"Multiple taxids found for class {class_name}: {taxids}")
                        
                        # At class level, we validate against phylum only
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                            if is_valid:
                                valid_taxids.append((taxid, lineage))
                        
                        if not valid_taxids:
                            logger.warning(f"None of the taxids for class {class_name} match the provided higher taxonomy")
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(f"Found single valid taxid ({taxid}) for class {class_name}")
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(valid_taxids, key=lambda x: x[1]['match_score'])[0]
                            logger.warning(f"Multiple valid taxids for class {class_name}, using best match: {best_taxid}")
                            return best_taxid
                    else:
                        # Single match - still validate against phylum
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(taxid, provided_taxonomy)
                        if is_valid:
                            logger.info(f"Confirmed taxid {taxid} for {class_name} is valid given input taxonomy")
                            return taxid
                        else:
                            logger.warning(f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy")
                            # Continue to next level
            except Exception as e:
                logger.error(f"Error searching for class taxid: {e}")

        # Try phylum level with minimal validation
        if phylum:
            search_term = f"{phylum}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching for phylum: {search_term}")
            try:
                result = self.search(db="taxonomy", term=search_term)
                if result and result.get("IdList") and len(result["IdList"]) > 0:
                    taxids = result["IdList"]
                    
                    # At phylum level we have less to validate against. Still check for multiple matches
                    if len(taxids) > 1:
                        logger.warning(f"Multiple taxids found for phylum {phylum}: {taxids}")
                        # Take the first one since we're already at phylum level
                        taxid = taxids[0]
                    else:
                        taxid = taxids[0]
                    
                    logger.info(f"Confirmed taxid {taxid} for {phylum} is valid given input taxonomy")
                    return taxid
            except Exception as e:
                logger.error(f"Error searching for phylum taxid: {e}")

        logger.warning(f"Could not find any valid taxid for {genus} {species}")
        return None

    # Calls fetch_taxonomy to access fetched NCBI taxonomy information
    # Validate the fetched taxonomy is consistent with the input taxonomy
    # Calcualte taxonomy match scores
    def validate_taxonomy_consistency(self, taxid, provided_taxonomy):
        logger.info(f"Validating taxonomy:")
        logger.info(f"Comparing provided taxonomy against NCBI taxonomy for taxid {taxid}")
        
        try:
            # Use the fetch_taxonomy method directly from EntrezHandler to get lineage information
            complete_lineage, rank_info, current_rank, taxid_info = self.fetch_taxonomy(taxid)
            
            if not complete_lineage:
                logger.warning(f"No taxonomy data available for comparison (taxid: {taxid})")
                return False, {'match_score': 0}
            
            # Build dictionary of the fetched taxonomy organized by rank
            fetched_taxonomy = {}
            for taxon in complete_lineage:
                if taxon in rank_info:
                    rank = rank_info[taxon]
                    fetched_taxonomy[rank.lower()] = taxon
            
            # Add the current taxon if it has a rank
            if current_rank != "no rank":
                current_taxon = complete_lineage[-1] if complete_lineage else ""
                fetched_taxonomy[current_rank.lower()] = current_taxon

            logger.debug(f"Taxonomy comparison details - NCBI: {fetched_taxonomy}, Provided: {provided_taxonomy}")
                
            # Check consistency between provided and fetched taxonomy
            # Start with higher taxonomic ranks which are less likely to have homonyms
            match_score = 0
            
            # Store details of what matched and what didn't
            match_details = {}
            
            # Check phylum
            if provided_taxonomy["phylum"] and "phylum" in fetched_taxonomy:
                phylum_match = provided_taxonomy["phylum"].lower() == fetched_taxonomy["phylum"].lower()
                match_details["phylum"] = "match" if phylum_match else "mismatch"
                if phylum_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(f"Phylum mismatch for taxid {taxid}: "
                                 f"provided '{provided_taxonomy['phylum']}' vs "
                                 f"fetched '{fetched_taxonomy['phylum']}'")
            
            # Check class
            if provided_taxonomy["class"] and "class" in fetched_taxonomy:
                class_match = provided_taxonomy["class"].lower() == fetched_taxonomy["class"].lower()
                match_details["class"] = "match" if class_match else "mismatch"
                if class_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(f"Class mismatch for taxid {taxid}: "
                                 f"provided '{provided_taxonomy['class']}' vs "
                                 f"fetched '{fetched_taxonomy['class']}'")
            
            # Check order
            if provided_taxonomy["order"] and "order" in fetched_taxonomy:
                order_match = provided_taxonomy["order"].lower() == fetched_taxonomy["order"].lower()
                match_details["order"] = "match" if order_match else "mismatch"
                if order_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(f"Order mismatch for taxid {taxid}: "
                                 f"provided '{provided_taxonomy['order']}' vs "
                                 f"fetched '{fetched_taxonomy['order']}'")
            
            # Check family
            if provided_taxonomy["family"] and "family" in fetched_taxonomy:
                family_match = provided_taxonomy["family"].lower() == fetched_taxonomy["family"].lower()
                match_details["family"] = "match" if family_match else "mismatch"
                if family_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(f"Family mismatch for taxid {taxid}: "
                                 f"provided '{provided_taxonomy['family']}' vs "
                                 f"fetched '{fetched_taxonomy['family']}'")
            
            # Check genus
            if provided_taxonomy["genus"] and "genus" in fetched_taxonomy:
                genus_match = provided_taxonomy["genus"].lower() == fetched_taxonomy["genus"].lower()
                match_details["genus"] = "match" if genus_match else "mismatch"
                if genus_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(f"Genus mismatch for taxid {taxid}: "
                                 f"provided '{provided_taxonomy['genus']}' vs "
                                 f"fetched '{fetched_taxonomy['genus']}'")
            
            # Determine if this taxid is valid based on matches
            # We consider it valid if:
            # 1. No conflicts in higher taxonomy (phylum, class)
            # 2. Or if no higher taxonomy was provided to check against
            
            # Check for higher taxonomy mismatches that would invalidate the match
            higher_taxonomy_conflict = False
            
            if (provided_taxonomy["phylum"] and 
                "phylum" in fetched_taxonomy and 
                match_details.get("phylum") == "mismatch"):
                higher_taxonomy_conflict = True
            
            if (provided_taxonomy["class"] and 
                "class" in fetched_taxonomy and 
                match_details.get("class") == "mismatch"):
                higher_taxonomy_conflict = True
            
            is_valid = not higher_taxonomy_conflict
            
            # Log validation result
            if is_valid:
                logger.info(f"Taxonomy validation: taxid {taxid} passed validation with match score: {match_score}")
            else:
                logger.warning(f"Taxonomy validation: taxid {taxid} failed validation due to higher taxonomy conflicts")
            
            return is_valid, {
                'match_score': match_score,
                'details': match_details,
                'lineage': complete_lineage,
                'fetched_taxonomy': fetched_taxonomy
            }
            
        except Exception as e:
            logger.error(f"Error validating taxonomy for taxid {taxid}: {e}")
            logger.error("Full error details:", exc_info=True)
            return False, {'match_score': 0}


# =============================================================================
# Processing and validation of fetched entries and sequences
# =============================================================================
class SequenceProcessor:
    def __init__(self, config: Config, entrez: EntrezHandler):
        self.config = config
        self.entrez = entrez

    # Handle 'contig' lines in feature table
    def parse_contig_line(self, contig_line: str) -> Optional[Tuple[str, int, int]]:
        try:
            # Remove 'join(' and ')' if present
            cleaned = contig_line.strip().replace("join(", "").replace(")", "")

            # Parse the contig reference. Format is typically: WVEN01000006.2:1..16118
            if ":" in cleaned:
                contig_id, coords = cleaned.split(":")
                if ".." in coords:
                    start, end = coords.split("..")
                    return contig_id, int(start), int(end)
        except Exception as e:
            logger.error(f"Error parsing CONTIG line '{contig_line}': {e}")
        return None

    # Handle WGS records that might not contain sequence data, fetching associated sequence if listed in 'contig' line
    def fetch_wgs_sequence(self, record: SeqRecord) -> Optional[SeqRecord]:
        try:
            # Check for CONTIG line in annotations
            contig_line = record.annotations.get("contig", None)
            if not contig_line:
                logger.warning(f"No CONTIG line found in WGS record {record.id}")
                return None

            # Parse the CONTIG line
            contig_info = self.parse_contig_line(contig_line)
            if not contig_info:
                logger.error(f"Could not parse CONTIG line: {contig_line}")
                return None

            contig_id, start, end = contig_info
            logger.info(
                f"Found WGS contig reference: {contig_id} positions {start}..{end}"
            )

            # Fetch the actual contig sequence
            try:
                handle = self.entrez.fetch(
                    db="nucleotide",
                    id=contig_id,
                    rettype="fasta",
                    retmode="text",
                )
                contig_record = next(SeqIO.parse(handle, "fasta"))
                handle.close()

                if not contig_record or not contig_record.seq:
                    logger.error(f"Failed to fetch sequence for contig {contig_id}")
                    return None

                # Extract the relevant portion
                start_idx = start - 1  # Convert to 0-based indexing
                sequence = contig_record.seq[start_idx:end]

                if not sequence:
                    logger.error(
                        f"Extracted sequence is empty for positions {start}..{end}"
                    )
                    return None

                # Create new record with the sequence
                new_record = record[:]
                new_record.seq = sequence
                logger.info(f"Successfully extracted {len(sequence)}bp from WGS contig")

                return new_record

            except Exception as e:
                logger.error(f"Error fetching WGS contig {contig_id}: {e}")
                return None

        except Exception as e:
            logger.error(f"Error processing WGS record: {e}")
            return None

    # Wrapper function to fetch nucleotide sequences, including WGS records if they have an available sequence
    def fetch_nucleotide_record(self, record_id: str) -> Optional[SeqRecord]:
        try:
            # Fetch the genbank record
            handle = self.entrez.fetch(
                db="nucleotide", id=record_id, rettype="gb", retmode="text"
            )
            record = next(SeqIO.parse(handle, "genbank"))
            handle.close()

            # Check if it's a WGS record
            is_wgs = False
            if hasattr(record, "annotations"):
                keywords = record.annotations.get("keywords", [])
                if "WGS" in keywords:
                    is_wgs = True
                    logger.info(f"WGS record detected for {record_id}")

            # For WGS records, check if there's a complete sequence available directly
            if is_wgs:
                if record.seq is not None and len(record.seq) > 0:
                    try:
                        # Verify sequence can be accessed
                        seq_str = str(record.seq)
                        if seq_str and not seq_str.startswith("?"):
                            logger.info(
                                f"WGS record {record_id} has a complete sequence "
                                f"of length {len(record.seq)}"
                            )
                            return record
                        else:
                            logger.info(
                                f"WGS record {record_id} has a placeholder/"
                                f"incomplete sequence"
                            )
                    except Exception as e:
                        logger.warning(
                            f"Unable to access sequence content directly from "
                            f"WGS record {record_id}: {e}"
                        )

                # If there is no sequence, check for 'contig' line
                if "contig" in record.annotations:
                    logger.info(
                        f"WGS record {record_id} has a CONTIG line, attempting "
                        f"to fetch underlying sequence"
                    )
                    wgs_record = self.fetch_wgs_sequence(record)
                    if wgs_record:
                        return wgs_record
                    else:
                        logger.warning(
                            f"Failed to fetch sequence from WGS CONTIG for "
                            f"{record_id}"
                        )

                # If no sequence, log and return None
                logger.info(
                    f"WGS record {record_id} does not have a usable sequence - "
                    f"skipping"
                )
                return None

            # Skip unverified sequences
            if (
                "unverified" in record.description.lower()
                or "UNVERIFIED" in record.description
            ):
                logger.info(f"Unverified sequence detected for {record_id} - skipping")
                return None

            # For non-WGS and verified records, verify sequence content
            if record.seq is not None and len(record.seq) > 0:
                try:
                    # Verify sequence can be accessed
                    _ = str(record.seq)
                    return record
                except Exception as e:
                    logger.error(f"Undefined sequence content for {record_id}: {e}")
                    return None

            return None

        except Exception as e:
            logger.error(f"Error fetching nucleotide sequence for {record_id}: {e}")
            return None

    # Extract CDS region from a sequence record with fallbacks for sequence types and name variations
    def extract_nucleotide(
        self, record: SeqRecord, gene_name: str, single_mode: bool = False
    ) -> Optional[SeqRecord]:
        # Prepare gene name variations for matching
        gene_variations = set()
        pattern_variations = []

        # Get minimum size threshold for single mode
        min_size = self.config.min_nucleotide_size_single_mode if single_mode else 100

        # Initialize base_gene at the beginning to ensure it's always defined
        base_gene = gene_name.lower()

        if base_gene in self.config._protein_coding_genes:
            # Get the gene variations from config
            variations = self.config._protein_coding_genes[base_gene]
            gene_variations = {v.split("[")[0].strip('"').lower() for v in variations}

            # Add common pattern variations for different writing styles
            if base_gene == "rbcl":
                pattern_variations = [
                    "rbcl",
                    "rbc-l",
                    "rbc l",
                    "rubisco",
                    "ribulose-1,5-bisphosphate",
                    "ribulose bisphosphate",
                ]
            elif base_gene.startswith("cox"):
                pattern_variations = ["cytochrome c oxidase"]
            elif base_gene == "cytb":
                pattern_variations = [
                    "cytb",
                    "cyt b",
                    "cyt-b",
                    "cytochrome b",
                    "cytochrome-b",
                ]
            elif base_gene == "matk":
                pattern_variations = [
                    "matk",
                    "mat-k",
                    "mat k",
                    "maturase k",
                    "maturase-k",
                    "maturase",
                    "chloroplast maturase k",
                    "trnk-matk",
                ]
            elif base_gene == "nd1":
                pattern_variations = [
                    "nd1",
                    "nd-1",
                    "nd 1",
                    "nadh1",
                    "nadh-1",
                    "nadh 1",
                    "nadh dehydrogenase 1",
                    "nadh dehydrogenase subunit 1",
                    "nadh-dehydrogenase 1",
                    "nad1",
                    "nad-1",
                ]
            elif base_gene == "nd2":
                pattern_variations = [
                    "nd2",
                    "nd-2",
                    "nd 2",
                    "nadh2",
                    "nadh-2",
                    "nadh 2",
                    "nadh dehydrogenase 2",
                    "nadh dehydrogenase subunit 2",
                    "nadh-dehydrogenase 2",
                    "nad2",
                    "nad-2",
                ]
            elif base_gene == "nd4":
                pattern_variations = [
                    "nd4",
                    "nd-4",
                    "nd 4",
                    "nadh4",
                    "nadh-4",
                    "nadh 4",
                    "nadh dehydrogenase 4",
                    "nadh dehydrogenase subunit 4",
                    "nadh-dehydrogenase 4",
                    "nad4",
                    "nad-4",
                ]
            elif base_gene == "nd5":
                pattern_variations = [
                    "nd5",
                    "nd-5",
                    "nd 5",
                    "nadh5",
                    "nadh-5",
                    "nadh 5",
                    "nadh dehydrogenase 5",
                    "nadh dehydrogenase subunit 5",
                    "nadh-dehydrogenase 5",
                    "nad5",
                    "nad-5",
                ]
            elif base_gene == "atp6":
                pattern_variations = [
                    "atp6",
                    "atp-6",
                    "atp 6",
                    "atpase6",
                    "atpase-6",
                    "atpase 6",
                    "atp synthase 6",
                    "atp synthase subunit 6",
                    "atp synthase f0 subunit 6",
                    "atpase subunit 6",
                    "atpase subunit a",
                ]
            elif base_gene == "atp8":
                pattern_variations = [
                    "atp8",
                    "atp-8",
                    "atp 8",
                    "atpase8",
                    "atpase-8",
                    "atpase 8",
                    "atp synthase 8",
                    "atp synthase subunit 8",
                    "atp synthase f0 subunit 8",
                    "atpase subunit 8",
                ]
        elif base_gene == "16s" or base_gene == "16s rrna" or base_gene == "rrn16":
            pattern_variations = [
                "16s",
                "16s rrna",
                "16s ribosomal rna",
                "16s ribosomal",
                "16 s rrna",
                "16 s",
                "rrn16",
                "rrn 16",
            ]
        elif base_gene == "18s" or base_gene == "18s rrna" or base_gene == "rrn18":
            pattern_variations = [
                "18s",
                "18s rrna",
                "18s ribosomal rna",
                "18s ribosomal",
                "18 s rrna",
                "18 s",
                "rrn18",
                "rrn 18",
                "small subunit ribosomal rna",
                "ssu rrna",
                "ssu",
            ]
        elif base_gene == "28s" or base_gene == "28s rrna" or base_gene == "rrn28":
            pattern_variations = [
                "28s",
                "28s rrna",
                "28s ribosomal rna",
                "28s ribosomal",
                "28 s rrna",
                "28 s",
                "rrn28",
                "rrn 28",
                "large subunit ribosomal rna",
                "lsu rrna",
                "lsu",
            ]
        elif base_gene == "its" or base_gene == "its1" or base_gene == "its2":
            pattern_variations = [
                "internal transcribed spacer",
                "its region",
                "its1-5.8s-its2",
                "its 1",
                "its 2",
                "its1",
                "its2",
                "its 1-5.8s-its 2",
                "ribosomal its",
                "rrna its",
            ]
        elif (
            base_gene == "trnh-psba" or base_gene == "psba-trnh" or base_gene == "psba"
        ):
            pattern_variations = [
                "trnh-psba",
                "psba-trnh",
                "trnh psba",
                "psba trnh",
                "trnh-psba spacer",
                "psba-trnh spacer",
                "trnh-psba intergenic spacer",
                "trnh psba intergenic",
                "psba trnh intergenic",
            ]
        else:
            logger.warning(f"No defined variations for gene {gene_name}")
            # For any other gene (including tetraspanin), add reasonable variations
            gene_variations = {base_gene}
            pattern_variations = [
                base_gene,
                f"{base_gene} gene",
                f"{base_gene} protein",
                f"{base_gene}-like",
            ]

        # If pattern variations are present, add them to regular variations
        if pattern_variations:
            gene_variations.update(pattern_variations)

        logger.info(f"Using gene variations for matching: {gene_variations}")

        # STEP 1: Try to find a CDS feature with EXACT match to target gene
        found_cds = None

        # First, look for exact matches to our target gene name
        target_gene = gene_name.lower()
        for feature in record.features:
            if feature.type != "CDS":
                continue

            qualifiers = []
            for field in ["gene", "product", "note"]:
                qualifiers.extend(feature.qualifiers.get(field, []))
            logger.info(f"Found CDS qualifiers: {qualifiers}")

            # Check for exact match first
            for qualifier in qualifiers:
                qualifier_lower = qualifier.lower()

                # Exact match to target gene (e.g., 'cox1', 'coi')
                if (
                    target_gene in qualifier_lower.split()
                    or f"{target_gene}" == qualifier_lower
                ):
                    logger.info(
                        f"Found exact match for {target_gene} in qualifier: "
                        f"{qualifier}"
                    )
                    found_cds = feature
                    break

                # For cox genes, check for the specific number match
                if target_gene.startswith("cox"):
                    if "cox" in qualifier_lower and target_gene[-1] in qualifier_lower:
                        if (
                            f"cox{target_gene[-1]}" in qualifier_lower
                            or f"cox {target_gene[-1]}" in qualifier_lower
                        ):
                            logger.info(
                                f"Found cox{target_gene[-1]} match in qualifier: "
                                f"{qualifier}"
                            )
                            found_cds = feature
                            break

                # For coi/cox1 specific matching
                if target_gene == "cox1":
                    if (
                        "coi" in qualifier_lower.split()
                        or "co1" in qualifier_lower.split()
                    ):
                        logger.info(
                            f"Found COI/CO1 match for cox1 in qualifier: {qualifier}"
                        )
                        found_cds = feature
                        break

            # If exact match found then break loop
            if found_cds:
                break

            # If no exact match found, try the more general matching for any variant
        if not found_cds:
            for feature in record.features:
                if feature.type != "CDS":
                    continue

                qualifiers = []
                for field in ["gene", "product", "note"]:
                    qualifiers.extend(feature.qualifiers.get(field, []))

                # Check if any variation matches in the qualifiers
                if any(
                    var in qualifier.lower()
                    for qualifier in qualifiers
                    for var in gene_variations
                ):
                    logger.info(f"Found match using variations for {gene_name}")
                    found_cds = feature
                    break

        # If matching CDS found then extract it
        if found_cds:
            try:
                cds_record = record[:]
                cds_record.seq = found_cds.extract(record.seq)

                if len(cds_record.seq) >= 100:  # Sanity check for minimum length
                    logger.info(
                        f"Successfully extracted CDS of length {len(cds_record.seq)} "
                        f"(accession {record.id})"
                    )
                    return cds_record
                else:
                    logger.warning(
                        f"Extracted CDS too short ({len(cds_record.seq)}bp) "
                        f"(accession {record.id})"
                    )
            except Exception as e:
                logger.error(f"CDS extraction error for {record.id}: {e}")
                logger.error("Full error details:", exc_info=True)

        # If not in single mode then don't use fallbacks
        if not single_mode:
            logger.debug(
                f"No valid CDS found for gene {gene_name} (accession {record.id})"
            )
            return None

        logger.info(
            f"No CDS feature found, trying fallbacks for single mode "
            f"(accession {record.id})"
        )

        # Define reasonable size limits for different genes - only used in single mode
        max_gene_sizes = {
            "rbcl": 2000,  # Typical rbcL is ~1400bp
            "cox1": 2000,  # Typical cox1 is ~1500bp
            "cox2": 2000,  # Typical cox2 is ~1500bp
            "cox3": 2000,  # Typical cox3 is ~1500bp
            "cytb": 1800,  # Typical cytb is ~1100bp
            "nd1": 1800,  # Typical nd1 is ~1000bp
            "nd2": 1800,  # Typical nd2 is ~1000bp
            "nd4": 1800,  # Typical nd4 is ~1300bp
            "nd5": 2000,  # Typical nd5 is ~1700bp
            "matk": 2000,  # Typical matK is ~1500bp
            "atp6": 1200,  # Typical atp6 is ~800bp
            "atp8": 1000,  # Typical atp8 is ~400bp
            "16s": 2000,  # Typical 16S is ~1600bp
            "18s": 2500,  # Typical 18S is ~1800bp
            "28s": 3500,  # Typical 28S can be ~3000bp
            "12s": 1500,  # Typical 12S is ~1000bp
            "its": 3000,  # Typical ITS region is variable in length
            "its1": 1500,  # Typical ITS1 is variable in length
            "its2": 1500,  # Typical ITS2 is variable in length
            "trnh-psba": 1000,  # Typical trnH-psbA is ~500-700bp
        }

        # Get maximum acceptable size for this gene
        max_size = max_gene_sizes.get(
            base_gene, 3000
        )  # Default to 3000 for unknown genes

        # FALLBACK 1: Check for gene feature with matching name but no CDS
        for feature in record.features:
            if feature.type == "gene":
                gene_qualifiers = feature.qualifiers.get("gene", [])
                gene_notes = feature.qualifiers.get("note", [])
                all_qualifiers = gene_qualifiers + gene_notes

                # Check for exact match first
                if target_gene in [q.lower() for q in all_qualifiers]:
                    logger.info(f"Found exact gene match: {target_gene}")
                    try:
                        gene_record = record[:]
                        gene_record.seq = feature.extract(record.seq)

                        if len(gene_record.seq) > max_size:
                            logger.warning(
                                f"Extracted gene region too large "
                                f"({len(gene_record.seq)} bp > {max_size} bp limit) "
                                f"- skipping"
                            )
                            continue

                        if len(gene_record.seq) >= min_size:
                            logger.info(
                                f"Successfully extracted gene region of length "
                                f"{len(gene_record.seq)}"
                            )
                            return gene_record
                        else:
                            logger.warning(
                                f"Extracted gene region too short "
                                f"({len(gene_record.seq)} bp < {min_size} bp)"
                            )
                    except Exception as e:
                        logger.error(f"Gene region extraction error: {e}")

                # More general matching
                qualifier_text = " ".join(all_qualifiers).lower()

                # Check if any variation matches in the qualifiers
                if any(var in qualifier_text for var in gene_variations):
                    try:
                        logger.info(
                            f"Found matching gene feature, using gene region "
                            f"(accession {record.id})"
                        )
                        gene_record = record[:]
                        gene_record.seq = feature.extract(record.seq)

                        # Check if the extracted sequence is too large
                        if len(gene_record.seq) > max_size:
                            logger.warning(
                                f"Extracted gene region too large "
                                f"({len(gene_record.seq)} bp > {max_size} bp limit) "
                                f"- skipping (accession {record.id})"
                            )
                            continue

                        if len(gene_record.seq) >= min_size:
                            logger.info(
                                f"Successfully extracted gene region of length "
                                f"{len(gene_record.seq)} (accession {record.id})"
                            )
                            return gene_record
                        else:
                            logger.warning(
                                f"Extracted gene region too short "
                                f"({len(gene_record.seq)} bp < {min_size} bp) "
                                f"(accession {record.id})"
                            )
                    except Exception as e:
                        logger.error(
                            f"Gene region extraction error for {record.id}: {e}"
                        )

        # FALLBACK 2: Check if it is an mRNA sequence with no CDS feature
        mol_type = ""
        for feature in record.features:
            if feature.type == "source" and "mol_type" in feature.qualifiers:
                mol_type = feature.qualifiers["mol_type"][0].lower()
                break

        if mol_type in ["mrna", "est"]:
            logger.info(
                f"Sequence is mRNA/EST, checking if description matches target gene "
                f"(accession {record.id})"
            )

            description_lower = record.description.lower()

            # Check if any of our variations appear in the description
            matching_vars = [var for var in gene_variations if var in description_lower]
            if matching_vars:
                logger.info(
                    f"mRNA/EST description matches gene variations: {matching_vars}"
                )

                # Check if the sequence is too large
                if len(record.seq) > max_size:
                    logger.warning(
                        f"mRNA/EST sequence too large ({len(record.seq)} bp > "
                        f"{max_size} bp limit) - skipping (accession {record.id})"
                    )
                    return None

                if len(record.seq) >= min_size:
                    logger.info(
                        f"Using complete mRNA/EST of length {len(record.seq)} as "
                        f"it matches gene variations (accession {record.id})"
                    )
                    return record
                else:
                    logger.warning(
                        f"mRNA/EST sequence too short ({len(record.seq)} bp < "
                        f"{min_size} bp) (accession {record.id})"
                    )
            else:
                logger.info("Description doesn't match any target gene variation")

        # FALLBACK 3: If it's a partial sequence entry, check if gene name appears in description
        description_lower = record.description.lower()
        if "partial" in description_lower:
            matching_vars = [var for var in gene_variations if var in description_lower]
            if matching_vars:
                logger.info(
                    f"Found partial sequence matching gene variations: {matching_vars}"
                )

                # Check if the sequence is too large
                if len(record.seq) > max_size:
                    logger.warning(
                        f"Partial sequence too large ({len(record.seq)} bp > "
                        f"{max_size} bp limit) - skipping (accession {record.id})"
                    )
                    return None

                if len(record.seq) >= min_size:
                    logger.info(
                        f"Using entire partial sequence of length {len(record.seq)} "
                        f"(accession {record.id})"
                    )
                    return record
                else:
                    logger.warning(
                        f"Partial sequence too short ({len(record.seq)} bp < "
                        f"{min_size} bp) (accession {record.id})"
                    )

        # FALLBACK 4: For all records, check if the target gene is in the organism name or sequence ID
        # This is a last resort when in single-taxid mode and are desperate for more sequences
        org_name = ""
        for feature in record.features:
            if feature.type == "source" and "organism" in feature.qualifiers:
                org_name = feature.qualifiers["organism"][0].lower()
                break

        if base_gene in record.id.lower() or base_gene in org_name:
            logger.info(
                f"Last resort: Gene name {base_gene} found in sequence ID or "
                f"organism name (accession {record.id})"
            )

            # Check if the sequence is too large
            if len(record.seq) > max_size:
                logger.warning(
                    f"Sequence too large ({len(record.seq)} bp > {max_size} bp "
                    f"limit) - skipping as last resort (accession {record.id})"
                )
                return None

            if len(record.seq) >= min_size:
                logger.info(
                    f"Using entire sequence of length {len(record.seq)} as a "
                    f"last resort (accession {record.id})"
                )
                return record
            else:
                logger.warning(
                    f"Sequence too short ({len(record.seq)} bp < {min_size} bp) "
                    f"- skipping as last resort (accession {record.id})"
                )

        logger.debug(
            f"No valid CDS or fallback found for gene {gene_name} "
            f"(accession {record.id})"
        )
        return None

    # Parse complex coded_by expressions, including complement() and join() statements
    def parse_coded_by(
        self, coded_by: str
    ) -> Tuple[Optional[List[Tuple[str, Optional[Tuple[int, int]]]]], bool]:
        logger.info(f"Parsing coded_by qualifier: {coded_by}")
        try:
            # Determine if complement first
            is_complement = coded_by.startswith("complement(")

            # Remove outer wrapper (complement or join)
            if is_complement:
                coded_by = coded_by[10:-1]  # Remove 'complement(' and final ')'
            elif coded_by.startswith("join("):
                coded_by = coded_by[5:-1]  # Remove 'join(' and final ')'

            # Split by comma while preserving full coordinates
            segments_raw = []
            current_segment = ""
            in_parentheses = 0

            # Careful splitting to preserve full coordinates
            for char in coded_by:
                if char == "(":
                    in_parentheses += 1
                elif char == ")":
                    in_parentheses -= 1
                elif char == "," and in_parentheses == 0:  # Only split at top level
                    segments_raw.append(current_segment.strip())
                    current_segment = ""
                    continue
                current_segment += char
            segments_raw.append(current_segment.strip())

            # Clean up the segments
            cleaned_segments = []
            for seg in segments_raw:
                seg = seg.strip().strip("()")
                cleaned_segments.append(seg)

            logger.debug(f"Cleaned segments: {cleaned_segments}")

            # Process all segments
            result = []
            all_coordinates_valid = True

            for segment in cleaned_segments:
                if not segment:  # Skip empty segments
                    continue

                logger.debug(f"Processing segment: '{segment}'")

                # Extract accession and coordinates
                if ":" in segment:
                    accession, coords = segment.split(":")
                    accession = accession.strip()
                    coords = coords.strip()
                    logger.debug(
                        f"Split into accession: '{accession}', coords: '{coords}'"
                    )

                    if ".." in coords:
                        coord_parts = coords.split("..")
                        if len(coord_parts) != 2:
                            logger.error(f"Invalid coordinate format: {coords}")
                            all_coordinates_valid = False
                            break

                        start_str, end_str = coord_parts
                        # Remove any non-digit characters
                        start_str = "".join(c for c in start_str if c.isdigit())
                        end_str = "".join(c for c in end_str if c.isdigit())

                        logger.debug(
                            f"Cleaned coordinate strings - Start: '{start_str}', "
                            f"End: '{end_str}'"
                        )

                        try:
                            start = int(start_str)
                            end = int(end_str)
                            logger.debug(
                                f"Parsed coordinates - Start: {start}, End: {end}"
                            )

                            if start <= 0 or end <= 0:
                                logger.error(
                                    "Invalid coordinates: must be positive numbers"
                                )
                                all_coordinates_valid = False
                                break

                            if start > end:
                                logger.error(
                                    f"Invalid coordinate range: {start}..{end}"
                                )
                                all_coordinates_valid = False
                                break

                            result.append((accession, (start, end)))
                            logger.info(
                                f"Successfully parsed coordinates: {start}-{end} "
                                f"for {accession}"
                            )

                        except ValueError as ve:
                            logger.error(
                                f"Failed to parse coordinates '{coords}': {ve}"
                            )
                            all_coordinates_valid = False
                            break
                    else:
                        logger.error(f"Missing coordinate separator '..' in {coords}")
                        all_coordinates_valid = False
                        break
                else:
                    logger.error(f"Missing accession separator ':' in {segment}")
                    all_coordinates_valid = False
                    break

            if not all_coordinates_valid or not result:
                logger.error("Failed to parse one or more segments")
                return None, False

            logger.debug(f"Successfully parsed {len(result)} segments")
            return result, is_complement

        except Exception as e:
            logger.error(f"Error parsing coded_by: {coded_by}, error: {e}")
            logger.error("Full error details:", exc_info=True)
            return None, False

    # Fetch nucleotide sequence corresponding to a protein record, handling both RefSeq coded_by qualifiers and UniProt xrefs.
    def fetch_nucleotide_from_protein(
        self, protein_record: SeqRecord, gene_name: str
    ) -> Optional[SeqRecord]:
        try:
            logger.info(
                f"Attempting to fetch nucleotide sequence from protein record "
                f"{protein_record.id}"
            )

            # Try coded_by qualifier for RefSeq records
            cds_feature = next(
                (f for f in protein_record.features if f.type == "CDS"), None
            )
            if cds_feature and "coded_by" in cds_feature.qualifiers:
                coded_by = cds_feature.qualifiers["coded_by"][0]

                parsed_result = self.parse_coded_by(coded_by)
                if parsed_result:
                    segments, is_complement = parsed_result

                    if not segments:
                        logger.error(
                            f"No valid segments found in coded_by qualifier for "
                            f"{protein_record.id}"
                        )
                        return None

                    # Fetch full sequence for first accession
                    first_accession = segments[0][0]
                    logger.info(
                        f"Fetching nucleotide sequence for accession: "
                        f"{first_accession}"
                    )

                    # Use enhanced nucleotide fetching
                    nucleotide_record = self.fetch_nucleotide_record(first_accession)
                    if not nucleotide_record:
                        logger.error(
                            f"Failed to fetch nucleotide sequence for accession: "
                            f"{first_accession}"
                        )
                        return None

                    # Extract and join all segments
                    complete_sequence = ""
                    for accession, coordinates in segments:
                        if accession != first_accession:
                            nucleotide_record = self.fetch_nucleotide_record(accession)
                            if not nucleotide_record:
                                logger.error(
                                    f"Failed to fetch additional sequence: "
                                    f"{accession}"
                                )
                                continue

                        if coordinates:
                            start, end = coordinates
                            segment_seq = str(nucleotide_record.seq[start - 1 : end])
                            if len(segment_seq) == 0:
                                logger.error(
                                    f"Zero-length sequence extracted using "
                                    f"coordinates {start}..{end} "
                                    f"(accession {accession})"
                                )
                                return None
                        else:
                            segment_seq = str(nucleotide_record.seq)

                        complete_sequence += segment_seq

                        # Handle complement if needed
                    if is_complement:
                        complete_sequence = str(
                            Seq(complete_sequence).reverse_complement()
                        )

                    # Create new record with complete sequence
                    new_record = nucleotide_record[:]
                    new_record.seq = Seq(complete_sequence)
                    logger.info(
                        f"***Successfully extracted nucleotide sequence: "
                        f"Length {len(complete_sequence)} "
                        f"(from protein record: {protein_record.id})"
                    )

                    return new_record

            logger.warning(
                f"No valid nucleotide reference found in protein record "
                f"{protein_record.id}"
            )
            return None

        except Exception as e:
            logger.error(
                f"Error in fetch_nucleotide_from_protein for "
                f"{protein_record.id}: {e}"
            )
            logger.error("Full error details:", exc_info=True)
            return None

    # Central nucleotide & protein search function using fetched taxid
    def try_fetch_at_taxid(
        self,
        current_taxid: str,
        rank_name: str,
        taxon_name: str,
        sequence_type: str,
        gene_name: str,
        protein_records: List[SeqRecord],
        nucleotide_records: List[SeqRecord],
        best_taxonomy: List[str],
        best_matched_rank: Optional[str],
        fetch_all: bool = False,
        progress_counters: Optional[Dict[str, int]] = None,
    ) -> Tuple[bool, bool, List[str], Optional[str], List[SeqRecord], List[SeqRecord]]:
        protein_found = False
        nucleotide_found = False

        # Initialise progress tracking if provided
        sequence_counter = (
            progress_counters.get("sequence_counter", 0) if progress_counters else 0
        )
        max_sequences = (
            progress_counters.get("max_sequences", None) if progress_counters else None
        )

        # Set minimum protein size for single mode
        min_protein_size = self.config.min_protein_size_single_mode if fetch_all else 0

        try:
            # Handle protein search for 'protein' or 'both' types
            if sequence_type in ["protein", "both"] and (
                not protein_records or fetch_all
            ):
                # Modify search string based on fetch_all mode
                if fetch_all and self.config.protein_length_threshold <= 0:
                    # No size filtering when fetch_all is True and threshold is 0 or negative
                    protein_search = f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp]"
                else:
                    protein_search = (
                        f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp] "
                        f"AND {self.config.protein_length_threshold}:10000[SLEN]"
                    )

                logger.info(
                    f"Searching protein database at rank {rank_name} ({taxon_name}) with term: {protein_search}"
                )

                try:
                    protein_results = self.entrez.search(
                        db="protein", term=protein_search
                    )
                    if protein_results and protein_results.get("IdList"):
                        id_list = protein_results.get("IdList")
                        logger.info(f"Found {len(id_list)} protein records")
                        if len(id_list) > 5:  # Only log IDs if there are not too many
                            logger.info(f"Protein IDs: {id_list}")
                        
                        # Update progress_counters with actual total if in fetch_all mode
                        if fetch_all and progress_counters:
                            # If max_sequences is specified, use min(max_sequences, len(id_list))
                            # Otherwise use actual number of sequences found
                            total_sequences = min(max_sequences, len(id_list)) if max_sequences else len(id_list)
                            progress_counters['total_sequences'] = total_sequences
                        
                        # For non-fetch_all mode, apply prefiltering if there are many IDs
                        processed_ids = id_list
                        if not fetch_all and len(id_list) > 10:
                            logger.info(
                                f"Prefiltering {len(id_list)} proteins based on length information"
                            )

                            # Get summaries and sort by length
                            try:
                                sorted_summaries = []
                                batch_size = 200

                                for i in range(0, len(id_list), batch_size):
                                    batch_ids = id_list[i : i + batch_size]
                                    id_string = ",".join(batch_ids)

                                    logger.debug(
                                        f"Fetching summary for batch of {len(batch_ids)} IDs"
                                    )
                                    try:
                                        handle = Entrez.esummary(
                                            db="protein", id=id_string
                                        )
                                        batch_summaries = Entrez.read(handle)
                                        handle.close()

                                        # Extract sequence lengths from summaries
                                        for summary in batch_summaries:
                                            seq_id = summary.get("Id", "")
                                            seq_length = int(summary.get("Length", 0))
                                            sorted_summaries.append(
                                                (seq_id, seq_length)
                                            )

                                        # Add delay between batches
                                        if i + batch_size < len(id_list):
                                            sleep(uniform(0.5, 1.0))
                                    except Exception as batch_e:
                                        logger.error(
                                            f"Error in batch summary fetch: {batch_e}"
                                        )
                                        continue

                                # Check if any summaries retrieved
                                if not sorted_summaries:
                                    logger.error(
                                        "Failed to fetch any sequence summaries, using all IDs"
                                    )
                                else:
                                    # Sort by length (descending)
                                    sorted_summaries.sort(
                                        key=lambda x: x[1], reverse=True
                                    )

                                    # Take only top 250 IDs by sequence length
                                    processed_ids = [
                                        item[0] for item in sorted_summaries[:250]
                                    ]
                                    logger.info(
                                        f"Successfully filtered to top proteins by length (longest: {sorted_summaries[0][1]} aa)"
                                    )

                            except Exception as e:
                                logger.error(f"Error in prefiltering: {e}")
                                logger.error("Full error details:", exc_info=True)
                                logger.warning("Using all IDs without length filtering")

                        # Log how many IDs are processed
                        logger.info(f"Processing {len(processed_ids)} protein record")

                        # Process the filtered or complete ID list
                        for protein_id in processed_ids:
                            # Check if reached the max_sequences limit
                            if max_sequences and sequence_counter >= max_sequences:
                                logger.info(
                                    f"Reached maximum sequence limit ({max_sequences}). Stopping search."
                                )
                                break

                            # Add logging for protein fetch attempt
                            logger.info(
                                f"Attempting to fetch protein sequence for {gene_name} (GI:{protein_id})"
                            )

                            handle = self.entrez.fetch(
                                db="protein",
                                id=protein_id,
                                rettype="gb",
                                retmode="text",
                            )
                            if handle:
                                temp_record = next(SeqIO.parse(handle, "genbank"))
                                handle.close()

                                # Add logging for successful protein fetch
                                logger.info(
                                    f"***Successfully fetched protein sequence: Length {len(temp_record.seq)} (accession {temp_record.id})"
                                )

                                # Filter out false matches like "16S rRNA methylase" for rRNA targets
                                is_target_gene = True
                                if gene_name.lower() in self.config._rRNA_genes:
                                    # Check for misleading annotations
                                    for feature in temp_record.features:
                                        if (
                                            feature.type == "CDS"
                                            and "product" in feature.qualifiers
                                        ):
                                            product = feature.qualifiers["product"][
                                                0
                                            ].lower()
                                            if (
                                                f"{gene_name.lower()} rrna" in product
                                                and any(
                                                    x in product
                                                    for x in [
                                                        "methylase",
                                                        "methyltransferase",
                                                        "pseudouridylate",
                                                        "synthase",
                                                    ]
                                                )
                                            ):
                                                logger.info(
                                                    f"Skipping false match with product: {product}"
                                                )
                                                is_target_gene = False
                                                break

                                if not is_target_gene:
                                    continue

                                # Only skip UniProt/Swiss-Prot protein accession numbers in non-single mode
                                if not fetch_all:
                                    # Skip problematic UniProt/Swiss-Prot protein accession numbers
                                    if re.match(
                                        r"^[A-Z]\d+", temp_record.id
                                    ) and not re.match(r"^[A-Z]{2,}", temp_record.id):
                                        logger.info(
                                            f"Skipping UniProtKB/Swiss-Prot protein record {temp_record.id}"
                                        )
                                        continue

                                # Check minimum protein size in single mode
                                if (
                                    fetch_all
                                    and len(temp_record.seq) < min_protein_size
                                ):
                                    logger.warning(
                                        f"Protein sequence too short ({len(temp_record.seq)} aa < {min_protein_size} aa) - skipping (accession {temp_record.id})"
                                    )
                                    continue

                                if fetch_all:
                                    protein_records.append(temp_record)
                                    protein_found = True
                                    if not best_taxonomy:
                                        best_taxonomy = temp_record.annotations.get(
                                            "taxonomy", []
                                        )

                                    # Update counter and log progress
                                    sequence_counter += 1
                                    if progress_counters:

                                        progress_counters['sequence_counter'] = sequence_counter
                                    
                                    # Log progress with actual total sequences
                                    total_to_use = progress_counters.get('total_sequences', max_sequences or len(id_list))
                                    logger.info(f"====>>> Progress: {sequence_counter}/{total_to_use} sequences processed")
                                    # Also print a direct progress line for the GUI to pick up
                                    percentage = int(100 * sequence_counter / total_to_use)
                                    print(f"progress: {percentage}%")
                                        progress_counters["sequence_counter"] = (
                                            sequence_counter
                                        )

                                    # Log progress
                                    if max_sequences:
                                        logger.info(
                                            f"====>>> Progress: {sequence_counter}/{max_sequences} sequences processed"
                                        )
                                    else:
                                        # If max_sequences is None, use the total found sequences
                                        total_found = len(id_list)
                                        logger.info(
                                            f"Progress: {sequence_counter}/{total_found} sequences processed"
                                        )

                                else:
                                    # Keep only longest sequence
                                    if not protein_records or len(
                                        temp_record.seq
                                    ) > len(protein_records[0].seq):
                                        protein_records.clear()
                                        protein_records.append(temp_record)
                                        protein_found = True
                                        best_taxonomy = temp_record.annotations.get(
                                            "taxonomy", []
                                        )

                            # For batch mode (--type both), try to fetch corresponding nucleotide
                            if (
                                protein_found
                                and not fetch_all
                                and sequence_type == "both"
                            ):
                                nucleotide_record = self.fetch_nucleotide_from_protein(
                                    protein_records[0], gene_name
                                )
                                if nucleotide_record:
                                    nucleotide_records.clear()
                                    nucleotide_records.append(nucleotide_record)
                                    nucleotide_found = True
                                    logger.debug(
                                        "Successfully fetched corresponding nucleotide sequence"
                                    )
                                    break  # Exit loop after finding the first valid protein and nucleotide pair
                                else:
                                    logger.warning(
                                        "Failed to fetch corresponding nucleotide sequence"
                                    )
                                    protein_records.clear()
                                    protein_found = False

                except Exception as e:
                    logger.error(f"Error searching protein database: {e}")

            # Handle nucleotide search
            if (
                (sequence_type == "nucleotide")
                or (sequence_type == "both" and fetch_all)  # Single taxid mode
                or (sequence_type == "both" and not nucleotide_found)
            ):  # Fallback for batch mode

                # Reset counter if switching to nucleotide search in 'both' mode
                if fetch_all and sequence_type == "both" and protein_records:
                    sequence_counter = 0
                    if progress_counters:
                        progress_counters["sequence_counter"] = sequence_counter

                # Modify search string based on fetch_all mode and add exclusion terms for rRNA genes
                search_exclusions = ""
                if gene_name.lower() in self.config._rRNA_genes:
                    search_exclusions = " NOT methylase[Title] NOT methyltransferase[Title] NOT pseudouridylate[Title] NOT synthase[Title]"

                if fetch_all:
                    nucleotide_search = f"{self.config.gene_search_term}{search_exclusions} AND txid{current_taxid}[Organism:exp]"
                else:
                    nucleotide_search = (
                        f"{self.config.gene_search_term}{search_exclusions} AND txid{current_taxid}[Organism:exp] "
                        f"AND {self.config.nucleotide_length_threshold}:60000[SLEN]"
                    )

                logger.info(
                    f"Searching nucleotide database at rank {rank_name} ({taxon_name}) with term: {nucleotide_search}"
                )

                try:
                    nucleotide_results = self.entrez.search(
                        db="nucleotide", term=nucleotide_search
                    )
                    if nucleotide_results and nucleotide_results.get("IdList"):
                        id_list = nucleotide_results.get("IdList")
                        logger.info(f"Found {len(id_list)} nucleotide sequence IDs")
                        if len(id_list) > 5:  # Only log IDs if there are not too many
                            logger.debug(f"Nucleotide IDs: {id_list}")
                        
                        # Update the progress_counters with the actual total in fetch_all mode
                        if fetch_all and progress_counters:
                            # If max_sequences is specified, use min(max_sequences, len(id_list))
                            # Otherwise use the actual number of sequences found
                            total_sequences = min(max_sequences, len(id_list)) if max_sequences else len(id_list)
                            progress_counters['total_sequences'] = total_sequences

                        # Apply the same prefiltering optimisation for nucleotide sequences
                        processed_ids = id_list
                        if not fetch_all and len(id_list) > 10:
                            logger.info(
                                f"Prefiltering {len(id_list)} nucleotide sequences based on length information"
                            )

                            # Get summaries and sort by length
                            try:
                                sorted_summaries = []
                                batch_size = 200  # Fetch in batches of 200

                                for i in range(0, len(id_list), batch_size):
                                    batch_ids = id_list[i : i + batch_size]
                                    id_string = ",".join(batch_ids)

                                    logger.debug(
                                        f"Fetching summary for batch of {len(batch_ids)} IDs"
                                    )
                                    try:
                                        handle = Entrez.esummary(
                                            db="nucleotide", id=id_string
                                        )
                                        batch_summaries = Entrez.read(handle)
                                        handle.close()

                                        # Extract sequence lengths from summaries
                                        for summary in batch_summaries:
                                            seq_id = summary.get("Id", "")
                                            seq_length = int(summary.get("Length", 0))
                                            sorted_summaries.append(
                                                (seq_id, seq_length)
                                            )

                                        # Add delay between batches
                                        if i + batch_size < len(id_list):
                                            sleep(uniform(0.5, 1.0))
                                    except Exception as batch_e:
                                        logger.error(
                                            f"Error in batch summary fetch: {batch_e}"
                                        )
                                        continue

                                # Check if any summaries retrieved
                                if not sorted_summaries:
                                    logger.error(
                                        "Failed to fetch any sequence summaries, using all IDs"
                                    )
                                else:
                                    # Sort by length (descending)
                                    sorted_summaries.sort(
                                        key=lambda x: x[1], reverse=True
                                    )

                                    # Take only top 250 IDs by sequence length
                                    processed_ids = [
                                        item[0] for item in sorted_summaries[:250]
                                    ]
                                    logger.info(
                                        f"Successfully filtered to top nucleotide sequences by length (longest: {sorted_summaries[0][1]} bp)"
                                    )

                            except Exception as e:
                                logger.error(f"Error in nucleotide prefiltering: {e}")
                                logger.error("Full error details:", exc_info=True)
                                logger.warning("Using all IDs without length filtering")

                        # Log how many IDs are processed
                        logger.info(f"Processing {len(processed_ids)} nucleotide IDs")

                        for seq_id in processed_ids:
                            # Check if reached the max_sequences limit
                            if max_sequences and sequence_counter >= max_sequences:
                                logger.info(
                                    f"Reached maximum sequence limit ({max_sequences}). Stopping search."
                                )
                                break

                            try:
                                logger.info(
                                    f"Attempting to fetch nucleotide sequence (GI:{seq_id})"
                                )
                                temp_record = self.fetch_nucleotide_record(seq_id)

                                if temp_record:
                                    logger.info(
                                        f"Successfully fetched nucleotide sequence of length {len(temp_record.seq)} (accession {temp_record.id})"
                                    )

                                    # Check for misleading annotations even after search filtering
                                    is_target_gene = True
                                    if gene_name.lower() in self.config._rRNA_genes:
                                        for feature in temp_record.features:
                                            if (
                                                feature.type == "CDS"
                                                and "product" in feature.qualifiers
                                            ):
                                                product = feature.qualifiers["product"][
                                                    0
                                                ].lower()
                                                if (
                                                    f"{gene_name.lower()} rrna"
                                                    in product
                                                    and any(
                                                        x in product
                                                        for x in [
                                                            "methylase",
                                                            "methyltransferase",
                                                            "pseudouridylate",
                                                            "synthase",
                                                        ]
                                                    )
                                                ):
                                                    logger.info(
                                                        f"Skipping record with misleading product: {product}"
                                                    )
                                                    is_target_gene = False
                                                    break

                                    if not is_target_gene:
                                        continue

                                    if gene_name.lower() in self.config._rRNA_genes:
                                        # For rRNA genes, extract the specific rRNA feature
                                        rRNA_record = self.extract_rRNA(
                                            temp_record, gene_name, fetch_all
                                        )

                                        if rRNA_record:
                                            # Only proceed if we successfully extracted the rRNA feature
                                            logger.info(
                                                f"Successfully extracted {gene_name} rRNA feature of length {len(rRNA_record.seq)} from {temp_record.id}"
                                            )
                                            if fetch_all:
                                                nucleotide_records.append(rRNA_record)
                                                nucleotide_found = True
                                                if not best_taxonomy:
                                                    best_taxonomy = (
                                                        temp_record.annotations.get(
                                                            "taxonomy", []
                                                        )
                                                    )

                                                # Update counter and log progress
                                                sequence_counter += 1
                                                if progress_counters:

                                                    progress_counters['sequence_counter'] = sequence_counter
                                                
                                                # Log progress with actual total sequences
                                                total_to_use = progress_counters.get('total_sequences', max_sequences or len(id_list))
                                                logger.info(f"Progress: {sequence_counter}/{total_to_use} sequences processed")
                                                # Also print a direct progress line for the GUI to pick up
                                                percentage = int(100 * sequence_counter / total_to_use)
                                                print(f"progress: {percentage}%")
                                                    progress_counters[
                                                        "sequence_counter"
                                                    ] = sequence_counter

                                                # Log progress
                                                if max_sequences:
                                                    logger.info(
                                                        f"Progress: {sequence_counter}/{max_sequences} sequences processed"
                                                    )
                                                else:
                                                    # If max_sequences is None, use the total found sequences
                                                    total_found = len(id_list)
                                                    logger.info(
                                                        f"Progress: {sequence_counter}/{total_found} sequences processed"
                                                    )

                                            else:
                                                # Keep only longest rRNA
                                                if not nucleotide_records or len(
                                                    rRNA_record.seq
                                                ) > len(nucleotide_records[0].seq):
                                                    nucleotide_records.clear()
                                                    nucleotide_records.append(
                                                        rRNA_record
                                                    )
                                                    nucleotide_found = True
                                                    best_taxonomy = (
                                                        temp_record.annotations.get(
                                                            "taxonomy", []
                                                        )
                                                    )
                                                    logger.info(
                                                        f"Found valid {gene_name} rRNA sequence. Stopping search since sequences are sorted by length."
                                                    )
                                                    break  # Exit the loop after finding the first valid sequence
                                        else:
                                            # Skip records where no rRNA feature was found
                                            logger.info(
                                                f"No {gene_name} rRNA feature found in {temp_record.id} - skipping"
                                            )
                                            continue
                                    else:
                                        # For protein-coding genes, extract CDS
                                        logger.info(
                                            f"Attempting to extract CDS from nucleotide sequence (accession {temp_record.id})"
                                        )
                                        cds_record = self.extract_nucleotide(
                                            temp_record, gene_name, fetch_all
                                        )
                                        if cds_record:
                                            logger.info(
                                                f"Using extracted CDS in search results (accession {temp_record.id})"
                                            )
                                            if fetch_all:
                                                nucleotide_records.append(cds_record)
                                                nucleotide_found = True
                                                if not best_taxonomy:
                                                    best_taxonomy = (
                                                        temp_record.annotations.get(
                                                            "taxonomy", []
                                                        )
                                                    )

                                                # Update counter and log progress
                                                sequence_counter += 1
                                                if progress_counters:

                                                    progress_counters['sequence_counter'] = sequence_counter
                                                
                                                # Log progress with actual total sequences
                                                total_to_use = progress_counters.get('total_sequences', max_sequences or len(id_list))
                                                logger.info(f"Progress: {sequence_counter}/{total_to_use} sequences processed")
                                                # Also print a direct progress line for the GUI to pick up
                                                percentage = int(100 * sequence_counter / total_to_use)
                                                print(f"progress: {percentage}%")
                                                    progress_counters[
                                                        "sequence_counter"
                                                    ] = sequence_counter

                                                # Log progress
                                                if max_sequences:
                                                    logger.info(
                                                        f"Progress: {sequence_counter}/{max_sequences} sequences processed"
                                                    )
                                                else:
                                                    # If max_sequences is None, use the total found sequences
                                                    total_found = len(id_list)
                                                    logger.info(
                                                        f"Progress: {sequence_counter}/{total_found} sequences processed"
                                                    )

                                            else:
                                                # Keep only longest CDS
                                                if not nucleotide_records or len(
                                                    cds_record.seq
                                                ) > len(nucleotide_records[0].seq):
                                                    nucleotide_records.clear()
                                                    nucleotide_records.append(
                                                        cds_record
                                                    )
                                                    nucleotide_found = True
                                                    best_taxonomy = (
                                                        temp_record.annotations.get(
                                                            "taxonomy", []
                                                        )
                                                    )
                                                    logger.info(
                                                        f"Found valid {gene_name} CDS sequence. Stopping search since sequences are sorted by length."
                                                    )
                                                    break  # Exit the loop after finding the first valid sequence
                                        else:
                                            logger.warning(
                                                f"Failed to extract CDS from nucleotide sequence (accession {temp_record.id})"
                                            )
                                            continue

                            except Exception as e:
                                logger.error(f"Error processing sequence {seq_id}: {e}")
                                continue
                except Exception as e:
                    logger.error(f"Error searching nucleotide database: {e}")
                    nucleotide_results = None

            if protein_found or nucleotide_found:
                current_match = (
                    f"{rank_name}:{taxon_name}"
                    if rank_name
                    else f"exact match:{taxon_name}"
                )
                if not best_matched_rank or (
                    rank_name and not best_matched_rank.startswith("exact")
                ):
                    best_matched_rank = current_match

        except Exception as e:
            logger.error(f"Error in try_fetch_at_taxid for taxid {current_taxid}: {e}")
            logger.error("Full error details:", exc_info=True)


        return protein_found, nucleotide_found, best_taxonomy, best_matched_rank, protein_records, nucleotide_records
    
        return (
            protein_found,
            nucleotide_found,
            best_taxonomy,
            best_matched_rank,
            protein_records,
            nucleotide_records,
        )

    # Handles taxonomy traversal (i.e. taxonomy walking) for target sequence searches
    def search_and_fetch_sequences(
        self,
        taxid: str,
        gene_name: str,
        sequence_type: str,
        fetch_all: bool = False,
        progress_counters: Optional[Dict[str, int]] = None,
    ) -> Tuple[List[SeqRecord], List[SeqRecord], List[str], str]:
        # Initialise empty lists for records
        protein_records = []
        nucleotide_records = []
        best_taxonomy = []
        best_matched_rank = None

        # Fetch taxonomy first (from cache if available)
        logger.debug(f"Starting sequence search for {gene_name} using taxid {taxid}")
        taxonomy, taxon_ranks, initial_rank, taxon_ids = self.entrez.fetch_taxonomy(taxid)
        if not taxonomy:
            logger.error(f"Could not fetch taxonomy for taxID ({taxid}), cannot search for sequences")
            return [], [], [], "No taxonomy found"

        # Get ordered list of ranks to traverse
        current_taxonomy = taxonomy[:]
        current_taxon = current_taxonomy.pop()  # Start with species
        current_rank = taxon_ranks.get(current_taxon, "unknown")
        current_taxid = taxid

        # Traverse taxonomy from species up
        while True:
            logger.info(
                f"Attempting search at {current_rank} level: {current_taxon} (taxid: {current_taxid})"
            )

            (
                protein_found,
                nucleotide_found,
                best_taxonomy,
                best_matched_rank,
                protein_records,
                nucleotide_records,
            ) = self.try_fetch_at_taxid(
                current_taxid,
                current_rank,
                current_taxon,
                sequence_type,
                gene_name,
                protein_records,
                nucleotide_records,
                best_taxonomy,
                best_matched_rank,
                fetch_all,
                progress_counters,
            )

            # For single-taxid mode with fetch_all, only search at the exact taxid level
            if fetch_all:
                break

            # For batch mode, continue searching up taxonomy if needed for different sequence 'types'
            found_required_sequences = False
            if sequence_type == "both":
                found_required_sequences = protein_records and nucleotide_records
            elif sequence_type == "protein":
                found_required_sequences = bool(protein_records)
            elif sequence_type == "nucleotide":
                found_required_sequences = bool(nucleotide_records)

            if found_required_sequences:
                break

            # Log and continue up taxonomy if not found
            logger.info(
                f"No sequences found at {current_rank} level ({current_taxon}), traversing up taxonomy"
            )

            # Stop if no more levels or rank is too high
            if (
                current_rank
                in ["class", "subphylum", "phylum", "kingdom", "superkingdom"]
                or current_taxon == "cellular organisms"
                or not current_taxonomy
            ):
                logger.info(
                    f"Reached {current_rank} rank, stopping taxonomic traversal"
                )
                break

            # Update current variables
            current_taxon = current_taxonomy.pop()
            current_rank = taxon_ranks.get(current_taxon, "unknown")
            current_taxid = taxon_ids.get(current_taxon)
            if not current_taxid:
                continue

            # Add delay between attempts
            sleep(uniform(1, 2))

        # Set final matched rank
        matched_rank = best_matched_rank if best_matched_rank else "No match"

        # Different return logic based on mode
        if fetch_all:
            # Single taxid mode: return what was found
            if not protein_records and not nucleotide_records:
                logger.warning("No sequences found")
                return [], [], [], "No match"
            logger.info(
                f"Single taxid mode: Found {len(protein_records)} protein and {len(nucleotide_records)} nucleotide sequences"
            )
            return (
                protein_records,
                nucleotide_records,
                best_taxonomy,
                matched_rank,
            )
        else:
            # Batch mode: require both for 'both' type
            if sequence_type == "both" and (
                not protein_records or not nucleotide_records
            ):
                logger.warning(
                    "Failed to find both protein and corresponding nucleotide sequence"
                )
                return [], [], [], "No match"
            elif sequence_type == "protein" and not protein_records:
                logger.warning("No protein sequence found")
                return [], [], [], "No match"
            elif sequence_type == "nucleotide" and not nucleotide_records:
                logger.warning("No nucleotide sequence found")
                return [], [], [], "No match"

        logger.info(f"Search completed! Matched at rank: {matched_rank}")
        return protein_records, nucleotide_records, best_taxonomy, matched_rank

    # Extract rRNA feature of specified type from record
    def extract_rRNA(self, record, gene_name, single_mode=False):
        # Get minimum size threshold for single-taxid mode or batch mode
        min_size = self.config.min_nucleotide_size_single_mode if single_mode else 100

        # Convert gene name to lowercase and normalise
        rRNA_type = (
            gene_name.lower()
            .replace("s rrna", "s")
            .replace(" rrna", "")
            .replace("rrn", "")
        )

        # Handle special gene name synonyms for different rRNA types
        if rRNA_type == "rrs":
            rRNA_type = "16s"
        elif rRNA_type == "rrl":
            rRNA_type = "23s"
        elif rRNA_type == "mt-rrn1":
            rRNA_type = "12s"
        elif rRNA_type == "ssu":
            rRNA_type = "18s"
        elif rRNA_type == "lsu" and gene_name.lower() != "5s":
            # LSU usually refers to 23S in bacteria or 28S in eukaryotes,
            # so bias toward the user's query if available
            if "23" in gene_name:
                rRNA_type = "23s"
            elif "28" in gene_name:
                rRNA_type = "28s"
            else:
                # LSU without specification - will match either
                rRNA_type = "lsu"

        logger.info(f"Looking for {rRNA_type} rRNA features in record {record.id}")

        # Define alternative names for different rRNA types
        rRNA_alternatives = {
            "16s": ["16s", "rrs", "rrn16", "ssu"],  # Small subunit bacterial
            "18s": ["18s", "rrn18", "ssu"],  # Small subunit eukaryotic
            "23s": ["23s", "rrl", "rrn23", "lsu"],  # Large subunit bacterial
            "28s": ["28s", "rrn28", "lsu"],  # Large subunit eukaryotic
            "12s": ["12s", "mt-rrn1", "mt 12s"],  # Mitochondrial SSU
            "5s": ["5s", "rrn5", "rrn5s", "rrna 5s"],  # 5S bacterial
        }

        # Get the set of alternative names for our target rRNA
        target_alternatives = set()
        for key, alternatives in rRNA_alternatives.items():
            if rRNA_type in alternatives:
                target_alternatives.update(alternatives)

        if not target_alternatives:
            # If we don't have a pre-defined set, just use the original type
            target_alternatives = {rRNA_type}

        logger.info(
            f"Searching for these rRNA name variations in feature table: {target_alternatives}"
        )

        # First look for rRNA feature type
        for feature in record.features:
            if feature.type == "rRNA":
                # Check product qualifier
                if "product" in feature.qualifiers:
                    product = feature.qualifiers["product"][0].lower()

                    # See if product matches any of our target alternatives
                    is_match = False
                    for alt in target_alternatives:
                        if alt in product and "ribosomal" in product:
                            is_match = True
                            break

                    if is_match:
                        # Skip if product also contains misleading terms
                        if any(
                            x in product
                            for x in [
                                "methylase",
                                "methyltransferase",
                                "pseudouridylate",
                                "synthase",
                            ]
                        ):
                            continue
                        logger.info(
                            f"Found matching rRNA feature with product: {product}"
                        )
                        try:
                            rRNA_record = record[:]
                            rRNA_record.seq = feature.extract(record.seq)
                            if len(rRNA_record.seq) >= min_size:
                                logger.info(
                                    f"Successfully extracted rRNA of length {len(rRNA_record.seq)}"
                                )
                                return rRNA_record
                        except Exception as e:
                            logger.error(f"Error extracting rRNA feature: {e}")

                # Check gene qualifier
                if "gene" in feature.qualifiers:
                    gene = feature.qualifiers["gene"][0].lower()

                    # Check if gene matches any of our target alternatives
                    is_match = False
                    for alt in target_alternatives:
                        if gene == alt or gene == f"{alt} rrna":
                            is_match = True
                            break

                    if is_match:
                        logger.info(f"Found matching rRNA feature with gene: {gene}")
                        try:
                            rRNA_record = record[:]
                            rRNA_record.seq = feature.extract(record.seq)
                            if len(rRNA_record.seq) >= min_size:
                                logger.info(
                                    f"Successfully extracted rRNA of length {len(rRNA_record.seq)}"
                                )
                                return rRNA_record
                        except Exception as e:
                            logger.error(f"Error extracting rRNA feature: {e}")

        # Look for gene features
        for feature in record.features:
            if feature.type == "gene":
                if "gene" in feature.qualifiers:
                    gene = feature.qualifiers["gene"][0].lower()

                    # Check if gene matches any target alternatives
                    is_match = False
                    for alt in target_alternatives:
                        if gene == alt or gene == f"{alt} rrna":
                            is_match = True
                            break

                    if is_match:
                        logger.info(f"Found matching gene feature: {gene}")
                        try:
                            rRNA_record = record[:]
                            rRNA_record.seq = feature.extract(record.seq)
                            if len(rRNA_record.seq) >= min_size:
                                logger.info(
                                    f"Successfully extracted gene region of length {len(rRNA_record.seq)}"
                                )
                                return rRNA_record
                        except Exception as e:
                            logger.error(f"Error extracting gene feature: {e}")

        logger.info(f"No specific {rRNA_type} rRNA feature found in record {record.id}")
        return None


# =============================================================================
# Output file and directory management
# =============================================================================
class OutputManager:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.nucleotide_dir = output_dir / "nucleotide"
        self.genbank_dir = output_dir / "genbank"
        self.failed_searches_path = output_dir / "failed_searches.csv"
        self.sequence_refs_path = output_dir / "sequence_references.csv"

        self._setup_directories()
        self._setup_files()

    # Create main output directories
    def _setup_directories(self):
        make_out_dir(self.output_dir)
        make_out_dir(self.nucleotide_dir)
        make_out_dir(self.genbank_dir)

    # Initialize output csv files and headers
    def _setup_files(self):
        if not self.failed_searches_path.exists():
            with open(self.failed_searches_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["process_id", "taxid", "error_type", "timestamp"])

        if not self.sequence_refs_path.exists():
            with open(self.sequence_refs_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "process_id",
                        "taxid",
                        "protein_accession",
                        "protein_length",
                        "nucleotide_accession",
                        "nucleotide_length",
                        "matched_rank",
                        "ncbi_taxonomy",
                        "reference_name",
                        "protein_reference_path",
                        "nucleotide_reference_path",
                    ]
                )

    # Log any failed sequence searches in csv
    def log_failure(self, process_id: str, taxid: str, error_type: str):
        with open(self.failed_searches_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    process_id,
                    taxid,
                    error_type,
                    time.strftime("%Y-%m-%d %H:%M:%S"),
                ]
            )

    # Write fetched sequence metadata to main output csv
    def write_sequence_reference(self, data: Dict[str, Any]):
        with open(self.sequence_refs_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    data["process_id"],
                    data["taxid"],
                    data.get("protein_accession", ""),
                    data.get("protein_length", ""),
                    data.get("nucleotide_accession", ""),
                    data.get("nucleotide_length", ""),
                    data.get("matched_rank", "unknown"),
                    data.get("taxonomy", ""),
                    data["process_id"],
                    data.get("protein_path", ""),
                    data.get("nucleotide_path", ""),
                ]
            )

    # Creates nucleotide and/or protein csv outputs for single-taxid mode
    def save_sequence_summary(self, sequences: List[SeqRecord], file_type: str):
        if not sequences:
            logger.info(f"No {file_type} sequences to summarize")
            return

        file_path = self.output_dir / f"fetched_{file_type}_sequences.csv"

        with open(file_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Accession", "Length", "Description"])

            for record in sequences:
                # Get the full description/name
                description = record.description
                accession = record.id
                length = len(record.seq)

                writer.writerow([accession, length, description])

        logger.info(f"Wrote {file_type} sequence summary to {file_path}")

# Fetch and save GenBank file for a specific record ID
def save_genbank_file(entrez: EntrezHandler, record_id: str, db: str, output_path: Path):
    try:
        logger.info(f"Fetching GenBank file for {db} record {record_id}")
        handle = entrez.fetch(
            db=db, 
            id=record_id, 
            rettype="gb", 
            retmode="text"
        )
  
        if handle:
            with open(output_path, 'w') as f:
                f.write(handle.read())
            logger.info(f"Successfully saved GenBank file to {output_path}")
            return True
        else:
             logger.warning(f"Failed to fetch GenBank file for {db} record {record_id}")
             return False
     
    except Exception as e:
        logger.error(f"Error saving GenBank file for {record_id}: {e}")
        logger.error("Full error details:", exc_info=True)
        return False

# =============================================================================
# Processing functions
# =============================================================================
# Process input per sample, retrieving and storing fetched sequences
def process_sample(
    process_id: str,
    taxid: str,
    sequence_type: str,
    processor: SequenceProcessor,
    output_manager: OutputManager,
    gene_name: str,
    save_genbank: bool = False,
) -> None:
    try:
        # Define output paths
        protein_path = output_manager.output_dir / f"{process_id}.fasta"
        nucleotide_path = output_manager.nucleotide_dir / f"{process_id}_dna.fasta"
        
        # Define GenBank paths if needed
        protein_gb_path = output_manager.genbank_dir / f"{process_id}.gb"
        nucleotide_gb_path = output_manager.genbank_dir / f"{process_id}_dna.gb"

        # Check if files already exist
        if (sequence_type in ["protein", "both"] and protein_path.exists()) or (
            sequence_type in ["nucleotide", "both"] and nucleotide_path.exists()
        ):
            logger.info(f"Sequence file(s) already exist for {process_id}. Skipping.")
            return

        # Fetch sequences (returns lists)
        protein_records, nucleotide_records, taxonomy, matched_rank = (
            processor.search_and_fetch_sequences(taxid, gene_name, sequence_type)
        )

        # Extract single records from lists for 'batch' mode
        protein_record = protein_records[0] if protein_records else None
        nucleotide_record = nucleotide_records[0] if nucleotide_records else None

        sequences_found = False
        result_data = {
            "process_id": process_id,
            "taxid": taxid,
            "matched_rank": matched_rank,
            "taxonomy": "; ".join(taxonomy) if taxonomy else "",
        }

        # Process protein sequence
        if protein_record and sequence_type in ["protein", "both"]:
            try:
                result_data["protein_accession"] = protein_record.id
                result_data["protein_length"] = len(protein_record.seq)
                result_data["protein_path"] = str(protein_path.absolute())

                logger.info(
                    f"SELECTED SEQUENCE: {protein_record.id}: Length {len(protein_record.seq)}aa"
                )

                # Store original ID for GenBank download
                original_id = protein_record.id
                
                # Write FASTA
                protein_record.id = process_id
                protein_record.description = ""
                SeqIO.write(protein_record, protein_path, "fasta")
                logger.info(f"Written protein sequence to '{protein_path}'")
                sequences_found = True
                
                # Download GenBank file if requested
                if save_genbank:
                    save_genbank_file(
                        processor.entrez,
                        original_id,
                        "protein",
                        protein_gb_path
                    )
                
            except Exception as e:
                logger.error(f"Error writing protein sequence: {e}")
                if protein_path.exists():
                    protein_path.unlink()

        # Process nucleotide sequence
        if nucleotide_record and sequence_type in ["nucleotide", "both"]:
            try:
                result_data["nucleotide_accession"] = nucleotide_record.id
                result_data["nucleotide_length"] = len(nucleotide_record.seq)
                result_data["nucleotide_path"] = str(nucleotide_path.absolute())

                logger.info(
                    f"SELECTED SEQUENCE: {nucleotide_record.id}: Length {len(nucleotide_record.seq)}bp"
                )
                
                # Store original ID for GenBank download
                original_id = nucleotide_record.id
                
                # Write FASTA
                nucleotide_record.id = process_id
                nucleotide_record.description = ""
                SeqIO.write(nucleotide_record, nucleotide_path, "fasta")
                logger.info(f"Written nucleotide sequence to '{nucleotide_path}'")
                sequences_found = True
                
                # Download GenBank file if requested
                if save_genbank:
                    save_genbank_file(
                        processor.entrez,
                        original_id,
                        "nucleotide",
                        nucleotide_gb_path
                    )
                
            except Exception as e:
                logger.error(f"Error writing nucleotide sequence: {e}")
                if nucleotide_path.exists():
                    nucleotide_path.unlink()

        if sequences_found:
            output_manager.write_sequence_reference(result_data)
        else:
            output_manager.log_failure(process_id, taxid, "No sequences found")
            logger.warning(f"No valid sequences found for taxID {taxid}")

    except Exception as e:
        logger.error(f"Error processing sample {process_id}: {e}")
        output_manager.log_failure(process_id, taxid, f"Processing error: {str(e)}")

# Process inputs for'single' taxid mode, fetching all or N available sequences
def process_single_taxid(
    taxid: str,
    gene_name: str,
    sequence_type: str,
    processor: SequenceProcessor,
    output_dir: Path,
    max_sequences: Optional[int] = None,
    save_genbank: bool = False,
) -> None:
    try:
        # Initialize progress counters
        progress_counters = {
            "sequence_counter": 0,
            "max_sequences": max_sequences,
        }

        # Fetch all sequences with progress tracking
        protein_records, nucleotide_records, taxonomy, matched_rank = (
            processor.search_and_fetch_sequences(
                taxid,
                gene_name,
                sequence_type,
                fetch_all=True,
                progress_counters=progress_counters,
            )
        )

        if not protein_records and not nucleotide_records:
            logger.warning(f"No sequences found for taxid {taxid}")
            return

        # Apply maximum sequence limit if specified
        if max_sequences is not None:
            if sequence_type in ["protein", "both"] and protein_records:
                if len(protein_records) > max_sequences:
                    logger.info(
                        f"Limiting protein records from {len(protein_records)} to {max_sequences} as specified"
                    )
                    protein_records = protein_records[:max_sequences]

            if sequence_type in ["nucleotide", "both"] and nucleotide_records:
                if len(nucleotide_records) > max_sequences:
                    logger.info(
                        f"Limiting nucleotide records from {len(nucleotide_records)} to {max_sequences} as specified"
                    )
                    nucleotide_records = nucleotide_records[:max_sequences]

        # Create output manager
        output_manager = OutputManager(output_dir)

        # Save protein sequences
        if sequence_type in ["protein", "both"] and protein_records:
            for i, record in enumerate(protein_records):
                # Store original ID for GenBank download
                original_id = record.id
                
                # Save FASTA
                filename = f"{record.id}.fasta"
                output_path = output_dir / filename
                SeqIO.write(record, output_path, "fasta")
                logger.info(
                    f"Written protein sequence {i+1}/{len(protein_records)} to '{output_path}'"
                )
                
                # Save GenBank if requested
                if save_genbank:
                    gb_path = output_manager.genbank_dir / f"{record.id}.gb"
                    save_genbank_file(
                        processor.entrez,
                        original_id,
                        "protein",
                        gb_path
                    )

            # Use output manager to save summary
            output_manager.save_sequence_summary(protein_records, "protein")
            logger.info(
                "======================================================================================="
            )
            logger.info(
                f"-----          Saved summary of {len(protein_records)} protein sequences          -----"
            )
            logger.info(
                "======================================================================================="
            )

        # Save nucleotide sequences
        if sequence_type in ["nucleotide", "both"] and nucleotide_records:
            nucleotide_dir = output_dir / "nucleotide"
            make_out_dir(nucleotide_dir)
            for i, record in enumerate(nucleotide_records):
                # Store original ID for GenBank download
                original_id = record.id
                
                # Save FASTA
                filename = f"{record.id}.fasta"
                output_path = nucleotide_dir / filename
                SeqIO.write(record, output_path, "fasta")
                logger.info(
                    f"Written nucleotide sequence {i+1}/{len(nucleotide_records)} to '{output_path}'"
                )
                
                # Save GenBank if requested
                if save_genbank:
                    gb_path = output_manager.genbank_dir / f"{record.id}.gb"
                    save_genbank_file(
                        processor.entrez,
                        original_id,
                        "nucleotide",
                        gb_path
                    )

            # Use output manager to save summary
            output_manager.save_sequence_summary(nucleotide_records, "nucleotide")
            logger.info(
                "======================================================================================="
            )
            logger.info(
                f"------      Saved summary of {len(nucleotide_records)} nucleotide sequences       -----"
            )
            logger.info(
                "======================================================================================="
            )

    except Exception as e:
        logger.error(f"Error processing taxid {taxid}: {e}")
        logger.error("Full error details:", exc_info=True)

# Process samples.csv
def process_taxid_csv(csv_path, gene_name, sequence_type, processor, output_manager, save_genbank=False):
    try:
        samples_csv = Path(csv_path)
        logger.info(f"Samples file: {samples_csv}")

        with open(samples_csv, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            process_id_col = get_process_id_column(reader.fieldnames)

            if not process_id_col:
                logger.error("Could not find process ID column in input CSV.")
                sys.exit(1)

            # Count total samples
            total_samples = sum(1 for _ in reader)
            f.seek(0)
            next(reader)

            # Initialize progress tracking
            log_progress(0, total_samples)

            # Process each sample
            for i, row in enumerate(reader, 1):
                try:
                    taxid = row["taxid"].strip()
                    process_id = row[process_id_col].strip()

                    logger.info("")
                    logger.info(
                        f"====== Processing sample {i}/{total_samples}: {process_id} (taxID: {taxid}) ======"
                    )

                    process_sample(
                        process_id=process_id,
                        taxid=taxid,
                        sequence_type=sequence_type,
                        processor=processor,
                        output_manager=output_manager,
                        gene_name=gene_name,
                        save_genbank=save_genbank,
                    )

                    # Log progress
                    log_progress(i, total_samples)

                    # Add a small delay between samples
                    sleep(uniform(0.5, 1.0))

                except Exception as e:
                    logger.error(f"Error processing row {i}: {e}")
                    continue

            # Log final progress
            log_progress(total_samples, total_samples)

    except Exception as e:
        logger.error(f"Fatal error processing taxid CSV: {e}")
        sys.exit(1)

# Process samples_taxonomy.csv
def process_taxonomy_csv(
    csv_path, gene_name, sequence_type, processor, output_manager, entrez, save_genbank=False
):
    try:
        taxonomy_csv = Path(csv_path)
        logger.info(f"Samples file: {taxonomy_csv}")

        # Read in csv file
        with open(taxonomy_csv, "r", newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)

            # Check for required columns
            required_columns = ["ID", "genus", "species"]
            missing_columns = [
                col
                for col in required_columns
                if col.lower() not in [field.lower() for field in reader.fieldnames]
            ]
            if missing_columns:
                logger.error(
                    f"Missing required columns in taxonomy CSV: {missing_columns}"
                )
                logger.error(
                    "CSV must contain at least 'ID', 'genus', and 'species' columns"
                )
                sys.exit(1)

            # Map actual column names to expected column names (case-insensitive)
            column_map = {}
            for expected_col in [
                "ID",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ]:
                for actual_col in reader.fieldnames:
                    if actual_col.lower() == expected_col.lower():
                        column_map[expected_col] = actual_col

            # Get ID column
            process_id_col = get_process_id_column(reader.fieldnames)
            if not process_id_col:
                if "ID" in column_map:
                    process_id_col = column_map["ID"]
                else:
                    logger.error("Could not find process ID column in input CSV.")
                    sys.exit(1)

            # Count total samples
            total_samples = sum(1 for _ in reader)
            f.seek(0)
            next(reader)  # Skip header

            # Initialize progress tracking
            log_progress(0, total_samples)

            # Process each sample
            for i, row in enumerate(reader, 1):
                try:
                    process_id = row[process_id_col].strip()

                    logger.info("")
                    logger.info(
                        f"====== Processing sample {i}/{total_samples}: {process_id} ======"
                    )

                    # Extract taxonomic information
                    phylum = (
                        row.get(column_map.get("phylum", ""), "").strip()
                        if "phylum" in column_map
                        else ""
                    )
                    class_name = (
                        row.get(column_map.get("class", ""), "").strip()
                        if "class" in column_map
                        else ""
                    )
                    order = (
                        row.get(column_map.get("order", ""), "").strip()
                        if "order" in column_map
                        else ""
                    )
                    family = (
                        row.get(column_map.get("family", ""), "").strip()
                        if "family" in column_map
                        else ""
                    )
                    genus = (
                        row.get(column_map.get("genus", ""), "").strip()
                        if "genus" in column_map
                        else ""
                    )
                    species = (
                        row.get(column_map.get("species", ""), "").strip()
                        if "species" in column_map
                        else ""
                    )

                    # Validate required fields
                    if not genus or not species:
                        logger.warning(
                            f"Missing required genus or species for {process_id}"
                        )
                        output_manager.log_failure(
                            process_id,
                            "unknown",
                            "Missing required taxonomy fields",
                        )
                        continue

                    # Fetch taxid from taxonomic information
                    taxid = entrez.fetch_taxid_from_taxonomy(
                        phylum, class_name, order, family, genus, species
                    )

                    if not taxid:
                        logger.warning(
                            f"Could not resolve taxid for {process_id} ({genus}, {species})"
                        )
                        output_manager.log_failure(
                            process_id, "unknown", "Could not resolve taxid"
                        )
                        continue

                    logger.info(
                        f"Starting sequence search:"
                    )
                    logger.info(
                        f"Using taxid {taxid} for sequence search of {process_id} ({genus}, {species})"
                    )

                    # Process the sample with resolved taxid
                    process_sample(
                        process_id=process_id,
                        taxid=taxid,
                        sequence_type=sequence_type,
                        processor=processor,
                        output_manager=output_manager,
                        gene_name=gene_name,
                        save_genbank=save_genbank,
                    )

                    # Log progress
                    log_progress(i, total_samples)

                    # Add a small delay between samples
                    sleep(uniform(0.5, 1.0))

                except Exception as e:
                    logger.error(f"Error processing row {i}: {e}")
                    continue

            # Log final progress
            log_progress(total_samples, total_samples)

    except Exception as e:
        logger.error(f"Fatal error processing taxonomy CSV: {e}")
        sys.exit(1)

def main():
    print("Starting gene_fetch.py")
    parser = setup_argument_parser()
    args = parser.parse_args()

    gene_name = args.gene.lower()
    output_dir = Path(args.out)
    sequence_type = args.type.lower()
    save_genbank = args.genbank  # Get genbank flag

    # Make sure output directory exists before setting up logging
    make_out_dir(output_dir)
    logger = setup_logging(output_dir)

    # Log if GenBank download is enabled
    if save_genbank:
        logger.info("GenBank download mode enabled - will save .gb files in genbank/ subdirectory")

    # Initialize components with required email/api_key
    try:
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
                save_genbank,  # Pass genbank flag
            )

    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)

    logger.info("**********************************************************")
    logger.info("           ? ?? ? Gene fetch complete ? ?? ?           ")
    logger.info("**********************************************************")


if __name__ == "__main__":
    main()
