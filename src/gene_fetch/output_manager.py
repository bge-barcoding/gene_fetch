# src/gene_fetch/output_manager.py
"""
OutputManager class for Gene Fetch.
Handles output file and directory management.
"""

import csv
import time
import logging
import traceback
from pathlib import Path
from typing import Dict, Any, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .core import logger, make_out_dir
from .entrez_handler import EntrezHandler


# =============================================================================
# Output file and directory management
# =============================================================================
class OutputManager:
    def __init__(self, output_dir: Path, save_genbank: bool = False, create_sequence_refs: bool = True):
        self.output_dir = output_dir
        self.protein_dir = output_dir / "protein"
        self.nucleotide_dir = output_dir / "nucleotide"
        self.save_genbank = save_genbank
        self.create_sequence_refs = create_sequence_refs
        
        # Only set up GenBank paths if needed
        if save_genbank:
            self.genbank_dir = output_dir / "genbank"
            self.protein_genbank_dir = self.genbank_dir / "protein"
            self.nucleotide_genbank_dir = self.genbank_dir / "nucleotide"
        else:
            self.genbank_dir = None
            self.protein_genbank_dir = None
            self.nucleotide_genbank_dir = None

        self.failed_searches_path = output_dir / "failed_searches.csv"
        
        # Only set up sequence references if needed
        if create_sequence_refs:
            self.sequence_refs_path = output_dir / "sequence_references.csv"
        else:
            self.sequence_refs_path = None

        self._setup_directories()
        self._setup_files()

    # Create main output directories
    def _setup_directories(self):
        make_out_dir(self.output_dir)
        make_out_dir(self.protein_dir)
        make_out_dir(self.nucleotide_dir)
        
        # Only create GenBank directories if needed
        if self.save_genbank:
            make_out_dir(self.genbank_dir)
            make_out_dir(self.protein_genbank_dir)
            make_out_dir(self.nucleotide_genbank_dir)

    # Initialize output csv files and headers
    def _setup_files(self):
        if not self.failed_searches_path.exists():
            with open(self.failed_searches_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["process_id", "taxid", "error_type", "timestamp"])

        # Only create sequence_references.csv if needed
        if self.create_sequence_refs and not self.sequence_refs_path.exists():
            with open(self.sequence_refs_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "process_id",
                        "input_taxa",
                        "first_matched_taxid",
                        "first_matched_taxid_rank",
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
        # Only write if sequence_references.csv was created
        if not self.create_sequence_refs:
            logger.debug("Skipping sequence reference write - not enabled for this mode")
            return
            
        with open(self.sequence_refs_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    data["process_id"],
                    data.get("input_taxa", ""),
                    data.get("first_matched_taxid", ""),
                    data.get("first_matched_taxid_rank", ""),
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
            writer.writerow(["Accession", "Length", "Description", "searched_taxid"])

            for record in sequences:
                # Get the full description/name
                description = record.description
                accession = record.id
                length = len(record.seq)
                searched_taxid = record.annotations.get("searched_taxid", "")

                writer.writerow([accession, length, description, searched_taxid])

        logger.info(f"Wrote {file_type} sequence summary to {file_path}")

    # Fetch and save GenBank file for a specific record ID
    def save_genbank_file(self, entrez: EntrezHandler, record_id: str, db: str, output_path: Path):
        """Fetch and save GenBank file for a specific record ID."""
        if not self.save_genbank:
            logger.warning("GenBank saving was not enabled - skipping file save")
            return False
        
        try:
            logger.info(f"Fetching GenBank file for {db} record {record_id}")
            handle = entrez.fetch(db=db, id=record_id, rettype="gb", retmode="text")

            if handle:
                with open(output_path, "w") as f:
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