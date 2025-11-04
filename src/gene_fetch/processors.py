# src/gene_fetch/processors.py
"""
Processing functions for Gene Fetch.
Contains functions for batch and single-taxid processing.
"""

import csv
import logging
import sys
import traceback
from pathlib import Path
from random import uniform
from time import sleep
from typing import Optional, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .core import logger, make_out_dir, log_progress, get_process_id_column
from .entrez_handler import EntrezHandler
from .sequence_processor import SequenceProcessor
from .output_manager import OutputManager


# =============================================================================
# Processing functions
# =============================================================================
# Create FASTA header based on the specified format.
def create_fasta_header(record, header_format, process_id, taxid=None, taxonomy=None):
    if header_format == "basic":
        return process_id
    elif header_format == "detailed":
        # Format: ID|taxid|accession_number|GenBank description|length
        header_parts = []
        
        # ID
        header_parts.append(process_id if process_id else "null")
        
        # taxid  
        header_parts.append(str(taxid) if taxid else "null")
        
        # accession_number
        header_parts.append(record.id if record.id else "null")
        
        # GenBank description (clean up the description)
        if hasattr(record, 'description') and record.description:
            description = record.description
            # Remove the accession number from description if it's at the start
            if record.id and description.startswith(record.id):
                description = description[len(record.id):].lstrip()
            # Remove any leading whitespace or punctuation
            description = description.lstrip(" .,;:")
            # Replace any pipe characters that would interfere with our delimiter
            description = description.replace("|", ";")
            header_parts.append(description if description else "null")
        else:
            header_parts.append("null")
            
        # length
        header_parts.append(str(len(record.seq)) if record.seq else "null")
        
        return "|".join(header_parts)
    else:
        return process_id  # fallback to basic

# Process input per sample, retrieving and storing fetched sequences
def process_sample(
    process_id: str,
    taxid: str,
    sequence_type: str,
    processor: SequenceProcessor,
    output_manager: OutputManager,
    gene_name: str,
    save_genbank: bool = False,
    header_format: str = "basic",
    input_taxa: str = "",
) -> None:
    try:
        # Define output paths
        protein_path = output_manager.protein_dir / f"{process_id}.fasta"
        nucleotide_path = output_manager.nucleotide_dir / f"{process_id}.fasta"

        # Define GenBank paths if needed
        if save_genbank:
            protein_gb_path = output_manager.protein_genbank_dir / f"{process_id}.gb"
            nucleotide_gb_path = output_manager.nucleotide_genbank_dir / f"{process_id}.gb"

        # Check if files already exist (simple resume logic)
        if (sequence_type in ["protein", "both"] and protein_path.exists()) or (
            sequence_type in ["nucleotide", "both"] and nucleotide_path.exists()
        ):
            logger.info(f"Sequence file(s) already exist for {process_id}. Skipping.")
            return

        # Fetch sequences (returns lists)
        protein_records, nucleotide_records, taxonomy, matched_rank, first_matched_taxid, first_matched_taxid_rank = (
            processor.search_and_fetch_sequences(taxid, gene_name, sequence_type)
        )

        # Extract single records from lists for 'batch' mode
        protein_record = protein_records[0] if protein_records else None
        nucleotide_record = nucleotide_records[0] if nucleotide_records else None

        sequences_found = False
        result_data = {
            "process_id": process_id,
            "input_taxa": input_taxa,
            "first_matched_taxid": first_matched_taxid if first_matched_taxid else taxid,
            "first_matched_taxid_rank": first_matched_taxid_rank if first_matched_taxid_rank else "unknown",
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

                # Create FASTA header based on format
                fasta_header = create_fasta_header(
                    protein_record, header_format, process_id, taxid, taxonomy
                )

                # Write FASTA with custom header
                protein_record.id = fasta_header
                protein_record.description = ""
                SeqIO.write(protein_record, protein_path, "fasta")
                logger.info(f"Written protein sequence to '{protein_path}'")
                sequences_found = True

                # Download GenBank file if requested
                if save_genbank:
                    output_manager.save_genbank_file(
                       processor.entrez, original_id, "protein", protein_gb_path
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

                # Create FASTA header based on format
                fasta_header = create_fasta_header(
                    nucleotide_record, header_format, process_id, taxid, taxonomy
                )

                # Write FASTA with custom header
                nucleotide_record.id = fasta_header
                nucleotide_record.description = ""
                SeqIO.write(nucleotide_record, nucleotide_path, "fasta")
                logger.info(f"Written nucleotide sequence to '{nucleotide_path}'")
                sequences_found = True

                # Download GenBank file if requested
                if save_genbank:
                    output_manager.save_genbank_file(
                       processor.entrez, original_id, "nucleotide", nucleotide_gb_path
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
    header_format: str = "basic",
) -> None:
    try:
        # Initialise progress counters
        progress_counters = {
            "sequence_counter": 0,
            "max_sequences": max_sequences,
        }

        # Fetch all sequences with progress tracking
        protein_records, nucleotide_records, taxonomy, matched_rank, first_matched_taxid, first_matched_taxid_rank = (  # UPDATED
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

        # Add searched_taxid to all records for CSV summary
        for record in protein_records:
            record.annotations["searched_taxid"] = taxid

        for record in nucleotide_records:
            record.annotations["searched_taxid"] = taxid

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

        # Create output manager (disable sequence_references.csv for single mode)
        output_manager = OutputManager(output_dir, save_genbank, create_sequence_refs=False)

        # Save protein sequences (only if protein type requested and sequences found)
        if sequence_type in ["protein", "both"] and protein_records:
            for i, record in enumerate(protein_records):
                # Store original ID for GenBank download and filename
                original_id = record.id
                
                # Create safe filename from original ID
                safe_filename = "".join(c for c in original_id if c.isalnum() or c in "._-")
                if not safe_filename:
                    safe_filename = f"sequence_{i+1}"

                # Create FASTA header based on format
                if header_format == "detailed":
                    # Use the taxonomy from the record annotations
                    record_taxonomy = record.annotations.get("taxonomy", [])
                    fasta_header = create_fasta_header(
                        record, header_format, original_id, taxid, record_taxonomy
                    )
                    record.id = fasta_header
                # For basic format, keep the original record.id (no changes needed)

                record.description = ""  # Clear description for clean FASTA

                # Save FASTA using safe filename
                filename = f"{safe_filename}.fasta"
                output_path = output_manager.protein_dir / filename
                SeqIO.write(record, output_path, "fasta")
                logger.info(
                    f"Written protein sequence {i+1}/{len(protein_records)} to '{output_path}'"
                )

                # Save GenBank if requested
                if save_genbank:
                    gb_path = output_manager.protein_genbank_dir / f"{safe_filename}.gb"
                    output_manager.save_genbank_file(processor.entrez, original_id, "protein", gb_path)

            # Save summary (only if sequences were found)
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

        # Save nucleotide sequences (only if nucleotide type requested and sequences found)
        if sequence_type in ["nucleotide", "both"] and nucleotide_records:
            for i, record in enumerate(nucleotide_records):
                # Store original ID for GenBank download and filename
                original_id = record.id
                original_description = record.description  # Save the record rescription as well
                
                # Create safe filename from original ID
                safe_filename = "".join(c for c in original_id if c.isalnum() or c in "._-")
                if not safe_filename:
                    safe_filename = f"sequence_{i+1}"

                # Create FASTA header based on format
                if header_format == "detailed":
                    # Use the taxonomy from the record annotations
                    record_taxonomy = record.annotations.get("taxonomy", [])
                    fasta_header = create_fasta_header(
                        record, header_format, original_id, taxid, record_taxonomy
                    )
                    record.id = fasta_header
                # For basic format, keep the original record.id (no changes needed)

                record.description = ""  # Clear description for a 'clean' FASTA header

                # Save FASTA using safe filename
                filename = f"{safe_filename}.fasta"
                output_path = output_manager.nucleotide_dir / filename
                SeqIO.write(record, output_path, "fasta")
                logger.info(
                    f"Written nucleotide sequence {i+1}/{len(nucleotide_records)} to '{output_path}'"
                )

                # RESTORE THE DESCRIPTION for the CSV summary
                record.description = original_description
                
                # Save GenBank if requested
                if save_genbank:
                    gb_path = output_manager.nucleotide_genbank_dir / f"{safe_filename}.gb"
                    output_manager.save_genbank_file(processor.entrez, original_id, "nucleotide", gb_path)

            # Save summary (only if sequences were found)
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
def process_taxid_csv(
    csv_path, gene_name, sequence_type, processor, output_manager, save_genbank=False, header_format="basic"
):
    try:
        samples_csv = Path(csv_path)
        logger.debug(f"Samples file: {samples_csv}")

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

            # Initialise progress tracking
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
                    
                    # Fetch taxon name for this taxid
                    input_taxa = ""
                    try:
                        taxonomy, taxon_ranks, initial_rank, taxon_ids = processor.entrez.fetch_taxonomy(taxid)
                        if taxonomy:
                            # Get the most specific taxon name (last item in taxonomy list)
                            input_taxa = taxonomy[-1] if taxonomy else ""
                    except Exception as e:
                        logger.warning(f"Could not fetch taxon name for taxid {taxid}: {e}")
                    
                    process_sample(
                        process_id=process_id,
                        taxid=taxid,
                        sequence_type=sequence_type,
                        processor=processor,
                        output_manager=output_manager,
                        gene_name=gene_name,
                        save_genbank=save_genbank,
                        header_format=header_format,
                        input_taxa=input_taxa,
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
    csv_path,
    gene_name,
    sequence_type,
    processor,
    output_manager,
    entrez,
    save_genbank=False,
    header_format="basic",
):
    try:
        taxonomy_csv = Path(csv_path)
        logger.info(f"Samples file: {taxonomy_csv}")

        # Read in csv file
        with open(taxonomy_csv, "r", newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)

            # Check for required ID column
            if "id" not in [field.lower() for field in reader.fieldnames]:
                logger.error("Missing required 'ID' column in taxonomy CSV")
                sys.exit(1)

            # Check for at least one taxonomic rank column
            taxonomic_columns = ["phylum", "class", "order", "family", "genus", "species"]
            available_taxonomic_columns = [
                col for col in taxonomic_columns 
                if col.lower() in [field.lower() for field in reader.fieldnames]
            ]

            if not available_taxonomic_columns:
                logger.error(
                    f"CSV must contain at least one taxonomic rank column from: {taxonomic_columns}"
                )
                logger.error(
                    f"Available columns: {list(reader.fieldnames)}"
                )
                sys.exit(1)

            logger.info(f"Found taxonomic columns: {available_taxonomic_columns}")

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

                    # Validate we have at least some taxonomic information
                    if not any([phylum, class_name, order, family, genus, species]):
                        logger.warning(f"No taxonomic information provided for {process_id}"
                        )
                        output_manager.log_failure(
                            process_id, 
                            "unknown", 
                            "No taxonomic information provided")
                        continue

                    # Fetch taxid from taxonomic information
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

                    # Construct input_taxa from most specific taxonomic rank
                    input_taxa = species if species else (genus if genus else (family if family else (order if order else (class_name if class_name else phylum))))

                    logger.info(f"Starting sequence search:")
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
                        header_format=header_format,
                        input_taxa=input_taxa,
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
