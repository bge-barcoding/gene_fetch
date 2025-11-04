# src/gene_fetch/sequence_processor.py
"""
SequenceProcessor class for Gene Fetch.
Handles the processing and validation of fetched sequences.
"""

import re
import logging
from time import sleep
from random import uniform
from typing import Optional, Tuple, List, Dict, Any

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .core import logger, Config
from .entrez_handler import EntrezHandler, check_ncbi_status


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
            # For any other gene add reasonable variations
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
                
                if len(cds_record.seq) >= 100:
                    # Try to preserve description
                    if not cds_record.description or cds_record.description == cds_record.id:
                        logger.info(f"Description empty or just ID, attempting to restore from: '{record.description}'")
                        cds_record.description = record.description
                    
                    logger.info(f"Successfully extracted CDS of length {len(cds_record.seq)} (accession {record.id})")
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

                        # Update the progress_counters with actual total if in fetch_all mode
                        if fetch_all and progress_counters:
                            # If max_sequences specified, use min(max_sequences, len(id_list))
                            # Otherwise use the actual number of sequences found
                            total_sequences = (
                                min(max_sequences, len(id_list))
                                if max_sequences
                                else len(id_list)
                            )
                            progress_counters["total_sequences"] = total_sequences

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

                if fetch_all and self.config.nucleotide_length_threshold <= 0:
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
                                                                                        
                                            ##### DEBUGGING
                                            logger.info(f"CDS record description BEFORE appending: '{cds_record.description}'")

                                            if fetch_all:
                                                nucleotide_records.append(cds_record)
                                                ##### DEBUGGING
                                                logger.info(f"CDS record description AFTER appending: '{nucleotide_records[-1].description}'")
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
    ) -> Tuple[List[SeqRecord], List[SeqRecord], List[str], str, Optional[str], Optional[str]]:
        
        # Initialise empty lists for records
        protein_records = []
        nucleotide_records = []
        best_taxonomy = []
        best_matched_rank = None
        
        # Capture first_matched_taxid_rank
        first_matched_taxid = taxid  # The input taxid is always the first matched taxid
        first_matched_taxid_rank = None  # Will be set below

        # Fetch taxonomy first (from cache if available)
        logger.debug(f"Starting sequence search for {gene_name} using taxid {taxid}")
        taxonomy, taxon_ranks, initial_rank, taxon_ids = self.entrez.fetch_taxonomy(
            taxid
        )
        if not taxonomy:
            logger.error(
                f"Could not fetch taxonomy for taxID ({taxid}), cannot search for sequences"
            )
            return [], [], [], "No taxonomy found", None, None

        # Get ordered list of ranks to traverse
        current_taxonomy = taxonomy[:]
        current_taxon = current_taxonomy.pop()  # Start with species
        current_rank = taxon_ranks.get(current_taxon, "unknown")
        current_taxid = taxid
        
        # Set the first_matched_taxid_rank
        first_matched_taxid_rank = f"{current_rank}:{current_taxon}"

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
                return [], [], [], "No match", None, None
            logger.info(
                f"Single taxid mode: Found {len(protein_records)} protein and {len(nucleotide_records)} nucleotide sequences"
            )
            return (
                protein_records,
                nucleotide_records,
                best_taxonomy,
                matched_rank,
                first_matched_taxid,
                first_matched_taxid_rank,
            )
        else:
            # Batch mode: require both for 'both' type
            if sequence_type == "both" and (
                not protein_records or not nucleotide_records
            ):
                logger.warning(
                    "Failed to find both protein and corresponding nucleotide sequence"
                )
                return [], [], [], "No match", None, None
            elif sequence_type == "protein" and not protein_records:
                logger.warning("No protein sequence found")
                return [], [], [], "No match", None, None
            elif sequence_type == "nucleotide" and not nucleotide_records:
                logger.warning("No nucleotide sequence found")
                return [], [], [], "No match", None, None

        logger.info(f"Search completed! Matched at rank: {matched_rank}")
        return protein_records, nucleotide_records, best_taxonomy, matched_rank, first_matched_taxid, first_matched_taxid_rank

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