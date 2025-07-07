# src/gene_fetch/entrez_handler.py
"""
EntrezHandler class for Gene Fetch.
Manages NCBI Entrez API calls, caching, and error recovery.
"""

import time
import logging
from functools import wraps
from ratelimit import limits, sleep_and_retry
from time import sleep
from random import uniform
from urllib.error import HTTPError
from http.client import IncompleteRead
from Bio import Entrez
from typing import Optional, Tuple, List, Dict, Any

from .core import logger, Config

# Function to check NCBI service status
def check_ncbi_status():
    try:
        handle = Entrez.einfo()
        _ = Entrez.read(handle)
        handle.close()
        # If einfo can be successfully queried then services are up
        return True
    except Exception as e:
        logger.warning(f"NCBI service check failed: {str(e)}")
        return False

# Retry decorator with NCBI status checking and exponential retry delay
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
        self.taxonomy_cache = (
            {}
        )  # taxid -> (complete_lineage, rank_info, current_rank, taxid_info)

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
                                f"{attempt + 1}/{max_retries}. Retrying in {delay} seconds... "
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

            logger.info(
                f"Successfully retrieved taxonomy information from NCBI for taxID: {taxid}"
            )
            logger.debug(f"Lineage: {complete_lineage}")
            logger.debug(f"Rank info: {rank_info}")

            # Store taxonomy in cache
            self.taxonomy_cache[taxid] = (
                complete_lineage,
                rank_info,
                current_rank,
                taxid_info,
            )

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
            f" Family: {family}, Order: {order}, Class: {class_name}, Phylum: {phylum})"
        )

        # Store the provided taxonomy for validation
        provided_taxonomy = {
            "phylum": phylum.strip() if phylum else "",
            "class": class_name.strip() if class_name else "",
            "order": order.strip() if order else "",
            "family": family.strip() if family else "",
            "genus": genus.strip() if genus else "",
            "species": species.strip() if species else "",
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
                        logger.warning(
                            f"Multiple taxids found for species {full_species}: {taxids}"
                        )

                        # Validate fetched taxids against the input taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(
                                taxid, provided_taxonomy
                            )
                            if is_valid:
                                valid_taxids.append((taxid, lineage))

                        if not valid_taxids:
                            logger.warning(
                                f"None of the taxids for {full_species} match the provided higher taxonomy"
                            )
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(
                                f"Found single valid taxid ({taxid}) for species {full_species}"
                            )
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(
                                valid_taxids, key=lambda x: x[1]["match_score"]
                            )[0]
                            logger.warning(
                                f"Multiple valid taxids for {full_species}, using best match: {best_taxid}"
                            )
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(
                            taxid, provided_taxonomy
                        )
                        if is_valid:
                            logger.info(
                                f"Found taxid ({taxid}) for species {full_species}"
                            )
                            return taxid
                        else:
                            logger.warning(
                                f"Taxid {taxid} for {full_species} does not match the provided higher taxonomy"
                            )
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
                        logger.warning(
                            f"Multiple taxids found for genus {genus}: {taxids}"
                        )

                        # Validate each taxid against the provided taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(
                                taxid, provided_taxonomy
                            )
                            if is_valid:
                                valid_taxids.append((taxid, lineage))

                        if not valid_taxids:
                            logger.warning(
                                f"None of the taxids for genus {genus} match the provided higher taxonomy"
                            )
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(
                                f"Found single valid taxid ({taxid}) for genus {genus}"
                            )
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(
                                valid_taxids, key=lambda x: x[1]["match_score"]
                            )[0]
                            logger.warning(
                                f"Multiple valid taxids for genus {genus}, using best match: {best_taxid}"
                            )
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(
                            taxid, provided_taxonomy
                        )
                        if is_valid:
                            logger.info(
                                f"Confirmed taxid {taxid} for {genus} is valid given input taxonomy"
                            )
                            return taxid
                        else:
                            logger.warning(
                                f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy"
                            )
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
                        logger.warning(
                            f"Multiple taxids found for family {family}: {taxids}"
                        )

                        # Validate each taxid against the provided taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(
                                taxid, provided_taxonomy
                            )
                            if is_valid:
                                valid_taxids.append((taxid, lineage))

                        if not valid_taxids:
                            logger.warning(
                                f"None of the taxids for family {family} match the provided higher taxonomy"
                            )
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(
                                f"Found single valid taxid ({taxid}) for family {family}"
                            )
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(
                                valid_taxids, key=lambda x: x[1]["match_score"]
                            )[0]
                            logger.warning(
                                f"Multiple valid taxids for family {family}, using best match: {best_taxid}"
                            )
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(
                            taxid, provided_taxonomy
                        )
                        if is_valid:
                            logger.info(
                                f"Confirmed taxid {taxid} for {family} is valid given input taxonomy"
                            )
                            return taxid
                        else:
                            logger.warning(
                                f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy"
                            )
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
                        logger.warning(
                            f"Multiple taxids found for order {order}: {taxids}"
                        )

                        # Validate each taxid against the provided taxonomy
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(
                                taxid, provided_taxonomy
                            )
                            if is_valid:
                                valid_taxids.append((taxid, lineage))

                        if not valid_taxids:
                            logger.warning(
                                f"None of the taxids for order {order} match the provided higher taxonomy"
                            )
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(
                                f"Found single valid taxid ({taxid}) for order {order}"
                            )
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(
                                valid_taxids, key=lambda x: x[1]["match_score"]
                            )[0]
                            logger.warning(
                                f"Multiple valid taxids for order {order}, using best match: {best_taxid}"
                            )
                            return best_taxid
                    else:
                        # Single match - still validate
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(
                            taxid, provided_taxonomy
                        )
                        if is_valid:
                            logger.info(
                                f"Confirmed taxid {taxid} for {order} is valid given input taxonomy"
                            )
                            return taxid
                        else:
                            logger.warning(
                                f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy"
                            )
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
                        logger.warning(
                            f"Multiple taxids found for class {class_name}: {taxids}"
                        )

                        # At class level, we validate against phylum only
                        valid_taxids = []
                        for taxid in taxids:
                            is_valid, lineage = self.validate_taxonomy_consistency(
                                taxid, provided_taxonomy
                            )
                            if is_valid:
                                valid_taxids.append((taxid, lineage))

                        if not valid_taxids:
                            logger.warning(
                                f"None of the taxids for class {class_name} match the provided higher taxonomy"
                            )
                            # Continue to next level
                        elif len(valid_taxids) == 1:
                            taxid = valid_taxids[0][0]
                            logger.info(
                                f"Found single valid taxid ({taxid}) for class {class_name}"
                            )
                            return taxid
                        else:
                            # Multiple valid taxids - take the one with most matching higher taxonomy
                            best_taxid = max(
                                valid_taxids, key=lambda x: x[1]["match_score"]
                            )[0]
                            logger.warning(
                                f"Multiple valid taxids for class {class_name}, using best match: {best_taxid}"
                            )
                            return best_taxid
                    else:
                        # Single match - still validate against phylum
                        taxid = taxids[0]
                        is_valid, lineage = self.validate_taxonomy_consistency(
                            taxid, provided_taxonomy
                        )
                        if is_valid:
                            logger.info(
                                f"Confirmed taxid {taxid} for {class_name} is valid given input taxonomy"
                            )
                            return taxid
                        else:
                            logger.warning(
                                f"Taxonomic mismatch for taxid {taxid}, it does not match the input higher taxonomy"
                            )
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
                        logger.warning(
                            f"Multiple taxids found for phylum {phylum}: {taxids}"
                        )
                        # Take the first one since we're already at phylum level
                        taxid = taxids[0]
                    else:
                        taxid = taxids[0]

                    logger.info(
                        f"Confirmed taxid {taxid} for {phylum} is valid given input taxonomy"
                    )
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
        logger.info(
            f"Comparing provided taxonomy against NCBI taxonomy for taxid {taxid}"
        )

        try:
            # Use the fetch_taxonomy method directly from EntrezHandler to get lineage information
            complete_lineage, rank_info, current_rank, taxid_info = self.fetch_taxonomy(
                taxid
            )

            if not complete_lineage:
                logger.warning(
                    f"No taxonomy data available for comparison (taxid: {taxid})"
                )
                return False, {"match_score": 0}

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

            logger.debug(
                f"Taxonomy comparison details - NCBI: {fetched_taxonomy}, Provided: {provided_taxonomy}"
            )

            # Check consistency between provided and fetched taxonomy
            # Start with higher taxonomic ranks which are less likely to have homonyms
            match_score = 0

            # Store details of what matched and what didn't
            match_details = {}

            # Check phylum
            if provided_taxonomy["phylum"] and "phylum" in fetched_taxonomy:
                phylum_match = (
                    provided_taxonomy["phylum"].lower()
                    == fetched_taxonomy["phylum"].lower()
                )
                match_details["phylum"] = "match" if phylum_match else "mismatch"
                if phylum_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(
                        f"Phylum mismatch for taxid {taxid}: "
                        f"provided '{provided_taxonomy['phylum']}' vs "
                        f"fetched '{fetched_taxonomy['phylum']}'"
                    )

            # Check class
            if provided_taxonomy["class"] and "class" in fetched_taxonomy:
                class_match = (
                    provided_taxonomy["class"].lower()
                    == fetched_taxonomy["class"].lower()
                )
                match_details["class"] = "match" if class_match else "mismatch"
                if class_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(
                        f"Class mismatch for taxid {taxid}: "
                        f"provided '{provided_taxonomy['class']}' vs "
                        f"fetched '{fetched_taxonomy['class']}'"
                    )

            # Check order
            if provided_taxonomy["order"] and "order" in fetched_taxonomy:
                order_match = (
                    provided_taxonomy["order"].lower()
                    == fetched_taxonomy["order"].lower()
                )
                match_details["order"] = "match" if order_match else "mismatch"
                if order_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(
                        f"Order mismatch for taxid {taxid}: "
                        f"provided '{provided_taxonomy['order']}' vs "
                        f"fetched '{fetched_taxonomy['order']}'"
                    )

            # Check family
            if provided_taxonomy["family"] and "family" in fetched_taxonomy:
                family_match = (
                    provided_taxonomy["family"].lower()
                    == fetched_taxonomy["family"].lower()
                )
                match_details["family"] = "match" if family_match else "mismatch"
                if family_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(
                        f"Family mismatch for taxid {taxid}: "
                        f"provided '{provided_taxonomy['family']}' vs "
                        f"fetched '{fetched_taxonomy['family']}'"
                    )

            # Check genus
            if provided_taxonomy["genus"] and "genus" in fetched_taxonomy:
                genus_match = (
                    provided_taxonomy["genus"].lower()
                    == fetched_taxonomy["genus"].lower()
                )
                match_details["genus"] = "match" if genus_match else "mismatch"
                if genus_match:
                    match_score += 1  # Changed to 1 point per match
                else:
                    logger.warning(
                        f"Genus mismatch for taxid {taxid}: "
                        f"provided '{provided_taxonomy['genus']}' vs "
                        f"fetched '{fetched_taxonomy['genus']}'"
                    )

            # Determine if this taxid is valid based on matches
            # We consider it valid if:
            # 1. No conflicts in higher taxonomy (phylum, class)
            # 2. Or if no higher taxonomy was provided to check against

            # Check for higher taxonomy mismatches that would invalidate the match
            higher_taxonomy_conflict = False

            if (
                provided_taxonomy["phylum"]
                and "phylum" in fetched_taxonomy
                and match_details.get("phylum") == "mismatch"
            ):
                higher_taxonomy_conflict = True

            if (
                provided_taxonomy["class"]
                and "class" in fetched_taxonomy
                and match_details.get("class") == "mismatch"
            ):
                higher_taxonomy_conflict = True

            is_valid = not higher_taxonomy_conflict

            # Log validation result
            if is_valid:
                logger.info(
                    f"Taxonomy validation: taxid {taxid} passed validation with match score: {match_score}"
                )
            else:
                logger.warning(
                    f"Taxonomy validation: taxid {taxid} failed validation due to higher taxonomy conflicts"
                )

            return is_valid, {
                "match_score": match_score,
                "details": match_details,
                "lineage": complete_lineage,
                "fetched_taxonomy": fetched_taxonomy,
            }

        except Exception as e:
            logger.error(f"Error validating taxonomy for taxid {taxid}: {e}")
            logger.error("Full error details:", exc_info=True)
            return False, {"match_score": 0}
