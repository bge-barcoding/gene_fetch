# tests/conftest.py

import pytest
from unittest.mock import patch, MagicMock
from io import StringIO

@pytest.fixture
def mock_ncbi_services():
    """Mock all NCBI API calls for tests."""
    
    # Sample GenBank format for a protein record
    protein_record = """LOCUS       P00395                 513 aa            linear   MAM 29-SEP-2021
DEFINITION  Cytochrome c oxidase subunit 1.
ACCESSION   P00395
VERSION     P00395.2
DBSOURCE    UniProtKB: locus COX1_HUMAN, accession P00395;
KEYWORDS    Copper; Electron transport; Heme; Iron; Metal-binding;
            Mitochondrion; Oxidoreductase; Transmembrane; Transport.
SOURCE      Homo sapiens (human)
...
"""
    
    # Sample GenBank format for a nucleotide record
    nucleotide_record = """LOCUS       NC_012920            16569 bp    DNA     circular PRI 06-APR-2023
DEFINITION  Homo sapiens mitochondrion, complete genome.
ACCESSION   NC_012920
VERSION     NC_012920.1
DBLINK      BioProject: PRJNA30353
KEYWORDS    RefSeq.
SOURCE      mitochondrion Homo sapiens (human)
...
"""
    
    # Sample search results in XML format
    search_results = {
        "Count": "5",
        "RetMax": "5",
        "RetStart": "0",
        "IdList": ["16754380", "16754378", "16754376", "16754374", "16754372"],
        "TranslationSet": [],
        "TranslationStack": [
            {"Term": "cox1[Gene]", "Field": "Gene", "Count": "16754", "Explode": "N"},
            {"Term": "Homo sapiens[Organism]", "Field": "Organism", "Count": "21550362", "Explode": "Y"},
            "AND"
        ],
        "QueryTranslation": "cox1[Gene] AND Homo sapiens[Organism]"
    }
    
    # Sample taxonomy information
    taxonomy_info = [
        {
            "ScientificName": "Homo sapiens",
            "Rank": "species",
            "TaxId": "9606",
            "Division": "Primates",
            "LineageEx": [
                {"TaxId": "131567", "ScientificName": "cellular organisms", "Rank": "no rank"},
                {"TaxId": "2759", "ScientificName": "Eukaryota", "Rank": "superkingdom"},
                {"TaxId": "33208", "ScientificName": "Metazoa", "Rank": "kingdom"},
                {"TaxId": "7711", "ScientificName": "Chordata", "Rank": "phylum"},
                {"TaxId": "40674", "ScientificName": "Mammalia", "Rank": "class"},
                {"TaxId": "9443", "ScientificName": "Primates", "Rank": "order"},
                {"TaxId": "9604", "ScientificName": "Hominidae", "Rank": "family"},
                {"TaxId": "9605", "ScientificName": "Homo", "Rank": "genus"}
            ]
        }
    ]
    
    # Create mock functions that return these sample data
    def mock_fetch(*args, **kwargs):
        db = kwargs.get('db', args[0] if args else None)
        if db == 'protein':
            return StringIO(protein_record)
        elif db == 'nucleotide':
            return StringIO(nucleotide_record)
        elif db == 'taxonomy':
            # Return taxonomy info in a format expected by the code
            return MagicMock(read=lambda: str(taxonomy_info))
        return StringIO("")
    
    def mock_search(*args, **kwargs):
        # Return the search results as is
        return search_results
    
    # Apply the patches for all Entrez functions
    with patch('gene_fetch.entrez_handler.Entrez.efetch', side_effect=mock_fetch) as mock_efetch, \
         patch('gene_fetch.entrez_handler.Entrez.esearch', side_effect=mock_search) as mock_esearch, \
         patch('gene_fetch.entrez_handler.Entrez.read', return_value=search_results) as mock_read:
            
            # Yield to allow the test to run
            yield {
                'mock_efetch': mock_efetch,
                'mock_esearch': mock_esearch,
                'mock_read': mock_read
            }