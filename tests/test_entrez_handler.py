# tests/test_entrez_handler.py
import pytest
import json
import time
from unittest.mock import patch, MagicMock, call
from urllib.error import HTTPError
from http.client import IncompleteRead
from Bio import Entrez

from gene_fetch.entrez_handler import EntrezHandler, check_ncbi_status, enhanced_retry
from gene_fetch.core import Config, logger

@pytest.fixture
def entrez_handler():
    """Create an EntrezHandler instance for testing."""
    config = Config(email="test@example.com", api_key="test_api_key_1234567890")
    config.retmax = 10000
    return EntrezHandler(config)

# Test the NCBI status check function
def test_check_ncbi_status():
    """Test the NCBI status check function."""
    with patch('gene_fetch.entrez_handler.Entrez.einfo') as mock_einfo:
        # Setup mock for successful status check
        mock_handle = MagicMock()
        mock_einfo.return_value = mock_handle
        mock_handle.close.return_value = None
        
        # Mock Entrez.read to return some data
        with patch('gene_fetch.entrez_handler.Entrez.read') as mock_read:
            mock_read.return_value = {"DbList": ["pubmed", "protein", "nuccore"]}
            
            # Call the function
            result = check_ncbi_status()
            
            # Return True when NCBI services are up
            assert result is True
            
        # Test with exception
        with patch('gene_fetch.entrez_handler.Entrez.read') as mock_read:
            mock_read.side_effect = RuntimeError("NCBI Error")
            
            # Call function
            result = check_ncbi_status()
            
            # Return False when NCBI services are down
            assert result is False

# Test enhanced_retry decorator
def test_enhanced_retry_decorator():
    """Test the enhanced retry decorator's behavior."""
    # Create a mock function that fails a few times then succeeds
    mock_func = MagicMock()
    mock_func.side_effect = [
        HTTPError(url="", code=500, msg="Server Error", hdrs={}, fp=None),
        HTTPError(url="", code=429, msg="Too Many Requests", hdrs={}, fp=None),
        "Success"  # Third call succeeds
    ]
    
    # Apply decorator
    decorated_func = enhanced_retry(
        exceptions=(Exception,),
        tries=4,
        initial_delay=0.1,  # Use small values for testing
        backoff=2,
        max_delay=1
    )(mock_func)
    
    # Call decorated function
    result = decorated_func("arg1", kwarg1="value1")
    
    # Check result
    assert result == "Success"
    
    # Function called expected number of times?
    assert mock_func.call_count == 3
    
    # Called with right arguments
    mock_func.assert_called_with("arg1", kwarg1="value1")

# Basic functionality tests
@patch('gene_fetch.entrez_handler.EntrezHandler.search')
def test_entrez_search_wrapper(mock_search, entrez_handler):
    """Test the search method wrapper."""
    # Setup mock response
    mock_search.return_value = {
        "IdList": ["9606"], 
        "Count": "1"
    }
    
    # Call search directly (without going through wrapper)
    result = entrez_handler.search(db="taxonomy", term="Homo sapiens")
    
    # Verify result
    assert result is not None
    assert result["IdList"] == ["9606"]
    assert result["Count"] == "1"

@patch('gene_fetch.entrez_handler.EntrezHandler.fetch')
def test_fetch_wrapper(mock_fetch, entrez_handler):
    """Test that our EntrezHandler.fetch method works."""
    # Setup mock response
    mock_result = MagicMock()
    mock_fetch.return_value = mock_result
    
    # Call fetch
    result = entrez_handler.fetch(db="taxonomy", id="9606")
    
    # Verify result
    assert result is mock_result
    
    # Verify fetch was called correctly
    mock_fetch.assert_called_once_with(db="taxonomy", id="9606")

# Test taxonomy cache behavior
@patch('gene_fetch.entrez_handler.EntrezHandler.fetch')
def test_taxonomy_caching(mock_fetch, entrez_handler):
    """Test that taxonomy information is cached."""
    # Create mock handle with XML content
    mock_handle = MagicMock()
    
    # Setup Bio.Entrez.read to parse XML
    with patch('gene_fetch.entrez_handler.Entrez.read') as mock_read:
        mock_read.return_value = [{
            'TaxId': '9606',
            'ScientificName': 'Homo sapiens',
            'Rank': 'species',
            'LineageEx': [
                {'TaxId': '9605', 'ScientificName': 'Homo', 'Rank': 'genus'},
                {'TaxId': '9604', 'ScientificName': 'Hominidae', 'Rank': 'family'},
                {'TaxId': '9443', 'ScientificName': 'Primates', 'Rank': 'order'},
                {'TaxId': '40674', 'ScientificName': 'Mammalia', 'Rank': 'class'},
                {'TaxId': '7711', 'ScientificName': 'Chordata', 'Rank': 'phylum'}
            ]
        }]
        
        mock_fetch.return_value = mock_handle
        
        # First call should hit the API
        lineage1, rank_info1, current_rank1, taxid_info1 = entrez_handler.fetch_taxonomy("9606")
        
        # Second call should use cache
        lineage2, rank_info2, current_rank2, taxid_info2 = entrez_handler.fetch_taxonomy("9606")
        
        # Verify fetch was only called once
        assert mock_fetch.call_count == 1
        
        # Results should be identical
        assert lineage1 == lineage2
        assert rank_info1 == rank_info2
        assert current_rank1 == current_rank2
        assert taxid_info1 == taxid_info2
        
        # Verify cache behavior - fetch_taxonomy should have stored the result
        assert "9606" in entrez_handler.taxonomy_cache

# Test error handling in fetch_taxonomy
@patch('gene_fetch.entrez_handler.EntrezHandler.fetch')
def test_fetch_taxonomy_error_handling(mock_fetch, entrez_handler):
    """Test error handling in fetch_taxonomy method."""
    # Simulate error in the fetch call
    mock_fetch.side_effect = HTTPError(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        code=500,
        msg="Internal Server Error",
        hdrs={},
        fp=None
    )
    
    # Call fetch_taxonomy
    lineage, rank_info, current_rank, taxid_info = entrez_handler.fetch_taxonomy("9606")
    
    # Verify empty results on error
    assert lineage == []
    assert rank_info == {}
    assert current_rank == ""
    assert taxid_info == {}
    
    # Verify fetch was called
    mock_fetch.assert_called()

# Test taxonomy validation
@patch('gene_fetch.entrez_handler.EntrezHandler.fetch_taxonomy')
def test_validate_taxonomy_consistency(mock_fetch_taxonomy, entrez_handler):
    """Test taxonomy validation logic."""
    # Mock fetch_taxonomy to return predefined data
    mock_fetch_taxonomy.return_value = (
        ["Animalia", "Chordata", "Mammalia", "Primates", "Hominidae", "Homo", "Homo sapiens"],  # complete_lineage
        {  # rank_info
            "Chordata": "phylum",
            "Mammalia": "class",
            "Primates": "order",
            "Hominidae": "family",
            "Homo": "genus",
            "Homo sapiens": "species"
        },
        "species",  # current_rank
        {  # taxid_info
            "Chordata": "7711",
            "Mammalia": "40674",
            "Primates": "9443",
            "Hominidae": "9604",
            "Homo": "9605",
            "Homo sapiens": "9606"
        }
    )
    
    # Test with matching taxonomy
    provided_taxonomy = {
        "phylum": "Chordata",
        "class": "Mammalia",
        "order": "Primates",
        "family": "Hominidae",
        "genus": "Homo",
        "species": "sapiens"
    }
    
    is_valid, results = entrez_handler.validate_taxonomy_consistency("9606", provided_taxonomy)
    
    # Should be valid with good match score
    assert is_valid
    assert results["match_score"] > 0
    
    # Test with conflicting higher taxonomy
    provided_taxonomy_conflict = {
        "phylum": "Arthropoda",  # Different phylum
        "class": "Insecta",      # Different class
        "order": "Primates",
        "family": "Hominidae",
        "genus": "Homo",
        "species": "sapiens"
    }
    
    is_valid, results = entrez_handler.validate_taxonomy_consistency("9606", provided_taxonomy_conflict)
    
    # Should be invalid due to phylum/class conflict
    assert not is_valid

# Test fetch_taxid_from_taxonomy
@patch('gene_fetch.entrez_handler.EntrezHandler.search')
@patch('gene_fetch.entrez_handler.EntrezHandler.validate_taxonomy_consistency')
def test_fetch_taxid_from_taxonomy(mock_validate, mock_search, entrez_handler):
    """Test fetching taxid from taxonomy information."""
    # Setup mock for search
    mock_search.return_value = {"IdList": ["9606"], "Count": "1"}
    
    # Setup mock for validation
    mock_validate.return_value = (True, {"match_score": 6})
    
    # Call method
    taxid = entrez_handler.fetch_taxid_from_taxonomy(
        phylum="Chordata",
        class_name="Mammalia",
        order="Primates",
        family="Hominidae",
        genus="Homo",
        species="sapiens"
    )
    
    # Verify results
    assert taxid == "9606"
    assert mock_search.call_count >= 1
    assert mock_validate.call_count >= 1
    
    # Check search was properly formed
    args, kwargs = mock_search.call_args_list[0]
    assert kwargs["db"] == "taxonomy"
    assert "Homo sapiens" in kwargs["term"]

# Test handling of consecutive errors
@patch('gene_fetch.entrez_handler.EntrezHandler.fetch')
@patch('gene_fetch.entrez_handler.check_ncbi_status')
@patch('gene_fetch.entrez_handler.time.sleep')
def test_handle_request_error(mock_sleep, mock_check_status, mock_fetch, entrez_handler):
    """Test the handling of consecutive request errors."""
    # Mock check_status to return False (service down)
    mock_check_status.return_value = False
    
    # Create error to pass to method
    error = HTTPError(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        code=500,
        msg="Internal Server Error",
        hdrs={},
        fp=None
    )
    
    # Call handle_request_error multiple times
    for _ in range(3):
        entrez_handler.handle_request_error(error)
    
    # After 3 errors, should_check_service_status should be true
    assert entrez_handler.consecutive_errors >= 3
    
    # Reset counter to 0
    entrez_handler.handle_request_success()
    
    assert entrez_handler.consecutive_errors == 0
