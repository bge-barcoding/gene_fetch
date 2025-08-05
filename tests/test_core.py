# tests/test_core.py
import pytest
import tempfile
import logging
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

from gene_fetch.core import (
    Config, 
    make_out_dir, 
    setup_logging, 
    log_progress, 
    get_process_id_column
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


def test_make_out_dir(temp_dir):
    """Test making output directories."""
    # Test creating a new directory
    test_dir = temp_dir / "test_dir"
    make_out_dir(test_dir)
    assert test_dir.exists()
    assert test_dir.is_dir()
    
    # Test with nested directories
    nested_dir = temp_dir / "parent" / "child" / "grandchild"
    make_out_dir(nested_dir)
    assert nested_dir.exists()
    assert nested_dir.is_dir()
    
    # Test with existing directory (should not raise an exception)
    make_out_dir(test_dir)  # This should not raise an exception


def test_setup_logging(temp_dir):
    """Test setting up logging."""
    # Setup logging
    logger = setup_logging(temp_dir)
    
    # Check that the log file was created
    log_file = temp_dir / "gene_fetch.log"
    assert log_file.exists()
    
    # Check logger settings
    assert logger.level == logging.INFO
    assert len(logger.handlers) == 2  # Console and file handlers
    
    # Test logging something and check if it appears in the log file
    test_message = "Test log message"
    logger.info(test_message)
    
    with open(log_file, "r") as f:
        log_content = f.read()
        assert test_message in log_content


@patch('gene_fetch.core.logger')
def test_log_progress(mock_logger):
    """Test logging progress."""
    # Test initial progress
    log_progress(0, 100)
    mock_logger.info.assert_any_call("")
    mock_logger.info.assert_any_call("Starting processing: 0/100 samples processed (0%)")
    
    # Test intermediate progress (multiple of interval)
    mock_logger.reset_mock()
    log_progress(10, 100)
    mock_logger.info.assert_called_once()
    # Extract the call arguments for verification
    call_args = mock_logger.info.call_args[0][0]
    assert "Progress: 10/100" in call_args
    assert "10.00%" in call_args  # Changed from "10%" to "10.00%"
    
    # Test intermediate progress (not multiple of interval)
    mock_logger.reset_mock()
    log_progress(11, 100)
    mock_logger.info.assert_not_called()
    
    # Test custom interval
    mock_logger.reset_mock()
    log_progress(5, 100, interval=5)
    mock_logger.info.assert_called_once()
    
    # Test final progress
    mock_logger.reset_mock()
    log_progress(100, 100)
    mock_logger.info.assert_called_once_with("Processed: 100/100 samples processed (100%)")


@patch('gene_fetch.core.logger')
def test_get_process_id_column(mock_logger):
    """Test identifying the process ID column."""
    # Test exact match
    header = ["ID", "taxid", "genus", "species"]
    assert get_process_id_column(header) == "ID"
    
    # Test with whitespace
    header = [" ID ", "taxid", "genus", "species"]
    assert get_process_id_column(header) == " ID "
    
    # Test case-insensitive match
    header = ["id", "taxid", "genus", "species"]
    assert get_process_id_column(header) == "id"
    
    # Test with alternative names
    header = ["process_id", "taxid", "genus", "species"]
    assert get_process_id_column(header) == "process_id"
    
    header = ["Process ID", "taxid", "genus", "species"]
    assert get_process_id_column(header) == "Process ID"
    
    header = ["sample", "taxid", "genus", "species"]
    assert get_process_id_column(header) == "sample"
    
    # Test when no match is found
    header = ["no_match", "taxid", "genus", "species"]
    assert get_process_id_column(header) is None
    mock_logger.error.assert_called_once()


def test_config_initialisation():
    """Test Config class initialisation with valid parameters."""
    # Test valid initialisation
    email = "test@example.com"
    api_key = "valid_test_key_12345"
    
    config = Config(email=email, api_key=api_key)
    
    # Check attributes
    assert config.email == email
    assert config.api_key == api_key
    assert config.max_calls_per_second == 10
    assert config.fetch_batch_size == 200
    assert config.batch_delay == (1, 2)
    assert "protein" in config.valid_sequence_types
    assert "nucleotide" in config.valid_sequence_types
    assert "both" in config.valid_sequence_types


def test_config_initialisation_missing_email():
    """Test Config class initialisation with missing email."""
    with pytest.raises(ValueError) as excinfo:
        Config(email="", api_key="valid_test_key_12345")
    assert "Email address required" in str(excinfo.value)


def test_config_initialisation_missing_api_key():
    """Test Config class initialisation with missing API key."""
    with pytest.raises(ValueError) as excinfo:
        Config(email="test@example.com", api_key="")
    assert "API key required" in str(excinfo.value)


def test_config_update_thresholds():
    """Test updating sequence length thresholds."""
    config = Config(email="test@example.com", api_key="valid_test_key_12345")
    
    # Initial values
    assert config.protein_length_threshold == 500
    assert config.nucleotide_length_threshold == 1000
    
    # Update thresholds
    config.update_thresholds(protein_size=200, nucleotide_size=500)
    
    # Check updated values
    assert config.protein_length_threshold == 200
    assert config.nucleotide_length_threshold == 500


def test_config_set_gene_search_term_rRNA():
    """Test setting gene search term for rRNA genes."""
    config = Config(email="test@example.com", api_key="valid_test_key_12345")
    
    # Set for 16S rRNA gene
    canonical_gene_name, search_type = config.set_gene_search_term("16s")
    
    # Check search term and type
    assert "16S ribosomal RNA[Gene]" in config.gene_search_term
    assert "16S rRNA[Gene]" in config.gene_search_term
    assert canonical_gene_name == "16s"
    assert search_type == "rRNA"


def test_config_set_gene_search_term_protein_coding():
    """Test setting gene search term for protein-coding genes."""
    config = Config(email="test@example.com", api_key="valid_test_key_12345")
    
    # Set for rbcL gene
    canonical_gene_name, search_type = config.set_gene_search_term("rbcl")
    
    # Check search term and type
    assert "rbcL[Gene]" in config.gene_search_term
    assert "RBCL[Gene]" in config.gene_search_term
    assert canonical_gene_name == "rbcl"
    assert search_type == "protein-coding"


def test_config_set_gene_search_term_generic():
    """Test setting gene search term for generic genes."""
    config = Config(email="test@example.com", api_key="valid_test_key_12345")
    
    # Set for a generic gene
    canonical_gene_name, search_type = config.set_gene_search_term("unknown_gene")
    
    # Check search term and type
    assert "unknown_gene[Title]" in config.gene_search_term
    assert "unknown_gene[Gene]" in config.gene_search_term
    assert '"unknown_gene"[Protein Name]' in config.gene_search_term  # Fixed: was [Text Word]
    assert canonical_gene_name == "unknown_gene"
    assert search_type == "generic"


def test_config_set_gene_search_term_alias_mapping():
    """Test that aliases are correctly mapped to canonical names."""
    config = Config(email="test@example.com", api_key="valid_test_key_12345")

    # Test COI -> cox1 mapping
    canonical_gene_name, search_type = config.set_gene_search_term("coi")
    assert canonical_gene_name == "cox1"
    assert search_type == "protein-coding"
    assert "cox1[Gene]" in config.gene_search_term

    # Test uppercase COI -> cox1 mapping
    canonical_gene_name, search_type = config.set_gene_search_term("COI")
    assert canonical_gene_name == "cox1"
    assert search_type == "protein-coding"

    # Test rbcL -> rbcl mapping
    canonical_gene_name, search_type = config.set_gene_search_term("rbcL")
    assert canonical_gene_name == "rbcl" 
    assert search_type == "protein-coding"

    # Test MATK -> matk mapping
    canonical_gene_name, search_type = config.set_gene_search_term("MATK")
    assert canonical_gene_name == "matk"
    assert search_type == "protein-coding"

    # Test that canonical names still work
    canonical_gene_name, search_type = config.set_gene_search_term("cox1")
    assert canonical_gene_name == "cox1"
    assert search_type == "protein-coding"


def test_config_validate_credentials_valid_email_and_api_key():
    """Test credential validation with valid email and API key."""
    # These should not raise exceptions
    Config.validate_credentials("test@example.com", "valid_test_key_12345")
    Config.validate_credentials("user@domain.co.uk", "test_api_key_1234567890")
    Config.validate_credentials("real.user@ncbi.gov", "1234567890abcdef")


def test_config_validate_credentials_invalid_email():
    """Test credential validation with invalid email formats."""
    with pytest.raises(ValueError, match="Invalid email format"):
        Config.validate_credentials("invalid-email", "valid_test_key_12345")
    
    with pytest.raises(ValueError, match="Invalid email format"):
        Config.validate_credentials("missing@domain", "valid_test_key_12345")
    
    with pytest.raises(ValueError, match="Invalid email format"):
        Config.validate_credentials("@domain.com", "valid_test_key_12345")


def test_config_validate_credentials_invalid_api_key():
    """Test credential validation with invalid API keys."""
    # Test explicitly fake keys
    with pytest.raises(ValueError, match="Invalid API key"):
        Config.validate_credentials("test@example.com", "fake_key")
    
    with pytest.raises(ValueError, match="Invalid API key"):
        Config.validate_credentials("test@example.com", "fake_api_key")
    
    # Test short keys (not in test list)
    with pytest.raises(ValueError, match="Invalid API key"):
        Config.validate_credentials("test@example.com", "short")
    
    with pytest.raises(ValueError, match="Invalid API key"):
        Config.validate_credentials("test@example.com", "abc123")


def test_config_validate_credentials_test_keys_allowed():
    """Test that specific test API keys are allowed."""
    # These test keys should be allowed regardless of length rules
    Config.validate_credentials("test@example.com", "test_api_key_1234567890")
    Config.validate_credentials("test@example.com", "valid_test_key_12345")