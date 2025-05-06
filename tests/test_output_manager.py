import pytest
import tempfile
import csv
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from gene_fetch.output_manager import OutputManager
from gene_fetch.entrez_handler import EntrezHandler
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


@pytest.fixture
def output_manager():
    """Create OutputManager instance for testing with temp dir."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create OutputManager with the temporary directory
        output_dir = Path(tmpdirname)
        manager = OutputManager(output_dir)
        yield manager


def test_initialization(output_manager):
    """OutputManager initialises correctly with required dirs and files."""
    # Check if dirs created
    assert output_manager.output_dir.exists()
    assert output_manager.nucleotide_dir.exists()
    assert output_manager.genbank_dir.exists()
    
    # Check if files created with correct headers
    assert output_manager.failed_searches_path.exists()
    assert output_manager.sequence_refs_path.exists()
    
    # Verify headers in failed_searches.csv
    with open(output_manager.failed_searches_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header == ["process_id", "taxid", "error_type", "timestamp"]
    
    # Verify headers in sequence_references.csv
    with open(output_manager.sequence_refs_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header == [
            "process_id", "taxid", "protein_accession", "protein_length", 
            "nucleotide_accession", "nucleotide_length", "matched_rank", 
            "ncbi_taxonomy", "reference_name", "protein_reference_path", 
            "nucleotide_reference_path"
        ]


def test_log_failure(output_manager):
    """Fails correctly logged to the failed_searches.csv file."""
    # Log test failure
    output_manager.log_failure("test_process", "9606", "test_error")
    
    # Verify failure was logged
    with open(output_manager.failed_searches_path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        row = next(reader)
        assert row[0] == "test_process"
        assert row[1] == "9606"
        assert row[2] == "test_error"


def test_write_sequence_reference(output_manager):
    """Write sequence metadata to the sequence_references.csv file."""
    # Create test data
    test_data = {
        "process_id": "test_process",
        "taxid": "9606",
        "protein_accession": "P12345",
        "protein_length": "500",
        "nucleotide_accession": "NM_12345",
        "nucleotide_length": "1500",
        "matched_rank": "species",
        "taxonomy": "Homo sapiens",
        "protein_path": "protein/P12345.fasta",
        "nucleotide_path": "nucleotide/NM_12345.fasta"
    }
    
    # Write reference
    output_manager.write_sequence_reference(test_data)
    
    # Verify written correctly
    with open(output_manager.sequence_refs_path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        row = next(reader)
        assert row[0] == "test_process"
        assert row[1] == "9606"
        assert row[2] == "P12345"
        assert row[3] == "500"
        assert row[4] == "NM_12345"
        assert row[5] == "1500"
        assert row[6] == "species"
        assert row[7] == "Homo sapiens"
        assert row[8] == "test_process"
        assert row[9] == "protein/P12345.fasta"
        assert row[10] == "nucleotide/NM_12345.fasta"


def test_save_sequence_summary(output_manager):
    """Save summary of sequences to a CSV file."""
    # Create test sequence records
    sequences = [
        SeqRecord(Seq("ACGT" * 25), id="SEQ1", description="Test sequence 1"),
        SeqRecord(Seq("TGCA" * 50), id="SEQ2", description="Test sequence 2")
    ]
    
    # Save summary
    output_manager.save_sequence_summary(sequences, "test")
    
    # Verify summary file was created
    summary_path = output_manager.output_dir / "fetched_test_sequences.csv"
    assert summary_path.exists()
    
    with open(summary_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header == ["Accession", "Length", "Description"]
        
        row1 = next(reader)
        assert row1[0] == "SEQ1"
        assert row1[1] == "100"  # 4 * 25 = 100
        assert row1[2] == "Test sequence 1"
        
        row2 = next(reader)
        assert row2[0] == "SEQ2"
        assert row2[1] == "200"  # 4 * 50 = 200
        assert row2[2] == "Test sequence 2"


@patch('gene_fetch.output_manager.EntrezHandler')
def test_save_genbank_file(mock_entrez, output_manager):
    """Save GenBank file fetched from NCBI."""
    # Create mock handle with test content
    mock_handle = MagicMock()
    mock_handle.read.return_value = "MOCK GENBANK CONTENT"
    
    # Mock EntrezHandler's fetch method
    mock_entrez_instance = MagicMock()
    mock_entrez_instance.fetch.return_value = mock_handle
    
    # Define test parameters
    record_id = "NM_12345"
    db = "nucleotide"
    output_path = output_manager.genbank_dir / "test_genbank.gb"
    
    # Call method
    result = output_manager.save_genbank_file(mock_entrez_instance, record_id, db, output_path)
    
    # Verify result and file
    assert result is True
    assert output_path.exists()
    
    # Check file content
    with open(output_path, 'r') as f:
        content = f.read()
        assert content == "MOCK GENBANK CONTENT"
    
    # Verify mock was called
    mock_entrez_instance.fetch.assert_called_once_with(
        db=db, id=record_id, rettype="gb", retmode="text"
    )


@patch('gene_fetch.output_manager.EntrezHandler')
def test_save_genbank_file_error(mock_entrez, output_manager):
    """Handle of errors when saving a GenBank file."""
    # Mock EntrezHandler's fetch method to raise an exception
    mock_entrez_instance = MagicMock()
    mock_entrez_instance.fetch.side_effect = Exception("Test error")
    
    # Define test parameters
    record_id = "NM_12345"
    db = "nucleotide"
    output_path = output_manager.genbank_dir / "test_error.gb"
    
    # Call method
    result = output_manager.save_genbank_file(mock_entrez_instance, record_id, db, output_path)
    
    # Verify result
    assert result is False
    assert not output_path.exists()


def test_standalone_save_genbank_file():
    """Test standalone save_genbank_file function."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Set up test parameters
        output_dir = Path(tmpdirname)
        output_path = output_dir / "test_genbank.gb"
        
        # Mock EntrezHandler
        mock_entrez = MagicMock()
        mock_handle = MagicMock()
        mock_handle.read.return_value = "STANDALONE MOCK CONTENT"
        mock_entrez.fetch.return_value = mock_handle
        
        # Import standalone function
        from gene_fetch.output_manager import save_genbank_file
        
        # Call function
        with patch('gene_fetch.output_manager.OutputManager') as mock_output_manager:
            # Set up mock OutputManager
            mock_instance = MagicMock()
            mock_instance.save_genbank_file.return_value = True
            mock_output_manager.return_value = mock_instance
            
            # Call standalone function
            result = save_genbank_file(mock_entrez, "NM_12345", "nucleotide", output_path)
            
            # Verify result
            assert result is True
            mock_instance.save_genbank_file.assert_called_once_with(
                mock_entrez, "NM_12345", "nucleotide", output_path
            )