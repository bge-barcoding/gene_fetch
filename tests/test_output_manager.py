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
    """Create OutputManager instance for testing with temp dir (no GenBank)."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create OutputManager with the temporary directory
        output_dir = Path(tmpdirname)
        manager = OutputManager(output_dir, save_genbank=False, create_sequence_refs=True)
        yield manager


@pytest.fixture
def output_manager_with_genbank():
    """Create OutputManager instance for testing with temp dir (with GenBank)."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create OutputManager with the temporary directory and GenBank enabled
        output_dir = Path(tmpdirname)
        manager = OutputManager(output_dir, save_genbank=True, create_sequence_refs=True)
        yield manager


def test_initialisation(output_manager):
    """OutputManager initialises correctly with required dirs and files (no GenBank)."""
    # Check if dirs created
    assert output_manager.output_dir.exists()
    assert output_manager.protein_dir.exists()
    assert output_manager.nucleotide_dir.exists()
    
    # GenBank directories should not be created when save_genbank=False
    assert output_manager.genbank_dir is None
    assert output_manager.protein_genbank_dir is None
    assert output_manager.nucleotide_genbank_dir is None
    
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
            "process_id", "input_taxa", "first_matched_taxid", "first_matched_taxid_rank",  # NEW COLUMNS
            "protein_accession", "protein_length", 
            "nucleotide_accession", "nucleotide_length", "matched_rank", 
            "ncbi_taxonomy", "reference_name", "protein_reference_path", 
            "nucleotide_reference_path"
        ]

def test_initialisation_with_genbank(output_manager_with_genbank):
    """OutputManager initialises correctly with GenBank dirs when enabled."""
    # Check if dirs created
    assert output_manager_with_genbank.output_dir.exists()
    assert output_manager_with_genbank.protein_dir.exists()
    assert output_manager_with_genbank.nucleotide_dir.exists()
    
    # GenBank directories should be created when save_genbank=True
    assert output_manager_with_genbank.genbank_dir.exists()
    assert output_manager_with_genbank.protein_genbank_dir.exists()
    assert output_manager_with_genbank.nucleotide_genbank_dir.exists()


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
        "input_taxa": "Homo sapiens",
        "first_matched_taxid": "9606",
        "first_matched_taxid_rank": "species:Homo sapiens",
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
        assert row[1] == "Homo sapiens"
        assert row[2] == "9606"
        assert row[3] == "species:Homo sapiens"
        assert row[4] == "P12345"
        assert row[5] == "500"
        assert row[6] == "NM_12345"
        assert row[7] == "1500"
        assert row[8] == "species"
        assert row[9] == "Homo sapiens"
        assert row[10] == "test_process"
        assert row[11] == "protein/P12345.fasta"
        assert row[12] == "nucleotide/NM_12345.fasta"


def test_save_sequence_summary(output_manager):
    """Save summary of sequences to a CSV file."""
    # Create test sequence records with annotations
    seq1 = SeqRecord(Seq("ACGT" * 25), id="SEQ1", description="Test sequence 1")
    seq1.annotations["searched_taxid"] = "9606"
    
    seq2 = SeqRecord(Seq("TGCA" * 50), id="SEQ2", description="Test sequence 2")
    seq2.annotations["searched_taxid"] = "9606"
    
    sequences = [seq1, seq2]
    
    # Save summary
    output_manager.save_sequence_summary(sequences, "test")
    
    # Verify summary file was created
    summary_path = output_manager.output_dir / "fetched_test_sequences.csv"
    assert summary_path.exists()
    
    with open(summary_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header == ["Accession", "Length", "Description", "searched_taxid"]
        
        row1 = next(reader)
        assert row1[0] == "SEQ1"
        assert row1[1] == "100"
        assert row1[2] == "Test sequence 1"
        assert row1[3] == "9606"
        
        row2 = next(reader)
        assert row2[0] == "SEQ2"
        assert row2[1] == "200"
        assert row2[2] == "Test sequence 2"
        assert row2[3] == "9606"


@patch('gene_fetch.output_manager.EntrezHandler')
def test_save_genbank_file(mock_entrez, output_manager_with_genbank):
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
    output_path = output_manager_with_genbank.nucleotide_genbank_dir / "test_genbank.gb"
    
    # Call method
    result = output_manager_with_genbank.save_genbank_file(mock_entrez_instance, record_id, db, output_path)
    
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
def test_save_genbank_file_error(mock_entrez, output_manager_with_genbank):
    """Handle of errors when saving a GenBank file."""
    # Mock EntrezHandler's fetch method to raise an exception
    mock_entrez_instance = MagicMock()
    mock_entrez_instance.fetch.side_effect = Exception("Test error")
    
    # Define test parameters
    record_id = "NM_12345"
    db = "nucleotide"
    output_path = output_manager_with_genbank.nucleotide_genbank_dir / "test_error.gb"
    
    # Call method
    result = output_manager_with_genbank.save_genbank_file(mock_entrez_instance, record_id, db, output_path)
    
    # Verify result
    assert result is False
    assert not output_path.exists()


def test_save_genbank_file_disabled(output_manager):
    """Test GenBank file saving when disabled."""
    # Mock EntrezHandler
    mock_entrez_instance = MagicMock()
    
    # Define test parameters
    record_id = "NM_12345"
    db = "nucleotide"
    output_path = Path("/tmp/test_genbank.gb")  # Dummy path
    
    # Call method with save_genbank=False
    result = output_manager.save_genbank_file(mock_entrez_instance, record_id, db, output_path)
    
    # Verify result is False and no fetch was called
    assert result is False
    mock_entrez_instance.fetch.assert_not_called()
    
@pytest.fixture
def output_manager_no_refs():
    """Create OutputManager instance for testing with temp dir (no sequence refs)."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create OutputManager with the temporary directory but no sequence references
        output_dir = Path(tmpdirname)
        manager = OutputManager(output_dir, save_genbank=False, create_sequence_refs=False)
        yield manager


def test_initialisation_no_sequence_refs(output_manager_no_refs):
    """OutputManager initialises correctly without sequence_references.csv."""
    # Check if dirs created
    assert output_manager_no_refs.output_dir.exists()
    assert output_manager_no_refs.protein_dir.exists()
    assert output_manager_no_refs.nucleotide_dir.exists()
    
    # sequence_refs_path should be None
    assert output_manager_no_refs.sequence_refs_path is None
    
    # Check if failed_searches.csv was still created
    assert output_manager_no_refs.failed_searches_path.exists()
    
    # Verify headers in failed_searches.csv
    with open(output_manager_no_refs.failed_searches_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header == ["process_id", "taxid", "error_type", "timestamp"]


def test_write_sequence_reference_disabled(output_manager_no_refs):
    """Test that write_sequence_reference does nothing when disabled."""
    # Create test data
    test_data = {
        "process_id": "test_process",
        "taxid": "9606",
        "protein_accession": "P12345",
        "protein_length": "500"
    }
    
    # This should not raise an error and should not create any file
    output_manager_no_refs.write_sequence_reference(test_data)
    
    # Verify no sequence_references.csv was created
    sequence_refs_path = output_manager_no_refs.output_dir / "sequence_references.csv"
    assert not sequence_refs_path.exists()