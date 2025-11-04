# tests/test_processors.py
import pytest
import tempfile
import csv
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from gene_fetch.processors import (
    process_sample,
    process_single_taxid,
    process_taxid_csv,
    process_taxonomy_csv
)
from gene_fetch.sequence_processor import SequenceProcessor
from gene_fetch.output_manager import OutputManager
from gene_fetch.entrez_handler import EntrezHandler

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


@pytest.fixture
def mock_processor():
    """Create mock SequenceProcessor for testing."""
    processor = MagicMock(spec=SequenceProcessor)
    
    # Mock search_and_fetch_sequences method
    protein_record = SeqRecord(Seq("ACGT" * 25), id="PROT123", description="Test protein")
    nucleotide_record = SeqRecord(Seq("ACGT" * 50), id="NM_12345", description="Test nucleotide")
    
    processor.search_and_fetch_sequences.return_value = (
        [protein_record],  # protein_records
        [nucleotide_record],  # nucleotide_records
        ["Homo", "sapiens"],  # taxonomy
        "species:Homo sapiens",  # matched_rank
        "9606",  # first_matched_taxid
        "species:Homo sapiens"  # first_matched_taxid_rank
    )
    
    # Mock entrez property
    processor.entrez = MagicMock(spec=EntrezHandler)
    processor.entrez.fetch_taxonomy.return_value = (
        ["Homo", "sapiens"],
        {"Homo sapiens": "species"},
        "species",
        {"Homo sapiens": "9606"}
    )
    
    return processor

@pytest.fixture
def mock_output_manager():
    """Create mock OutputManager for testing."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dir = Path(tmpdirname)
        # Create OutputManager with GenBank enabled and sequence refs enabled for testing
        manager = OutputManager(output_dir, save_genbank=True, create_sequence_refs=True)
        
        # Spy on write methods using MagicMock
        manager.write_sequence_reference = MagicMock(wraps=manager.write_sequence_reference)
        manager.log_failure = MagicMock(wraps=manager.log_failure)
        manager.save_sequence_summary = MagicMock(wraps=manager.save_sequence_summary)
        manager.save_genbank_file = MagicMock(wraps=manager.save_genbank_file)
        
        yield manager


@patch('gene_fetch.processors.SeqIO')
def test_process_sample_both_types(mock_seqio, mock_processor, mock_output_manager):
    """Process sample with both protein and nucleotide types."""
    # Set up test parameters
    process_id = "TEST123"
    taxid = "9606"
    sequence_type = "both"
    gene_name = "rbcl"
    
    # Call function
    process_sample(
        process_id=process_id,
        taxid=taxid,
        sequence_type=sequence_type,
        processor=mock_processor,
        output_manager=mock_output_manager,
        gene_name=gene_name,
        save_genbank=True
    )
    
    # Check search was performed correctly
    mock_processor.search_and_fetch_sequences.assert_called_once_with(
        taxid, gene_name, sequence_type
    )
    
    # Check sequences were written
    assert mock_seqio.write.call_count == 2  # Once for protein, once for nucleotide
    
    # Check GenBank files were saved
    assert mock_output_manager.save_genbank_file.call_count == 2  # Once for protein, once for nucleotide
    
    # Check sequence reference was written
    mock_output_manager.write_sequence_reference.assert_called_once()
    
    # Check no failure was logged
    mock_output_manager.log_failure.assert_not_called()


@patch('gene_fetch.processors.SeqIO')
def test_process_sample_protein_only(mock_seqio, mock_processor, mock_output_manager):
    """Process sample with protein type only."""
    # Set up test parameters
    process_id = "TEST123"
    taxid = "9606"
    sequence_type = "protein"
    gene_name = "rbcl"
    
    # Call function
    process_sample(
        process_id=process_id,
        taxid=taxid,
        sequence_type=sequence_type,
        processor=mock_processor,
        output_manager=mock_output_manager,
        gene_name=gene_name,
        save_genbank=False
    )
    
    # Check search was performed correctly
    mock_processor.search_and_fetch_sequences.assert_called_once_with(
        taxid, gene_name, sequence_type
    )
    
    # Check only protein sequence was written
    assert mock_seqio.write.call_count == 1
    
    # Check no GenBank files were saved (save_genbank=False)
    mock_output_manager.save_genbank_file.assert_not_called()
    
    # Check sequence reference was written
    mock_output_manager.write_sequence_reference.assert_called_once()


@patch('gene_fetch.processors.SeqIO')
def test_process_sample_no_sequences(mock_seqio, mock_processor, mock_output_manager):
    """Process sample when no sequences are found."""
    # Set up test parameters
    process_id = "TEST123"
    taxid = "9606"
    sequence_type = "both"
    gene_name = "rbcl"
    
    # Mock search_and_fetch_sequences method to return no sequences
    mock_processor.search_and_fetch_sequences.return_value = (
        [],  # protein_records
        [],  # nucleotide_records
        [],  # taxonomy
        "No match",  # matched_rank
        None, # first_matched_taxid
        None # first_matched_taxid_rank
    )
    
    # Call function
    process_sample(
        process_id=process_id,
        taxid=taxid,
        sequence_type=sequence_type,
        processor=mock_processor,
        output_manager=mock_output_manager,
        gene_name=gene_name
    )
    
    # Check search was performed correctly
    mock_processor.search_and_fetch_sequences.assert_called_once_with(
        taxid, gene_name, sequence_type
    )
    
    # Check no sequences were written
    mock_seqio.write.assert_not_called()
    
    # Check failure was logged
    mock_output_manager.log_failure.assert_called_once_with(
        process_id, taxid, "No sequences found"
    )
    
    # Check no sequence reference was written
    mock_output_manager.write_sequence_reference.assert_not_called()


@patch('gene_fetch.processors.SeqIO')
@patch('gene_fetch.processors.OutputManager')
def test_process_single_taxid(mock_output_manager_class, mock_seqio, mock_processor):
    """Process single taxid to fetch all available sequences."""
    # Set up test parameters
    taxid = "9606"
    gene_name = "rbcl"
    sequence_type = "both"
    
    # Create a mock output manager instance
    mock_output_manager_instance = MagicMock()
    mock_output_manager_instance.protein_dir = Path("/tmp/protein")
    mock_output_manager_instance.nucleotide_dir = Path("/tmp/nucleotide")
    mock_output_manager_instance.protein_genbank_dir = Path("/tmp/genbank/protein")
    mock_output_manager_instance.nucleotide_genbank_dir = Path("/tmp/genbank/nucleotide")
    mock_output_manager_class.return_value = mock_output_manager_instance
    
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dir = Path(tmpdirname)
        
        # Mock multiple records return for single-taxid mode
        protein_records = [
            SeqRecord(Seq("ACGT" * 25), id=f"PROT{i}", description=f"Test protein {i}")
            for i in range(1, 4)  # 3 protein records
        ]
        
        nucleotide_records = [
            SeqRecord(Seq("ACGT" * 50), id=f"NM_{i}", description=f"Test nucleotide {i}")
            for i in range(1, 5)  # 4 nucleotide records
        ]
        
        mock_processor.search_and_fetch_sequences.return_value = (
            protein_records,
            nucleotide_records,
            ["Homo", "sapiens"],
            "species:Homo sapiens",
            "9606",
            "species:Homo sapiens"
        )
        
        # Call function
        process_single_taxid(
            taxid=taxid,
            gene_name=gene_name,
            sequence_type=sequence_type,
            processor=mock_processor,
            output_dir=output_dir,
            max_sequences=5,
            save_genbank=True
        )
        
        # Check OutputManager was created with correct parameters
        mock_output_manager_class.assert_called_once_with(output_dir, True, create_sequence_refs=False)
        
        # Check search was performed correctly with fetch_all=True
        mock_processor.search_and_fetch_sequences.assert_called_once()
        args, kwargs = mock_processor.search_and_fetch_sequences.call_args
        assert args[0] == taxid
        assert args[1] == gene_name
        assert args[2] == sequence_type
        assert kwargs['fetch_all'] is True
        assert 'progress_counters' in kwargs
        
        # Check sequences were written (3 proteins + 4 nucleotides)
        assert mock_seqio.write.call_count == 7
        
        # Check GenBank files were saved
        assert mock_output_manager_instance.save_genbank_file.call_count == 7
        
        # Check summaries were saved
        assert mock_output_manager_instance.save_sequence_summary.call_count == 2


@patch('gene_fetch.processors.SeqIO')
@patch('gene_fetch.processors.OutputManager')
def test_process_single_taxid_max_limit(mock_output_manager_class, mock_seqio, mock_processor):
    """Process single taxid with max limit that actually limits the results."""
    # Set up test parameters
    taxid = "9606"
    gene_name = "rbcl"
    sequence_type = "both"
    max_sequences = 2  # Limit to 2 sequences of each type
    
    # Create a mock output manager instance
    mock_output_manager_instance = MagicMock()
    mock_output_manager_instance.protein_dir = Path("/tmp/protein")
    mock_output_manager_instance.nucleotide_dir = Path("/tmp/nucleotide")
    mock_output_manager_class.return_value = mock_output_manager_instance
    
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dir = Path(tmpdirname)
        
        # Mock multiple records return for single-taxid mode
        protein_records = [
            SeqRecord(Seq("ACGT" * 25), id=f"PROT{i}", description=f"Test protein {i}")
            for i in range(1, 6)  # 5 records
        ]
        
        nucleotide_records = [
            SeqRecord(Seq("ACGT" * 50), id=f"NM_{i}", description=f"Test nucleotide {i}")
            for i in range(1, 7)  # 6 records
        ]
        
        mock_processor.search_and_fetch_sequences.return_value = (
            protein_records,
            nucleotide_records,
            ["Homo", "sapiens"],
            "species:Homo sapiens",
            "9606",
            "species:Homo sapiens"
        )
        
        # Call function
        process_single_taxid(
            taxid=taxid,
            gene_name=gene_name,
            sequence_type=sequence_type,
            processor=mock_processor,
            output_dir=output_dir,
            max_sequences=max_sequences,
            save_genbank=False
        )
        
        # Check only max_sequences (2) of each type were processed
        assert mock_seqio.write.call_count == 4  # 2 proteins + 2 nucleotides
        
        # Check summaries were saved with the limited sequences
        call_args_list = mock_output_manager_instance.save_sequence_summary.call_args_list
        assert len(call_args_list) == 2
        
        # Verify first call (protein)
        protein_args = call_args_list[0][0]
        assert len(protein_args[0]) == 2  # Only 2 sequences
        assert protein_args[1] == "protein"
        
        # Verify second call (nucleotide)
        nucleotide_args = call_args_list[1][0]
        assert len(nucleotide_args[0]) == 2  # Only 2 sequences
        assert nucleotide_args[1] == "nucleotide"


@patch('gene_fetch.processors.csv.DictReader')
@patch('gene_fetch.processors.get_process_id_column')
@patch('gene_fetch.processors.log_progress')
@patch('gene_fetch.processors.sleep')
@patch('gene_fetch.processors.process_sample')
def test_process_taxid_csv(
    mock_process_sample, mock_sleep, mock_log_progress, 
    mock_get_process_id, mock_dict_reader, 
    mock_processor, mock_output_manager
):
    """Process CSV file with taxids."""
    # Set up test parameters
    csv_path = "test_samples.csv"
    gene_name = "rbcl"
    sequence_type = "both"
    
    # Mock CSV reader
    mock_get_process_id.return_value = "sample_id"
    
    # Create mock rows
    mock_rows = [
        {"sample_id": "SAMPLE1", "taxid": "9606"},
        {"sample_id": "SAMPLE2", "taxid": "9913"},
        {"sample_id": "SAMPLE3", "taxid": "10090"}
    ]
    
    # Set up DictReader mock to return rows
    mock_dict_reader.return_value.__iter__.return_value = mock_rows
    mock_dict_reader.return_value.fieldnames = ["sample_id", "taxid"]
    
    # Mock file opening and reading
    with patch('builtins.open', MagicMock()):
        # Call function
        process_taxid_csv(
            csv_path=csv_path,
            gene_name=gene_name,
            sequence_type=sequence_type,
            processor=mock_processor,
            output_manager=mock_output_manager,
            save_genbank=True
        )
        
        # Check process_sample was called for each row
        assert mock_process_sample.call_count == 3
        
        # Check parameters for first call
        args, kwargs = mock_process_sample.call_args_list[0]
        assert kwargs["process_id"] == "SAMPLE1"
        assert kwargs["taxid"] == "9606"
        assert kwargs["gene_name"] == gene_name
        assert kwargs["sequence_type"] == sequence_type
        assert kwargs["processor"] == mock_processor
        assert kwargs["output_manager"] == mock_output_manager
        assert kwargs["save_genbank"] is True
        
        # Check progress logging happened, at least at start and end
        assert mock_log_progress.call_count >= 2 
        
        # Check sleep was called between samples
        assert mock_sleep.call_count == 3


@patch('gene_fetch.processors.csv.DictReader')
@patch('gene_fetch.processors.get_process_id_column')
@patch('gene_fetch.processors.log_progress')
@patch('gene_fetch.processors.process_sample')
def test_process_taxonomy_csv(
    mock_process_sample, mock_log_progress, 
    mock_get_process_id, mock_dict_reader, 
    mock_processor, mock_output_manager
):
    """Process CSV file with taxonomy information."""
    # Set up test parameters
    csv_path = "test_taxonomy.csv"
    gene_name = "rbcl"
    sequence_type = "both"
    
    # Mock CSV reader
    mock_get_process_id.return_value = "ID"
    
    # Create mock rows with taxonomy information
    mock_rows = [
        {"ID": "SAMPLE1", "genus": "Homo", "species": "sapiens", "family": "Hominidae"},
        {"ID": "SAMPLE2", "genus": "Bos", "species": "taurus", "family": "Bovidae"},
        {"ID": "SAMPLE3", "genus": "Mus", "species": "musculus", "family": "Muridae"}
    ]
    
    # Set up DictReader mock to return rows
    mock_dict_reader.return_value.__iter__.return_value = mock_rows
    mock_dict_reader.return_value.fieldnames = ["ID", "genus", "species", "family"]
    
    # Mock EntrezHandler
    mock_entrez = MagicMock(spec=EntrezHandler)
    # Return different taxids for different genus-species combinations
    mock_entrez.fetch_taxid_from_taxonomy.side_effect = ["9606", "9913", "10090"]
    
    # Mock file opening and reading
    with patch('builtins.open', MagicMock()):
        # Call function
        process_taxonomy_csv(
            csv_path=csv_path,
            gene_name=gene_name,
            sequence_type=sequence_type,
            processor=mock_processor,
            output_manager=mock_output_manager,
            entrez=mock_entrez,
            save_genbank=True
        )
        
        # Check taxids were fetched for each row
        assert mock_entrez.fetch_taxid_from_taxonomy.call_count == 3
        
        # Check parameters for the first taxid fetch
        args, kwargs = mock_entrez.fetch_taxid_from_taxonomy.call_args_list[0]
        assert args[4] == "Homo"  # genus
        assert args[5] == "sapiens"  # species
        
        # Check process_sample was called for each row
        assert mock_process_sample.call_count == 3
        
        # Check parameters for first process_sample call
        args, kwargs = mock_process_sample.call_args_list[0]
        assert kwargs["process_id"] == "SAMPLE1"
        assert kwargs["taxid"] == "9606"  # Resolved taxid for Homo sapiens
        assert kwargs["gene_name"] == gene_name
        assert kwargs["sequence_type"] == sequence_type
        assert kwargs["processor"] == mock_processor
        assert kwargs["output_manager"] == mock_output_manager
        assert kwargs["save_genbank"] is True


@patch('gene_fetch.processors.csv.DictReader')
@patch('gene_fetch.processors.get_process_id_column')
def test_process_taxonomy_csv_missing_columns(
    mock_get_process_id, mock_dict_reader,
    mock_processor, mock_output_manager
):
    """Process taxonomy CSV with missing required columns."""
    # Set up test parameters
    csv_path = "test_taxonomy.csv"
    gene_name = "rbcl"
    sequence_type = "both"
    
    # Mock CSV reader - missing all taxonomic columns
    mock_get_process_id.return_value = "ID"
    mock_dict_reader.return_value.fieldnames = ["ID", "description", "notes"]
    
    # Mock EntrezHandler
    mock_entrez = MagicMock(spec=EntrezHandler)
    
    # Mock file opening and reading
    with patch('builtins.open', MagicMock()):
        # Call function - should exit with error
        with pytest.raises(SystemExit):
            process_taxonomy_csv(
                csv_path=csv_path,
                gene_name=gene_name,
                sequence_type=sequence_type,
                processor=mock_processor,
                output_manager=mock_output_manager,
                entrez=mock_entrez,
                save_genbank=True
            )


@patch('gene_fetch.processors.csv.DictReader')
@patch('gene_fetch.processors.get_process_id_column')
@patch('gene_fetch.processors.log_progress')
def test_process_taxonomy_csv_taxid_not_found(
    mock_log_progress, mock_get_process_id, mock_dict_reader,
    mock_processor, mock_output_manager
):
    """Process taxonomy CSV when taxid can't be resolved."""
    # Set up test parameters
    csv_path = "test_taxonomy.csv"
    gene_name = "rbcl"
    sequence_type = "both"
    
    # Mock CSV reader
    mock_get_process_id.return_value = "ID"
    
    # Create mock rows with taxonomy information
    mock_rows = [
        {"ID": "SAMPLE1", "genus": "Unknown", "species": "species", "family": ""}
    ]
    
    # Set up DictReader mock to return rows
    mock_dict_reader.return_value.__iter__.return_value = mock_rows
    mock_dict_reader.return_value.fieldnames = ["ID", "genus", "species", "family"]
    
    # Mock EntrezHandler to return None (no taxid found)
    mock_entrez = MagicMock(spec=EntrezHandler)
    mock_entrez.fetch_taxid_from_taxonomy.return_value = None
    
    # Mock file opening and reading
    with patch('builtins.open', MagicMock()):
        # Call function
        process_taxonomy_csv(
            csv_path=csv_path,
            gene_name=gene_name,
            sequence_type=sequence_type,
            processor=mock_processor,
            output_manager=mock_output_manager,
            entrez=mock_entrez,
            save_genbank=True
        )
        
        # Check failure was logged
        mock_output_manager.log_failure.assert_called_once_with(
            "SAMPLE1", "unknown", "Could not resolve taxid"
        )
        
        # Check process_sample was not called
        assert mock_processor.search_and_fetch_sequences.call_count == 0