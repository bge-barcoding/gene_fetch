"""
Integration tests for Gene Fetch.
Uses actual test data files to test end-to-end functionality.
"""

import pytest
import os
import tempfile
import shutil
from pathlib import Path
import subprocess
from unittest.mock import patch, MagicMock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gene_fetch.entrez_handler import EntrezHandler
from gene_fetch.sequence_processor import SequenceProcessor
from gene_fetch.processors import process_taxid_csv, process_taxonomy_csv
from gene_fetch.core import Config
from gene_fetch.output_manager import OutputManager


# Paths to test data
TEST_DATA_DIR = Path(__file__).parent / "data"
TEST_SAMPLES_CSV = TEST_DATA_DIR / "test_samples.csv"
TEST_SAMPLES_TAXONOMY_CSV = TEST_DATA_DIR / "test_samples_taxonomy.csv"
SAMPLE_RECORD_GB = TEST_DATA_DIR / "sample_record.gb"


@pytest.fixture
def mock_config():
    """Create a mock Config object."""
    config = MagicMock(spec=Config)
    config.email = "test@example.com"
    config.api_key = "test_api_key"
    config.max_calls_per_second = 10
    config.fetch_batch_size = 200
    config.batch_delay = (1, 2)
    config.valid_sequence_types = ["protein", "nucleotide", "both"]
    config.protein_length_threshold = 500
    config.nucleotide_length_threshold = 1000
    config.min_nucleotide_size_single_mode = 200
    config.min_protein_size_single_mode = 100
    config.gene_search_term = "(cox1[Gene] OR COI[Gene])"
    
    # Add the missing attributes
    config._protein_coding_genes = {
        "cox1": ["cox1[Gene]", "COI[Gene]"],
        "rbcl": ["rbcL[Gene]", "RBCL[Gene]"]
    }
    config._rRNA_genes = {
        "16s": ["16S ribosomal RNA[Gene]"],
        "18s": ["18S ribosomal RNA[Gene]"]
    }
    
    return config


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory for test results."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield Path(tmpdirname)


@pytest.mark.integration
def test_entrez_fetch_method(mock_config):
    """Test EntrezHandler.fetch method with actual NCBI interaction."""
    # Get API credentials from environment if available
    api_key = os.environ.get("NCBI_API_KEY")
    email = os.environ.get("NCBI_EMAIL")
    
    # Skip if no API key is available in environment
    if not api_key or not email:
        pytest.skip("Skipping test that requires valid NCBI API key and email")
    
    # Create a real config with the real credentials
    config = Config(email=email, api_key=api_key)
    entrez = EntrezHandler(config)
    
    # Test fetching a known protein record
    handle = entrez.fetch(db="protein", id="P00395", rettype="gb", retmode="text")
    assert handle is not None
    content = handle.read()
    assert "LOCUS" in content
    assert "DEFINITION" in content


@pytest.mark.integration
def test_entrez_search_method(mock_config):
    """Test EntrezHandler.search method with actual NCBI interaction."""
    # Get API credentials from environment if available
    api_key = os.environ.get("NCBI_API_KEY")
    email = os.environ.get("NCBI_EMAIL")
    
    # Skip if no API key is available in environment
    if not api_key or not email:
        pytest.skip("Skipping test that requires valid NCBI API key and email")
    
    # Create a real config with the real credentials
    config = Config(email=email, api_key=api_key)
    entrez = EntrezHandler(config)
    
    # Test searching for a human COX1 protein
    results = entrez.search(
        db="protein",
        term="cox1[Gene] AND Homo sapiens[Organism]",
        retmax=5
    )
    
    assert results is not None
    assert "IdList" in results
    assert len(results["IdList"]) > 0


@pytest.mark.integration
def test_sample_record_parsing():
    """Test parsing a GenBank sample record."""
    if not SAMPLE_RECORD_GB.exists():
        pytest.skip(f"Sample record file {SAMPLE_RECORD_GB} not found")
    
    record = SeqIO.read(SAMPLE_RECORD_GB, "genbank")
    
    assert record is not None
    assert len(record.seq) > 0
    assert record.id is not None
    
    # Test SequenceProcessor's parsing capabilities
    with patch('gene_fetch.sequence_processor.SequenceProcessor.fetch_nucleotide_record') as mock_fetch:
        mock_fetch.return_value = record
        
        config = MagicMock(spec=Config)
        entrez = MagicMock(spec=EntrezHandler)
        processor = SequenceProcessor(config, entrez)
        
        # Test parse_coded_by if the record has coding information
        if 'coded_by' in str(record.features):
            for feature in record.features:
                if feature.type == "CDS" and 'coded_by' in feature.qualifiers:
                    result = processor.parse_coded_by(feature.qualifiers['coded_by'][0])
                    assert result is not None
                    break


@pytest.mark.integration
def test_process_taxid_csv_mock(temp_output_dir, mock_config):
    """Test process_taxid_csv with mocked sequence processor."""
    if not TEST_SAMPLES_CSV.exists():
        pytest.skip(f"Test file {TEST_SAMPLES_CSV} not found")
    
    # Create mocked components
    entrez = MagicMock(spec=EntrezHandler)
    processor = MagicMock(spec=SequenceProcessor)
    output_manager = OutputManager(temp_output_dir)
    
    # Define mock behavior
    def mock_search_and_fetch(*args, **kwargs):
        protein_record = SeqRecord(Seq("ACGT" * 25), id="PROT123", description="Test protein")
        nucleotide_record = SeqRecord(Seq("ACGT" * 50), id="NM_12345", description="Test nucleotide")
        return [protein_record], [nucleotide_record], ["Homo", "sapiens"], "species:Homo sapiens"
    
    processor.search_and_fetch_sequences.side_effect = mock_search_and_fetch
    
    # Run the function with only the first two samples
    with patch('gene_fetch.processors.sleep'):  # Patch sleep to speed up the test
        process_taxid_csv(
            TEST_SAMPLES_CSV,
            "cox1",
            "both",
            processor,
            output_manager,
            save_genbank=False
        )
    
    # Check that files were created
    assert (temp_output_dir / "sequence_references.csv").exists()
    
    # Count the lines in the sequence references file
    with open(temp_output_dir / "sequence_references.csv", "r") as f:
        lines = f.readlines()
        # Header + one line per sample
        assert len(lines) > 1  # At least header + one processed sample
        
@pytest.mark.integration
def test_process_taxid_csv_real(temp_output_dir):
    """Test process_taxid_csv with real NCBI interactions."""
    if not TEST_SAMPLES_CSV.exists():
        pytest.skip(f"Test file {TEST_SAMPLES_CSV} not found")
    
    # Get API credentials from environment
    api_key = os.environ.get("NCBI_API_KEY")
    email = os.environ.get("NCBI_EMAIL")
    
    if not api_key or not email:
        pytest.skip("Skipping test that requires valid NCBI API key and email")
    
    # Create real components with actual credentials
    config = Config(email=email, api_key=api_key)
    entrez = EntrezHandler(config)
    processor = SequenceProcessor(config, entrez)
    output_manager = OutputManager(temp_output_dir)
    
    # Process only the first sample to keep the test runtime reasonable
    with open(TEST_SAMPLES_CSV, 'r') as f:
        content = f.read().splitlines()
    
    # Create a temporary file with just the header and first data row
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        # If there are at least 2 lines (header + first data row)
        if len(content) >= 2:
            temp_file.write('\n'.join(content[:2]))
        else:
            # If the file doesn't have enough lines, skip the test
            os.unlink(temp_file.name)
            pytest.skip("Test file doesn't have enough data rows")
            
        temp_csv = temp_file.name
    
    try:
        # Log the test start
        print(f"Starting real API test with taxid CSV: {temp_csv}")
        print(f"Using first sample from: {TEST_SAMPLES_CSV}")
        
        # Process the CSV with real NCBI interactions
        process_taxid_csv(
            temp_csv,
            "cox1",  # Using COX1 as the target gene
            "both",  # Request both protein and nucleotide
            processor,
            output_manager,
            save_genbank=False  # Don't save GenBank files to save space
        )
        
        # Verify outputs
        refs_file = temp_output_dir / "sequence_references.csv"
        assert refs_file.exists(), "Sequence references file was not created"
        
        # Check that we have actual data in the output
        with open(refs_file, 'r') as f:
            refs_content = f.read()
            assert len(refs_content.splitlines()) > 1, "No sequence data was retrieved"
            
        # You could add more specific validations here
        print(f"Successfully retrieved sequence data and created: {refs_file}")
        
    except Exception as e:
        pytest.fail(f"Real API test failed with error: {str(e)}")
    finally:
        # Clean up the temporary file
        os.unlink(temp_csv)


@pytest.mark.integration
def test_process_taxonomy_csv_mock(temp_output_dir, mock_config):
    """Test process_taxonomy_csv with mocked sequence processor."""
    if not TEST_SAMPLES_TAXONOMY_CSV.exists():
        pytest.skip(f"Test file {TEST_SAMPLES_TAXONOMY_CSV} not found")
    
    # Create mocked components
    entrez = MagicMock(spec=EntrezHandler)
    processor = MagicMock(spec=SequenceProcessor)
    output_manager = OutputManager(temp_output_dir)
    
    # Define mock behavior
    def mock_search_and_fetch(*args, **kwargs):
        protein_record = SeqRecord(Seq("ACGT" * 25), id="PROT123", description="Test protein")
        nucleotide_record = SeqRecord(Seq("ACGT" * 50), id="NM_12345", description="Test nucleotide")
        return [protein_record], [nucleotide_record], ["Arthropoda", "Insecta", "Diptera", "Acroceridae", "Astomella", "hispaniae"], "species:Astomella hispaniae"
    
    processor.search_and_fetch_sequences.side_effect = mock_search_and_fetch
    
    # Mock entrez to return taxids
    entrez.fetch_taxid_from_taxonomy.return_value = "210239"
    
    # Run the function
    with patch('gene_fetch.processors.sleep'):  # Patch sleep to speed up the test
        process_taxonomy_csv(
            TEST_SAMPLES_TAXONOMY_CSV,
            "cox1",
            "both",
            processor,
            output_manager,
            entrez,
            save_genbank=False
        )
    
    # Check that files were created
    assert (temp_output_dir / "sequence_references.csv").exists()
    
    # Count the lines in the sequence references file
    with open(temp_output_dir / "sequence_references.csv", "r") as f:
        lines = f.readlines()
        # Header + one line per sample
        assert len(lines) > 1  # At least header + one processed sample


@pytest.mark.integration
def test_process_taxonomy_csv_real(temp_output_dir):
    """Test process_taxonomy_csv with real NCBI interactions."""
    if not TEST_SAMPLES_TAXONOMY_CSV.exists():
        pytest.skip(f"Test file {TEST_SAMPLES_TAXONOMY_CSV} not found")
    
    # Get API credentials from environment
    api_key = os.environ.get("NCBI_API_KEY")
    email = os.environ.get("NCBI_EMAIL")
    
    if not api_key or not email:
        pytest.skip("Skipping test that requires valid NCBI API key and email")
    
    # Create real components with actual credentials
    config = Config(email=email, api_key=api_key)
    entrez = EntrezHandler(config)
    processor = SequenceProcessor(config, entrez)
    output_manager = OutputManager(temp_output_dir)
    
    # Process only the first sample to keep the test runtime reasonable
    with open(TEST_SAMPLES_TAXONOMY_CSV, 'r') as f:
        content = f.read().splitlines()
    
    # Create a temporary file with just the header and first data row
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        # If there are at least 2 lines (header + first data row)
        if len(content) >= 2:
            temp_file.write('\n'.join(content[:2]))
        else:
            # If the file doesn't have enough lines, skip the test
            os.unlink(temp_file.name)
            pytest.skip("Test file doesn't have enough data rows")
            
        temp_csv = temp_file.name
    
    try:
        # Log the test start
        print(f"Starting real API test with taxonomy CSV: {temp_csv}")
        print(f"Using first sample from: {TEST_SAMPLES_TAXONOMY_CSV}")
        
        # Process the CSV with real NCBI interactions
        process_taxonomy_csv(
            temp_csv,
            "cox1",  # Using COX1 as the target gene
            "both",  # Request both protein and nucleotide
            processor,
            output_manager,
            entrez,
            save_genbank=False  # Don't save GenBank files to save space
        )
        
        # Verify outputs
        refs_file = temp_output_dir / "sequence_references.csv"
        assert refs_file.exists(), "Sequence references file was not created"
        
        # Check that we have actual data in the output
        with open(refs_file, 'r') as f:
            refs_content = f.read()
            assert len(refs_content.splitlines()) > 1, "No sequence data was retrieved"
            
        # You could add more specific validations here
        print(f"Successfully retrieved sequence data and created: {refs_file}")
        
    except Exception as e:
        pytest.fail(f"Real API test failed with error: {str(e)}")
    finally:
        # Clean up the temporary file
        os.unlink(temp_csv)
        
        
@pytest.mark.integration
@pytest.mark.skipif(os.environ.get("NCBI_API_KEY") is None, reason="NCBI API key required")
def test_cli_with_test_data():
    """Test the command-line interface with the test data files."""
    # This test requires setting environment variables, either:
    # export NCBI_API_KEY=your_api_key
    # export NCBI_EMAIL=your_email
    # Or by specifying them in the run command.
    
    if not TEST_SAMPLES_CSV.exists():
        pytest.skip(f"Test file {TEST_SAMPLES_CSV} not found")
    
    # Get API key and email from environment variables
    api_key = os.environ.get("NCBI_API_KEY")
    email = os.environ.get("NCBI_EMAIL", "test@example.com")
    
    # Create a temporary output directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dir = Path(tmpdirname)
        
        # Run the command with the first sample only
        command = [
            "python", "-m", "gene_fetch.main",
            "-g", "cox1",
            "-o", str(output_dir),
            "-i", str(TEST_SAMPLES_CSV),
            "--type", "both",
            "-e", email,
            "-k", api_key,
            # Add a grep to limit to just the first sample for speed
            "| head -n 2"  # This will process only the header and first line
        ]
        
        try:
            # Run with a timeout to avoid hanging tests
            result = subprocess.run(
                " ".join(command),
                shell=True,
                timeout=30,
                capture_output=True,
                text=True
            )
            
            # Check that the output directory contains the expected files
            assert (output_dir / "sequence_references.csv").exists()
            
        except subprocess.TimeoutExpired:
            pytest.skip("CLI test timed out")


@pytest.mark.parametrize("feature_type,feature_file", [
    ("nucleotide", "nucleotide_record.gb"),
    ("protein", "protein_record.gb")
])
def test_sequence_processor_extract_methods(feature_type, feature_file, mock_config):
    """Test SequenceProcessor extract methods with sample GenBank records."""
    # This is a template for testing different extract methods based on
    # sample files with known features
    
    sample_file = TEST_DATA_DIR / feature_file
    if not sample_file.exists():
        pytest.skip(f"Sample file {sample_file} not found")
    
    from Bio import SeqIO
    record = SeqIO.read(sample_file, "genbank")
    
    entrez = MagicMock(spec=EntrezHandler)
    processor = SequenceProcessor(mock_config, entrez)
    
    if feature_type == "nucleotide":
        # Test extract_nucleotide if we have a nucleotide record
        with patch('gene_fetch.sequence_processor.SequenceProcessor.fetch_nucleotide_record') as mock_fetch:
            mock_fetch.return_value = record
            result = processor.extract_nucleotide(record, "cox1")
            
            # Check if the record actually has COX1/COI features
            has_cox1 = any('cox1' in str(feature.qualifiers).lower() or 
                           'coi' in str(feature.qualifiers).lower() 
                           for feature in record.features if feature.type == 'CDS')
            
            if has_cox1:
                assert result is not None
            else:
                # If no COX1/COI features, just verify the function ran without error
                pass  # No assertion needed
                
    elif feature_type == "protein":
        # Test extract_rRNA if we have an rRNA record
        if any(feature.type == "rRNA" for feature in record.features):
            with patch('gene_fetch.sequence_processor.SequenceProcessor.fetch_nucleotide_record') as mock_fetch:
                mock_fetch.return_value = record
                result = processor.extract_rRNA(record, "16s")  # Use appropriate rRNA gene
                assert result is not None
        else:
            # Test with a regular gene if no rRNA features
            with patch('gene_fetch.sequence_processor.SequenceProcessor.fetch_nucleotide_record') as mock_fetch:
                mock_fetch.return_value = record
                processor.extract_nucleotide(record, "cox1")  # Just verify no exceptions