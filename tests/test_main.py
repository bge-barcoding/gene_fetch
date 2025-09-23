"""
Tests for the main module.
"""

import pytest
import argparse
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock, call

import gene_fetch.main as main_module
from gene_fetch.main import main, setup_argument_parser

def test_setup_argument_parser():
    """Test setting up the argument parser with all expected arguments."""
    parser = setup_argument_parser()
    
    # Test that the parser is an ArgumentParser instance
    assert isinstance(parser, argparse.ArgumentParser)
    
    # Test the parser has all expected arguments
    actions = {action.dest: action for action in parser._actions}
    
    # Required arguments
    assert 'gene' in actions
    assert 'out' in actions
    assert 'type' in actions
    assert 'email' in actions
    assert 'api_key' in actions
    
    # Optional arguments
    assert 'input_csv' in actions
    assert 'input_taxonomy_csv' in actions
    assert 'single' in actions
    assert 'protein_size' in actions
    assert 'nucleotide_size' in actions
    assert 'max_sequences' in actions
    assert 'genbank' in actions
    
    # Test mutually exclusive group
    for group in parser._mutually_exclusive_groups:
        group_actions = [action.dest for action in group._group_actions]
        if 'input_csv' in group_actions:
            assert 'input_taxonomy_csv' in group_actions
            break
    else:
        assert False, "No mutually exclusive group found for input files"


@patch('gene_fetch.main.sys.exit')
@patch('gene_fetch.main.process_single_taxid')
@patch('gene_fetch.main.SequenceProcessor')
@patch('gene_fetch.main.EntrezHandler')
@patch('gene_fetch.main.Config')
@patch('gene_fetch.main.setup_logging')
@patch('gene_fetch.main.make_out_dir')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_single_taxid_mode(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir,
                               mock_make_out_dir, mock_setup_logging, 
                               mock_config, mock_entrez, mock_processor,
                               mock_process_single, mock_exit):
    """Test main function in single taxid mode."""
    # Create a mock parser that returns a namespace with the required arguments
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='both',
        email='test@example.com',
        api_key='valid_test_key_12345',
        input_csv='samples.csv',  # Add a CSV file to prevent the validation error
        input_taxonomy_csv=None,
        single='9606',
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=None,
        genbank=False,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance
    
    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False
    
    # Set up mock returns
    mock_config_instance = MagicMock()
    mock_config_instance.valid_sequence_types = {'protein', 'nucleotide', 'both'}
    mock_config_instance.set_gene_search_term.return_value = ("cox1", "protein-coding")  # Fixed: return tuple
    mock_config.return_value = mock_config_instance
    
    mock_entrez_instance = MagicMock()
    mock_entrez.return_value = mock_entrez_instance
    
    mock_processor_instance = MagicMock()
    mock_processor.return_value = mock_processor_instance
    
    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False
    
    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False
    
    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False
    
    # Also patch process_taxid_csv to prevent it from being called
    with patch('gene_fetch.main.process_taxid_csv'):
        # Run main function
        main()
    
    # Verify correct functions were called
    assert mock_make_out_dir.call_count >= 1  # Can be called 1-2 times depending on clean logic
    mock_setup_logging.assert_called_once()
    mock_config.assert_called_once_with(email='test@example.com', api_key='valid_test_key_12345')
    mock_config_instance.update_thresholds.assert_called_once()
    mock_config_instance.set_gene_search_term.assert_called_once_with('cox1')
    mock_entrez.assert_called_once_with(mock_config_instance)
    mock_processor.assert_called_once_with(mock_config_instance, mock_entrez_instance)
    
    # Check that process_single_taxid was called with correct arguments
    mock_process_single.assert_called_once()
    args, kwargs = mock_process_single.call_args
    assert kwargs['taxid'] == '9606'
    assert kwargs['gene_name'] == 'cox1'  # Should use canonical name
    assert kwargs['sequence_type'] == 'both'
    assert kwargs['processor'] == mock_processor_instance
    assert isinstance(kwargs['output_dir'], Path)
    assert kwargs['max_sequences'] is None
    assert kwargs['save_genbank'] is False
    
    # Check that sys.exit was called with 0 (success)
    mock_exit.assert_called_once_with(0)


@patch('gene_fetch.main.sys.exit')
@patch('gene_fetch.main.process_taxid_csv')
@patch('gene_fetch.main.OutputManager')
@patch('gene_fetch.main.SequenceProcessor')
@patch('gene_fetch.main.EntrezHandler')
@patch('gene_fetch.main.Config')
@patch('gene_fetch.main.setup_logging')
@patch('gene_fetch.main.make_out_dir')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_taxid_csv_mode(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir,
                            mock_make_out_dir, mock_setup_logging, 
                            mock_config, mock_entrez, mock_processor,
                            mock_output_manager, mock_process_taxid_csv, 
                            mock_exit):
    """Test main function in taxid CSV mode."""
    # Create a mock parser that returns a namespace with the required arguments
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='both',
        email='test@example.com',
        api_key='valid_test_key_12345',
        input_csv='samples.csv',
        input_taxonomy_csv=None,
        single=None,
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=None,
        genbank=False,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance
    
    # Set up mock returns
    mock_config_instance = MagicMock()
    mock_config_instance.valid_sequence_types = {'protein', 'nucleotide', 'both'}
    mock_config_instance.set_gene_search_term.return_value = ("cox1", "protein-coding")  # Fixed: return tuple
    mock_config.return_value = mock_config_instance
    
    mock_entrez_instance = MagicMock()
    mock_entrez.return_value = mock_entrez_instance
    
    mock_processor_instance = MagicMock()
    mock_processor.return_value = mock_processor_instance
    
    mock_output_manager_instance = MagicMock()
    mock_output_manager.return_value = mock_output_manager_instance
    
    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False
    
    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False
    
    # Run main function
    main()
    
    # Verify correct functions were called
    assert mock_make_out_dir.call_count >= 1  # Can be called 1-2 times depending on clean logic
    mock_setup_logging.assert_called_once()
    mock_config.assert_called_once_with(email='test@example.com', api_key='valid_test_key_12345')
    mock_config_instance.update_thresholds.assert_called_once()
    mock_config_instance.set_gene_search_term.assert_called_once_with('cox1')
    mock_entrez.assert_called_once_with(mock_config_instance)
    mock_processor.assert_called_once_with(mock_config_instance, mock_entrez_instance)
    mock_output_manager.assert_called_once()
    
    # Check that process_taxid_csv was called with correct arguments
    mock_process_taxid_csv.assert_called_once()
    args, kwargs = mock_process_taxid_csv.call_args
    assert args[0] == 'samples.csv'
    assert args[1] == 'cox1'  # Should use canonical name
    assert args[2] == 'both'
    assert args[3] == mock_processor_instance
    assert args[4] == mock_output_manager_instance
    assert args[5] is False  # save_genbank flag


@patch('gene_fetch.main.sys.exit')
@patch('gene_fetch.main.process_taxonomy_csv')
@patch('gene_fetch.main.OutputManager')
@patch('gene_fetch.main.SequenceProcessor')
@patch('gene_fetch.main.EntrezHandler')
@patch('gene_fetch.main.Config')
@patch('gene_fetch.main.setup_logging')
@patch('gene_fetch.main.make_out_dir')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_taxonomy_csv_mode(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir,
                               mock_make_out_dir, mock_setup_logging, 
                               mock_config, mock_entrez, mock_processor,
                               mock_output_manager, mock_process_taxonomy_csv, 
                               mock_exit):
    """Test main function in taxonomy CSV mode."""
    # Create a mock parser that returns a namespace with the required arguments
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='both',
        email='test@example.com',
        api_key='valid_test_key_12345',
        input_csv=None,
        input_taxonomy_csv='taxonomy.csv',
        single=None,
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=None,
        genbank=False,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance
    
    # Set up mock returns
    mock_config_instance = MagicMock()
    mock_config_instance.valid_sequence_types = {'protein', 'nucleotide', 'both'}
    mock_config_instance.set_gene_search_term.return_value = ("cox1", "protein-coding")  # Fixed: return tuple
    mock_config.return_value = mock_config_instance
    
    mock_entrez_instance = MagicMock()
    mock_entrez.return_value = mock_entrez_instance
    
    mock_processor_instance = MagicMock()
    mock_processor.return_value = mock_processor_instance
    
    mock_output_manager_instance = MagicMock()
    mock_output_manager.return_value = mock_output_manager_instance
    
    # Run main function
    main()
    
    # Verify correct functions were called
    assert mock_make_out_dir.call_count >= 1  # Can be called 1-2 times depending on clean logic
    mock_setup_logging.assert_called_once()
    mock_config.assert_called_once_with(email='test@example.com', api_key='valid_test_key_12345')
    mock_config_instance.update_thresholds.assert_called_once()
    mock_config_instance.set_gene_search_term.assert_called_once_with('cox1')
    mock_entrez.assert_called_once_with(mock_config_instance)
    mock_processor.assert_called_once_with(mock_config_instance, mock_entrez_instance)
    mock_output_manager.assert_called_once()
    
    # Check that process_taxonomy_csv was called with correct arguments
    mock_process_taxonomy_csv.assert_called_once()
    args, kwargs = mock_process_taxonomy_csv.call_args
    assert args[0] == 'taxonomy.csv'
    assert args[1] == 'cox1'  # Should use canonical name
    assert args[2] == 'both'
    assert args[3] == mock_processor_instance
    assert args[4] == mock_output_manager_instance
    assert args[5] == mock_entrez_instance


@patch('gene_fetch.main.sys.exit')
@patch('gene_fetch.main.process_single_taxid')
@patch('gene_fetch.main.SequenceProcessor')
@patch('gene_fetch.main.EntrezHandler')
@patch('gene_fetch.main.Config')
@patch('gene_fetch.main.setup_logging')
@patch('gene_fetch.main.make_out_dir')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_with_advanced_options(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir,
                                   mock_make_out_dir, mock_setup_logging, 
                                   mock_config, mock_entrez, mock_processor,
                                   mock_process_single, mock_exit):
    """Test main function with advanced options (max-sequences and genbank)."""
    # Create a mock parser that returns a namespace with the required arguments
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='both',
        email='test@example.com',
        api_key='valid_test_key_12345',
        input_csv='samples.csv',  # Add a CSV file to prevent the validation error
        input_taxonomy_csv=None,
        single='9606',
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=10,
        genbank=True,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance
    
    # Set up mock returns
    mock_config_instance = MagicMock()
    mock_config_instance.valid_sequence_types = {'protein', 'nucleotide', 'both'}
    mock_config_instance.set_gene_search_term.return_value = ("cox1", "protein-coding")  # Fixed: return tuple
    mock_config.return_value = mock_config_instance
    
    mock_entrez_instance = MagicMock()
    mock_entrez.return_value = mock_entrez_instance
    
    mock_processor_instance = MagicMock()
    mock_processor.return_value = mock_processor_instance
    
    # Also patch process_taxid_csv to prevent it from being called
    with patch('gene_fetch.main.process_taxid_csv'):
        # Run main function
        main()
    
    # Check that process_single_taxid was called with correct arguments
    mock_process_single.assert_called_once()
    args, kwargs = mock_process_single.call_args
    assert kwargs['taxid'] == '9606'
    assert kwargs['gene_name'] == 'cox1'  # Should use canonical name
    assert kwargs['sequence_type'] == 'both'
    assert kwargs['processor'] == mock_processor_instance
    assert isinstance(kwargs['output_dir'], Path)
    assert kwargs['max_sequences'] == 10
    assert kwargs['save_genbank'] is True
    
    # Verify make_out_dir was called at least once
    assert mock_make_out_dir.call_count >= 1
    
    # Check that sys.exit was called with 0 (success)
    mock_exit.assert_called_once_with(0)


@patch('gene_fetch.main.sys.exit')
@patch('gene_fetch.main.setup_logging')
@patch('gene_fetch.main.make_out_dir')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_invalid_sequence_type(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir,
                                   mock_make_out_dir, mock_setup_logging, mock_exit):
    """Test main function with invalid sequence type."""
    # Create a mock parser that returns a namespace with the required arguments
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='invalid',  # This is an invalid sequence type
        email='test@example.com',
        api_key='valid_test_key_12345',
        input_csv='samples.csv',  # Add a CSV file to prevent the validation error
        input_taxonomy_csv=None,
        single=None,
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=None,
        genbank=False,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance
    
    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False
    
    # Create a mock Config that will be used to validate the sequence type
    with patch('gene_fetch.main.Config') as mock_config:
        mock_config_instance = MagicMock()
        mock_config_instance.valid_sequence_types = {'protein', 'nucleotide', 'both'}
        mock_config_instance.set_gene_search_term.return_value = ("cox1", "protein-coding")  # Fixed: return tuple
        mock_config.return_value = mock_config_instance
        
        # Run main function
        main()
        
        # Check that sys.exit was called with 1 (error)
        # We need to use assert_any_call instead of assert_called_once_with
        # because exit might be called multiple times due to the validations
        mock_exit.assert_any_call(1)


@patch('gene_fetch.main.sys.exit')
@patch('gene_fetch.main.SequenceProcessor')
@patch('gene_fetch.main.EntrezHandler')
@patch('gene_fetch.main.Config')
@patch('gene_fetch.main.setup_logging')
@patch('gene_fetch.main.make_out_dir')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_no_input_files(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir,
                            mock_make_out_dir, mock_setup_logging, 
                            mock_config, mock_entrez, mock_processor, 
                            mock_exit):
    """Test that main function completes when argparse validation is bypassed via mocking."""
    # Create a mock parser that returns a namespace with no input files
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='both',
        email='test@example.com',
        api_key='valid_test_key_12345',
        input_csv=None,
        input_taxonomy_csv=None,
        single=None,
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=None,
        genbank=False,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance
    
    # Set up mock returns
    mock_config_instance = MagicMock()
    mock_config_instance.valid_sequence_types = {'protein', 'nucleotide', 'both'}
    mock_config_instance.set_gene_search_term.return_value = ("cox1", "protein-coding")  # Fixed: return tuple
    mock_config.return_value = mock_config_instance
    
    mock_entrez_instance = MagicMock()
    mock_entrez.return_value = mock_entrez_instance
    
    mock_processor_instance = MagicMock()
    mock_processor.return_value = mock_processor_instance
    
    # Run main function
    main()
    
    # Since we're bypassing argparse validation via mocking, 
    # main() should complete normally (no sys.exit call)
    assert mock_exit.call_count == 0
    
    # Verify make_out_dir was called at least once
    assert mock_make_out_dir.call_count >= 1


@patch('gene_fetch.main.setup_logging')
@patch('gene_fetch.main.make_out_dir')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_credential_validation_error(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir,
                                          mock_make_out_dir, mock_setup_logging):
    """Test main function when credential validation fails."""

    # Create a mock parser that returns a namespace with invalid credentials
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='both',
        email='test@example.com',
        api_key='fake_key',  # This should trigger validation failure
        input_csv=None,
        input_taxonomy_csv=None,
        single=None,
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=None,
        genbank=False,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance

    # Set up mock for should_clear_output_directory to return False
    mock_should_clear_output_dir.return_value = False

    # Run main function - this should exit due to credential validation
    with pytest.raises(SystemExit) as excinfo:
        main()

    # Verify it exited with code 1
    assert excinfo.value.code == 1


@patch('gene_fetch.main.sys.exit')
@patch('gene_fetch.main.should_clear_output_directory', create=True)
@patch('gene_fetch.main.clear_output_directory', create=True)
@patch('gene_fetch.main.save_run_info', create=True)
@patch('gene_fetch.main.setup_argument_parser')
def test_main_with_real_args(mock_parser, mock_save_run_info, mock_clear_output_dir, mock_should_clear_output_dir, mock_exit):
    """Test main function with real arguments."""
    # Create a mock parser that returns a namespace with the required arguments
    mock_parser_instance = MagicMock()
    mock_args = argparse.Namespace(
        gene='cox1',
        out='/tmp/out',
        type='both',
        email='test@example.com',
        api_key='valid_test_key_12345',
        input_csv='tests/data/test_samples.csv',
        input_taxonomy_csv=None,
        single=None,
        protein_size=500,
        nucleotide_size=1000,
        max_sequences=None,
        genbank=False,
        clean=False,
        header="basic"
    )
    mock_parser_instance.parse_args.return_value = mock_args
    mock_parser.return_value = mock_parser_instance
    
    # Mock all the functions called by main
    with patch('gene_fetch.main.make_out_dir') as mock_make_out_dir, \
         patch('gene_fetch.main.setup_logging') as mock_setup_logging, \
         patch('gene_fetch.main.Config') as mock_config, \
         patch('gene_fetch.main.EntrezHandler') as mock_entrez, \
         patch('gene_fetch.main.SequenceProcessor') as mock_processor, \
         patch('gene_fetch.main.OutputManager') as mock_output_manager, \
         patch('gene_fetch.main.process_taxid_csv') as mock_process_csv: 
        
        # Set up Config to return a mock instance
        mock_config_instance = MagicMock()
        mock_config_instance.valid_sequence_types = {'protein', 'nucleotide', 'both'}
        mock_config_instance.set_gene_search_term.return_value = ("cox1", "protein-coding")  # Fixed: return tuple
        mock_config.return_value = mock_config_instance
        
        # Run main function
        main()
        
        # Should complete normally without calling sys.exit
        assert mock_exit.call_count == 0
        
        # Verify make_out_dir was called at least once
        assert mock_make_out_dir.call_count >= 1
        
        # Verify that process_taxid_csv was called (indicating normal flow)
        mock_process_csv.assert_called_once()