import sys
import os
from pathlib import Path
import pytest

# Explicitly add src directory to path if needed
# src_path = Path(__file__).parent.parent / "src"
# sys.path.insert(0, str(src_path))

def test_import():
    """Simple test to verify imports work."""
    from gene_fetch.entrez_handler import EntrezHandler
    from gene_fetch.core import Config
    
    # Create Config with required arguments
    config = Config(email="test@example.com", api_key="valid_test_key_12345")
    
    # Create the EntrezHandler
    handler = EntrezHandler(config)
    
    # Verify objects are created correctly
    assert handler is not None
    assert config.email == "test@example.com"
    assert config.api_key == "valid_test_key_12345"
    
    print("Successfully imported and created EntrezHandler!")
