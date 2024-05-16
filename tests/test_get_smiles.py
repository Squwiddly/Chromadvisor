#To run the test, put in the terminal : pytest -s <path/to/Chromadvisor/test_name_files.py
import pytest
from Chromadvisor_without_interface import get_smiles

# Test the function

def test_get_smiles_found():
    molecule_name = "water"
    result = get_smiles(molecule_name)
    assert result == "O"

def test_get_smiles_not_found():
    molecule_name = "nonexistent_molecule"
    result = get_smiles(molecule_name)
    assert result == "Molecule not found. Please try another name."

def test_get_smiles_nothing():
    molecule_name = ""
    result = get_smiles(molecule_name)
    assert "An error occurred" in result
