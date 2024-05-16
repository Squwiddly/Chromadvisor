
from Chromadvisor_without_interface import on_submit

def test_on_submit_molecule_found():
    # Mocking the get_smiles function to return a valid SMILES string
    molecule_name = "Acetone"
    expected_result = {'ketone': {'count': 1, }} #'positions': [(1, 2)]}}
    assert on_submit(molecule_name) == expected_result

def test_on_submit_molecule_not_found():
    # Mocking the get_smiles function to return None
    molecule_name = "UnknownMolecule"
    assert on_submit(molecule_name) == None
