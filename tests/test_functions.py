import pytest

#test the get_smiles function
from Chromadvisor_without_interface import get_smiles

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

########################################################################################################
#test the find_functional_groups function
from Chromadvisor_without_interface import find_functional_groups

def test_find_functional_groups_carboxylic_acid():
    smiles = "CC(=O)O"
    result = find_functional_groups(smiles)
    assert result.get("carboxylic acid") == {'count': 1} # 'positions': [(1, 2, 3)]}

def test_find_functional_groups_amide():
    smiles = "CC(=O)NC"
    result = find_functional_groups(smiles)
    assert result.get("amide") == {'count': 1} # 'positions': [(2, 3, 4)]}

def test_find_functional_groups_no_functional_groups():
    smiles = "CCC"
    result = find_functional_groups(smiles)
    assert result == {}

def test_find_functional_groups_invalid_smiles():
    smiles = "InvalidSMILES"
    result = find_functional_groups(smiles)
    assert result is None

########################################################################################################
#test the on_submit function
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

########################################################################################################
#test the log and recommendation function
from Chromadvisor import calculate_logp_and_recommend_solvent

#test on chlorphenamine (apolar)
def test_recommend_solvent_apolaire():
    molecule_smiles = "c1([C@@H](c2ccccn2)CCN(C)C)ccc(Cl)cc1"
    result = recommend_solvent(molecule_smiles)
    assert result == "Use an apolar eluent such as : hexane or toluene with a bit of EtOH/acetone."

#test on ortho-Cresol
def test_recommend_solvent_apolar_polar():
    molecule_smiles = "Cc1ccccc1O"
    result = recommend_solvent(molecule_smiles)
    assert result == "Use a mix of DCM/MeOH or a mix of ethyl acetate/hexane."

#test on glucose (polar)
def test_recommend_solvent_polar():
    molecule_smiles = "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1"
    result = recommend_solvent(molecule_smiles)
    assert result == "Use a polar eluent such as : MeOH or acetone, even water."


