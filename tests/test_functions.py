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
#test the function display_molecule_2d
from Code_without_interface import display_molecule_2d

def test_display_molecule_2d_found():
    name_molecule = "123"
    result = display_molecule_2d(name_molecule)
    assert result == None

########################################################################################################
#test the find_functional_groups function
from Chromadvisor_without_interface import find_functional_groups

def test_find_functional_groups_salicylic_acid():
    smiles = "C1=CC=C(C(=C1)C(=O)O)O"
    result = find_functional_groups(smiles)
    assert result == {'carboxylic acid': {'count': 1, 'positions': (((6, 7, 8),))},
                      'phenol': {'count': 2, 'positions': ((6, 3, 2), (6, 3, 4))}}

def test_find_functional_groups_propionamide():
    smiles = "CCC(=O)N"
    result = find_functional_groups(smiles)
    assert result.get("amide") == {'count': 1, 'positions': ((4, 2, 3, 1),)}

def test_find_functional_groups_no_functional_groups():
    smiles = "CCC"
    result = find_functional_groups(smiles)
    assert result == {}

def test_find_functional_groups_invalid_smiles():
    smiles = "InvalidSMILES"
    result = find_functional_groups(smiles)
    assert result is None


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


