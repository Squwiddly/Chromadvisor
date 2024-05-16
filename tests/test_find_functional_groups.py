
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
