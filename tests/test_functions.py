import pytest

#test the get_smiles function
from chromadvisor_pack import get_smiles

def test_get_smiles_found():
    molecule_name = "water"
    result = get_smiles(molecule_name)
    assert result == "O"
  
def test_get_smiles_not_found():
    molecule_name = "nonexistent_molecule"
    result = get_smiles(molecule_name)
    assert result is None

def test_get_smiles_error(mocker):
    # Mock the get_compounds function to raise an exception
    mock_get_compounds = mocker.patch('pubchempy.get_compounds', side_effect=Exception("Test Exception"))
    result = get_smiles("invalid_molecule_name")
    assert result is None

########################################################################################################
#test the find_functional_groups function
from chromadvisor_pack import find_functional_groups

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
#test the display_molecule_2d function
from chromadvisor_pack import display_molecule_2d
from unittest.mock import patch
def test_display_molecule_2d():
    # Mock the functions used by display_molecule_2d
    with patch('your_module.Chem.MolFromSmiles') as mock_MolFromSmiles, \
         patch('your_module.Draw.MolToImage') as mock_MolToImage, \
         patch('your_module.ImageTk.PhotoImage') as mock_PhotoImage, \
         patch('your_module.ttk.Label') as mock_Label:

        # Set the return values for the mocked functions
        mock_MolFromSmiles.return_value = "mock_molecule"
        mock_MolToImage.return_value = "mock_image"
        mock_PhotoImage.return_value = "mock_photo_image"

        # Call the function with a valid SMILES string
        display_molecule_2d("CCO", "mock_parent_window")

        # Assert that the functions were called with the correct arguments
        mock_MolFromSmiles.assert_called_once_with("CCO")
        mock_MolToImage.assert_called_once_with("mock_molecule")
        mock_PhotoImage.assert_called_once_with("mock_image")
        mock_Label.assert_called_once_with("mock_parent_window", image="mock_photo_image")


########################################################################################################
#test the log and recommendation function
from chromadvisor_pack import calculate_logp_and_recommend_solvent

#test on chlorphenamine (apolar)
def test_recommend_solvent_apolaire():
    molecule_smiles = "c1([C@@H](c2ccccn2)CCN(C)C)ccc(Cl)cc1"
    result = calculate_logp_and_recommend_solvent(molecule_smiles)
    assert result == (3.8186000000000027, "Use an apolar eluent such as : hexane or toluene with a bit of EtOH/acetone.")

#test on ortho-Cresol
def test_recommend_solvent_apolar_polar():
    molecule_smiles = "Cc1ccccc1O"
    result = calculate_logp_and_recommend_solvent(molecule_smiles)
    assert result == (1.70062, "Use a mix of DCM/MeOH or a mix of ethyl acetate/hexane.")

#test on glucose (polar)
def test_recommend_solvent_polar():
    molecule_smiles = "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1"
    result = calculate_logp_and_recommend_solvent(molecule_smiles)
    assert result == (-3.2214000000000005, 'Use a polar eluent such as : MeOH or acetone, even water.')

def test_recommend_solvant_invalid():
    molecule_smiles = "InvalidSMILES"
    result = calculate_logp_and_recommend_solvent(molecule_smiles)
    assert result == (None, "Invalid SMILES")



