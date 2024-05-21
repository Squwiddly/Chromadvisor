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
#test the display_molecul_2d function
from tkinter import Tk
@pytest.fixture
def parent_window():
    # Create a parent window for testing
    parent = Tk()
    yield parent
    parent.destroy()

def test_display_molecule_2d_valid_smiles(parent_window, mocker):
    smiles = "CCO"
    with mocker.patch('tkinter.messagebox.showerror') as mock_showerror:
        display_molecule_2d(smiles, parent_window)
        assert mock_showerror.call_count == 0  # Ensure messagebox.showerror is not called
        # Ensure that the molecule image label is created and packed
        assert len(parent_window.winfo_children()) == 1  # Assuming no other widgets are present

def test_display_molecule_2d_invalid_smiles(parent_window, mocker):
    smiles = "invalid_smiles"
    with mocker.patch('tkinter.messagebox.showerror') as mock_showerror:
        display_molecule_2d(smiles, parent_window)
        # Ensure that messagebox.showerror is called with the appropriate error message
        mock_showerror.assert_called_once_with("Error", "Impossible to convert the SMILES into a molecule.")
        # Ensure that no molecule image label is created
        assert len(parent_window.winfo_children()) == 0


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



