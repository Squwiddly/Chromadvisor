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
from unittest import mock
from unittest.mock import MagicMock
from chromadvisor_pack.functions import display_molecule_2d
from tkinter import Tk
from PIL import Image, ImageTk
import io
# Importez les modules nécessaires pour le test
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except ImportError:
    Chem = None
    Draw = None

try:
    from PIL import ImageTk
except ImportError:
    ImageTk = None

def create_test_image():
    # Créer une petite image de test en mémoire
    img = Image.new('RGB', (10, 10), color = 'red')
    return ImageTk.PhotoImage(img)

def test_display_molecule_2d_success(mocker):
    smiles = 'CCO'
    parent_window = Tk()

    # Mocking Chem.MolFromSmiles
    mock_mol = MagicMock()
    mocker.patch('rdkit.Chem.MolFromSmiles', return_value=mock_mol)

    # Mocking Draw.MolToImage to return an image
    mock_image = create_test_image()
    mocker.patch('rdkit.Chem.Draw.MolToImage', return_value=mock_image)

    # Mocking ImageTk.PhotoImage to use the real image
    mocker.patch('PIL.ImageTk.PhotoImage', return_value=mock_image)

    display_molecule_2d(smiles, parent_window)

    Chem.MolFromSmiles.assert_called_once_with(smiles)
    Draw.MolToImage.assert_called_once_with(mock_mol)
    ImageTk.PhotoImage.assert_called_once_with(mock_image)
    
    # Cleanup
    parent_window.destroy()

def test_display_molecule_2d_invalid_smiles(mocker):
    smiles = 'invalid_smiles'
    parent_window = Tk()

    # Mocking Chem.MolFromSmiles to return None
    mocker.patch('rdkit.Chem.MolFromSmiles', return_value=None)

    # Mocking messagebox.showerror
    mock_showerror = mocker.patch('tkinter.messagebox.showerror')

    display_molecule_2d(smiles, parent_window)

    Chem.MolFromSmiles.assert_called_once_with(smiles)
    mock_showerror.assert_called_once_with("Error", "Impossible to convert the SMILES into a molecule.")

    # Cleanup
    parent_window.destroy()


########################################################################################################
#test the generate_3d_structure function
from unittest import mock
from chromadvisor_pack.functions import generate_3d_structure
from rdkit import Chem

def test_generate_3d_structure_success(mocker):
    smiles = 'CCO'
    # Mocking Chem.MolFromSmiles to return a molecule
    mocker.patch('rdkit.Chem.MolFromSmiles', return_value=Chem.MolFromSmiles(smiles))
    # Mocking py3Dmol.view
    mock_viewer = mocker.Mock()
    mocker.patch('py3Dmol.view', return_value=mock_viewer)
    result = generate_3d_structure(smiles)
    assert result == mock_viewer.show.return_value

def test_generate_3d_structure_invalid_smiles(mocker):
    smiles = 'invalid_smiles'
    # Mocking Chem.MolFromSmiles to return None
    mocker.patch('rdkit.Chem.MolFromSmiles', return_value=None)
    result = generate_3d_structure(smiles)
    assert result is None


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



