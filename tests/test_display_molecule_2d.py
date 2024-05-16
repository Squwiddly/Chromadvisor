
from Code_without_interface import display_molecule_2d

def test_display_molecule_2d_found():
    name_molecule = "123"
    result = display_molecule_2d(name_molecule)
    assert result == "Error", "Impossible to convert the SMILES into a molecule."
