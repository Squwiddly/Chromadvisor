"""Recommending an eluant mixture based of the log(P) of the desired molecule"""
from .functions_without_interface, .functions import (
    get_smiles,
    find_functional_groups,
    calculate_logp_and_recommend_solvent,
    display_molecule_2d,
    generate_3d_structure
)
__version__ = "0.0.1"
