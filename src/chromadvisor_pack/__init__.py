"""Recommending an eluant mixture based of the log(P) of the desired molecule"""
from .functions import (
    get_smiles,
    find_functional_groups,
    calculate_logp_and_recommend_solvent,
    display_molecule_2d,
    on_submit,
    generate_3d_structure
)
__version__ = "0.0.1"
