"""structure 包 - 晶体结构处理层"""
from .cif_io import read_cif, structure_from_params, get_lattice_params, get_atom_info
from .symmetry_validator import SymmetryValidator
from .site_mapping import SiteMapper
from .coordinate_transform import (
    frac_to_cart, cart_to_frac,
    apply_affine_transform, wrap_to_unit_cell,
    build_supercell, transform_lattice, coordinates_are_equal,
)

__all__ = [
    "read_cif", "structure_from_params", "get_lattice_params", "get_atom_info",
    "SymmetryValidator",
    "SiteMapper",
    "frac_to_cart", "cart_to_frac",
    "apply_affine_transform", "wrap_to_unit_cell",
    "build_supercell", "transform_lattice", "coordinates_are_equal",
]
