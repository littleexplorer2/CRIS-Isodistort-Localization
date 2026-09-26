from __future__ import annotations

from pymatgen.core import Lattice, Structure

from isocore.utils.lattice import (
    centering_primitive_matrix,
    identity_matrix,
    same_lattice,
)
from tests_dev.manual.method3_affine_equivalence import parent_operations


def _i4mmm_structure() -> Structure:
    return Structure.from_spacegroup(
        "I4/mmm",
        Lattice.tetragonal(4.0, 10.0),
        ["Eu"],
        [[0, 0, 0]],
    )


def test_independent_validator_derives_conventional_centering_from_operations() -> None:
    parent = parent_operations(_i4mmm_structure(), symprec=1e-5)

    assert parent.space_group_number == 139
    assert same_lattice(parent.lattice, centering_primitive_matrix("I"))


def test_independent_validator_respects_uploaded_primitive_nonstandard_axes() -> None:
    primitive = _i4mmm_structure().get_primitive_structure()
    # A determinant-one shear changes the uploaded axes without changing the
    # physical primitive cell.  Its integer translations must therefore be
    # Z^3 even though spglib's standardized HM symbol still begins with I.
    primitive.make_supercell([[1, 1, 0], [0, 1, 0], [0, 0, 1]])

    parent = parent_operations(primitive, symprec=1e-5)

    assert parent.space_group_number == 139
    assert same_lattice(parent.lattice, identity_matrix())
