"""Tests for test_method3_affine_equivalence; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.utils.lattice import centering_primitive_matrix
from backend.utils.lattice import identity_matrix
from backend.utils.lattice import same_lattice
from tests.method3_affine_equivalence import parent_operations
from tests.support_backend_method3 import _i4mmm_parent


# Original coverage: test_method3_affine_validator.py.

def test_independent_validator_derives_conventional_centering_from_operations() -> None:
    parent = parent_operations(_i4mmm_parent(), symprec=1e-5)

    assert parent.space_group_number == 139
    assert same_lattice(parent.lattice, centering_primitive_matrix("I"))

def test_independent_validator_respects_uploaded_primitive_nonstandard_axes() -> None:
    primitive = _i4mmm_parent().get_primitive_structure()
    # A determinant-one shear changes the uploaded axes without changing the
    # physical primitive cell.  Its integer translations must therefore be
    # Z^3 even though spglib's standardized HM symbol still begins with I.
    primitive.make_supercell([[1, 1, 0], [0, 1, 0], [0, 0, 1]])

    parent = parent_operations(primitive, symprec=1e-5)

    assert parent.space_group_number == 139
    assert same_lattice(parent.lattice, identity_matrix())
