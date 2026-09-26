from __future__ import annotations

from fractions import Fraction

import pytest
import spglib

from isocore.distortion.affine_embeddings import (
    AffineEmbedding,
    ParentAffineGroup,
    _hall_database_operations,
    _hall_setting_transforms,
    affine_equivalence,
    embedding_from_identity,
    standardized_embedding_identity,
)
from isocore.utils.lattice import identity_matrix


@pytest.mark.parametrize(
    ("source_hall", "preferred_hall"),
    ((408, 409), (402, 403)),
)
def test_standardization_converts_spglib_origin_choice_to_preferred_hall(
    source_hall: int,
    preferred_hall: int,
) -> None:
    """SG 129/125 choice 1 must become exact ISODISTORT choice 2."""

    source_type = spglib.get_spacegroup_type(source_hall)
    assert source_type is not None
    lattice = identity_matrix()
    operations = _hall_database_operations(source_hall)
    parent = ParentAffineGroup(
        lattice=lattice,
        operations=operations,
        space_group_number=int(source_type.number),
        hall_number=source_hall,
        symbol=str(source_type.international_short),
        standard_basis=lattice,
        standard_origin=(Fraction(0), Fraction(0), Fraction(0)),
    )
    embedding = AffineEmbedding(lattice, operations, source_hall)

    identity = standardized_embedding_identity(
        embedding,
        parent,
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 8.0)),
        int(source_type.number),
    )

    assert identity is not None
    assert "hall_number" not in identity
    reconstructed = embedding_from_identity(identity, parent)
    assert reconstructed.hall_number == preferred_hall
    assert (
        affine_equivalence(
            reconstructed,
            embedding,
            parent,
            allow_parent_conjugacy=False,
        )
        == "equal"
    )


def test_hall_setting_transform_supports_conventional_axis_relabelling() -> None:
    """The setting bridge is not restricted to origin-choice-only examples."""

    transforms = _hall_setting_transforms(4, 3)  # P 2 (c) -> P 2y (b)

    assert transforms
    assert any(matrix != identity_matrix() for matrix, _shift in transforms)


def test_explicit_hall_number_must_belong_to_requested_space_group() -> None:
    lattice = identity_matrix()
    operations = _hall_database_operations(408)
    parent = ParentAffineGroup(
        lattice=lattice,
        operations=operations,
        space_group_number=129,
        hall_number=408,
        symbol="P4/nmm",
    )

    with pytest.raises(ValueError, match="belongs to SG 125, not SG 129"):
        embedding_from_identity(
            {
                "space_group_number": 129,
                "symbol": "P4/nmm",
                "basis": lattice,
                "origin": (0, 0, 0),
                "hall_number": 403,
            },
            parent,
        )
