from __future__ import annotations

from fractions import Fraction

import pytest
from data_dir import experiment_data_dir
from pymatgen.core import Structure

from isocore.distortion.affine_embeddings import (
    AffineEmbedding,
    AffineOperation,
    ParentAffineGroup,
    build_affine_quotient,
    locate_embedding_subgroup,
    parent_affine_group,
)
from isocore.distortion.inverse_landau import FixedSpaceFeasibilityAnalyzer
from isocore.distortion.inverse_landau_adapter import (
    build_selected_character_representation,
    build_selected_representation,
    diagnose_embedding_feasibility,
    stable_embedding_id,
)
from isocore.utils.lattice import identity_matrix, multiply, rational_matrix


def _p2_parent() -> ParentAffineGroup:
    identity = identity_matrix()
    twofold = rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    return ParentAffineGroup(
        lattice=identity,
        operations=(
            AffineOperation(identity, zero),
            AffineOperation(twofold, zero),
        ),
        space_group_number=3,
        hall_number=3,
        symbol="P2",
    )


def test_affine_quotient_locates_and_validates_embedding_subgroup() -> None:
    parent = _p2_parent()
    child = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    quotient = build_affine_quotient(parent, child)
    embedding = AffineEmbedding(
        child,
        (
            AffineOperation(
                identity_matrix(),
                (Fraction(0), Fraction(0), Fraction(0)),
            ),
        ),
    )

    subgroup = locate_embedding_subgroup(embedding, quotient)

    assert len(quotient.operations) == 4
    assert quotient.translation_index == 2
    assert subgroup == (quotient.identity_index,)


def test_embedding_id_ignores_basis_representative_and_operation_order() -> None:
    parent = _p2_parent()
    child = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    operations = parent.operations
    first = AffineEmbedding(child, operations)
    expected = stable_embedding_id(parent, first)
    unimodular_transforms = (
        rational_matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]]),
        rational_matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
        rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        rational_matrix([[1, 1, 0], [1, 2, 0], [0, 0, 1]]),
    )

    for unimodular in unimodular_transforms:
        alternative_child = multiply(unimodular, child)
        second = AffineEmbedding(alternative_child, tuple(reversed(operations)))
        assert stable_embedding_id(parent, second) == expected


def test_embedding_id_keeps_distinct_origin_and_orientation_embeddings() -> None:
    parent = _p2_parent()
    child_x = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    child_z = rational_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    identity = AffineOperation(
        identity_matrix(),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    twofold = parent.operations[1]
    shifted_twofold = AffineOperation(
        twofold.rotation,
        (Fraction(1), Fraction(0), Fraction(0)),
    )
    origin_zero = AffineEmbedding(child_x, (identity, twofold))
    origin_shifted = AffineEmbedding(child_x, (identity, shifted_twofold))
    oriented = AffineEmbedding(child_z, (identity,))
    reference_orientation = AffineEmbedding(child_x, (identity,))

    assert stable_embedding_id(parent, origin_zero) != stable_embedding_id(
        parent, origin_shifted
    )
    assert stable_embedding_id(parent, reference_orientation) != stable_embedding_id(
        parent, oriented
    )


@pytest.mark.parametrize("filename", ["EuAl4 Parent.cif", "NdNiO2 own.cif"])
def test_real_parent_embedding_builds_exact_nonselectable_diagnostic(
    filename: str,
) -> None:
    structure = Structure.from_file(experiment_data_dir() / filename)
    parent = parent_affine_group(structure)
    embedding = AffineEmbedding(
        parent.lattice,
        parent.operations,
        parent.hall_number,
    )

    result = diagnose_embedding_feasibility(
        parent,
        embedding,
        structure,
        ("strain", "displacive"),
    )

    assert result.status == "embedding_feasible"
    assert result.route_status == "unresolved_irrep_decomposition"
    assert not result.selectable
    assert result.site_count is not None and result.site_count > 1
    assert result.max_site_match_error_angstrom is not None
    assert result.max_site_match_error_angstrom < 1e-7
    assert result.fixed_space is not None
    assert result.fixed_space.is_reachable


def test_real_ndnio2_strain_cannot_break_translation_but_displacement_can() -> None:
    structure = Structure.from_file(experiment_data_dir() / "NdNiO2 own.cif")
    parent = parent_affine_group(structure)
    child = multiply(
        rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]]),
        parent.lattice,
    )
    identity_embedding = AffineEmbedding(
        child,
        (
            AffineOperation(
                identity_matrix(),
                (Fraction(0), Fraction(0), Fraction(0)),
            ),
        ),
    )

    strain = diagnose_embedding_feasibility(
        parent,
        identity_embedding,
        structure,
        "strain",
    )
    displacive = diagnose_embedding_feasibility(
        parent,
        identity_embedding,
        structure,
        "displacive",
    )

    assert strain.status == "embedding_infeasible"
    assert strain.fixed_space is not None
    assert not strain.fixed_space.is_reachable
    assert displacive.status == "embedding_feasible"
    assert displacive.fixed_space is not None
    assert displacive.fixed_space.is_reachable


def test_real_selected_character_bundle_matches_dense_representation() -> None:
    structure = Structure.from_file(experiment_data_dir() / "NdNiO2 own.cif")
    parent = parent_affine_group(structure)
    child = multiply(
        rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]]),
        parent.lattice,
    )
    embedding = AffineEmbedding(
        child,
        (
            AffineOperation(
                identity_matrix(),
                (Fraction(0), Fraction(0), Fraction(0)),
            ),
        ),
    )
    quotient = build_affine_quotient(parent, child)
    subgroup = locate_embedding_subgroup(embedding, quotient)
    dense = build_selected_representation(
        parent, quotient, structure, ("strain", "displacive")
    )
    character = build_selected_character_representation(
        parent, quotient, structure, ("strain", "displacive")
    )

    assert dense.representation is not None
    assert character.representation is not None
    assert character.site_count == dense.site_count
    assert character.component_names == dense.component_names
    dense_result = FixedSpaceFeasibilityAnalyzer(dense.representation).analyze(subgroup)
    character_result = FixedSpaceFeasibilityAnalyzer(
        character.representation
    ).analyze(subgroup)
    assert character_result == dense_result


def test_rotational_and_magnetic_are_explicitly_unsupported() -> None:
    structure = Structure.from_file(experiment_data_dir() / "NdNiO2 own.cif")
    parent = parent_affine_group(structure)
    embedding = AffineEmbedding(parent.lattice, parent.operations)

    result = diagnose_embedding_feasibility(
        parent,
        embedding,
        structure,
        ("rotational", "magnetic"),
    )

    assert result.status == "unsupported_distortion_types"
    assert result.is_symmetry_reachable is None
    assert {item[0] for item in result.unsupported_types} == {
        "rotational",
        "magnetic",
    }
    assert not result.selectable
