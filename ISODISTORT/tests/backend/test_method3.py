"""Tests for backend.method3; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.method3 import Method3Query
from backend.method3 import Method3Solver
from backend.method3 import search as search_methods_module
from backend.method3.coupled_routes import CoupledRouteCandidate
from backend.method3.coupled_routes import lift_embedding_subgroup
from backend.method3.coupled_routes import resolve_coupled_route_witness
from backend.method3.inverse_landau import FiniteGroup
from backend.method3.inverse_landau import FixedSpaceFeasibilityAnalyzer
from backend.method3.inverse_landau import RationalRepresentation
from backend.method3.inverse_landau import analyze_fixed_space
from backend.method3.inverse_landau import direct_sum_representations
from backend.method3.inverse_landau import homogeneous_strain_representation
from backend.method3.inverse_landau import polar_vector_representation
from backend.method3.inverse_landau import remove_uniform_site_vectors
from backend.method3.inverse_landau import site_displacement_representation
from backend.method3.inverse_landau import site_vector_character_representation
from backend.method3.inverse_landau_adapter import build_selected_character_representation
from backend.method3.inverse_landau_adapter import build_selected_representation
from backend.method3.inverse_landau_adapter import diagnose_embedding_feasibility
from backend.method3.inverse_landau_adapter import diagnose_embeddings_feasibility
from backend.method3.inverse_landau_adapter import stable_embedding_id
from backend.method3.search import _commensurate_line_parameters
from backend.method3.search import _eval_simple_k_component
from backend.method3.search import _kvec_from_template
from backend.method3.search import _line_kpoint_template
from backend.method3.search import _parse_affine_k_component
from backend.method3.search import _parse_basis_rows
from backend.method3.search import _reciprocal_vectors_equivalent
from backend.method3.search import _same_method3_route
from backend.method3.search import _selected_primitive_lattice
from backend.method3.search import _subgroup_primitive_lattice
from backend.method3.search import _validate_centering
from backend.utils.affine_embeddings import AffineEmbedding
from backend.utils.affine_embeddings import AffineOperation
from backend.utils.affine_embeddings import ParentAffineGroup
from backend.utils.affine_embeddings import build_affine_quotient
from backend.utils.affine_embeddings import locate_embedding_subgroup
from backend.utils.affine_embeddings import parent_affine_group
from backend.utils.config_loader import Config
from backend.utils.k_lattice import k_compatible_with_supercell
from backend.utils.k_lattice import param_value_candidates
from backend.utils.lattice import as_fraction
from backend.utils.lattice import centering_primitive_matrix
from backend.utils.lattice import identity_matrix
from backend.utils.lattice import multiply
from backend.utils.lattice import rational_matrix
from backend.utils.lattice import same_lattice
from backend.wrappers import IrrepInfo
from backend.wrappers import KPointInfo
from backend.wrappers import SubgroupInfo
from backend.wrappers import SubgroupInfo as _search_methods_SubgroupInfo
from data_dir import experiment_data_dir
from dataclasses import asdict
from dataclasses import replace
from fractions import Fraction
from pymatgen.core import Lattice
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from tests.support_backend_method3 import _DummyIsoWrapper
from tests.support_backend_method3 import _i4mmm_parent
from tests.support_backend_method3 import _method3_route
from types import SimpleNamespace
import numpy as np
import pytest


# Original coverage: test_coupled_routes.py.

ZERO = (Fraction(0), Fraction(0), Fraction(0))

def _operation(rotation):
    return AffineOperation(rational_matrix(rotation), ZERO)

def _route(label: str, space_group_number: int) -> SubgroupInfo:
    return SubgroupInfo(
        index=0,
        space_group_number=space_group_number,
        irrep_label=label,
        opd_symbol="P1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )

def _klein_parent():
    identity = _operation([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    first = _operation([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
    second = _operation([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    product = _operation([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    parent = ParentAffineGroup(
        lattice=identity_matrix(),
        operations=(identity, first, second, product),
        space_group_number=47,
        hall_number=227,
        symbol="Pmmm",
    )
    target = AffineEmbedding(identity_matrix(), (identity,))
    candidates = [
        CoupledRouteCandidate(
            _route("GM2+", 10),
            AffineEmbedding(identity_matrix(), (identity, first)),
        ),
        CoupledRouteCandidate(
            _route("GM4+", 10),
            AffineEmbedding(identity_matrix(), (identity, second)),
        ),
    ]
    return parent, target, candidates

def test_resolve_coupled_route_witness_uses_exact_stabilizer_intersection():
    parent, target, candidates = _klein_parent()

    witness = resolve_coupled_route_witness(parent, target, candidates)

    assert witness is not None
    assert [route.irrep_label for route in witness.components] == ["GM2+", "GM4+"]
    assert witness.target_subgroup_order == 1
    assert witness.quotient_order == 4
    assert witness.verification_method == "exact_affine_stabilizer_intersection"

def test_lift_embedding_subgroup_adds_translation_cosets():
    parent, _, _ = _klein_parent()
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    quotient = build_affine_quotient(parent, child_lattice)
    identity_route = AffineEmbedding(identity_matrix(), (parent.operations[0],))

    lifted = lift_embedding_subgroup(identity_route, quotient)

    assert len(lifted) == 2
    assert quotient.identity_index in lifted

def test_coupled_route_state_budget_fails_without_truncation():
    parent, target, candidates = _klein_parent()

    with pytest.raises(RuntimeError, match="No candidate result was truncated"):
        resolve_coupled_route_witness(
            parent,
            target,
            candidates,
            max_intersection_states=1,
        )

# Original coverage: test_inverse_landau.py.

def _cyclic_two() -> FiniteGroup:
    return FiniteGroup.from_operation(
        ("e", "a"),
        "e",
        lambda left, right: "e" if left == right else "a",
        name="C2",
    )

def _klein_four() -> FiniteGroup:
    elements = ((0, 0), (1, 0), (0, 1), (1, 1))
    return FiniteGroup.from_operation(
        elements,
        (0, 0),
        lambda left, right: (left[0] ^ right[0], left[1] ^ right[1]),
        name="C2 x C2",
    )

def test_trivial_representation_rejects_a_spurious_broken_subgroup() -> None:
    group = _cyclic_two()
    trivial = RationalRepresentation.from_matrices(
        group,
        {"e": ((1,),), "a": ((1,),)},
        name="trivial",
    )

    broken = analyze_fixed_space(trivial, ("e",))
    unbroken = analyze_fixed_space(trivial, ("e", "a"))

    assert broken.fixed_dimension == 1
    assert broken.pointwise_stabilizer == ("e", "a")
    assert not broken.is_reachable
    assert "outside H" in (broken.obstruction or "")
    assert unbroken.is_reachable

def test_arbitrary_matrices_still_require_full_exact_homomorphism() -> None:
    group = _cyclic_two()

    with pytest.raises(ValueError, match="do not form a representation"):
        RationalRepresentation.from_matrices(
            group,
            {"e": ((1,),), "a": ((2,),)},
            name="invalid arbitrary matrices",
        )

def test_homogeneous_strain_cannot_break_a_pure_translation() -> None:
    group = FiniteGroup.from_operation(
        ("identity", "translation"),
        "identity",
        lambda left, right: "identity" if left == right else "translation",
    )
    identity_rotation = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    strain = homogeneous_strain_representation(
        group,
        {
            "identity": identity_rotation,
            "translation": identity_rotation,
        },
    )

    analysis = analyze_fixed_space(strain, ("identity",))

    assert strain.dimension == 6
    assert analysis.fixed_dimension == 6
    assert analysis.pointwise_stabilizer == ("identity", "translation")
    assert not analysis.is_reachable

def test_polar_vector_has_no_inversion_even_fixed_space() -> None:
    group = FiniteGroup.from_operation(
        ("e", "inversion"),
        "e",
        lambda left, right: "e" if left == right else "inversion",
    )
    polar = polar_vector_representation(
        group,
        {
            "e": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            "inversion": ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
        },
    )

    polar_subgroup = analyze_fixed_space(polar, ("e",))
    centrosymmetric = analyze_fixed_space(polar, ("e", "inversion"))

    assert polar_subgroup.is_reachable
    assert polar_subgroup.pointwise_stabilizer == ("e",)
    assert centrosymmetric.projector == (
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert centrosymmetric.fixed_dimension == 0
    assert not centrosymmetric.is_reachable

def test_coupled_direct_sum_intersects_component_stabilizers() -> None:
    group = _klein_four()
    first = RationalRepresentation.from_matrices(
        group,
        {
            element: ((-1 if element[1] else 1,),)
            for element in group.elements
        },
        name="first character",
    )
    second = RationalRepresentation.from_matrices(
        group,
        {
            element: ((-1 if element[0] else 1,),)
            for element in group.elements
        },
        name="second character",
    )

    first_only = analyze_fixed_space(first, ((0, 0),))
    second_only = analyze_fixed_space(second, ((0, 0),))
    coupled = analyze_fixed_space(
        direct_sum_representations((first, second)),
        ((0, 0),),
    )

    assert set(first_only.pointwise_stabilizer) == {(0, 0), (1, 0)}
    assert set(second_only.pointwise_stabilizer) == {(0, 0), (0, 1)}
    assert not first_only.is_reachable
    assert not second_only.is_reachable
    assert coupled.fixed_dimension == 2
    assert coupled.pointwise_stabilizer == ((0, 0),)
    assert coupled.is_reachable

def test_character_analyzer_is_exactly_equivalent_to_dense_projectors() -> None:
    group = _klein_four()
    first = RationalRepresentation.from_matrices(
        group,
        {element: ((-1 if element[1] else 1,),) for element in group.elements},
    )
    second = RationalRepresentation.from_matrices(
        group,
        {element: ((-1 if element[0] else 1,),) for element in group.elements},
    )
    representation = direct_sum_representations((first, second))
    analyzer = FixedSpaceFeasibilityAnalyzer(representation)

    for subgroup in (
        ((0, 0),),
        ((0, 0), (1, 0)),
        group.elements,
    ):
        dense = analyze_fixed_space(representation, subgroup)
        summary = analyzer.analyze(subgroup)
        assert summary.fixed_dimension == dense.fixed_dimension
        assert summary.pointwise_stabilizer == dense.pointwise_stabilizer
        assert summary.is_reachable == dense.is_reachable
        assert summary.obstruction == dense.obstruction

def test_commensurate_site_displacements_can_break_a_translation() -> None:
    group = FiniteGroup.from_operation(
        ("e", "translation"),
        "e",
        lambda left, right: "e" if left == right else "translation",
    )
    identity_rotation = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    displacements = site_displacement_representation(
        group,
        {"e": (0, 1), "translation": (1, 0)},
        {"e": identity_rotation, "translation": identity_rotation},
    )

    broken = analyze_fixed_space(displacements, ("e",))
    translation_even = analyze_fixed_space(displacements, ("e", "translation"))

    assert broken.fixed_dimension == 6
    assert broken.pointwise_stabilizer == ("e",)
    assert broken.is_reachable
    assert translation_even.fixed_dimension == 3
    assert translation_even.projector[0][0] == Fraction(1, 2)
    assert translation_even.projector[0][3] == Fraction(1, 2)
    assert translation_even.is_reachable

def test_internal_site_character_matches_dense_quotient_exactly() -> None:
    group = FiniteGroup.from_operation(
        ("e", "translation"),
        "e",
        lambda left, right: "e" if left == right else "translation",
    )
    identity_rotation = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    permutations = {"e": (0, 1), "translation": (1, 0)}
    vector = polar_vector_representation(
        group,
        {"e": identity_rotation, "translation": identity_rotation},
    )
    full = site_displacement_representation(
        group,
        permutations,
        {"e": identity_rotation, "translation": identity_rotation},
    )
    dense_internal = remove_uniform_site_vectors(full, site_count=2)
    character_internal = site_vector_character_representation(
        group,
        permutations,
        vector,
        remove_uniform_vectors=True,
    )

    dense_analyzer = FixedSpaceFeasibilityAnalyzer(dense_internal)
    character_analyzer = FixedSpaceFeasibilityAnalyzer(character_internal)
    for subgroup in (("e",), group.elements):
        assert character_analyzer.analyze(subgroup) == dense_analyzer.analyze(subgroup)

# Original coverage: test_inverse_landau_adapter.py.

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

def test_batch_feasibility_matches_individual_character_criterion() -> None:
    structure = Structure.from_file(experiment_data_dir() / "NdNiO2 own.cif")
    parent = parent_affine_group(structure)
    identity = AffineOperation(
        identity_matrix(),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    embeddings = [
        AffineEmbedding(parent.lattice, parent.operations, parent.hall_number),
        AffineEmbedding(parent.lattice, (identity,)),
    ]

    batch = diagnose_embeddings_feasibility(
        parent,
        embeddings,
        structure,
        ("strain", "displacive"),
    )
    individual = [
        diagnose_embedding_feasibility(
            parent,
            embedding,
            structure,
            ("strain", "displacive"),
        )
        for embedding in embeddings
    ]

    assert [item.status for item in batch] == [item.status for item in individual]
    assert [item.embedding_id for item in batch] == [
        item.embedding_id for item in individual
    ]
    assert [item.fixed_space.is_reachable for item in batch] == [
        item.fixed_space.is_reachable for item in individual
    ]

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

def test_selected_species_displacement_keeps_relative_uniform_vector() -> None:
    structure = Structure.from_file(experiment_data_dir() / "NdNiO2 own.cif")
    parent = parent_affine_group(structure)
    quotient = build_affine_quotient(parent, parent.lattice)

    dense = build_selected_representation(
        parent,
        quotient,
        structure,
        "displacive",
        displacive_species=("Nd",),
    )
    character = build_selected_character_representation(
        parent,
        quotient,
        structure,
        "displacive",
        displacive_species=("Nd",),
    )

    assert dense.site_count == 1
    assert dense.representation is not None
    assert dense.representation.dimension == 3
    assert character.representation is not None
    assert character.representation.dimension == 3
    assert character.representation.characters == tuple(
        sum(matrix[index][index] for index in range(len(matrix)))
        for matrix in dense.representation.matrices
    )

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

# Original coverage: test_search_methods.py.

def test_method3_search_merges_same_operator_subgroup_across_settings() -> None:
    parent_structure = _i4mmm_parent()
    rotations = [
        operation.rotation_matrix.tolist()
        for operation in SpacegroupAnalyzer(parent_structure).get_symmetry_operations(
            cartesian=False
        )
    ]

    class _AlternativeSettingsIso(_DummyIsoWrapper):
        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            common = dict(
                space_group_number=71,
                space_group_symbol="Immm",
                subgroup_index=2,
                size=1,
                origin=[0, 0, 0],
                k_point_label="GM",
                irrep_label="GM1+",
                opd_symbol="P1",
            )
            return [
                _search_methods_SubgroupInfo(
                    index=0,
                    basis_vectors=[[1, 0, 0], [0, 0, -1], [0, 1, 0]],
                    **common,
                ),
                _search_methods_SubgroupInfo(
                    index=1,
                    basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    **common,
                ),
            ]

    items = Method3Solver(_AlternativeSettingsIso()).search(
        139,
        Method3Query(
            space_group_type=71,
            direct_sublattice_centering="d",
            parent_structure=parent_structure,
            parent_rotations=rotations,
            symmetry_tolerance=1e-5,
        ),
    )

    assert len(items) == 1
    assert len(items[0].routes or []) == 1

def test_method3_first_page_includes_affine_only_embeddings() -> None:
    parent_structure = Structure.from_spacegroup(
        "P2",
        Lattice.monoclinic(3.0, 4.0, 5.0, 100.0),
        ["H", "He"],
        [[0.123, 0.234, 0.345], [0.173, 0.272, 0.414]],
    )
    rotations = [
        operation.rotation_matrix.tolist()
        for operation in SpacegroupAnalyzer(parent_structure).get_symmetry_operations(
            cartesian=False
        )
    ]

    class _NoSingleIrRoutes(_DummyIsoWrapper):
        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            return []

    items = Method3Solver(_NoSingleIrRoutes()).search(
        3,
        Method3Query(
            space_group_type=3,
            supercell_basis=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            direct_sublattice_centering="P",
            parent_structure=parent_structure,
            parent_rotations=rotations,
            symmetry_tolerance=1e-5,
            include_affine_only_diagnostics=True,
        ),
    )

    assert len(items) == 2
    assert all(item.routes == [] for item in items)
    assert all(
        item.route_resolution == "affine_only_unresolved_coupled_route"
        for item in items
    )
    assert all(item.affine_embedding is not None for item in items)
    assert all(str(item.embedding_id).startswith("m3e1-") for item in items)
    assert len({item.embedding_id for item in items}) == len(items)

    default_items = Method3Solver(_NoSingleIrRoutes()).search(
        3,
        Method3Query(
            space_group_type=3,
            supercell_basis=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            direct_sublattice_centering="P",
            parent_structure=parent_structure,
            parent_rotations=rotations,
            symmetry_tolerance=1e-5,
        ),
    )
    assert default_items == []

def test_method3_route_dedup_recognizes_equivalent_k_star_parameters() -> None:
    first = _method3_route("L1", k_parameters=["1/3"])
    second = _method3_route("L1", k_parameters=["2/3"])
    first.k_coordinates = ["1/3", "0", "0"]
    second.k_coordinates = ["-1/3", "0", "0"]
    inversion = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]

    assert _same_method3_route(
        first, second, identity_matrix(), [inversion]
    )

def _i4mmm_coordinate_presentations() -> list[Structure]:
    conventional = _i4mmm_parent()
    primitive = conventional.get_primitive_structure()
    nonstandard = primitive.copy()
    nonstandard.make_supercell([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    return [conventional, primitive, nonstandard]

def _structure_rotations(structure: Structure) -> list[list[list[float]]]:
    return [
        operation.rotation_matrix.tolist()
        for operation in SpacegroupAnalyzer(
            structure, symprec=1e-5
        ).get_symmetry_operations(cartesian=False)
    ]

def test_method3_full_entry_uses_uploaded_parent_lattice_for_stage_a() -> None:
    class _NoRoutes(_DummyIsoWrapper):
        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            return []

    signatures: list[list[tuple[int, tuple[tuple[Fraction, ...], ...]]]] = []
    for structure in _i4mmm_coordinate_presentations():
        per_structure = []
        for target, centering, expected_resolution in (
            (139, "d", "affine_only_unresolved_coupled_route"),
            (123, "P", "affine_only_infeasible"),
        ):
            items = Method3Solver(_NoRoutes()).search(
                139,
                Method3Query(
                    space_group_type=target,
                    supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    direct_sublattice_centering=centering,
                    parent_structure=structure,
                    parent_rotations=_structure_rotations(structure),
                    symmetry_tolerance=1e-5,
                    include_affine_only_diagnostics=True,
                ),
            )
            assert len(items) == 1
            assert items[0].routes == []
            assert items[0].route_resolution == expected_resolution
            per_structure.append((
                items[0].subgroup.space_group_number,
                tuple(tuple(row) for row in items[0].subgroup.basis_vectors),
            ))
        signatures.append(per_structure)

    # Default I with identity basis is the full parent lattice; explicit P
    # with identity basis is its index-two primitive subgroup lattice.  Both
    # standard-frame selections must be independent of the uploaded cell axes.
    assert signatures[0] == signatures[1] == signatures[2]

def test_method3_parametric_congruence_uses_uploaded_parent_coordinates() -> None:
    class _LineIso(_DummyIsoWrapper):
        def __init__(self) -> None:
            self.requested_parameters: list[tuple[str, ...]] = []

        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            return [
                KPointInfo(
                    label="GM",
                    coordinates=["0", "0", "0"],
                    is_special=True,
                    parameters=[],
                ),
                KPointInfo(
                    label="LINE",
                    coordinates=["a", "a", "0"],
                    is_special=False,
                    parameters=["a"],
                ),
            ]

        def list_irreps(self, parent_sg, k_point, k_parameters=None):
            self.requested_parameters.append(tuple(k_parameters or []))
            return [IrrepInfo(label="L1", dimension=1)]

        @staticmethod
        def _include_irrep(irrep, distortion_types):
            return True

        def list_subgroups(self, *args, **kwargs):
            return [
                _search_methods_SubgroupInfo(
                    index=0,
                    space_group_number=123,
                    space_group_symbol="P4/mmm",
                    subgroup_index=18,
                    size=18,
                    opd_symbol="P1",
                    opd_dir_raw="(a)",
                    basis_vectors=[[3, 0, 0], [0, 3, 0], [0, 0, 1]],
                    origin=[0, 0, 0],
                    k_point_label="LINE",
                    irrep_label="L1",
                    parent_sg=139,
                )
            ]

    for structure in _i4mmm_coordinate_presentations():
        backend = _LineIso()
        items = Method3Solver(backend).search(
            139,
            Method3Query(
                space_group_type=123,
                supercell_basis=[[3, 0, 0], [0, 3, 0], [0, 0, 1]],
                direct_sublattice_centering="P",
                parent_structure=structure,
                parent_rotations=_structure_rotations(structure),
                symmetry_tolerance=1e-5,
            ),
        )
        # a=0 is the special Gamma endpoint; a=1/3 and 2/3 are one star under
        # inversion.  The same canonical backend query must result in every
        # uploaded coordinate presentation.
        assert backend.requested_parameters == [("1/3",)]
        assert len(items) == 1
        assert items[0].subgroup.k_parameters == ["1/3"]

def test_method_3_space_group_precedence():
    engine = Method3Solver(_DummyIsoWrapper())

    items = engine.search(
        225,
        Method3Query(
            point_group="m-3m",
            space_group_type=62,
            supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        ),
    )

    assert len(items) == 1
    assert items[0].subgroup.space_group_number == 62

def test_method_3_reuses_special_k_enumeration_for_repeated_filters():
    class _CountingIso(_DummyIsoWrapper):
        calls = 0

        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            self.calls += 1
            return super().enumerate_all_special_subgroups(parent_sg, distortion_types)

    backend = _CountingIso()
    engine = Method3Solver(backend)
    results = []
    for sg in (225, 62):
        results.append(engine.search(
            225,
            Method3Query(
                distortion_types=["strain", "displacive"],
                space_group_type=sg,
                supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            ),
        ))
    assert backend.calls == 1
    repeated = engine.search(
        225,
        Method3Query(
            distortion_types=["strain", "displacive"],
            space_group_type=225,
            supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        ),
    )
    assert backend.calls == 1
    assert results[0][0].subgroup is not repeated[0].subgroup

def test_method_3_preserves_all_official_centerings():
    assert _validate_centering("P") == "P"
    assert _validate_centering("p") == "P"
    assert _validate_centering("d") == "d"
    for letter in "ABCIFR":
        assert _validate_centering(letter) == letter

def test_method_3_default_and_p_define_distinct_centered_lattices():
    basis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    default, _ = _selected_primitive_lattice(
        139,
        Method3Query(space_group_type=71, direct_sublattice_centering="d"),
        basis,
    )
    primitive, _ = _selected_primitive_lattice(
        139,
        Method3Query(space_group_type=71, direct_sublattice_centering="P"),
        basis,
    )
    assert not same_lattice(default, primitive)

def test_method_3_identity_basis_is_exact_not_containment():
    class _IdentityIso(_DummyIsoWrapper):
        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            base = super().enumerate_all_special_subgroups(parent_sg, distortion_types)[2]
            larger = replace(
                base,
                index=3,
                basis_vectors=[[2, 0, 0], [0, 2, 0], [0, 0, 2]],
            )
            return [base, larger]

    items = Method3Solver(_IdentityIso()).search(
        225,
        Method3Query(
            space_group_type=62,
            supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            direct_sublattice_centering="P",
        ),
    )
    assert len(items) == 1
    assert items[0].subgroup.basis_vectors == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

def test_method_3_commensurate_ld_inference_helpers():
    """(0,0,6) includes every [0,1) lock-in, not only reciprocal divisors."""
    M = np.diag([1, 1, 6])
    assert M is not None
    cands = param_value_candidates(M)
    assert set(cands) == {
        Fraction(1, 6),
        Fraction(1, 3),
        Fraction(1, 2),
        Fraction(2, 3),
        Fraction(5, 6),
    }
    assert _commensurate_line_parameters(M, ["0", "0", "g"], "g") == [
        Fraction(0),
        Fraction(1, 6),
        Fraction(1, 3),
        Fraction(1, 2),
        Fraction(2, 3),
        Fraction(5, 6),
    ]
    k = _kvec_from_template(["0", "0", "g"], "g", Fraction(1, 6))
    assert k is not None
    assert k_compatible_with_supercell(k, M)
    assert not k_compatible_with_supercell(
        _kvec_from_template(["0", "0", "g"], "g", Fraction(1, 5)), M
    )

def test_exact_fraction_parsing_never_collapses_close_or_large_values():
    assert as_fraction("0.5000000001") == Fraction(5_000_000_001, 10_000_000_000)
    assert as_fraction("1/2000001") == Fraction(1, 2_000_001)
    assert as_fraction("2000001/2000000") == Fraction(2_000_001, 2_000_000)
    assert as_fraction("9876543210123456789/2000001") == Fraction(
        9_876_543_210_123_456_789,
        2_000_001,
    )
    basis = _parse_basis_rows([
        ["2000001/2000000", 0, 0],
        [0, "1/2000001", 0],
        [0, 0, "0.5000000001"],
    ])
    assert basis[0][0] == Fraction(2_000_001, 2_000_000)
    assert basis[1][1] == Fraction(1, 2_000_001)
    assert basis[2][2] == Fraction(5_000_000_001, 10_000_000_000)
    # A tolerance-based implementation incorrectly treats this as 1/2.
    assert not k_compatible_with_supercell(
        ("0.5000000001", "0", "0"),
        [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    assert k_compatible_with_supercell(
        ("1/2000001", "0", "0"),
        [[2_000_001, 0, 0], [0, 1, 0], [0, 0, 1]],
    )

def test_affine_k_parser_supports_website_minus_parameter_plus_one():
    assert _parse_affine_k_component("-a+1", "a") == (
        Fraction(-1),
        Fraction(1),
    )
    assert _parse_affine_k_component("1/2*a-3/7", "a") == (
        Fraction(1, 2),
        Fraction(-3, 7),
    )
    assert _eval_simple_k_component("-a+1", "a", Fraction(1, 3)) == Fraction(2, 3)
    assert _kvec_from_template(
        ["a", "-a+1", "0"],
        "a",
        Fraction(1, 3),
    ) == (Fraction(1, 3), Fraction(2, 3), Fraction(0))
    assert _parse_affine_k_component("a*b", "a") is None
    assert _parse_affine_k_component("b+1", "a") is None

    fallback = KPointInfo(
        label="Y",
        coordinates=["a", "-a+1", "0"],
        parameters=[],
        is_special=False,
    )
    assert _line_kpoint_template(1, fallback) == (["a", "-a+1", "0"], "a")

def test_affine_y_line_lockins_are_solved_exactly():
    values = _commensurate_line_parameters(
        [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
        ["a", "-a+1", "0"],
        "a",
    )
    assert values == [Fraction(0), Fraction(1, 2)]

def test_subgroup_primitive_lattice_prefers_actual_hm_centering():
    subgroup = _search_methods_SubgroupInfo(
        index=0,
        space_group_number=139,  # Standard setting is I, deliberately contradictory.
        space_group_symbol="P4/mmm",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    assert _subgroup_primitive_lattice(subgroup) == identity_matrix()

def test_method_3_parametric_search_uses_nonrepresentative_star_arm():
    class _StarArmIso(_DummyIsoWrapper):
        def __init__(self):
            self.requested_parameters = []

        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            return [
                KPointInfo(
                    label="LINE",
                    coordinates=["0", "a", "0"],
                    is_special=False,
                    parameters=["a"],
                )
            ]

        def list_irreps(self, parent_sg, k_point, k_parameters=None):
            self.requested_parameters.append(tuple(k_parameters or []))
            return [IrrepInfo(label="L1", dimension=1)]

        @staticmethod
        def _include_irrep(irrep, distortion_types):
            return True

        def list_subgroups(self, *args, **kwargs):
            return [
                _search_methods_SubgroupInfo(
                    index=0,
                    space_group_number=99,
                    space_group_symbol="P4mm",
                    subgroup_index=2,
                    size=2,
                    opd_symbol="P1",
                    opd_dir_raw="(a)",
                    basis_vectors=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                    origin=[0, 0, 0],
                    k_point_label="LINE",
                    irrep_label="L1",
                    parent_sg=123,
                )
            ]

    backend = _StarArmIso()
    rotation_swapping_x_y = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
    items = Method3Solver(backend).search(
        123,
        Method3Query(
            space_group_type=99,
            supercell_basis=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            parent_rotations=[rotation_swapping_x_y],
        ),
    )

    assert backend.requested_parameters == [("1/2",)]
    assert len(items) == 1
    assert items[0].subgroup.k_parameters == ["1/2"]

def test_method_3_deduplicates_symmetry_related_parameter_values_before_backend():
    class _EquivalentParameterIso(_DummyIsoWrapper):
        def __init__(self):
            self.requested_parameters = []

        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            return [
                KPointInfo(
                    label="X",
                    coordinates=["0", "1/2", "0"],
                    is_special=True,
                    parameters=[],
                ),
                KPointInfo(
                    label="Y",
                    coordinates=["a", "1/2", "0"],
                    is_special=False,
                    parameters=["a"],
                ),
            ]

        def list_irreps(self, parent_sg, k_point, k_parameters=None):
            self.requested_parameters.append(tuple(k_parameters or []))
            return [IrrepInfo(label="Y2", dimension=1)]

        @staticmethod
        def _include_irrep(irrep, distortion_types):
            return True

        def list_subgroups(self, *args, **kwargs):
            return [
                _search_methods_SubgroupInfo(
                    index=0,
                    space_group_number=51,
                    space_group_symbol="Pmma",
                    subgroup_index=12,
                    size=6,
                    opd_symbol="P2",
                    opd_dir_raw="(a)",
                    basis_vectors=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
                    origin=[0, 0, 0],
                    k_point_label="Y",
                    irrep_label="Y2",
                    parent_sg=123,
                )
            ]

    backend = _EquivalentParameterIso()
    inversion = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
    items = Method3Solver(backend).search(
        123,
        Method3Query(
            space_group_type=51,
            supercell_basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
            direct_sublattice_centering="P",
            parent_rotations=[inversion],
        ),
    )

    # Y(a=1/3) and Y(a=2/3) are related by inversion plus a reciprocal
    # lattice vector.  Querying both gives duplicate physical embeddings with
    # different backend origin representatives, so only the canonical first
    # exact value may reach list_irreps/list_subgroups.
    assert len(backend.requested_parameters) == 1
    assert len(items) == 1
    assert items[0].subgroup.k_parameters == ["1/3"]

def test_method_3_skips_parameter_endpoint_on_special_k_star():
    class _SpecialEndpointIso(_DummyIsoWrapper):
        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            return [
                KPointInfo(
                    label="SP",
                    coordinates=["0", "1/2", "0"],
                    is_special=True,
                    parameters=[],
                ),
                KPointInfo(
                    label="LINE",
                    coordinates=["a", "0", "0"],
                    is_special=False,
                    parameters=["a"],
                ),
            ]

        def list_irreps(self, parent_sg, k_point, k_parameters=None):
            raise AssertionError("special endpoint must not use the parametric backend")

    rotation_swapping_x_y = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
    items = Method3Solver(_SpecialEndpointIso()).search(
        123,
        Method3Query(
            space_group_type=99,
            supercell_basis=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            parent_rotations=[rotation_swapping_x_y],
        ),
    )

    assert items == []

def test_reciprocal_equivalence_respects_centered_parent_translations():
    parent_primitive = centering_primitive_matrix("I")
    gamma = (Fraction(0), Fraction(0), Fraction(0))
    # Odd h+k+l is not a reciprocal vector of an I-centered direct lattice.
    assert not _reciprocal_vectors_equivalent(
        (Fraction(1), Fraction(1), Fraction(1)),
        gamma,
        parent_primitive,
    )
    # Even h+k+l has integer phase on every I-parent primitive translation.
    assert _reciprocal_vectors_equivalent(
        (Fraction(1), Fraction(1), Fraction(0)),
        gamma,
        parent_primitive,
    )

def test_method_3_parametric_backend_failure_is_not_reported_as_empty():
    class _BrokenParametricIso(_DummyIsoWrapper):
        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            raise OSError("backend unavailable")

    try:
        Method3Solver(_BrokenParametricIso()).search(
            123,
            Method3Query(
                space_group_type=99,
                supercell_basis=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            ),
        )
        raise AssertionError("expected RuntimeError")
    except RuntimeError as exc:
        assert "failed to list k points" in str(exc)

def test_method_3_genuine_parametric_subgroup_failure_still_raises():
    class _BrokenSubgroupIso(_DummyIsoWrapper):
        def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
            return []

        def list_k_points(self, parent_sg):
            return [
                KPointInfo(
                    label="GM",
                    coordinates=["0", "0", "0"],
                    is_special=True,
                    parameters=[],
                ),
                KPointInfo(
                    label="LINE",
                    coordinates=["a", "0", "0"],
                    is_special=False,
                    parameters=["a"],
                ),
            ]

        def list_irreps(self, parent_sg, k_point, k_parameters=None):
            return [IrrepInfo(label="L1", dimension=1)]

        @staticmethod
        def _include_irrep(irrep, distortion_types):
            return True

        def list_subgroups(self, *args, **kwargs):
            raise OSError("genuine parametric backend failure")

    try:
        Method3Solver(_BrokenSubgroupIso()).search(
            123,
            Method3Query(
                space_group_type=99,
                supercell_basis=[[3, 0, 0], [0, 1, 0], [0, 0, 1]],
            ),
        )
        raise AssertionError("expected RuntimeError")
    except RuntimeError as exc:
        assert "failed to list subgroups" in str(exc)
        assert "1/3" in str(exc)

class _Method3BudgetIso(_DummyIsoWrapper):
    """Three exact parameter searches, each with two active IRs."""

    def __init__(self) -> None:
        self.irrep_calls = 0
        self.subgroup_calls = 0

    def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
        return []

    def list_k_points(self, parent_sg):
        return [
            KPointInfo(
                label="GM",
                coordinates=["0", "0", "0"],
                is_special=True,
                parameters=[],
            ),
            KPointInfo(
                label="LINE",
                coordinates=["a", "0", "0"],
                is_special=False,
                parameters=["a"],
            ),
        ]

    def list_irreps(self, parent_sg, k_point, k_parameters=None):
        self.irrep_calls += 1
        return [
            IrrepInfo(label="L1", dimension=1),
            IrrepInfo(label="L2", dimension=1),
        ]

    @staticmethod
    def _include_irrep(irrep, distortion_types):
        return True

    def list_subgroups(self, *args, **kwargs):
        self.subgroup_calls += 1
        return []

def _run_method3_budget_search(backend: _Method3BudgetIso):
    return Method3Solver(backend).search(
        123,
        Method3Query(
            space_group_type=99,
            supercell_basis=[[4, 0, 0], [0, 1, 0], [0, 0, 1]],
            direct_sublattice_centering="P",
        ),
    )

def test_method_3_parametric_value_budget_fails_before_list_irreps(monkeypatch):
    backend = _Method3BudgetIso()
    monkeypatch.setattr(
        search_methods_module,
        "get_config",
        lambda: SimpleNamespace(
            method3_max_parametric_values=2,
            method3_max_backend_queries=0,
        ),
    )

    with pytest.raises(RuntimeError, match=r"3 parametric values.*method3_max_parametric_values=2"):
        _run_method3_budget_search(backend)

    assert backend.irrep_calls == 0
    assert backend.subgroup_calls == 0

def test_method_3_million_index_budget_stops_lazy_candidates_before_backend(
    monkeypatch,
):
    backend = _Method3BudgetIso()
    monkeypatch.setattr(
        search_methods_module,
        "get_config",
        lambda: SimpleNamespace(
            method3_max_parametric_values=2,
            method3_max_backend_queries=0,
        ),
    )
    original = search_methods_module._iter_line_congruence_solutions
    yielded = 0

    def observed_solutions(*args, **kwargs):
        nonlocal yielded
        for value in original(*args, **kwargs):
            yielded += 1
            if yielded > 8:
                raise AssertionError("candidate generator consumed past the proven limit")
            yield value

    monkeypatch.setattr(
        search_methods_module,
        "_iter_line_congruence_solutions",
        observed_solutions,
    )

    with pytest.raises(
        RuntimeError,
        match=r"at least 3 parametric values.*method3_max_parametric_values=2",
    ):
        Method3Solver(backend).search(
            123,
            Method3Query(
                space_group_type=99,
                supercell_basis=[[1_000_000, 0, 0], [0, 1, 0], [0, 0, 1]],
                direct_sublattice_centering="P",
            ),
        )

    # Gamma is skipped, then the first three non-special values prove excess.
    assert yielded == 4
    assert backend.irrep_calls == 0
    assert backend.subgroup_calls == 0

def test_method_3_backend_query_lower_bound_fails_before_list_irreps(monkeypatch):
    backend = _Method3BudgetIso()
    monkeypatch.setattr(
        search_methods_module,
        "get_config",
        lambda: SimpleNamespace(
            method3_max_parametric_values=0,
            method3_max_backend_queries=2,
        ),
    )

    with pytest.raises(RuntimeError, match=r"at least 3 backend queries.*method3_max_backend_queries=2"):
        _run_method3_budget_search(backend)

    assert backend.irrep_calls == 0
    assert backend.subgroup_calls == 0

def test_method_3_exact_query_budget_fails_before_list_subgroups(monkeypatch):
    backend = _Method3BudgetIso()
    monkeypatch.setattr(
        search_methods_module,
        "get_config",
        lambda: SimpleNamespace(
            method3_max_parametric_values=0,
            method3_max_backend_queries=8,
        ),
    )

    with pytest.raises(RuntimeError, match=r"9 backend queries \(3 list_irreps \+ 6 list_subgroups\)"):
        _run_method3_budget_search(backend)

    assert backend.irrep_calls == 3
    assert backend.subgroup_calls == 0

def test_method_3_zero_cost_limits_are_unlimited(monkeypatch):
    backend = _Method3BudgetIso()
    monkeypatch.setattr(
        search_methods_module,
        "get_config",
        lambda: SimpleNamespace(
            method3_max_parametric_values=0,
            method3_max_backend_queries=0,
        ),
    )

    assert _run_method3_budget_search(backend) == []
    assert backend.irrep_calls == 3
    assert backend.subgroup_calls == 6

@pytest.mark.parametrize(
    "key",
    ["method3_max_parametric_values", "method3_max_backend_queries"],
)
def test_method_3_negative_cost_limit_is_configuration_error(key):
    config = object.__new__(Config)
    config._cfg = {"runtime": {key: -1}}

    with pytest.raises(ValueError, match=rf"runtime\.{key}=-1.*0 = unlimited"):
        getattr(config, key)

def test_method_3_restores_official_k_active_parameter_scale() -> None:
    class _ParametricIso(_DummyIsoWrapper):
        def list_k_points(self, parent_sg):
            return [
                KPointInfo(
                    label="LD",
                    coordinates=["0", "0", "a"],
                    is_special=False,
                    parameters=["a"],
                )
            ]

        def list_irreps(self, parent_sg, k_point, k_parameters=None):
            return [IrrepInfo(label="LD1", dimension=2)]

        @staticmethod
        def _include_irrep(irrep, distortion_types):
            return True

        def list_subgroups(self, *args, **kwargs):
            return [
                _search_methods_SubgroupInfo(
                    index=0,
                    space_group_number=99,
                    space_group_symbol="P4mm",
                    subgroup_index=24,
                    size=12,
                    opd_symbol="C1",
                    opd_dir_raw="(a,b)",
                    basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
                    origin=[0, 0, 0],
                    k_point_label="LD",
                    irrep_label="LD1",
                    parent_sg=139,
                    # Simulate the bundled iso's internal LD parameter scale.
                    k_active_raw=" (0,0,1/12)",
                )
            ]

    items = Method3Solver(_ParametricIso()).search(
        139,
        Method3Query(
            space_group_type=99,
            supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
        ),
    )

    routes = [route for item in items for route in (item.routes or [item.subgroup])]
    by_parameter = {
        tuple(route.k_parameters): route.official_fields()["k_active"]
        for route in routes
    }
    assert by_parameter[("1/6",)] == "(0,0,1/6)"
    assert "(0,0,1/12)" not in by_parameter.values()

# Original coverage: test_session_result_contract.py.

def test_method3_identity_has_one_serializable_source():
    from backend.method3 import Method3ResultItem
    subgroup = SubgroupInfo(index=0, space_group_number=1,
                            method3_route_resolution="exact_fixed_space",
                            method3_embedding_id="m3e1-" + "a" * 64)
    item = Method3ResultItem(subgroup, "1", [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    restored = SubgroupInfo(**asdict(subgroup))
    assert restored.method3_embedding_id == item.embedding_id
    assert item.route_resolution == restored.method3_route_resolution
    with pytest.raises(AttributeError):
        item.embedding_id = "conflicting copy"
    assert "embedding_id" not in vars(item)
