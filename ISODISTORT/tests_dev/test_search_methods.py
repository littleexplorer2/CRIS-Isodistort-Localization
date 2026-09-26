"""
测试：ISODISTORT Method 1-4 搜索流程（使用桩对象，不依赖真实二进制）。
"""

from copy import deepcopy
from dataclasses import dataclass, replace
from fractions import Fraction
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.api import IsoDistort
from isocore.api import core_api as core_api_module
from isocore.backend import (
    BushMode,
    DistortionMode,
    IrrepInfo,
    IsoWrapper,
    KPointInfo,
    SubgroupInfo,
)
from isocore.distortion import (
    AffineEmbedding,
    AffineOperation,
    IsoSearchEngine,
    Method1Query,
    Method1ResultItem,
    Method2Query,
    Method3Query,
    Method3ResultItem,
    Method4Query,
    ParentAffineGroup,
    affine_equivalence,
    embedding_from_identity,
    enumerate_lifted_point_subgroups,
    enumerate_target_affine_embeddings,
    parent_affine_group,
)
from isocore.distortion import affine_embeddings as affine_embeddings_module
from isocore.distortion import search_methods as search_methods_module
from isocore.distortion import superspace as superspace_module
from isocore.distortion.search_methods import (
    _basis_is_same_lattice,
    _basis_is_sublattice_of,
    _commensurate_line_parameters,
    _eval_simple_k_component,
    _integer_basis_matrix,
    _k_compatible_with_supercell,
    _kvec_from_template,
    _line_kpoint_template,
    _param_value_candidates,
    _parse_affine_k_component,
    _parse_basis_rows,
    _reciprocal_vectors_equivalent,
    _same_method3_route,
    _selected_primitive_lattice,
    _subgroup_primitive_lattice,
    _validate_centering,
)
from isocore.utils.config_loader import Config
from isocore.utils.lattice import (
    as_fraction,
    centering_primitive_matrix,
    identity_matrix,
    inverse,
    multiply,
    rational_matrix,
    same_lattice,
)


def _i4mmm_parent() -> Structure:
    return Structure.from_spacegroup(
        "I4/mmm",
        Lattice.tetragonal(4.0, 10.0),
        ["Eu"],
        [[0, 0, 0]],
    )


def test_method3_api_propagates_physical_tolerance_model(monkeypatch) -> None:
    # This test exercises only the API/configuration boundary; no external
    # ISO executable or WSL session is part of the tolerance model.
    monkeypatch.setattr(core_api_module, "IsoWrapper", SimpleNamespace)
    monkeypatch.setattr(core_api_module, "SmodesWrapper", SimpleNamespace)
    api = IsoDistort()
    api.structure = _i4mmm_parent()
    api.symmetry_info = {"space_group_number": 139}
    captured = {}

    def _capture(parent_sg, query):
        captured["parent_sg"] = parent_sg
        captured["query"] = query
        return []

    monkeypatch.setattr(api._search, "method_3_search", _capture)
    assert api.search_method_3(space_group_type=139) == []

    cfg = api.cfg
    query = captured["query"]
    assert captured["parent_sg"] == 139
    assert query.symmetry_tolerance == cfg.symmetry_cartesian_tolerance_angstrom
    assert (
        query.symmetry_angle_tolerance_degrees
        == cfg.symmetry_angle_tolerance_degrees
    )
    assert (
        query.affine_exact_tolerance
        == cfg.affine_exact_cartesian_tolerance_angstrom
    )
    assert (
        query.fractional_coordinate_tolerance
        == cfg.fractional_coordinate_tolerance
    )


def test_exact_affine_spglib_boundary_uses_configured_angstrom_tolerance(
    monkeypatch,
) -> None:
    captured = {}

    def _identify(rotations, translations, lattice=None, symprec=None):
        captured["symprec"] = symprec
        return SimpleNamespace(number=1)

    monkeypatch.setattr(
        affine_embeddings_module.spglib,
        "get_spacegroup_type_from_symmetry",
        _identify,
    )
    identity = identity_matrix()
    embedding = AffineEmbedding(
        lattice=identity,
        operations=(
            AffineOperation(
                identity,
                (Fraction(0), Fraction(0), Fraction(0)),
            ),
        ),
    )

    assert (
        affine_embeddings_module.identify_embedding_space_group(
            embedding, np.eye(3)
        )
        == 1
    )
    assert captured["symprec"] == (
        Config().affine_exact_cartesian_tolerance_angstrom
    )

    assert (
        affine_embeddings_module.identify_embedding_space_group(
            embedding, np.eye(3) * 1e10
        )
        == 1
    )
    expected_roundoff_guard = 256 * np.finfo(float).eps * 1e10
    assert captured["symprec"] == pytest.approx(expected_roundoff_guard)


def test_method3_exact_affine_identity_accepts_alternative_setting() -> None:
    parent = parent_affine_group(_i4mmm_parent(), 1e-5)
    first = {
        "space_group_number": 71,
        "symbol": "Immm",
        "basis": [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
        "origin": [0, 0, 0],
    }
    second = {
        "space_group_number": 71,
        "symbol": "Immm",
        "basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "origin": [0, 0, 0],
    }

    assert affine_equivalence(
        embedding_from_identity(first, parent),
        embedding_from_identity(second, parent),
        parent,
        allow_parent_conjugacy=False,
    ) == "equal"


def test_parent_affine_group_uses_uploaded_cell_translation_lattice() -> None:
    conventional = _i4mmm_parent()
    conventional_group = parent_affine_group(conventional, 1e-5)
    assert same_lattice(
        conventional_group.lattice,
        centering_primitive_matrix("I"),
    )

    primitive = conventional.get_primitive_structure()
    primitive_group = parent_affine_group(primitive, 1e-5)
    assert primitive_group.space_group_number == 139
    assert same_lattice(primitive_group.lattice, identity_matrix())

    nonstandard = conventional.copy()
    shear = rational_matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    nonstandard.make_supercell([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    nonstandard_group = parent_affine_group(nonstandard, 1e-5)
    expected = multiply(centering_primitive_matrix("I"), inverse(shear))
    assert nonstandard_group.space_group_number == 139
    assert same_lattice(nonstandard_group.lattice, expected)


def test_parent_conjugacy_does_not_erase_embedding_multiplicity() -> None:
    parent = parent_affine_group(_i4mmm_parent(), 1e-5)
    first = embedding_from_identity(
        {
            "space_group_number": 140,
            "symbol": "I4/mcm",
            "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
            "origin": ["-1/2", "-1/2", 0],
        },
        parent,
    )
    second = embedding_from_identity(
        {
            "space_group_number": 140,
            "symbol": "I4/mcm",
            "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
            "origin": [0, 0, "1/2"],
        },
        parent,
    )

    assert affine_equivalence(
        first, second, parent, allow_parent_conjugacy=False
    ) is None
    assert affine_equivalence(first, second, parent) == "parent_conjugate"
    embeddings = [first, second]
    assert len(embeddings) == 2


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
                SubgroupInfo(
                    index=0,
                    basis_vectors=[[1, 0, 0], [0, 0, -1], [0, 1, 0]],
                    **common,
                ),
                SubgroupInfo(
                    index=1,
                    basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    **common,
                ),
            ]

    items = IsoSearchEngine(_AlternativeSettingsIso()).method_3_search(
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


def test_method3_exact_affine_identity_finds_parent_conjugacy() -> None:
    parent = parent_affine_group(_i4mmm_parent(), 1e-5)
    first = {
        "space_group_number": 140,
        "symbol": "I4/mcm",
        "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
        "origin": ["-1/2", "-1/2", 0],
    }
    second = {
        "space_group_number": 140,
        "symbol": "I4/mcm",
        "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
        "origin": [0, 0, "1/2"],
    }

    assert affine_equivalence(
        embedding_from_identity(first, parent),
        embedding_from_identity(second, parent),
        parent,
    ) == "parent_conjugate"


def test_method3_finite_quotient_finds_embeddings_without_ir_routes() -> None:
    """Stage A is group-theoretic: no single-IR route is needed to find lifts."""

    identity = identity_matrix()
    twofold = rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(
            AffineOperation(identity, zero),
            AffineOperation(twofold, zero),
        ),
        space_group_number=3,
        hall_number=3,
        symbol="P2",
    )
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])

    embeddings = enumerate_lifted_point_subgroups(
        parent, child_lattice, [[identity, twofold]]
    )

    assert len(embeddings) == 2
    assert all(len(embedding.operations) == 2 for embedding in embeddings)
    assert all(isinstance(embedding, AffineEmbedding) for embedding in embeddings)


def test_method3_affine_lifts_are_solved_from_point_group_generators() -> None:
    """A cyclic point group has one lift variable, not one per rotation."""

    identity = identity_matrix()
    quarter_turn = rational_matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    half_turn = multiply(quarter_turn, quarter_turn)
    inverse_quarter = multiply(half_turn, quarter_turn)
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=tuple(
            AffineOperation(rotation, zero)
            for rotation in (identity, quarter_turn, half_turn, inverse_quarter)
        ),
        space_group_number=75,
        hall_number=349,
        symbol="P4",
    )
    child_lattice = rational_matrix([[4, 0, 0], [0, 4, 0], [0, 0, 1]])

    # There are 16 lifts of the single C4 generator.  The former per-rotation
    # Cartesian estimate was 16**3 and incorrectly rejected this exact search
    # at the same boundary before solving the factor-set constraints.
    embeddings = enumerate_lifted_point_subgroups(
        parent,
        child_lattice,
        [[identity, quarter_turn, half_turn, inverse_quarter]],
        max_lift_assignments=16,
    )

    assert len(embeddings) == 16
    assert all(len(embedding.operations) == 4 for embedding in embeddings)


def test_method3_finite_quotient_restricts_to_lattice_normalizer() -> None:
    identity = identity_matrix()
    quarter_turn = rational_matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    half_turn = rational_matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    inverse_quarter = rational_matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=tuple(
            AffineOperation(rotation, zero)
            for rotation in (identity, quarter_turn, half_turn, inverse_quarter)
        ),
        space_group_number=75,
        hall_number=349,
        symbol="P4",
    )
    # The 2a x b x c sublattice is not invariant under a 90-degree rotation.
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])

    with pytest.raises(ValueError, match="outside the parent group"):
        enumerate_lifted_point_subgroups(
            parent,
            child_lattice,
            [[identity, quarter_turn, half_turn, inverse_quarter]],
        )


def test_method3_affine_quotient_limit_precedes_translation_coset_allocation(
    monkeypatch,
) -> None:
    identity = identity_matrix()
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(AffineOperation(identity, zero),),
        space_group_number=1,
        hall_number=1,
        symbol="P1",
    )
    child_lattice = rational_matrix(
        [[1_000_000, 0, 0], [0, 1, 0], [0, 0, 1]]
    )

    def unexpected_coset_allocation(*_args, **_kwargs):
        raise AssertionError("translation cosets must not be allocated after preflight")

    monkeypatch.setattr(
        affine_embeddings_module,
        "translation_cosets",
        unexpected_coset_allocation,
    )
    with pytest.raises(
        RuntimeError,
        match=r"quotient order 1000000.*no translation cosets",
    ):
        enumerate_lifted_point_subgroups(
            parent,
            child_lattice,
            [[identity]],
            max_quotient_order=512,
        )


def test_method3_affine_lift_limit_precedes_translation_coset_allocation(
    monkeypatch,
) -> None:
    identity = identity_matrix()
    twofold = rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(
            AffineOperation(identity, zero),
            AffineOperation(twofold, zero),
        ),
        space_group_number=3,
        hall_number=3,
        symbol="P2",
    )
    child_lattice = rational_matrix(
        [[1_000_000, 0, 0], [0, 1, 0], [0, 0, 1]]
    )

    def unexpected_coset_allocation(*_args, **_kwargs):
        raise AssertionError("translation cosets must not be allocated after preflight")

    monkeypatch.setattr(
        affine_embeddings_module,
        "translation_cosets",
        unexpected_coset_allocation,
    )
    with pytest.raises(
        RuntimeError,
        match=r"1000000 generator-lift assignments.*no translation cosets",
    ):
        enumerate_lifted_point_subgroups(
            parent,
            child_lattice,
            [[identity, twofold]],
            max_quotient_order=2_000_000,
            max_lift_assignments=10,
        )


def test_method3_target_type_enumeration_identifies_synthetic_p2_lifts() -> None:
    identity = identity_matrix()
    twofold = rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(
            AffineOperation(identity, zero),
            AffineOperation(twofold, zero),
        ),
        space_group_number=3,
        hall_number=3,
        symbol="P2",
    )
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])

    embeddings = enumerate_target_affine_embeddings(
        parent,
        child_lattice,
        3,
        [[3.0, 0, 0], [0, 4.0, 0], [0, 0, 5.0]],
    )

    assert len(embeddings) == 2


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

    items = IsoSearchEngine(_NoSingleIrRoutes()).method_3_search(
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

    default_items = IsoSearchEngine(_NoSingleIrRoutes()).method_3_search(
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


def test_sublattice_check_general_matrix():
    """官网 Conventional/Primitive lattice 过滤的子格判定（一般 3x3 矩阵）。"""
    ident = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    doubled = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    fmmm = [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]
    assert _basis_is_sublattice_of(doubled, doubled) is True
    assert _basis_is_sublattice_of(ident, ident) is True
    # 倍格是单位格的子格；单位格不是倍格的子格
    assert _basis_is_sublattice_of(doubled, ident) is True
    assert _basis_is_sublattice_of(ident, doubled) is False
    # Fmmm 型非对角超胞（L_fmmm = {x+y 为偶数}）：倍格(全偶坐标)是其子格
    assert _basis_is_sublattice_of(fmmm, fmmm) is True
    assert _basis_is_sublattice_of(doubled, fmmm) is True
    assert _basis_is_sublattice_of(fmmm, doubled) is False


def test_method1_selector_matches_one_lattice_class():
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    doubled = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    assert _basis_is_same_lattice(identity, identity)
    assert not _basis_is_same_lattice(doubled, identity)


def test_sublattice_check_parent_point_group_orbit():
    """官网 lattice 选项含母相点群旋转：B 相对 S@R 亦算命中。"""
    s = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    # 右乘点群旋转：把 (2,1,1) 型子格转到 (1,2,1) 型（一般非 GL 等价）
    r = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    b = (np.asarray(s, dtype=float) @ r.T).tolist()
    assert _basis_is_sublattice_of(b, s) is False
    assert _basis_is_sublattice_of(b, s, parent_rotations=[r]) is True


@dataclass
class _DummyIsoWrapper:
    """最小桩对象：提供 search engine 所需接口。"""

    def list_k_points(self, parent_sg):
        return []

    def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
        return [
            SubgroupInfo(index=0, space_group_number=225, space_group_symbol="Fm-3m",
                         subgroup_index=1, size=1, is_maximal=True,
                         opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
                         basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         k_point_label="GM", irrep_label="GM1+"),
            SubgroupInfo(index=1, space_group_number=139, space_group_symbol="I4/mmm",
                         subgroup_index=4, size=1, is_maximal=True,
                         opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
                         basis_vectors=[[0, 1, 0], [0, 0, 1], [1, 0, 0]],
                         k_point_label="X", irrep_label="X1"),
            SubgroupInfo(index=2, space_group_number=62, space_group_symbol="Pnma",
                         subgroup_index=8, size=1, is_maximal=False,
                         opd_symbol="C1", opd_vector=[1.0, 1.0, 0.0],
                         basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         k_point_label="L", irrep_label="L2"),
        ]

    def calc_distortion_modes(self, parent_sg, subgroup, wyckoff_letters):
        return [
            DistortionMode(
                irrep_label=subgroup.irrep_label,
                dimension=1,
                opd_symbol=subgroup.opd_symbol,
                bush_modes=[
                    BushMode(irrep_label=subgroup.irrep_label,
                             opd_symbol=subgroup.opd_symbol,
                             wyckoff_letter="a",
                             point=[0.0, 0.0, 0.0],
                             displacements=[[1.0, 0.0, 0.0]]),
                ],
            ),
        ]


def _make_engine() -> IsoSearchEngine:
    return IsoSearchEngine(_DummyIsoWrapper())


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
        for target, centering in ((139, "d"), (123, "P")):
            items = IsoSearchEngine(_NoRoutes()).method_3_search(
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
            assert items[0].route_resolution == "affine_only_unresolved_coupled_route"
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
                SubgroupInfo(
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
        items = IsoSearchEngine(backend).method_3_search(
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


def test_method_1_filters():
    engine = _make_engine()

    all_items = engine.method_1_search(225, Method1Query())
    assert len(all_items) == 3

    cubic_only = engine.method_1_search(225, Method1Query(crystal_system="cubic"))
    assert len(cubic_only) == 1
    assert cubic_only[0].subgroup.space_group_number == 225

    sg_filtered = engine.method_1_search(
        225, Method1Query(subgroup_space_group=62)
    )
    assert len(sg_filtered) == 1
    assert sg_filtered[0].subgroup.space_group_number == 62

    maximal_only = engine.method_1_search(
        225, Method1Query(maximal_subgroup_only=True)
    )
    assert all(item.is_maximal for item in maximal_only)
    assert len(maximal_only) == 2


def test_method_1_lattice_filter():
    engine = _make_engine()
    # 请求 (2,2,2) 子格（对角阵）：仅基矢可被 2 整除的子群保留
    items = engine.method_1_search(
        225, Method1Query(lattice=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    )
    # 桩数据中 index=2 (Pnma, 单位基矢) 被过滤
    assert len(items) == 0 or all(
        item.subgroup.index != 2 for item in items
    )


def test_method_2_metadata_and_modes():
    engine = _make_engine()
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)

    result = engine.method_2_search(
        225,
        subs,
        Method2Query(subgroup_idx=1, distortion_type="displacement"),
        wyckoff_letters=["a", "b"],
    )

    assert result.subgroup.index == 1
    assert result.subgroup.space_group_number == 139
    assert len(result.modes) == 1
    assert result.modes[0].bush_modes[0].displacements == [[1.0, 0.0, 0.0]]
    # k / IR / OPD 元数据来自所选子群（k 点坐标等不再重复传参）
    assert result.metadata["k_point_label"] == "X"
    assert result.metadata["irrep_label"] == "X1"
    assert result.metadata["number_of_independent_modulations"] == 0


def test_iso_bush_queries_every_wyckoff_orbit_before_merging(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = object.__new__(IsoWrapper)
    seen: list[str] = []

    def fake_session(commands: str, timeout=None) -> str:
        del timeout
        seen.append(commands)
        return """
Irrep (ML) Dir(dom) Wyckoff Point            Displacement
GM5+       C1(1)    d       (0,1/2,1/4)      (1,0,0), (0,1,0)
Irrep (ML) Dir(dom) Wyckoff Point            Displacement
GM5+       C1(1)    e       (0,0,z)          (1,0,0)
"""

    monkeypatch.setattr(wrapper, "_run_session", fake_session)
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=2,
        k_point_label="GM",
        irrep_label="GM5+",
        opd_symbol="C1",
    )

    modes = wrapper.calc_distortion_modes(139, subgroup, ["d", "e", "d"])

    assert len(seen) == 1
    assert seen[0].count("VALUE WYCKOFF d\nDISPLAY BUSH") == 1
    assert seen[0].count("VALUE WYCKOFF e\nDISPLAY BUSH") == 1
    assert len(modes) == 3
    assert [mode.wyckoff_site for mode in modes] == ["d", "d", "e"]
    assert [mode.opd_component for mode in modes] == ["a", "b", "a"]
    assert [mode.dimension for mode in modes] == [2, 2, 1]
    assert [mode.amplitude_key for mode in modes] == [
        "GM5+__d__C1(1)__a",
        "GM5+__d__C1(1)__b",
        "GM5+__e__C1(1)__a",
    ]
    assert all(
        len(row.displacements) == 1
        for mode in modes
        for row in mode.bush_modes
    )
    assert modes[0].bush_modes[0].displacements == [[1.0, 0.0, 0.0]]
    assert modes[1].bush_modes[0].displacements == [[0.0, 1.0, 0.0]]


def test_method_2_nmod_accepted():
    engine = _make_engine()
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)
    result = engine.method_2_search(
        225,
        subs,
        Method2Query(subgroup_idx=1, number_of_independent_modulations=1),
        wyckoff_letters=["a"],
    )
    assert result.metadata["number_of_independent_modulations"] == 1
    assert result.modes


def test_method_2_negative_nmod_rejected():
    engine = _make_engine()
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)
    try:
        engine.method_2_search(
            225,
            subs,
            Method2Query(subgroup_idx=1, number_of_independent_modulations=-1),
            wyckoff_letters=["a"],
        )
        raise AssertionError("expected ValueError")
    except ValueError as exc:
        assert "number_of_independent_modulations" in str(exc)


def test_method_2_missing_subgroup_rejected():
    engine = _make_engine()
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)
    try:
        engine.method_2_search(
            225, subs, Method2Query(subgroup_idx=99),
            wyckoff_letters=["a"],
        )
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


def test_method_3_space_group_precedence():
    engine = _make_engine()

    items = engine.method_3_search(
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
    engine = IsoSearchEngine(backend)
    results = []
    for sg in (225, 62):
        results.append(engine.method_3_search(
            225,
            Method3Query(
                distortion_types=["strain", "displacive"],
                space_group_type=sg,
                supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            ),
        ))
    assert backend.calls == 1
    repeated = engine.method_3_search(
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

    items = IsoSearchEngine(_IdentityIso()).method_3_search(
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
    M = _integer_basis_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 6]])
    assert M is not None
    cands = _param_value_candidates(M)
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
    assert _k_compatible_with_supercell(k, M)
    assert not _k_compatible_with_supercell(
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
    assert not _k_compatible_with_supercell(
        ("0.5000000001", "0", "0"),
        [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    assert _k_compatible_with_supercell(
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
    subgroup = SubgroupInfo(
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
                SubgroupInfo(
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
    items = IsoSearchEngine(backend).method_3_search(
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
                SubgroupInfo(
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
    items = IsoSearchEngine(backend).method_3_search(
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
    items = IsoSearchEngine(_SpecialEndpointIso()).method_3_search(
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
        IsoSearchEngine(_BrokenParametricIso()).method_3_search(
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
        IsoSearchEngine(_BrokenSubgroupIso()).method_3_search(
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
    return IsoSearchEngine(backend).method_3_search(
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
        IsoSearchEngine(backend).method_3_search(
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
                SubgroupInfo(
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

    items = IsoSearchEngine(_ParametricIso()).method_3_search(
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


def _method3_route(
    irrep: str,
    *,
    index: int = 0,
    origin: list[float] | None = None,
    k_parameters: list[str] | None = None,
) -> SubgroupInfo:
    return SubgroupInfo(
        index=index,
        space_group_number=99,
        space_group_symbol="P4mm",
        subgroup_index=24,
        size=12,
        opd_symbol="C1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
        origin=origin or [0, 0, 0],
        k_point_label="LD",
        irrep_label=irrep,
        k_parameters=k_parameters or ["1/6"],
        parent_sg=139,
    )


def _method3_filter_api(smodes) -> IsoDistort:
    api = object.__new__(IsoDistort)
    api.structure = [SimpleNamespace(species_string="Eu")]
    api.symmetry_info = {"space_group_number": 139, "wyckoff_sites": []}
    api.distortion_types = ["strain", "displacive"]
    api.distortion_scope = {
        "displacive": ["Eu"],
        "rotational": [],
        "occupational": ["Eu"],
        "magnetic": [],
        "strain": [],
    }
    api._smodes = smodes
    return api


def test_method_3_filters_routes_with_parameter_specific_activity() -> None:
    class _RecordingSmodes:
        def __init__(self) -> None:
            self.calls = []

        def active_irreps(
            self, structure, parent_sg, wyckoff_sites, k_label,
            k_parameters=None, species_filter=None,
        ):
            self.calls.append((k_label, k_parameters, species_filter))
            return {"LD1"}

    smodes = _RecordingSmodes()
    api = _method3_filter_api(smodes)
    # Route identity stores the website's physical g=1/6, while the bundled
    # backend represents this LD parameter on its internal half scale.
    api._resolve_iso_kparams = lambda label, params: ["1/12"]
    inactive = _method3_route("LD2", index=7)
    active = _method3_route("LD1", index=8)
    missing = _method3_route("LD3", index=9, origin=[0.5, 0, 0])
    items = [
        Method3ResultItem(
            subgroup=inactive,
            point_group="4mm",
            basis=inactive.basis_vectors,
            routes=[inactive, active],
        ),
        Method3ResultItem(
            subgroup=missing,
            point_group="4mm",
            basis=missing.basis_vectors,
            routes=[missing],
        ),
        Method3ResultItem(
            subgroup=_method3_route("LD1", index=10, origin=[0, 0.5, 0]),
            point_group="4mm",
            basis=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
            routes=[],
        ),
    ]

    filtered = api._filter_method3_routes_by_types(items, ["displacive"])

    assert len(filtered) == 1
    assert filtered[0].routes == [active]
    assert filtered[0].subgroup is active
    assert filtered[0].subgroup.index == 0
    assert smodes.calls == [("LD", ["1/12"], {"Eu"})]


def test_method_2_rejects_affine_only_method3_diagnostic_row() -> None:
    unresolved = _method3_route("", index=4)
    api = _method3_filter_api(SimpleNamespace())
    api.subgroups = [unresolved]
    api._unresolved_method3_embedding_keys = {
        api._method3_embedding_guard_key(unresolved)
    }

    # The guard is a scientific path identity, not a Python object address:
    # an equivalent result reconstructed by another API/UI layer remains
    # non-selectable even without the transient marker attribute.
    equivalent = deepcopy(unresolved)
    assert api._method3_embedding_guard_key(equivalent) == (
        api._method3_embedding_guard_key(unresolved)
    )

    close_but_distinct = deepcopy(unresolved)
    close_but_distinct.origin = [1e-11, 0.0, 0.0]
    assert api._method3_embedding_guard_key(close_but_distinct) != (
        api._method3_embedding_guard_key(unresolved)
    )

    with pytest.raises(RuntimeError, match="no resolved single-IR or coupled-IR"):
        api.search_method_2(
            subgroup_idx=4,
            candidates=[equivalent],
        )


def test_method_3_route_filter_uses_existing_magnetic_and_occupational_semantics(
    monkeypatch,
) -> None:
    class _NoSmodes:
        def active_irreps(self, *args, **kwargs):
            raise AssertionError("smodes must not run for magnetic/occupational-only filters")

    class _OccupationalGenerator:
        def generate(self, structure, wyckoff_sites, subgroup, species_scope):
            assert species_scope == {"Eu"}
            return (
                [SimpleNamespace(validated=True)]
                if subgroup.irrep_label == "Xocc" else []
            )

    monkeypatch.setattr(
        core_api_module, "OccupationalModeGenerator", _OccupationalGenerator
    )
    api = _method3_filter_api(_NoSmodes())
    magnetic = _method3_route("mLD1")
    occupational = _method3_route("Xocc", origin=[0.5, 0, 0])
    ordinary = _method3_route("LD1", origin=[0, 0.5, 0])
    items = [
        Method3ResultItem(magnetic, "4mm", magnetic.basis_vectors, [magnetic]),
        Method3ResultItem(
            occupational, "4mm", occupational.basis_vectors, [occupational]
        ),
        Method3ResultItem(ordinary, "4mm", ordinary.basis_vectors, [ordinary]),
    ]

    magnetic_only = api._filter_method3_routes_by_types(items, ["magnetic"])
    assert [item.subgroup.irrep_label for item in magnetic_only] == ["mLD1"]

    items = [
        Method3ResultItem(magnetic, "4mm", magnetic.basis_vectors, [magnetic]),
        Method3ResultItem(
            occupational, "4mm", occupational.basis_vectors, [occupational]
        ),
        Method3ResultItem(ordinary, "4mm", ordinary.basis_vectors, [ordinary]),
    ]
    occupational_only = api._filter_method3_routes_by_types(items, ["occupational"])
    assert [item.subgroup.irrep_label for item in occupational_only] == ["Xocc"]


def test_method_3_route_filter_rejects_unvalidated_occupational_guess(
    monkeypatch,
) -> None:
    class _NoSmodes:
        def active_irreps(self, *args, **kwargs):
            raise AssertionError("smodes must not run for occupational-only filters")

    class _UnvalidatedOccupationalGenerator:
        def generate(self, structure, wyckoff_sites, subgroup, species_scope):
            return [SimpleNamespace(validated=False)]

    monkeypatch.setattr(
        core_api_module,
        "OccupationalModeGenerator",
        _UnvalidatedOccupationalGenerator,
    )
    api = _method3_filter_api(_NoSmodes())
    route = _method3_route("Xocc")
    item = Method3ResultItem(route, "4mm", route.basis_vectors, [route])

    assert api._filter_method3_routes_by_types([item], ["occupational"]) == []


def test_method_1_strain_filter_uses_symmetric_tensor_fixed_space() -> None:
    api = IsoDistort()
    api.structure = _i4mmm_parent()
    api.symmetry_info = {
        "space_group_number": 139,
        "wyckoff_sites": [],
    }

    def _item(
        irrep: str,
        child_sg: int,
        child_symbol: str,
        basis: list[list[int]],
    ) -> Method1ResultItem:
        route = SubgroupInfo(
            index=0,
            parent_sg=139,
            k_point_label="GM",
            irrep_label=irrep,
            size=1,
            subgroup_index=1,
            space_group_number=child_sg,
            space_group_symbol=child_symbol,
            basis_vectors=basis,
            origin=[0, 0, 0],
        )
        return Method1ResultItem(
            subgroup=route,
            crystal_system="tetragonal",
            is_maximal=False,
        )

    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    gm1 = _item("GM1+", 139, "I4/mmm", identity)
    gm4 = _item(
        "GM4+",
        69,
        "Fmmm",
        [[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
    )
    gm3 = _item("GM3+", 87, "I4/m", identity)
    polar = _item("GM5-", 25, "Pmm2", identity)

    kept = api._filter_method1_by_types(
        [gm1, gm4, gm3, polar], ["strain"]
    )

    assert [item.subgroup.irrep_label for item in kept] == ["GM1+", "GM4+"]


def test_special_modes_complete_only_wholly_rootless_orbits(monkeypatch) -> None:
    subgroup = _method3_route("P2")
    subgroup.k_parameters = []
    bush = DistortionMode(
        irrep_label="P2",
        wyckoff_site="d",
        amplitude_key="P2__d__C1__a",
        bush_modes=[
            BushMode(
                irrep_label="P2",
                opd_symbol="C1",
                wyckoff_letter="d",
                point=[0, 0, 0],
                displacements=[[1, 0, 0]],
            )
        ],
    )
    secondary = DistortionMode(
        irrep_label="GM1+",
        wyckoff_site="e",
        amplitude_key="GM1+[0,0,0]__e__A1(a)",
    )

    class _Iso:
        def calc_distortion_modes(self, parent_sg, selected, wyckoff_letters):
            assert wyckoff_letters == ["d", "e"]
            return [bush]

    def _fixed_space(parent, symmetry_info, selected, letters, smodes, **kwargs):
        assert letters == ["e"]
        return superspace_module.ParametricModeResult(
            modes=[secondary],
            supercell_displacements={secondary.amplitude_key: np.zeros((2, 3))},
            labels={secondary.amplitude_key: "secondary"},
            nmod=0,
        )

    monkeypatch.setattr(superspace_module, "compute_parametric_modes", _fixed_space)
    result = superspace_module.compute_special_modes_with_rootless_supplement(
        object(),
        {"space_group_number": 139},
        subgroup,
        ["d", "e"],
        _Iso(),
        object(),
    )

    assert result.modes == [bush, secondary]
    assert set(result.supercell_displacements) == {secondary.amplitude_key}
    assert "rootless orbit(s) e" in result.note


def test_method_4_decomposition_recovery():
    engine = _make_engine()

    parent = Structure(Lattice.cubic(5.0), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    distorted = parent.copy()
    distorted.translate_sites([0], [0.03, 0.0, 0.0], frac_coords=True)

    mode_displacements = {
        "GM1+": np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        "R5-": np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
    }

    result = engine.method_4_decompose(
        parent,
        distorted,
        mode_displacements,
        Method4Query(atom_matching_method="nearest-site"),
    )

    assert abs(result.amplitudes["GM1+"] - 0.03) < 1e-8
    assert abs(result.amplitudes["R5-"]) < 1e-8
    assert result.rms_residual < 1e-10
