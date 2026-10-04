"""Fail-closed provenance tests for ISO microscopic displacement columns."""

from __future__ import annotations

import hashlib
from dataclasses import replace
from fractions import Fraction
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from isocore.backend.iso_mode_models import (
    InvariantDirection,
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
    microscopic_atom_order_id,
)
from isocore.backend.iso_wrapper import (
    BushMode,
    DistortionMode,
    IsoWrapper,
    KPointInfo,
    SubgroupInfo,
    microscopic_extension_digest,
)
from isocore.distortion.distortion_mapper import (
    DistortionMapper,
    MappedMicroscopicColumn,
    verify_mapped_microscopic_columns,
)
from isocore.distortion.superspace import (
    _extend_microscopic_modes_to_reference_domain,
    _select_verified_microscopic_modes,
    _symbolic_subspace_proof,
    _unresolved_modes,
)
from isocore.io.cif_displacive_model import (
    ChildAtom,
    DisplaciveModelError,
    ExactSeitzOperation,
    ModeScaleProvenance,
    RawModeColumn,
    project_modes_to_free_rows,
)
from isocore.utils.exceptions import OutputParseError

QUERY_DIGEST = hashlib.sha256(b"synthetic microscopic query").hexdigest()


def _block(
    *,
    irrep: str = "X1+",
    source_order: int = 0,
) -> MicroscopicVectorBlock:
    return MicroscopicVectorBlock(
        global_irrep=irrep,
        wyckoff_letter="a",
        site_irrep="A",
        rows=(
            MicroscopicVectorRow(
                point_raw=("0", "0", "0"),
                displacements=(
                    (Fraction(1), Fraction(0), Fraction(0)),
                    (Fraction(0), Fraction(1), Fraction(0)),
                ),
            ),
            MicroscopicVectorRow(
                point_raw=("1", "0", "0"),
                displacements=(
                    (Fraction(-1), Fraction(0), Fraction(0)),
                    (Fraction(0), Fraction(-1), Fraction(0)),
                ),
            ),
        ),
        source_order=source_order,
        direction_symbol="VECTOR,A,B",
    )


def _provenance(column: int, *, block: MicroscopicVectorBlock | None = None):
    return MicroscopicColumnProvenance.from_exact_column(
        query_digest=QUERY_DIGEST,
        query_order=0,
        query_irrep_label=(block or _block()).global_irrep,
        source_parent_sg=62,
        source_k_coordinates=("0", "1/2", "0"),
        source_subgroup_space_group_number=1,
        source_subgroup_basis=np.eye(3, dtype=int).tolist(),
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block or _block(),
        column_index=column,
    )


def _subgroup(*, irrep: str = "GM1+") -> SubgroupInfo:
    return SubgroupInfo(
        index=0,
        space_group_number=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label=irrep,
        parent_sg=62,
        k_coordinates=["0", "0", "0"],
        basis_vectors=np.eye(3).tolist(),
        origin=[0.0, 0.0, 0.0],
    )


def _mode_from_block(block: MicroscopicVectorBlock, column: int) -> DistortionMode:
    component = chr(ord("a") + column)
    provenance = _provenance(column, block=block)
    bushes = [
        BushMode(
            irrep_label=block.global_irrep,
            opd_symbol="P1",
            wyckoff_letter=block.wyckoff_letter,
            point=[float(Fraction(value)) for value in row.point_raw],
            point_raw=list(row.point_raw),
            displacements=[
                [float(value) for value in row.displacements[column]]
            ],
        )
        for row in block.rows
    ]
    return DistortionMode(
        irrep_label=block.global_irrep,
        dimension=block.column_count,
        wyckoff_site=block.wyckoff_letter,
        wyckoff_orbit_id="orbit-a",
        opd_symbol="P1",
        opd_dir_raw=(
            "(" + ",".join(
                chr(ord("a") + index) for index in range(block.column_count)
            ) + ")"
        ),
        bush_modes=bushes,
        amplitude_key=f"{block.global_irrep}__{component}",
        site_irrep=block.site_irrep,
        k_coords_label="0,1/2,0",
        opd_component=component,
        mode_identity=ModeIdentity(
            parent_sg=62,
            global_irrep=block.global_irrep,
            k_coordinates=("0", "1/2", "0"),
            wyckoff_letter=block.wyckoff_letter,
            orbit_id="orbit-a",
            site_irrep=block.site_irrep,
            component_index=column,
            component_label=component,
            source="iso_microscopic",
            status="partial",
            reason="awaiting_bush_or_smodes_subspace_validation",
        ),
        microscopic_provenance=provenance,
    )


def _current_mode(vector: tuple[float, float, float], index: int) -> DistortionMode:
    return DistortionMode(
        irrep_label="X1+",
        dimension=2,
        wyckoff_site="a",
        bush_modes=[
            BushMode(
                irrep_label="X1+",
                opd_symbol="P1",
                wyckoff_letter="a",
                point=[0.0, 0.0, 0.0],
                point_raw=["0", "0", "0"],
                displacements=[list(vector)],
            ),
            BushMode(
                irrep_label="X1+",
                opd_symbol="P1",
                wyckoff_letter="a",
                point=[1.0, 0.0, 0.0],
                point_raw=["1", "0", "0"],
                displacements=[[-value for value in vector]],
            ),
        ],
        amplitude_key=f"current-{index}",
    )


def _mapped(
    component_label: str,
    provenance: MicroscopicColumnProvenance,
    displacements: np.ndarray,
    *,
    canonical: bool,
    frame_id: str = "child-frame",
    atom_ids: tuple[str, ...] = ("a", "b"),
    amplitude_key: str | None = None,
    mode_identity: ModeIdentity | None = None,
) -> MappedMicroscopicColumn:
    bound_provenance = provenance.bind_mapping(
        frame_id=frame_id,
        atom_ids=atom_ids,
    )
    identity = mode_identity or ModeIdentity(
        parent_sg=62,
        global_irrep=provenance.source_global_irrep,
        k_coordinates=("0", "1/2", "0"),
        wyckoff_letter=provenance.source_wyckoff_letter,
        orbit_id="orbit-a",
        site_irrep=provenance.source_site_irrep,
        component_index=provenance.source_column_index,
        component_label=component_label,
        source="iso_microscopic",
        status="verified",
    )
    return MappedMicroscopicColumn(
        mode_id=identity.stable_token,
        amplitude_key=amplitude_key or f"amp-{component_label}",
        frame_id=frame_id,
        atom_order_id=microscopic_atom_order_id(atom_ids),
        atom_ids=atom_ids,
        displacements=displacements,
        mode_identity=identity,
        provenance=bound_provenance,
        relation_status="canonical" if canonical else "unverified",
        current_over_canonical_signed_scale=1.0 if canonical else None,
        cartesian_residual_angstrom=0.0 if canonical else None,
    )


def test_direction_resolved_microscopic_query_binds_exact_route_and_opd():
    subgroup = _subgroup(irrep="X1+")
    irreps, letters, commands = IsoWrapper._microscopic_query_commands(
        62,
        subgroup,
        ["X1+", "GM1+"],
        ["a"],
        direction_symbols={
            "X1+": "VECTOR,0,A",
            "GM1+": "VECTOR,A",
        },
    )

    assert irreps == ("X1+", "GM1+")
    assert letters == ("a",)
    text = "\n".join(commands)
    assert "VALUE KPOINT X\nVALUE IRREP X1+\nVALUE DIRECTION VECTOR,0,A" in text
    assert "VALUE KPOINT GM\nVALUE IRREP GM1+\nVALUE DIRECTION VECTOR,A" in text
    assert text.count("VALUE SUBGROUP 1") == 2
    assert text.count("DISPLAY DISTORTION") == 2


def test_microscopic_vector_query_rejects_syntax_error_with_parseable_full_block():
    wrapper = object.__new__(IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: """
        ******Syntax error: VECTOR,A/2,0
        Irrep (ML) Wyckoff Irrep Point Projected Vectors
        X1+ a A (0,0,0) (1,0,0)
        *
    """

    with pytest.raises(OutputParseError, match="VECTOR"):
        wrapper.get_microscopic_vector_blocks(
            62,
            _subgroup(irrep="X1+"),
            ["X1+"],
            ["a"],
            direction_symbols={"X1+": "VECTOR,A"},
        )


def test_microscopic_vector_query_rejects_unconditioned_column_width():
    wrapper = object.__new__(IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: """
        Irrep (ML) Wyckoff Irrep Point Projected Vectors
        X1+ a A (0,0,0) (1,0,0), (0,1,0)
        *
    """

    with pytest.raises(OutputParseError, match="错误列数"):
        wrapper.get_microscopic_vector_blocks(
            62,
            _subgroup(irrep="X1+"),
            ["X1+"],
            ["a"],
            direction_symbols={"X1+": "VECTOR,A"},
        )


@pytest.mark.parametrize(
    "unexpected_row",
    (
        "Y1+ a A (0,0,0) (1,0,0)",
        "X1+ b A (0,0,0) (1,0,0)",
    ),
)
def test_microscopic_vector_query_rejects_unrequested_output_blocks(
    unexpected_row,
):
    wrapper = object.__new__(IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: f"""
        Irrep (ML) Wyckoff Irrep Point Projected Vectors
        {unexpected_row}
        *
    """

    with pytest.raises(OutputParseError, match="未请求"):
        wrapper.get_microscopic_vector_blocks(
            62,
            _subgroup(irrep="X1+"),
            ["X1+"],
            ["a"],
            direction_symbols={"X1+": "VECTOR,A"},
        )


def test_exact_invariant_vector_resolves_one_cached_direction_symbol():
    wrapper = object.__new__(IsoWrapper)
    calls = []

    def list_subgroups(parent_sg, k_point, irrep_label):
        calls.append((parent_sg, k_point, irrep_label))
        return [
            SimpleNamespace(
                opd_dir_raw="(a, 0)", opd_symbol="P1",
                space_group_number=2, size=2,
            ),
            SimpleNamespace(
                opd_dir_raw="(a; a)", opd_symbol="P3",
                space_group_number=1, size=1,
            ),
        ]

    wrapper.list_subgroups = list_subgroups
    invariant = InvariantDirection("X1+", "(a,a)", 1, "P1", 1)

    assert wrapper._direction_symbol_for_invariant(62, invariant) == "P3"
    assert wrapper._direction_symbol_for_invariant(62, invariant) == "P3"
    assert calls == [(62, "X", "X1+")]

    wrapper._direction_symbol_cache = {}
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(a,a)", opd_symbol="P3",
            space_group_number=1, size=1,
        ),
        SimpleNamespace(
            opd_dir_raw="(b,b)", opd_symbol="C1",
            space_group_number=1, size=1,
        ),
    ]
    with pytest.raises(ValueError, match="one ISO OPD token"):
        wrapper._direction_symbol_for_invariant(62, invariant)


def test_exact_direction_vector_serializer_preserves_component_order():
    assert IsoWrapper._direction_vector_selector("(0,a)") == "VECTOR,0,A"
    assert (
        IsoWrapper._direction_vector_selector("(a;b;b;-a)")
        == "VECTOR,A,B,B,-A"
    )
    assert (
        IsoWrapper._direction_vector_selector("(2*a+b;4*a-3*b)")
        == "VECTOR,A+B,2A-3B"
    )
    assert IsoWrapper._direction_vector_selector("(3*a,a)") == "VECTOR,3A,A"
    assert IsoWrapper._direction_vector_selector("(a/2,a/3)") == "VECTOR,3A,2A"


@pytest.mark.parametrize(
    "value",
    [
        "(0,0)",
        "(a*a,0)",
        "(sqrt(2)*a,0)",
        "(0.5*a,0)",
        "(a+1,0)",
        "(a+b,2*a+2*b)",
        "(a+A,0)",
        "a,0",
    ],
)
def test_exact_direction_vector_serializer_rejects_unrepresentable_forms(value):
    with pytest.raises(ValueError):
        IsoWrapper._direction_vector_selector(value)


def test_primary_direction_symbol_rejects_permuted_star_arm():
    wrapper = object.__new__(IsoWrapper)
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(a;a)", opd_symbol="P1",
            space_group_number=127, size=4,
        ),
        SimpleNamespace(
            opd_dir_raw="(2*b;0)", opd_symbol="P3",
            space_group_number=64, size=2,
        ),
        SimpleNamespace(
            opd_dir_raw="(a;b)", opd_symbol="C1",
            space_group_number=55, size=4,
        ),
    ]

    invariant = InvariantDirection("X2+", "(0,a)", 64, "Cmce", 2)

    with pytest.raises(ValueError, match=r"tokens=\[\]"):
        wrapper._direction_symbol_for_invariant(139, invariant)


def test_direction_symbol_rejects_wrong_isotropy_identity_and_ambiguity():
    wrapper = object.__new__(IsoWrapper)
    invariant = InvariantDirection("X2+", "(0,a)", 64, "Cmce", 2)
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(a;0)", opd_symbol="P3",
            space_group_number=65, size=2,
        )
    ]
    with pytest.raises(ValueError, match=r"tokens=\[\]"):
        wrapper._direction_symbol_for_invariant(139, invariant)

    wrapper._direction_symbol_cache = {}
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(0;b)", opd_symbol="P3",
            space_group_number=64, size=2,
        ),
        SimpleNamespace(
            opd_dir_raw="(0;c)", opd_symbol="C9",
            space_group_number=64, size=2,
        ),
    ]
    with pytest.raises(ValueError, match="one ISO OPD token"):
        wrapper._direction_symbol_for_invariant(139, invariant)


def test_parser_modes_retain_exact_column_and_unresolved_source_facts():
    block = replace(
        _block(irrep="X1+", source_order=7),
        direction_symbol="VECTOR,A,B",
    )
    wrapper = object.__new__(IsoWrapper)
    wrapper.get_microscopic_vector_blocks = lambda *_args, **_kwargs: [block]
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1),
        InvariantDirection("X1+", "(a,a)", 1, "P1", 1),
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P1"
    wrapper.list_k_points = lambda _parent: [
        KPointInfo("X", ["0", "1/2", "0"], [], True)
    ]

    modes = IsoWrapper.calc_microscopic_distortion_modes(
        wrapper,
        62,
        _subgroup(),
        ["a"],
        irrep_labels=["X1+"],
    )

    assert [mode.microscopic_provenance.source_block_order for mode in modes] == [7, 7]
    assert [mode.microscopic_provenance.source_column_index for mode in modes] == [0, 1]
    assert all(mode.microscopic_provenance.query_order == 0 for mode in modes)
    assert all(mode.microscopic_provenance.source_row_count == 2 for mode in modes)
    assert modes[0].microscopic_provenance.exact_row_digest == modes[1].microscopic_provenance.exact_row_digest
    assert modes[0].microscopic_provenance.exact_vector_digest != modes[1].microscopic_provenance.exact_vector_digest
    assert all(
        mode.microscopic_provenance.unresolved_fields
        == ("source_frame_id", "source_atom_order_id")
        for mode in modes
    )
    assert [mode.mode_identity.k_coordinates for mode in modes] == [
        ("0", "1/2", "0"),
        ("0", "1/2", "0"),
    ]
    assert all(mode.opd_dir_raw == "(a,a)" for mode in modes)
    assert all(
        mode.microscopic_provenance.source_direction_symbol == "VECTOR,A,B"
        for mode in modes
    )


@pytest.mark.parametrize(
    "listed",
    [
        [KPointInfo("X", ["a", "0", "0"], ["a"], False)],
        [
            KPointInfo("X", ["0", "1/2", "0"], [], True),
            KPointInfo("X", ["1/2", "0", "0"], [], True),
        ],
    ],
)
def test_parametric_or_ambiguous_listed_secondary_k_remains_unresolved(listed):
    block = _block(irrep="X1+")
    wrapper = object.__new__(IsoWrapper)
    wrapper.get_microscopic_vector_blocks = lambda *_args, **_kwargs: [block]
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1),
        InvariantDirection("X1+", "(a,a)", 1, "P1", 1),
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P1"
    wrapper.list_k_points = lambda _parent: listed
    modes = IsoWrapper.calc_microscopic_distortion_modes(
        wrapper, 62, _subgroup(), ["a"], irrep_labels=["X1+"]
    )
    assert all(mode.mode_identity.status == "unresolved" for mode in modes)
    assert all(
        mode.mode_identity.reason == "exact_special_k_identity_unresolved"
        for mode in modes
    )
    assert all(mode.mode_identity.k_coordinates == () for mode in modes)


def test_microscopic_calc_requires_exact_direction_for_explicit_irreps():
    wrapper = object.__new__(IsoWrapper)
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1)
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P1"

    with pytest.raises(ValueError, match="requested irrep"):
        IsoWrapper.calc_microscopic_distortion_modes(
            wrapper,
            62,
            _subgroup(),
            ["a"],
            irrep_labels=["X1+"],
        )


def test_microscopic_calc_rejects_primary_opd_direction_mismatch():
    wrapper = object.__new__(IsoWrapper)
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1)
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P2"

    with pytest.raises(ValueError, match="selected subgroup OPD"):
        IsoWrapper.calc_microscopic_distortion_modes(
            wrapper,
            62,
            _subgroup(),
            ["a"],
            irrep_labels=["GM1+"],
        )


def test_superspace_installs_canonical_iso_basis_after_nondiagonal_change():
    block = _block()
    canonical = [_mode_from_block(block, 0), _mode_from_block(block, 1)]
    mixed = [
        _current_mode((1.0, 1.0, 0.0), 0),
        _current_mode((1.0, -1.0, 0.0), 1),
    ]
    parent = Structure(Lattice.cubic(4.0), ["Fe", "Fe"], [[0, 0, 0], [0, 0, 0]])
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, mixed, list(reversed(canonical)), parent_sg=62
    )
    assert diagnostics[("X1+", "a")].matched
    assert [mode.amplitude_key for mode in selected] == [
        canonical[0].amplitude_key,
        canonical[1].amplitude_key,
    ]
    assert [mode.microscopic_provenance for mode in selected] == [
        canonical[0].microscopic_provenance,
        canonical[1].microscopic_provenance,
    ]
    assert all(mode.mode_identity.status == "verified" for mode in selected)

    selected, diagnostics = _select_verified_microscopic_modes(
        parent, [_current_mode((1.0, 0.0, 0.0), 0)], [canonical[0]], parent_sg=62
    )
    assert diagnostics[("X1+", "a")].reason == "microscopic_source_block_column_missing"
    assert selected[0].mode_identity.status == "unresolved"


def _sparse_domain_extension_fixture():
    block = MicroscopicVectorBlock(
        global_irrep="X1+",
        wyckoff_letter="a",
        site_irrep="A",
        rows=(MicroscopicVectorRow(
            point_raw=("0", "0", "0"),
            displacements=(
                (Fraction(1), Fraction(0), Fraction(0)),
                (Fraction(0), Fraction(1), Fraction(0)),
            ),
        ),),
        source_order=0,
        direction_symbol="VECTOR,A,B",
    )
    canonical = [_mode_from_block(block, 0), _mode_from_block(block, 1)]
    reference = [
        _current_mode((1.0, 1.0, 0.0), 0),
        _current_mode((1.0, -1.0, 0.0), 1),
    ]
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, canonical, parent_sg=62,
    )
    return parent, reference, canonical, selected, diagnostics


def test_sparse_canonical_rows_are_uniquely_extended_with_separate_evidence():
    _parent, _reference, canonical, selected, diagnostics = (
        _sparse_domain_extension_fixture()
    )

    assert diagnostics[("X1+", "a")].matched
    assert len(selected) == 2
    assert all(len(mode.bush_modes) == 2 for mode in selected)
    assert all(mode.mode_identity.status == "verified" for mode in selected)
    assert [mode.microscopic_provenance for mode in selected] == [
        mode.microscopic_provenance for mode in canonical
    ]
    extensions = [mode.microscopic_domain_extension for mode in selected]
    assert all(extension is not None for extension in extensions)
    assert [extension.coefficient_column for extension in extensions] == [
        pytest.approx((0.5, 0.5)),
        pytest.approx((0.5, -0.5)),
    ]
    assert all(
        extension.relation_status == "verified_bush_domain_extension"
        and extension.comparison_domain_kind == "exact_common_points"
        and extension.reference_rank == 2
        and extension.canonical_condition_number <= extension.maximum_condition_number
        and extension.matrix_convention
        == "R_common@T=C_common;C_full=R_full@T"
        for extension in extensions
    )
    assert np.allclose(
        [mode.bush_modes[1].displacements[0] for mode in selected],
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
    )
    assert all(len(mode.bush_modes) == 1 for mode in canonical)
    unresolved = _unresolved_modes(selected, parent_sg=62, reason="forced_failure")
    assert all(mode.microscopic_domain_extension is None for mode in unresolved)


def test_sparse_extension_is_invariant_to_bush_basis_column_order():
    parent, reference, canonical, baseline, _diagnostics = (
        _sparse_domain_extension_fixture()
    )

    selected, diagnostics = _select_verified_microscopic_modes(
        parent,
        list(reversed(reference)),
        list(reversed(canonical)),
        parent_sg=62,
    )

    assert diagnostics[("X1+", "a")].matched
    assert [mode.amplitude_key for mode in selected] == [
        mode.amplitude_key for mode in canonical
    ]
    assert [mode.microscopic_provenance for mode in selected] == [
        mode.microscopic_provenance for mode in canonical
    ]
    for actual, expected in zip(selected, baseline, strict=True):
        assert np.allclose(
            [bush.displacements[0] for bush in actual.bush_modes],
            [bush.displacements[0] for bush in expected.bush_modes],
        )
        assert actual.mode_identity == expected.mode_identity
        assert actual.microscopic_domain_extension is not None
        assert (
            actual.microscopic_domain_extension.canonical_source_order_key
            == actual.microscopic_provenance.source_order_key
        )


def test_domain_extension_proof_rejects_rebound_columns_and_basis():
    parent, reference, canonical, _selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    lattice = np.asarray(parent.lattice.matrix, dtype=float)
    proof = _symbolic_subspace_proof(reference, canonical, lattice)
    assert proof.validation.matched

    with pytest.raises(ValueError, match="proof_candidate_binding_mismatch"):
        _extend_microscopic_modes_to_reference_domain(
            reference, list(reversed(canonical)), lattice, proof,
        )

    rebound = replace(canonical[0], bush_modes=list(canonical[1].bush_modes))
    with pytest.raises(
        ValueError, match="proof_candidate_source_vector_digest_mismatch",
    ):
        _extend_microscopic_modes_to_reference_domain(
            reference, [rebound, canonical[1]], lattice, proof,
        )

    with pytest.raises(ValueError, match="proof_reference_basis_mismatch"):
        _extend_microscopic_modes_to_reference_domain(
            list(reversed(reference)), canonical, lattice, proof,
        )


def test_domain_extension_evidence_survives_mapping_and_raw_contract():
    parent, _reference, _canonical, selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    mapper = DistortionMapper()
    mapped_arrays = {
        mode.amplitude_key: np.asarray(
            [[float(index + 1), 0.0, 0.0], [0.0, float(index + 1), 0.0]]
        )
        for index, mode in enumerate(selected)
    }
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: mapped_arrays
    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        [],
        selected,
        np.eye(3).tolist(),
        subgroup_context=_subgroup(),
        frame_id="child-frame",
        atom_ids=("atom-0", "atom-1"),
    )

    assert [column.domain_extension_evidence for column in mapped] == [
        mode.microscopic_domain_extension for mode in selected
    ]
    raw = RawModeColumn.from_values(
        mapped[0].mode_id,
        0,
        mapped[0].frame_id,
        mapped[0].atom_ids,
        mapped[0].displacements,
        ModeScaleProvenance(
            source="iso_microscopic_with_verified_bush_domain_extension",
            convention="mapped_unmixed_verified_full_domain_extension",
            direction_resolved=True,
        ),
        mode_identity=mapped[0].mode_identity,
        microscopic_provenance=mapped[0].provenance,
        microscopic_domain_extension=mapped[0].domain_extension_evidence,
    )
    assert raw.source_binding_issue() is None
    assert raw.microscopic_domain_extension == mapped[0].domain_extension_evidence


def test_mapper_rejects_tampered_or_incomplete_domain_extension_evidence():
    parent, _reference, _canonical, selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    mapper = DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.ones((2, 3), dtype=float) for mode in selected
    }

    def _map(modes, subgroup=None):
        return mapper.map_microscopic_columns_to_supercell(
            parent,
            [],
            modes,
            np.eye(3).tolist(),
            subgroup_context=subgroup or _subgroup(),
            frame_id="child-frame",
            atom_ids=("atom-0", "atom-1"),
        )

    first_extension = selected[0].microscopic_domain_extension
    second_extension = selected[1].microscopic_domain_extension
    assert first_extension is not None and second_extension is not None
    wrong_coefficient = replace(
        selected[0],
        microscopic_domain_extension=replace(
            first_extension,
            coefficient_column=(0.25, 0.5),
        ),
    )
    with pytest.raises(ValueError, match="change-of-basis digest"):
        _map([wrong_coefficient, selected[1]])

    wrong_output_order = replace(
        selected[1],
        microscopic_domain_extension=replace(
            second_extension,
            canonical_output_column_index=0,
        ),
    )
    with pytest.raises(ValueError, match="canonical column order"):
        _map([selected[0], wrong_output_order])

    swapped_output_order = [
        replace(
            selected[column],
            microscopic_domain_extension=replace(
                selected[column].microscopic_domain_extension,
                canonical_output_column_index=1 - column,
            ),
        )
        for column in range(2)
    ]
    with pytest.raises(ValueError, match="differs from exact source order"):
        _map(swapped_output_order)

    wrong_source_index = replace(
        selected[0],
        microscopic_domain_extension=replace(
            first_extension,
            canonical_source_order_key=(0, 9, 0),
        ),
    )
    with pytest.raises(ValueError, match="another source index"):
        _map([wrong_source_index, selected[1]])

    damaged_bush = list(selected[0].bush_modes)
    damaged_bush[0] = replace(
        damaged_bush[0], displacements=[[9.0, 0.0, 0.0]],
    )
    with pytest.raises(ValueError, match="current full column"):
        _map([replace(selected[0], bush_modes=damaged_bush), selected[1]])

    with pytest.raises(ValueError, match="source basis is incomplete"):
        _map([selected[0]])

    wrong_transport = replace(
        selected[0],
        microscopic_domain_extension=replace(
            first_extension,
            transport_digest="0" * 64,
        ),
    )
    with pytest.raises(ValueError, match="group evidence is inconsistent"):
        _map([wrong_transport, selected[1]])

    with pytest.raises(ValueError, match="numerical proof is invalid"):
        replace(first_extension, maximum_normalized_residual=1.01)

    singular_change = ((0.5, 0.5), (0.5, 0.5))
    singular_digest = microscopic_extension_digest(singular_change)
    singular_modes = [
        replace(
            mode,
            microscopic_domain_extension=replace(
                mode.microscopic_domain_extension,
                coefficient_column=tuple(
                    singular_change[row][column] for row in range(2)
                ),
                change_of_basis_digest=singular_digest,
            ),
        )
        for column, mode in enumerate(selected)
    ]
    with pytest.raises(ValueError, match="change-of-basis rank loss"):
        _map(singular_modes)

    wrong_condition = [
        replace(
            mode,
            microscopic_domain_extension=replace(
                mode.microscopic_domain_extension,
                change_of_basis_condition_number=(
                    mode.microscopic_domain_extension
                    .change_of_basis_condition_number
                    * 1.01
                ),
            ),
        )
        for mode in selected
    ]
    with pytest.raises(ValueError, match="change-of-basis condition mismatch"):
        _map(wrong_condition)

    shifted_context = replace(_subgroup(), origin=[1.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="query context differs"):
        _map(selected, shifted_context)


def test_domain_extension_groups_are_independent_per_physical_orbit():
    parent, _reference, _canonical, first_orbit, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    second_orbit = [
        replace(
            mode,
            amplitude_key=f"{mode.amplitude_key}__orbit-b",
            wyckoff_orbit_id="orbit-b",
            mode_identity=replace(mode.mode_identity, orbit_id="orbit-b"),
        )
        for mode in first_orbit
    ]
    modes = [*first_orbit, *second_orbit]
    mapper = DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.ones((2, 3), dtype=float) for mode in modes
    }

    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        [],
        modes,
        np.eye(3).tolist(),
        subgroup_context=_subgroup(),
        frame_id="child-frame",
        atom_ids=("atom-0", "atom-1"),
    )
    assert len(mapped) == 4
    assert len({column.mode_id for column in mapped}) == 4
    assert len({
        column.provenance.exact_source_token for column in mapped
    }) == 2

    with pytest.raises(ValueError, match="source basis is incomplete"):
        mapper.map_microscopic_columns_to_supercell(
            parent,
            [],
            [*first_orbit, second_orbit[0]],
            np.eye(3).tolist(),
            subgroup_context=_subgroup(),
            frame_id="child-frame",
            atom_ids=("atom-0", "atom-1"),
        )


def test_public_raw_contract_revalidates_complete_extension_groups():
    parent, _reference, _canonical, selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    mapper = DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.asarray(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float,
        )
        for mode in selected
    }
    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        [],
        selected,
        np.eye(3).tolist(),
        subgroup_context=_subgroup(),
        frame_id="child-frame",
        atom_ids=("atom-0", "atom-1"),
    )

    def _raw(column, index, *, extension=None, source=None):
        return RawModeColumn.from_values(
            column.mode_id,
            index,
            column.frame_id,
            column.atom_ids,
            column.displacements,
            ModeScaleProvenance(
                source=source or (
                    "iso_microscopic_with_verified_bush_domain_extension"
                ),
                convention="mapped_unmixed_verified_full_domain_extension",
                direction_resolved=True,
            ),
            mode_identity=column.mode_identity,
            microscopic_provenance=column.provenance,
            microscopic_domain_extension=(
                column.domain_extension_evidence
                if extension is None
                else extension
            ),
            label=f"mode-{index}",
        )

    raw = [_raw(column, index) for index, column in enumerate(mapped)]
    operation = ExactSeitzOperation.from_values(np.eye(3), (0, 0, 0))
    atoms = (
        ChildAtom.from_values("atom-0", (0, 0, 0), "Fe"),
        ChildAtom.from_values("atom-1", (Fraction(1, 2), 0, 0), "Fe"),
    )
    with pytest.raises(
        DisplaciveModelError,
        match=r"invalid_domain_extension_evidence.*source basis is incomplete",
    ):
        project_modes_to_free_rows(
            (),
            atoms,
            (operation,),
            [replace(raw[0], column_index=0)],
            frame_id="child-frame",
            reference_lattice=np.eye(3),
        )

    first_extension = raw[0].microscopic_domain_extension
    assert first_extension is not None
    tampered = replace(
        raw[0],
        microscopic_domain_extension=replace(
            first_extension,
            coefficient_column=(0.25, 0.5),
        ),
    )
    with pytest.raises(
        DisplaciveModelError,
        match=r"invalid_domain_extension_evidence.*change-of-basis digest",
    ):
        project_modes_to_free_rows(
            (),
            atoms,
            (operation,),
            [tampered, raw[1]],
            frame_id="child-frame",
            reference_lattice=np.eye(3),
        )

    assert replace(
        raw[0],
        microscopic_domain_extension=None,
    ).source_binding_issue() == "microscopic scale source claims missing domain extension"
    assert replace(
        raw[0],
        scale_provenance=ModeScaleProvenance(
            source="iso_microscopic_display_distortion",
            convention="mapped_unmixed_fractional_source_column",
            direction_resolved=True,
        ),
    ).source_binding_issue() == "microscopic scale source omits verified domain extension"


def test_domain_extension_rejects_rank_loss_residual_and_ill_conditioned_change():
    parent, reference, canonical, _selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    lattice = np.asarray(parent.lattice.matrix, dtype=float)
    proof = _symbolic_subspace_proof(reference, canonical, lattice)
    assert proof.validation.matched

    rank_lost = replace(
        proof,
        reference_comparison=(
            np.asarray([[1.0, 0.0, 0.0]]),
            np.asarray([[0.0, 0.0, 0.0]]),
        ),
    )
    with pytest.raises(ValueError, match="reference_common_rank_loss"):
        _extend_microscopic_modes_to_reference_domain(
            reference, canonical, lattice, rank_lost,
        )

    residual = replace(
        proof,
        candidate_comparison=(
            np.asarray([[1.0, 0.0, 1.0]]),
            np.asarray([[0.0, 1.0, 0.0]]),
        ),
    )
    with pytest.raises(ValueError, match="residual_exceeds_tolerance"):
        _extend_microscopic_modes_to_reference_domain(
            reference, canonical, lattice, residual,
        )

    ill_conditioned = replace(
        proof,
        reference_comparison=(
            np.asarray([[1.0, 0.0, 0.0]]),
            np.asarray([[0.0, 1.0e-5, 0.0]]),
        ),
        candidate_comparison=(
            np.asarray([[1.0e-5, 0.0, 0.0]]),
            np.asarray([[0.0, 1.0, 0.0]]),
        ),
    )
    with pytest.raises(ValueError, match="change_of_basis_ill_conditioned"):
        _extend_microscopic_modes_to_reference_domain(
            reference, canonical, lattice, ill_conditioned,
        )


def test_mapper_preserves_source_provenance_and_signed_scalars():
    block = _block()
    modes = [_mode_from_block(block, 0), _mode_from_block(block, 1)]
    modes = [
        replace(mode, mode_identity=replace(mode.mode_identity, status="verified", reason=None))
        for mode in modes
    ]
    first = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    second = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    mapper = DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        modes[0].amplitude_key: first,
        modes[1].amplitude_key: second,
    }
    source = mapper.map_microscopic_columns_to_supercell(
        Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]]),
        [],
        modes,
        np.eye(3).tolist(),
        subgroup_context=_subgroup(),
        frame_id="child-frame",
        atom_ids=("a", "b"),
    )
    assert [column.mode_identity for column in source] == [
        mode.mode_identity for mode in modes
    ]
    assert [column.provenance.exact_vector_digest for column in source] == [
        mode.microscopic_provenance.exact_vector_digest for mode in modes
    ]
    assert [column.provenance.source_float_column_digest for column in source] == [
        mode.microscopic_provenance.source_float_column_digest for mode in modes
    ]
    assert all(column.provenance.unresolved_fields == () for column in source)
    assert all(
        column.provenance.source_frame_id == "child-frame"
        and column.provenance.source_atom_order_id
        == microscopic_atom_order_id(("a", "b"))
        for column in source
    )

    current = [
        _mapped(
            "b",
            source[1].provenance,
            -3.0 * second,
            canonical=False,
            amplitude_key=source[1].amplitude_key,
            mode_identity=source[1].mode_identity,
        ),
        _mapped(
            "a",
            source[0].provenance,
            2.0 * first,
            canonical=False,
            amplitude_key=source[0].amplitude_key,
            mode_identity=source[0].mode_identity,
        ),
    ]
    verified = verify_mapped_microscopic_columns(
        source,
        current,
        [[2.0, 0.3, 0.0], [0.0, 3.0, 0.2], [0.0, 0.0, 4.0]],
    )
    assert [column.mode_id for column in verified] == [
        source[0].mode_id,
        source[1].mode_id,
    ]
    assert [column.current_over_canonical_signed_scale for column in verified] == pytest.approx([2.0, -3.0])
    assert all(column.relation_status == "scalar_verified" for column in verified)


def test_mapper_rejects_tampered_direct_displacement_and_outer_routing_fields():
    block = _block()
    modes = [
        replace(
            _mode_from_block(block, column),
            mode_identity=replace(
                _mode_from_block(block, column).mode_identity,
                status="verified",
                reason=None,
            ),
        )
        for column in range(2)
    ]
    mapper = DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.ones((2, 3), dtype=float) for mode in modes
    }

    def _map(candidate_modes, subgroup=None):
        return mapper.map_microscopic_columns_to_supercell(
            Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]]),
            [],
            candidate_modes,
            np.eye(3).tolist(),
            subgroup_context=subgroup or _subgroup(),
            frame_id="child-frame",
            atom_ids=("atom-0", "atom-1"),
        )

    damaged_rows = list(modes[0].bush_modes)
    damaged_rows[0] = replace(
        damaged_rows[0],
        displacements=[[9.0, 0.0, 0.0]],
    )
    with pytest.raises(ValueError, match="direct microscopic displacements"):
        _map([replace(modes[0], bush_modes=damaged_rows), modes[1]])

    with pytest.raises(ValueError, match="direction differs"):
        _map([replace(modes[0], opd_dir_raw="(a,a)"), modes[1]])

    with pytest.raises(ValueError, match="mode OPD differs"):
        _map([replace(modes[0], opd_symbol="P2"), modes[1]])

    wrong_bush_opd = list(modes[0].bush_modes)
    wrong_bush_opd[0] = replace(wrong_bush_opd[0], opd_symbol="P2")
    with pytest.raises(ValueError, match="row routing fields differ"):
        _map([replace(modes[0], bush_modes=wrong_bush_opd), modes[1]])

    routing_changes = (
        {"irrep_label": "WRONG"},
        {"wyckoff_site": "b"},
        {"wyckoff_orbit_id": "orbit-b"},
        {"site_irrep": "B"},
        {"opd_component": "z"},
        {"dimension": 9},
        {"k_coords_label": "0,0,0"},
    )
    for changes in routing_changes:
        with pytest.raises(ValueError, match="routing fields differ"):
            _map([replace(modes[0], **changes), modes[1]])

    forged_parent_identity = replace(modes[0].mode_identity, parent_sg=63)
    with pytest.raises(ValueError, match="source provenance disagree"):
        _map([
            replace(modes[0], mode_identity=forged_parent_identity),
            modes[1],
        ])

    forged_k_identity = replace(
        modes[0].mode_identity,
        k_coordinates=("0", "0", "0"),
    )
    with pytest.raises(ValueError, match="source provenance disagree"):
        _map([
            replace(
                modes[0],
                k_coords_label="0,0,0",
                mode_identity=forged_k_identity,
            ),
            modes[1],
        ])

    damaged_row_identity = list(modes[0].bush_modes)
    damaged_row_identity[0] = replace(damaged_row_identity[0], wyckoff_letter="b")
    with pytest.raises(ValueError, match="row routing fields differ"):
        _map([replace(modes[0], bush_modes=damaged_row_identity), modes[1]])

    context_changes = (
        {"space_group_number": 2},
        {"basis_vectors": [[2, 0, 0], [0, 1, 0], [0, 0, 1]]},
        {"origin": [1, 0, 0]},
        {"irrep_label": "GM2+"},
        {"opd_symbol": "P2"},
        {"opd_dir_raw": "(0,a)"},
    )
    for changes in context_changes:
        with pytest.raises(ValueError, match=r"query context|primary direction"):
            _map(modes, replace(_subgroup(), **changes))


def test_constant_point_raw_is_authoritative_and_numeric_tamper_fails_closed():
    block = MicroscopicVectorBlock(
        global_irrep="GM1+",
        wyckoff_letter="a",
        site_irrep="A",
        rows=(
            MicroscopicVectorRow(
                point_raw=("0", "0", "0"),
                displacements=((Fraction(1), Fraction(0), Fraction(0)),),
            ),
        ),
        source_order=0,
        direction_symbol="VECTOR,A",
    )
    mode = _mode_from_block(block, 0)
    gamma_provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=QUERY_DIGEST,
        query_order=0,
        query_irrep_label="GM1+",
        source_parent_sg=62,
        source_k_coordinates=("0", "0", "0"),
        source_subgroup_space_group_number=1,
        source_subgroup_basis=np.eye(3, dtype=int).tolist(),
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block,
        column_index=0,
    )
    mode = replace(
        mode,
        k_coords_label="0,0,0",
        mode_identity=replace(
            mode.mode_identity,
            k_coordinates=("0", "0", "0"),
            status="verified",
            reason=None,
        ),
        microscopic_provenance=gamma_provenance,
        bush_modes=[replace(mode.bush_modes[0], point=[0.25, 0.0, 0.0])],
    )
    mapper = DistortionMapper()
    with pytest.raises(ValueError, match="differs from authoritative point_raw"):
        mapper.map_microscopic_columns_to_supercell(
            Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]]),
            [{
                "equivalent_indices": [0],
                "wyckoff_letter": "a",
                "orbit_id": "orbit-a",
            }],
            [mode],
            np.eye(3).tolist(),
            subgroup_context=_subgroup(irrep="GM1+"),
            frame_id="child-frame",
            atom_ids=("atom-0",),
        )


def test_mapper_rejects_mixing_missing_frame_and_atom_order_mismatches():
    p0, p1 = _provenance(0), _provenance(1)
    first = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    second = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    canonical = [
        _mapped("m0", p0, first, canonical=True),
        _mapped("m1", p1, second, canonical=True),
    ]
    lattice = np.eye(3)

    mixed = [
        _mapped("m0", p0, first + second, canonical=False),
        _mapped("m1", p1, first - second, canonical=False),
    ]
    with pytest.raises(ValueError, match="non_diagonal_column_mixing"):
        verify_mapped_microscopic_columns(canonical, mixed, lattice)
    with pytest.raises(ValueError, match="missing or unexpected"):
        verify_mapped_microscopic_columns(canonical, mixed[:1], lattice)

    wrong_frame = [
        _mapped("m0", p0, first, canonical=False, frame_id="other-frame"),
        _mapped("m1", p1, second, canonical=False),
    ]
    with pytest.raises(ValueError, match="frame mismatch"):
        verify_mapped_microscopic_columns(canonical, wrong_frame, lattice)

    wrong_order = [
        _mapped("m0", p0, first[::-1], canonical=False, atom_ids=("b", "a")),
        _mapped("m1", p1, second, canonical=False),
    ]
    with pytest.raises(ValueError, match="atom-order mismatch"):
        verify_mapped_microscopic_columns(canonical, wrong_order, lattice)


def test_mapped_column_rejects_incomplete_or_mismatched_mode_identity():
    provenance = _provenance(0)
    vector = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    verified = ModeIdentity(
        parent_sg=62,
        global_irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        wyckoff_letter="a",
        orbit_id="orbit-a",
        site_irrep="A",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="verified",
    )

    with pytest.raises(ValueError, match="not verified"):
        _mapped(
            "a",
            provenance,
            vector,
            canonical=True,
            mode_identity=replace(verified, status="unresolved", reason="missing"),
        )
    with pytest.raises(ValueError, match="source provenance disagree"):
        _mapped(
            "a",
            provenance,
            vector,
            canonical=True,
            mode_identity=replace(verified, global_irrep="WRONG"),
        )


def test_mapped_column_keeps_legacy_positional_fields_and_checks_actual_atom_order():
    original = _mapped(
        "a",
        _provenance(0),
        np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        canonical=True,
    )
    positional = MappedMicroscopicColumn(
        original.mode_id,
        original.amplitude_key,
        original.frame_id,
        original.atom_order_id,
        original.displacements,
        original.provenance,
        original.relation_status,
        original.current_over_canonical_signed_scale,
        original.cartesian_residual_angstrom,
        atom_ids=original.atom_ids,
        mode_identity=original.mode_identity,
    )

    assert positional.mode_identity == original.mode_identity
    assert positional.atom_ids == original.atom_ids
    with pytest.raises(ValueError, match="atom-order ID does not match"):
        replace(original, atom_order_id="forged-order")


def test_verifier_rejects_display_identity_change_even_when_stable_token_matches():
    provenance = _provenance(0)
    vector = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    canonical = _mapped("a", provenance, vector, canonical=True)
    changed_identity = replace(canonical.mode_identity, component_label="renamed")
    assert changed_identity.stable_token == canonical.mode_id
    current = _mapped(
        "renamed",
        provenance,
        vector,
        canonical=False,
        mode_identity=changed_identity,
        amplitude_key=canonical.amplitude_key,
    )

    with pytest.raises(ValueError, match="mode identity mismatch"):
        verify_mapped_microscopic_columns((canonical,), (current,), np.eye(3))
