"""Audited ISO microscopic/site-irrep and macroscopic-tensor evidence tests."""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from collections import Counter
from dataclasses import replace
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from backend.models.iso_mode_models import (
    InvariantDirection,
    IsoPointOperationCatalog,
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
    OrientedSiteCharacterTable,
    SiteIrrepCharacter,
    SiteOperation,
    build_site_irrep_projectors,
    validate_cartesian_subspaces,
)
from backend.utils.exceptions import WrapperRunError
from backend.utils.opd_format import (
    _g_allowed,
    _k_equivalent_exact,
    k_star_fraction_vectors,
)
from backend.utils.parent_header import parent_wyckoff_display
from backend.utils.text_parser import (
    parse_invariant_direction_table,
    parse_macroscopic_tensor_blocks,
    parse_microscopic_vector_blocks,
    parse_symbolic_tensor_expression,
    parse_wyckoff_character_tables,
)
from backend.wrappers.iso_wrapper import (
    BushMode,
    DistortionMode,
    IsoWrapper,
    SubgroupInfo,
)
from features.input_cif import read_cif
from features.input_cif.symmetry_validator import SymmetryValidator
from features.method2.superspace import (
    ParametricModeResult,
    _common_point_rows,
    _display_site_label,
    _instantiate_bush_modes_by_orbit,
    _match_candidate_points_to_reference,
    _ordered_complete_microscopic_space,
    _proven_real_transport_phase,
    _replace_rootless_supplement_basis,
    _select_verified_microscopic_modes,
    _site_irrep_for_disp,
    validate_symbolic_mode_subspaces,
)

ISODISTORT_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ISODISTORT_ROOT.parent


def _character_table(row: dict) -> OrientedSiteCharacterTable:
    return OrientedSiteCharacterTable(
        wyckoff_letter=row["wyckoff_letter"],
        point_group=row["point_group"],
        operations=tuple(SiteOperation(**item) for item in row["operations"]),
        irreps=tuple(SiteIrrepCharacter(**item) for item in row["irreps"]),
    )


def _matrix(projectors, label: str) -> np.ndarray:
    item = next(projector for projector in projectors if projector.label == label)
    return np.asarray(item.matrix, dtype=float)


def _mode(
    points: list[tuple[str, str, str]],
    vectors: list[tuple[float, float, float]],
    *,
    irrep: str = "N1+",
    letter: str = "e",
    site_irrep: str = "A1",
    k_coordinates: tuple[str, str, str] = ("1/2", "0", "1/2"),
    source_block_order: int = 0,
) -> DistortionMode:
    bushes = [
        BushMode(
            irrep_label=irrep,
            opd_symbol="4D1",
            wyckoff_letter=letter,
            point=[0.0, 0.0, 0.0],
            point_raw=list(point),
            displacements=[list(vector)],
        )
        for point, vector in zip(points, vectors, strict=True)
    ]
    exact_block = MicroscopicVectorBlock(
        irrep,
        letter,
        site_irrep,
        tuple(
            MicroscopicVectorRow(
                point_raw=point,
                displacements=(tuple(Fraction(str(value)) for value in vector),),
            )
            for point, vector in zip(points, vectors, strict=True)
        ),
        source_order=source_block_order,
        direction_symbol="VECTOR,A",
    )
    try:
        exact_source_k = tuple(str(Fraction(value)) for value in k_coordinates)
    except (ValueError, ZeroDivisionError):
        exact_source_k = ()
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=hashlib.sha256(f"synthetic:{irrep}".encode()).hexdigest(),
        query_order=0,
        query_irrep_label=irrep,
        source_parent_sg=139,
        source_k_coordinates=exact_source_k,
        source_subgroup_space_group_number=2,
        source_subgroup_basis=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="N1+",
        source_subgroup_opd_symbol="4D1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=exact_block,
        column_index=0,
    )
    return DistortionMode(
        irrep_label=irrep,
        dimension=1,
        wyckoff_site=letter,
        bush_modes=bushes,
        amplitude_key=f"{irrep}__{letter}__test",
        site_irrep=site_irrep,
        k_coords_label=",".join(k_coordinates),
        mode_identity=ModeIdentity(
            parent_sg=139,
            global_irrep=irrep,
            k_coordinates=k_coordinates,
            wyckoff_letter=letter,
            site_irrep=site_irrep,
            component_index=0,
            component_label="a",
            source="iso_microscopic",
            status="verified",
        ),
        microscopic_provenance=provenance,
    )


def test_parse_exact_embedding_directions():
    rows = parse_invariant_direction_table(
        """
        *********Irrep (ML) Dir       Subgroup     Size
        GM1+       (a)       139 I4/mmm   1
        GM5+       (a,b)     2 P-1        1
        N1+        (a,b,c,d) 2 P-1        8
        *
        """
    )
    assert rows == [
        {
            "irrep_label": "GM1+",
            "direction_raw": "(a)",
            "subgroup_number": 139,
            "subgroup_symbol": "I4/mmm",
            "size": 1,
        },
        {
            "irrep_label": "GM5+",
            "direction_raw": "(a,b)",
            "subgroup_number": 2,
            "subgroup_symbol": "P-1",
            "size": 1,
        },
        {
            "irrep_label": "N1+",
            "direction_raw": "(a,b,c,d)",
            "subgroup_number": 2,
            "subgroup_symbol": "P-1",
            "size": 8,
        },
    ]


def test_parse_parametric_embedding_directions_retains_exact_k_context():
    rows = parse_invariant_direction_table(
        """
        *************Irrep (ML) k params k vector  Dir   Subgroup   Size
        GM1+                (0,0,0)   (a)   139 I4/mmm 1
        LD1        1/12     (0,0,1/6) (a,b) 99 P4mm    12
        LD1        1/6      (0,0,1/3) (a,b) 99 P4mm    6
        *
        """
    )

    assert rows[0]["k_coordinates"] == ("0", "0", "0")
    assert rows[0]["k_parameters"] == ()
    assert rows[1]["k_parameters"] == ("1/12",)
    assert rows[1]["k_coordinates"] == ("0", "0", "1/6")
    assert rows[2]["k_parameters"] == ("1/6",)
    assert rows[2]["k_coordinates"] == ("0", "0", "1/3")


def test_microscopic_parser_preserves_fractions_rows_and_repeated_blocks():
    rows = parse_microscopic_vector_blocks(
        """
        Irrep (ML) Wyckoff Irrep Point        Projected Vectors
        N1+        g       A1    (0,1/2,z)    (0,0,1/3), (1,0,0)
                                 (1,1/2,z+1)  (0,0,-1/3), (-1,0,0)
        N1+        g       A1    (0,1/2,z)    (0,1/2,0), (0,-1/2,0)
                                 (1,1/2,z+1)  (0,-1/2,0), (0,1/2,0)
        *
        """
    )
    assert len(rows) == 2
    assert [row["source_order"] for row in rows] == [0, 1]
    assert all(row["site_irrep"] == "A1" for row in rows)
    assert len(rows[0]["rows"]) == 2
    assert rows[0]["rows"][0]["point_raw"] == ("0", "1/2", "z")
    assert rows[0]["rows"][0]["displacements"][0] == (
        Fraction(0), Fraction(0), Fraction(1, 3)
    )
    assert rows[1]["rows"][1]["displacements"][1] == (
        Fraction(0), Fraction(1, 2), Fraction(0)
    )


def test_microscopic_parser_preserves_decimal_quantization_intervals():
    rows = parse_microscopic_vector_blocks(
        "GM1+ a A1 (0,0,0) (0.577,-1.155,1/2)"
    )

    assert rows[0]["rows"][0]["displacements"][0] == (
        Fraction(577, 1000),
        Fraction(-231, 200),
        Fraction(1, 2),
    )
    assert rows[0]["rows"][0]["displacement_half_steps"][0] == (
        Fraction(1, 2000),
        Fraction(1, 2000),
        Fraction(0),
    )


def test_iso_microscopic_block_order_defines_copy_and_component_identity():
    row = MicroscopicVectorRow(
        point_raw=("0", "0", "z"),
        displacements=(
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
        ),
    )
    blocks = [
        MicroscopicVectorBlock(
            "GM5+", "g", "A1", (row,), source_order=index,
            direction_symbol="VECTOR,A,B",
        )
        for index in range(2)
    ]
    wrapper = object.__new__(IsoWrapper)
    wrapper.get_microscopic_vector_blocks = lambda *_args, **_kwargs: blocks
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM5+", "(a;b)", 2, "P-1", 8)
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "4D1"
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=2,
        opd_symbol="4D1",
        irrep_label="GM5+",
        parent_sg=139,
        k_coordinates=["0", "0", "0"],
        basis_vectors=np.eye(3).tolist(),
        origin=[0.0, 0.0, 0.0],
    )
    modes = IsoWrapper.calc_microscopic_distortion_modes(
        wrapper,
        139,
        subgroup,
        ["g"],
        irrep_labels=["GM5+"],
    )
    assert [mode.site_irrep for mode in modes] == ["A1_1", "A1_1", "A1_2", "A1_2"]
    assert [mode.opd_component for mode in modes] == ["a", "b", "a", "b"]
    assert [mode.mode_identity.copy_index for mode in modes] == [1, 1, 2, 2]
    assert [mode.mode_identity.component_index for mode in modes] == [0, 1, 0, 1]
    assert all(mode.mode_identity.status == "partial" for mode in modes)
    assert all(
        mode.mode_identity.reason == "awaiting_bush_or_smodes_subspace_validation"
        for mode in modes
    )
    assert len({mode.mode_identity.stable_token for mode in modes}) == 4


def test_oriented_c2v_and_d2h_polar_vector_projectors():
    text = (
        "c D2h (E|0,0,0) Ag 1, B3g 1, B1g 1, B2g 1, Au 1, B3u 1, "
        "B1u 1, B2u 1, (C2x|0,1,0) Ag 1, B3g 1, B1g -1, B2g -1, "
        "Au 1, B3u 1, B1u -1, B2u -1, (C2y|0,0,0) Ag 1, B3g -1, "
        "B1g -1, B2g\n"
        "             1, Au 1, B3u -1, B1u -1, B2u 1, (C2z|0,1,0) "
        "Ag 1, B3g -1, B1g 1, B2g -1, Au 1, B3u -1, B1u 1, B2u -1, "
        "(I|0,1,0) Ag 1, B3g 1, B1g 1, B2g 1, Au -1, B3u -1, B1u -1, "
        "B2u -1,\n"
        "             (SGx|0,0,0) Ag 1, B3g 1, B1g -1, B2g -1, Au -1, "
        "B3u -1, B1u 1, B2u 1, (SGy|0,1,0) Ag 1, B3g -1, B1g -1, "
        "B2g 1, Au -1, B3u 1, B1u 1, B2u -1, (SGz|0,0,0) Ag 1, "
        "B3g -1, B1g 1,\n"
        "             B2g -1, Au -1, B3u 1, B1u -1, B2u 1\n"
        "g C2v (E|0,0,0) A1 1, B2 1, A2 1, B1 1, (C2z|0,1,0) A1 1, "
        "B2 -1, A2 1, B1 -1, (SGx|0,0,0) A1 1, B2 1, A2 -1, B1 -1, "
        "(SGy|0,1,0) A1 1, B2 -1, A2 -1, B1 1\n*\n"
    )
    tables = parse_wyckoff_character_tables(text)
    assert [(row["wyckoff_letter"], row["point_group"]) for row in tables] == [
        ("c", "D2h"), ("g", "C2v")
    ]
    catalog = IsoPointOperationCatalog.from_data_space(
        ISODISTORT_ROOT / "resources" / "isobyu" / "data_space.txt"
    )
    d2h = build_site_irrep_projectors(_character_table(tables[0]), catalog)
    c2v = build_site_irrep_projectors(_character_table(tables[1]), catalog)
    assert _matrix(c2v, "A1") == pytest.approx(np.diag([0.0, 0.0, 1.0]))
    assert _matrix(c2v, "B1") == pytest.approx(np.diag([1.0, 0.0, 0.0]))
    assert _matrix(c2v, "B2") == pytest.approx(np.diag([0.0, 1.0, 0.0]))
    assert _matrix(c2v, "A2") == pytest.approx(np.zeros((3, 3)))
    assert _matrix(d2h, "B3u") == pytest.approx(np.diag([1.0, 0.0, 0.0]))
    assert _matrix(d2h, "B2u") == pytest.approx(np.diag([0.0, 1.0, 0.0]))
    assert _matrix(d2h, "B1u") == pytest.approx(np.diag([0.0, 0.0, 1.0]))


def test_macroscopic_parser_preserves_copies_components_and_decimal_evidence():
    tetragonal = """
    *******Irrep (ML) Basis Functions
    GM1+       xx+yy
               zz
    GM2+       xx-yy
    GM4+       xy
    GM5+       yz,-xz
    *
    """
    rows = parse_macroscopic_tensor_blocks(tetragonal)
    assert [(row["global_irrep"], len(row["components"])) for row in rows] == [
        ("GM1+", 1), ("GM1+", 1), ("GM2+", 1), ("GM4+", 1), ("GM5+", 2)
    ]
    assert rows[0]["components"][0]["coefficients"] == (
        Fraction(1), Fraction(1), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),
    )
    assert rows[-1]["components"][1]["coefficients"] == (
        Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(-1), Fraction(0),
    )

    cubic = parse_symbolic_tensor_expression("1.732xx-1.732yy")
    assert cubic["quality"] == "printed_decimal"
    assert cubic["coefficients"][:2] == (Decimal("1.732"), Decimal("-1.732"))
    radical = parse_symbolic_tensor_expression("sqrt(3)xx")
    assert radical["quality"] == "unresolved"
    assert radical["coefficients"] is None


def test_macroscopic_wrapper_assigns_one_based_copy_indices_and_fails_closed():
    output = """
    Irrep (ML) Basis Functions
    GM1+       xx+yy
               zz
    GM2+       xx-yy
    GM4+       xy
    GM5+       yz,-xz
    *
    """
    wrapper = object.__new__(IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: output
    basis = IsoWrapper.get_macroscopic_tensor_basis(wrapper, 139)
    assert basis.status == "verified"
    assert [(block.global_irrep, block.copy_index) for block in basis.blocks] == [
        ("GM1+", 1), ("GM1+", 2), ("GM2+", None),
        ("GM4+", None), ("GM5+", None),
    ]

    wrapper._run_session = lambda *_args, **_kwargs: "no table"
    missing = IsoWrapper.get_macroscopic_tensor_basis(wrapper, 139)
    assert missing.status == "unresolved"
    assert missing.blocks == ()
    assert missing.reason == "iso_macroscopic_tensor_output_empty"

    def _unavailable(*_args, **_kwargs):
        raise WrapperRunError("iso", 1, "unavailable")

    wrapper._run_session = _unavailable
    unavailable = IsoWrapper.get_macroscopic_tensor_basis(wrapper, 139)
    assert unavailable.status == "unresolved"
    assert unavailable.blocks == ()
    assert unavailable.reason == "iso_macroscopic_tensor_output_unavailable"


def test_cartesian_subspace_validation_accepts_basis_rotation_and_rejects_mismatch():
    lattice = np.asarray([[2.0, 0.3, 0.0], [0.0, 3.0, 0.2], [0.0, 0.0, 4.0]])
    x = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    y = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    rotated = [x + y, x - y]
    matched = validate_cartesian_subspaces([x, y], rotated, lattice)
    assert matched.matched
    assert matched.reference_rank == matched.candidate_rank == 2
    mismatch = validate_cartesian_subspaces([x, y], [x], lattice)
    assert not mismatch.matched
    assert mismatch.reason == "rank_mismatch"


def _rootless_basis_fixture():
    orbit_id = "physical-orbit-a"
    x = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    y = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    reference_modes = [
        DistortionMode(
            irrep_label="GM1+",
            wyckoff_site="e",
            wyckoff_orbit_id=orbit_id,
            amplitude_key=key,
            mode_identity=ModeIdentity.unresolved(
                parent_sg=139,
                global_irrep="GM1+",
                k_coordinates=("0", "0", "0"),
                wyckoff_letter="e",
                orbit_id=orbit_id,
                reason="smodes_basis_has_no_iso_microscopic_identity",
            ),
        )
        for key in ("smodes-1", "smodes-2")
    ]
    supplement = ParametricModeResult(
        modes=reference_modes,
        supercell_displacements={"smodes-1": x + y, "smodes-2": x - y},
        labels={"smodes-1": "unresolved 1", "smodes-2": "unresolved 2"},
        nmod=0,
    )
    microscopic = []
    for key, site_irrep, source_order in (
        ("canonical-x", "A1", 0),
        ("canonical-y", "B1", 1),
    ):
        mode = _mode(
            [("x", "0", "z")],
            [(1.0, 0.0, 0.0)],
            irrep="GM1+",
            site_irrep=site_irrep,
            k_coordinates=("0", "0", "0"),
            source_block_order=source_order,
        )
        microscopic.append(replace(
            mode,
            amplitude_key=key,
            wyckoff_orbit_id=orbit_id,
            mode_identity=replace(mode.mode_identity, orbit_id=orbit_id),
        ))
    return supplement, microscopic, {"canonical-x": x, "canonical-y": y}


def test_rootless_whole_space_replacement_accepts_nondiagonal_basis_change():
    supplement, microscopic, arrays = _rootless_basis_fixture()
    lattice = np.asarray([[2.0, 0.4, 0.0], [0.0, 3.0, 0.3], [0.0, 0.0, 4.0]])

    replaced = _replace_rootless_supplement_basis(
        supplement,
        microscopic,
        arrays,
        lattice,
    )

    assert [mode.amplitude_key for mode in replaced.modes] == [
        "canonical-x",
        "canonical-y",
    ]
    assert set(replaced.supercell_displacements) == {"canonical-x", "canonical-y"}
    assert np.array_equal(replaced.supercell_displacements["canonical-x"], arrays["canonical-x"])
    assert all(mode.mode_identity.status == "verified" for mode in replaced.modes)


def test_rootless_whole_space_replacement_rejects_cross_orbit_and_missing_column():
    supplement, microscopic, arrays = _rootless_basis_fixture()
    crossed = [
        microscopic[0],
        replace(
            microscopic[1],
            wyckoff_orbit_id="physical-orbit-b",
            mode_identity=replace(
                microscopic[1].mode_identity,
                orbit_id="physical-orbit-b",
            ),
        ),
    ]
    with pytest.raises(ValueError, match="cross or omit physical orbits"):
        _replace_rootless_supplement_basis(
            supplement,
            crossed,
            arrays,
            np.eye(3),
        )
    with pytest.raises(ValueError, match="complete child-fixed space"):
        _replace_rootless_supplement_basis(
            supplement,
            microscopic[:1],
            {"canonical-x": arrays["canonical-x"]},
            np.eye(3),
        )


def test_rootless_microscopic_queries_split_same_digest_by_query_order():
    _supplement, microscopic, _arrays = _rootless_basis_fixture()
    first = microscopic[0]
    second = microscopic[1]
    shared_digest = first.microscopic_provenance.query_digest
    second_provenance = replace(
        second.microscopic_provenance,
        query_digest=shared_digest,
        query_order=1,
        query_irrep_label="GM2+",
        source_global_irrep="GM2+",
    )
    second_identity = replace(second.mode_identity, global_irrep="GM2+")
    second = replace(
        second,
        irrep_label="GM2+",
        microscopic_provenance=second_provenance,
        mode_identity=second_identity,
    )

    ordered = _ordered_complete_microscopic_space([second, first])

    assert [mode.irrep_label for mode in ordered] == ["GM1+", "GM2+"]

    duplicate_order = replace(
        second,
        microscopic_provenance=replace(second_provenance, query_order=0),
    )
    with pytest.raises(ValueError, match="multiple irreps"):
        _ordered_complete_microscopic_space([first, duplicate_order])

    missing_source_column = replace(
        first,
        microscopic_provenance=replace(
            first.microscopic_provenance,
            source_block_column_count=2,
        ),
    )
    with pytest.raises(ValueError, match="source_block_column_missing"):
        _ordered_complete_microscopic_space([missing_source_column])


def test_symbolic_point_alignment_handles_permutation_and_gamma_integer_shifts():
    reference = [
        _mode(
            [("x", "0", "-z"), ("y", "0", "z")],
            [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        )
    ]
    permuted = [
        _mode(
            [("y", "0", "z"), ("x", "0", "-z")],
            [(0.0, 1.0, 0.0), (1.0, 0.0, 0.0)],
        )
    ]
    assert validate_symbolic_mode_subspaces(
        reference, permuted, np.eye(3)
    ).matched

    shifted_reference = [
        _mode(
            [("x", "0", "-z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    equivalent = [
        _mode(
            [("x+1", "0", "1-z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    assert validate_symbolic_mode_subspaces(
        shifted_reference, equivalent, np.eye(3)
    ).matched


def test_non_gamma_exact_common_rows_prove_nondiagonal_whole_space_match():
    reference_points = [
        ("x", "0", "z"),
        ("y", "0", "z"),
        ("y+1", "0", "z"),
    ]
    candidate_points = [
        ("x", "0", "z"),
        ("y", "0", "z"),
        ("y+2", "0", "z"),
    ]
    reference = [
        _mode(
            reference_points,
            [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
        ),
        _mode(
            reference_points,
            [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
            site_irrep="B1",
        ),
    ]
    candidate = [
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        ),
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)],
            site_irrep="B1",
        ),
    ]

    validation = validate_symbolic_mode_subspaces(
        reference, candidate, np.eye(3)
    )

    assert validation.matched
    assert validation.reference_rank == validation.candidate_rank == 2


def test_non_gamma_rank_losing_common_rows_still_require_phase_proof():
    reference_points = [("x", "0", "z"), ("y", "0", "z")]
    candidate_points = [("x", "0", "z"), ("y+1", "0", "z")]
    reference = [
        _mode(reference_points, [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)]),
        _mode(
            reference_points,
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            site_irrep="B1",
        ),
    ]
    candidate = [
        _mode(candidate_points, [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)]),
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            site_irrep="B1",
        ),
    ]

    with pytest.raises(ValueError, match="star-arm phase evidence"):
        validate_symbolic_mode_subspaces(reference, candidate, np.eye(3))


def test_x_star_common_minus_phase_extends_sparse_microscopic_column():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z")],
            [(1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x+1", "0", "z")],
            [(-1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]

    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139,
    )

    assert diagnostics[("X1+", "e")].matched
    assert selected[0].mode_identity.status == "verified"
    assert selected[0].microscopic_domain_extension is not None
    assert (
        selected[0].microscopic_domain_extension.comparison_domain_kind
        == "phase_proven_transported_points"
    )
    assert selected[0].bush_modes[0].point_raw == ["x", "0", "z"]
    assert selected[0].bush_modes[0].displacements == [[1.0, 0.0, 0.0]]


def test_star_phase_matching_assigns_distinct_equivalent_bush_images():
    reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x-1", "0", "z"), ("x", "0", "z")],
            [(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]

    validation = validate_symbolic_mode_subspaces(
        reference, candidate, np.eye(3)
    )
    assert validation.matched
    assert validation.reference_rank == validation.candidate_rank == 1
    reference_points, reference_rows = _common_point_rows(reference)
    candidate_points, _candidate_rows = _common_point_rows(candidate)
    transports = _match_candidate_points_to_reference(
        reference_points,
        reference_rows,
        candidate_points,
        parent_sg=139,
        k_coordinates=(Fraction(1, 2), Fraction(1, 2), Fraction(0)),
    )
    assert [point for point, _translation in transports] == [
        reference_points[1],
        reference_points[0],
    ]
    assert [translation for _point, translation in transports] == [
        (-2, 0, 0),
        (0, 0, 0),
    ]


def test_symbolic_point_matching_fails_when_no_complete_injection_exists():
    reference_points = [("x", "0", "z"), ("x+1", "0", "z")]
    candidate_points = [
        ("x-1", "0", "z"),
        ("x+2", "0", "z"),
        ("x+3", "0", "z"),
    ]
    reference = [
        _mode(
            reference_points,
            [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        ),
        _mode(
            reference_points,
            [(0.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
            site_irrep="B1",
            source_block_order=1,
            k_coordinates=("0", "0", "0"),
        ),
    ]
    candidate = [
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0)] * 3,
            k_coordinates=("0", "0", "0"),
        ),
        _mode(
            candidate_points,
            [(0.0, 1.0, 0.0)] * 3,
            site_irrep="B1",
            source_block_order=1,
            k_coordinates=("0", "0", "0"),
        ),
    ]

    with pytest.raises(ValueError, match="non-injectively"):
        validate_symbolic_mode_subspaces(reference, candidate, np.eye(3))


@pytest.mark.parametrize(
    ("irrep", "k_coordinates"),
    [
        ("N1+", ("1/2", "0", "1/2")),
        ("DT1", ("1/4", "0", "0")),
    ],
)
def test_non_gamma_integer_translation_without_star_arm_identity_fails_closed(
    irrep, k_coordinates,
):
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z"), ("y", "0", "z")],
            [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            irrep=irrep,
            k_coordinates=k_coordinates,
        )
    ]
    candidate = [
        _mode(
            [("x+1", "0", "z"), ("y", "0", "z")],
            [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            irrep=irrep,
            k_coordinates=k_coordinates,
        )
    ]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert diagnostics[(irrep, "e")].reason == "subspace_validation_error"
    assert selected[0].mode_identity.status == "unresolved"


def test_r_hexagonal_reciprocal_centering_keeps_exact_star_arms_distinct():
    zero = Fraction(0)
    half = Fraction(1, 2)

    assert _g_allowed(1, 1, 0, "R")
    assert not _g_allowed(1, 0, 0, "R")
    assert _k_equivalent_exact((zero, zero, zero), (Fraction(1), Fraction(1), zero), "R")
    assert not _k_equivalent_exact(
        (zero, zero, zero), (Fraction(1), zero, zero), "R"
    )
    assert k_star_fraction_vectors((zero, half, zero), 146) == [
        (zero, half, zero),
        (half, zero, Fraction(1)),
        (Fraction(3, 2), half, zero),
    ]


def test_r_hexagonal_star_phase_proof_rejects_mixed_real_arm_phases():
    with pytest.raises(ValueError, match="star-arm phase evidence"):
        _proven_real_transport_phase(
            (-1, 0, 0),
            parent_sg=146,
            k_coordinates=(Fraction(0), Fraction(1, 2), Fraction(0)),
        )


def test_unknown_reciprocal_centering_fails_closed():
    with pytest.raises(ValueError, match="unsupported reciprocal-lattice centering"):
        _g_allowed(0, 0, 0, "Q")


def test_unparseable_symbolic_transport_k_fails_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("a", "0", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x+1", "0", "z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("a", "0", "0"),
        )
    ]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert diagnostics[("N1+", "e")].reason == "microscopic_k_identity_not_exact"
    assert selected[0].mode_identity.status == "unresolved"


def test_gamma_multiple_integer_transporters_require_column_consistency():
    consistent_reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x+2", "0", "z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    assert validate_symbolic_mode_subspaces(
        consistent_reference, candidate, np.eye(3)
    ).matched

    conflicting_reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    with pytest.raises(ValueError, match="unresolved phase"):
        validate_symbolic_mode_subspaces(
            conflicting_reference, candidate, np.eye(3)
        )


def test_unproven_symbolic_point_equivalence_fails_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [_mode([("x", "0", "z")], [(1.0, 0.0, 0.0)])]
    candidate = [_mode([("y", "0", "z")], [(1.0, 0.0, 0.0)])]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert diagnostics[("N1+", "e")].reason == "subspace_validation_error"
    assert len(selected) == len(reference)
    assert all(mode.site_irrep == "" for mode in selected)
    assert all(mode.mode_identity.status == "unresolved" for mode in selected)


def test_one_proven_star_phase_edge_is_sufficient_for_sparse_extension():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
        )
    ]
    candidate = [_mode([("x+2", "0", "z")], [(1.0, 0.0, 0.0)])]

    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )

    assert diagnostics[("N1+", "e")].matched
    assert selected[0].mode_identity.status == "verified"
    assert selected[0].microscopic_domain_extension is not None
    assert [bush.displacements for bush in selected[0].bush_modes] == [
        [[1.0, 0.0, 0.0]],
        [[-1.0, 0.0, 0.0]],
    ]


def test_subspace_rank_mismatch_and_missing_copy_fail_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    point = [("x", "0", "z")]
    reference = [
        _mode(point, [(1.0, 0.0, 0.0)]),
        _mode(point, [(0.0, 1.0, 0.0)], site_irrep="B1"),
    ]
    candidate = [_mode(point, [(1.0, 0.0, 0.0)])]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert not diagnostics[("N1+", "e")].matched
    assert diagnostics[("N1+", "e")].reason == "rank_mismatch"
    assert len(selected) == 2
    assert all(mode.site_irrep == "" for mode in selected)
    assert all(mode.mode_identity.reason == "rank_mismatch" for mode in selected)


def test_duplicate_microscopic_column_fails_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    point = [("x", "0", "z")]
    reference = [_mode(point, [(1.0, 0.0, 0.0)])]
    duplicate = [
        _mode(point, [(1.0, 0.0, 0.0)]),
        _mode(point, [(1.0, 0.0, 0.0)], source_block_order=1),
    ]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, duplicate, parent_sg=139
    )
    assert not diagnostics[("N1+", "e")].matched
    assert (
        diagnostics[("N1+", "e")].reason
        == "candidate_columns_not_independent"
    )
    assert len(selected) == 1
    assert selected[0].mode_identity.status == "unresolved"
    assert selected[0].mode_identity.reason == "candidate_columns_not_independent"


def test_five_physical_e_orbits_receive_distinct_stable_identities():
    mode = _mode([("0", "0", "z")], [(0.0, 0.0, 1.0)])
    sites = [
        {
            "wyckoff_letter": "e",
            "species": species,
            "orbit_id": f"physical-e-{index}",
            "display_label": label,
        }
        for index, (species, label) in enumerate(
            [("La", "La1"), ("La", "La2"), ("Ni", "Ni2"), ("O", "O4"), ("O", "O2")]
        )
    ]
    clones = _instantiate_bush_modes_by_orbit([mode], sites, ["e"])
    assert len(clones) == 5
    assert {item.wyckoff_orbit_id for item in clones} == {
        f"physical-e-{index}" for index in range(5)
    }
    assert {item.mode_identity.orbit_id for item in clones} == {
        f"physical-e-{index}" for index in range(5)
    }
    assert len({item.mode_identity.stable_token for item in clones}) == 5
    assert len({item.amplitude_key for item in clones}) == 5


def test_missing_or_duplicate_physical_orbit_metadata_fails_closed():
    mode = _mode([("0", "0", "z")], [(0.0, 0.0, 1.0)])
    duplicate = [
        {"wyckoff_letter": "e", "orbit_id": "same"},
        {"wyckoff_letter": "e", "orbit_id": "same"},
    ]
    with pytest.raises(ValueError, match="not unique"):
        _instantiate_bush_modes_by_orbit([mode], duplicate, ["e"])

    sites = [{"wyckoff_letter": "e", "orbit_id": "present"}]
    with pytest.raises(ValueError, match="metadata is missing"):
        _instantiate_bush_modes_by_orbit(
            [mode], sites, ["e"], requested_orbit_ids=["absent"]
        )
    with pytest.raises(ValueError, match="contain duplicates"):
        _instantiate_bush_modes_by_orbit(
            [mode], sites, ["e"], requested_orbit_ids=["present", "present"]
        )

    unresolved = _instantiate_bush_modes_by_orbit([mode], [], ["e"])
    assert len(unresolved) == 1
    assert unresolved[0].site_irrep == ""
    assert unresolved[0].mode_identity.status == "unresolved"
    assert unresolved[0].mode_identity.reason == "physical_orbit_metadata_unavailable"


def test_legacy_single_vector_site_irrep_entry_is_explicitly_unresolved():
    assert _site_irrep_for_disp("C2v", np.asarray([0.0, 0.0, 1.0])) == ""
    assert _site_irrep_for_disp("D2h", np.asarray([1.0, 0.0, 0.0])) == ""


def _wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        result = subprocess.run(
            ["wsl.exe", "--status"],  # noqa: S607 - platform executable
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return result.returncode == 0


@pytest.mark.skipif(not _wsl_available(), reason="WSL/ISO binary unavailable")
def test_live_4310_n1plus_4d1_all_192_site_irrep_labels_match_official():
    parent_path = WORKSPACE / "experiment_data" / "4310_tetra.cif"
    official_path = (
        WORKSPACE / "output_compare" / "4310_tetra.cif" / "官网" / "Method1"
        / "N1+_4D1_SG2" / "subgroup.cif"
    )
    if not parent_path.is_file() or not official_path.is_file():
        pytest.skip("validated 4310 parent/official output unavailable")

    subgroup = SubgroupInfo(
        index=0,
        space_group_number=2,
        space_group_symbol="P-1",
        subgroup_index=64,
        size=8,
        opd_symbol="4D1",
        opd_vector=[1.0, 1.0, 1.0, 1.0],
        basis_vectors=[[2, 0, 0], [0, 2, 0], [-1, -1, 1]],
        origin=[0.0, 0.0, 0.0],
        k_point_label="N",
        irrep_label="N1+",
        k_coordinates=["1/2", "0", "1/2"],
        parent_sg=139,
        opd_dir_raw="(a;b;c;d)",
        basis_raw="(2,0,0),(0,2,0),(-1,-1,1)",
        origin_raw="(0,0,0)",
    )
    parent = read_cif(parent_path)
    symmetry_info = SymmetryValidator().validate(parent)
    sites = symmetry_info["wyckoff_sites"]
    parent_wyckoff_display(parent, sites, parent_path)
    assert sum(site["wyckoff_letter"] == "e" for site in sites) == 5

    wrapper = IsoWrapper()
    bush = wrapper.calc_distortion_modes(139, subgroup, ["e", "g", "c"])
    microscopic = wrapper.calc_microscopic_distortion_modes(
        139, subgroup, ["e", "g", "c"]
    )
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, bush, microscopic, parent_sg=139
    )
    assert diagnostics
    assert all(item.matched for item in diagnostics.values())
    assert len(selected) == 96
    assert all(mode.mode_identity.status == "verified" for mode in selected)
    modes = _instantiate_bush_modes_by_orbit(selected, sites, ["e", "g", "c"])
    assert len(modes) == 192
    assert len({mode.wyckoff_orbit_id for mode in modes if mode.wyckoff_site == "e"}) == 5

    official_pattern = re.compile(
        r"^\s*\d+\s+\S+\[[^]]+\](?P<irrep>\S+?)\([^)]*\)"
        r"\[(?P<site>[^:]+):(?P<letter>[a-z]):dsp\]"
        r"(?P<site_irrep>[^\s(]+)\((?P<component>[^)]+)\)"
    )
    expected: Counter[tuple[str, str, str, str, str]] = Counter()
    for line in official_path.read_text(encoding="utf-8-sig").splitlines():
        match = official_pattern.match(line)
        if match is not None:
            expected[(
                match.group("irrep"),
                match.group("site"),
                match.group("letter"),
                match.group("site_irrep"),
                match.group("component"),
            )] += 1
    assert sum(expected.values()) == 192

    actual: Counter[tuple[str, str, str, str, str]] = Counter()
    for mode in modes:
        site_label, _species = _display_site_label(
            sites, mode.wyckoff_site, mode.wyckoff_orbit_id
        )
        actual[(
            mode.irrep_label,
            site_label,
            mode.wyckoff_site,
            mode.site_irrep,
            mode.opd_component,
        )] += 1
    assert actual == expected
