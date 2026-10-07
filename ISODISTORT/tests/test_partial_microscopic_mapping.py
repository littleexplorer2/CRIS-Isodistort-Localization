"""Pre-verification mapping tests for exact ISO microscopic columns."""

from __future__ import annotations

import contextlib
import hashlib
import io
import shutil
import subprocess
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from backend.api import IsoDistort
from backend.models.iso_mode_models import (
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
)
from backend.utils import get_config
from backend.utils.lattice import centering_primitive_matrix
from backend.wrappers.iso_wrapper import BushMode, DistortionMode, SubgroupInfo
from features.export.distortion_formats import render_cif
from features.input_cif import build_supercell
from features.method2.distortion_mapper import (
    DistortionMapper,
    _has_exact_unit_phase,
    _make_exact_translation_lattice,
    _snap_to_exact_translation,
)

WORKSPACE = Path(__file__).resolve().parents[2]
QUERY_DIGEST = hashlib.sha256(
    b"synthetic exact secondary-irrep mapping query"
).hexdigest()


def test_parent_translation_snap_uses_metric_threshold_before_exact_phase():
    lattice = _make_exact_translation_lattice(
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        np.diag([4.0, 5.0, 6.0]),
    )

    assert _snap_to_exact_translation(
        np.asarray([1.0 + 2e-5, 0.0, 0.0]), lattice, 1e-4,
    ) == (Fraction(1), Fraction(0), Fraction(0))
    assert _snap_to_exact_translation(
        np.asarray([1.0 + 3e-5, 0.0, 0.0]), lattice, 1e-4,
    ) is None
    assert _snap_to_exact_translation(
        np.asarray([0.25, 0.0, 0.0]), lattice, 1e-4,
    ) is None


def test_parent_translation_snap_is_gl_invariant_and_exact_for_i_centering():
    metric = np.diag([4.0, 4.0, 6.0])
    conventional = _make_exact_translation_lattice(
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)), metric,
    )
    sheared = _make_exact_translation_lattice(
        ((1, 0, 0), (10, 1, 0), (0, 0, 1)), metric,
    )
    delta = np.asarray([1.0 + 1e-6, 1.0 - 1e-6, 0.0])
    expected = (Fraction(1), Fraction(1), Fraction(0))
    assert _snap_to_exact_translation(delta, conventional, 1e-4) == expected
    assert _snap_to_exact_translation(delta, sheared, 1e-4) == expected

    body_centered = _make_exact_translation_lattice(
        centering_primitive_matrix("I"), metric,
    )
    centered = _snap_to_exact_translation(
        np.asarray([0.5 + 1e-6, 0.5 - 1e-6, 0.5]),
        body_centered,
        1e-4,
    )
    assert centered == (Fraction(1, 2),) * 3
    assert _has_exact_unit_phase(((Fraction(1), Fraction(1), Fraction(0)),), centered)
    assert not _has_exact_unit_phase(
        ((Fraction(1), Fraction(0), Fraction(0)),), centered,
    )


def _partial_secondary_mapping_fixture(
    *,
    k_coordinates: tuple[str, str, str] = ("0", "0", "0"),
) -> tuple[Structure, list[dict], SubgroupInfo, DistortionMode, list[list[float]]]:
    """Return one sparse secondary column bound to a complete query context."""

    parent = Structure(
        Lattice(
            [
                [3.0, 0.2, 0.1],
                [0.1, 4.0, 0.3],
                [0.2, 0.4, 5.0],
            ]
        ),
        ["H", "He"],
        [[0.0, 0.0, 0.0], [0.23, 0.31, 0.41]],
    )
    orbit_id = "hydrogen-a"
    sites = [{
        "wyckoff_letter": "a",
        "orbit_id": orbit_id,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    basis = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=1,
        subgroup_index=2,
        size=2,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label="X1-",
        k_point_label="X",
        k_coordinates=["1/2", "0", "0"],
        parent_sg=1,
        basis_vectors=basis,
        basis_raw="(2,0,0),(0,1,0),(0,0,1)",
        origin=[0.0, 0.0, 0.0],
        origin_raw="(0,0,0)",
    )
    block = MicroscopicVectorBlock(
        global_irrep="GM1+",
        direction_symbol="VECTOR,A",
        wyckoff_letter="a",
        site_irrep="A1",
        source_order=0,
        rows=(MicroscopicVectorRow(
            point_raw=("0", "0", "0"),
            displacements=((Fraction(1), Fraction(0), Fraction(0)),),
        ),),
    )
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=QUERY_DIGEST,
        query_order=0,
        query_irrep_label="GM1+",
        source_parent_sg=1,
        source_k_coordinates=k_coordinates,
        source_subgroup_space_group_number=1,
        source_subgroup_basis=basis,
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="X1-",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block,
        column_index=0,
    )
    bush = BushMode(
        irrep_label="GM1+",
        opd_symbol="P1",
        wyckoff_letter="a",
        point=[0.0, 0.0, 0.0],
        point_raw=["0", "0", "0"],
        displacements=[[1.0, 0.0, 0.0]],
    )
    identity = ModeIdentity(
        parent_sg=1,
        global_irrep="GM1+",
        k_coordinates=k_coordinates,
        wyckoff_letter="a",
        orbit_id=orbit_id,
        site_irrep="A1",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="partial",
        reason="awaiting_bush_or_smodes_subspace_validation",
    )
    mode = DistortionMode(
        irrep_label="GM1+",
        dimension=1,
        mode_type="displacive",
        basis_vectors=[[1.0, 0.0, 0.0]],
        wyckoff_site="a",
        wyckoff_orbit_id=orbit_id,
        k_point_label="GM",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        bush_modes=[bush],
        amplitude_key="GM1+[0,0,0]__a@hydrogen-a__A1(a)",
        site_irrep="A1",
        k_coords_label=",".join(k_coordinates),
        opd_component="a",
        mode_identity=identity,
        microscopic_provenance=provenance,
    )
    return parent, sites, subgroup, mode, basis


def _map_partial_fixture(
    parent: Structure,
    sites: list[dict],
    subgroup: SubgroupInfo,
    mode: DistortionMode,
    basis: list[list[float]],
) -> np.ndarray:
    return DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        [mode],
        basis,
        subgroup_translation_lattice=basis,
        subgroup_context=subgroup,
    )[mode.amplitude_key]


def test_exact_partial_secondary_column_can_cover_unit_phase_parent_translations():
    parent, sites, subgroup, mode, basis = _partial_secondary_mapping_fixture()

    mapped = _map_partial_fixture(parent, sites, subgroup, mode, basis)

    child_species = [
        site.species_string for site in build_supercell(parent, basis)
    ]
    hydrogen_rows = mapped[
        np.asarray([species == "H" for species in child_species], dtype=bool)
    ]
    assert hydrogen_rows.shape == (2, 3)
    assert np.allclose(hydrogen_rows, [[0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
    assert mode.mode_identity.status == "partial"


@pytest.mark.parametrize(
    "damage",
    ["missing-provenance", "identity-k", "context-q", "raw-row"],
)
def test_partial_mapping_rejects_missing_or_tampered_source_evidence(damage):
    parent, sites, subgroup, mode, basis = _partial_secondary_mapping_fixture()
    if damage == "missing-provenance":
        mode = replace(mode, microscopic_provenance=None)
    elif damage == "identity-k":
        mode = replace(
            mode,
            mode_identity=replace(
                mode.mode_identity,
                k_coordinates=("1/2", "0", "0"),
            ),
            k_coords_label="1/2,0,0",
        )
    elif damage == "context-q":
        subgroup = replace(
            subgroup,
            origin=[0.5, 0.0, 0.0],
            origin_raw="(1/2,0,0)",
        )
    else:
        mode = replace(
            mode,
            bush_modes=[
                replace(mode.bush_modes[0], point_raw=["1", "0", "0"])
            ],
        )

    with pytest.raises(ValueError, match="partial microscopic mode"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)


def test_partial_mapping_does_not_repeat_a_nonunit_k_phase():
    parent, sites, subgroup, mode, basis = _partial_secondary_mapping_fixture(
        k_coordinates=("1/2", "0", "0"),
    )

    with pytest.raises(ValueError, match="does not cover child atom"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)


def _active_arm_mapping_fixture(
    *,
    irrep: str,
    k_coordinates: tuple[str, str, str],
    direction: str,
    active_k_coordinates: tuple[tuple[str, str, str], ...],
    basis: list[list[float]],
    points: tuple[tuple[str, str, str], ...],
    arrows: tuple[tuple[int, int, int], ...],
) -> tuple[Structure, list[dict], SubgroupInfo, DistortionMode, list[list[float]]]:
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    orbit_id = "hydrogen-a"
    sites = [{
        "wyckoff_letter": "a",
        "orbit_id": orbit_id,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=1,
        subgroup_index=round(abs(np.linalg.det(np.asarray(basis)))),
        size=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label=irrep,
        k_point_label=irrep.rstrip("0123456789+-"),
        k_coordinates=list(k_coordinates),
        parent_sg=123,
        basis_vectors=basis,
        basis_raw=",".join(
            "(" + ",".join(str(int(value)) for value in row) + ")"
            for row in basis
        ),
        origin=[0.0, 0.0, 0.0],
        origin_raw="(0,0,0)",
    )
    rows = tuple(
        MicroscopicVectorRow(
            point_raw=point,
            displacements=(tuple(Fraction(value) for value in arrow),),
        )
        for point, arrow in zip(points, arrows, strict=True)
    )
    block = MicroscopicVectorBlock(
        global_irrep=irrep,
        direction_symbol=direction,
        wyckoff_letter="a",
        site_irrep="A1",
        source_order=0,
        rows=rows,
    )
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=QUERY_DIGEST,
        query_order=0,
        query_irrep_label=irrep,
        source_parent_sg=123,
        source_k_coordinates=k_coordinates,
        source_active_k_coordinates=active_k_coordinates,
        source_subgroup_space_group_number=1,
        source_subgroup_basis=basis,
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label=irrep,
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block,
        column_index=0,
    )
    bushes = [
        BushMode(
            irrep_label=irrep,
            opd_symbol="P1",
            wyckoff_letter="a",
            point=[float(Fraction(value)) for value in point],
            point_raw=list(point),
            displacements=[[float(value) for value in arrow]],
        )
        for point, arrow in zip(points, arrows, strict=True)
    ]
    identity = ModeIdentity(
        parent_sg=123,
        global_irrep=irrep,
        k_coordinates=k_coordinates,
        wyckoff_letter="a",
        orbit_id=orbit_id,
        site_irrep="A1",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="partial",
        reason="awaiting_bush_or_smodes_subspace_validation",
    )
    mode = DistortionMode(
        irrep_label=irrep,
        dimension=1,
        mode_type="displacive",
        basis_vectors=[list(arrows[0])],
        wyckoff_site="a",
        wyckoff_orbit_id=orbit_id,
        k_point_label=subgroup.k_point_label,
        opd_symbol="P1",
        opd_dir_raw="(" + direction.split(",", 1)[1] + ")",
        bush_modes=bushes,
        amplitude_key=f"{irrep}__a__A1(a)",
        site_irrep="A1",
        k_coords_label=",".join(k_coordinates),
        opd_component="a",
        mode_identity=identity,
        microscopic_provenance=provenance,
    )
    return parent, sites, subgroup, mode, basis


def test_x_active_arm_provenance_allows_only_inactive_axis_repetition():
    fixture = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )

    mapped = _map_partial_fixture(*fixture)

    assert np.allclose(
        sorted(mapped[:, 0]),
        [-1 / 3, -1 / 3, -1 / 3, 1 / 3, 1 / 3, 1 / 3],
    )


def test_dt_active_conjugate_pair_allows_transverse_repetition():
    fixture = _active_arm_mapping_fixture(
        irrep="DT1",
        k_coordinates=("0", "1/3", "0"),
        direction="VECTOR,0,0,A,0",
        active_k_coordinates=(("1/3", "0", "0"), ("2/3", "0", "0")),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("1", "0", "0"), ("2", "0", "0")),
        arrows=((0, 0, 0), (1, 0, 0), (-1, 0, 0)),
    )

    mapped = _map_partial_fixture(*fixture)

    assert np.allclose(
        sorted(mapped[:, 0]),
        [-1 / 3, -1 / 3, 0.0, 0.0, 1 / 3, 1 / 3],
    )


def test_sparse_x_column_without_active_arm_provenance_stays_fail_closed():
    parent, sites, subgroup, mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_coordinates=(),
            source_active_k_resolution_kind=None,
        ),
    )

    with pytest.raises(ValueError, match="does not cover child atom"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)


def test_tampered_active_arm_provenance_is_rejected_before_repetition():
    parent, sites, subgroup, mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_coordinates=(("1/2", "0", "0"),),
        ),
    )

    with pytest.raises(ValueError, match="disagrees with its exact VECTOR"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)


def test_tampered_active_arm_resolution_kind_is_rejected_before_repetition():
    parent, sites, subgroup, mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_resolution_kind="compatible_envelope",
        ),
    )

    with pytest.raises(ValueError, match="disagrees with its exact VECTOR"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)


def test_method1_candidate_method2_flow_forwards_active_arm_context():
    """The selected Method-1 candidate remains the mapping proof context."""

    parent, sites, subgroup, mode, _basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mapper = DistortionMapper()
    mapper.map_modes_to_atoms = lambda *_args, **_kwargs: {}

    api = object.__new__(IsoDistort)
    api.structure = parent
    api.symmetry_info = {
        "space_group_number": 123,
        "space_group_symbol": "P4/mmm",
        "wyckoff_sites": sites,
    }
    api.cfg = get_config()
    api.subgroups = []
    api._unresolved_method3_embedding_keys = set()
    api.distortion_types = ["displacive"]
    api.distortion_scope = {"displacive": ["*"]}
    api._selected_subgroup = None
    api._dist_mapper = mapper
    api._iso = SimpleNamespace(list_k_points=lambda *_args, **_kwargs: [])
    api._smodes = None
    api._search = SimpleNamespace(
        method_2_search=lambda *_args, **_kwargs: SimpleNamespace(
            subgroup=subgroup,
            modes=[mode],
            metadata={},
        )
    )
    api._compute_scoped_modes = lambda *_args, **_kwargs: [mode]
    api._resolve_k_vector = lambda *_args, **_kwargs: [0.0, 0.5, 0.0]
    api._mode_context_key = lambda *_args, **_kwargs: ("test-context",)
    api.mode_occupancies = {}
    api.mode_displacements_sc = {}
    api._mode_label_overrides = {}

    result = api.search_method_2(
        subgroup.index,
        distortion_type=["displacive"],
        candidates=[subgroup],
    )

    assert result.subgroup is subgroup
    assert api._selected_subgroup is subgroup
    assert set(api.mode_displacements_sc) == {mode.amplitude_key}
    assert np.allclose(
        sorted(api.mode_displacements_sc[mode.amplitude_key][:, 0]),
        [-1 / 3, -1 / 3, -1 / 3, 1 / 3, 1 / 3, 1 / 3],
    )


def test_microscopic_column_mapping_forwards_and_revalidates_active_arm_context():
    """The high-level mapping adapter must not drop or trust active-arm evidence."""

    parent, sites, subgroup, partial_mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        partial_mode,
        mode_identity=replace(
            partial_mode.mode_identity,
            status="verified",
            reason=None,
        ),
    )
    child = build_supercell(parent, basis)
    kwargs = {
        "subgroup_context": subgroup,
        "frame_id": "active-arm-child-frame",
        "atom_ids": tuple(f"atom-{index}" for index in range(len(child))),
        "subgroup_translation_lattice": basis,
    }
    mapper = DistortionMapper()

    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        sites,
        [mode],
        basis,
        **kwargs,
    )

    assert len(mapped) == 1
    assert mapped[0].mode_id == mode.mode_identity.stable_token
    assert mapped[0].provenance.source_active_k_coordinates == (
        ("0", "1/2", "0"),
    )

    tampered = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_coordinates=(("1/2", "0", "0"),),
        ),
    )
    with pytest.raises(ValueError, match="disagrees with its exact VECTOR"):
        mapper.map_microscopic_columns_to_supercell(
            parent,
            sites,
            [tampered],
            basis,
            **kwargs,
        )


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
def test_live_4310_secondary_e_modes_export_for_p_and_x_primary_routes():
    parent_path = WORKSPACE / "experiment_data" / "4310_tetra.cif"
    if not parent_path.is_file():
        pytest.skip("4310 parent input unavailable")

    api = IsoDistort(language="en")
    api.set_distortion_scope({
        "displacive": ["*"],
        "occupational": [],
        "strain": ["*"],
        "magnetic": [],
        "rotational": [],
    })
    with contextlib.redirect_stdout(io.StringIO()):
        api.load_structure(parent_path)
        rows = api.search_method_1(distortion_types=["strain", "displacive"])
    candidates = [row.subgroup for row in rows]
    target_specs = (
        ("P2", "C1", 120, "(1/2,1/2,0)", 18),
        ("X1-", "C1", 50, "(0,1/2,0)", 16),
    )
    for irrep, opd, space_group, origin_raw, expected_count in target_specs:
        target = next(
            subgroup
            for subgroup in candidates
            if (
                subgroup.irrep_label == irrep
                and subgroup.opd_symbol == opd
                and subgroup.space_group_number == space_group
                and subgroup.origin_raw == origin_raw
            )
        )
        with contextlib.redirect_stdout(io.StringIO()):
            api.search_method_2(
                subgroup_idx=target.index,
                distortion_type=["displacive"],
                number_of_independent_modulations=0,
                candidates=candidates,
            )
        assert len(api.distortion_modes) == expected_count
        assert all(
            mode.mode_identity is not None
            and mode.mode_identity.status == "verified"
            and mode.mode_identity.source == "iso_microscopic"
            for mode in api.distortion_modes
        )
        spec = api._spec_for_subgroup(
            target,
            use_current_modes=True,
            use_generated_structure=False,
            folder_name=f"test_{irrep}_{opd}",
        )
        cif_text = render_cif(spec.structure, spec)
        assert "_iso_displacivemode_ID" in cif_text
