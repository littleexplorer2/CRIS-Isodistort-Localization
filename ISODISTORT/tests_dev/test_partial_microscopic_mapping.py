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

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from isocore.api import IsoDistort
from isocore.backend.iso_mode_models import (
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
)
from isocore.backend.iso_wrapper import BushMode, DistortionMode, SubgroupInfo
from isocore.distortion.distortion_mapper import DistortionMapper
from isocore.io.distortion_formats import render_cif
from isocore.structure import build_supercell

WORKSPACE = Path(__file__).resolve().parents[2]
QUERY_DIGEST = hashlib.sha256(
    b"synthetic exact secondary-irrep mapping query"
).hexdigest()


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

