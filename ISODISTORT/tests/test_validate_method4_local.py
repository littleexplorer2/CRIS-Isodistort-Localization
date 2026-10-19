"""Tests for test_validate_method4_local; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.models.iso_mode_models import MicroscopicColumnProvenance
from backend.models.iso_mode_models import MicroscopicVectorBlock
from backend.models.iso_mode_models import MicroscopicVectorRow
from backend.models.iso_mode_models import ModeIdentity
from backend.wrappers import BushMode
from backend.wrappers import DistortionMode
from fractions import Fraction
from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests.validate_method4_local import _identity_readiness
from tests.validate_method4_local import _mode_source_quantization_bound
from tests.validate_method4_local import _raw_coordinate_roundtrip
from tests.validate_method4_local import _reconstructed_max_displacement
from types import SimpleNamespace
import numpy as np
import pytest


# Original coverage: test_method4_official_audit.py.

def test_source_quantization_bound_uses_actual_mode_matrix_sigma_min() -> None:
    half = Fraction(1, 2000)
    vector_row = MicroscopicVectorRow(
        point_raw=("0", "0", "0"),
        displacements=((Fraction(577, 1000), Fraction(0), Fraction(0)),),
        displacement_half_steps=((half, Fraction(0), Fraction(0)),),
    )
    block = MicroscopicVectorBlock(
        global_irrep="GM1+",
        wyckoff_letter="a",
        site_irrep="A",
        rows=(vector_row,),
        source_order=0,
        direction_symbol="VECTOR,A",
    )
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest="a" * 64,
        query_order=0,
        query_irrep_label="GM1+",
        source_parent_sg=1,
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
    mode = DistortionMode(
        irrep_label="GM1+",
        amplitude_key="mode",
        bush_modes=[
            BushMode(
                irrep_label="GM1+",
                opd_symbol="P1",
                wyckoff_letter="a",
                point=[0.0, 0.0, 0.0],
                displacements=[[0.577, 0.0, 0.0]],
            )
        ],
        microscopic_provenance=provenance,
    )
    structure = Structure(Lattice.cubic(4.0), ["Na"], [[0.0, 0.0, 0.0]])
    iso = SimpleNamespace(structure=structure, distortion_modes=[mode])
    result = SimpleNamespace(
        metadata={"primitive_cell_multiplicity": 1, "supercell_size": 1.0},
        rms_residual=0.0,
        amplitudes={"mode": 0.2},
    )

    evidence = _mode_source_quantization_bound(
        iso,
        result,
        {"mode": np.array([[0.5, 0.0, 0.0]])},
        structure.lattice.matrix,
    )

    assert evidence["status"] == "available"
    assert evidence["unit_mode_matrix_smallest_singular_value"] == pytest.approx(1.0)
    assert evidence["unit_mode_matrix_l2_error_bound"] > 0.0
    assert evidence["As_vector_l2_error_bound_angstrom"] > 0.0
    assert (
        evidence["mode_records"]["mode"]
        ["normfactor_error_bound_inverse_angstrom"]
        > 0.0
    )

def _validation_mode(label: str, identity: ModeIdentity) -> DistortionMode:
    return DistortionMode(
        irrep_label=identity.global_irrep,
        amplitude_key=label,
        wyckoff_site=identity.wyckoff_letter,
        site_irrep=identity.display_site_irrep or "",
        opd_component=identity.component_label or "a",
        mode_identity=identity,
    )

def test_local_identity_readiness_requires_verified_unique_mode_identities() -> None:
    verified = ModeIdentity(
        parent_sg=123,
        global_irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        wyckoff_letter="a",
        orbit_id="orbit-a",
        site_irrep="Eu",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="verified",
    )
    label = "X1+[0,1/2,0]__a__Eu(a)"
    ready = _identity_readiness(
        SimpleNamespace(distortion_modes=[_validation_mode(label, verified)]),
        [label],
    )
    assert ready["status"] == "pass"

    unresolved = ModeIdentity.unresolved(
        parent_sg=123,
        global_irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        wyckoff_letter="a",
        orbit_id="orbit-a",
        reason="missing exact microscopic source",
    )
    blocked = _identity_readiness(
        SimpleNamespace(distortion_modes=[_validation_mode(label, unresolved)]),
        [label],
    )
    assert blocked["status"] == "inconclusive"
    assert "missing exact microscopic source" in blocked["unresolved"][0]

def test_local_raw_coordinates_are_not_compared_without_frozen_basis() -> None:
    row = {
        "expect": "success",
        "contributions": {"old-label": 0.25},
    }
    result = SimpleNamespace(raw_coefficients={"new-label": 1.5})
    tolerances = {"exact_amplitude_abs": 1e-8, "noise_amplitude_abs": 1e-4}
    current = {"schema": 1, "sha256": "a" * 64}

    comparison = _raw_coordinate_roundtrip(
        row,
        result,
        tolerances,
        frozen_contract=None,
        current_contract=current,
    )

    assert comparison["status"] == "not_comparable"
    assert comparison["reason"] == "generation_basis_fingerprint_absent"

def test_local_raw_coordinates_fail_only_in_identical_basis() -> None:
    row = {"expect": "success", "contributions": {"mode": 0.25}}
    result = SimpleNamespace(raw_coefficients={"mode": 1.5})
    tolerances = {"exact_amplitude_abs": 1e-8, "noise_amplitude_abs": 1e-4}
    contract = {"schema": 1, "sha256": "b" * 64}

    comparison = _raw_coordinate_roundtrip(
        row,
        result,
        tolerances,
        frozen_contract=contract,
        current_contract=contract,
    )

    assert comparison["status"] == "fail"
    assert comparison["reason"] == "identical_basis_fingerprint"

def test_reconstructed_dmax_uses_current_physical_mode_vectors() -> None:
    modes = {
        "a": np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        "b": np.array([[0.0, 0.0, 0.0], [0.0, -2.0, 0.0]]),
    }
    result = SimpleNamespace(raw_coefficients={"a": 0.3, "b": -0.1})

    assert np.isclose(
        _reconstructed_max_displacement(result, modes, np.eye(3)),
        0.3,
    )
