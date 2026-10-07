import json
import sys
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from backend.models.iso_mode_models import (
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
)
from backend.wrappers import BushMode, DistortionMode
from tests.manual import audit_method4_official as official_audit
from tests.manual import prepare_method4_inputs as preparation
from tests.manual.audit_method4_official import (
    _amplitude_comparison_tolerance,
    _complete_gamma_site_subspaces,
    _complete_origin_phase_amplitude_subspaces,
    _complete_required_amplitude_subspaces,
    _expected_failure_matches,
    _isoviz_applied_strain,
    _isoviz_applied_strain_with_rounding,
    _local_mode_identity,
    _official_mode_identity,
    _parameter_mode_subspace_key,
    _parameter_origin_phase_turns,
    _subspace_amplitude_tolerances,
)
from tests.manual.method4_provenance import (
    PROVENANCE_KIND,
    PROVENANCE_SCHEMA,
    build_method4_provenance_state,
)
from tests.manual.validate_method4_local import (
    _identity_readiness,
    _mode_source_quantization_bound,
    _raw_coordinate_roundtrip,
    _reconstructed_max_displacement,
)


def test_isoviz_strain_parser_applies_archived_mode_vectors() -> None:
    text = """
!strainmodelist
  1    0.06016   0.10000    1 GM1+strain_1(a)
   0.70711   0.70711   0.00000   0.00000   0.00000   0.00000
  2   -0.06016   0.10000    2 GM2+strain(a)
   0.70711  -0.70711   0.00000   0.00000   0.00000   0.00000
  3    0.09683   0.10000    4 GM5+strain(a)
   0.00000   0.00000   0.00000   1.41421   0.00000   0.00000
!displacivemodelist
"""

    np.testing.assert_allclose(
        _isoviz_applied_strain(text),
        [0.0, 0.06016 * 2 * 0.70711, 0.0, 0.09683 * 1.41421, 0.0, 0.0],
        atol=1e-12,
    )
    applied, bound = _isoviz_applied_strain_with_rounding(text)
    exact = np.array(
        [0.0, 0.06016 * np.sqrt(2), 0.0, 0.09683 * np.sqrt(2), 0.0, 0.0]
    )
    assert np.all(np.abs(applied - exact) <= bound + 1e-12)


def test_expected_failure_matching_rejects_generic_server_errors() -> None:
    assert _expected_failure_matches(
        "reject_species_or_stoichiometry",
        "Types of atoms in subgroup do not match types of atoms in parent.",
    )
    assert _expected_failure_matches(
        "reject_unmatched_atom",
        "Failed to find match. Try using a larger value for dmax.",
    )
    assert not _expected_failure_matches(
        "reject_incompatible_lattice",
        "Internal server error. Please try again later.",
    )


def test_gamma_mode_identity_uses_stable_component() -> None:
    assert _local_mode_identity("GM5+__d__P1(1)__a") == _official_mode_identity(
        "[0,0,0]GM5+(a,b)[Al1:d:dsp]E(a)"
    )


def test_complete_gamma_site_subspace_accepts_repeated_official_components() -> None:
    local_rows = [
        {
            "local_label": f"GM5-__f__C1(1)__{component}",
            "identity": list(_local_mode_identity(
                f"GM5-__f__C1(1)__{component}"
            )),
            "official_label": None,
        }
        for component in "abcd"
    ]
    official_labels = [
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(a)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(b)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B2u(a)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B2u(b)",
    ]
    official_rows = [
        {
            "label": label,
            "identity": list(_official_mode_identity(label)),
        }
        for label in official_labels
    ]

    groups = _complete_gamma_site_subspaces(local_rows, official_rows)

    assert len(groups) == 1
    key, grouped_local, grouped_official = groups[0]
    assert key == ("GM5-", "f")
    assert grouped_local == local_rows
    assert grouped_official == official_rows


def test_incomplete_gamma_site_subspace_is_not_compared_by_norm() -> None:
    local_rows = [
        {
            "local_label": f"GM5-__f__C1(1)__{component}",
            "identity": list(_local_mode_identity(
                f"GM5-__f__C1(1)__{component}"
            )),
            "official_label": None,
        }
        for component in "abcd"
    ]
    official_labels = [
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(a)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(b)",
    ]
    official_rows = [
        {
            "label": label,
            "identity": list(_official_mode_identity(label)),
        }
        for label in official_labels
    ]

    assert _complete_gamma_site_subspaces(local_rows, official_rows) == []


def test_exact_amplitude_tolerance_does_not_absorb_structure_mismatch() -> None:
    tolerance, evidence = _amplitude_comparison_tolerance(
        {},
        {"origin_aligned_root_sum_squared_distance_angstrom": 0.25},
    )

    assert np.isclose(tolerance, 5.1e-6)
    assert evidence["aligned_structure_rss_bound_angstrom"] == 0.0


def _noisy_matrix_evidence(
    sigma_min: float = 1.0,
    matrix_error: float = 0.0,
    supercell_size: float = 4.0,
    source_as: float = 0.0,
) -> dict:
    return {
        "metadata": {"supercell_size": supercell_size},
        "mode_source_quantization_bound": {
            "status": "available",
            "unit_mode_matrix_smallest_singular_value": sigma_min,
            "unit_mode_matrix_l2_error_bound": matrix_error,
            "As_vector_l2_error_bound_angstrom": source_as,
            "Ap_vector_l2_error_bound_angstrom": source_as / 2.0,
        },
    }


def test_noisy_amplitude_tolerance_uses_matrix_sensitivity_and_Ap_scale() -> None:
    tolerance, evidence = _amplitude_comparison_tolerance(
        {"noise_sigma_fractional": 2e-5},
        {"origin_aligned_root_sum_squared_distance_angstrom": 4.2e-5},
        _noisy_matrix_evidence(),
    )

    assert np.isclose(tolerance, 4.71e-5)
    assert np.isclose(evidence["aligned_structure_rss_bound_angstrom"], 4.2e-5)
    assert evidence["As_total_angstrom"] == pytest.approx(4.71e-5)
    assert evidence["Ap_total_angstrom"] == pytest.approx(2.61e-5)


def test_noisy_nonorthogonal_modes_need_more_than_the_structure_rss() -> None:
    # Both columns have unit norm, yet their coefficient map amplifies a y
    # perturbation by about 100. Cauchy--Schwarz on the columns is insufficient.
    unit_matrix = np.column_stack((
        np.array([1.0, 0.0]),
        np.array([1.0, 0.01]) / np.sqrt(1.0001),
    ))
    delta = np.array([0.0, 1e-5])
    coefficient_change = np.linalg.lstsq(unit_matrix, delta, rcond=None)[0]
    rss = float(np.linalg.norm(delta))
    sigma_min = float(np.linalg.svd(unit_matrix, compute_uv=False)[-1])
    _, evidence = _amplitude_comparison_tolerance(
        {"noise_sigma_fractional": 2e-5},
        {"origin_aligned_root_sum_squared_distance_angstrom": rss},
        _noisy_matrix_evidence(sigma_min=sigma_min),
    )

    as_bound = evidence["aligned_structure_As_vector_l2_bound_angstrom"]
    assert np.linalg.norm(coefficient_change) > 100.0 * rss
    assert as_bound == pytest.approx(rss / sigma_min)
    assert np.linalg.norm(coefficient_change) <= as_bound
    assert evidence["aligned_structure_Ap_vector_l2_bound_angstrom"] == pytest.approx(
        as_bound / 2.0
    )


def test_noisy_structure_and_source_bounds_share_weyl_denominator_once() -> None:
    _, evidence = _amplitude_comparison_tolerance(
        {"noise_sigma_fractional": 2e-5},
        {"origin_aligned_root_sum_squared_distance_angstrom": 4e-5},
        _noisy_matrix_evidence(sigma_min=0.5, matrix_error=0.1, source_as=2e-4),
    )

    assert evidence["amplitude_sensitivity_denominator"] == pytest.approx(0.4)
    assert evidence["As_total_angstrom"] == pytest.approx(5.1e-6 + 2e-4 + 4e-5 / 0.4)
    assert evidence["Ap_total_angstrom"] == pytest.approx(5.1e-6 + (2e-4 + 4e-5 / 0.4) / 2.0)
    assert _subspace_amplitude_tolerances(evidence, 4) == pytest.approx(
        (1.02e-5 + 3e-4, 1.02e-5 + 1.5e-4)
    )


@pytest.mark.parametrize("invalid_evidence", [
    None,
    {"mode_source_quantization_bound": {"status": "unavailable"}},
    _noisy_matrix_evidence(sigma_min=0.1, matrix_error=0.1),
    _noisy_matrix_evidence(sigma_min=0.1, matrix_error=0.2),
    _noisy_matrix_evidence(sigma_min=float("nan")),
    _noisy_matrix_evidence(matrix_error=-0.1),
    _noisy_matrix_evidence(supercell_size=0.0),
    {"mode_source_quantization_bound": _noisy_matrix_evidence()["mode_source_quantization_bound"]},
    {"metadata": {"supercell_size": 4.0}, "mode_source_quantization_bound": {"status": "available"}},
])
def test_noisy_amplitude_tolerance_without_valid_matrix_evidence_is_inconclusive(
    invalid_evidence,
) -> None:
    tolerance, evidence = _amplitude_comparison_tolerance(
        {"noise_sigma_fractional": 2e-5},
        {"origin_aligned_root_sum_squared_distance_angstrom": 4e-5},
        invalid_evidence,
    )

    assert tolerance is None
    assert evidence["status"] == "inconclusive"
    assert evidence["reason"]
    assert evidence["As_total_angstrom"] is None
    assert evidence["Ap_total_angstrom"] is None
    with pytest.raises(ValueError, match="unavailable"):
        _subspace_amplitude_tolerances(evidence, 1)


def test_noisy_amplitude_tolerance_requires_structure_difference_evidence() -> None:
    tolerance, evidence = _amplitude_comparison_tolerance(
        {"noise_sigma_fractional": 2e-5}, {}, _noisy_matrix_evidence(),
    )

    assert tolerance is None
    assert evidence["status"] == "inconclusive"
    assert "structure RSS evidence" in evidence["reason"]


def test_amplitude_tolerance_adds_propagated_source_quantization_bound() -> None:
    tolerance, evidence = _amplitude_comparison_tolerance(
        {},
        {},
        {
            "mode_source_quantization_bound": {
                "status": "available",
                "As_vector_l2_error_bound_angstrom": 2.0e-4,
                "Ap_vector_l2_error_bound_angstrom": 1.0e-4,
            }
        },
    )

    assert tolerance == pytest.approx(2.051e-4)
    assert _subspace_amplitude_tolerances(evidence, 4) == pytest.approx(
        (2.102e-4, 1.102e-4)
    )


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


def test_parameter_mode_identity_retains_wavevector_and_site_irrep() -> None:
    assert _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_2(b)"
    ) == _official_mode_identity(
        "[0,0,1/6]LD1(a,b)[Al2:e:dsp]A1_2(b)"
    )
    assert _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_2(b)"
    ) != _local_mode_identity(
        "LD1[0,0,1/3]__e__A1_2(b)"
    )
    assert _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_1(b)"
    ) != _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_2(b)"
    )


def test_complete_parameter_subspace_uses_basis_independent_group() -> None:
    identity_a = _local_mode_identity("LD1[0,0,1/6]__a__A2u(a)")
    identity_b = _local_mode_identity("LD1[0,0,1/6]__a__A2u(b)")
    assert _parameter_mode_subspace_key(identity_a) == (
        "LD1[0,0,1/6]", "a", "A2u",
    )
    rows = [
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(a)",
            "identity": list(identity_a),
            "official_label": "official-a",
        },
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(b)",
            "identity": list(identity_b),
            "official_label": "official-b",
        },
    ]
    required = {row["local_label"] for row in rows}
    groups = _complete_required_amplitude_subspaces(
        rows, required, {identity_a, identity_b},
    )
    assert groups == [rows]


def test_partial_parameter_subspace_keeps_componentwise_comparison() -> None:
    identity_a = _local_mode_identity("LD1[0,0,1/6]__a__A2u(a)")
    identity_b = _local_mode_identity("LD1[0,0,1/6]__a__A2u(b)")
    rows = [
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(a)",
            "identity": list(identity_a),
            "official_label": "official-a",
        },
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(b)",
            "identity": list(identity_b),
            "official_label": "official-b",
        },
    ]
    groups = _complete_required_amplitude_subspaces(
        rows,
        {"LD1[0,0,1/6]__a__A2u(a)"},
        {identity_a, identity_b},
    )
    assert groups == []


def test_parameter_origin_phase_is_exact_and_periodic() -> None:
    key = ("Y1[1/3,1/2,0]", "d", "Eu")

    assert _parameter_origin_phase_turns(key, [1.0, 0.0, 0.0]) == Fraction(1, 3)
    assert _parameter_origin_phase_turns(key, [0.0, -1.0, 0.0]) == Fraction(1, 2)
    assert _parameter_origin_phase_turns(key, [3.0, 0.0, 0.0]) == 0


def test_origin_phase_subspace_accepts_full_local_branch_in_official_expansion() -> None:
    identity_a = _local_mode_identity("Y1[1/3,1/2,0]__d__Eu(a)")
    identity_b = _official_mode_identity(
        "[1/3,1/2,0]Y1(a,b;0,0)[ND:d:dsp]Eu(b)"
    )
    local_rows = [{
        "local_label": "Y1[1/3,1/2,0]__d__Eu(a)",
        "identity": list(identity_a),
        "official_label": "official-a",
    }]
    official_rows = [
        {"label": "official-a", "identity": list(identity_a)},
        {"label": "official-b", "identity": list(identity_b)},
    ]

    groups = _complete_origin_phase_amplitude_subspaces(
        local_rows,
        {local_rows[0]["local_label"]},
        official_rows,
        [1.0, 0.0, 0.0],
    )

    assert len(groups) == 1
    key, grouped_local, grouped_official, phase_turns = groups[0]
    assert key == ("Y1[1/3,1/2,0]", "d", "Eu")
    assert grouped_local == local_rows
    assert grouped_official == official_rows
    assert phase_turns == Fraction(1, 3)


def test_origin_phase_subspace_rejects_trivial_phase_or_partial_local_branch() -> None:
    identity_a = _local_mode_identity("Y1[1/3,1/2,0]__d__Eu(a)")
    local_rows = [{
        "local_label": "Y1[1/3,1/2,0]__d__Eu(a)",
        "identity": list(identity_a),
        "official_label": "official-a",
    }]
    official_rows = [{"label": "official-a", "identity": list(identity_a)}]

    assert _complete_origin_phase_amplitude_subspaces(
        local_rows,
        {local_rows[0]["local_label"]},
        official_rows,
        [3.0, 0.0, 0.0],
    ) == []
    assert _complete_origin_phase_amplitude_subspaces(
        local_rows,
        set(),
        official_rows,
        [1.0, 0.0, 0.0],
    ) == []


def test_dmax_scaling_probes_inside_unwrapped_linear_regime(monkeypatch) -> None:
    lattice = Lattice.cubic(10.0)
    zero = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])

    def wrapped_generator(_iso, contributions):
        coefficient = contributions.get("mode", 0.0)
        return Structure(
            lattice,
            ["H"],
            [[0.0, 0.0, (2.0 * coefficient) % 1.0]],
        )

    monkeypatch.setattr(preparation, "_generated", wrapped_generator)
    contributions, generated = preparation._scaled_contributions_for_dmax(
        object(), zero, {"mode": 1.0}, 0.30,
    )

    assert np.isclose(contributions["mode"], 0.015)
    displacement = np.linalg.norm(generated.cart_coords[0] - zero.cart_coords[0])
    assert np.isclose(displacement, 0.30)


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


def test_method4_provenance_covers_code_config_runtime_backend_and_inputs(
    tmp_path: Path,
) -> None:
    package_root = tmp_path / "ISODISTORT"
    files = {
        "tests/manual/validate_method4_local.py": "validator-v1",
        "tests/manual/prepare_method4_inputs.py": "preparer-v1",
        "backend/api/core_api.py": "core-v1",
        "features/method4/strain_modes.py": "strain-modes-v1",
        "pyproject.toml": "[project]\nname='fixture'\n",
        "requirements.txt": "numpy\nscipy\n",
        "runtime/iso": "iso-v1",
        "runtime/findsym": "findsym-v1",
        "runtime/smodes": "smodes-v1",
        "runtime/const.dat": "const-v1",
        "runtime/data_isotropy.txt": "data-v1",
    }
    for relative, content in files.items():
        path = package_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    config_path = package_root / "resources" / "config" / "settings.yaml"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(
        json.dumps({
            "isobyu": {
                "bin_dir": "../../runtime",
                "data_dir": "../../runtime",
                "iso_bin": "iso",
                "findsym_bin": "findsym",
                "smodes_bin": "smodes",
            }
        }),
        encoding="utf-8",
    )
    input_root = tmp_path / "inputs"
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("parent-v1", encoding="utf-8")
    (input_root / "Daughter.cif").write_text("daughter-v1", encoding="utf-8")
    manifest_path = package_root / "docs" / "manifests" / "method4.json"
    manifest_path.parent.mkdir(parents=True)
    manifest = {
        "signature": "manifest-source-v1",
        "parents": {
            "fixture": {
                "parent_cif": "inputs/Parent.cif",
                "cases": {
                    "C01": {"daughter_cif": "inputs/Daughter.cif"},
                },
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    selection = {"parent": None, "context": None, "case_id": None}

    def signature_state() -> dict:
        return build_method4_provenance_state(
            manifest_path,
            selection,
            workspace=tmp_path,
            package_root=package_root,
        )

    state = signature_state()
    labels = {record["label"] for record in state["inputs"]}
    assert {
        "runtime-config",
        "input-manifest",
        "production:backend/api/core_api.py",
        "production:features/method4/strain_modes.py",
        "backend:iso-binary",
        "backend:findsym-binary",
        "backend:smodes-binary",
        "backend:iso-data:const.dat",
        "backend:iso-data:data_isotropy.txt",
        "parent:fixture",
        "daughter:fixture:C01",
    } <= labels
    assert state["runtime"]["python"]["executable"]["sha256"]
    assert set(state["runtime"]["distributions"]) == {
        "numpy", "scipy", "pymatgen", "spglib", "PyYAML",
    }

    signatures = [state["signature"]]
    mutation_path = package_root / "features" / "method4" / "strain_modes.py"
    mutation_path.write_text("strain-modes-v2", encoding="utf-8")
    signatures.append(signature_state()["signature"])
    config_path.write_text(config_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    signatures.append(signature_state()["signature"])
    manifest["marker"] = "manifest-v2"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    signatures.append(signature_state()["signature"])
    (package_root / "runtime" / "findsym").unlink()
    signatures.append(signature_state()["signature"])

    assert len(set(signatures)) == len(signatures)


def _valid_local_report(
    state: dict,
    manifest: dict,
    manifest_path: Path,
) -> dict:
    selection = {"parent": None, "context": None, "case_id": None}
    return {
        "schema": 3,
        "kind": "method4_local_inverse_consistency",
        "status": "pass",
        "failed": 0,
        "provenance_failure_count": 0,
        "manifest_sha256": official_audit._sha256(manifest_path),
        "manifest_source_signature": manifest["signature"],
        "tolerances": manifest["local_validation_tolerances"],
        "selection": selection,
        "case_count": 1,
        "parents": {"fixture": {"cases": {"C01": {"status": "pass"}}}},
        "provenance": {
            "schema": PROVENANCE_SCHEMA,
            "kind": PROVENANCE_KIND,
            "signature": state["signature"],
            "state": state,
            "run_start_signature": state["signature"],
            "run_end_signature_before_manifest_update": state["signature"],
            "stable_during_run": True,
            "changed_surfaces_during_run": [],
            "manifest_self_update": {
                "requested": False,
                "performed": False,
                "changed_surfaces": [],
                "unexpected_changed_surfaces": [],
            },
        },
    }


def test_official_gate_independently_accepts_only_current_local_provenance(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = {
        "signature": "manifest-source-v1",
        "local_validation_tolerances": {"exact": 1e-8},
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    selection = {"parent": None, "context": None, "case_id": None}
    current_state = {
        "schema": PROVENANCE_SCHEMA,
        "kind": PROVENANCE_KIND,
        "selection": selection,
        "signature": "current-generation-signature",
    }
    calls = []

    def independently_recompute(path, selected):
        calls.append((path, selected))
        return current_state

    monkeypatch.setattr(
        official_audit,
        "build_method4_provenance_state",
        independently_recompute,
    )
    local = _valid_local_report(current_state, manifest, manifest_path)
    accepted = official_audit._local_report_provenance_gate(
        local,
        manifest,
        manifest_path,
    )
    assert accepted["status"] == "pass"

    stale = json.loads(json.dumps(local))
    stale["provenance"]["signature"] = "old-generation-signature"
    stale["provenance"]["state"]["signature"] = "old-generation-signature"
    stale["provenance"]["run_start_signature"] = "old-generation-signature"
    stale["provenance"]["run_end_signature_before_manifest_update"] = (
        "old-generation-signature"
    )
    rejected = official_audit._local_report_provenance_gate(
        stale,
        manifest,
        manifest_path,
    )

    assert rejected["status"] == "fail"
    assert rejected["recorded_signature"] == "old-generation-signature"
    assert rejected["recomputed_signature"] == "current-generation-signature"
    assert any("independent recomputation" in issue for issue in rejected["issues"])
    assert len(calls) == 2


def test_official_main_blocks_legacy_report_before_case_audit(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = {
        "signature": "manifest-source-v1",
        "local_validation_tolerances": {},
        "readable_case_folders": {"C01": "case-folder"},
        "parents": {"fixture": {"cases": {"C01": {}}}},
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    local_path = tmp_path / "legacy-local.json"
    local_path.write_text(json.dumps({"schema": 1}), encoding="utf-8")
    report_path = tmp_path / "official-report.json"
    recomputed = {
        "schema": PROVENANCE_SCHEMA,
        "kind": PROVENANCE_KIND,
        "selection": {"parent": None, "context": None, "case_id": None},
        "signature": "current-generation-signature",
    }
    monkeypatch.setattr(
        official_audit,
        "build_method4_provenance_state",
        lambda _path, _selection: recomputed,
    )

    def forbidden_case_audit(*_args, **_kwargs):
        raise AssertionError("case audit must not run after a provenance gate failure")

    monkeypatch.setattr(official_audit, "_audit_case", forbidden_case_audit)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "audit_method4_official.py",
            "--manifest",
            str(manifest_path),
            "--local-report",
            str(local_path),
            "--json-output",
            str(report_path),
        ],
    )

    assert official_audit.main() == 1
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["summary"]["status"] == "fail"
    assert report["summary"]["audited_case_count"] == 0
    assert report["summary"]["blocked_case_count"] == 1
    assert report["local_provenance_gate"]["status"] == "fail"
    assert any(
        "schema predates" in issue
        for issue in report["local_provenance_gate"]["issues"]
    )


@pytest.mark.parametrize("status,failed", [("inconclusive", 0), ("fail", 2)])
def test_official_gate_accepts_current_provenance_despite_local_overall_status(
    tmp_path: Path, monkeypatch, status: str, failed: int,
) -> None:
    manifest = {"signature": "current-manifest", "local_validation_tolerances": {}}
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    state = {"signature": "current-generation"}
    monkeypatch.setattr(
        official_audit, "build_method4_provenance_state", lambda *_args: state,
    )
    local = _valid_local_report(state, manifest, manifest_path)
    local.update(status=status, failed=failed, inconclusive=1)
    assert official_audit._local_report_provenance_gate(
        local, manifest, manifest_path,
    )["status"] == "pass"


def test_official_gate_rejects_schema2_even_with_current_provenance(
    tmp_path: Path, monkeypatch,
) -> None:
    manifest = {"signature": "current-manifest", "local_validation_tolerances": {}}
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    state = {"signature": "current-generation"}
    monkeypatch.setattr(
        official_audit, "build_method4_provenance_state", lambda *_args: state,
    )
    local = _valid_local_report(state, manifest, manifest_path)
    local["schema"] = 2
    gate = official_audit._local_report_provenance_gate(local, manifest, manifest_path)
    assert gate["status"] == "fail"
    assert len(gate["issues"]) == 1
    assert "schema 3" in gate["issues"][0]


def _html_page(title: str, body: str) -> str:
    return (
        f"<html><head><title>{title}</title></head>"
        f"<body><h1>{title}</h1>{body}</body></html>"
    )


def _method4_result_html(body: str, *, include_form: bool = True) -> str:
    form = (
        '<form method="POST"><input name="input" value="displaydistort"></form>'
        if include_form else ""
    )
    return _html_page("ISODISTORT: distortion", body + form)


def _method4_details_html(body: str) -> str:
    return _html_page(
        "ISODISTORT: complete modes details",
        "<strong>Subgroup details</strong>"
        "<strong>Undistorted superstructure</strong>"
        "<strong>Distorted superstructure</strong>"
        f"<pre>{body}</pre>",
    )


def _method4_basis_html() -> str:
    return _html_page(
        "ISODISTORT: distorted structure (basis)",
        '<form method="POST">'
        '<input name="input" value="distort">'
        '<input name="origintype" value="method4">'
        '<input name="chooseorigin" value="false" checked>'
        '<input name="trynearest" value="true" checked>'
        "</form>",
    )


def _official_success_fixture(tmp_path: Path, monkeypatch):
    """Synthetic archive exercises audit control flow, not scientific equivalence."""
    directory = tmp_path / "archive" / "fixture" / "官网" / "Method4" / "case-folder"
    directory.mkdir(parents=True)
    subgroup = (
        "Subgroup: 1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
        "origin=(0,0,0), s=1, i=1"
    )
    official_label = "[0,0,0]GM1+(a)[X:a:dsp]A1(a)"
    result_html = _method4_result_html(
        subgroup,
    )
    details_html = _method4_details_html(
        f"{subgroup}\nDisplacive mode definitions\n"
        f"a1{official_label} normfactor = 1.0\n"
        f"Displacive mode amplitudes\n{official_label} 0.1 0.1 0.1\n"
        "Parent-cell strain mode definitions\n"
        "Parent-cell strain mode amplitudes\n",
    )
    basis_html = _method4_basis_html()
    files = {
        "arbitrary result capture.HTM": result_html,
        "renamed details.HTML": details_html,
        "manual basis name.hTm": basis_html,
        "subgroup.cif": "_iso_displacivemode_number 1\n_iso_strainmode_number 0\n",
        "topas.str": "prm !a1 0.1\n",
        "data.isoviz": (
            "!strainmodelist\n1 0.00000 0 1 E\n"
            "0.00000 0.00000 0.00000 0.00000 0.00000 0.00000\n"
            "!displacivemodelist\n"
        ),
    }
    for name, content in files.items():
        (directory / name).write_text(content, encoding="utf-8")
    daughter = tmp_path / "daughter.cif"
    daughter.write_text("synthetic frozen input", encoding="utf-8")
    monkeypatch.setattr(official_audit, "OFFICIAL_ROOT", tmp_path / "archive")
    monkeypatch.setattr(official_audit, "WORKSPACE", tmp_path)
    comparison = {
        "expected_atom_count": 1, "actual_atom_count": 1,
        "expected_species": {"X": 1}, "actual_species": {"X": 1},
        "lattice_metric_relative_residual": 0.0,
        "max_species_matched_distance_angstrom": 0.0,
    }
    monkeypatch.setattr(
        official_audit, "_compare_exported_cif", lambda *_args: (comparison, 1e-4),
    )
    case = {"daughter_cif": "daughter.cif", "sha256": official_audit._sha256(daughter),
            "context": "gamma", "expect": "success"}
    parent = {"resolved_contexts": {"gamma": {"basis_vectors": np.eye(3).tolist()}}}
    local = {
        "status": "inconclusive",
        "inverse_consistency": {"status": "pass", "issues": []},
        "identity_readiness": {"status": "inconclusive", "unresolved": ["missing identity"]},
        "strain_voigt_engineering": dict.fromkeys(("xx", "yy", "zz", "2yz", "2xz", "2xy"), 0),
        "mode_normfactors_inverse_angstrom": {"GM1+__a__A1(a)": 999.0},
        "actual_As_angstrom": {"GM1+__a__A1(a)": 999.0},
        "actual_Ap_angstrom": {"GM1+__a__A1(a)": 999.0},
    }
    return case, parent, local, comparison


def test_method4_html_resolution_accepts_arbitrary_stems_and_htm_suffix(
    tmp_path: Path,
    monkeypatch,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)

    result = official_audit._audit_case(
        "fixture", "C01", case, parent, "case-folder", local,
    )

    selected = result["selected_html_pages"]
    assert Path(selected["result_page"]["path"]).name == "arbitrary result capture.HTM"
    assert Path(selected["details_page"]["path"]).name == "renamed details.HTML"
    assert Path(selected["basis_page"]["path"]).name == "manual basis name.hTm"


@pytest.mark.parametrize("role", ["result", "details", "basis"])
def test_method4_html_resolution_rejects_duplicate_roles(
    tmp_path: Path,
    monkeypatch,
    role: str,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)
    directory = tmp_path / "archive" / "fixture" / "官网" / "Method4" / "case-folder"
    source = {
        "result": directory / "arbitrary result capture.HTM",
        "details": directory / "renamed details.HTML",
        "basis": directory / "manual basis name.hTm",
    }[role]
    (directory / f"duplicate-{role}.html").write_text(
        source.read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    result = official_audit._audit_case(
        "fixture", "C01", case, parent, "case-folder", local,
    )

    assert result["status"] == "fail"
    assert any("found 2" in issue for issue in result["issues"])


@pytest.mark.parametrize("missing_role", ["result", "details"])
def test_method4_html_resolution_rejects_missing_or_link_only_pages(
    tmp_path: Path,
    monkeypatch,
    missing_role: str,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)
    directory = tmp_path / "archive" / "fixture" / "官网" / "Method4" / "case-folder"
    target = {
        "result": directory / "arbitrary result capture.HTM",
        "details": directory / "renamed details.HTML",
    }[missing_role]
    target.write_text(
        _html_page(
            "Unrelated page",
            '<a href="official.html">ISODISTORT: distortion</a>'
            '<a href="details.html">ISODISTORT: complete modes details</a>',
        ),
        encoding="utf-8",
    )

    result = official_audit._audit_case(
        "fixture", "C01", case, parent, "case-folder", local,
    )

    assert result["status"] == "fail"
    assert any("missing required" in issue for issue in result["issues"])
    assert result["html_evidence"]["unknown_pages"]


def test_method4_html_parsing_uses_the_captured_inventory(
    tmp_path: Path,
    monkeypatch,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)
    directory = tmp_path / "archive" / "fixture" / "官网" / "Method4" / "case-folder"
    inventory = official_audit.scan_official_html(directory, recursive=True)
    result_page = directory / "arbitrary result capture.HTM"
    result_page.write_text(_html_page("Changed after capture", ""), encoding="utf-8")

    result = official_audit._audit_case(
        "fixture",
        "C01",
        case,
        parent,
        "case-folder",
        local,
        html_inventory=inventory,
    )

    assert result["official_subgroup"]["space_group_number"] == 1
    assert result["status"] == "inconclusive"


def test_nested_diagnostic_run_is_resolved_without_polluting_main_case(
    tmp_path: Path,
    monkeypatch,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)
    directory = tmp_path / "archive" / "fixture" / "官网" / "Method4" / "case-folder"
    diagnostic = directory / "large-displacement-nearest-failed"
    diagnostic.mkdir()
    diagnostic_page = diagnostic / "renamed rejection.HTM"
    diagnostic_page.write_text(
        _method4_result_html(
            "Failed to find match. Try using a larger value for dmax.",
            include_form=False,
        ),
        encoding="utf-8",
    )

    result = official_audit._audit_case(
        "fixture", "C01", case, parent, "case-folder", local,
    )

    assert result["status"] == "inconclusive"
    assert Path(result["selected_html_pages"]["result_page"]["path"]).parent == directory
    diagnostics = result["diagnostic_html_subdirectories"]
    assert len(diagnostics) == 1
    assert diagnostics[0]["status"] == "pass"
    assert (
        Path(diagnostics[0]["selected_pages"]["result_page"]["path"])
        == diagnostic_page
    )


def test_diagnostic_subdirectory_without_result_page_fails_closed(
    tmp_path: Path,
    monkeypatch,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)
    directory = tmp_path / "archive" / "fixture" / "官网" / "Method4" / "case-folder"
    diagnostic = directory / "incomplete-diagnostic"
    diagnostic.mkdir()
    (diagnostic / "capture.png").write_bytes(b"not-a-real-image")

    result = official_audit._audit_case(
        "fixture", "C01", case, parent, "case-folder", local,
    )

    assert result["status"] == "fail"
    assert any("diagnostic HTML directory" in issue for issue in result["issues"])
    assert any("missing required distortion_result" in issue for issue in result["issues"])


def test_selected_case_html_signature_is_rename_stable_and_detects_drift(
    tmp_path: Path,
    monkeypatch,
) -> None:
    _official_success_fixture(tmp_path, monkeypatch)
    manifest = {
        "readable_case_folders": {"C01": "case-folder"},
        "parents": {"fixture": {"cases": {"C01": {}}}},
    }
    directory = tmp_path / "archive" / "fixture" / "官网" / "Method4" / "case-folder"

    def signature() -> str:
        evidence, _inventories, _errors = (
            official_audit._capture_selected_html_evidence(
                manifest,
                parent_filter=None,
                case_filter=None,
            )
        )
        return evidence["sha256"]

    original = signature()
    result_page = directory / "arbitrary result capture.HTM"
    renamed = directory / "manually renamed result.htm"
    result_page.rename(renamed)
    assert signature() == original

    original_html = renamed.read_text(encoding="utf-8")
    renamed.write_text(original_html + "<!-- changed -->", encoding="utf-8")
    assert signature() != original
    renamed.write_text(original_html, encoding="utf-8")

    nested = directory / "diagnostic"
    nested.mkdir()
    moved = nested / renamed.name
    renamed.rename(moved)
    assert signature() != original
    moved.rename(renamed)

    (directory / "unknown.html").write_text(
        _html_page("Unrelated", "new evidence"), encoding="utf-8",
    )
    with_unknown = signature()
    assert with_unknown != original
    (directory / "duplicate-result.html").write_text(
        original_html,
        encoding="utf-8",
    )
    assert signature() != with_unknown


@pytest.mark.parametrize("identity_status", ["inconclusive", "unresolved"])
@pytest.mark.parametrize("overall_status", ["pass", "inconclusive"])
def test_official_case_audits_exports_but_skips_unresolved_local_amplitudes(
    tmp_path: Path, monkeypatch, identity_status: str, overall_status: str,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)
    local["identity_readiness"]["status"] = identity_status
    local["status"] = overall_status
    result = official_audit._audit_case("fixture", "C01", case, parent, "case-folder", local)
    assert result["status"] == "inconclusive"
    assert result["issues"] == []
    assert result["inconclusive_reasons"]
    assert result["export_consistency"]["cif_displacive_mode_count"] == 1
    assert result["exported_cif_vs_frozen_daughter"]["actual_atom_count"] == 1
    assert result["local_vs_official_applied_strain"]["max_absolute_difference"] == 0
    assert result["local_identity_amplitude_comparison"]["status"] == "skipped"
    assert "overlapping_mode_normalization" not in result
    assert "overlapping_mode_amplitudes" not in result


@pytest.mark.parametrize("failed_layer", [
    "inverse", "missing_inverse", "identity", "official", "local_case",
])
def test_official_case_failure_overrides_unresolved_identity(
    tmp_path: Path, monkeypatch, failed_layer: str,
) -> None:
    case, parent, local, comparison = _official_success_fixture(tmp_path, monkeypatch)
    if failed_layer == "inverse":
        local["inverse_consistency"]["status"] = "fail"
    elif failed_layer == "missing_inverse":
        local.pop("inverse_consistency")
        local["status"] = "pass"
    elif failed_layer == "identity":
        local["identity_readiness"]["status"] = "fail"
    elif failed_layer == "official":
        comparison["actual_atom_count"] = 2
    else:
        local["status"] = "fail"
    result = official_audit._audit_case("fixture", "C01", case, parent, "case-folder", local)
    assert result["status"] == "fail"
    assert result["issues"]


def test_official_case_compares_amplitudes_when_local_identity_is_ready(
    tmp_path: Path, monkeypatch,
) -> None:
    case, parent, local, _comparison = _official_success_fixture(tmp_path, monkeypatch)
    local["status"] = "pass"
    local["identity_readiness"] = {"status": "pass", "unresolved": []}
    result = official_audit._audit_case("fixture", "C01", case, parent, "case-folder", local)
    assert result["status"] == "fail"
    assert result["local_identity_amplitude_comparison"]["status"] == "ready"
    assert result["overlapping_mode_amplitudes"]
    assert any("normfactor mismatch" in issue for issue in result["issues"])
    assert any("As/Ap mismatch" in issue for issue in result["issues"])


def test_noisy_official_case_without_matrix_sensitivity_is_inconclusive(
    tmp_path: Path, monkeypatch,
) -> None:
    case, parent, local, comparison = _official_success_fixture(tmp_path, monkeypatch)
    case["noise_sigma_fractional"] = 2e-5
    comparison["origin_aligned_root_sum_squared_distance_angstrom"] = 4e-5
    local["status"] = "pass"
    local["identity_readiness"] = {"status": "pass", "unresolved": []}

    result = official_audit._audit_case("fixture", "C01", case, parent, "case-folder", local)

    assert result["status"] == "inconclusive"
    assert result["issues"] == []
    assert result["shared_mode_amplitude_tolerance"]["status"] == "inconclusive"
    assert result["local_identity_amplitude_comparison"]["status"] == "skipped"
    assert "overlapping_mode_amplitudes" not in result
    assert result["export_consistency"]["cif_displacive_mode_count"] == 1
    assert result["local_vs_official_applied_strain"]["max_absolute_difference"] == 0


@pytest.mark.parametrize("statuses,expected_status,exit_code", [
    (["inconclusive", "pass_with_warnings"], "inconclusive", 2),
    (["inconclusive", "fail"], "fail", 1),
])
def test_official_main_summarizes_inconclusive_and_failure_priority(
    tmp_path: Path, monkeypatch, statuses: list[str], expected_status: str, exit_code: int,
) -> None:
    manifest = {
        "signature": "current-manifest", "local_validation_tolerances": {},
        "readable_case_folders": {"C01": "first", "C02": "second"},
        "parents": {"fixture": {"cases": {"C01": {}, "C02": {}}}},
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    state = {"signature": "current-generation"}
    monkeypatch.setattr(
        official_audit, "build_method4_provenance_state", lambda *_args: state,
    )
    local = _valid_local_report(state, manifest, manifest_path)
    local.update(status="inconclusive", inconclusive=1)
    local_path = tmp_path / "local.json"
    local_path.write_text(json.dumps(local), encoding="utf-8")
    report_path = tmp_path / "official.json"
    observed = []

    def audit_case(*args, **_kwargs):
        observed.append(args[1])
        return {"status": statuses[len(observed) - 1]}

    monkeypatch.setattr(official_audit, "_audit_case", audit_case)
    monkeypatch.setattr(sys, "argv", [
        "audit_method4_official.py", "--manifest", str(manifest_path),
        "--local-report", str(local_path), "--json-output", str(report_path),
    ])
    assert official_audit.main() == exit_code
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["schema"] == 5
    assert report["local_provenance_gate"]["status"] == "pass"
    assert observed == ["C01", "C02"]
    assert report["summary"]["status"] == expected_status
    assert report["summary"]["inconclusive_case_count"] == 1
    assert report["summary"]["failed_case_count"] == statuses.count("fail")
    assert report["summary"]["blocked_case_count"] == 0


@pytest.mark.parametrize(
    ("flag", "value", "expected_issue"),
    [
        ("--parent", "", "--parent must not be empty"),
        ("--parent", "unknown", "unknown --parent value"),
        ("--case-id", "", "--case-id must not be empty"),
        ("--case-id", "unknown", "unknown --case-id value"),
    ],
)
def test_official_main_reports_empty_or_unknown_filter_as_inconclusive(
    tmp_path: Path,
    monkeypatch,
    flag: str,
    value: str,
    expected_issue: str,
) -> None:
    manifest = {
        "signature": "current-manifest",
        "local_validation_tolerances": {},
        "readable_case_folders": {"C01": "case-folder"},
        "parents": {"fixture": {"cases": {"C01": {}}}},
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    state = {"signature": "current-generation"}
    monkeypatch.setattr(
        official_audit, "build_method4_provenance_state", lambda *_args: state,
    )
    local_path = tmp_path / "local.json"
    local_path.write_text(
        json.dumps(_valid_local_report(state, manifest, manifest_path)),
        encoding="utf-8",
    )
    report_path = tmp_path / "official.json"

    def forbidden_case_audit(*_args, **_kwargs):
        raise AssertionError("invalid filters must not audit any case")

    monkeypatch.setattr(official_audit, "_audit_case", forbidden_case_audit)
    monkeypatch.setattr(sys, "argv", [
        "audit_method4_official.py",
        "--manifest", str(manifest_path),
        "--local-report", str(local_path),
        "--json-output", str(report_path),
        flag, value,
    ])

    assert official_audit.main() == 2
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["schema"] == 5
    assert report["selection_gate"]["status"] == "inconclusive"
    assert any(
        expected_issue in issue for issue in report["selection_gate"]["issues"]
    )
    assert report["summary"]["selected_case_count"] == 0
    assert report["summary"]["audited_case_count"] == 0
    assert report["summary"]["downloaded_case_count"] == 0
    assert report["summary"]["status"] == "inconclusive"
    assert report["cases"] == []


@pytest.mark.parametrize(
    "statuses",
    [
        ["not_downloaded"],
        ["pass", "not_downloaded"],
    ],
)
def test_official_main_never_passes_with_not_downloaded_case(
    tmp_path: Path,
    monkeypatch,
    statuses: list[str],
) -> None:
    case_ids = [f"C{index:02d}" for index in range(1, len(statuses) + 1)]
    manifest = {
        "signature": "current-manifest",
        "local_validation_tolerances": {},
        "readable_case_folders": {
            case_id: f"case-{case_id}" for case_id in case_ids
        },
        "parents": {
            "fixture": {"cases": {case_id: {} for case_id in case_ids}},
        },
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    state = {"signature": "current-generation"}
    monkeypatch.setattr(
        official_audit, "build_method4_provenance_state", lambda *_args: state,
    )
    local = _valid_local_report(state, manifest, manifest_path)
    local["case_count"] = len(case_ids)
    local["parents"] = {
        "fixture": {
            "cases": {case_id: {"status": "pass"} for case_id in case_ids},
        },
    }
    local_path = tmp_path / "local.json"
    local_path.write_text(json.dumps(local), encoding="utf-8")
    report_path = tmp_path / "official.json"
    observed: list[str] = []

    def audit_case(*args, **_kwargs):
        observed.append(args[1])
        return {"status": statuses[len(observed) - 1]}

    monkeypatch.setattr(official_audit, "_audit_case", audit_case)
    monkeypatch.setattr(sys, "argv", [
        "audit_method4_official.py",
        "--manifest", str(manifest_path),
        "--local-report", str(local_path),
        "--json-output", str(report_path),
    ])

    assert official_audit.main() == 2
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["schema"] == 5
    assert observed == case_ids
    assert report["summary"]["selected_case_count"] == len(statuses)
    assert report["summary"]["audited_case_count"] == len(statuses)
    assert report["summary"]["downloaded_case_count"] == len(statuses) - 1
    assert report["summary"]["not_downloaded_case_count"] == 1
    assert report["summary"]["status"] == "inconclusive"


def test_official_main_fails_when_html_inventory_changes_during_run(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = {
        "signature": "current-manifest",
        "local_validation_tolerances": {},
        "readable_case_folders": {"C01": "case-folder"},
        "parents": {"fixture": {"cases": {"C01": {}}}},
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    state = {"signature": "current-generation"}
    monkeypatch.setattr(
        official_audit,
        "build_method4_provenance_state",
        lambda *_args: state,
    )
    local = _valid_local_report(state, manifest, manifest_path)
    local_path = tmp_path / "local.json"
    local_path.write_text(json.dumps(local), encoding="utf-8")
    report_path = tmp_path / "official.json"
    archive_root = tmp_path / "archive"
    case_directory = archive_root / "fixture" / "官网" / "Method4" / "case-folder"
    case_directory.mkdir(parents=True)
    (case_directory / "initial.html").write_text(
        _html_page("Unrelated", "initial"),
        encoding="utf-8",
    )
    monkeypatch.setattr(official_audit, "OFFICIAL_ROOT", archive_root)

    def mutate_archive(*_args, **_kwargs):
        (case_directory / "late-added.htm").write_text(
            _html_page("Unrelated", "late"),
            encoding="utf-8",
        )
        return {"status": "pass"}

    monkeypatch.setattr(official_audit, "_audit_case", mutate_archive)
    monkeypatch.setattr(sys, "argv", [
        "audit_method4_official.py",
        "--manifest", str(manifest_path),
        "--local-report", str(local_path),
        "--json-output", str(report_path),
    ])

    assert official_audit.main() == 1
    report = json.loads(report_path.read_text(encoding="utf-8"))
    gate = report["official_html_evidence_gate"]
    assert report["schema"] == 5
    assert gate["status"] == "fail"
    assert gate["stable_during_run"] is False
    assert gate["at_start"]["sha256"] != gate["at_completion"]["sha256"]
    assert report["summary"]["html_evidence_gate_failure_count"] == 1
    assert report["summary"]["status"] == "fail"
