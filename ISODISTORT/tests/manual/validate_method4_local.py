"""Run the frozen Method 4 local generate/decompose validation matrix.

This is an inverse-consistency and input-rejection audit.  It deliberately
does not claim official ISODISTORT amplitude normalization equivalence; that
requires the separately archived official Method 4 output for the same CIFs.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import sys
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from backend.api import IsoDistort  # noqa: E402
from backend.wrappers import SubgroupInfo  # noqa: E402
from features.input_cif.cif_io import read_cif  # noqa: E402
from tests.manual.method4_provenance import (  # noqa: E402
    PROVENANCE_KIND,
    PROVENANCE_SCHEMA,
    build_method4_provenance_state,
    changed_provenance_surfaces,
    sha256_file,
)

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method4_download_manifest.json"
DEFAULT_REPORT = PACKAGE_ROOT / "output" / "validation" / "method4_local_validation.json"


def _sha256(path: Path) -> str:
    return sha256_file(path)


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _float_hex_rows(values: Any) -> list[list[str]]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 2:
        raise ValueError("Method 4 basis arrays must be two-dimensional")
    if not bool(np.all(np.isfinite(array))):
        raise ValueError("Method 4 basis arrays must be finite")
    return [[float(value).hex() for value in row] for row in array]


def _mode_key(mode: Any) -> str:
    return str(mode.amplitude_key or mode.irrep_label)


def _identity_record(mode: Any) -> dict[str, Any]:
    identity = getattr(mode, "mode_identity", None)
    if identity is None:
        return {
            "status": "unresolved",
            "reason": "mode_identity_missing",
        }
    return {
        "stable_token": identity.stable_token,
        "parent_sg": int(identity.parent_sg),
        "global_irrep": str(identity.global_irrep),
        "k_coordinates": [str(value) for value in identity.k_coordinates],
        "wyckoff_letter": str(identity.wyckoff_letter),
        "orbit_id": str(identity.orbit_id),
        "site_irrep": identity.site_irrep,
        "copy_index": identity.copy_index,
        "component_index": identity.component_index,
        "component_label": identity.component_label,
        "source": str(identity.source),
        "status": str(identity.status),
        "reason": identity.reason,
    }


def _identity_readiness(iso: IsoDistort, labels: list[str]) -> dict[str, Any]:
    """Require one complete scientific identity for every installed column."""

    by_key: dict[str, list[Any]] = {}
    for mode in iso.distortion_modes:
        by_key.setdefault(_mode_key(mode), []).append(mode)
    records: dict[str, Any] = {}
    failures: list[str] = []
    unresolved: list[str] = []
    verified_tokens: dict[str, str] = {}
    for label in labels:
        matches = by_key.get(label, [])
        if len(matches) != 1:
            failures.append(
                f"mode key {label!r} resolves to {len(matches)} mode definitions"
            )
            records[label] = {
                "status": "invalid",
                "reason": "mode_key_not_unique",
                "match_count": len(matches),
            }
            continue
        mode = matches[0]
        record = _identity_record(mode)
        records[label] = record
        if record["status"] != "verified":
            unresolved.append(
                f"mode {label!r} identity is {record['status']}: {record.get('reason')}"
            )
            continue
        required = (
            record["global_irrep"],
            record["k_coordinates"],
            record["wyckoff_letter"],
            record["orbit_id"],
            record["site_irrep"],
            record["component_index"],
            record["component_label"],
        )
        if (
            not required[0]
            or len(required[1]) != 3
            or not required[2]
            or not required[3]
            or not required[4]
            or required[5] is None
            or not required[6]
        ):
            failures.append(f"verified mode {label!r} has incomplete identity metadata")
            continue
        token = str(record["stable_token"])
        previous = verified_tokens.get(token)
        if previous is not None:
            failures.append(
                f"modes {previous!r} and {label!r} share identity token {token}"
            )
        else:
            verified_tokens[token] = label
    extra_keys = sorted(set(by_key).difference(labels))
    if extra_keys:
        failures.append(f"installed mode definitions have no displacement columns: {extra_keys}")
    status = "fail" if failures else "inconclusive" if unresolved else "pass"
    return {
        "status": status,
        "issues": failures,
        "unresolved": unresolved,
        "mode_count": len(labels),
        "modes": records,
    }


def _raw_basis_contract(
    iso: IsoDistort,
    reference,
    mode_displacements: dict[str, Any],
) -> dict[str, Any]:
    """Fingerprint the signed/scaled raw coordinate basis and its atom frame."""

    subgroup = getattr(iso, "_selected_subgroup", None)
    ordered_labels = list(mode_displacements)
    identities = {
        _mode_key(mode): _identity_record(mode)
        for mode in iso.distortion_modes
    }
    payload = {
        "schema": 1,
        "embedding": {
            "basis": getattr(subgroup, "basis_vectors", None),
            "origin": getattr(subgroup, "origin", None),
            "space_group_number": getattr(subgroup, "space_group_number", None),
        },
        "atom_frame": {
            "lattice": _float_hex_rows(reference.lattice.matrix),
            "species": [str(site.specie) for site in reference],
            "fractional_coordinates": _float_hex_rows(reference.frac_coords),
        },
        "ordered_columns": [
            {
                "label": label,
                "shape": list(np.asarray(mode_displacements[label]).shape),
                "fractional_vectors": _float_hex_rows(mode_displacements[label]),
                "identity": identities.get(label),
            }
            for label in ordered_labels
        ],
    }
    return {
        "schema": 1,
        "sha256": _payload_sha256(payload),
        "mode_count": len(ordered_labels),
        "atom_count": len(reference),
        "ordered_labels": ordered_labels,
    }


def _prepare_context(
    iso: IsoDistort,
    gamma: dict[str, Any],
    parameter: dict[str, Any],
    context: str,
    resolved_context: dict[str, Any] | None = None,
) -> list[str]:
    target = gamma if context == "gamma" else parameter
    if context == "gamma":
        rows = iso.search_method_1(
            distortion_types=["strain", "displacive"],
            subgroup_space_group=int(target["sg"]),
        )
    elif context == "parameter":
        if resolved_context is None:
            rows = iso.search_method_3(
                distortion_types=["strain", "displacive"],
                space_group_type=int(target["sg"]),
                supercell_basis=target["basis"],
                direct_sublattice_centering="P",
                lattice_type="direct",
                generate_if_missing=False,
            )
        else:
            identity = resolved_context["candidate_identity"]
            candidate = SubgroupInfo(
                index=0,
                space_group_number=int(resolved_context["space_group_number"]),
                space_group_symbol=str(resolved_context["space_group_symbol"]),
                subgroup_index=int(resolved_context["subgroup_index"]),
                size=int(resolved_context["size"]),
                opd_symbol=str(identity["opd"]),
                basis_vectors=resolved_context["basis_vectors"],
                origin=resolved_context["origin"],
                k_point_label=str(resolved_context["k_point_label"]),
                irrep_label=str(identity["irrep"]),
                k_parameters=list(resolved_context["k_parameters"]),
                k_coordinates=list(resolved_context["k_coordinates"]),
                parent_sg=int((iso.symmetry_info or {})["space_group_number"]),
                opd_dir_raw=str(identity["dir"]),
                k_active_raw=str(identity["k_active"]),
            )
            rows = [SimpleNamespace(subgroup=candidate)]
    else:
        raise ValueError(f"Unknown Method 4 context {context!r}")
    candidate = next(
        item.subgroup
        for item in rows
        if item.subgroup.irrep_label == target["irrep"]
        and item.subgroup.opd_symbol == target["opd"]
        and item.subgroup.space_group_number == int(target["sg"])
    )
    iso.search_method_2(
        candidate.index,
        distortion_type=["displacive"],
        number_of_independent_modulations=0,
        candidates=[item.subgroup for item in rows],
    )
    return list(iso.mode_displacements_sc or iso.mode_displacements)


def _failure_matches(expectation: str, exc: Exception) -> bool:
    message = str(exc).lower()
    if expectation == "reject_species_or_stoichiometry":
        return "species" in message or "atom count" in message
    if expectation == "reject_incompatible_lattice":
        return "lattice" in message
    if expectation == "reject_unmatched_atom":
        return "cannot match" in message
    return False


def _raw_coordinate_roundtrip(
    row: dict[str, Any],
    result,
    tolerances: dict[str, float],
    *,
    frozen_contract: dict[str, Any] | None,
    current_contract: dict[str, Any],
) -> dict[str, Any]:
    """Compare raw coefficients only in the identical signed/scaled basis.

    A low physical residual proves inverse consistency, but it does not make
    coefficients from two independently normalized bases comparable.  Older
    Method 4 manifests intentionally lack a complete basis fingerprint, so
    their raw values remain diagnostic rather than pass/fail evidence.
    """

    expected = {str(key): float(value) for key, value in row["contributions"].items()}
    actual = {str(key): float(value) for key, value in result.raw_coefficients.items()}
    report: dict[str, Any] = {
        "status": "not_comparable",
        "reason": "generation_basis_fingerprint_absent",
        "expected_generation_coordinates": expected,
        "current_coordinates": actual,
        "frozen_basis_contract": frozen_contract,
        "current_basis_contract": current_contract,
    }
    if not isinstance(frozen_contract, dict):
        return report
    if (
        frozen_contract.get("schema") != current_contract.get("schema")
        or frozen_contract.get("sha256") != current_contract.get("sha256")
    ):
        report["reason"] = "generation_and_current_basis_fingerprints_differ"
        return report

    noisy = row["expect"] == "success_with_nonzero_residual"
    tolerance = float(
        tolerances["noise_amplitude_abs" if noisy else "exact_amplitude_abs"]
    )
    issues: list[str] = []
    for label, value in expected.items():
        error = abs(actual.get(label, float("inf")) - value)
        if error > tolerance:
            issues.append(
                f"raw coordinate {label!r}: |{actual.get(label)!r} - {value}| "
                f"> {tolerance}"
            )
    unexpected = {
        label: value
        for label, value in actual.items()
        if label not in expected and abs(value) > tolerance
    }
    if unexpected:
        issues.append(f"unexpected nonzero raw coordinates: {unexpected}")
    report.update(
        status="fail" if issues else "pass",
        reason="identical_basis_fingerprint",
        tolerance=tolerance,
        issues=issues,
    )
    return report


def _reconstructed_max_displacement(
    result,
    mode_displacements: dict[str, Any],
    reference_lattice: Any,
) -> float:
    if set(result.raw_coefficients) != set(mode_displacements):
        raise ValueError("Method 4 coefficients and basis columns use different keys")
    displacement = np.zeros_like(
        np.asarray(next(iter(mode_displacements.values())), dtype=float)
    )
    for label, coefficient in result.raw_coefficients.items():
        displacement += float(coefficient) * np.asarray(
            mode_displacements[label], dtype=float
        )
    cartesian = displacement @ np.asarray(reference_lattice, dtype=float)
    return float(np.max(np.linalg.norm(cartesian, axis=1)))


def _mode_source_quantization_bound(
    iso: IsoDistort,
    result,
    mode_displacements: dict[str, Any],
    reference_lattice: Any,
) -> dict[str, Any]:
    """Propagate ISO source-token intervals through the actual mode matrix.

    ``DISPLAY DISTORTION`` coefficients are lexical measurements: ``0.577``
    denotes a central value with half-step 0.0005, whereas ``1/2`` is exact.
    For a mapped source row the parent-space-group rotation is Cartesian
    orthogonal, so the row's Cartesian interval radius is unchanged.  The
    bounds below then account for child max-component normalization, primitive
    RSS normalization, and the measured smallest singular value of the full
    unit-mode matrix.  If direct-source evidence is unavailable, the audit
    fails closed instead of inventing a tolerance.
    """

    labels = list(mode_displacements)
    by_key = {_mode_key(mode): mode for mode in iso.distortion_modes}
    parent_lattice = np.asarray(iso.structure.lattice.matrix, dtype=float)
    child_lattice = np.asarray(reference_lattice, dtype=float)
    inverse_child_norm = float(np.linalg.norm(np.linalg.inv(child_lattice), ord=2))
    primitive_multiplicity = int(result.metadata["primitive_cell_multiplicity"])
    if primitive_multiplicity <= 0:
        return {"status": "unavailable", "reason": "invalid_primitive_multiplicity"}
    sqrt_multiplicity = float(np.sqrt(primitive_multiplicity))
    mode_records: dict[str, Any] = {}
    unit_columns: list[np.ndarray] = []
    unit_column_error_bounds: list[float] = []

    for label in labels:
        mode = by_key.get(label)
        if mode is None:
            return {
                "status": "unavailable",
                "reason": f"mode_definition_missing:{label}",
            }
        provenance = getattr(mode, "microscopic_provenance", None)
        if provenance is None:
            return {
                "status": "unavailable",
                "reason": f"microscopic_source_provenance_missing:{label}",
            }
        if getattr(mode, "microscopic_domain_extension", None) is not None:
            return {
                "status": "unavailable",
                "reason": f"extended_source_precision_not_proven:{label}",
            }
        half_steps = tuple(provenance.source_component_half_steps)
        if len(half_steps) != len(mode.bush_modes):
            return {
                "status": "unavailable",
                "reason": f"source_precision_row_mismatch:{label}",
            }
        source_vectors = np.asarray(
            [bush.displacements[0] for bush in mode.bush_modes], dtype=float,
        )
        source_scale_lower = float(
            np.max(np.linalg.norm(source_vectors, axis=1))
        )
        if source_scale_lower <= 0:
            return {
                "status": "unavailable",
                "reason": f"zero_source_column:{label}",
            }
        row_cartesian_error = max(
            (
                float(sum(
                    float(component) * float(np.linalg.norm(parent_lattice[index]))
                    for index, component in enumerate(vector)
                ))
                for vector in half_steps
            ),
            default=0.0,
        )
        fractional = np.asarray(mode_displacements[label], dtype=float)
        internal_max = float(np.max(np.abs(fractional)))
        if internal_max <= 0:
            return {
                "status": "unavailable",
                "reason": f"zero_mapped_column:{label}",
            }
        canonical_fractional = fractional / internal_max
        canonical_cartesian = canonical_fractional @ child_lattice
        canonical_cartesian_norm = float(np.linalg.norm(canonical_cartesian))
        if canonical_cartesian_norm <= 0:
            return {
                "status": "unavailable",
                "reason": f"zero_canonical_column:{label}",
            }

        # Every target atom can inherit the same representative-row rounding
        # error under exact group transport.  This deliberately uses the full
        # child atom count, yielding a conservative Frobenius bound without
        # guessing which zero-looking entries were printed as exact zeros.
        unnormalized_cartesian_error = (
            float(np.sqrt(len(fractional))) * row_cartesian_error
        )
        child_fractional_component_error = (
            row_cartesian_error * inverse_child_norm
        )
        max_component_lower = source_scale_lower * internal_max
        denominator = max_component_lower - child_fractional_component_error
        if denominator <= 0:
            return {
                "status": "unavailable",
                "reason": f"source_interval_crosses_zero_scale:{label}",
            }
        canonical_cartesian_error = float(
            (
                unnormalized_cartesian_error
                + canonical_cartesian_norm * child_fractional_component_error
            )
            / denominator
        )
        primitive_norm = canonical_cartesian_norm / sqrt_multiplicity
        primitive_norm_error = canonical_cartesian_error / sqrt_multiplicity
        if primitive_norm_error >= primitive_norm:
            return {
                "status": "unavailable",
                "reason": f"source_interval_crosses_zero_norm:{label}",
            }
        normfactor_error = float(
            primitive_norm_error
            / (primitive_norm * (primitive_norm - primitive_norm_error))
        )
        if canonical_cartesian_error >= canonical_cartesian_norm:
            return {
                "status": "unavailable",
                "reason": f"source_interval_destroys_mode_direction:{label}",
            }
        unit_column_error = float(
            sqrt_multiplicity
            * 2.0
            * canonical_cartesian_error
            / (canonical_cartesian_norm - canonical_cartesian_error)
        )
        unit_columns.append(
            (canonical_cartesian / primitive_norm).reshape(-1)
        )
        unit_column_error_bounds.append(unit_column_error)
        mode_records[label] = {
            "source_component_half_steps": [
                [str(value) for value in vector] for vector in half_steps
            ],
            "source_scale_lower_parent_fractional": source_scale_lower,
            "source_row_cartesian_error_bound_angstrom": row_cartesian_error,
            "canonical_column_cartesian_error_bound_angstrom": (
                canonical_cartesian_error
            ),
            "unit_column_l2_error_bound": unit_column_error,
            "normfactor_error_bound_inverse_angstrom": normfactor_error,
        }

    unit_matrix = np.column_stack(unit_columns)
    singular = np.linalg.svd(unit_matrix, compute_uv=False)
    sigma_min = float(singular[-1]) if singular.size else 0.0
    matrix_error = float(np.linalg.norm(unit_column_error_bounds))
    residual_rss = float(
        result.rms_residual * np.sqrt(unit_matrix.shape[0])
    )
    amplitude_norm = float(np.linalg.norm([
        float(result.amplitudes[label]) for label in labels
    ]))
    denominator = sigma_min - matrix_error
    if denominator <= 0:
        return {
            "status": "unavailable",
            "reason": "source_interval_exceeds_mode_matrix_sigma_min",
            "mode_records": mode_records,
            "unit_mode_matrix_smallest_singular_value": sigma_min,
            "unit_mode_matrix_l2_error_bound": matrix_error,
        }
    as_error = float(
        (residual_rss + matrix_error * amplitude_norm) / denominator
    )
    supercell_size = float(result.metadata["supercell_size"])
    ap_error = as_error / float(np.sqrt(supercell_size))
    return {
        "status": "available",
        "source": "DISPLAY DISTORTION lexical component half-steps",
        "matrix_norm": "Cartesian spectral/Frobenius upper bound",
        "mode_records": mode_records,
        "unit_mode_matrix_smallest_singular_value": sigma_min,
        "unit_mode_matrix_l2_error_bound": matrix_error,
        "local_residual_rss_angstrom": residual_rss,
        "local_As_vector_l2_norm_angstrom": amplitude_norm,
        "As_vector_l2_error_bound_angstrom": as_error,
        "Ap_vector_l2_error_bound_angstrom": ap_error,
        "derivation": (
            "Weyl/least-squares bound: (residual_rss + epsilon_A*||As||) / "
            "(sigma_min(A)-epsilon_A); Ap bound divides by sqrt(supercell_size)"
        ),
    }


def _check_success(
    row: dict[str, Any],
    result,
    tolerances: dict[str, float],
    *,
    iso: IsoDistort,
    mode_displacements: dict[str, Any],
    reference_lattice: Any,
) -> tuple[bool, list[str], dict[str, Any]]:
    noisy = row["expect"] == "success_with_nonzero_residual"
    expected = {str(key): float(value) for key, value in row["contributions"].items()}
    actual = {str(key): float(value) for key, value in result.raw_coefficients.items()}
    issues = []
    if not all(np.isfinite(value) for value in actual.values()):
        issues.append("raw mode coordinates contain non-finite values")
    matrix_rank = int(result.metadata["mode_matrix_rank"])
    matrix_columns = int(result.metadata["mode_matrix_columns"])
    if matrix_rank != matrix_columns or matrix_columns != len(actual):
        issues.append(
            "mode matrix is not one full-rank column basis: "
            f"rank={matrix_rank}, columns={matrix_columns}, coefficients={len(actual)}"
        )
    rms = float(result.rms_residual)
    if noisy:
        lower = float(tolerances["noise_rms_residual_angstrom_min"])
        upper = float(tolerances["noise_rms_residual_angstrom_max"])
        if not lower < rms < upper:
            issues.append(f"noise RMS {rms} angstrom is outside ({lower}, {upper})")
    else:
        upper = float(tolerances["exact_rms_residual_angstrom_max"])
        if rms > upper:
            issues.append(f"exact-case RMS {rms} angstrom exceeds {upper}")
    strain = {
        str(key): float(value)
        for key, value in (
            result.strain_applied_engineering_q_parent_basis.items()
        )
    }
    strain_norm = float(np.linalg.norm(list(strain.values())))
    strain_tolerance = float(tolerances.get("strain_component_abs", 5e-7))
    if row.get("expect_nonzero_strain"):
        if strain_norm <= strain_tolerance:
            issues.append(
                "homogeneous-strain case returned a zero applied strain tensor"
            )
    elif max((abs(value) for value in strain.values()), default=0.0) > strain_tolerance:
        issues.append(f"unstrained case has nonzero applied strain: {strain}")
    strain_metric_residual = float(
        result.metadata["strain_reconstruction_relative_metric_residual"]
    )
    if strain_metric_residual > float(
        tolerances.get("strain_metric_relative_residual", 1e-10)
    ):
        issues.append(
            "strain tensor does not reconstruct the daughter lattice metric: "
            f"{strain_metric_residual}"
        )
    reconstructed_dmax = _reconstructed_max_displacement(
        result,
        mode_displacements,
        reference_lattice,
    )
    target_dmax = row.get("target_max_displacement_angstrom")
    dmax_tolerance = max(
        5.0 * float(tolerances["exact_rms_residual_angstrom_max"]),
        10.0 * float(result.max_abs_residual),
    )
    if target_dmax is not None and abs(reconstructed_dmax - float(target_dmax)) > dmax_tolerance:
        issues.append(
            "reconstructed maximum displacement differs from the frozen physical "
            f"target: |{reconstructed_dmax} - {float(target_dmax)}| > {dmax_tolerance}"
        )
    details = {
        "expected_amplitudes": expected,
        "actual_raw_coefficients": actual,
        "actual_As_angstrom": {
            str(key): float(value) for key, value in result.amplitudes.items()
        },
        "actual_Ap_angstrom": {
            str(key): float(value)
            for key, value in result.parent_cell_amplitudes.items()
        },
        "mode_normfactors_inverse_angstrom": {
            str(key): float(value) for key, value in result.mode_normfactors.items()
        },
        "rms_residual_angstrom": rms,
        "max_abs_residual_angstrom": float(result.max_abs_residual),
        "reconstructed_max_displacement_angstrom": reconstructed_dmax,
        "target_max_displacement_angstrom": target_dmax,
        "target_max_displacement_tolerance_angstrom": dmax_tolerance,
        "assignments": [int(value) for value in result.assignments],
        "strain_mode_amplitudes": result.strain_mode_amplitudes,
        "strain_modes": result.strain_modes,
        "strain_raw_coordinate_sum_parent_basis": (
            result.strain_raw_coordinate_sum_parent_basis
        ),
        "strain_applied_engineering_q_parent_basis": strain,
        "strain_tensor_parent_basis": result.strain_tensor_parent_basis,
        "strain_multiplier_parent_basis": result.strain_multiplier_parent_basis,
        "metadata": result.metadata,
        "mode_source_quantization_bound": _mode_source_quantization_bound(
            iso,
            result,
            mode_displacements,
            reference_lattice,
        ),
    }
    return not issues, issues, details


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--update-manifest-status",
        action="store_true",
        help="record the completed local run in the tracked input manifest",
    )
    parser.add_argument("--parent")
    parser.add_argument("--context", choices=("gamma", "parameter"))
    parser.add_argument("--case-id")
    args = parser.parse_args()
    if args.update_manifest_status and any(
        value is not None for value in (args.parent, args.context, args.case_id)
    ):
        raise ValueError(
            "--update-manifest-status requires an unfiltered full validation run"
        )
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    tolerances = manifest["local_validation_tolerances"]
    selection = {
        "parent": args.parent,
        "context": args.context,
        "case_id": args.case_id,
    }
    provenance_at_start = build_method4_provenance_state(
        args.manifest,
        selection,
    )
    report: dict[str, Any] = {
        "schema": 3,
        "kind": "method4_local_inverse_consistency",
        "started_at": datetime.now(UTC).isoformat(),
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": _sha256(args.manifest),
        "manifest_source_signature": manifest.get("signature"),
        "tolerances": tolerances,
        "scope_note": (
            "Local inverse consistency, rejection, and scientific identity readiness. "
            "Raw generation coordinates are compared only when a complete signed/scaled "
            "basis fingerprint is identical. Official As/Ap amplitude normalization "
            "equivalence remains a separate archived-site differential."
        ),
        "selection": selection,
        "parents": {},
    }
    failures = 0
    inconclusive = 0
    for parent_name, parent_row in manifest["parents"].items():
        if args.parent and parent_name != args.parent:
            continue
        iso = IsoDistort(language="en")
        parent_path = WORKSPACE / parent_row["parent_cif"]
        with contextlib.redirect_stdout(io.StringIO()):
            iso.load_structure(parent_path)
        parent_report = {"contexts": {}, "cases": {}}
        report["parents"][parent_name] = parent_report
        selected_rows = {
            case_id: row
            for case_id, row in parent_row["cases"].items()
            if (args.context is None or row["context"] == args.context)
            and (args.case_id is None or case_id == args.case_id)
        }
        contexts = {row["context"] for row in selected_rows.values()}
        for context in sorted(contexts):
            with contextlib.redirect_stdout(io.StringIO()):
                labels = _prepare_context(
                    iso,
                    parent_row["gamma_context"],
                    parent_row["parameter_context"],
                    context,
                    parent_row.get("resolved_contexts", {}).get(context),
                )
            identity_readiness = _identity_readiness(iso, labels)
            frozen_context = parent_row.get("resolved_contexts", {}).get(context, {})
            parent_report["contexts"][context] = {
                "mode_labels": labels,
                "identity_readiness": identity_readiness,
            }
            for case_id, row in selected_rows.items():
                if row["context"] != context:
                    continue
                path = WORKSPACE / row["daughter_cif"]
                case_report: dict[str, Any] = {
                    "expect": row["expect"],
                    "daughter_cif": str(path.resolve()),
                    "sha256_expected": row["sha256"],
                    "sha256_actual": _sha256(path) if path.is_file() else None,
                }
                if case_report["sha256_actual"] != row["sha256"]:
                    case_report.update(status="fail", issues=["daughter CIF hash mismatch"])
                    failures += 1
                    parent_report["cases"][case_id] = case_report
                    continue
                kwargs = {
                    key: row[key]
                    for key in (
                        "atom_matching_method",
                        "robust_distance_threshold",
                        "provided_origin_shift",
                    )
                    if key in row
                }
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        result = iso.search_method_4(path, **kwargs)
                except Exception as exc:  # noqa: BLE001 - validation records exact failure
                    passed = row["expect"].startswith("reject_") and _failure_matches(
                        row["expect"], exc,
                    )
                    case_report.update(
                        status="pass" if passed else "fail",
                        exception_type=type(exc).__name__,
                        exception=str(exc),
                    )
                    failures += int(not passed)
                else:
                    if row["expect"].startswith("reject_"):
                        case_report.update(
                            status="fail",
                            issues=["expected rejection but decomposition succeeded"],
                        )
                        failures += 1
                    else:
                        distorted = read_cif(path)
                        reference, mode_displacements = iso._method4_reference_modes(
                            distorted
                        )
                        current_contract = _raw_basis_contract(
                            iso,
                            reference,
                            mode_displacements,
                        )
                        inverse_passed, inverse_issues, details = _check_success(
                            row,
                            result,
                            tolerances,
                            iso=iso,
                            mode_displacements=mode_displacements,
                            reference_lattice=reference.lattice.matrix,
                        )
                        raw_roundtrip = _raw_coordinate_roundtrip(
                            row,
                            result,
                            tolerances,
                            frozen_contract=frozen_context.get("raw_basis_contract"),
                            current_contract=current_contract,
                        )
                        issues = list(inverse_issues)
                        if raw_roundtrip["status"] == "fail":
                            issues.extend(raw_roundtrip.get("issues", []))
                        if not inverse_passed or raw_roundtrip["status"] == "fail":
                            status = "fail"
                        elif identity_readiness["status"] == "fail":
                            status = "fail"
                            issues = [
                                *issues,
                                *identity_readiness["issues"],
                            ]
                        elif identity_readiness["status"] == "inconclusive":
                            status = "inconclusive"
                        else:
                            status = "pass"
                        case_report.update(
                            status=status,
                            issues=issues,
                            inverse_consistency={
                                "status": "pass" if inverse_passed else "fail",
                                "issues": inverse_issues,
                            },
                            identity_readiness=identity_readiness,
                            raw_coordinate_roundtrip=raw_roundtrip,
                            **details,
                        )
                        failures += int(status == "fail")
                        inconclusive += int(status == "inconclusive")
                parent_report["cases"][case_id] = case_report

        unprepared = IsoDistort(language="en")
        with contextlib.redirect_stdout(io.StringIO()):
            unprepared.load_structure(parent_path)
        if not selected_rows:
            continue
        probe = next(iter(selected_rows.values()))
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                unprepared.search_method_4(WORKSPACE / probe["daughter_cif"])
        except RuntimeError as exc:
            parent_report["unprepared_modes"] = {
                "status": "pass",
                "exception": str(exc),
            }
        else:
            parent_report["unprepared_modes"] = {"status": "fail"}
            failures += 1

    all_cases = [
        case
        for parent in report["parents"].values()
        for case in parent["cases"].values()
    ]
    finished_at = datetime.now(UTC).isoformat()
    provenance_before_manifest_update = build_method4_provenance_state(
        args.manifest,
        selection,
    )
    changed_during_run = changed_provenance_surfaces(
        provenance_at_start,
        provenance_before_manifest_update,
    )
    provenance_stable = not changed_during_run
    report.update(
        finished_at=finished_at,
        case_count=len(all_cases),
        passed=sum(case["status"] == "pass" for case in all_cases),
        inconclusive=inconclusive,
        failed=failures,
        status=(
            "fail"
            if failures or not provenance_stable
            else "inconclusive"
            if inconclusive
            else "pass"
        ),
    )
    manifest_self_update: dict[str, Any] = {
        "requested": bool(args.update_manifest_status),
        "performed": False,
        "allowed_changed_surface": "input-manifest",
    }
    if args.update_manifest_status:
        manifest["method4_executed"] = True
        manifest["local_validation_status"] = report["status"]
        manifest["local_validation_report"] = str(
            args.json_output.relative_to(WORKSPACE),
        ).replace("\\", "/")
        manifest["local_validation_completed_at"] = report["finished_at"]
        args.manifest.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        manifest_self_update["performed"] = True
        report["manifest_status_updated"] = True
    provenance_at_completion = build_method4_provenance_state(
        args.manifest,
        selection,
    )
    post_validation_changes = changed_provenance_surfaces(
        provenance_before_manifest_update,
        provenance_at_completion,
    )
    allowed_post_validation = {"input-manifest"} if args.update_manifest_status else set()
    unexpected_post_validation = [
        label for label in post_validation_changes if label not in allowed_post_validation
    ]
    provenance_stable = provenance_stable and not unexpected_post_validation
    if not provenance_stable:
        report["status"] = "fail"
        if args.update_manifest_status and manifest.get("local_validation_status") != "fail":
            manifest["local_validation_status"] = "fail"
            args.manifest.write_text(
                json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            provenance_at_completion = build_method4_provenance_state(
                args.manifest,
                selection,
            )
            post_validation_changes = changed_provenance_surfaces(
                provenance_before_manifest_update,
                provenance_at_completion,
            )
            unexpected_post_validation = [
                label
                for label in post_validation_changes
                if label not in allowed_post_validation
            ]
    manifest_self_update["changed_surfaces"] = post_validation_changes
    manifest_self_update["unexpected_changed_surfaces"] = unexpected_post_validation
    report["manifest_sha256"] = _sha256(args.manifest)
    report["provenance_failure_count"] = int(not provenance_stable)
    report["provenance"] = {
        "schema": PROVENANCE_SCHEMA,
        "kind": PROVENANCE_KIND,
        "signature": provenance_at_completion["signature"],
        "state": provenance_at_completion,
        "run_start_signature": provenance_at_start["signature"],
        "run_end_signature_before_manifest_update": (
            provenance_before_manifest_update["signature"]
        ),
        "stable_during_run": provenance_stable,
        "changed_surfaces_during_run": changed_during_run,
        "manifest_self_update": manifest_self_update,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "passed": report["passed"],
        "inconclusive": report["inconclusive"],
        "failed": report["failed"],
        "provenance_failed": bool(report["provenance_failure_count"]),
        "report": str(args.json_output.resolve()),
    }, ensure_ascii=False))
    return {"pass": 0, "fail": 1, "inconclusive": 2}[report["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
