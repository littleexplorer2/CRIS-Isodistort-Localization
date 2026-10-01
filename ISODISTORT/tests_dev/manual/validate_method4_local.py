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

from isocore.api import IsoDistort  # noqa: E402
from isocore.backend import SubgroupInfo  # noqa: E402

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method4_download_manifest.json"
DEFAULT_REPORT = PACKAGE_ROOT / "output" / "validation" / "method4_local_validation.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def _check_success(
    row: dict[str, Any],
    result,
    tolerances: dict[str, float],
) -> tuple[bool, list[str], dict[str, Any]]:
    noisy = row["expect"] == "success_with_nonzero_residual"
    amplitude_tolerance = float(
        tolerances["noise_amplitude_abs" if noisy else "exact_amplitude_abs"]
    )
    expected = {str(key): float(value) for key, value in row["contributions"].items()}
    actual = {str(key): float(value) for key, value in result.raw_coefficients.items()}
    issues = []
    for label, value in expected.items():
        error = abs(actual.get(label, float("inf")) - value)
        if error > amplitude_tolerance:
            issues.append(
                f"amplitude {label!r}: |{actual.get(label)!r} - {value}| "
                f"> {amplitude_tolerance}"
            )
    unexpected = {
        label: value
        for label, value in actual.items()
        if label not in expected and abs(value) > amplitude_tolerance
    }
    if unexpected:
        issues.append(f"unexpected nonzero amplitudes: {unexpected}")
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
        for key, value in result.strain_voigt_engineering.items()
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
        "assignments": [int(value) for value in result.assignments],
        "strain_voigt_engineering": strain,
        "strain_tensor": result.strain_tensor,
        "metadata": result.metadata,
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
    report: dict[str, Any] = {
        "schema": 1,
        "kind": "method4_local_inverse_consistency",
        "started_at": datetime.now(UTC).isoformat(),
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": _sha256(args.manifest),
        "manifest_source_signature": manifest.get("signature"),
        "tolerances": tolerances,
        "scope_note": (
            "Local inverse consistency and rejection only; official As/Ap amplitude "
            "normalization equivalence is not established by this report."
        ),
        "selection": {
            "parent": args.parent,
            "context": args.context,
            "case_id": args.case_id,
        },
        "parents": {},
    }
    failures = 0
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
            parent_report["contexts"][context] = {"mode_labels": labels}
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
                        passed, issues, details = _check_success(row, result, tolerances)
                        case_report.update(
                            status="pass" if passed else "fail",
                            issues=issues,
                            **details,
                        )
                        failures += int(not passed)
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
    report.update(
        finished_at=datetime.now(UTC).isoformat(),
        case_count=len(all_cases),
        passed=sum(case["status"] == "pass" for case in all_cases),
        failed=failures,
        status="pass" if failures == 0 else "fail",
    )
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
        report["manifest_sha256"] = _sha256(args.manifest)
        report["manifest_status_updated"] = True
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "passed": report["passed"],
        "failed": report["failed"],
        "report": str(args.json_output.resolve()),
    }, ensure_ascii=False))
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
