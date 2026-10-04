"""Checkpointed current-source audit of every saved 4310 Method 1 candidate.

The official tree is a read-only reference.  Generated CIFs, the checkpoint,
and the final report are written only under ``ISODISTORT/output/validation``.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import platform
import sys
import time
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[3]
PROJECT_ROOT = WORKSPACE / "ISODISTORT"
MANUAL_ROOT = PROJECT_ROOT / "tests_dev" / "manual"
VALIDATOR_ROOT = WORKSPACE / "ISODISTORT_VALIDATE"
for entry in (PROJECT_ROOT, MANUAL_ROOT, VALIDATOR_ROOT):
    if str(entry) not in sys.path:
        sys.path.insert(0, str(entry))

import validate_method_outputs as validator  # noqa: E402
from isodistort_validate.compare_cif import compare_cif  # noqa: E402

from isocore.api import IsoDistort  # noqa: E402
from isocore.io.distortion_formats import (  # noqa: E402
    render_cif,
    unique_folder_name,
)

PARENT = WORKSPACE / "experiment_data" / "4310_tetra.cif"
OFFICIAL_ROOT = (
    WORKSPACE / "output_compare" / "4310_tetra.cif" / "官网" / "Method1"
)
VALIDATION_ROOT = PROJECT_ROOT / "output" / "validation"
DEFAULT_DEST_ROOT = VALIDATION_ROOT / "method1_4310_live_postfix_cifs_20261003"
DEFAULT_CHECKPOINT = VALIDATION_ROOT / (
    "method1_4310_live_postfix_checkpoint_20261003.json"
)
DEFAULT_REPORT = VALIDATION_ROOT / (
    "method1_4310_live_postfix_audit_20261003.json"
)
DEFAULT_CIF_REANALYSIS_REPORT = VALIDATION_ROOT / (
    "method1_4310_live_postfix_cif_precision_reanalysis_v2_20261004.json"
)
SUPERSEDED_CIF_REANALYSIS_REPORT = VALIDATION_ROOT / (
    "method1_4310_live_postfix_cif_precision_reanalysis_20261003.json"
)
SUPERSEDED_CIF_REANALYSIS_SHA256 = (
    "59bff2a8a7ede1051e8355ed0267b696e62a57e5b882d2631b0419d6f3c5d2cd"
)
SCHEMA = 1


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256(usedforsecurity=False)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_or_none(value: Any) -> Any:
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _json_default(value: Any) -> Any:
    """Convert NumPy scalar results without weakening JSON validation."""
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(
        f"Object of type {value.__class__.__name__} is not JSON serializable"
    )


def _json_safe(value: Any) -> Any:
    """Recursively replace non-finite diagnostic scalars with JSON ``null``.

    Comparison diagnostics can legitimately use infinity for an absent bound.
    JSON has no portable representation for NaN or infinity, so retain that
    distinction as ``null`` while leaving unknown object types for the strict
    encoder to reject.
    """

    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(
            _json_safe(payload),
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            allow_nan=False,
            default=_json_default,
        ),
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(path)


def _classify_cif_result(
    setting: Any,
    *,
    direct_structure_equal: bool,
    equivalent_setting_label: str,
) -> str:
    """Make the precision-aware result authoritative; direct is only a label."""
    if setting.equivalent:
        return "direct" if direct_structure_equal else equivalent_setting_label
    return "inconclusive" if not setting.conclusive else "failed"


def _source_record_set_integrity(
    source: dict[str, Any],
    source_records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], int | None, list[Any], list[Any]]:
    """Validate completeness and uniqueness before trusting a frozen record set."""
    failures: list[dict[str, Any]] = []
    summary = source.get("summary", {})
    identity = source.get("identity", {})
    declared_counts = {
        "summary.completed_records": summary.get("completed_records"),
        "identity.paired_candidates": identity.get("paired_candidates"),
        "identity.official_candidates": identity.get("official_candidates"),
        "identity.live_candidates": identity.get("live_candidates"),
    }
    concrete_counts = {
        label: int(value)
        for label, value in declared_counts.items()
        if value is not None
    }
    expected_record_count = next(iter(concrete_counts.values()), None)
    for label, expected in concrete_counts.items():
        if expected != len(source_records):
            failures.append(
                {
                    "code": "source-record-count-mismatch",
                    "field": label,
                    "expected": expected,
                    "actual": len(source_records),
                }
            )
    if not concrete_counts:
        failures.append({"code": "source-record-count-not-declared"})

    identities = [record.get("identity") for record in source_records]
    positions = [record.get("position") for record in source_records]
    if len(set(identities)) != len(identities) or any(
        value is None for value in identities
    ):
        failures.append({"code": "source-identities-not-unique-and-complete"})
    if len(set(positions)) != len(positions) or any(
        value is None for value in positions
    ):
        failures.append({"code": "source-positions-not-unique-and-complete"})
    if positions and all(isinstance(value, int) for value in positions):
        expected_positions = set(range(1, len(positions) + 1))
        if set(positions) != expected_positions:
            failures.append({"code": "source-positions-not-contiguous"})

    non_cif_source_errors = [
        {
            "identity": record.get("identity"),
            "error": error,
        }
        for record in source_records
        for error in record.get("errors", [])
        if error.get("code") != "cif-semantic-mismatch"
    ]
    if non_cif_source_errors:
        failures.append(
            {
                "code": "source-report-has-non-cif-failures",
                "count": len(non_cif_source_errors),
                "samples": non_cif_source_errors[:8],
            }
        )
    if source.get("signature_unchanged_at_completion") is not True:
        failures.append(
            {"code": "source-signature-was-not-stable-at-completion"}
        )
    if any(
        identity.get(key)
        for key in ("duplicate_live", "duplicate_official", "extra_live", "missing_live")
    ):
        failures.append({"code": "source-candidate-pairing-not-complete"})
    return failures, expected_record_count, identities, positions


def _frozen_report_status(
    records: list[dict[str, Any]],
    integrity_failures: list[dict[str, Any]],
) -> str:
    accepted = {"direct", "precision-equivalent-setting"}
    return (
        "complete-passed"
        if records
        and all(record.get("status") in accepted for record in records)
        and not integrity_failures
        else "complete-failed"
    )


def _frozen_pair_preflight(local_cif: Path, official_cif: Path) -> dict[str, Any]:
    """Strictly parse every frozen input needed to define one U,q comparison."""
    headers = {
        "local": validator.parse_candidate_header(local_cif),
        "official": validator.parse_candidate_header(official_cif),
    }
    official_inverse = validator._inverse(headers["official"].basis)
    transform = validator._multiply(headers["local"].basis, official_inverse)
    origin_shift = validator._vector_multiply(
        tuple(
            headers["local"].origin[index] - headers["official"].origin[index]
            for index in range(3)
        ),
        official_inverse,
    )
    determinant = validator._determinant(transform)
    integer_transform = all(
        value.denominator == 1 for row in transform for value in row
    )
    sides: dict[str, Any] = {}
    for side, path in (("local", local_cif), ("official", official_cif)):
        block = validator._cif_block(path)
        rows = validator._atom_site_rows(block)
        operations = validator._space_group_operations(block)
        operation_tags = [
            tag for tag in validator._SPACE_GROUP_OPERATION_TAGS if tag in block
        ]
        sides[side] = {
            "atom_site_raw_row_count": len(rows),
            "space_group_operation_tags": operation_tags,
            "space_group_declared_operation_count": len(operations),
            "space_group_unique_operation_count": len(
                {validator._operation_key(operation) for operation in operations}
            ),
            "magnetic_moment_tags": validator._magnetic_moment_tags(block),
        }
    return {
        **sides,
        "transform": {
            "unique": True,
            "matrix": validator._fraction_matrix_json(transform),
            "origin_shift": [
                validator._fraction_json(value) for value in origin_shift
            ],
            "determinant": validator._fraction_json(determinant),
            "integer": integer_transform,
            "unimodular": integer_transform and abs(determinant) == 1,
        },
    }


def _frozen_preflight_summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate strict-parser evidence without hiding per-record failures."""
    preflights = [
        record["preflight"] for record in records if "preflight" in record
    ]
    operation_counts: Counter[tuple[str, int, int]] = Counter()
    operation_tags: Counter[tuple[str, tuple[str, ...]]] = Counter()
    row_counts: defaultdict[str, list[int]] = defaultdict(list)
    magnetic_cifs = 0
    for preflight in preflights:
        for side in ("local", "official"):
            details = preflight[side]
            operation_counts[
                (
                    side,
                    details["space_group_declared_operation_count"],
                    details["space_group_unique_operation_count"],
                )
            ] += 1
            operation_tags[
                (side, tuple(details["space_group_operation_tags"]))
            ] += 1
            row_counts[side].append(details["atom_site_raw_row_count"])
            magnetic_cifs += bool(details["magnetic_moment_tags"])
    return {
        "record_pairs_expected": len(records),
        "record_pairs_parsed": len(preflights),
        "parse_error_count": len(records) - len(preflights),
        "unique_transform_count": sum(
            bool(item["transform"]["unique"]) for item in preflights
        ),
        "non_unimodular_transform_count": sum(
            not bool(item["transform"]["unimodular"]) for item in preflights
        ),
        "magnetic_cif_count": magnetic_cifs,
        "atom_site_raw_row_count_ranges": {
            side: [min(values), max(values)]
            for side, values in sorted(row_counts.items())
            if values
        },
        "space_group_operation_tag_counts": [
            {"side": side, "tags": list(tags), "count": count}
            for (side, tags), count in sorted(operation_tags.items())
        ],
        "space_group_operation_count_distribution": [
            {
                "side": side,
                "declared": declared,
                "unique": unique,
                "count": count,
            }
            for (side, declared, unique), count in sorted(operation_counts.items())
        ],
    }


def _signature_state() -> dict[str, Any]:
    paths: list[tuple[str, Path]] = [
        ("audit-script", Path(__file__).resolve()),
        ("saved-output-validator", MANUAL_ROOT / "validate_method_outputs.py"),
        (
            "cif-semantic-validator",
            VALIDATOR_ROOT / "isodistort_validate" / "compare_cif.py",
        ),
        ("parent", PARENT),
        ("config", PROJECT_ROOT / "config" / "settings.yaml"),
    ]
    paths.extend(
        (
            f"source:{path.relative_to(PROJECT_ROOT).as_posix()}",
            path,
        )
        for path in sorted(
            (PROJECT_ROOT / "isocore").rglob("*.py"),
            key=lambda item: item.as_posix().casefold(),
        )
    )
    paths.extend(
        (f"backend:{path.name}", path)
        for path in sorted(
            (PROJECT_ROOT / "isobyu").iterdir(),
            key=lambda item: item.name.casefold(),
        )
        if path.is_file()
        and (
            path.name in {"iso", "smodes", "findsym", "const.dat"}
            or path.name.startswith("data_")
        )
    )
    official_files = sorted(
        (
            path
            for path in OFFICIAL_ROOT.rglob("*")
            if path.is_file()
            and (
                path.name.casefold() == "subgroup.cif"
                or "complete modes details" in path.name.casefold()
            )
        ),
        key=lambda item: item.relative_to(OFFICIAL_ROOT).as_posix().casefold(),
    )
    paths.extend(
        (f"official:{path.relative_to(OFFICIAL_ROOT).as_posix()}", path)
        for path in official_files
    )
    fingerprints = {
        label: {
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for label, path in paths
    }
    versions = {}
    for name in ("numpy", "pymatgen", "spglib", "PyYAML"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "<missing>"
    runtime = {
        "python": sys.version,
        "implementation": sys.implementation.name,
        "platform": platform.platform(),
        "packages": versions,
    }
    digest = hashlib.sha256(usedforsecurity=False)
    digest.update(
        json.dumps(
            {"files": fingerprints, "runtime": runtime},
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    return {
        "signature": digest.hexdigest(),
        "runtime": runtime,
        "files": fingerprints,
        "official_reference_file_count": len(official_files),
    }


def _candidate_fingerprint(subgroup: Any) -> str:
    return validator._live12_candidate_fingerprint(subgroup)


def _completed_record_reusable(record: Any, candidate_fingerprint: str) -> bool:
    if not isinstance(record, dict):
        return False
    if record.get("status") != "passed":
        return False
    if record.get("candidate_fingerprint") != candidate_fingerprint:
        return False
    expected_hash = record.get("generated_cif_sha256")
    generated_folder = record.get("generated_folder")
    if not isinstance(expected_hash, str) or not isinstance(generated_folder, str):
        return False
    local_cif = Path(generated_folder) / "subgroup.cif"
    return local_cif.is_file() and _sha256(local_cif) == expected_hash


def _mode_matrix_diagnostics(
    displacements: dict[str, np.ndarray],
    expected_atoms: int,
) -> dict[str, Any]:
    invalid: list[dict[str, Any]] = []
    flattened: list[np.ndarray] = []
    for label, raw in displacements.items():
        values = np.asarray(raw, dtype=float)
        if values.shape != (expected_atoms, 3):
            invalid.append(
                {
                    "label": label,
                    "code": "shape",
                    "shape": list(values.shape),
                    "expected": [expected_atoms, 3],
                }
            )
            continue
        if not np.all(np.isfinite(values)):
            invalid.append({"label": label, "code": "nonfinite"})
            continue
        norm = float(np.linalg.norm(values))
        if norm <= 1e-12:
            invalid.append({"label": label, "code": "zero"})
            continue
        flattened.append(values.reshape(-1))
    if flattened:
        matrix = np.vstack(flattened)
        singular = np.linalg.svd(matrix, compute_uv=False)
        tolerance = (
            max(matrix.shape) * np.finfo(float).eps * float(singular[0])
            if singular.size
            else 0.0
        )
        rank = int(np.sum(singular > tolerance))
    else:
        tolerance = 0.0
        rank = 0
    return {
        "invalid_vectors": invalid,
        "matrix_rank": rank,
        "rank_tolerance": tolerance,
        "linearly_independent": rank == len(displacements),
    }


def _report_from_checkpoint(checkpoint: dict[str, Any]) -> dict[str, Any]:
    records = list(checkpoint["records"].values())
    status_counts = Counter(record["status"] for record in records)
    mode_matches = sum(
        record.get("official_mode_count") == record.get("live_mode_count")
        for record in records
        if record.get("official_mode_count") is not None
        and record.get("live_mode_count") is not None
    )
    cif_matches = sum(
        record.get("cif_semantic_status") in {"direct", "equivalent-setting"}
        for record in records
    )
    bush_passes = sum(record.get("bush_coverage_passed") is True for record in records)
    bush_na = sum(record.get("bush_coverage_status") == "not-applicable" for record in records)
    failures = [record for record in records if record["status"] != "passed"]
    return {
        "schema": SCHEMA,
        "scope": (
            "Current-source live generation of all saved official 4310 Method 1 "
            "candidates. Official files are read-only; generated CIFs are audit artifacts."
        ),
        "criteria": [
            "Official and live candidate identities pair one-to-one by the complete basis-independent OPD identity.",
            "Every live displacive mode count equals the official Complete modes details count.",
            (
                "Each rooted DISPLAY BUSH mode maps without uncovered/conflicting "
                "child atoms; all final vectors have the child-cell shape, finite "
                "nonzero entries, and full row rank."
            ),
            (
                "Each generated zero-amplitude CIF is consistent with the paired "
                "official CIF at their serialized output precision, directly or "
                "after an exact GL(3,Z) setting transform."
            ),
        ],
        "limitations": [
            "This finite 4310 batch is regression evidence, not proof for every crystal.",
            (
                "CIF comparison checks the zero-amplitude structure and metadata "
                "handled by the validator; it does not independently prove every "
                "numerical mode vector or normalization."
            ),
            "No VESTA/IsoVIZ GUI validation and no official website resubmission are performed.",
            (
                "A no-displacive-mode candidate has no BUSH mapping to test and is "
                "recorded as not-applicable rather than a BUSH pass."
            ),
            (
                "Mutable WSL working/cache files are excluded from the signature; "
                "the 4310 Method 1 set is re-enumerated live and every paired "
                "candidate is recomputed in the same invocation."
            ),
        ],
        "signature_state": checkpoint["signature_state"],
        "signature_unchanged_at_completion": checkpoint.get(
            "signature_unchanged_at_completion", False
        ),
        "started_at": checkpoint["started_at"],
        "updated_at": checkpoint["updated_at"],
        "completed_at": checkpoint.get("completed_at"),
        "elapsed_seconds": checkpoint["elapsed_seconds"],
        "status": checkpoint["status"],
        "identity": checkpoint["identity"],
        "summary": {
            "official_candidates": checkpoint["identity"]["official_candidates"],
            "live_candidates": checkpoint["identity"]["live_candidates"],
            "paired_candidates": checkpoint["identity"]["paired_candidates"],
            "completed_records": len(records),
            "passed_records": status_counts.get("passed", 0),
            "failed_records": len(failures),
            "mode_count_matches": mode_matches,
            "cif_semantic_matches": cif_matches,
            "bush_coverage_passes": bush_passes,
            "bush_not_applicable": bush_na,
            "status_counts": dict(status_counts),
        },
        "failure_samples": failures[:20],
        "records": records,
        "generated_cif_root": checkpoint["generated_cif_root"],
        "checkpoint": checkpoint["checkpoint_path"],
    }


def _reanalyze_saved_cifs(source_report: Path, output_report: Path) -> dict[str, Any]:
    """Rejudge frozen CIF pairs without rerunning candidate or mode generation."""
    if source_report.resolve() == output_report.resolve():
        raise ValueError("reanalysis output must not overwrite its source report")
    source_hash_before = _sha256(source_report)
    source = json.loads(source_report.read_text(encoding="utf-8"))
    code_paths = {
        "audit_script": Path(__file__).resolve(),
        "precision_validator": (MANUAL_ROOT / "validate_method_outputs.py").resolve(),
        "direct_validator": (
            VALIDATOR_ROOT / "isodistort_validate" / "compare_cif.py"
        ).resolve(),
    }
    code_hashes_before = {name: _sha256(path) for name, path in code_paths.items()}
    raw_records = source.get("records")
    if isinstance(raw_records, dict):
        source_records = list(raw_records.values())
    elif isinstance(raw_records, list):
        source_records = raw_records
    else:
        raise ValueError("source report has no records list or mapping")
    source_records.sort(key=lambda record: int(record.get("position", 0)))
    (
        integrity_failures,
        expected_record_count,
        identities,
        positions,
    ) = _source_record_set_integrity(source, source_records)
    superseded_hash_before = (
        _sha256(SUPERSEDED_CIF_REANALYSIS_REPORT)
        if SUPERSEDED_CIF_REANALYSIS_REPORT.is_file()
        else None
    )
    if superseded_hash_before != SUPERSEDED_CIF_REANALYSIS_SHA256:
        integrity_failures.append(
            {
                "code": "superseded-report-missing-or-changed",
                "path": str(SUPERSEDED_CIF_REANALYSIS_REPORT.resolve()),
                "expected_sha256": SUPERSEDED_CIF_REANALYSIS_SHA256,
                "actual_sha256": superseded_hash_before,
            }
        )

    expected_official_hashes = {
        str(Path(item["path"]).resolve()).casefold(): item.get("sha256")
        for label, item in source.get("signature_state", {}).get("files", {}).items()
        if str(label).startswith("official:") and isinstance(item, dict) and item.get("path")
    }

    records: list[dict[str, Any]] = []
    corpus_digest = hashlib.sha256(usedforsecurity=False)
    frozen_hashes_before: dict[str, tuple[Path, str]] = {}
    for source_record in source_records:
        local_cif = Path(source_record["generated_folder"]) / "subgroup.cif"
        official_cif = Path(source_record["official_folder"]) / "subgroup.cif"
        local_hash = _sha256(local_cif)
        official_hash = _sha256(official_cif)
        frozen_hashes_before[str(local_cif.resolve()).casefold()] = (
            local_cif.resolve(),
            local_hash,
        )
        frozen_hashes_before[str(official_cif.resolve()).casefold()] = (
            official_cif.resolve(),
            official_hash,
        )
        record_integrity_issues: list[dict[str, Any]] = []
        expected_local_hash = source_record.get("generated_cif_sha256")
        if expected_local_hash != local_hash:
            record_integrity_issues.append(
                {
                    "code": "frozen-generated-cif-hash-mismatch",
                    "expected": expected_local_hash,
                    "actual": local_hash,
                }
            )
        expected_official_hash = expected_official_hashes.get(
            str(official_cif.resolve()).casefold()
        )
        if expected_official_hash is None:
            record_integrity_issues.append(
                {"code": "official-cif-hash-missing-from-source-signature"}
            )
        elif expected_official_hash != official_hash:
            record_integrity_issues.append(
                {
                    "code": "frozen-official-cif-hash-mismatch",
                    "expected": expected_official_hash,
                    "actual": official_hash,
                }
            )
        corpus_digest.update(
            json.dumps(
                {
                    "identity": source_record["identity"],
                    "local": local_hash,
                    "official": official_hash,
                },
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        )
        record: dict[str, Any] = {
            "position": source_record.get("position"),
            "identity": source_record["identity"],
            "generated_cif": str(local_cif.resolve()),
            "official_cif": str(official_cif.resolve()),
            "generated_cif_sha256": local_hash,
            "official_cif_sha256": official_hash,
            "original_record_status": source_record.get("status"),
            "original_cif_semantic_status": source_record.get(
                "cif_semantic_status"
            ),
            "integrity_issues": record_integrity_issues,
        }
        try:
            record["preflight"] = _frozen_pair_preflight(local_cif, official_cif)
            direct = None
            direct_error = None
            try:
                direct = compare_cif(local_cif, official_cif, ignore_atom_order=True)
            except (OSError, AttributeError, TypeError, ValueError) as exc:
                direct_error = {
                    "type": f"{type(exc).__module__}.{type(exc).__qualname__}",
                    "message": str(exc),
                }
            setting = validator.compare_cif_alternate_settings(
                local_cif,
                official_cif,
                direct_details=direct.details if direct is not None else {},
            )
            status = _classify_cif_result(
                setting,
                direct_structure_equal=bool(
                    direct is not None and direct.structure_equal
                ),
                equivalent_setting_label="precision-equivalent-setting",
            )
            if record_integrity_issues:
                status = "error"
                integrity_failures.append(
                    {
                        "code": "frozen-record-integrity-failure",
                        "identity": source_record.get("identity"),
                        "issues": record_integrity_issues,
                    }
                )
            record.update(
                {
                    "status": status,
                    "direct_comparison": (
                        {
                            "structure_equal": direct.structure_equal,
                            "lattice_equal": direct.details.get("lattice_equal"),
                            "coordinates_equal": direct.details.get("coordinates_equal"),
                            "occupancies_equal": direct.details.get("occupancies_equal"),
                            "declared_space_group_equal": direct.details.get(
                                "declared_space_group_equal"
                            ),
                            "inferred_space_group_equal": direct.details.get(
                                "inferred_space_group_equal"
                            ),
                            "lattice_max_abs_difference": _finite_or_none(
                                direct.details.get("lattice_max_abs_difference")
                            ),
                            "coordinate_max_periodic_difference": _finite_or_none(
                                direct.details.get("coordinate_max_periodic_difference")
                            ),
                            "issues": direct.details.get("issues", []),
                        }
                        if direct is not None
                        else {"error": direct_error}
                    ),
                    "precision_setting_comparison": {
                        "equivalent": setting.equivalent,
                        "conclusive": setting.conclusive,
                        "details": setting.details,
                    },
                }
            )
        except (OSError, TypeError, ValueError) as exc:
            record.update(
                {
                    "status": "error",
                    "error": {
                        "type": f"{type(exc).__module__}.{type(exc).__qualname__}",
                        "message": str(exc),
                    },
                }
            )
        records.append(record)

    preflight_summary = _frozen_preflight_summary(records)
    if preflight_summary["parse_error_count"]:
        integrity_failures.append(
            {
                "code": "frozen-preflight-incomplete",
                "count": preflight_summary["parse_error_count"],
            }
        )
    if preflight_summary["unique_transform_count"] != len(records):
        integrity_failures.append(
            {
                "code": "frozen-setting-transforms-not-all-unique",
                "expected": len(records),
                "actual": preflight_summary["unique_transform_count"],
            }
        )
    if preflight_summary["non_unimodular_transform_count"]:
        integrity_failures.append(
            {
                "code": "frozen-setting-transform-not-unimodular",
                "count": preflight_summary["non_unimodular_transform_count"],
            }
        )

    status_counts = Counter(record["status"] for record in records)
    failed = (
        status_counts.get("failed", 0)
        + status_counts.get("error", 0)
        + status_counts.get("inconclusive", 0)
    )
    source_hash_after = _sha256(source_report)
    if source_hash_after != source_hash_before:
        integrity_failures.append(
            {
                "code": "source-report-changed-during-reanalysis",
                "sha256_before": source_hash_before,
                "sha256_after": source_hash_after,
            }
        )
    code_hashes_after = {name: _sha256(path) for name, path in code_paths.items()}
    code_changes = {
        name: {
            "sha256_before": code_hashes_before[name],
            "sha256_after": code_hashes_after[name],
        }
        for name in code_paths
        if code_hashes_before[name] != code_hashes_after[name]
    }
    if code_changes:
        integrity_failures.append(
            {
                "code": "validator-code-changed-during-reanalysis",
                "files": code_changes,
            }
        )
    frozen_changes: list[dict[str, Any]] = []
    for path_text, (path, expected_hash) in frozen_hashes_before.items():
        actual_hash = _sha256(path) if path.is_file() else None
        if actual_hash != expected_hash:
            frozen_changes.append(
                {
                    "path": path_text,
                    "sha256_before": expected_hash,
                    "sha256_after": actual_hash,
                }
            )
    if frozen_changes:
        integrity_failures.append(
            {
                "code": "frozen-cif-corpus-changed-during-reanalysis",
                "count": len(frozen_changes),
                "samples": frozen_changes[:8],
            }
        )
    superseded_hash_after = (
        _sha256(SUPERSEDED_CIF_REANALYSIS_REPORT)
        if SUPERSEDED_CIF_REANALYSIS_REPORT.is_file()
        else None
    )
    if superseded_hash_after != superseded_hash_before:
        integrity_failures.append(
            {
                "code": "superseded-report-changed-during-reanalysis",
                "sha256_before": superseded_hash_before,
                "sha256_after": superseded_hash_after,
            }
        )
    integrity_ok = not integrity_failures
    report = {
        "schema": 2,
        "kind": "method1-frozen-cif-precision-reanalysis",
        "created_at": _now(),
        "status": _frozen_report_status(records, integrity_failures),
        "scope": (
            "Independent serialized-output-precision consistency reanalysis of the "
            "already generated 4310 Method 1 CIF pairs; no candidate search or mode "
            "generation is run."
        ),
        "criterion": {
            "basis": "U = B_local @ inverse(B_official) is exact integer with abs(det(U)) = 1",
            "origin": "q = (origin_local - origin_official) @ inverse(B_official)",
            "lattice": (
                "a concrete joint cell-parameter witness lies inside both printed "
                "half-last-digit boxes and satisfies G_local = U @ G_official @ U.T"
            ),
            "sites": (
                "raw atom-site rows retain per-row coordinate and occupancy precision; "
                "x_official = x_local @ U + q (mod 1), followed by an exact-species, "
                "occupancy-constrained perfect orbit match"
            ),
            "space_group": (
                "declared Seitz operation sets are exactly conjugated by U,q and must "
                "match; IT-number disagreement fails and spglib inference is diagnostic"
            ),
            "magnetic_moments": (
                "CIFs carrying atom-site moment tags are inconclusive until axial-vector "
                "and setting transformations are implemented"
            ),
        },
        "limitations": [
            (
                "This checks zero-amplitude CIF serialized structural consistency, "
                "not every mode-vector component or normalization."
            ),
            (
                "Passing means the serialized values are mutually consistent within "
                "their printed precision; it does not recover or prove identical "
                "unrounded structures."
            ),
            "This finite 4310 batch is regression evidence, not proof for every crystal.",
        ],
        "source_report": {
            "path": str(source_report.resolve()),
            "sha256_before": source_hash_before,
            "sha256_after": source_hash_after,
            "unchanged": source_hash_before == source_hash_after,
            "status": source.get("status"),
            "signature": source.get("signature_state", {}).get("signature")
            or source.get("signature"),
        },
        "supersedes": {
            "path": str(SUPERSEDED_CIF_REANALYSIS_REPORT.resolve()),
            "sha256_before": superseded_hash_before,
            "sha256_after": superseded_hash_after,
            "unchanged": superseded_hash_before == superseded_hash_after,
            "status": "invalidated",
            "reasons": [
                "componentwise Gram-interval overlap did not prove a joint cell witness",
                "the direct comparator could bypass the precision-aware decision",
                "atom-site occupancy and precision were not preserved per raw loop row",
                "declared Seitz operation equivalence was not proved under U,q",
                "magnetic moments and report-integrity failures were not hard gates",
            ],
        },
        "code": {
            name: {
                "path": str(path),
                "sha256": code_hashes_after[name],
                "sha256_before": code_hashes_before[name],
                "sha256_after": code_hashes_after[name],
                "unchanged": code_hashes_before[name] == code_hashes_after[name],
            }
            for name, path in code_paths.items()
        },
        "integrity": {
            "passed": integrity_ok,
            "expected_record_count": expected_record_count,
            "actual_record_count": len(records),
            "unique_identity_count": len(set(identities)),
            "unique_position_count": len(set(positions)),
            "failures": integrity_failures,
        },
        "preflight": preflight_summary,
        "cif_corpus_sha256": corpus_digest.hexdigest(),
        "summary": {
            "records": len(records),
            "direct_matches": status_counts.get("direct", 0),
            "precision_equivalent_setting_matches": status_counts.get(
                "precision-equivalent-setting", 0
            ),
            "failed": status_counts.get("failed", 0),
            "errors": status_counts.get("error", 0),
            "inconclusive": status_counts.get("inconclusive", 0),
            "consistent_at_serialized_output_precision": len(records) - failed,
            "status_counts": dict(status_counts),
        },
        "records": records,
    }
    _atomic_json(output_report, report)
    return report


def _save(
    checkpoint: dict[str, Any],
    started: float,
    checkpoint_path: Path,
    report_path: Path,
) -> None:
    checkpoint["updated_at"] = _now()
    checkpoint["elapsed_seconds"] = round(time.perf_counter() - started, 6)
    _atomic_json(checkpoint_path, checkpoint)
    _atomic_json(report_path, _report_from_checkpoint(checkpoint))


def _stale_path(path: Path, signature: str) -> Path:
    return path.with_name(f"{path.stem}.stale-{signature[:16]}{path.suffix}")


def _archive_incompatible_state(
    previous: dict[str, Any],
    checkpoint_path: Path,
    report_path: Path,
    replacement_signature: str,
) -> None:
    old_signature = str(previous.get("signature") or "unknown")
    previous["status"] = "stale-source-signature"
    previous["stale_at"] = _now()
    previous["stale_reason"] = "source/input/runtime signature changed before completion"
    previous["superseded_by_signature"] = replacement_signature
    _atomic_json(_stale_path(checkpoint_path, old_signature), previous)
    try:
        old_report = _report_from_checkpoint(previous)
    except (KeyError, TypeError, ValueError):
        old_report = {
            "schema": SCHEMA,
            "signature": old_signature,
            "status": "stale-source-signature",
            "stale_reason": previous["stale_reason"],
            "superseded_by_signature": replacement_signature,
            "raw_checkpoint": previous,
        }
    _atomic_json(_stale_path(report_path, old_signature), old_report)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--dest-root", type=Path, default=DEFAULT_DEST_ROOT)
    parser.add_argument(
        "--reanalyze-existing-report",
        type=Path,
        help="rejudge CIFs recorded in a frozen report without rerunning modes",
    )
    parser.add_argument(
        "--reanalyze-output",
        type=Path,
        default=DEFAULT_CIF_REANALYSIS_REPORT,
        help="new machine report written by --reanalyze-existing-report",
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        help="run only the first N paired identities as an explicitly partial smoke",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        help="ignore a compatible checkpoint and recompute requested candidates",
    )
    args = parser.parse_args(argv)
    if args.max_candidates is not None and args.max_candidates < 1:
        parser.error("--max-candidates must be at least 1")
    if args.reanalyze_existing_report is not None:
        source_report = args.reanalyze_existing_report.expanduser().resolve()
        output_report = args.reanalyze_output.expanduser().resolve()
        result = _reanalyze_saved_cifs(source_report, output_report)
        print(json.dumps(result["summary"], ensure_ascii=False, indent=2), flush=True)
        print(f"REPORT {output_report}", flush=True)
        print(f"REPORT_SHA256 {_sha256(output_report)}", flush=True)
        return 0 if result["status"] == "complete-passed" else 1
    checkpoint_path = args.checkpoint.expanduser().resolve()
    report_path = args.report.expanduser().resolve()
    dest_base = args.dest_root.expanduser().resolve()
    started = time.perf_counter()
    signature_state = _signature_state()
    signature = signature_state["signature"]
    dest_root = dest_base / signature[:16]
    checkpoint: dict[str, Any] = {
        "schema": SCHEMA,
        "signature": signature,
        "signature_state": signature_state,
        "status": "running",
        "started_at": _now(),
        "updated_at": _now(),
        "elapsed_seconds": 0.0,
        "generated_cif_root": str(dest_root),
        "checkpoint_path": str(checkpoint_path),
        "identity": {},
        "records": {},
    }
    if checkpoint_path.is_file():
        try:
            previous = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = None
        if (
            not args.restart
            and
            isinstance(previous, dict)
            and previous.get("schema") == SCHEMA
            and previous.get("signature") == signature
            and isinstance(previous.get("records"), dict)
        ):
            checkpoint = previous
            checkpoint["signature_state"] = signature_state
            checkpoint["status"] = "running"
            checkpoint["generated_cif_root"] = str(dest_root)
            checkpoint["checkpoint_path"] = str(checkpoint_path)
            print(
                f"Resuming {len(checkpoint['records'])} checkpointed candidate(s)",
                flush=True,
            )
        elif isinstance(previous, dict):
            _archive_incompatible_state(
                previous,
                checkpoint_path,
                report_path,
                signature,
            )

    official, official_duplicates = validator._index_candidates(OFFICIAL_ROOT)

    api = IsoDistort(language="en")
    api.set_distortion_scope(
        {
            "displacive": ["*"],
            "occupational": [],
            "strain": ["*"],
            "magnetic": [],
            "rotational": [],
        }
    )
    with contextlib.redirect_stdout(io.StringIO()):
        api.load_structure(PARENT)
        rows = api.search_method_1(distortion_types=["strain", "displacive"])
    candidates = [row.subgroup for row in rows]
    live_grouped: dict[str, list[Any]] = defaultdict(list)
    for subgroup in candidates:
        live_grouped[validator._live_subgroup_identity(subgroup)].append(subgroup)
    live = {
        identity: values[0]
        for identity, values in live_grouped.items()
        if len(values) == 1
    }
    live_duplicates = {
        identity: [_candidate_fingerprint(value) for value in values]
        for identity, values in live_grouped.items()
        if len(values) > 1
    }
    paired = sorted(official.keys() & live.keys())
    checkpoint["identity"] = {
        "official_candidates": len(official),
        "live_candidates": len(live_grouped),
        "paired_candidates": len(paired),
        "missing_live": sorted(official.keys() - live.keys()),
        "extra_live": sorted(live.keys() - official.keys()),
        "duplicate_official": [
            {
                "identity": identity,
                "folders": [str(path) for path in paths],
            }
            for identity, paths in official_duplicates
        ],
        "duplicate_live": live_duplicates,
        "parametric_live_candidates": sum(
            bool(getattr(subgroup, "k_parameters", None))
            for subgroup in candidates
        ),
    }
    _save(checkpoint, started, checkpoint_path, report_path)

    used_names: set[str] = set()
    folder_names = {
        identity: unique_folder_name(
            live[identity],
            used_names,
            export_method=1,
            sequence=position,
        )
        for position, identity in enumerate(paired, start=1)
    }
    dest_root.mkdir(parents=True, exist_ok=True)

    requested = paired
    if args.max_candidates is not None:
        requested = paired[: args.max_candidates]

    for position, identity in enumerate(requested, start=1):
        subgroup = live[identity]
        fingerprint = _candidate_fingerprint(subgroup)
        previous = checkpoint["records"].get(identity)
        if _completed_record_reusable(previous, fingerprint):
            if position % 10 == 0 or position == len(requested):
                print(
                    f"4310 Method1 {position}/{len(requested)} (checkpoint)",
                    flush=True,
                )
            continue

        candidate_started = time.perf_counter()
        official_folder = official[identity]
        official_count = validator._count_official_displacive_modes(official_folder)
        folder = dest_root / folder_names[identity]
        folder.mkdir(parents=True, exist_ok=True)
        local_cif = folder / "subgroup.cif"
        errors: list[dict[str, Any]] = []
        record: dict[str, Any] = {
            "position": position,
            "identity": identity,
            "candidate_fingerprint": fingerprint,
            "official_folder": str(official_folder.resolve()),
            "generated_folder": str(folder.resolve()),
            "official_mode_count": official_count,
            "live_mode_count": None,
            "rooted_bush_mode_count": None,
            "rootless_supplement_mode_count": None,
            "bush_coverage_status": "not-run",
            "bush_coverage_passed": False,
            "cif_semantic_status": "not-run",
            "errors": errors,
        }
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                api.search_method_2(
                    subgroup_idx=subgroup.index,
                    distortion_type=["displacive"],
                    number_of_independent_modulations=0,
                    candidates=candidates,
                )
            child = api._supercell_for_subgroup(subgroup)
            displacements = api._lifted_mode_displacements(subgroup)
            live_count = len(displacements)
            record["live_mode_count"] = live_count
            rooted = [mode for mode in api.distortion_modes if mode.bush_modes]
            rooted_keys = {api._mode_session_key(mode) for mode in rooted}
            final_keys = set(displacements)
            missing_rooted = sorted(rooted_keys - final_keys)
            vector_diagnostics = _mode_matrix_diagnostics(displacements, len(child))
            record.update(
                {
                    "child_atom_count": len(child),
                    "rooted_bush_mode_count": len(rooted_keys),
                    "rootless_supplement_mode_count": live_count - len(rooted_keys),
                    "missing_rooted_mode_keys": missing_rooted,
                    "vector_diagnostics": vector_diagnostics,
                }
            )
            all_mode_keys = {
                api._mode_session_key(mode) for mode in api.distortion_modes
            }
            missing_mode_keys = sorted(all_mode_keys - final_keys)
            unexpected_mode_keys = sorted(final_keys - all_mode_keys)
            record["missing_distortion_mode_keys"] = missing_mode_keys
            record["unexpected_final_mode_keys"] = unexpected_mode_keys
            if missing_mode_keys or unexpected_mode_keys:
                errors.append(
                    {
                        "code": "mode-key-set-mismatch",
                        "missing": missing_mode_keys,
                        "unexpected": unexpected_mode_keys,
                    }
                )
            if not rooted_keys:
                record["bush_coverage_status"] = "not-applicable"
                record["bush_coverage_passed"] = None
            elif not missing_rooted and not vector_diagnostics["invalid_vectors"]:
                # map_bush_modes_to_supercell raises on every uncovered child atom
                # or conflicting arrow.  Reaching this branch proves all rooted
                # BUSH modes traversed that coverage gate.
                record["bush_coverage_status"] = "passed"
                record["bush_coverage_passed"] = True
            else:
                record["bush_coverage_status"] = "failed"
                errors.append(
                    {
                        "code": "bush-or-vector-coverage",
                        "missing_rooted_mode_keys": missing_rooted,
                        "invalid_vectors": vector_diagnostics["invalid_vectors"],
                    }
                )
            if not vector_diagnostics["linearly_independent"]:
                errors.append(
                    {
                        "code": "rank-deficient-mode-basis",
                        "mode_count": live_count,
                        "matrix_rank": vector_diagnostics["matrix_rank"],
                    }
                )
            if official_count != live_count:
                errors.append(
                    {
                        "code": "mode-count-mismatch",
                        "official": official_count,
                        "live": live_count,
                    }
                )

            spec = api._spec_for_subgroup(
                subgroup,
                use_current_modes=True,
                use_generated_structure=False,
                folder_name=folder.name,
            )
            cif_text = render_cif(spec.structure, spec)
            temporary = local_cif.with_suffix(".cif.tmp")
            temporary.write_text(cif_text, encoding="utf-8", newline="\n")
            temporary.replace(local_cif)
            local_header = validator.parse_candidate_header(local_cif)
            official_header = validator.parse_candidate_header(
                official_folder / "subgroup.cif"
            )
            record["generated_cif_sha256"] = _sha256(local_cif)
            record["generated_identity"] = local_header.identity
            record["basis_same_lattice"] = validator.bases_generate_same_lattice(
                local_header.basis, official_header.basis
            )
            if not record["basis_same_lattice"]:
                errors.append({"code": "non-equivalent-basis"})
            if local_header.identity != identity:
                errors.append(
                    {
                        "code": "generated-identity-mismatch",
                        "expected": identity,
                        "generated": local_header.identity,
                    }
                )
            comparison = None
            direct_comparison_error = None
            try:
                comparison = compare_cif(
                    local_cif,
                    official_folder / "subgroup.cif",
                    ignore_atom_order=True,
                )
            except (OSError, AttributeError, TypeError, ValueError) as exc:
                direct_comparison_error = {
                    "type": f"{type(exc).__module__}.{type(exc).__qualname__}",
                    "message": str(exc),
                }
            record["cif_comparison"] = (
                {
                    "direct_structure_equal": comparison.structure_equal,
                    "details": comparison.details,
                }
                if comparison is not None
                else {"error": direct_comparison_error}
            )
            setting_comparison = validator.compare_cif_alternate_settings(
                local_cif,
                official_folder / "subgroup.cif",
                direct_details=comparison.details if comparison is not None else {},
            )
            record["precision_setting_comparison"] = {
                "equivalent": setting_comparison.equivalent,
                "conclusive": setting_comparison.conclusive,
                "details": setting_comparison.details,
            }
            record["cif_semantic_status"] = _classify_cif_result(
                setting_comparison,
                direct_structure_equal=bool(
                    comparison is not None and comparison.structure_equal
                ),
                equivalent_setting_label="equivalent-setting",
            )
            if record["cif_semantic_status"] == "inconclusive":
                errors.append(
                    {
                        "code": "cif-semantic-inconclusive",
                        "reasons": setting_comparison.details.get(
                            "inconclusive_reasons", []
                        ),
                    }
                )
            elif record["cif_semantic_status"] == "failed":
                errors.append(
                    {
                        "code": "cif-semantic-mismatch",
                        "issues": setting_comparison.details.get("issues", []),
                    }
                )
        except Exception as exc:  # noqa: BLE001 - record and continue full audit
            errors.append(
                {
                    "code": "candidate-exception",
                    "type": f"{type(exc).__module__}.{type(exc).__qualname__}",
                    "message": str(exc),
                }
            )
        record["elapsed_seconds"] = round(
            time.perf_counter() - candidate_started, 6
        )
        record["status"] = "passed" if not errors else "failed"
        checkpoint["records"][identity] = record
        _save(checkpoint, started, checkpoint_path, report_path)
        print(
            f"4310 Method1 {position}/{len(requested)} {record['status']} "
            f"modes={record.get('live_mode_count')}/{official_count} "
            f"cif={record.get('cif_semantic_status')} "
            f"elapsed={record['elapsed_seconds']:.2f}s",
            flush=True,
        )

    identity_ok = not any(
        checkpoint["identity"][key]
        for key in (
            "missing_live",
            "extra_live",
            "duplicate_official",
            "duplicate_live",
        )
    )
    full_run = args.max_candidates is None or args.max_candidates >= len(paired)
    records_complete = len(checkpoint["records"]) == len(paired)
    records_pass = records_complete and all(
        record.get("status") == "passed"
        for record in checkpoint["records"].values()
    )
    end_signature_state = _signature_state()
    signature_unchanged = end_signature_state["signature"] == signature
    checkpoint["signature_unchanged_at_completion"] = signature_unchanged
    checkpoint["completion_signature"] = end_signature_state["signature"]
    checkpoint["completed_at"] = _now()
    checkpoint["status"] = (
        "complete-passed"
        if full_run and identity_ok and records_pass and signature_unchanged
        else "complete-source-drift"
        if not signature_unchanged
        else "partial-smoke-passed"
        if not full_run
        and all(
            checkpoint["records"][identity].get("status") == "passed"
            for identity in requested
        )
        else "complete-failed"
    )
    _save(checkpoint, started, checkpoint_path, report_path)
    summary = _report_from_checkpoint(checkpoint)["summary"]
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    print(f"REPORT {report_path}", flush=True)
    return (
        0
        if checkpoint["status"] in {"complete-passed", "partial-smoke-passed"}
        else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
