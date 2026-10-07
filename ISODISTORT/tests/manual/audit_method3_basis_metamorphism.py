"""Audit Method 3 invariance under exact GL(3,Z) query-basis changes.

The official Method 3 input is a representative basis of a direct lattice.
Left multiplication by an integer unimodular matrix changes that representative
without changing the lattice.  Each frozen official query is rerun with two
such representatives; the exact content-addressed embedding multiset must be
unchanged.  This is a metamorphic generalization test, not new official gold.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from backend.api import IsoDistort  # noqa: E402
from backend.utils.lattice import (  # noqa: E402
    centering_primitive_matrix,
    multiply,
    rational_matrix,
    same_lattice,
)
from backend.utils.opd_format import _centering_letter  # noqa: E402

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method3_download_manifest.json"
REPORT_PATH = PACKAGE_ROOT / "output" / "validation" / "method3_basis_metamorphic_audit.json"
TRANSFORMS = {
    "cyclic_rows": [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    "signed_ab_rotation": [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _source_signature(manifest: dict[str, Any], manifest_path: Path) -> str:
    """Hash every local input that can change this audit's scientific result."""
    paths = [
        Path(__file__),
        manifest_path,
        PACKAGE_ROOT / "backend" / "api" / "core_api.py",
        PACKAGE_ROOT / "features" / "method1" / "affine_embeddings.py",
        PACKAGE_ROOT / "features" / "method1" / "search_methods.py",
        PACKAGE_ROOT / "backend" / "utils" / "lattice.py",
        PACKAGE_ROOT / "backend" / "utils" / "opd_format.py",
    ]
    paths.extend(
        WORKSPACE / "experiment_data" / parent["parent_cif"]
        for parent in manifest["parents"]
    )
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.resolve().relative_to(WORKSPACE)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _transform_basis(basis: list[list[Any]], transform: list[list[int]]) -> list[list[int]]:
    result = np.asarray(transform, dtype=object) @ np.asarray(basis, dtype=object)
    return [[int(value) for value in row] for row in result.tolist()]


def _identity(item) -> str:
    embedding_id = str(getattr(item, "embedding_id", "") or "").strip()
    if embedding_id:
        return embedding_id
    subgroup = item.subgroup
    fields = subgroup.official_fields()
    return json.dumps(
        {
            key: fields[key]
            for key in ("space_group_number", "basis", "origin", "s", "i")
        },
        sort_keys=True,
        separators=(",", ":"),
    )


def _run(
    iso: IsoDistort,
    case: dict[str, Any],
    basis: list[list[int]],
    centering: str,
) -> list:
    with contextlib.redirect_stdout(io.StringIO()):
        return iso.search_method_3(
            distortion_types=case["distortion_types"],
            space_group_type=case["space_group_type"],
            supercell_basis=basis,
            direct_sublattice_centering=centering,
            lattice_type=case["lattice_type"],
            generate_if_missing=False,
        )


def _equivalent_centering(
    case: dict[str, Any],
    transformed_basis: list[list[int]],
) -> str:
    original_letter = _centering_letter(int(case["space_group_type"]))
    original = multiply(
        centering_primitive_matrix(original_letter),
        rational_matrix(case["supercell_basis"]),
    )
    letters = [original_letter, "P", "A", "B", "C", "I", "F", "R"]
    for letter in dict.fromkeys(letters):
        candidate = multiply(
            centering_primitive_matrix(letter),
            rational_matrix(transformed_basis),
        )
        if same_lattice(candidate, original):
            return "d" if letter == original_letter else letter
    raise ValueError(
        "No website centering option makes the transformed conventional basis "
        "primitive-lattice equivalent to the original query"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--json-output", type=Path, default=REPORT_PATH)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    report: dict[str, Any] = {
        "schema": 1,
        "kind": "method3_exact_basis_metamorphism",
        "started_at": datetime.now(UTC).isoformat(),
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": _sha256(args.manifest),
        "source_signature": _source_signature(manifest, args.manifest),
        "transforms": TRANSFORMS,
        "scope_note": (
            "GL(3,Z) representative-basis invariance over the existing 40-query "
            "official matrix; this does not add new parent crystal systems or prove "
            "the unresolved coupled-IR product path."
        ),
        "parents": [],
    }
    failures = 0
    case_count = 0
    variant_count = 0
    for parent in manifest["parents"]:
        iso = IsoDistort(language="en")
        parent_path = WORKSPACE / "experiment_data" / parent["parent_cif"]
        with contextlib.redirect_stdout(io.StringIO()):
            iso.load_structure(parent_path)
        parent_report = {"parent_cif": parent["parent_cif"], "cases": []}
        report["parents"].append(parent_report)
        for case in parent["cases"]:
            baseline = _run(
                iso,
                case,
                case["supercell_basis"],
                case["direct_sublattice_centering"],
            )
            expected = Counter(_identity(item) for item in baseline)
            case_report = {
                "id": case["id"],
                "baseline_count": len(baseline),
                "baseline_embedding_ids": sorted(expected.elements()),
                "variants": [],
            }
            case_count += 1
            for name, transform in TRANSFORMS.items():
                basis = _transform_basis(case["supercell_basis"], transform)
                try:
                    centering = _equivalent_centering(case, basis)
                    actual_rows = _run(iso, case, basis, centering)
                except Exception as exc:  # noqa: BLE001 - audit records exact failure
                    failures += 1
                    variant_count += 1
                    case_report["variants"].append({
                        "name": name,
                        "basis": basis,
                        "status": "fail",
                        "exception_type": type(exc).__name__,
                        "exception": str(exc),
                    })
                    continue
                actual = Counter(_identity(item) for item in actual_rows)
                missing = list((expected - actual).elements())
                extra = list((actual - expected).elements())
                passed = not missing and not extra
                failures += int(not passed)
                variant_count += 1
                case_report["variants"].append({
                    "name": name,
                    "basis": basis,
                    "direct_sublattice_centering": centering,
                    "candidate_count": len(actual_rows),
                    "missing_embedding_ids": sorted(missing),
                    "extra_embedding_ids": sorted(extra),
                    "status": "pass" if passed else "fail",
                })
            parent_report["cases"].append(case_report)
            case_failures = sum(
                row["status"] != "pass" for row in case_report["variants"]
            )
            print(
                f"Method 3 basis metamorphism {parent['parent_cif']} "
                f"{case['id']}: {len(case_report['variants']) - case_failures}/"
                f"{len(case_report['variants'])}",
                flush=True,
            )
    report.update(
        finished_at=datetime.now(UTC).isoformat(),
        case_count=case_count,
        variant_count=variant_count,
        passed=variant_count - failures,
        failed=failures,
        status="pass" if failures == 0 else "fail",
    )
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "cases": case_count,
        "variants": variant_count,
        "failed": failures,
        "report": str(args.json_output.resolve()),
    }, ensure_ascii=False))
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
