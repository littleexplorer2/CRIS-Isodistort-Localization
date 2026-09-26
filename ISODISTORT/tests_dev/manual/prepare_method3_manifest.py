"""Freeze and preflight the 40-query Method 3 official comparison batch.

The manifest is deliberately independent of directory ordering in the saved
Method 1/2 references.  Each case is an explicit website query and names one
candidate that must occur in the local result.  Preflight progress is written
atomically after every query so an interrupted WSL run can resume safely.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import os
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from isocore.api import IsoDistort  # noqa: E402

OUTPUT_DIR = PACKAGE_ROOT / "output" / "validation"
MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method3_download_manifest.json"
CHECKPOINT_PATH = OUTPUT_DIR / "method3_preflight_checkpoint.json"
IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
PARENT_SPACE_GROUP_TYPES = {
    "EuAl4 Parent.cif": 139,
    "NdNiO2 own.cif": 123,
}


def _case(
    case_id: str,
    sg: int,
    basis: list[list[int]],
    expected_ir: str,
    expected_opd: str,
    centering: str,
    category: str,
) -> dict[str, Any]:
    return {
        "id": case_id,
        "readable_case_folder": f"{case_id} - SG{sg} - {expected_ir} {expected_opd}",
        "space_group_type": sg,
        "supercell_basis": basis,
        "direct_sublattice_centering": centering,
        "lattice_type": "direct",
        "distortion_types": ["strain", "displacive"],
        "expected_candidate": {"irrep": expected_ir, "opd": expected_opd, "sg": sg},
        "category": category,
    }


# Eight unit-volume queries (identity or an axis-oriented representative), eight
# special-k supercells and four commensurate parameter-k queries per parent.
# The two seed queries already smoke-tested in docs/DEVELOPMENT_PLAN.md are
# M3-EU-17 and M3-ND-17.
CASES: dict[str, list[dict[str, Any]]] = {
    "EuAl4 Parent.cif": [
        _case("M3-EU-01", 139, IDENTITY, "GM1+", "P1", "d", "identity"),
        _case("M3-EU-02", 71, IDENTITY, "GM2+", "P1", "d", "identity"),
        _case("M3-EU-03", 107, IDENTITY, "GM3-", "P1", "d", "identity"),
        _case("M3-EU-04", 12, [[0, 1, 1], [-1, 0, 0], [0, -1, 0]], "GM5+", "P1", "d", "oriented-cell"),
        _case("M3-EU-05", 8, [[-1, 1, 0], [0, 0, 1], [1, 0, 0]], "GM5-", "C1", "d", "oriented-cell"),
        _case("M3-EU-06", 44, [[0, 1, 0], [0, 0, 1], [1, 0, 0]], "GM5-", "P1", "d", "oriented-cell"),
        _case("M3-EU-07", 123, IDENTITY, "M1+", "P1", "d", "identity"),
        _case("M3-EU-08", 129, IDENTITY, "M3-", "P1", "d", "identity"),
        _case("M3-EU-09", 47, [[1, 1, 0], [-1, 1, 0], [0, 0, 1]], "X1+", "C1", "d", "special-supercell"),
        _case("M3-EU-10", 125, [[1, -1, 0], [1, 1, 0], [0, 0, 1]], "X1-", "P1", "d", "special-supercell"),
        _case("M3-EU-11", 119, [[1, 1, 0], [-1, 1, 0], [0, 0, 2]], "P1", "C1", "d", "special-supercell"),
        _case("M3-EU-12", 140, [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]], "P2", "P3", "d", "special-supercell"),
        _case("M3-EU-13", 2, [[2, 0, 0], [0, 2, 0], [-1, -1, 1]], "N1+", "4D1", "d", "special-supercell"),
        _case("M3-EU-14", 65, [[2, 0, 0], [0, 0, -2], [0, 1, 0]], "N1+", "P1", "d", "special-supercell"),
        _case("M3-EU-15", 122, [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]], "P5", "P1", "d", "special-supercell"),
        _case("M3-EU-16", 24, [[0, 0, -2], [-1, -1, 0], [-1, 1, 0]], "P5", "C3", "d", "special-supercell"),
        _case("M3-EU-17", 99, [[1, 0, 0], [0, 1, 0], [0, 0, 6]], "LD1", "C1", "d", "parameter-k"),
        _case("M3-EU-18", 123, [[1, 0, 0], [0, 1, 0], [0, 0, 6]], "LD1", "P1", "d", "parameter-k"),
        _case("M3-EU-19", 105, [[1, 0, 0], [0, 1, 0], [0, 0, 6]], "LD2", "C1", "d", "parameter-k"),
        _case("M3-EU-20", 20, [[-1, 1, 0], [-1, -1, 0], [0, 0, 6]], "LD5", "C4", "d", "parameter-k"),
    ],
    "NdNiO2 own.cif": [
        _case("M3-ND-01", 123, IDENTITY, "GM1+", "P1", "d", "identity"),
        _case("M3-ND-02", 47, IDENTITY, "GM2+", "P1", "d", "identity"),
        _case("M3-ND-03", 99, IDENTITY, "GM3-", "P1", "d", "identity"),
        _case("M3-ND-04", 115, IDENTITY, "GM4-", "P1", "d", "identity"),
        _case("M3-ND-05", 2, [[0, 0, 1], [1, 0, 0], [0, 1, 0]], "GM5+", "C1", "d", "oriented-cell"),
        _case("M3-ND-06", 10, [[0, 0, -1], [-1, 0, 0], [0, 1, 0]], "GM5+", "P1", "d", "oriented-cell"),
        _case("M3-ND-07", 6, [[0, 1, 0], [0, 0, 1], [1, 0, 0]], "GM5-", "C1", "d", "oriented-cell"),
        _case("M3-ND-08", 25, [[0, 1, 0], [0, 0, 1], [1, 0, 0]], "GM5-", "P1", "d", "oriented-cell"),
        _case("M3-ND-09", 123, [[1, 1, 0], [-1, 1, 0], [0, 0, 1]], "M1+", "P1", "d", "special-supercell"),
        _case("M3-ND-10", 127, [[1, 1, 0], [-1, 1, 0], [0, 0, 1]], "M2+", "P1", "d", "special-supercell"),
        _case("M3-ND-11", 129, [[1, 1, 0], [-1, 1, 0], [0, 0, 1]], "M2-", "P1", "d", "special-supercell"),
        _case("M3-ND-12", 139, [[1, 1, 0], [-1, 1, 0], [0, 0, 2]], "A1+", "P1", "d", "special-supercell"),
        _case("M3-ND-13", 140, [[1, 1, 0], [-1, 1, 0], [0, 0, 2]], "A2+", "P1", "d", "special-supercell"),
        _case("M3-ND-14", 47, [[2, 0, 0], [0, 2, 0], [0, 0, 1]], "X1+", "C1", "d", "special-supercell"),
        _case("M3-ND-15", 123, [[1, 0, 0], [0, 1, 0], [0, 0, 2]], "Z1+", "P1", "d", "special-supercell"),
        _case("M3-ND-16", 11, [[0, 1, 0], [0, 0, 2], [1, 0, 0]], "Z5+", "C1", "d", "special-supercell"),
        _case("M3-ND-17", 47, [[-3, 0, 0], [0, 0, 1], [0, 2, 0]], "Y1", "P1", "d", "parameter-k"),
        _case("M3-ND-18", 25, [[0, 2, 0], [0, 0, 1], [3, 0, 0]], "Y1", "C1", "d", "parameter-k"),
        _case("M3-ND-19", 49, [[-3, 0, 0], [0, 0, 1], [0, 2, 0]], "Y2", "P1", "d", "parameter-k"),
        _case("M3-ND-20", 51, [[0, 2, 0], [-3, 0, 0], [0, 0, 1]], "Y2", "P2", "d", "parameter-k"),
    ],
}


def _validate_case_matrix() -> None:
    ids: set[str] = set()
    for parent, cases in CASES.items():
        if len(cases) != 20:
            raise AssertionError(f"{parent}: expected 20 Method 3 cases, got {len(cases)}")
        categories: dict[str, int] = {}
        centerings: set[str] = set()
        for case in cases:
            if case["id"] in ids:
                raise AssertionError(f"duplicate Method 3 case id: {case['id']}")
            ids.add(case["id"])
            categories[case["category"]] = categories.get(case["category"], 0) + 1
            centerings.add(case["direct_sublattice_centering"])
            if case["expected_candidate"]["sg"] != case["space_group_type"]:
                raise AssertionError(f"{case['id']}: expected SG does not match query")
        expected_categories = {
            "EuAl4 Parent.cif": {"identity": 5, "oriented-cell": 3, "special-supercell": 8, "parameter-k": 4},
            "NdNiO2 own.cif": {"identity": 4, "oriented-cell": 4, "special-supercell": 8, "parameter-k": 4},
        }
        if categories != expected_categories[parent]:
            raise AssertionError(f"{parent}: invalid coverage {categories}")
        if centerings != {"d"}:
            raise AssertionError(f"{parent}: official download cases must all use Default centering")


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temp, path)


def _signature() -> str:
    digest = hashlib.sha256()
    for path in [
        Path(__file__),
        PACKAGE_ROOT / "isocore" / "api" / "core_api.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "search_methods.py",
        PACKAGE_ROOT / "isocore" / "backend" / "iso_wrapper.py",
        WORKSPACE / "experiment_data" / "EuAl4 Parent.cif",
        WORKSPACE / "experiment_data" / "NdNiO2 own.cif",
    ]:
        digest.update(path.relative_to(WORKSPACE).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _base_manifest(signature: str) -> dict[str, Any]:
    return {
        "schema": 2,
        "purpose": "Exact official-site download batch for Method 3 differential validation",
        "signature": signature,
        # The manifest freezes scientific inputs only.  Download/comparison
        # status is intentionally kept in the dated audit reports so editing
        # evidence cannot invalidate or recursively rewrite this input file.
        "download_status": "see_output_validation_audit_reports",
        "official_download_root": "output_compare/<parent>/官网/Method3/<readable-case-folder>/",
        "local_download_root": "output_compare/<parent>/现有网页版交互/Method3/<readable-case-folder>/",
        "official_centering_instruction": (
            "Select the radio button immediately after Default for every case. "
            "Do not select P: it means primitive/no centering on the official site."
        ),
        "per_case_downloads": [
            "complete Method 3 result page (HTML)",
            (
                "one subgroup ZIP for every returned candidate; the expected candidate alone "
                "is not sufficient when the table has additional rows"
            ),
        ],
        "preflight_semantics": (
            "Local single-IR subset presence check only. local_subset_candidate_count is not "
            "an official expected Method 3 embedding count or proof of completeness."
        ),
        "parents": [
            {
                "parent_cif": parent,
                "parent_space_group_type": PARENT_SPACE_GROUP_TYPES[parent],
                "case_count": len(cases),
                "cases": cases,
            }
            for parent, cases in CASES.items()
        ],
        "preflight": {"status": "not_run", "completed": 0, "total": 40, "results": {}},
    }


def _candidate_record(item: Any) -> dict[str, Any]:
    sg = item.subgroup
    return {
        "irrep": sg.irrep_label,
        "opd": sg.opd_symbol,
        "sg": sg.space_group_number,
        "symbol": sg.space_group_symbol,
        "basis": sg.basis_vectors,
        "origin": sg.origin,
        "k_active": sg.official_fields()["k_active"],
    }


def _run_preflight(manifest: dict[str, Any]) -> int:
    signature = manifest["signature"]
    checkpoint: dict[str, Any] = {"signature": signature, "results": {}}
    if CHECKPOINT_PATH.exists():
        saved = json.loads(CHECKPOINT_PATH.read_text(encoding="utf-8"))
        if saved.get("signature") == signature:
            checkpoint = saved
        else:
            print("Ignoring incompatible Method 3 checkpoint", flush=True)

    results: dict[str, Any] = checkpoint["results"]
    total = sum(len(cases) for cases in CASES.values())
    done = len(results)
    for parent_name, cases in CASES.items():
        iso = IsoDistort()
        with contextlib.redirect_stdout(io.StringIO()):
            iso.load_structure(WORKSPACE / "experiment_data" / parent_name)
        for case in cases:
            case_id = case["id"]
            if results.get(case_id, {}).get("status") == "pass":
                continue
            expected = case["expected_candidate"]
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    items = iso.search_method_3(
                        distortion_types=case["distortion_types"],
                        space_group_type=case["space_group_type"],
                        supercell_basis=case["supercell_basis"],
                        direct_sublattice_centering=case["direct_sublattice_centering"],
                        lattice_type=case["lattice_type"],
                        generate_if_missing=case["category"] == "parameter-k",
                    )
                candidates = [_candidate_record(item) for item in items]
                matched = any(
                    row["irrep"] == expected["irrep"]
                    and row["opd"] == expected["opd"]
                    and row["sg"] == expected["sg"]
                    for row in candidates
                )
                if not matched:
                    raise AssertionError(
                        f"expected {expected['irrep']} {expected['opd']} SG {expected['sg']} not returned"
                    )
                result = {
                    "status": "pass",
                    "local_subset_candidate_count": len(candidates),
                    "expected_candidate_found": True,
                    "candidates": candidates,
                }
            except Exception as exc:  # noqa: BLE001 - checkpoint the scientific failure
                result = {"status": "fail", "error": f"{type(exc).__name__}: {exc}"}
            results[case_id] = result
            done += 1
            checkpoint["updated_at"] = datetime.now(UTC).isoformat()
            _atomic_json(CHECKPOINT_PATH, checkpoint)
            print(f"Method 3 preflight {done}/{total}: {case_id} {result['status']}", flush=True)

    failures = [case_id for case_id, result in results.items() if result["status"] != "pass"]
    manifest["preflight"] = {
        "status": "pass" if not failures and len(results) == total else "fail",
        "completed": len(results),
        "total": total,
        "failures": failures,
        "results": results,
        "finished_at": datetime.now(UTC).isoformat(),
    }
    _atomic_json(MANIFEST_PATH, manifest)
    return 0 if manifest["preflight"]["status"] == "pass" else 1


def main() -> int:
    _validate_case_matrix()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-only",
        action="store_true",
        help="write the frozen manifest without invoking WSL/iso",
    )
    args = parser.parse_args()
    manifest = _base_manifest(_signature())
    if args.write_only:
        _atomic_json(MANIFEST_PATH, manifest)
        print(MANIFEST_PATH)
        return 0
    return _run_preflight(manifest)


if __name__ == "__main__":
    raise SystemExit(main())
