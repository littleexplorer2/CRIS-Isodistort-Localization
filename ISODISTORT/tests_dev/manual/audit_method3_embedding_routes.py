"""Audit single-IR routes for saved official Method 3 embeddings.

This validation tool is intentionally data-driven and read-only with respect
to the manifest and ``output_compare``.  For every embedding parsed from a
saved official Method 3 result page it asks the bundled ISOTROPY ``iso``
binary to ``DISPLAY DIRECTION`` for the exact parent/subgroup/basis/origin.

An embedding has an exact single-IR route when the returned direction table
contains the requested subgroup type with the same primitive-cell size ``s``.
Because ``DISPLAY DIRECTION`` is conditioned on the supplied basis and origin,
equal subgroup type and size also imply equal subgroup order (and therefore
the selected embedded subgroup rather than a proper supergroup).  When the
command succeeds but no such row exists, the embedding requires an
intersection of two or more single-IR stabilizers; recovering the actual
coupled route is a separate ``DISPLAY ISOTROPY COUPLED``/fixed-subspace task.

No crystal, case, basis, origin, IR, or expected classification is embedded in
this script.  Inputs come exclusively from the manifest and downloaded
official artifacts.  By default only authoritative result-table HTML is used.
The explicit ``--allow-candidate-inventory`` mode can audit identities recovered
from downloaded ``subgroup.cif`` files, but marks them provisional because the
downloaded set alone cannot prove that no website row was omitted.

Examples::

    python tests_dev/manual/audit_method3_embedding_routes.py
    python tests_dev/manual/audit_method3_embedding_routes.py \
        --case M3-EU-04 --json-output output/validation/route-audit.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
from collections.abc import Callable
from dataclasses import asdict
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from isocore.backend import IsoWrapper  # noqa: E402

try:  # Package import in pytest; direct import when this file is executed.
    from .audit_method3_downloads import (
        SubgroupIdentity,
        audit_case,
        manifest_case_expectation,
        parse_result_page,
    )
except ImportError:  # pragma: no cover - exercised by direct CLI invocation
    from audit_method3_downloads import (
        SubgroupIdentity,
        audit_case,
        manifest_case_expectation,
        parse_result_page,
    )

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method3_download_manifest.json"
OUTPUT_COMPARE_ROOT = WORKSPACE / "output_compare"
REPORT_PATH = PACKAGE_ROOT / "output" / "validation" / "method3_embedding_route_audit.json"
SCHEMA = 2


def _identity_from_json(item: dict[str, Any]) -> SubgroupIdentity:
    return SubgroupIdentity(
        space_group_number=int(item["space_group_number"]),
        symbol=str(item["symbol"]),
        basis=tuple(tuple(Fraction(value) for value in row) for row in item["basis"]),
        origin=tuple(Fraction(value) for value in item["origin"]),
        sublattice_index=int(item["s"]),
        group_index=int(item["i"]),
    )

_PARENT_SG_RE = re.compile(r"Space\s+Group\s*:\s*(\d+)\b", re.IGNORECASE)
_DIRECTION_HEADER_RE = re.compile(
    r"Irrep\s+\(ML\).*?k\s+vector.*?Dir.*?Subgroup.*?Size",
    re.IGNORECASE,
)
_DIRECTION_ROW_RE = re.compile(
    r"^\s*\*?\s*"
    r"(?P<irrep>\S+)\s+"
    r"(?:(?P<k_parameters>(?!\()[^\s]+)\s+)?"
    r"(?P<k_vector>\([^)]*\))\s+"
    r"(?P<direction>\([^)]*\))\s+"
    r"(?P<space_group_number>\d+)\s+"
    r"(?P<symbol>\S+)\s+"
    r"(?P<size>\d+)\s*$"
)


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_text(path: Path) -> str:
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "cp1252"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _vector_text(values: tuple[Fraction, ...]) -> str:
    return ",".join(_fraction_text(value) for value in values)


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _find_result_pages(case_dir: Path) -> list[tuple[Path, list[SubgroupIdentity]]]:
    """Return only genuine result pages; search/details pages are ignored."""

    result: list[tuple[Path, list[SubgroupIdentity]]] = []
    if not case_dir.is_dir():
        return result
    for path in sorted(case_dir.glob("*.htm*")):
        try:
            identities = parse_result_page(path)
        except ValueError:
            continue
        result.append((path, identities))
    return result


def parse_parent_space_group(path: Path) -> int:
    """Read the parent SG number printed in an official Method 3 result page."""

    match = _PARENT_SG_RE.search(_read_text(path))
    if match is None:
        raise ValueError("official result page does not contain a parent 'Space Group:' header")
    return int(match.group(1))


def build_display_direction_session(parent_sg: int, identity: SubgroupIdentity) -> list[str]:
    """Build the exact, non-interactive ISO command stream for one embedding."""

    basis = " ".join(_vector_text(row) for row in identity.basis)
    origin = _vector_text(identity.origin)
    return [
        "PAGE NOBREAK",
        "SCREEN 300",
        f"VALUE PARENT {parent_sg}",
        f"VALUE SUBGROUP {identity.space_group_number}",
        f"VALUE BASIS {basis}",
        f"VALUE ORIGIN {origin}",
        "SHOW KPOINT",
        "SHOW SUBGROUP",
        "SHOW SIZE",
        "DISPLAY DIRECTION",
        "QUIT",
    ]


def parse_display_direction_output(stdout: str) -> tuple[list[dict[str, Any]], list[str]]:
    """Parse ISO's fixed-width ``DISPLAY DIRECTION`` table.

    The parser does not infer routes from prose.  A missing table header is an
    operational/unsupported result and is therefore reported as indeterminate,
    never as evidence for a coupled-only embedding.
    """

    lines = stdout.splitlines()
    header_seen = any(_DIRECTION_HEADER_RE.search(line) for line in lines)
    errors: list[str] = []
    if not header_seen:
        diagnostic = next(
            (
                line.strip(" *")
                for line in lines
                if "error" in line.casefold()
                or "not been implemented" in line.casefold()
                or "not been selected" in line.casefold()
            ),
            "DISPLAY DIRECTION table header was not found",
        )
        errors.append(diagnostic)
        return [], errors

    routes: list[dict[str, Any]] = []
    for line in lines:
        match = _DIRECTION_ROW_RE.fullmatch(line)
        if match is None:
            continue
        routes.append(
            {
                "irrep": match.group("irrep"),
                "k_parameters": match.group("k_parameters") or "",
                "k_vector": match.group("k_vector"),
                "direction": match.group("direction"),
                "space_group_number": int(match.group("space_group_number")),
                "symbol": match.group("symbol"),
                "size": int(match.group("size")),
                "raw": line.strip(),
            }
        )
    return routes, errors


def classify_embedding_routes(
    identity: SubgroupIdentity,
    routes: list[dict[str, Any]],
    parse_errors: list[str],
) -> tuple[str, list[dict[str, Any]]]:
    """Classify one exact embedding without consulting an expected-answer list."""

    if parse_errors:
        return "indeterminate", []
    exact = [
        route
        for route in routes
        if route["space_group_number"] == identity.space_group_number
        and route["size"] == identity.sublattice_index
    ]
    return ("single_ir_exact" if exact else "coupled_ir_required"), exact


def audit_embedding(
    parent_sg: int,
    identity: SubgroupIdentity,
    run_iso: Callable[[str], str],
    *,
    include_raw_output: bool = True,
) -> dict[str, Any]:
    """Run and record one exact DISPLAY DIRECTION query."""

    commands = build_display_direction_session(parent_sg, identity)
    command_text = "\n".join(commands) + "\n"
    started = time.monotonic()
    try:
        stdout = run_iso(command_text)
    except Exception as exc:  # noqa: BLE001 - report each external-tool failure and continue
        return {
            "official_identity": identity.to_json(),
            "iso_input": {
                "parent_space_group_number": parent_sg,
                "subgroup_space_group_number": identity.space_group_number,
                "basis": [[_fraction_text(value) for value in row] for row in identity.basis],
                "origin": [_fraction_text(value) for value in identity.origin],
                "target_size": identity.sublattice_index,
                "target_index": identity.group_index,
            },
            "commands": commands,
            "classification": "indeterminate",
            "exact_single_ir_routes": [],
            "display_direction_rows": [],
            "errors": [f"{type(exc).__name__}: {exc}"],
            "elapsed_seconds": time.monotonic() - started,
        }

    routes, errors = parse_display_direction_output(stdout)
    classification, exact = classify_embedding_routes(identity, routes, errors)
    record: dict[str, Any] = {
        "official_identity": identity.to_json(),
        "iso_input": {
            "parent_space_group_number": parent_sg,
            "subgroup_space_group_number": identity.space_group_number,
            "basis": [[_fraction_text(value) for value in row] for row in identity.basis],
            "origin": [_fraction_text(value) for value in identity.origin],
            "target_size": identity.sublattice_index,
            "target_index": identity.group_index,
        },
        "commands": commands,
        "classification": classification,
        "exact_single_ir_routes": exact,
        "display_direction_rows": routes,
        "errors": errors,
        "elapsed_seconds": time.monotonic() - started,
        "stdout_sha256": hashlib.sha256(stdout.encode("utf-8")).hexdigest(),
    }
    if include_raw_output:
        record["raw_stdout"] = stdout
    return record


def _manifest_cases(payload: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    cases: list[tuple[str, dict[str, Any]]] = []
    for parent in payload.get("parents", []):
        parent_cif = str(parent.get("parent_cif", "")).strip()
        for case in parent.get("cases", []):
            cases.append((parent_cif, case))
    return cases


def _iso_provenance(wrapper: IsoWrapper, version_stdout: str) -> dict[str, Any]:
    version_match = re.search(r"^Isotropy,\s*Version\s+([^\r\n]+)", version_stdout, re.MULTILINE)
    binary = Path(wrapper.binary).resolve()
    data_dir = binary.parent
    data_files = sorted(
        path
        for path in data_dir.iterdir()
        if path.is_file() and (path.name == "const.dat" or path.name.startswith("data_"))
    )
    return {
        "version": version_match.group(1).strip() if version_match else None,
        "version_stdout": version_stdout,
        "binary": str(binary),
        "binary_sha256": _sha256(binary),
        "data_directory": str(data_dir),
        "data_file_sha256": {path.name: _sha256(path) for path in data_files},
    }


def run_audit(
    manifest_path: Path,
    output_compare_root: Path,
    *,
    parent_filter: str | None = None,
    case_filter: set[str] | None = None,
    timeout: float | None = None,
    include_raw_output: bool = True,
    allow_candidate_inventory: bool = False,
) -> dict[str, Any]:
    """Audit every eligible official row selected by the manifest."""

    manifest_path = manifest_path.resolve()
    output_compare_root = output_compare_root.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    selected = [
        (parent_cif, case)
        for parent_cif, case in _manifest_cases(manifest)
        if (parent_filter is None or parent_cif == parent_filter)
        and (not case_filter or str(case.get("id")) in case_filter)
    ]
    parent_space_groups = {
        str(parent.get("parent_cif", "")).strip(): int(
            parent.get("parent_space_group_type")
        )
        for parent in manifest.get("parents", [])
    }

    wrapper: IsoWrapper | None = None
    iso_error: str | None = None
    version_stdout = ""
    try:
        wrapper = IsoWrapper()
        version_stdout = wrapper.run_stdin(wrapper.binary, "QUIT\n", timeout=timeout)
    except Exception as exc:  # noqa: BLE001 - retain a machine-readable environment failure
        iso_error = f"{type(exc).__name__}: {exc}"

    provenance: dict[str, Any] = {
        "manifest": {
            "path": str(manifest_path),
            "sha256": _sha256(manifest_path),
            "declared_signature": manifest.get("signature"),
        },
        "output_compare_root": str(output_compare_root),
        "classification_rule": (
            "single_ir_exact iff DISPLAY DIRECTION succeeds and contains a row whose subgroup "
            "space-group number and primitive-cell Size equal the official embedding's SG and s; "
            "otherwise a successful table is coupled_ir_required. Missing/unsupported output is indeterminate."
        ),
        "scientific_scope": (
            "DISPLAY DIRECTION tests exact single-IR stabilizers for a supplied embedding. "
            "It does not enumerate the actual multi-IR intersection; that requires DISPLAY ISOTROPY COUPLED "
            "or an equivalent fixed-subspace algorithm."
        ),
    }
    if wrapper is not None and not iso_error:
        provenance["iso"] = _iso_provenance(wrapper, version_stdout)
    else:
        provenance["iso"] = {"error": iso_error}

    case_records: list[dict[str, Any]] = []
    for parent_cif, case in selected:
        case_id = str(case.get("id", ""))
        folder = str(case.get("readable_case_folder", ""))
        case_dir = output_compare_root / parent_cif / "官网" / "Method3" / folder
        record: dict[str, Any] = {
            "parent_cif": parent_cif,
            "case_id": case_id,
            "readable_case_folder": folder,
            "official_case_directory": str(case_dir),
            "manifest_space_group_type": case.get("space_group_type"),
            "status": "pending",
            "issues": [],
            "embeddings": [],
        }
        parent_sg = parent_space_groups[parent_cif]
        expectation = manifest_case_expectation(
            parent_cif=parent_cif,
            parent_space_group_number=parent_sg,
            case=case,
        )
        inventory_audit = audit_case(case_dir, expectation)
        record["official_download_audit_status"] = inventory_audit.status
        record["official_download_audit_issues"] = [
            asdict(issue) for issue in inventory_audit.issues
        ] + [
            asdict(issue)
            for candidate in inventory_audit.candidates
            for issue in candidate.issues
        ]
        pages = [Path(path) for path in inventory_audit.result_pages]
        page_errors = [
            issue
            for issue in inventory_audit.issues
            if issue.severity == "error"
            and (
                issue.code in {"invalid_result_context", "invalid_result_row"}
                or issue.code.startswith("result_")
            )
        ]
        if len(pages) == 1 and page_errors:
            record["status"] = "error"
            record["issues"].append(
                "saved result page is invalid or does not match the manifest query context"
            )
            case_records.append(record)
            continue
        if len(pages) != 1:
            identities = [
                _identity_from_json(item)
                for item in inventory_audit.downloaded_candidates
            ]
            duplicate_count = len(identities) - len(set(identities))
            invalid_candidate_count = sum(
                issue.code == "invalid_subgroup_cif"
                for candidate in inventory_audit.candidates
                for issue in candidate.issues
            )
            wrong_types = sorted(
                {
                    identity.space_group_number
                    for identity in identities
                    if identity.space_group_number != int(case.get("space_group_type"))
                }
            )
            if (
                not allow_candidate_inventory
                or len(pages) > 1
                or not identities
                or duplicate_count
                or invalid_candidate_count
                or wrong_types
            ):
                record["status"] = "skipped_no_unique_official_result"
                record["issues"].append(
                    f"expected exactly one genuine result page, found {len(pages)}"
                )
                if duplicate_count:
                    record["issues"].append(
                        f"candidate inventory contains {duplicate_count} duplicate identities"
                    )
                if invalid_candidate_count:
                    record["issues"].append(
                        "candidate inventory contains "
                        f"{invalid_candidate_count} invalid subgroup.cif identities"
                    )
                if wrong_types:
                    record["issues"].append(
                        f"candidate inventory contains wrong target SG types: {wrong_types}"
                    )
                record["inventory_candidate_count"] = len(identities)
                case_records.append(record)
                continue
            record["official_source"] = "downloaded_candidate_inventory"
            record["official_completeness_proven"] = False
            record["issues"].append(
                "provisional: result-table HTML is missing; downloaded candidates do not prove completeness"
            )
        else:
            result_page = pages[0]
            try:
                identities = parse_result_page(result_page)
                parsed_parent_sg = parse_parent_space_group(result_page)
            except (OSError, UnicodeError, ValueError) as exc:
                record["status"] = "error"
                record["issues"].append(
                    f"invalid official result page: {type(exc).__name__}: {exc}"
                )
                case_records.append(record)
                continue
            record["official_result_page"] = str(result_page)
            record["official_result_page_sha256"] = _sha256(result_page)
            record["official_source"] = "result_table_html"
            record["official_completeness_proven"] = True
            if parsed_parent_sg != parent_sg:
                record["status"] = "error"
                record["issues"].append(
                    f"result page parent SG {parsed_parent_sg} does not match manifest {parent_sg}"
                )
                case_records.append(record)
                continue
        record["official_embedding_count"] = len(identities)
        record["parent_space_group_number"] = parent_sg

        if not identities:
            record["status"] = "verified_empty_official_result"
            case_records.append(record)
            continue
        if wrapper is None or iso_error:
            record["status"] = "error"
            record["issues"].append(iso_error or "ISO wrapper is unavailable")
            case_records.append(record)
            continue

        def run_iso(commands: str) -> str:
            return wrapper.run_stdin(wrapper.binary, commands, timeout=timeout)

        record["embeddings"] = [
            audit_embedding(
                parent_sg,
                identity,
                run_iso,
                include_raw_output=include_raw_output,
            )
            for identity in identities
        ]
        if any(item["classification"] == "indeterminate" for item in record["embeddings"]):
            record["status"] = "error"
        else:
            record["status"] = (
                "audited"
                if record.get("official_completeness_proven")
                else "audited_provisional"
            )
        case_records.append(record)

    embeddings = [embedding for case in case_records for embedding in case["embeddings"]]
    authoritative_embeddings = [
        embedding
        for case in case_records
        if case.get("official_completeness_proven") is True
        for embedding in case["embeddings"]
    ]
    provisional_embeddings = [
        embedding
        for case in case_records
        if case.get("official_source") == "downloaded_candidate_inventory"
        for embedding in case["embeddings"]
    ]
    counts = {
        classification: sum(item["classification"] == classification for item in embeddings)
        for classification in ("single_ir_exact", "coupled_ir_required", "indeterminate")
    }
    authoritative_counts = {
        f"authoritative_{classification}": sum(
            item["classification"] == classification
            for item in authoritative_embeddings
        )
        for classification in ("single_ir_exact", "coupled_ir_required", "indeterminate")
    }
    provisional_counts = {
        f"provisional_{classification}": sum(
            item["classification"] == classification
            for item in provisional_embeddings
        )
        for classification in ("single_ir_exact", "coupled_ir_required", "indeterminate")
    }
    return {
        "schema": SCHEMA,
        "generated_at": _now(),
        "tool": "tests_dev/manual/audit_method3_embedding_routes.py",
        "provenance": provenance,
        "filters": {
            "parent": parent_filter,
            "cases": sorted(case_filter or []),
            "raw_iso_output_included": include_raw_output,
            "allow_candidate_inventory": allow_candidate_inventory,
        },
        "summary": {
            "selected_case_count": len(selected),
            "audited_case_count": sum(case["status"] == "audited" for case in case_records),
            "provisional_audited_case_count": sum(
                case["status"] == "audited_provisional" for case in case_records
            ),
            "verified_empty_case_count": sum(
                case["status"] == "verified_empty_official_result" for case in case_records
            ),
            "skipped_case_count": sum(
                case["status"] == "skipped_no_unique_official_result" for case in case_records
            ),
            "error_case_count": sum(case["status"] == "error" for case in case_records),
            "official_embedding_count": len(authoritative_embeddings),
            "provisional_inventory_embedding_count": len(provisional_embeddings),
            "audited_embedding_count": len(embeddings),
            **counts,
            **authoritative_counts,
            **provisional_counts,
        },
        "cases": case_records,
    }


def _print_summary(report: dict[str, Any]) -> None:
    summary = report["summary"]
    print(
        "Method 3 embedding-route audit: "
        f"{summary['audited_case_count']} authoritative + "
        f"{summary['provisional_audited_case_count']} provisional cases / "
        f"{summary['audited_embedding_count']} embeddings audited "
        f"({summary['official_embedding_count']} authoritative + "
        f"{summary['provisional_inventory_embedding_count']} provisional); "
        f"single-IR={summary['single_ir_exact']}, coupled-required={summary['coupled_ir_required']}, "
        f"indeterminate={summary['indeterminate']}; "
        f"skipped cases={summary['skipped_case_count']}, error cases={summary['error_case_count']}"
    )
    for case in report["cases"]:
        classifications: dict[str, int] = {}
        for embedding in case["embeddings"]:
            key = embedding["classification"]
            classifications[key] = classifications.get(key, 0) + 1
        suffix = ", ".join(f"{key}={value}" for key, value in sorted(classifications.items()))
        if suffix:
            suffix = "; " + suffix
        print(f"[{case['status']}] {case['case_id']}{suffix}")


def _audit_exit_code(
    report: dict[str, Any], *, accept_provisional: bool = False
) -> int:
    """Return 0 only for complete, accepted route evidence."""

    summary = report["summary"]
    if summary["error_case_count"]:
        return 1
    if summary["selected_case_count"] == 0 or summary["skipped_case_count"]:
        return 2
    if summary["provisional_audited_case_count"] and not accept_provisional:
        return 2
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Read-only ISO DISPLAY DIRECTION audit of saved official Method 3 embeddings."
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output-compare-root", type=Path, default=OUTPUT_COMPARE_ROOT)
    parser.add_argument("--json-output", type=Path, default=REPORT_PATH)
    parser.add_argument("--parent", help="exact parent CIF name from the manifest")
    parser.add_argument("--case", action="append", dest="cases", help="case ID; repeat to select several")
    parser.add_argument("--timeout", type=float, default=None, help="per-ISO-call timeout in seconds")
    parser.add_argument(
        "--omit-raw-output",
        action="store_true",
        help="retain parsed rows and stdout hashes but omit raw ISO stdout from JSON",
    )
    parser.add_argument(
        "--allow-candidate-inventory",
        action="store_true",
        help=(
            "provisionally audit downloaded subgroup.cif identities when the authoritative "
            "result-table HTML is missing; this does not prove completeness"
        ),
    )
    parser.add_argument(
        "--accept-provisional",
        action="store_true",
        help=(
            "allow otherwise successful provisional candidate-inventory audits to exit 0; "
            "errors and genuinely skipped cases still exit nonzero"
        ),
    )
    args = parser.parse_args(argv)

    report = run_audit(
        args.manifest,
        args.output_compare_root,
        parent_filter=args.parent,
        case_filter=set(args.cases or []),
        timeout=args.timeout,
        include_raw_output=not args.omit_raw_output,
        allow_candidate_inventory=args.allow_candidate_inventory,
    )
    _atomic_json(args.json_output, report)
    _print_summary(report)
    print(f"JSON: {args.json_output.resolve()}")
    return _audit_exit_code(report, accept_provisional=args.accept_provisional)


if __name__ == "__main__":
    raise SystemExit(main())
