"""Compare local Method 3 embeddings with saved official result tables.

This is a resumable, read-only validation runner.  It obtains every query from
``method3_download_manifest.json`` and every authoritative reference embedding
from the saved official result-table HTML; it contains no crystal-specific
answer table.  Cases without exactly one valid result page are reported but
are not run locally by default.  An explicit ``--allow-candidate-inventory``
mode can additionally compare the identities recovered from downloaded
``subgroup.cif`` files.  Such cases remain marked provisional because a set of
downloaded candidates cannot prove that no website result row was omitted.

Examples::

    python tests_dev/manual/compare_method3_official_local.py
    python tests_dev/manual/compare_method3_official_local.py --parent "EuAl4 Parent.cif"
    python tests_dev/manual/compare_method3_official_local.py --case M3-EU-04 --restart

The report keeps two deliberately separate verdicts.  A literal comparison
checks the website's printed ``(SG, basis, origin, s, i)`` fields.  An exact
affine comparison then reconstructs subgroup operations and proves equality
or parent-space-group conjugacy modulo the child primitive translation
lattice.  Thus a normalized equivalence is never mislabeled as a literal
website representation match.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import os
import platform
import sys
import time
from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from copy import deepcopy
from dataclasses import asdict
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

import yaml

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from isocore.api import IsoDistort  # noqa: E402

try:  # Package import in pytest; direct import when this file is executed.
    from .audit_method3_downloads import (
        SubgroupIdentity,
        audit_case,
        manifest_case_expectation,
        manifest_parent_space_group_number,
        parse_result_page,
    )
    from .method3_affine_equivalence import compare_affine_embedding_sets
except ImportError:  # pragma: no cover - exercised by direct CLI invocation
    from audit_method3_downloads import (
        SubgroupIdentity,
        audit_case,
        manifest_case_expectation,
        manifest_parent_space_group_number,
        parse_result_page,
    )
    from method3_affine_equivalence import compare_affine_embedding_sets

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method3_download_manifest.json"
OUTPUT_COMPARE_ROOT = WORKSPACE / "output_compare"
VALIDATION_DIR = PACKAGE_ROOT / "output" / "validation"
REPORT_PATH = VALIDATION_DIR / "method3_official_local_comparison.json"
CHECKPOINT_PATH = VALIDATION_DIR / "method3_official_local_checkpoint.json"
INPUT_ROOT = WORKSPACE / "experiment_data"
CHECKPOINT_SCHEMA = 3
ISO_DATA_BASENAMES = (
    "const.dat",
    "data_diperiodic.txt",
    "data_images.txt",
    "data_irreps.txt",
    "data_isotropy.txt",
    "data_little.txt",
    "data_magnetic.txt",
    "data_space.txt",
    "data_ssg.txt",
    "data_ssgmag.txt",
    "data_wyckoff.txt",
)
RUNTIME_DISTRIBUTIONS = ("numpy", "pymatgen", "spglib", "PyYAML")


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _as_fraction(value: object) -> Fraction:
    """Preserve exact integer/fraction strings and stabilize float display values."""

    if isinstance(value, Fraction):
        return value
    if isinstance(value, bool):
        return Fraction(int(value))
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, float):
        return Fraction(str(value)).limit_denominator(1_000_000)
    return Fraction(str(value).strip())


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _json_value(value: object) -> object:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, Iterable):
        return [_json_value(item) for item in value]
    return repr(value)


def _identity_from_local(item: Any) -> SubgroupIdentity:
    subgroup = item.subgroup
    basis_source = getattr(item, "basis", None) or subgroup.basis_vectors
    basis = tuple(tuple(_as_fraction(value) for value in row) for row in basis_source)
    if len(basis) != 3 or any(len(row) != 3 for row in basis):
        raise ValueError(f"local Method 3 basis is not 3x3: {basis_source!r}")
    origin = tuple(_as_fraction(value) for value in (subgroup.origin or [0, 0, 0]))
    if len(origin) != 3:
        raise ValueError(f"local Method 3 origin does not contain three values: {subgroup.origin!r}")
    return SubgroupIdentity(
        space_group_number=int(subgroup.space_group_number),
        symbol=" ".join(str(subgroup.space_group_symbol).split()),
        basis=basis,  # type: ignore[arg-type]
        origin=origin,  # type: ignore[arg-type]
        sublattice_index=int(subgroup.size),
        group_index=int(subgroup.subgroup_index),
    )


def _route_record(route: Any) -> dict[str, Any]:
    fields: dict[str, Any] = {}
    official_fields = getattr(route, "official_fields", None)
    if callable(official_fields):
        try:
            fields = dict(official_fields())
        except (AttributeError, TypeError, ValueError):
            fields = {}
    return {
        "irrep": getattr(route, "irrep_label", ""),
        "opd": getattr(route, "opd_symbol", ""),
        "k_point": getattr(route, "k_point_label", ""),
        "k_parameters": _json_value(getattr(route, "k_parameters", [])),
        "k_active": fields.get("k_active", getattr(route, "k_active_raw", "")),
    }


def _local_record(item: Any, identity: SubgroupIdentity) -> dict[str, Any]:
    subgroup = item.subgroup
    item_routes = getattr(item, "routes", None)
    routes = item_routes if item_routes is not None else [subgroup]
    return {
        "identity": identity.to_json(),
        "raw_fields": {
            "space_group_number": getattr(subgroup, "space_group_number", None),
            "space_group_symbol": getattr(subgroup, "space_group_symbol", ""),
            "basis": _json_value(getattr(item, "basis", None)),
            "subgroup_basis_vectors": _json_value(getattr(subgroup, "basis_vectors", [])),
            "basis_raw": getattr(subgroup, "basis_raw", ""),
            "origin": _json_value(getattr(subgroup, "origin", [])),
            "origin_raw": getattr(subgroup, "origin_raw", ""),
            "s": getattr(subgroup, "size", None),
            "i": getattr(subgroup, "subgroup_index", None),
        },
        "routes": [_route_record(route) for route in routes],
    }


def _counter_json(counter: Counter[SubgroupIdentity]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for identity, multiplicity in sorted(counter.items(), key=lambda pair: pair[0].text):
        record = identity.to_json()
        record["multiplicity"] = multiplicity
        records.append(record)
    return records


def _sg_counts(identities: Iterable[SubgroupIdentity]) -> dict[str, int]:
    counts = Counter(identity.space_group_number for identity in identities)
    return {str(number): count for number, count in sorted(counts.items())}


def compare_embedding_sets(
    official: Sequence[SubgroupIdentity],
    local_items: Sequence[Any],
    *,
    parent_structure: Any | None = None,
    symmetry_tolerance: float | None = None,
) -> dict[str, Any]:
    """Compare literal fields and, when possible, exact affine embeddings."""

    local_identities: list[SubgroupIdentity] = []
    local_records: list[dict[str, Any]] = []
    conversion_errors: list[dict[str, Any]] = []
    for index, item in enumerate(local_items):
        try:
            identity = _identity_from_local(item)
        except (AttributeError, TypeError, ValueError, ZeroDivisionError) as exc:
            conversion_errors.append(
                {
                    "index": index,
                    "error": f"{type(exc).__name__}: {exc}",
                    "raw_item": repr(item),
                }
            )
            continue
        local_identities.append(identity)
        local_records.append(_local_record(item, identity))

    official_counter = Counter(official)
    local_counter = Counter(local_identities)
    missing = official_counter - local_counter
    unexpected = local_counter - official_counter
    official_sg_counts = _sg_counts(official)
    local_sg_counts = _sg_counts(local_identities)
    count_match = len(official) == len(local_items)
    exact_match = not conversion_errors and not missing and not unexpected and count_match
    result: dict[str, Any] = {
        "status": "exact_match" if exact_match else "different",
        "comparison_semantics": (
            "Literal comparison of SG number, normalized symbol whitespace, 3x3 basis, "
            "origin, s and i. Alternative settings are not treated as equal."
        ),
        "official_count": len(official),
        "local_embedding_count": len(local_items),
        "count_match": count_match,
        "official_space_group_counts": official_sg_counts,
        "local_space_group_counts": local_sg_counts,
        "space_group_multiset_match": official_sg_counts == local_sg_counts,
        "exact_representation_match": exact_match,
        "official_embeddings": [identity.to_json() for identity in official],
        "local_embeddings": local_records,
        "missing_exact_official_embeddings": _counter_json(missing),
        "unexpected_exact_local_embeddings": _counter_json(unexpected),
        "local_conversion_errors": conversion_errors,
        "interpretation_note": (
            "Literal and crystallographic verdicts are separate. The affine_equivalence "
            "record, when present, is the proof-level normalized comparison."
        ),
    }
    if parent_structure is not None:
        tolerance = float(symmetry_tolerance if symmetry_tolerance is not None else 1e-3)
        affine = compare_affine_embedding_sets(
            official,
            local_identities,
            parent_structure,
            symprec=tolerance,
        )
        result["affine_equivalence"] = affine
        result["missing_affine_official_embeddings"] = [
            official[index].to_json()
            for index in affine["unmatched_official_indices"]
        ]
        result["unexpected_affine_local_embeddings"] = [
            local_records[index]
            for index in affine["unmatched_local_indices"]
        ]
        if affine["complete_match"]:
            result["status"] = "exact_match" if exact_match else "affine_equivalent"
    return result


def _manifest_cases(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for parent in manifest.get("parents", []):
        parent_name = parent.get("parent_cif")
        if not isinstance(parent_name, str) or not parent_name:
            raise ValueError("each Method 3 manifest parent requires a non-empty parent_cif")
        parent_space_group_number = manifest_parent_space_group_number(parent)
        for case in parent.get("cases", []):
            if not isinstance(case, dict) or not case.get("id") or not case.get("readable_case_folder"):
                raise ValueError(f"invalid Method 3 case beneath {parent_name!r}")
            result.append(
                {
                    "parent_cif": parent_name,
                    "parent_space_group_number": parent_space_group_number,
                    "case": case,
                }
            )
    return result


def select_manifest_cases(
    manifest: dict[str, Any],
    parents: Sequence[str] = (),
    case_ids: Sequence[str] = (),
) -> list[dict[str, Any]]:
    """Select manifest cases using case-insensitive exact parent/case names."""

    parent_filter = {item.casefold() for item in parents}
    case_filter = {item.casefold() for item in case_ids}
    selected = [
        spec
        for spec in _manifest_cases(manifest)
        if (not parent_filter or spec["parent_cif"].casefold() in parent_filter)
        and (not case_filter or str(spec["case"]["id"]).casefold() in case_filter)
    ]
    if parents:
        known = {spec["parent_cif"].casefold() for spec in _manifest_cases(manifest)}
        unknown = sorted(item for item in parents if item.casefold() not in known)
        if unknown:
            raise ValueError("unknown parent filter(s): " + ", ".join(unknown))
    if case_ids:
        known = {str(spec["case"]["id"]).casefold() for spec in _manifest_cases(manifest)}
        unknown = sorted(item for item in case_ids if item.casefold() not in known)
        if unknown:
            raise ValueError("unknown case filter(s): " + ", ".join(unknown))
    return selected


def _case_root(
    template: str,
    output_compare_root: Path,
    parent_cif: str,
    readable_case_folder: str,
) -> Path:
    expanded = (
        template.replace("\\", "/")
        .replace("<parent>", parent_cif)
        .replace("<readable-case-folder>", readable_case_folder)
    )
    parts = [part for part in expanded.split("/") if part and part != "."]
    if any("<" in part or ">" in part for part in parts):
        raise ValueError(f"unresolved placeholder in official_download_root: {template!r}")
    if parts and parts[0].casefold() == output_compare_root.name.casefold():
        parts = parts[1:]
    return output_compare_root.joinpath(*parts)


def _official_input(
    case_dir: Path,
    *,
    parent_cif: str,
    parent_space_group_number: int,
    manifest_case: Mapping[str, Any],
    allow_candidate_inventory: bool = False,
) -> dict[str, Any]:
    expectation = manifest_case_expectation(
        parent_cif=parent_cif,
        parent_space_group_number=parent_space_group_number,
        case=manifest_case,
    )
    audit = audit_case(case_dir, expectation)
    if len(audit.result_pages) != 1:
        reason = "missing_result_page" if not audit.result_pages else "ambiguous_multiple_result_pages"
        inventory_identities = [
            SubgroupIdentity(
                space_group_number=int(item["space_group_number"]),
                symbol=str(item["symbol"]),
                basis=tuple(
                    tuple(Fraction(value) for value in row)
                    for row in item["basis"]
                ),
                origin=tuple(Fraction(value) for value in item["origin"]),
                sublattice_index=int(item["s"]),
                group_index=int(item["i"]),
            )
            for item in audit.downloaded_candidates
        ]
        target_sg = int(manifest_case["space_group_type"])
        inventory_errors = [
            asdict(issue)
            for issue in audit.issues
            if issue.code == "duplicate_downloaded_identity"
        ] + [
            asdict(issue)
            for candidate in audit.candidates
            for issue in candidate.issues
            if issue.code == "invalid_subgroup_cif"
        ]
        wrong_types = sorted(
            {
                identity.space_group_number
                for identity in inventory_identities
                if identity.space_group_number != target_sg
            }
        )
        if wrong_types:
            inventory_errors.append(
                {
                    "severity": "error",
                    "code": "candidate_target_space_group_mismatch",
                    "message": (
                        f"downloaded candidate SG types {wrong_types} do not match "
                        f"manifest target SG {target_sg}"
                    ),
                    "path": str(case_dir),
                }
            )
        if allow_candidate_inventory and inventory_identities and not inventory_errors:
            return {
                "eligible": True,
                "case_directory": str(case_dir),
                "result_page": None,
                "result_page_path": None,
                "identities": inventory_identities,
                "download_audit_status": audit.status,
                "official_source": "downloaded_candidate_inventory",
                "official_completeness_proven": False,
                "provisional_reason": reason,
                "audit_issues": [asdict(issue) for issue in audit.issues]
                + [
                    asdict(issue)
                    for candidate in audit.candidates
                    for issue in candidate.issues
                ],
            }
        return {
            "eligible": False,
            "reason": reason,
            "case_directory": str(case_dir),
            "result_pages": audit.result_pages,
            "html_pages": audit.html_pages,
            "download_audit_status": audit.status,
            "inventory_candidate_count": len(inventory_identities),
            "inventory_errors": inventory_errors,
        }
    context_errors = [
        issue
        for issue in audit.issues
        if issue.code == "invalid_result_context" or issue.code.startswith("result_")
    ]
    if context_errors:
        return {
            "eligible": False,
            "reason": "result_context_mismatch",
            "case_directory": str(case_dir),
            "result_pages": audit.result_pages,
            "html_pages": audit.html_pages,
            "download_audit_status": audit.status,
            "context_errors": [asdict(issue) for issue in context_errors],
        }
    result_path = Path(audit.result_pages[0])
    try:
        identities = parse_result_page(result_path)
    except ValueError as exc:
        return {
            "eligible": False,
            "reason": "invalid_result_page",
            "case_directory": str(case_dir),
            "result_pages": audit.result_pages,
            "html_pages": audit.html_pages,
            "download_audit_status": audit.status,
            "parse_error": str(exc),
        }
    return {
        "eligible": True,
        "case_directory": str(case_dir),
        "result_page": str(result_path),
        "result_page_path": result_path,
        "identities": identities,
        "download_audit_status": audit.status,
        "official_source": "result_table_html",
        "official_completeness_proven": True,
    }


def _hash_file(digest: Any, path: Path, label: str) -> None:
    digest.update(label.encode("utf-8"))
    digest.update(b"\0")
    if path.is_file():
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
    else:
        digest.update(b"<missing>")
    digest.update(b"\0")


def _configured_path(config_path: Path, value: object, fallback: str) -> Path:
    raw = value if isinstance(value, str) and value.strip() else fallback
    path = Path(raw)
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def _iso_runtime_dependencies() -> list[tuple[str, Path]]:
    """Discover immutable backend inputs without executing the Linux binary.

    Method 3 results depend on the exact ``iso`` executable and its bundled
    ISODATA tables, not only on Python source.  The configured paths are safe
    to inspect directly on both Windows and Linux.  Fixed expected table names
    ensure a missing file contributes an explicit ``<missing>`` marker; extra
    future ``data_*.txt`` tables are discovered and signed as well.  Mutable
    generated parameter-k caches are handled separately by disabling unsafe
    cross-process reuse for the cases that can consume them.
    """

    config_path = PACKAGE_ROOT / "config" / "settings.yaml"
    isobyu: dict[str, Any] = {}
    try:
        loaded = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict) and isinstance(loaded.get("isobyu"), dict):
            isobyu = loaded["isobyu"]
    except (OSError, UnicodeError, yaml.YAMLError):
        # The config file itself is also signed.  Falling back keeps checkpoint
        # invalidation deterministic even when it is temporarily absent/bad.
        pass

    bin_dir = _configured_path(config_path, isobyu.get("bin_dir"), "../isobyu")
    data_dir = _configured_path(config_path, isobyu.get("data_dir"), "../isobyu")
    iso_name = isobyu.get("iso_bin")
    if not isinstance(iso_name, str) or not iso_name.strip():
        iso_name = "iso"

    dependencies = [("iso-binary", bin_dir / iso_name)]
    known_names = set(ISO_DATA_BASENAMES)
    try:
        known_names.update(path.name for path in data_dir.glob("data_*.txt"))
    except OSError:
        pass
    dependencies.extend(
        (f"iso-data:{name}", data_dir / name)
        for name in sorted(known_names, key=str.casefold)
    )
    return dependencies


def _runtime_dependency_versions() -> dict[str, Any]:
    """Return a stable description of the interpreter and scientific stack."""

    distributions: dict[str, str] = {}
    for name in RUNTIME_DISTRIBUTIONS:
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "<missing>"
        except Exception as exc:  # noqa: BLE001 - signature must survive broken metadata
            distributions[name] = f"<unavailable:{type(exc).__module__}.{type(exc).__qualname__}>"
    return {
        "python": {
            "implementation": sys.implementation.name,
            "version": list(sys.version_info[:5]),
            "cache_tag": sys.implementation.cache_tag,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "distributions": distributions,
    }


def _case_uses_mutable_isotropy_cache(case: Mapping[str, Any]) -> bool:
    """Whether a case may read or create generated ``i*.iso`` parameter-k data."""

    return (
        str(case.get("category", "")).casefold() == "parameter-k"
        or bool(case.get("generate_if_missing"))
    )


def _isotropy_cache_signature_policy(
    selected: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    """Describe the safe checkpoint policy for mutable generated ISO data.

    The public cache inventory exposes only name/size/second-resolution mtime
    and header, not a content SHA-256.  In addition, ``generate_if_missing`` can
    change that inventory *during this very batch*.  Treating the pre-run list
    as a batch signature would therefore invalidate an interrupted run because
    of its own work, while still failing to prove cache content identity.

    Parameter-k records are consequently never restored across invocations.
    This stable marker documents why no pseudo-fingerprint is claimed.  Cases
    that cannot consume this cache retain normal content-signed checkpointing.
    """

    affected = sorted(
        _case_key(str(spec.get("parent_cif", "")), str(spec.get("case", {}).get("id", "")))
        for spec in selected
        if isinstance(spec.get("case"), Mapping)
        and _case_uses_mutable_isotropy_cache(spec["case"])
    )
    if not affected:
        return {"status": "not_required", "affected_case_keys": []}
    return {
        "status": "content_fingerprint_unavailable",
        "failure_marker": "public_cache_api_has_no_content_sha256",
        "resume_policy": "disable_cross_process_reuse_for_affected_cases",
        "affected_case_keys": affected,
    }


def _run_signature(
    manifest_path: Path,
    input_root: Path,
    selected: Sequence[dict[str, Any]],
    official_inputs: dict[str, dict[str, Any]],
    *,
    allow_candidate_inventory: bool = False,
) -> str:
    """Invalidate checkpoints when any scientific or official input changes.

    The comparison calls through a broad slice of :mod:`isocore` (including
    smodes/type filtering and superspace helpers), and its official side is
    produced by ``audit_method3_downloads``.  Signing only the immediately
    imported modules or only the result-table HTML can therefore reuse a stale
    checkpoint after either the parser/filtering logic or a downloaded core
    file changes.  Hash the complete Python implementation surface plus every
    validation-relevant file in each official case directory instead.
    """

    digest = hashlib.sha256(usedforsecurity=False)
    digest.update(
        f"allow-candidate-inventory:{int(allow_candidate_inventory)}\0".encode()
    )
    runtime_state = {
        "runtime_dependencies": _runtime_dependency_versions(),
        "generated_isotropy_cache": _isotropy_cache_signature_policy(selected),
    }
    digest.update(
        json.dumps(
            runtime_state,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    digest.update(b"\0")
    sources = [
        Path(__file__),
        PACKAGE_ROOT / "tests_dev" / "manual" / "audit_method3_downloads.py",
        PACKAGE_ROOT / "tests_dev" / "manual" / "method3_affine_equivalence.py",
        manifest_path,
        PACKAGE_ROOT / "config" / "settings.yaml",
    ]
    sources.extend(
        sorted(
            (PACKAGE_ROOT / "isocore").rglob("*.py"),
            key=lambda path: path.as_posix().casefold(),
        )
    )
    for path in sources:
        _hash_file(digest, path, str(path))
    for label, path in _iso_runtime_dependencies():
        _hash_file(digest, path, label)
    for parent_name in sorted({spec["parent_cif"] for spec in selected}):
        _hash_file(digest, input_root / parent_name, f"parent:{parent_name}")
    for case_key in sorted(official_inputs):
        record = official_inputs[case_key]
        case_directory = Path(str(record.get("case_directory", "")))
        digest.update(f"official-case:{case_key}\0".encode())
        if not case_directory.is_dir():
            digest.update(b"<missing-directory>\0")
            continue

        dependencies = [
            path
            for path in case_directory.rglob("*")
            if path.is_file()
            and not any(parent.name.casefold().endswith("_files") for parent in path.parents)
            and (
                path.suffix.casefold() in {".html", ".htm"}
                or path.name.casefold() in {"subgroup.cif", "data.isoviz", "topas.str"}
            )
        ]
        for path in sorted(
            dependencies,
            key=lambda item: item.relative_to(case_directory).as_posix().casefold(),
        ):
            label = path.relative_to(case_directory).as_posix()
            _hash_file(digest, path, f"official:{case_key}:{label}")
        if not dependencies:
            digest.update(b"<no-validation-files>\0")
    return digest.hexdigest()


def _case_key(parent_cif: str, case_id: str) -> str:
    return f"{parent_cif}::{case_id}"


def _load_checkpoint(path: Path, signature: str, *, resume: bool) -> dict[str, Any]:
    empty = {
        "schema": CHECKPOINT_SCHEMA,
        "signature": signature,
        "created_at": _now(),
        "results": {},
    }
    if not resume or not path.is_file():
        return empty
    try:
        existing = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return empty
    if existing.get("schema") != CHECKPOINT_SCHEMA or existing.get("signature") != signature:
        return empty
    if not isinstance(existing.get("results"), dict):
        return empty
    return existing


def _manifest_query(case: dict[str, Any]) -> dict[str, Any]:
    return {
        "distortion_types": case.get("distortion_types"),
        "point_group": case.get("point_group"),
        "space_group_type": case.get("space_group_type"),
        "supercell_basis": case.get("supercell_basis"),
        "direct_sublattice_centering": case.get("direct_sublattice_centering"),
        "lattice_type": case.get("lattice_type", "direct"),
        "generate_if_missing": case.get("generate_if_missing", case.get("category") == "parameter-k"),
    }


def _summary(cases: Sequence[dict[str, Any]], *, selected_count: int, complete: bool) -> dict[str, Any]:
    comparable = [
        case for case in cases if case["status"] != "skipped_no_unique_official_result"
    ]
    authoritative = [
        case for case in comparable if case.get("official_completeness_proven") is True
    ]
    provisional = [
        case
        for case in comparable
        if case.get("official_source") == "downloaded_candidate_inventory"
    ]

    def count_status(records: Sequence[dict[str, Any]], status: str) -> int:
        return sum(case["status"] == status for case in records)

    return {
        "complete": complete,
        "selected_case_count": selected_count,
        "comparable_case_count": len(comparable),
        "official_eligible_case_count": len(authoritative),
        "authoritative_result_table_case_count": len(authoritative),
        "provisional_candidate_inventory_case_count": len(provisional),
        "skipped_no_unique_official_result_count": sum(
            case["status"] == "skipped_no_unique_official_result" for case in cases
        ),
        "exact_match_count": count_status(comparable, "exact_match"),
        "affine_equivalent_count": count_status(comparable, "affine_equivalent"),
        "difference_count": count_status(comparable, "different"),
        "local_error_count": count_status(comparable, "local_error"),
        "authoritative_exact_match_count": count_status(authoritative, "exact_match"),
        "authoritative_affine_equivalent_count": count_status(
            authoritative, "affine_equivalent"
        ),
        "authoritative_difference_count": count_status(authoritative, "different"),
        "authoritative_local_error_count": count_status(authoritative, "local_error"),
        "provisional_exact_match_count": count_status(provisional, "exact_match"),
        "provisional_affine_equivalent_count": count_status(
            provisional, "affine_equivalent"
        ),
        "provisional_difference_count": count_status(provisional, "different"),
        "provisional_local_error_count": count_status(provisional, "local_error"),
        "resumed_case_count": sum(bool(case.get("resumed_from_checkpoint")) for case in cases),
    }


def _report_payload(
    *,
    signature: str,
    manifest_path: Path,
    output_compare_root: Path,
    filters: dict[str, Any],
    cases: list[dict[str, Any]],
    selected_count: int,
    complete: bool,
) -> dict[str, Any]:
    return {
        "schema": 3,
        "generated_at": _now(),
        "run_signature": signature,
        "manifest": str(manifest_path.resolve()),
        "output_compare_root": str(output_compare_root.resolve()),
        "filters": filters,
        "summary": _summary(cases, selected_count=selected_count, complete=complete),
        "cases": cases,
    }


def run_comparison(
    *,
    manifest_path: Path = MANIFEST_PATH,
    output_compare_root: Path = OUTPUT_COMPARE_ROOT,
    input_root: Path = INPUT_ROOT,
    report_path: Path = REPORT_PATH,
    checkpoint_path: Path = CHECKPOINT_PATH,
    parents: Sequence[str] = (),
    case_ids: Sequence[str] = (),
    resume: bool = True,
    allow_candidate_inventory: bool = False,
    iso_factory: Callable[[], Any] = IsoDistort,
) -> dict[str, Any]:
    """Run eligible manifest cases, checkpointing one complete comparison at a time."""

    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    all_specs = _manifest_cases(manifest)
    selected = select_manifest_cases(manifest, parents, case_ids)
    official_template = manifest.get("official_download_root")
    if not isinstance(official_template, str) or not official_template:
        raise ValueError("manifest is missing official_download_root")

    # The signature deliberately covers the whole manifest batch, not only the
    # current filter.  Thus ``--case`` runs can safely accumulate into and
    # resume from the same checkpoint as a later full run.
    official_inputs: dict[str, dict[str, Any]] = {}
    for spec in all_specs:
        case = spec["case"]
        key = _case_key(spec["parent_cif"], str(case["id"]))
        directory = _case_root(
            official_template,
            output_compare_root,
            spec["parent_cif"],
            str(case["readable_case_folder"]),
        )
        official_inputs[key] = _official_input(
            directory,
            parent_cif=spec["parent_cif"],
            parent_space_group_number=int(spec["parent_space_group_number"]),
            manifest_case=case,
            allow_candidate_inventory=allow_candidate_inventory,
        )

    signature = _run_signature(
        manifest_path,
        input_root,
        all_specs,
        official_inputs,
        allow_candidate_inventory=allow_candidate_inventory,
    )
    checkpoint = _load_checkpoint(checkpoint_path, signature, resume=resume)
    checkpoint_results: dict[str, Any] = checkpoint["results"]
    mutable_cache_keys = {
        _case_key(spec["parent_cif"], str(spec["case"]["id"]))
        for spec in all_specs
        if _case_uses_mutable_isotropy_cache(spec["case"])
    }
    # We cannot content-sign the generated cache through the public API.  Drop
    # any parameter-k results loaded from a prior invocation before considering
    # checkpoint reuse; ordinary cases remain safely resumable.
    for mutable_key in mutable_cache_keys:
        checkpoint_results.pop(mutable_key, None)
    parent_apis: dict[str, Any] = {}
    parent_setup_errors: dict[str, str] = {}
    report_cases: list[dict[str, Any]] = []
    # Also replace an incompatible/stale checkpoint when every selected case is
    # skipped, so ``--restart`` always has observable clearing semantics.
    _atomic_json(checkpoint_path, checkpoint)

    for spec in selected:
        parent_name = spec["parent_cif"]
        case = spec["case"]
        case_id = str(case["id"])
        key = _case_key(parent_name, case_id)
        official = official_inputs[key]
        base = {
            "key": key,
            "parent_cif": parent_name,
            "case_id": case_id,
            "readable_case_folder": case["readable_case_folder"],
            "manifest_query": _manifest_query(case),
            "official_case_directory": official["case_directory"],
            "official_source": official.get("official_source"),
            "official_completeness_proven": official.get(
                "official_completeness_proven", False
            ),
            "checkpoint_reuse_policy": (
                "disabled_mutable_generated_isotropy_cache"
                if key in mutable_cache_keys
                else "content_signed"
            ),
        }
        if not official["eligible"]:
            report_cases.append(
                {
                    **base,
                    "status": "skipped_no_unique_official_result",
                    "reason": official["reason"],
                    "result_pages": official["result_pages"],
                    "html_pages": official["html_pages"],
                    **({"parse_error": official["parse_error"]} if official.get("parse_error") else {}),
                    **(
                        {"context_errors": official["context_errors"]}
                        if official.get("context_errors")
                        else {}
                    ),
                }
            )
            continue

        if key in checkpoint_results:
            resumed_record = deepcopy(checkpoint_results[key])
            resumed_record["resumed_from_checkpoint"] = True
            report_cases.append(resumed_record)
            continue

        started = time.perf_counter()
        try:
            if parent_name in parent_setup_errors:
                raise RuntimeError(parent_setup_errors[parent_name])
            if parent_name not in parent_apis:
                try:
                    api = iso_factory()
                    with contextlib.redirect_stdout(io.StringIO()):
                        api.load_structure(input_root / parent_name)
                    parent_apis[parent_name] = api
                except Exception as exc:
                    message = f"parent backend setup failed: {type(exc).__name__}: {exc}"
                    parent_setup_errors[parent_name] = message
                    raise RuntimeError(message) from exc
            api = parent_apis[parent_name]
            query = _manifest_query(case)
            with contextlib.redirect_stdout(io.StringIO()):
                items = api.search_method_3(**query)
            validator = getattr(api, "_sym_val", None)
            comparison = compare_embedding_sets(
                official["identities"],
                items,
                parent_structure=getattr(api, "structure", None),
                symmetry_tolerance=getattr(validator, "tolerance", None),
            )
            record = {
                **base,
                "status": comparison["status"],
                "official_result_page": official.get("result_page"),
                "official_download_audit_status": official["download_audit_status"],
                **(
                    {
                        "provisional_reason": official["provisional_reason"],
                        "official_audit_issues": official["audit_issues"],
                    }
                    if not official.get("official_completeness_proven", False)
                    else {}
                ),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "comparison": comparison,
            }
        except Exception as exc:  # noqa: BLE001 - persist scientific/backend failures and continue the batch
            record = {
                **base,
                "status": "local_error",
                "official_result_page": official.get("result_page"),
                "official_download_audit_status": official["download_audit_status"],
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "error": f"{type(exc).__name__}: {exc}",
            }
        checkpoint_results[key] = record
        checkpoint["updated_at"] = _now()
        _atomic_json(checkpoint_path, checkpoint)
        report_cases.append(record)
        partial = _report_payload(
            signature=signature,
            manifest_path=manifest_path,
            output_compare_root=output_compare_root,
            filters={
                "parents": list(parents),
                "case_ids": list(case_ids),
                "allow_candidate_inventory": allow_candidate_inventory,
            },
            cases=report_cases,
            selected_count=len(selected),
            complete=False,
        )
        _atomic_json(report_path, partial)

    report = _report_payload(
        signature=signature,
        manifest_path=manifest_path,
        output_compare_root=output_compare_root,
        filters={
            "parents": list(parents),
            "case_ids": list(case_ids),
            "allow_candidate_inventory": allow_candidate_inventory,
        },
        cases=report_cases,
        selected_count=len(selected),
        complete=True,
    )
    _atomic_json(report_path, report)
    return report


def _print_summary(report: dict[str, Any], report_path: Path, checkpoint_path: Path) -> None:
    summary = report["summary"]
    print(
        "Method 3 official/local comparison: "
        f"comparable={summary['comparable_case_count']} "
        f"({summary['official_eligible_case_count']} authoritative + "
        f"{summary['provisional_candidate_inventory_case_count']} provisional), "
        f"exact={summary['exact_match_count']}, "
        f"affine-equivalent={summary['affine_equivalent_count']}, "
        f"differences={summary['difference_count']}, "
        f"errors={summary['local_error_count']}, "
        f"skipped={summary['skipped_no_unique_official_result_count']}"
    )
    for case in report["cases"]:
        suffix = " (checkpoint)" if case.get("resumed_from_checkpoint") else ""
        if "comparison" in case:
            comparison = case["comparison"]
            detail = f"official={comparison['official_count']} local={comparison['local_embedding_count']}"
        else:
            detail = case.get("reason", case.get("error", ""))
        print(f"[{case['status']}] {case['case_id']}: {detail}{suffix}")
    print(f"Report: {report_path.resolve()}")
    print(f"Checkpoint: {checkpoint_path.resolve()}")


def _comparison_exit_code(
    report: Mapping[str, Any], *, accept_provisional: bool = False
) -> int:
    """Return 0 only for complete, accepted evidence with no scientific mismatch."""

    summary = report["summary"]
    if summary["difference_count"] or summary["local_error_count"]:
        return 1
    if summary["comparable_case_count"] == 0:
        return 2
    if summary["skipped_no_unique_official_result_count"]:
        return 2
    if summary["provisional_candidate_inventory_case_count"] and not accept_provisional:
        return 2
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run manifest Method 3 queries and compare both literal and exact-affine "
            "local embeddings with official HTML."
        )
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output-compare-root", type=Path, default=OUTPUT_COMPARE_ROOT)
    parser.add_argument("--input-root", type=Path, default=INPUT_ROOT)
    parser.add_argument("--report", type=Path, default=REPORT_PATH)
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT_PATH)
    parser.add_argument("--parent", action="append", default=[], help="exact parent CIF name; repeatable")
    parser.add_argument("--case", dest="case_ids", action="append", default=[], help="exact case id; repeatable")
    parser.add_argument("--restart", action="store_true", help="ignore an otherwise compatible checkpoint")
    parser.add_argument(
        "--allow-candidate-inventory",
        action="store_true",
        help=(
            "provisionally compare identities recovered from downloaded subgroup.cif files "
            "when the authoritative result-table HTML is missing; this does not prove completeness"
        ),
    )
    parser.add_argument(
        "--accept-provisional",
        action="store_true",
        help=(
            "allow an otherwise successful provisional candidate-inventory comparison "
            "to exit 0; differences, errors and skipped cases still exit nonzero"
        ),
    )
    args = parser.parse_args(argv)

    try:
        report = run_comparison(
            manifest_path=args.manifest,
            output_compare_root=args.output_compare_root,
            input_root=args.input_root,
            report_path=args.report,
            checkpoint_path=args.checkpoint,
            parents=args.parent,
            case_ids=args.case_ids,
            resume=not args.restart,
            allow_candidate_inventory=args.allow_candidate_inventory,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    _print_summary(report, args.report, args.checkpoint)
    return _comparison_exit_code(report, accept_provisional=args.accept_provisional)


if __name__ == "__main__":
    raise SystemExit(main())
