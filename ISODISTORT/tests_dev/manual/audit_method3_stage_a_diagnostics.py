"""Audit the opt-in Method 3 affine-only Stage-A enumerator.

This diagnostic deliberately isolates the exact affine first stage from the
single-IR backend.  It runs the production ``IsoSearchEngine.method_3_search``
path with ``include_affine_only_diagnostics=True`` and a null route provider,
then compares the resulting affine subgroup embeddings with the identities in
the saved official route-audit report.

By default, only cases containing at least one official
``coupled_ir_required`` identity are selected; ``--all-cases`` extends the
same audit to the whole route report.  Exact Fraction Seitz reconstruction is
used for every edge in the official/candidate bipartite graph.  Parent-affine
conjugacy is accepted in the same way as the main Method-3 comparison tool,
while multiplicity is still preserved by maximum one-to-one matching.

Passing this audit proves affine membership only.  It does *not* prove that an
embedding has a nonzero fixed subspace for the selected displacive/strain
representations, and it does not supply coupled-IR amplitudes or modes.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import sys
import time
import warnings
from copy import deepcopy
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from isocore.distortion.search_methods import (  # noqa: E402
    IsoSearchEngine,
    Method3Query,
)
from isocore.structure import read_cif  # noqa: E402
from isocore.utils import get_config  # noqa: E402

try:  # Package import in pytest; direct import when this file is executed.
    from .audit_method3_downloads import SubgroupIdentity
    from .compare_method3_official_local import _identity_from_local
    from .method3_affine_equivalence import (
        affine_equivalence_witness,
        embedding_operations,
        parent_operations,
    )
except ImportError:  # pragma: no cover - direct CLI invocation
    from audit_method3_downloads import SubgroupIdentity
    from compare_method3_official_local import _identity_from_local
    from method3_affine_equivalence import (
        affine_equivalence_witness,
        embedding_operations,
        parent_operations,
    )


MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method3_download_manifest.json"
ROUTE_AUDIT_PATH = PACKAGE_ROOT / "output" / "validation" / "method3_embedding_route_audit.json"
REPORT_PATH = PACKAGE_ROOT / "output" / "validation" / "method3_stage_a_diagnostic_audit.json"
CHECKPOINT_PATH = (
    PACKAGE_ROOT
    / "output"
    / "validation"
    / "method3_stage_a_diagnostic_audit_checkpoint.json"
)
INPUT_ROOT = WORKSPACE / "experiment_data"
SCHEMA = 2
CHECKPOINT_SCHEMA = 1
RUNTIME_DISTRIBUTIONS = ("numpy", "spglib", "pymatgen")


class _NullRouteBackend:
    """Make the production search path return Stage-A candidates only."""

    def enumerate_all_special_subgroups(self, *_args: object) -> list[object]:
        return []

    def list_k_points(self, *_args: object) -> list[object]:
        return []


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _runtime_versions() -> dict[str, Any]:
    """Describe the interpreter and scientific stack that affect the audit."""

    distributions: dict[str, str] = {}
    for name in RUNTIME_DISTRIBUTIONS:
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "<missing>"
        except Exception as exc:  # noqa: BLE001 - provenance must remain serializable
            error_type = f"{type(exc).__module__}.{type(exc).__qualname__}"
            distributions[name] = f"<unavailable:{error_type}>"
    return {
        "python": {
            "implementation": sys.implementation.name,
            "version": list(sys.version_info[:5]),
            "cache_tag": sys.implementation.cache_tag,
            "build": sys.version,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "distributions": distributions,
    }


def _scientific_tolerances() -> dict[str, Any]:
    """Return configured tolerances together with their physical semantics."""

    cfg = get_config()
    return {
        "configured": {
            "symmetry_cartesian_tolerance_angstrom": float(
                cfg.symmetry_cartesian_tolerance_angstrom
            ),
            "symmetry_angle_tolerance_degrees": float(
                cfg.symmetry_angle_tolerance_degrees
            ),
            "affine_exact_cartesian_tolerance_angstrom": float(
                cfg.affine_exact_cartesian_tolerance_angstrom
            ),
            "fractional_coordinate_tolerance": float(
                cfg.fractional_coordinate_tolerance
            ),
        },
        "semantics": {
            "symmetry_cartesian_tolerance_angstrom": (
                "Absolute Cartesian distance passed to pymatgen/spglib when "
                "detecting the uploaded parent symmetry."
            ),
            "symmetry_angle_tolerance_degrees": (
                "Lattice-angle tolerance passed to pymatgen/spglib for parent "
                "symmetry detection."
            ),
            "affine_exact_cartesian_tolerance_angstrom": (
                "Configured Cartesian floor at synthetic spglib boundaries in "
                "Stage A; the implementation also applies its documented "
                "machine-roundoff floor before exact Seitz reconstruction."
            ),
            "fractional_coordinate_tolerance": (
                "Dimensionless float-to-rational residual accepted only at the "
                "standard-setting bridge; final affine equality is exact."
            ),
        },
    }


def _signature_inputs(
    manifest: dict[str, Any],
    manifest_path: Path,
    route_audit_path: Path,
    input_root: Path,
    *,
    package_root: Path,
    tool_path: Path,
) -> list[tuple[str, Path]]:
    """List every immutable scientific input needed for safe checkpoint reuse."""

    inputs = [
        ("manifest", manifest_path),
        ("official-route-audit", route_audit_path),
        ("source:stage-a-audit", tool_path),
        (
            "source:download-audit-parser",
            package_root / "tests_dev" / "manual" / "audit_method3_downloads.py",
        ),
        (
            "source:official-local-identity-adapter",
            package_root
            / "tests_dev"
            / "manual"
            / "compare_method3_official_local.py",
        ),
        (
            "source:affine-comparator",
            package_root
            / "tests_dev"
            / "manual"
            / "method3_affine_equivalence.py",
        ),
        (
            "source:search-methods",
            package_root / "isocore" / "distortion" / "search_methods.py",
        ),
        (
            "source:affine-embeddings",
            package_root / "isocore" / "distortion" / "affine_embeddings.py",
        ),
        ("source:lattice", package_root / "isocore" / "utils" / "lattice.py"),
        (
            "source:config-loader",
            package_root / "isocore" / "utils" / "config_loader.py",
        ),
        (
            "source:cif-reader",
            package_root / "isocore" / "structure" / "cif_io.py",
        ),
        (
            "source:structure-api",
            package_root / "isocore" / "structure" / "__init__.py",
        ),
        (
            "source:hall-setting-map",
            package_root / "isocore" / "io" / "isodistort_cif.py",
        ),
        ("config", package_root / "config" / "settings.yaml"),
    ]
    parent_names = sorted(
        {
            str(parent.get("parent_cif", "")).strip()
            for parent in manifest.get("parents", [])
            if str(parent.get("parent_cif", "")).strip()
        },
        key=str.casefold,
    )
    inputs.extend(
        (f"parent-cif:{parent_name}", input_root / parent_name)
        for parent_name in parent_names
    )
    return inputs


def _signature_state(
    manifest_path: Path,
    route_audit_path: Path,
    input_root: Path,
    *,
    manifest: dict[str, Any] | None = None,
    package_root: Path | None = None,
    tool_path: Path | None = None,
    runtime_versions: dict[str, Any] | None = None,
    tolerances: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the aggregate content/runtime signature for an interrupted run."""

    package_root = (PACKAGE_ROOT if package_root is None else package_root).resolve()
    tool_path = (Path(__file__) if tool_path is None else tool_path).resolve()
    if manifest is None:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if runtime_versions is None:
        runtime_versions = _runtime_versions()
    if tolerances is None:
        tolerances = _scientific_tolerances()
    inputs = _signature_inputs(
        manifest,
        manifest_path.resolve(),
        route_audit_path.resolve(),
        input_root.resolve(),
        package_root=package_root,
        tool_path=tool_path,
    )

    digest = hashlib.sha256(usedforsecurity=False)
    digest.update(
        json.dumps(
            {
                "checkpoint_schema": CHECKPOINT_SCHEMA,
                "report_schema": SCHEMA,
                "runtime_versions": runtime_versions,
                "scientific_tolerances": tolerances,
            },
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    digest.update(b"\0")

    records: list[dict[str, Any]] = []
    for label, path in inputs:
        resolved = path.resolve()
        exists = resolved.is_file()
        content_hash = _sha256(resolved) if exists else None
        digest.update(label.encode("utf-8"))
        digest.update(b"\0")
        digest.update((content_hash or "<missing>").encode("ascii"))
        digest.update(b"\0")
        records.append(
            {
                "label": label,
                "path": str(resolved),
                "exists": exists,
                "sha256": content_hash,
            }
        )
    return {
        "signature": digest.hexdigest(),
        "algorithm": "sha256(label, file-content-sha256, runtime, tolerances)",
        "runtime_versions": runtime_versions,
        "scientific_tolerances": tolerances,
        "inputs": records,
    }


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _case_key(parent_cif: str, case_id: str) -> str:
    return f"{parent_cif}::{case_id}"


def _new_checkpoint(signature_state: dict[str, Any]) -> dict[str, Any]:
    now = _now()
    return {
        "schema": CHECKPOINT_SCHEMA,
        "run_signature": signature_state["signature"],
        "signature_state": signature_state,
        "created_at": now,
        "updated_at": now,
        "results": {},
    }


def _load_checkpoint(
    path: Path,
    signature_state: dict[str, Any],
    *,
    resume: bool,
) -> dict[str, Any]:
    """Load a compatible checkpoint or reject stale scientific state.

    A signature mismatch is not silently treated as a fresh run: doing so can
    overwrite the only evidence from an interrupted calculation.  Callers
    must explicitly request a restart before replacing an incompatible file.
    """

    fresh = _new_checkpoint(signature_state)
    if not resume or not path.is_file():
        return fresh
    try:
        existing = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"cannot resume unreadable Stage-A checkpoint {path}; use --restart "
            "only after preserving it for diagnosis"
        ) from exc
    if existing.get("schema") != CHECKPOINT_SCHEMA:
        raise RuntimeError(
            "Stage-A checkpoint schema mismatch; use --restart to begin a new "
            "content-signed run"
        )
    if existing.get("run_signature") != signature_state["signature"]:
        raise RuntimeError(
            "Stage-A checkpoint signature mismatch: source, configuration, "
            "official evidence, parent CIF, tolerance, or runtime changed; "
            "use --restart to recompute instead of mixing runs"
        )
    if not isinstance(existing.get("results"), dict):
        raise RuntimeError(
            "Stage-A checkpoint has an invalid results object; use --restart "
            "only after preserving it for diagnosis"
        )
    return existing


def _identity_from_json(item: dict[str, Any]) -> SubgroupIdentity:
    return SubgroupIdentity(
        space_group_number=int(item["space_group_number"]),
        symbol=str(item["symbol"]),
        basis=tuple(
            tuple(Fraction(value) for value in row) for row in item["basis"]
        ),
        origin=tuple(Fraction(value) for value in item["origin"]),
        sublattice_index=int(item["s"]),
        group_index=int(item["i"]),
    )


def _manifest_specs(manifest: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    specs: dict[tuple[str, str], dict[str, Any]] = {}
    for parent in manifest.get("parents", []):
        parent_cif = str(parent["parent_cif"])
        parent_sg = int(parent["parent_space_group_type"])
        for case in parent.get("cases", []):
            specs[(parent_cif, str(case["id"]))] = {
                "parent_cif": parent_cif,
                "parent_space_group_number": parent_sg,
                "case": case,
            }
    return specs


def _maximum_matching(
    edges: dict[int, list[dict[str, Any]]],
    official_count: int,
    candidate_count: int,
) -> dict[str, Any]:
    candidate_matches: dict[int, tuple[int, dict[str, Any]]] = {}

    def assign(official_index: int, visited: set[int]) -> bool:
        ordered = sorted(
            edges.get(official_index, []),
            key=lambda edge: (
                edge["relationship"] != "same_affine_subgroup",
                edge["candidate_index"],
            ),
        )
        for edge in ordered:
            candidate_index = int(edge["candidate_index"])
            if candidate_index in visited:
                continue
            visited.add(candidate_index)
            previous = candidate_matches.get(candidate_index)
            if previous is None or assign(previous[0], visited):
                candidate_matches[candidate_index] = (official_index, edge)
                return True
        return False

    for official_index in range(official_count):
        assign(official_index, set())

    matches = [
        {
            "official_index": official_index,
            "candidate_index": candidate_index,
            **{key: value for key, value in edge.items() if key != "candidate_index"},
        }
        for candidate_index, (official_index, edge) in sorted(
            candidate_matches.items(), key=lambda item: item[1][0]
        )
    ]
    matched_official = {item["official_index"] for item in matches}
    matched_candidates = {item["candidate_index"] for item in matches}
    return {
        "match_count": len(matches),
        "matches": matches,
        "unmatched_official_indices": [
            index for index in range(official_count) if index not in matched_official
        ],
        "unmatched_candidate_indices": [
            index for index in range(candidate_count) if index not in matched_candidates
        ],
    }


def _equivalence_classes(
    models: list[Any], parent_model: Any
) -> tuple[list[list[int]], list[dict[str, Any]]]:
    """Partition embeddings by exact parent-affine conjugacy."""

    parents = list(range(len(models)))

    def root(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = root(left)
        right_root = root(right)
        if left_root != right_root:
            parents[right_root] = left_root

    witnesses: list[dict[str, Any]] = []
    for left in range(len(models)):
        for right in range(left + 1, len(models)):
            witness = affine_equivalence_witness(
                models[left], models[right], parent_model
            )
            if witness is None:
                continue
            union(left, right)
            witnesses.append({"left": left, "right": right, **witness})
    classes: dict[int, list[int]] = {}
    for index in range(len(models)):
        classes.setdefault(root(index), []).append(index)
    return list(classes.values()), witnesses


def _summarize(cases: list[dict[str, Any]], expected_case_count: int) -> dict[str, Any]:
    completed = [case for case in cases if case.get("status") == "completed"]
    authoritative = [case for case in completed if case["official_completeness_proven"]]
    provisional = [case for case in completed if not case["official_completeness_proven"]]

    def total(records: list[dict[str, Any]], field: str) -> int:
        return sum(int(record.get(field, 0)) for record in records)

    return {
        "complete": len(cases) == expected_case_count
        and all(case.get("status") == "completed" for case in cases),
        "selected_case_count": expected_case_count,
        "completed_case_count": len(completed),
        "error_case_count": sum(case.get("status") == "error" for case in cases),
        "resumed_case_count": sum(
            bool(case.get("resumed_from_checkpoint")) for case in cases
        ),
        "authoritative_case_count": len(authoritative),
        "provisional_case_count": len(provisional),
        "official_embedding_count": total(completed, "official_embedding_count"),
        "coupled_official_embedding_count": total(
            completed, "coupled_official_embedding_count"
        ),
        "stage_a_candidate_count": total(completed, "stage_a_candidate_count"),
        "matched_official_embedding_count": total(
            completed, "matched_official_embedding_count"
        ),
        "matched_coupled_official_embedding_count": total(
            completed, "matched_coupled_official_embedding_count"
        ),
        "unmatched_official_embedding_count": total(
            completed, "unmatched_official_embedding_count"
        ),
        "unmatched_coupled_official_embedding_count": total(
            completed, "unmatched_coupled_official_embedding_count"
        ),
        "unmatched_stage_a_candidate_count": total(
            completed, "unmatched_stage_a_candidate_count"
        ),
        "stage_a_candidate_without_any_official_equivalent_count": total(
            completed, "stage_a_candidate_without_any_official_equivalent_count"
        ),
        "stage_a_surplus_multiplicity_count": total(
            completed, "stage_a_surplus_multiplicity_count"
        ),
        "stage_a_surplus_parent_conjugate_only_count": total(
            completed, "stage_a_surplus_parent_conjugate_only_count"
        ),
        "stage_a_surplus_with_exact_edge_count": total(
            completed, "stage_a_surplus_with_exact_edge_count"
        ),
        "stage_a_parent_affine_orbit_count": total(
            completed, "stage_a_parent_affine_orbit_count"
        ),
        "stage_a_orbit_without_exact_official_representative_count": total(
            completed,
            "stage_a_orbit_without_exact_official_representative_count",
        ),
        "stage_a_orbit_with_multiple_exact_official_representatives_count": total(
            completed,
            "stage_a_orbit_with_multiple_exact_official_representatives_count",
        ),
        "official_parent_conjugate_pair_count": total(
            completed, "official_parent_conjugate_pair_count"
        ),
    }


def run_audit(
    *,
    manifest_path: Path = MANIFEST_PATH,
    route_audit_path: Path = ROUTE_AUDIT_PATH,
    report_path: Path = REPORT_PATH,
    checkpoint_path: Path | None = None,
    input_root: Path = INPUT_ROOT,
    case_ids: set[str] | None = None,
    coupled_only: bool = True,
    resume: bool = True,
) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    route_audit = json.loads(route_audit_path.read_text(encoding="utf-8"))
    specs = _manifest_specs(manifest)
    selected_route_cases = []
    for case in route_audit.get("cases", []):
        has_coupled = any(
            embedding.get("classification") == "coupled_ir_required"
            for embedding in case.get("embeddings", [])
        )
        if coupled_only and not has_coupled:
            continue
        if case_ids is not None and str(case["case_id"]) not in case_ids:
            continue
        selected_route_cases.append(case)

    checkpoint_path = checkpoint_path or report_path.with_name(
        f"{report_path.stem}_checkpoint.json"
    )
    tolerances = _scientific_tolerances()
    signature_state = _signature_state(
        manifest_path,
        route_audit_path,
        input_root,
        manifest=manifest,
        runtime_versions=_runtime_versions(),
        tolerances=tolerances,
    )
    checkpoint = _load_checkpoint(
        checkpoint_path,
        signature_state,
        resume=resume,
    )
    checkpoint_results: dict[str, Any] = checkpoint["results"]
    checkpoint["updated_at"] = _now()
    _atomic_json(checkpoint_path, checkpoint)

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "run_signature": signature_state["signature"],
        "generated_at": _now(),
        "tool": str(Path(__file__).resolve()),
        "scientific_scope": {
            "proves": (
                "Exact affine subgroup membership/equivalence modulo the child primitive "
                "translation lattice and under parent-space-group affine conjugacy."
            ),
            "does_not_prove": (
                "A nonzero fixed subspace for the selected strain/displacive "
                "representations, coupled-IR order-parameter activity, mode amplitudes, "
                "energetic accessibility, or a continuous phase-transition path."
            ),
            "stage_a_isolation": (
                "The production method_3_search path is called with "
                "include_affine_only_diagnostics=True while a null route backend "
                "suppresses all single-IR candidates."
            ),
        },
        "inputs": {
            "manifest": str(manifest_path.resolve()),
            "route_audit": str(route_audit_path.resolve()),
            "input_root": str(input_root.resolve()),
        },
        "provenance": signature_state,
        "checkpoint": {
            "path": str(checkpoint_path.resolve()),
            "resume_requested": resume,
            "compatible_saved_result_count": len(checkpoint_results),
        },
        "filters": {
            "case_ids": sorted(case_ids) if case_ids else [],
            "coupled_only": coupled_only,
        },
        "summary": {},
        "cases": [],
    }

    structure_cache: dict[str, Any] = {}
    parent_model_cache: dict[str, Any] = {}
    rotations_cache: dict[str, list[list[list[int]]]] = {}
    configured = tolerances["configured"]
    # These names and units are deliberately repeated at the call boundary:
    # a Cartesian Å tolerance must never become a fractional-coordinate or
    # dimensionless lattice residual by accident.
    symmetry_tolerance_angstrom = float(
        configured["symmetry_cartesian_tolerance_angstrom"]
    )
    symmetry_angle_tolerance_degrees = float(
        configured["symmetry_angle_tolerance_degrees"]
    )
    affine_exact_tolerance_angstrom = float(
        configured["affine_exact_cartesian_tolerance_angstrom"]
    )
    fractional_coordinate_tolerance = float(
        configured["fractional_coordinate_tolerance"]
    )

    report["summary"] = _summarize([], len(selected_route_cases))
    _atomic_json(report_path, report)

    for route_case in selected_route_cases:
        parent_cif = str(route_case["parent_cif"])
        case_id = str(route_case["case_id"])
        spec = specs[(parent_cif, case_id)]
        case = spec["case"]
        key = _case_key(parent_cif, case_id)
        if key in checkpoint_results:
            record = deepcopy(checkpoint_results[key])
            record["resumed_from_checkpoint"] = True
            report["cases"].append(record)
            report["generated_at"] = _now()
            report["summary"] = _summarize(
                report["cases"], len(selected_route_cases)
            )
            _atomic_json(report_path, report)
            continue

        started = time.perf_counter()
        base = {
            "parent_cif": parent_cif,
            "case_id": case_id,
            "official_source": route_case.get("official_source"),
            "official_completeness_proven": bool(
                route_case.get("official_completeness_proven")
            ),
        }
        try:
            if parent_cif not in structure_cache:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                structure = read_cif(input_root / parent_cif)
                dataset = SpacegroupAnalyzer(
                    structure,
                    symprec=symmetry_tolerance_angstrom,
                    angle_tolerance=symmetry_angle_tolerance_degrees,
                ).get_symmetry_dataset()
                if dataset is None:
                    raise ValueError(f"spglib failed for parent {parent_cif}")
                structure_cache[parent_cif] = structure
                rotations_cache[parent_cif] = [
                    rotation.tolist() for rotation in dataset.rotations
                ]
                parent_model_cache[parent_cif] = parent_operations(
                    structure,
                    symmetry_tolerance_angstrom,
                    symmetry_angle_tolerance_degrees,
                )
            structure = structure_cache[parent_cif]
            parent_model = parent_model_cache[parent_cif]

            embeddings = route_case.get("embeddings", [])
            official = [
                _identity_from_json(item["official_identity"])
                for item in embeddings
            ]
            coupled_indices = [
                index
                for index, item in enumerate(embeddings)
                if item.get("classification") == "coupled_ir_required"
            ]

            query = Method3Query(
                distortion_types=case.get("distortion_types"),
                point_group=case.get("point_group"),
                space_group_type=case.get("space_group_type"),
                supercell_basis=case.get("supercell_basis"),
                direct_sublattice_centering=case.get(
                    "direct_sublattice_centering"
                ),
                lattice_type=case.get("lattice_type", "direct"),
                parent_rotations=rotations_cache[parent_cif],
                parent_structure=structure,
                symmetry_tolerance=symmetry_tolerance_angstrom,
                symmetry_angle_tolerance_degrees=(
                    symmetry_angle_tolerance_degrees
                ),
                affine_exact_tolerance=affine_exact_tolerance_angstrom,
                fractional_coordinate_tolerance=(
                    fractional_coordinate_tolerance
                ),
                generate_if_missing=False,
                include_affine_only_diagnostics=True,
            )
            stage_a_items = IsoSearchEngine(_NullRouteBackend()).method_3_search(
                int(spec["parent_space_group_number"]), query
            )
            stage_a = [_identity_from_local(item) for item in stage_a_items]
            route_contract_errors = [
                {
                    "candidate_index": index,
                    "route_resolution": item.route_resolution,
                    "route_count": len(item.routes or []),
                }
                for index, item in enumerate(stage_a_items)
                if item.route_resolution != "affine_only_unresolved_coupled_route"
                or item.routes
            ]

            official_models = [
                embedding_operations(identity, parent_model) for identity in official
            ]
            candidate_models = [
                embedding_operations(identity, parent_model) for identity in stage_a
            ]
            candidate_classes, candidate_class_witnesses = _equivalence_classes(
                candidate_models, parent_model
            )
            _official_classes, official_class_witnesses = _equivalence_classes(
                official_models, parent_model
            )
            edges: dict[int, list[dict[str, Any]]] = {}
            candidate_edge_counts = [0 for _ in stage_a]
            for official_index, official_model in enumerate(official_models):
                for candidate_index, candidate_model in enumerate(candidate_models):
                    witness = affine_equivalence_witness(
                        candidate_model, official_model, parent_model
                    )
                    if witness is None:
                        continue
                    edges.setdefault(official_index, []).append(
                        {"candidate_index": candidate_index, **witness}
                    )
                    candidate_edge_counts[candidate_index] += 1

            matching = _maximum_matching(edges, len(official), len(stage_a))
            unmatched_candidate_indices = matching["unmatched_candidate_indices"]
            novel_indices = [
                index for index in unmatched_candidate_indices if not candidate_edge_counts[index]
            ]
            surplus_indices = [
                index for index in unmatched_candidate_indices if candidate_edge_counts[index]
            ]
            exact_edge_candidate_indices = {
                int(edge["candidate_index"])
                for official_edges in edges.values()
                for edge in official_edges
                if edge["relationship"] == "same_affine_subgroup"
            }
            surplus_with_exact_edge = [
                index for index in surplus_indices if index in exact_edge_candidate_indices
            ]
            surplus_parent_conjugate_only = [
                index for index in surplus_indices if index not in exact_edge_candidate_indices
            ]
            unmatched_official = set(matching["unmatched_official_indices"])
            unmatched_coupled = [
                index for index in coupled_indices if index in unmatched_official
            ]
            coupled_membership = [
                {
                    "official_index": index,
                    "official_identity": official[index].to_json(),
                    "matched": bool(edges.get(index)),
                    "candidate_matches": edges.get(index, []),
                }
                for index in coupled_indices
            ]
            orbit_records: list[dict[str, Any]] = []
            for orbit_index, members in enumerate(candidate_classes):
                exact_official_indices = sorted(
                    {
                        official_index
                        for official_index, official_edges in edges.items()
                        for edge in official_edges
                        if edge["candidate_index"] in members
                        and edge["relationship"] == "same_affine_subgroup"
                    }
                )
                affine_official_indices = sorted(
                    {
                        official_index
                        for official_index, official_edges in edges.items()
                        for edge in official_edges
                        if edge["candidate_index"] in members
                    }
                )
                orbit_records.append(
                    {
                        "orbit_index": orbit_index,
                        "candidate_indices": members,
                        "exact_official_indices": exact_official_indices,
                        "affine_equivalent_official_indices": affine_official_indices,
                    }
                )
            orbits_without_exact_official = [
                orbit["orbit_index"]
                for orbit in orbit_records
                if not orbit["exact_official_indices"]
            ]
            orbits_with_multiple_exact_official = [
                orbit["orbit_index"]
                for orbit in orbit_records
                if len(orbit["exact_official_indices"]) > 1
            ]

            record = {
                **base,
                "status": "completed",
                "resumed_from_checkpoint": False,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "query": {
                    "space_group_type": case.get("space_group_type"),
                    "supercell_basis": case.get("supercell_basis"),
                    "direct_sublattice_centering": case.get(
                        "direct_sublattice_centering"
                    ),
                    "distortion_types": case.get("distortion_types"),
                    "include_affine_only_diagnostics": True,
                    "scientific_tolerances": deepcopy(configured),
                },
                "official_embedding_count": len(official),
                "coupled_official_embedding_count": len(coupled_indices),
                "stage_a_candidate_count": len(stage_a),
                "matched_official_embedding_count": matching["match_count"],
                "matched_coupled_official_embedding_count": (
                    len(coupled_indices) - len(unmatched_coupled)
                ),
                "unmatched_official_embedding_count": len(
                    matching["unmatched_official_indices"]
                ),
                "unmatched_coupled_official_embedding_count": len(unmatched_coupled),
                "unmatched_stage_a_candidate_count": len(unmatched_candidate_indices),
                "stage_a_candidate_without_any_official_equivalent_count": len(
                    novel_indices
                ),
                "stage_a_surplus_multiplicity_count": len(surplus_indices),
                "stage_a_surplus_parent_conjugate_only_count": len(
                    surplus_parent_conjugate_only
                ),
                "stage_a_surplus_with_exact_edge_count": len(
                    surplus_with_exact_edge
                ),
                "stage_a_parent_affine_orbit_count": len(candidate_classes),
                "stage_a_orbit_without_exact_official_representative_count": len(
                    orbits_without_exact_official
                ),
                "stage_a_orbit_with_multiple_exact_official_representatives_count": len(
                    orbits_with_multiple_exact_official
                ),
                "official_parent_conjugate_pair_count": len(
                    official_class_witnesses
                ),
                "official_embeddings": [identity.to_json() for identity in official],
                "stage_a_candidates": [
                    {
                        "identity": identity.to_json(),
                        "route_resolution": item.route_resolution,
                        "route_count": len(item.routes or []),
                        "has_any_official_equivalent": bool(candidate_edge_counts[index]),
                    }
                    for index, (identity, item) in enumerate(
                        zip(stage_a, stage_a_items, strict=True)
                    )
                ],
                "coupled_membership": coupled_membership,
                "maximum_one_to_one_matching": matching,
                "all_affine_equivalence_edges": {
                    str(index): official_edges
                    for index, official_edges in sorted(edges.items())
                },
                "stage_a_parent_affine_orbits": orbit_records,
                "stage_a_candidate_pair_witnesses": candidate_class_witnesses,
                "official_pair_witnesses": official_class_witnesses,
                "candidate_without_any_official_equivalent_indices": novel_indices,
                "surplus_multiplicity_candidate_indices": surplus_indices,
                "surplus_parent_conjugate_only_candidate_indices": (
                    surplus_parent_conjugate_only
                ),
                "surplus_with_exact_edge_candidate_indices": surplus_with_exact_edge,
                "route_contract_errors": route_contract_errors,
            }
        except Exception as exc:  # noqa: BLE001 - preserve every diagnostic failure
            record = {
                **base,
                "status": "error",
                "resumed_from_checkpoint": False,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "error": f"{type(exc).__name__}: {exc}",
            }
        checkpoint_results[key] = deepcopy(record)
        checkpoint["updated_at"] = _now()
        _atomic_json(checkpoint_path, checkpoint)
        report["cases"].append(record)
        report["generated_at"] = _now()
        report["summary"] = _summarize(
            report["cases"], len(selected_route_cases)
        )
        _atomic_json(report_path, report)

    checkpoint["updated_at"] = _now()
    checkpoint["last_selection"] = deepcopy(report["filters"])
    checkpoint["last_summary"] = deepcopy(report["summary"])
    _atomic_json(checkpoint_path, checkpoint)
    report["checkpoint"]["saved_result_count"] = len(checkpoint_results)
    report["checkpoint"]["updated_at"] = checkpoint["updated_at"]
    _atomic_json(report_path, report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--report", type=Path, default=REPORT_PATH)
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT_PATH)
    parser.add_argument(
        "--restart",
        action="store_true",
        help="replace any checkpoint instead of resuming compatible signed cases",
    )
    parser.add_argument(
        "--all-cases",
        action="store_true",
        help="audit every case in the route report, not only cases with coupled-only rows",
    )
    args = parser.parse_args(argv)
    report = run_audit(
        report_path=args.report,
        checkpoint_path=args.checkpoint,
        case_ids=set(args.case) if args.case else None,
        coupled_only=not args.all_cases,
        resume=not args.restart,
    )
    print(json.dumps(report["summary"], ensure_ascii=False, indent=2))
    return 0 if report["summary"].get("complete") else 1


if __name__ == "__main__":
    raise SystemExit(main())
