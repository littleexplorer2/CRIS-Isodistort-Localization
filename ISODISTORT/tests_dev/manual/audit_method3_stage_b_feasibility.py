"""Batch-read Stage-A embeddings through the exact Stage-B feasibility core.

This is a diagnostic audit, not a product search path.  For each Method-3
case it builds ``N_G(T_s)/T_s`` and the selected physical representation once,
then locates and analyses every saved Stage-A candidate in that shared finite
group.  Results prove only embedding-level symmetry feasibility:

``P_H = |H|^-1 sum(D(h), h in H)``

and an embedding is feasible iff ``Fix(H)`` is non-zero and its pointwise
stabilizer is exactly ``H``.  IR decomposition, coupled-route labels, mode
amplitudes and selectable rows remain unresolved.

For scale, this batch uses the exactly equivalent character corollary
``dim Fix(H) = |H|^-1 sum_h trace(D(h))``.  An element ``g`` fixes all of
``Fix(H)`` iff ``dim Fix(<H,g>) == dim Fix(H)``.  No numerical rank or dense
projector multiplication is used.
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
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from isocore.distortion.affine_embeddings import (  # noqa: E402
    build_affine_quotient,
    embedding_from_identity,
    locate_embedding_subgroup,
    parent_affine_group,
)
from isocore.distortion.inverse_landau import (  # noqa: E402
    FixedSpaceFeasibilityAnalyzer,
)
from isocore.distortion.inverse_landau_adapter import (  # noqa: E402
    build_selected_character_representation,
    stable_embedding_id,
)
from isocore.structure import read_cif  # noqa: E402
from isocore.utils import get_config  # noqa: E402
from isocore.utils.lattice import same_lattice  # noqa: E402

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method3_download_manifest.json"
ROUTE_REPORT_PATH = (
    PACKAGE_ROOT / "output" / "validation" / "method3_embedding_route_audit.json"
)
STAGE_A_REPORT_PATH = (
    PACKAGE_ROOT / "output" / "validation" / "method3_stage_a_all_audit.json"
)
REPORT_PATH = (
    PACKAGE_ROOT / "output" / "validation" / "method3_stage_b_feasibility_audit.json"
)
CHECKPOINT_PATH = (
    PACKAGE_ROOT
    / "output"
    / "validation"
    / "method3_stage_b_feasibility_checkpoint.json"
)
INPUT_ROOT = WORKSPACE / "experiment_data"

SCHEMA = 1
CHECKPOINT_SCHEMA = 1
EXPECTED_FULL_COUNTS = {
    "case_count": 40,
    "official_embedding_count": 77,
    "stage_a_candidate_count": 142,
    "parent_affine_orbit_count": 80,
    "novel_orbit_count": 3,
}
EXPECTED_NOVEL_CASES = {"M3-ND-12": 1, "M3-ND-13": 2}
RUNTIME_DISTRIBUTIONS = ("numpy", "spglib", "pymatgen")


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256(usedforsecurity=False)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _runtime_versions() -> dict[str, Any]:
    distributions: dict[str, str] = {}
    for name in RUNTIME_DISTRIBUTIONS:
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "<missing>"
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


def _scientific_tolerances() -> dict[str, Any]:
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
            "site_matching": (
                "Configured Cartesian Å tolerance is used only to recover a "
                "unique finite-supercell site permutation from floating CIF coordinates."
            ),
            "finite_group": (
                "All quotient multiplication, representation homomorphism, "
                "character averages, generated-subgroup closures and stabilizers "
                "are then exact Fraction/integer arithmetic."
            ),
        },
    }


def _load_json(path: Path, description: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"cannot read {description}: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{description} must contain a JSON object")
    return payload


def _signature_state(
    *,
    manifest_path: Path,
    route_report_path: Path,
    stage_a_report_path: Path,
    input_root: Path,
    package_root: Path,
    tool_path: Path,
    runtime_versions: dict[str, Any],
    tolerances: dict[str, Any],
    max_quotient_order: int,
    case_ids: tuple[str, ...] | None,
) -> dict[str, Any]:
    manifest = _load_json(manifest_path, "Method-3 manifest")
    parent_names = sorted(
        {
            str(parent.get("parent_cif", "")).strip()
            for parent in manifest.get("parents", [])
            if str(parent.get("parent_cif", "")).strip()
        },
        key=str.casefold,
    )
    inputs = [
        ("manifest", manifest_path),
        ("official-route-report", route_report_path),
        ("stage-a-all-report", stage_a_report_path),
        ("source:stage-b-audit", tool_path),
        (
            "source:affine-embeddings",
            package_root / "isocore" / "distortion" / "affine_embeddings.py",
        ),
        (
            "source:inverse-landau",
            package_root / "isocore" / "distortion" / "inverse_landau.py",
        ),
        (
            "source:inverse-landau-adapter",
            package_root / "isocore" / "distortion" / "inverse_landau_adapter.py",
        ),
        ("source:lattice", package_root / "isocore" / "utils" / "lattice.py"),
        (
            "source:config-loader",
            package_root / "isocore" / "utils" / "config_loader.py",
        ),
        ("source:cif-reader", package_root / "isocore" / "structure" / "cif_io.py"),
        ("config", package_root / "config" / "settings.yaml"),
        *[(f"parent-cif:{name}", input_root / name) for name in parent_names],
    ]
    records: list[dict[str, Any]] = []
    for label, path in inputs:
        resolved = path.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"missing signed audit input {label}: {resolved}")
        records.append(
            {
                "label": label,
                "path": str(resolved),
                "sha256": _sha256(resolved),
                "size_bytes": resolved.stat().st_size,
            }
        )
    signature_payload = {
        "report_schema": SCHEMA,
        "checkpoint_schema": CHECKPOINT_SCHEMA,
        "algorithm": "inverse-landau-fixed-space-v1",
        "max_quotient_order": max_quotient_order,
        "case_ids": case_ids,
        "runtime_versions": runtime_versions,
        "scientific_tolerances": tolerances,
        "inputs": records,
    }
    signature = hashlib.sha256(
        json.dumps(
            signature_payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8"),
        usedforsecurity=False,
    ).hexdigest()
    return {"signature": signature, **signature_payload}


def _identity_key(identity: dict[str, Any]) -> tuple[object, ...]:
    return (
        int(identity["space_group_number"]),
        tuple(tuple(str(value) for value in row) for row in identity["basis"]),
        tuple(str(value) for value in identity.get("origin", (0, 0, 0))),
    )


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
        "cases": {},
    }


def _load_checkpoint(
    path: Path,
    signature_state: dict[str, Any],
    *,
    restart: bool,
) -> dict[str, Any]:
    if restart or not path.is_file():
        return _new_checkpoint(signature_state)
    payload = _load_json(path, "Stage-B checkpoint")
    if (
        payload.get("schema") != CHECKPOINT_SCHEMA
        or payload.get("run_signature") != signature_state["signature"]
    ):
        raise RuntimeError(
            "Stage-B checkpoint signature is stale; rerun with --restart "
            "instead of mixing results from different sources/configuration"
        )
    if not isinstance(payload.get("cases"), dict):
        raise ValueError("Stage-B checkpoint cases must be a JSON object")
    return payload


def _route_case_map(route_report: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    result: dict[tuple[str, str], dict[str, Any]] = {}
    for case in route_report.get("cases", []):
        key = (str(case.get("parent_cif", "")), str(case.get("case_id", "")))
        if key in result:
            raise ValueError(f"duplicate route-report case {key!r}")
        result[key] = case
    return result


def _audit_case(
    stage_case: dict[str, Any],
    route_case: dict[str, Any],
    *,
    structure: Any,
    parent: Any,
    max_quotient_order: int,
) -> dict[str, Any]:
    started = time.perf_counter()
    parent_cif = str(stage_case["parent_cif"])
    case_id = str(stage_case["case_id"])
    if stage_case.get("status") != "completed":
        raise ValueError(f"Stage-A case {case_id} is not completed")
    if route_case.get("status") != "audited":
        raise ValueError(f"route case {case_id} is not audited")

    candidate_rows = list(stage_case.get("stage_a_candidates", []))
    if not candidate_rows:
        raise ValueError(f"Stage-A case {case_id} has no candidates")
    embeddings = [
        embedding_from_identity(candidate["identity"], parent)
        for candidate in candidate_rows
    ]
    reference_lattice = embeddings[0].lattice
    if any(not same_lattice(item.lattice, reference_lattice) for item in embeddings[1:]):
        raise ValueError(f"Stage-A case {case_id} spans multiple translation lattices")

    build_started = time.perf_counter()
    quotient = build_affine_quotient(
        parent,
        reference_lattice,
        max_quotient_order=max_quotient_order,
    )
    selected_types = tuple(str(item) for item in stage_case["query"]["distortion_types"])
    bundle = build_selected_character_representation(
        parent,
        quotient,
        structure,
        selected_types,
    )
    if bundle.unsupported_types:
        raise ValueError(
            f"case {case_id} contains unsupported types: {bundle.unsupported_types!r}"
        )
    if bundle.representation is None:
        raise ValueError(f"case {case_id} produced no physical representation")
    feasibility_analyzer = FixedSpaceFeasibilityAnalyzer(bundle.representation)
    build_seconds = time.perf_counter() - build_started

    matching = stage_case["maximum_one_to_one_matching"]
    exact_candidate_indices = {
        int(match["candidate_index"]) for match in matching.get("matches", [])
    }
    candidate_to_orbit: dict[int, int] = {}
    for orbit in stage_case.get("stage_a_parent_affine_orbits", []):
        for candidate_index in orbit["candidate_indices"]:
            candidate_to_orbit[int(candidate_index)] = int(orbit["orbit_index"])

    candidate_results: list[dict[str, Any]] = []
    analyze_started = time.perf_counter()
    for candidate_index, (candidate, embedding) in enumerate(
        zip(candidate_rows, embeddings, strict=True)
    ):
        subgroup = locate_embedding_subgroup(embedding, quotient)
        fixed_space = feasibility_analyzer.analyze(subgroup)
        has_official_equivalent = bool(candidate["has_any_official_equivalent"])
        if candidate_index in exact_candidate_indices:
            relationship = "official_exact_representative"
        elif has_official_equivalent:
            relationship = "official_parent_conjugate_multiplicity"
        else:
            relationship = "no_official_affine_equivalent"
        expected_reachable = has_official_equivalent
        candidate_results.append(
            {
                "candidate_index": candidate_index,
                "orbit_index": candidate_to_orbit[candidate_index],
                "identity": candidate["identity"],
                "embedding_id": stable_embedding_id(
                    parent,
                    embedding,
                    quotient=quotient,
                ),
                "relationship_to_official": relationship,
                "subgroup_order": len(subgroup),
                "fixed_dimension": fixed_space.fixed_dimension,
                "pointwise_stabilizer_indices": list(
                    fixed_space.pointwise_stabilizer
                ),
                "is_symmetry_reachable": fixed_space.is_reachable,
                "obstruction": fixed_space.obstruction,
                "expected_reachable": expected_reachable,
                "expectation_pass": (
                    fixed_space.is_reachable == expected_reachable
                ),
                "route_status": "unresolved_irrep_decomposition",
                "selectable": False,
            }
        )
    analyze_seconds = time.perf_counter() - analyze_started

    route_embeddings = list(route_case.get("embeddings", []))
    official_rows = list(stage_case.get("official_embeddings", []))
    if len(route_embeddings) != len(official_rows):
        raise ValueError(f"route/Stage-A official count mismatch for {case_id}")
    for official_index, (route_embedding, official) in enumerate(
        zip(route_embeddings, official_rows, strict=True)
    ):
        if _identity_key(route_embedding["official_identity"]) != _identity_key(official):
            raise ValueError(
                f"route/Stage-A official identity mismatch in {case_id} index {official_index}"
            )
    match_by_official = {
        int(match["official_index"]): int(match["candidate_index"])
        for match in matching.get("matches", [])
    }
    official_results = []
    for official_index, (route_embedding, official) in enumerate(
        zip(route_embeddings, official_rows, strict=True)
    ):
        candidate_index = match_by_official[official_index]
        candidate_result = candidate_results[candidate_index]
        official_results.append(
            {
                "official_index": official_index,
                "official_identity": official,
                "route_classification": route_embedding["classification"],
                "matched_candidate_index": candidate_index,
                "matched_embedding_id": candidate_result["embedding_id"],
                "is_symmetry_reachable": candidate_result["is_symmetry_reachable"],
                "route_status": "unresolved_irrep_decomposition",
                "selectable": False,
            }
        )

    orbit_results: list[dict[str, Any]] = []
    for orbit in stage_case.get("stage_a_parent_affine_orbits", []):
        candidate_indices = [int(value) for value in orbit["candidate_indices"]]
        reachable_values = {
            bool(candidate_results[index]["is_symmetry_reachable"])
            for index in candidate_indices
        }
        official_indices = [int(value) for value in orbit["exact_official_indices"]]
        expected_reachable = bool(official_indices)
        orbit_results.append(
            {
                "orbit_index": int(orbit["orbit_index"]),
                "candidate_indices": candidate_indices,
                "embedding_ids": [
                    candidate_results[index]["embedding_id"]
                    for index in candidate_indices
                ],
                "official_indices": official_indices,
                "candidate_count": len(candidate_indices),
                "reachability_uniform_under_parent_affine_conjugacy": (
                    len(reachable_values) == 1
                ),
                "is_symmetry_reachable": (
                    next(iter(reachable_values)) if len(reachable_values) == 1 else None
                ),
                "expected_reachable": expected_reachable,
                "expectation_pass": (
                    len(reachable_values) == 1
                    and next(iter(reachable_values)) == expected_reachable
                ),
            }
        )

    errors: list[str] = []
    if any(not result["expectation_pass"] for result in candidate_results):
        errors.append("one or more candidate feasibility results contradict official equivalence")
    if any(not result["is_symmetry_reachable"] for result in official_results):
        errors.append("one or more official embeddings are infeasible")
    if any(not result["expectation_pass"] for result in orbit_results):
        errors.append("parent-affine orbit feasibility is nonuniform or unexpected")
    if len({result["embedding_id"] for result in candidate_results}) != len(
        candidate_results
    ):
        errors.append("distinct Stage-A candidates collapsed to one stable embedding_id")

    return {
        "parent_cif": parent_cif,
        "case_id": case_id,
        "status": "completed" if not errors else "failed_acceptance",
        "errors": errors,
        "selected_types": list(selected_types),
        "cache_contract": {
            "quotient_build_count": 1,
            "representation_build_count": 1,
            "feasibility_analyzer_build_count": 1,
            "candidate_analysis_count": len(candidate_results),
        },
        "quotient_order": len(quotient.operations),
        "translation_index": quotient.translation_index,
        "representation_dimension": bundle.representation.dimension,
        "representation_verification_method": (
            bundle.representation.verification_method
        ),
        "feasibility_verification_method": (
            "exact_character_average_and_generated_subgroup_dimension"
        ),
        "component_names": list(bundle.component_names),
        "site_count": bundle.site_count,
        "max_site_match_error_angstrom": bundle.max_site_match_error_angstrom,
        "official_embedding_count": len(official_results),
        "stage_a_candidate_count": len(candidate_results),
        "parent_affine_orbit_count": len(orbit_results),
        "candidate_results": candidate_results,
        "official_results": official_results,
        "parent_affine_orbits": orbit_results,
        "timing_seconds": {
            "build_shared_quotient_and_representation": build_seconds,
            "analyze_candidates": analyze_seconds,
            "total": time.perf_counter() - started,
        },
    }


def _summary(cases: list[dict[str, Any]], *, full_scope: bool) -> dict[str, Any]:
    completed = [case for case in cases if case.get("status") == "completed"]
    candidate_results = [
        candidate for case in completed for candidate in case["candidate_results"]
    ]
    official_results = [
        official for case in completed for official in case["official_results"]
    ]
    orbit_results = [
        orbit for case in completed for orbit in case["parent_affine_orbits"]
    ]
    novel_orbits = [orbit for orbit in orbit_results if not orbit["official_indices"]]
    surplus = [
        candidate
        for candidate in candidate_results
        if candidate["relationship_to_official"]
        == "official_parent_conjugate_multiplicity"
    ]
    novel_by_case = {
        case["case_id"]: sum(
            1 for orbit in case.get("parent_affine_orbits", []) if not orbit["official_indices"]
        )
        for case in completed
        if any(not orbit["official_indices"] for orbit in case.get("parent_affine_orbits", []))
    }
    actual_counts = {
        "case_count": len(cases),
        "official_embedding_count": len(official_results),
        "stage_a_candidate_count": len(candidate_results),
        "parent_affine_orbit_count": len(orbit_results),
        "novel_orbit_count": len(novel_orbits),
    }
    full_count_contract = (
        actual_counts == EXPECTED_FULL_COUNTS if full_scope else None
    )
    full_novel_case_contract = (
        novel_by_case == EXPECTED_NOVEL_CASES if full_scope else None
    )
    official_reachable = sum(
        bool(item["is_symmetry_reachable"]) for item in official_results
    )
    novel_infeasible = sum(
        item["is_symmetry_reachable"] is False for item in novel_orbits
    )
    surplus_reachable = sum(
        bool(item["is_symmetry_reachable"]) for item in surplus
    )
    complete = (
        len(completed) == len(cases)
        and official_reachable == len(official_results)
        and novel_infeasible == len(novel_orbits)
        and surplus_reachable == len(surplus)
        and all(item["expectation_pass"] for item in candidate_results)
        and all(item["expectation_pass"] for item in orbit_results)
        and (not full_scope or bool(full_count_contract and full_novel_case_contract))
    )
    return {
        "complete": complete,
        "scope": "full-40-case" if full_scope else "filtered-partial",
        "selected_case_count": len(cases),
        "completed_case_count": len(completed),
        "error_or_failed_case_count": len(cases) - len(completed),
        **actual_counts,
        "official_reachable_count": official_reachable,
        "official_infeasible_count": len(official_results) - official_reachable,
        "candidate_reachable_count": sum(
            bool(item["is_symmetry_reachable"]) for item in candidate_results
        ),
        "candidate_infeasible_count": sum(
            item["is_symmetry_reachable"] is False for item in candidate_results
        ),
        "novel_orbit_infeasible_count": novel_infeasible,
        "surplus_parent_conjugate_candidate_count": len(surplus),
        "surplus_parent_conjugate_reachable_count": surplus_reachable,
        "novel_orbits_by_case": novel_by_case,
        "full_count_contract": full_count_contract,
        "full_novel_case_contract": full_novel_case_contract,
        "expected_full_counts": EXPECTED_FULL_COUNTS,
        "expected_novel_cases": EXPECTED_NOVEL_CASES,
    }


def run_audit(
    *,
    manifest_path: Path = MANIFEST_PATH,
    route_report_path: Path = ROUTE_REPORT_PATH,
    stage_a_report_path: Path = STAGE_A_REPORT_PATH,
    input_root: Path = INPUT_ROOT,
    report_path: Path = REPORT_PATH,
    checkpoint_path: Path = CHECKPOINT_PATH,
    package_root: Path = PACKAGE_ROOT,
    max_quotient_order: int = 512,
    case_ids: tuple[str, ...] | None = None,
    restart: bool = False,
) -> dict[str, Any]:
    """Run or resume the signed, case-cached Stage-B diagnostic audit."""
    started = time.perf_counter()
    manifest_path = manifest_path.resolve()
    route_report_path = route_report_path.resolve()
    stage_a_report_path = stage_a_report_path.resolve()
    input_root = input_root.resolve()
    package_root = package_root.resolve()
    runtime_versions = _runtime_versions()
    tolerances = _scientific_tolerances()
    signature_state = _signature_state(
        manifest_path=manifest_path,
        route_report_path=route_report_path,
        stage_a_report_path=stage_a_report_path,
        input_root=input_root,
        package_root=package_root,
        tool_path=Path(__file__).resolve(),
        runtime_versions=runtime_versions,
        tolerances=tolerances,
        max_quotient_order=max_quotient_order,
        case_ids=case_ids,
    )
    route_report = _load_json(route_report_path, "official route report")
    stage_a_report = _load_json(stage_a_report_path, "Stage-A all report")
    if not bool(stage_a_report.get("summary", {}).get("complete")):
        raise ValueError("Stage-A all report is not complete")
    if int(route_report.get("summary", {}).get("audited_embedding_count", -1)) != 77:
        raise ValueError("official route report does not contain the frozen 77 embeddings")
    route_cases = _route_case_map(route_report)
    selected = [
        case
        for case in stage_a_report.get("cases", [])
        if case_ids is None or str(case.get("case_id")) in case_ids
    ]
    if case_ids is not None:
        missing = set(case_ids).difference(str(case.get("case_id")) for case in selected)
        if missing:
            raise ValueError(f"unknown requested Stage-A cases: {sorted(missing)!r}")
    checkpoint = _load_checkpoint(
        checkpoint_path,
        signature_state,
        restart=restart,
    )
    cached_cases: dict[str, Any] = checkpoint["cases"]
    initially_cached_keys = frozenset(cached_cases)
    parent_cache: dict[str, tuple[Any, Any]] = {}

    for index, stage_case in enumerate(selected, start=1):
        parent_cif = str(stage_case["parent_cif"])
        case_id = str(stage_case["case_id"])
        key = _case_key(parent_cif, case_id)
        if key in cached_cases:
            print(f"[{index}/{len(selected)}] resume {case_id}", flush=True)
            continue
        print(f"[{index}/{len(selected)}] audit {case_id}", flush=True)
        if parent_cif not in parent_cache:
            structure = read_cif(input_root / parent_cif)
            parent_cache[parent_cif] = (structure, parent_affine_group(structure))
        structure, parent = parent_cache[parent_cif]
        route_key = (parent_cif, case_id)
        try:
            route_case = route_cases[route_key]
            result = _audit_case(
                stage_case,
                route_case,
                structure=structure,
                parent=parent,
                max_quotient_order=max_quotient_order,
            )
        except Exception as exc:  # noqa: BLE001 - preserve batch evidence per case
            result = {
                "parent_cif": parent_cif,
                "case_id": case_id,
                "status": "error",
                "errors": [f"{type(exc).__module__}.{type(exc).__qualname__}: {exc}"],
            }
        cached_cases[key] = result
        checkpoint["updated_at"] = _now()
        _atomic_json(checkpoint_path, checkpoint)

    ordered_cases = [
        cached_cases[_case_key(str(case["parent_cif"]), str(case["case_id"]))]
        for case in selected
    ]
    summary = _summary(ordered_cases, full_scope=case_ids is None)
    report = {
        "schema": SCHEMA,
        "run_signature": signature_state["signature"],
        "generated_at": _now(),
        "tool": str(Path(__file__).resolve()),
        "scientific_scope": {
            "finite_group": "Q = N_G(T_s) / T_s",
            "projector_formula": "P_H = (1/|H|) * sum_{h in H} D(h)",
            "fixed_dimension_formula": (
                "dim Fix(H) = trace(P_H) = (1/|H|) * "
                "sum_{h in H} trace(D(h))"
            ),
            "pointwise_stabilizer_test": (
                "g fixes every vector in Fix(H) iff "
                "dim Fix(<H,g>) == dim Fix(H)"
            ),
            "feasibility_criterion": (
                "Fix(H) != {0} and no q outside H acts identically on all Fix(H)"
            ),
            "selected_representation": (
                "homogeneous strain direct-summed with finite-supercell internal "
                "atomic displacements; uniform rigid translation is quotiented out"
            ),
            "exact_representation_verification": (
                "Arbitrary supplied matrices receive a full exact multiplication-table "
                "homomorphism check. Structured representations use exact verified site "
                "permutations/vector actions plus functorial Sym^2 and exact character "
                "identities for tensor products, direct sums, and the uniform-vector "
                "invariant quotient."
            ),
            "claim_boundary": (
                "Embedding feasibility only. IR decomposition and coupled routes "
                "remain unresolved; every row is selectable=false."
            ),
            "parent_affine_multiplicity": (
                "All 142 candidates are analysed. Parent-affine orbit members must "
                "have uniform feasibility, but remain distinct embedding IDs."
            ),
        },
        "signature_state": signature_state,
        "runtime_versions": runtime_versions,
        "scientific_tolerances": tolerances,
        "checkpoint": {
            "path": str(checkpoint_path.resolve()),
            "resumed_case_count": sum(
                1
                for case in selected
                if _case_key(str(case["parent_cif"]), str(case["case_id"]))
                in initially_cached_keys
            ),
        },
        "summary": summary,
        "cases": ordered_cases,
        "elapsed_seconds": time.perf_counter() - started,
    }
    _atomic_json(report_path, report)
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--route-report", type=Path, default=ROUTE_REPORT_PATH)
    parser.add_argument("--stage-a-report", type=Path, default=STAGE_A_REPORT_PATH)
    parser.add_argument("--input-root", type=Path, default=INPUT_ROOT)
    parser.add_argument("--report", type=Path, default=REPORT_PATH)
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT_PATH)
    parser.add_argument("--max-quotient-order", type=int, default=512)
    parser.add_argument("--case", action="append", dest="case_ids")
    parser.add_argument(
        "--restart",
        action="store_true",
        help="discard a compatible/incompatible checkpoint and recompute selected cases",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    report = run_audit(
        manifest_path=args.manifest,
        route_report_path=args.route_report,
        stage_a_report_path=args.stage_a_report,
        input_root=args.input_root,
        report_path=args.report,
        checkpoint_path=args.checkpoint,
        max_quotient_order=args.max_quotient_order,
        case_ids=None if args.case_ids is None else tuple(args.case_ids),
        restart=args.restart,
    )
    print(json.dumps(report["summary"], ensure_ascii=False, indent=2))
    return 0 if report["summary"]["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
