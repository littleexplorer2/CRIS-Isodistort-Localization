"""Freeze reproducible Method 4 daughter-CIF inputs without running Method 4.

This is a preparation step for the future official-site differential.  It
creates twelve upload cases per parent from two real mode contexts and records
hashes after every file, making the run safe to resume after interruption.
It intentionally does not call ``search_method_4`` or compare amplitudes.
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

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from isocore.api import IsoDistort  # noqa: E402

OUTPUT_ROOT = PACKAGE_ROOT / "output" / "validation" / "method4_inputs"
MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method4_download_manifest.json"

READABLE_CASE_FOLDERS = {
    "G01-zero": "G01 - gamma - zero",
    "G02-positive-reordered": "G02 - gamma - positive - reordered",
    "G03-negative-wrapped": "G03 - gamma - negative - wrapped",
    "G04-two-mode": "G04 - gamma - two-mode",
    "G05-origin-shift": "G05 - gamma - origin-shift",
    "G06-near-bound-noise": "G06 - gamma - near-bound-noise",
    "P01-zero-supercell": "P01 - parameter-k - zero-supercell",
    "P02-positive-supercell": "P02 - parameter-k - positive-supercell",
    "P03-mixed-reordered-wrapped": "P03 - parameter-k - mixed-reordered-wrapped",
    "F01-species-mismatch": "F01 - expected-failure - species-mismatch",
    "F02-invalid-lattice": "F02 - expected-failure - invalid-lattice",
    "F03-distance-threshold": "F03 - expected-failure - distance-threshold",
}

PARENTS = {
    "EuAl4 Parent.cif": {
        "slug": "EuAl4",
        "gamma": {"sg": 12, "irrep": "GM5+", "opd": "P1"},
        "parameter": {
            "sg": 99,
            "irrep": "LD1",
            "opd": "C1",
            "basis": [[1, 0, 0], [0, 1, 0], [0, 0, 6]],
        },
    },
    "NdNiO2 own.cif": {
        "slug": "NdNiO2",
        "gamma": {"sg": 6, "irrep": "GM5-", "opd": "C1"},
        "parameter": {
            "sg": 47,
            "irrep": "Y1",
            "opd": "P1",
            "basis": [[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        },
    },
}


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temp, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _signature() -> str:
    digest = hashlib.sha256()
    for path in [
        Path(__file__),
        PACKAGE_ROOT / "isocore" / "api" / "core_api.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "distortion_mapper.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "search_methods.py",
        WORKSPACE / "experiment_data" / "EuAl4 Parent.cif",
        WORKSPACE / "experiment_data" / "NdNiO2 own.cif",
    ]:
        digest.update(path.relative_to(WORKSPACE).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _copy_structure(
    source: Structure,
    *,
    coords: np.ndarray | None = None,
    species: list[Any] | None = None,
    lattice: Lattice | None = None,
) -> Structure:
    return Structure(
        lattice or source.lattice,
        species or list(source.species),
        np.asarray(source.frac_coords if coords is None else coords, dtype=float),
        coords_are_cartesian=False,
        to_unit_cell=False,
    )


def _reordered(source: Structure) -> Structure:
    order = list(reversed(range(len(source))))
    return _copy_structure(
        source,
        species=[source.species[i] for i in order],
        coords=np.asarray([source.frac_coords[i] for i in order]),
    )


def _wrapped(source: Structure) -> Structure:
    coords = np.asarray(source.frac_coords, dtype=float).copy()
    offsets = np.zeros_like(coords)
    offsets[::2, 0] = 1.0
    offsets[1::2, 1] = -1.0
    return _copy_structure(source, coords=coords + offsets)


def _shifted(source: Structure, shift: list[float]) -> Structure:
    return _copy_structure(source, coords=np.asarray(source.frac_coords) + np.asarray(shift))


def _noisy(source: Structure, seed: int, sigma: float = 2e-5) -> Structure:
    rng = np.random.default_rng(seed)
    coords = np.asarray(source.frac_coords) + rng.normal(0.0, sigma, (len(source), 3))
    return _copy_structure(source, coords=coords)


def _species_mismatch(source: Structure) -> Structure:
    species = list(source.species)
    species[0] = "H"
    return _copy_structure(source, species=species)


def _invalid_lattice(source: Structure) -> Structure:
    matrix = np.asarray(source.lattice.matrix, dtype=float).copy()
    matrix[0] *= 1.08
    return _copy_structure(source, lattice=Lattice(matrix))


def _threshold_failure(source: Structure) -> Structure:
    coords = np.asarray(source.frac_coords, dtype=float).copy()
    coords[0] += [0.36, 0.31, 0.0]
    return _copy_structure(source, coords=coords)


def _prepare_context(iso: IsoDistort, spec: dict[str, Any], kind: str) -> list[str]:
    target = spec[kind]
    if kind == "gamma":
        rows = iso.search_method_1(
            distortion_types=["strain", "displacive"],
            subgroup_space_group=target["sg"],
        )
    else:
        rows = iso.search_method_3(
            distortion_types=["strain", "displacive"],
            space_group_type=target["sg"],
            supercell_basis=target["basis"],
            direct_sublattice_centering="P",
            lattice_type="direct",
            generate_if_missing=False,
        )
    candidate = next(
        item.subgroup
        for item in rows
        if item.subgroup.irrep_label == target["irrep"]
        and item.subgroup.opd_symbol == target["opd"]
        and item.subgroup.space_group_number == target["sg"]
    )
    iso.search_method_2(
        candidate.index,
        distortion_type=["displacive"],
        number_of_independent_modulations=0,
        candidates=[item.subgroup for item in rows],
    )
    labels = list(iso.mode_displacements_sc or iso.mode_displacements)
    if len(labels) < 2:
        raise RuntimeError(f"{kind} context returned fewer than two displacement modes")
    return labels


def _generated(iso: IsoDistort, contributions: dict[str, float]) -> Structure:
    with contextlib.redirect_stdout(io.StringIO()):
        result = iso.generate_mixed_distortion(contributions)
    return result.copy()


def _cases(iso: IsoDistort, spec: dict[str, Any], seed: int) -> list[dict[str, Any]]:
    with contextlib.redirect_stdout(io.StringIO()):
        gamma_labels = _prepare_context(iso, spec, "gamma")
    g1, g2 = gamma_labels[:2]
    gamma_zero = _generated(iso, {g1: 0.0})
    gamma_pos = _generated(iso, {g1: 0.08})
    gamma_neg = _generated(iso, {g1: -0.08})
    gamma_mix = _generated(iso, {g1: 0.06, g2: -0.035})
    gamma_bound = _generated(iso, {g1: 0.18, g2: 0.11})

    with contextlib.redirect_stdout(io.StringIO()):
        parameter_labels = _prepare_context(iso, spec, "parameter")
    p1, p2 = parameter_labels[:2]
    parameter_zero = _generated(iso, {p1: 0.0})
    parameter_pos = _generated(iso, {p1: 0.07})
    parameter_mix = _generated(iso, {p1: 0.05, p2: -0.03})

    return [
        {
            "id": "G01-zero",
            "structure": gamma_zero,
            "context": "gamma",
            "contributions": {g1: 0.0},
            "expect": "success",
        },
        {
            "id": "G02-positive-reordered",
            "structure": _reordered(gamma_pos),
            "context": "gamma",
            "contributions": {g1: 0.08},
            "transform": "reverse atom rows",
            "expect": "success",
        },
        {
            "id": "G03-negative-wrapped",
            "structure": _wrapped(gamma_neg),
            "context": "gamma",
            "contributions": {g1: -0.08},
            "transform": "integer periodic images",
            "expect": "success",
        },
        {
            "id": "G04-two-mode",
            "structure": gamma_mix,
            "context": "gamma",
            "contributions": {g1: 0.06, g2: -0.035},
            "expect": "success",
        },
        {
            "id": "G05-origin-shift",
            "structure": _shifted(gamma_mix, [0.125, 0.25, 0.0]),
            "context": "gamma",
            "contributions": {g1: 0.06, g2: -0.035},
            "provided_origin_shift": [0.125, 0.25, 0.0],
            "expect": "success",
        },
        {
            "id": "G06-near-bound-noise",
            "structure": _noisy(gamma_bound, seed),
            "context": "gamma",
            "contributions": {g1: 0.18, g2: 0.11},
            "noise_sigma_fractional": 2e-5,
            "expect": "success_with_nonzero_residual",
        },
        {
            "id": "P01-zero-supercell",
            "structure": parameter_zero,
            "context": "parameter",
            "contributions": {p1: 0.0},
            "expect": "success",
        },
        {
            "id": "P02-positive-supercell",
            "structure": parameter_pos,
            "context": "parameter",
            "contributions": {p1: 0.07},
            "expect": "success",
        },
        {
            "id": "P03-mixed-reordered-wrapped",
            "structure": _wrapped(_reordered(parameter_mix)),
            "context": "parameter",
            "contributions": {p1: 0.05, p2: -0.03},
            "transform": "reverse atom rows and integer periodic images",
            "expect": "success",
        },
        {
            "id": "F01-species-mismatch",
            "structure": _species_mismatch(gamma_zero),
            "context": "gamma",
            "expect": "reject_species_or_stoichiometry",
        },
        {
            "id": "F02-invalid-lattice",
            "structure": _invalid_lattice(gamma_zero),
            "context": "gamma",
            "expect": "reject_incompatible_lattice",
        },
        {
            "id": "F03-distance-threshold",
            "structure": _threshold_failure(gamma_zero),
            "context": "gamma",
            "atom_matching_method": "robust",
            "robust_distance_threshold": 0.1,
            "expect": "reject_unmatched_atom",
        },
    ]


def _write_cif(path: Path, structure: Structure) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".tmp.cif")
    CifWriter(structure, symprec=None).write_file(str(temp))
    os.replace(temp, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing batch whose signature differs",
    )
    args = parser.parse_args()
    signature = _signature()
    manifest: dict[str, Any] = {
        "schema": 1,
        "purpose": "Frozen daughter-CIF upload batch for the future Method 4 differential",
        "signature": signature,
        "method4_executed": False,
        "official_download_status": "not_downloaded",
        "official_download_root": "output_compare/<parent>/官网/Method4/<readable-case-folder>/",
        "local_download_root": "output_compare/<parent>/现有网页版交互/Method4/<readable-case-folder>/",
        "readable_case_folders": READABLE_CASE_FOLDERS,
        "per_case_downloads": [
            "the exact uploaded daughter.cif whose SHA-256 matches this manifest",
            "parent-CIF identity and Complete modes details for the selected path/subgroup context",
            "complete Method 4 result page (HTML)",
            "amplitude table txt/csv when offered; otherwise preserve the full HTML table",
            "official error page/text for expected failures",
        ],
        "local_only_checks": [
            "calling Method 4 before selecting/calculating modes must be rejected",
        ],
        "parents": {},
    }
    if MANIFEST_PATH.exists():
        previous = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        if previous.get("signature") != signature and not args.force:
            raise RuntimeError(
                "An input batch from a different source signature already exists; "
                "review it and pass --force only when deliberate replacement is intended"
            )
        if previous.get("signature") == signature:
            manifest = previous

    for parent_index, (parent_name, spec) in enumerate(PARENTS.items(), start=1):
        parent_entry = manifest["parents"].setdefault(
            parent_name,
            {
                "parent_cif": str((WORKSPACE / "experiment_data" / parent_name).relative_to(WORKSPACE)),
                "gamma_context": spec["gamma"],
                "parameter_context": spec["parameter"],
                "cases": {},
            },
        )
        existing = parent_entry["cases"]
        if len(existing) == 12 and all(
            (WORKSPACE / row["daughter_cif"]).is_file()
            and _sha256(WORKSPACE / row["daughter_cif"]) == row["sha256"]
            for row in existing.values()
        ):
            print(f"Method 4 inputs: {parent_name} already complete", flush=True)
            continue

        iso = IsoDistort()
        with contextlib.redirect_stdout(io.StringIO()):
            iso.load_structure(WORKSPACE / "experiment_data" / parent_name)
        for row in _cases(iso, spec, 20260924 + parent_index):
            case_id = row.pop("id")
            structure = row.pop("structure")
            destination = OUTPUT_ROOT / spec["slug"] / case_id / "daughter.cif"
            _write_cif(destination, structure)
            record = {
                **row,
                "daughter_cif": str(destination.relative_to(WORKSPACE)).replace("\\", "/"),
                "sha256": _sha256(destination),
                "atom_count": len(structure),
                "formula": structure.composition.formula,
            }
            existing[case_id] = record
            manifest["updated_at"] = datetime.now(UTC).isoformat()
            _atomic_json(MANIFEST_PATH, manifest)
            print(f"Method 4 input {parent_name}: {case_id}", flush=True)

    manifest["status"] = "inputs_frozen"
    manifest["case_count"] = sum(len(row["cases"]) for row in manifest["parents"].values())
    manifest["finished_at"] = datetime.now(UTC).isoformat()
    _atomic_json(MANIFEST_PATH, manifest)
    print(MANIFEST_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
