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
    # Keep the legacy archive folder stable: the official run inside it is
    # authoritative evidence that this case is a homogeneous-strain success.
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
        PACKAGE_ROOT / "isocore" / "backend" / "iso_wrapper.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "distortion_engine.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "distortion_mapper.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "search_methods.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "strain.py",
        PACKAGE_ROOT / "isocore" / "distortion" / "superspace.py",
        PACKAGE_ROOT / "isocore" / "structure" / "cif_io.py",
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


def _homogeneous_strain(source: Structure) -> Structure:
    matrix = np.asarray(source.lattice.matrix, dtype=float).copy()
    matrix[0] *= 1.08
    return _copy_structure(source, lattice=Lattice(matrix))


def _threshold_failure(source: Structure) -> Structure:
    coords = np.asarray(source.frac_coords, dtype=float).copy()
    coords[0] += [0.36, 0.31, 0.0]
    return _copy_structure(source, coords=coords)


def _prepare_context(
    iso: IsoDistort,
    spec: dict[str, Any],
    kind: str,
) -> tuple[list[str], dict[str, Any]]:
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
    fields = candidate.official_fields()
    context = {
        "candidate_identity": {
            key: fields[key]
            for key in (
                "irrep", "opd", "dir", "sg", "basis", "origin", "s", "i",
                "k_active",
            )
        },
        "basis_vectors": candidate.basis_vectors,
        "origin": candidate.origin,
        "space_group_number": candidate.space_group_number,
        "space_group_symbol": candidate.space_group_symbol,
        "subgroup_index": candidate.subgroup_index,
        "size": candidate.size,
        "k_point_label": candidate.k_point_label,
        "k_parameters": list(candidate.k_parameters or []),
        "k_coordinates": list(candidate.k_coordinates or []),
        "mode_labels": labels,
    }
    return labels, context


def _generated(iso: IsoDistort, contributions: dict[str, float]) -> Structure:
    with contextlib.redirect_stdout(io.StringIO()):
        result = iso.generate_mixed_distortion(contributions)
    return result.copy()


def _scaled_contributions_for_dmax(
    iso: IsoDistort,
    zero: Structure,
    relative_contributions: dict[str, float],
    target_dmax_angstrom: float,
) -> tuple[dict[str, float], Structure]:
    """Scale a mode mixture to a material-independent maximum displacement."""
    # Probe inside the linear, unwrapped regime.  A unit raw coefficient can
    # move atoms across periodic boundaries; subtracting those wrapped
    # coordinates would then underestimate the true per-coefficient vector.
    probe_scale = 1e-4
    probe_contributions = {
        label: float(value * probe_scale)
        for label, value in relative_contributions.items()
    }
    probe = _generated(iso, probe_contributions)
    displacements = np.asarray(probe.cart_coords) - np.asarray(zero.cart_coords)
    probe_dmax = float(np.max(np.linalg.norm(displacements, axis=1)))
    if not np.isfinite(probe_dmax) or probe_dmax <= 0:
        raise RuntimeError("parameter-mode mixture has no finite nonzero displacement")
    scale = target_dmax_angstrom / probe_dmax
    contributions = {
        label: float(value * scale)
        for label, value in probe_contributions.items()
    }
    return contributions, _generated(iso, contributions)


def _cases(
    iso: IsoDistort,
    spec: dict[str, Any],
    seed: int,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    with contextlib.redirect_stdout(io.StringIO()):
        gamma_labels, gamma_context = _prepare_context(iso, spec, "gamma")
    g1, g2 = gamma_labels[:2]
    gamma_zero = _generated(iso, {g1: 0.0})
    gamma_pos = _generated(iso, {g1: 0.08})
    gamma_neg = _generated(iso, {g1: -0.08})
    gamma_mix = _generated(iso, {g1: 0.06, g2: -0.035})
    gamma_bound = _generated(iso, {g1: 0.18, g2: 0.11})
    origin_shift_child = [0.125, 0.25, 0.0]
    # The generated daughter coordinates are shifted by +t_child.  The
    # website field asks for the conventional child-cell origin relative to
    # the parent origin, so the explicit origin is the inverse translation
    # expressed in parent fractional coordinates: -t_child @ B.
    origin_shift_parent = (
        -np.asarray(origin_shift_child, dtype=float)
        @ np.asarray(gamma_context["basis_vectors"], dtype=float)
    ).tolist()

    with contextlib.redirect_stdout(io.StringIO()):
        parameter_labels, parameter_context = _prepare_context(iso, spec, "parameter")
    p1, p2 = parameter_labels[:2]
    parameter_zero = _generated(iso, {p1: 0.0})
    parameter_pos_contributions, parameter_pos = _scaled_contributions_for_dmax(
        iso, parameter_zero, {p1: 1.0}, 0.30,
    )
    parameter_mix_contributions, parameter_mix = _scaled_contributions_for_dmax(
        iso, parameter_zero, {p1: 5.0, p2: -3.0}, 0.35,
    )

    cases = [
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
            "structure": _shifted(gamma_mix, origin_shift_child),
            "context": "gamma",
            "contributions": {g1: 0.06, g2: -0.035},
            "provided_origin_shift": origin_shift_child,
            "provided_origin_shift_coordinate_system": "daughter fractional",
            "official_origin_shift_parent_fractional": origin_shift_parent,
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
            "contributions": parameter_pos_contributions,
            "target_max_displacement_angstrom": 0.30,
            "expect": "success",
        },
        {
            "id": "P03-mixed-reordered-wrapped",
            "structure": _wrapped(_reordered(parameter_mix)),
            "context": "parameter",
            "contributions": parameter_mix_contributions,
            "target_max_displacement_angstrom": 0.35,
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
            "structure": _homogeneous_strain(gamma_zero),
            "context": "gamma",
            "contributions": {g1: 0.0},
            "lattice_transform": "scale child direct-lattice vector 1 by 1.08",
            "expect_nonzero_strain": True,
            "expect": "success",
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
    return cases, {
        "gamma": gamma_context,
        "parameter": parameter_context,
    }


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
        "schema": 2,
        "purpose": "Frozen daughter-CIF upload batch for the future Method 4 differential",
        "signature": signature,
        "method4_executed": False,
        "official_download_status": "not_downloaded",
        "official_download_root": "output_compare/<parent>/官网/Method4/<readable-case-folder>/",
        "local_download_root": "output_compare/<parent>/现有网页版交互/Method4/<readable-case-folder>/",
        "readable_case_folders": READABLE_CASE_FOLDERS,
        "per_case_downloads": [
            "the exact uploaded daughter.cif whose SHA-256 matches this manifest",
            "parent-CIF identity",
            (
                "basis/origin/matching page HTML plus a screenshot or print-to-PDF "
                "captured after dynamic controls were set"
            ),
            "complete Method 4 result page (HTML)",
            "Complete modes details page (HTML) with the full official mode labels and order",
            "the individually offered subgroup.cif, topas.str, and data.isoviz exports",
            "amplitude table txt/csv when offered; otherwise preserve the full HTML table",
            "official error page/text for expected failures",
        ],
        "case_specific_downloads": {
            "G05-origin-shift": [
                (
                    "case root: explicit-origin run using "
                    "official_origin_shift_parent_fractional, with the six standard "
                    "artifacts, uploaded daughter copy, and completed-controls screenshot/PDF"
                ),
                (
                    "auto-origin/: automatic-origin run with the same six standard "
                    "artifact names and a completed-controls screenshot/PDF"
                ),
            ],
            "G06-near-bound-noise": [
                "compare amplitudes for uniquely matched shared target modes",
                (
                    "record nonzero extra official P1 modes; do not require official "
                    "full-P1 residual to equal the local subgroup-restricted residual"
                ),
            ],
            "F02-invalid-lattice": [
                "legacy folder name retained so the completed official archive is not moved or redownloaded",
                (
                    "positive homogeneous-strain case: compare the applied six-component "
                    "IsoVIZ tensor, not displayed symmetry-mode amplitudes alone"
                ),
            ],
        },
        "local_only_checks": [
            "calling Method 4 before selecting/calculating modes must be rejected",
        ],
        "local_validation_tolerances": {
            "exact_amplitude_abs": 2e-7,
            "exact_rms_residual_angstrom_max": 1e-6,
            "noise_amplitude_abs": 5e-4,
            "noise_rms_residual_angstrom_min": 1e-7,
            "noise_rms_residual_angstrom_max": 5e-3,
            "strain_component_abs": 5e-7,
            "strain_metric_relative_residual": 1e-10,
        },
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
                "parent_cif": (
                    WORKSPACE / "experiment_data" / parent_name
                ).relative_to(WORKSPACE).as_posix(),
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
        cases, resolved_contexts = _cases(iso, spec, 20260924 + parent_index)
        parent_entry["resolved_contexts"] = resolved_contexts
        for row in cases:
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
