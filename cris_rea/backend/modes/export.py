"""JSON, Cartesian mode CSV, and H-asymmetric-unit CIF exports."""
import csv
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from .core import seitz_expression, ModeError
from ..structure import parse_cif, expand_structure, StructureError


def json_value(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    return value


def write_json(path, data):
    Path(path).write_text(json.dumps(json_value(data), ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def _quote(value):
    # Source labels are not used as CIF identifiers. Species from CIF may need quoting.
    if "'" not in value:
        return "'" + value + "'"
    if '"' not in value:
        return '"' + value + '"'
    raise ValueError("CIF value contains both quote kinds and is not supported for export.")


def cif_text(model, state, subgroup, title, *, full_p1=False):
    a = state["lattice"]
    lengths = np.linalg.norm(a, axis=0)
    angles = [math.degrees(math.acos(float(np.clip(
        np.dot(a[:, i], a[:, j]) / (lengths[i]*lengths[j]), -1, 1))))
        for i, j in ((1, 2), (0, 2), (0, 1))]
    lines = [f"data_{title}", "# H-preserving structure; special amplitudes may have higher symmetry.",
             "# Cartesian amplitude norm refers to the unstrained conventional subgroup cell.",
             f"_cris_rea_preserved_subgroup_number {subgroup.subgroup_number}",
             f"_space_group_IT_number {1 if full_p1 else subgroup.subgroup_number}",
             "_space_group_name_H-M_alt " + _quote("P 1" if full_p1 else subgroup.subgroup_symbol)]
    if full_p1:
        lines += ["# Full-cell P1 serialization avoids decimal rounding of H-special-site constraints.",
                  "# The verified H operations and exact reference positions remain in modes.json."]
    lines += [f"_cell_length_{axis} {value:.15g}" for axis, value in zip("abc", lengths)]
    lines += [f"_cell_angle_{name} {value:.15g}"
              for name, value in zip(("alpha", "beta", "gamma"), angles)]
    lines += ["loop_", "_space_group_symop_id", "_space_group_symop_operation_xyz"]
    lines += (["1 'x,y,z'"] if full_p1 else
              [f"{i+1} '{seitz_expression(op)}'" for i, op in enumerate(model["operations"])])
    lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol",
              "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
              "_atom_site_occupancy", "_atom_site_site_symmetry_multiplicity"]
    visited = set()
    for index, atom in enumerate(model["atoms"]):
        if index in visited:
            continue
        orbit = {index} if full_p1 else {p[index] for p in model["permutations"]}
        visited.update(orbit)
        xyz = " ".join(f"{v:.15g}" for v in state["fractional"][index])
        lines.append(f"{atom['label']} {_quote(atom['symbol'])} {xyz} 1 {len(orbit)}")
    return "\n".join(lines) + "\n"


def checked_cif(model, state, subgroup, title):
    """Validate the serialization; use full P1 if decimal ASU cannot preserve H constraints."""
    def check(text):
        parsed = parse_cif(text)
        expanded = expand_structure(parsed)
        if len(expanded) != len(model["atoms"]):
            raise ValueError("CIF expansion changed the atom count.")
        # 15 significant decimal digits -> <= 5e-15 fractional component rounding.
        # 1e-12 permits accumulation through at most crystallographic integer ops.
        # This checks serialization only; never changes coordinates or scientific tolerances.
        remaining = set(range(len(expanded)))
        error = 0.
        for atom, position in zip(model["atoms"], state["fractional"]):
            candidates = []
            for j in remaining:
                if expanded[j].type_symbol != atom["symbol"]:
                    continue
                delta = np.asarray(expanded[j].fractional, float) - position
                delta -= np.rint(delta)
                candidates.append((float(np.max(np.abs(delta))), j))
            if not candidates:
                raise ValueError("CIF expansion changed species.")
            distance, j = min(candidates)
            if distance > 1e-12:
                raise ValueError("CIF serialization exceeds its fractional rounding bound.")
            remaining.remove(j)
            error = max(error, distance)
        source_metric = state["lattice"].T@state["lattice"]
        metric = np.asarray(parsed.lattice.metric)
        scale = np.sqrt(np.outer(np.diag(source_metric), np.diag(source_metric)))
        metric_error = float(np.max(np.abs(metric-source_metric)/scale))
        if metric_error > 1e-12:
            raise ValueError("CIF cell roundtrip exceeds its rounding bound.")
        return {"atom_count": len(expanded), "fractional_rounding_residual": error,
                "metric_rounding_residual": metric_error, "rounding_bound": 1e-12}
    text = cif_text(model, state, subgroup, title)
    try:
        info = check(text)
        info["representation"] = "H_asymmetric_unit"
        return text, info
    except (StructureError, ValueError) as exc:
        reason = str(exc)
    text = cif_text(model, state, subgroup, title, full_p1=True)
    try:
        info = check(text)
    except (StructureError, ValueError) as exc:
        raise ModeError("cif_roundtrip_failure", str(exc)) from exc
    info.update({"representation": "P1_full_cell", "reason": reason,
                 "preserved_subgroup_number": subgroup.subgroup_number})
    return text, info


def export_modes(directory, model, state, subgroup, provenance):
    directory = Path(directory)
    directory.mkdir()
    names = [f"internal_{i+1}" for i in range(model["internal_modes"])]
    names += [f"translation_{i+1}" for i in range(model["rigid_translations"])]
    payload = {
        "schema_version": 1, "basis_convention": model["basis_convention"],
        "scope": "Complete H-invariant polar displacement and homogeneous symmetric strain spaces "
                 "on the stated finite cell; columns are not parent-IR-resolved.",
        "coordinates": {"lattice_columns": "Cartesian angstrom", "fractional": "subgroup conventional cell",
                        "basis_columns": "ISO vectors in parent conventional fractional coordinates",
                        "origin": "parent conventional fractional coordinates"},
        "normalization": {
            "displacement": "sum over all atoms in subgroup conventional cell of |e_i|^2 = 1",
            "amplitude_unit": "angstrom; unstrained Cartesian RSS",
            "mass_weighted": False,
            "parent_primitive_cells_in_export_cell": model["parent_primitive_cells"],
            "rss_per_parent_primitive_cell": "norm(amplitudes)/sqrt(parent_primitive_cells_in_export_cell)",
            "strain": "Kelvin/Frobenius orthonormal; dimensionless; F=I+E positive definite",
            "website_column_equivalence": "not established"},
        "subgroup": {"number": subgroup.subgroup_number, "symbol": subgroup.subgroup_symbol,
                     "basis": subgroup.basis, "origin": subgroup.origin,
                     "index": subgroup.index, "primitive_size": subgroup.size},
        "lattice": model["lattice"], "atoms": model["atoms"],
        "operation_xyz": [seitz_expression(op) for op in model["operations"]],
        "mode_names": names, "displacement_basis": model["basis"],
        "basis_shape": "rows = 3*atom_index + Cartesian component, columns = mode",
        "strain_basis_kelvin": model["strain_basis"],
        "kelvin_order": ["xx", "yy", "zz", "sqrt(2)*yz", "sqrt(2)*xz", "sqrt(2)*xy"],
        "selected_state": state, "validation": model["validation"], "provenance": provenance}
    write_json(directory / "modes.json", payload)
    with (directory / "modes.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mode", "atom", "species", "source_label", "dx", "dy", "dz",
                         "normalization", "amplitude_unit"])
        for j, name in enumerate(names):
            for i, atom in enumerate(model["atoms"]):
                writer.writerow([name, atom["label"], atom["symbol"], atom["source_label"],
                                 *[f"{v:.17g}" for v in model["basis"][3*i:3*i+3, j]],
                                 "supercell_Cartesian_RSS_1", "angstrom"])
    reference = {"lattice": model["lattice"],
                 "fractional": np.array([a["fractional_exact"] for a in model["atoms"]], float)}
    serialization = {}
    for name, data in (("reference", reference), ("distorted", state)):
        text, serialization[name] = checked_cif(model, data, subgroup, name)
        (directory / f"{name}.cif").write_text(text, encoding="utf-8")
    write_json(directory / "report.json", {
        "status": "validated_local_subspace", "validation": model["validation"],
        "scope": payload["scope"], "normalization": payload["normalization"],
        "provenance": provenance, "cif_serialization": serialization,
        "limitations": ["No parent-IR-resolved columns or ISODISTORT column matching.",
                        "No phonon frequencies, force constants, or mass-weighted normal modes.",
                        "No magnetic, occupational, incommensurate or modulated input.",
                        "Subgroup setting must match the official ISO conventional setting.",
                        "Declared H is preserved; maximal symmetry at special amplitudes is not identified."]})
