"""Write ISODISTORT-style ``data.isoviz`` (!tag layout for IsoVIZ).

Matches the official Save interactive distortion ascii layout closely enough
for IsoVIZ / ISOVIZ_INPUT amplitude patching. Atom orbits and mode vectors are
built from the local supercell + mode displacements when available; they are
not a byte-for-byte copy of the website file.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SpaceGroup

from ..utils.schoenflies import hm_symbol

_DEFAULT_ATOMIC_RADIUS = 0.4
_DEFAULT_BOND_MIN = 0.0
_DEFAULT_BOND_MAX = 2.5
_DEFAULT_APPLET_WIDTH = 1024


def render_isodistort_isoviz(spec: Any) -> str:
    """Return official-style IsoVIZ text for ``spec``."""
    sg = spec.subgroup
    parent = spec.parent_structure
    child = spec.structure
    parent_sg = int(spec.parent_sg or 0)
    parent_short = (spec.parent_symbol or hm_symbol(parent_sg) or "").strip()
    child_n = int(sg.space_group_number or 1)
    child_short = (sg.space_group_symbol or hm_symbol(child_n) or str(child_n)).strip()

    lines: list[str] = [
        "#isodistort_version_number ",
        "!isoversion 6.12",
        "",
        "#atom_sphere_radius_in_angstroms ",
        f"!atommaxradius    {_DEFAULT_ATOMIC_RADIUS:.5f}",
        "",
        "#angstroms_per_magneton ",
        "!angstromspermagneton    0.50000",
        "",
        "#angstroms_per_radian ",
        "!angstromsperradian    4.00000",
        "",
        "#default_isotropy_uiso_parameter ",
        "!defaultuiso    0.04000",
        "",
        "#maximum_bond_length_in_angstroms ",
        f"!maxbondlength    {_DEFAULT_BOND_MAX:.5f}",
        "",
        "#minimum_bond_length_in_angstroms ",
        f"!minbondlength    {_DEFAULT_BOND_MIN:.5f}",
        "",
        "#view_width_in_pixels ",
        f"!appletwidth   {_DEFAULT_APPLET_WIDTH}",
        "",
        "#number_of_incommensurate_modulations ",
        "!numberOfModulations      0",
        "",
    ]
    lines.extend(_space_group_block("parent", parent_sg, parent_short, parent))
    lines.extend(_basis_block(spec, parent, child))
    lines.extend(_space_group_block("child", child_n, child_short, child))
    lines.extend(_ops_block("child", child_n, child_short))
    lines.extend(_ops_block("parent", parent_sg, parent_short))
    lines.extend(_atom_blocks(spec, parent, child))
    lines.extend(_mode_blocks(spec))
    if spec.note:
        lines.append(f"# note: {spec.note}")
        lines.append("")
    return "\n".join(lines) + "\n"


def _space_group_block(
    role: str,
    number: int,
    short: str,
    structure: Structure | None,
) -> list[str]:
    full = _hm_quoted(number, short)
    hall = _hall_quoted(number)
    prefix = role
    lines = [
        f"#ITA number of {role} space group ",
        f"!{prefix}_space_group_number ",
        f"{number}",
        "",
        f"#ITA name of {role} space group ",
        f"!{prefix}_space_group_name_short ",
        f"{short:<10s}",
        "",
        f"#Hermann-Mauguin name of {role} space group ",
        f"!{prefix}_space_group_name ",
        full,
        "",
        f"#Hall symbol of {role} space group ",
        f"!{prefix}_space_group_Hall_symbol ",
        hall,
        "",
    ]
    if role == "parent" and structure is not None:
        lat = structure.lattice
        lines.extend(
            [
                "#parentcell_parameters_in_angstroms_and_radians ",
                "!parentcell ",
                "        "
                + "  ".join(
                    f"{v:10.5f}"
                    for v in (
                        lat.a,
                        lat.b,
                        lat.c,
                        math.radians(lat.alpha),
                        math.radians(lat.beta),
                        math.radians(lat.gamma),
                    )
                ),
                "",
            ]
        )
    return lines


def _basis_block(spec: Any, parent: Structure | None, child: Structure) -> list[str]:
    sg = spec.subgroup
    origin = np.asarray(sg.origin or [0.0, 0.0, 0.0], dtype=float)
    basis = np.asarray(sg.basis_vectors or np.eye(3), dtype=float)
    # parentbasis: parent vectors in child coordinates ≈ inv(basis)
    try:
        parent_in_child = np.linalg.inv(basis)
    except np.linalg.LinAlgError:
        parent_in_child = np.eye(3)
    lines = [
        "#parent_origin_in_supercell_units ",
        "!parentorigin ",
        "   " + "  ".join(f"{float(v):8.5f}" for v in origin),
        "",
        "#rows_are_conventional_parent_basis_vectors_relative_to_conventional_child_basis  ",
        "!parentbasis ",
    ]
    for row in parent_in_child:
        lines.append("  " + "  ".join(f"{float(v):8.5f}" for v in row))
    lines.append("")
    for tag, comment in (
        ("conv2primparentbasis", "primitive_parent"),
        ("conv2primchildbasis", "primitive_child"),
    ):
        lines.append(
            f"#rows_are_primitive_{comment}_basis_vectors_relativeto_conventional_"
            f"{'parent' if 'parent' in tag else 'child'}_basis "
        )
        lines.append(f"!{tag} ")
        # Identity fallback when centering unknown — IsoVIZ still loads structure.
        for row in np.eye(3):
            lines.append("  " + "  ".join(f"{float(v):8.5f}" for v in row))
        lines.append("")
    _ = parent, child
    return lines


def _ops_block(role: str, number: int, short: str) -> list[str]:
    ops, centering = _ops_and_centering(number, short)
    lines = [
        f"#rows_are_representative_4x4_affine_conventional_{role}_space_group_operations ",
        f"!conv{role}spacegroupops",
    ]
    for op in ops:
        lines.append(_affine16(op))
    lines.append("")
    lines.append(f"#rows_are_conventional_{role}_space_group_centering_translations ")
    lines.append(f"!conv{role}centeringvecs ")
    for vec in centering:
        lines.append("  " + "  ".join(f"{float(v):8.5f}" for v in vec))
    lines.append("")
    return lines


def _atom_blocks(spec: Any, parent: Structure | None, child: Structure) -> list[str]:
    src = parent if parent is not None else child
    # Unique element stems from parent ASU when possible.
    types: list[tuple[str, str]] = []
    seen: set[str] = set()
    for site in src:
        stem = _stem(site.label, site.species_string)
        if stem in seen:
            continue
        seen.add(stem)
        types.append((stem, site.species_string))
    lines = [
        "#parentatom/label/element ",
        "!atomtypelist ",
    ]
    for i, (stem, el) in enumerate(types, start=1):
        lines.append(f"   {i} {stem} {el} ")
    lines.append("")
    lines.append("#parentatom/subatom/label ")
    lines.append("!atomsubtypelist ")
    for i, (stem, _el) in enumerate(types, start=1):
        lines.append(f"   {i}   1 {stem}_1 ")
    lines.append("")
    lines.append("#parentatom/type/subatom/x/y/z/_for_each_subatom ")
    lines.append("!atomcoordlist ")
    # Emit every atom in the child supercell, mapped onto parent-type index by element.
    el_to_type = {el: i for i, (_stem, el) in enumerate(types, start=1)}
    counters: dict[int, int] = {}
    for site in child:
        tidx = el_to_type.get(site.species_string, 1)
        counters[tidx] = counters.get(tidx, 0) + 1
        x, y, z = (float(c) for c in site.frac_coords)
        lines.append(
            f"    {tidx:d}    1    {counters[tidx]}  {x:8.5f}  {y:8.5f}  {z:8.5f} "
        )
    lines.append("")
    n_atoms = len(child)
    lines.append("!atomsinunitcell ")
    for i in range(n_atoms):
        lines.append(f"  {1 if i < n_atoms else 0}")
    lines.append("")
    lines.append("!atomocclist ")
    for _ in range(n_atoms):
        lines.append("  1.00000")
    lines.append("")
    for empty in ("!atommaglist ", "!atomrotlist ", "!bondlist "):
        lines.append(empty)
        lines.append("")
    return lines


def _mode_blocks(spec: Any) -> list[str]:
    sg = spec.subgroup
    irrep = (sg.irrep_label or "IR").strip() or "IR"
    lines = [
        "!irreplist ",
        f"  1 {irrep:<10s}",
        "",
        "#strainmodenum/amp/maxamp/irrepnum/modelabel/modevector_for_each_mode ",
        "!strainmodelist ",
        "",
        "#parentatom/dispmodenum/amp/maxamp/irrepnum/modelabel/(modevector_for_each_subatom)_for_each_mode ",
        "!displacivemodelist ",
    ]
    disp = spec.mode_displacements_sc or {}
    labels = spec.mode_labels or {}
    amps = spec.amplitudes or {}
    if not disp:
        lines.append("")
        return lines
    for n, (key, arr) in enumerate(disp.items(), start=1):
        pretty = labels.get(key, key)
        amp = float(amps.get(key, 0.0))
        mat = np.asarray(arr, dtype=float)
        maxamp = float(np.max(np.linalg.norm(mat, axis=1))) if mat.size else 1.0
        if maxamp < 1e-12:
            maxamp = 1.0
        # Official header: parentatom_index, mode_index, amp, maxamp, irrepnum, label
        lines.append(f"    1    {n}  {amp:8.5f}  {maxamp:8.5f}    1 {pretty} ")
        for vec in mat:
            lines.append(
                "  "
                + "  ".join(f"{float(v):8.5f}" for v in vec)
            )
    lines.append("")
    return lines


def _hm_quoted(number: int, short: str) -> str:
    from .isodistort_cif import _HM_FULL_1992, _space_hm_symbol

    full = _HM_FULL_1992.get(int(number)) or _space_hm_symbol(short) or short
    return f'"{full}"'


def _hall_quoted(number: int) -> str:
    try:
        from .isodistort_cif import _hall_for_number

        _hall_n, hall = _hall_for_number(int(number))
        return f'"{hall}"'
    except Exception:  # noqa: BLE001
        return f'"{int(number):03d}"'


def _ops_and_centering(number: int, short: str) -> tuple[list[SymmOp], list[np.ndarray]]:
    try:
        from .isodistort_cif import _hall_for_number, _ops_from_hall

        hall_n, _hall = _hall_for_number(int(number))
        ops = _ops_from_hall(hall_n)
        if ops:
            return ops, _centering_from_symbol(short)
    except Exception:  # noqa: BLE001
        pass
    try:
        pmg = SpaceGroup(short.split()[0]) if short else SpaceGroup.from_int_number(int(number))
        ops = list(pmg.symmetry_ops)
        return ops, _centering_from_symbol(short or pmg.symbol)
    except Exception:  # noqa: BLE001
        op = SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3))
        return [op], [np.zeros(3)]


def _centering_from_symbol(symbol: str) -> list[np.ndarray]:
    letter = (symbol or "P")[:1].upper()
    zero = np.zeros(3)
    if letter == "I":
        return [zero, np.array([0.5, 0.5, 0.5])]
    if letter == "C":
        return [zero, np.array([0.5, 0.5, 0.0])]
    if letter == "A":
        return [zero, np.array([0.0, 0.5, 0.5])]
    if letter == "B":
        return [zero, np.array([0.5, 0.0, 0.5])]
    if letter == "F":
        return [
            zero,
            np.array([0.0, 0.5, 0.5]),
            np.array([0.5, 0.0, 0.5]),
            np.array([0.5, 0.5, 0.0]),
        ]
    return [zero]


def _affine16(op: SymmOp) -> str:
    rot = np.asarray(op.rotation_matrix, dtype=float)
    trans = np.mod(np.asarray(op.translation_vector, dtype=float), 1.0)
    mat = np.eye(4)
    mat[:3, :3] = rot
    mat[:3, 3] = trans
    vals = [f"{float(v):8.5f}" for v in mat.reshape(-1)]
    return "  " + "  ".join(vals)


def _stem(label: str | None, species: str) -> str:
    text = (label or "").strip()
    if text:
        return text.split("_")[0]
    return species
