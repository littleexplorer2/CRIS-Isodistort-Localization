"""Write ISODISTORT-style ``data.isoviz`` (!tag layout for IsoVIZ).

Must open in the ISOTROPY Suite **IsoVIZ** program (CRIS root ``ISOViz``
shortcut). Layout follows the official Save interactive distortion ascii
closely enough for IsoVIZ / ISOVIZ_INPUT; byte-identical copy of the website
file is not required (see repo ``agent.md``).

IsoVIZ requires that each ``!displacivemodelist`` entry's mode-vector count
equals the number of ``!atomcoordlist`` rows for that ``parentatom`` type.
Types must be parent Wyckoff stems (Eu1/Al1/Al2), not bare elements.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
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

# Parent-site stem in mode labels, e.g. [Al2:e:dsp]
_MODE_SITE_RE = re.compile(r"\[([A-Za-z]+\d*)\s*:")


@dataclass
class _AtomLayout:
    """Atom-type layout shared by ``!atom*`` and ``!displacivemodelist``."""

    types: list[tuple[str, str]]  # (stem, element), 1-based index = position+1
    # For each child site (in structure order): 1-based parent-type index
    child_type_idx: list[int]


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
    lines.extend(_basis_block(spec, parent_short, child_short))
    lines.extend(_space_group_block("child", child_n, child_short, child))
    lines.extend(_ops_block("child", child_n, child_short))
    lines.extend(_ops_block("parent", parent_sg, parent_short))
    atom_lines, layout = _atom_blocks(spec, parent, child)
    lines.extend(atom_lines)
    lines.extend(_mode_blocks(spec, layout))
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


def _basis_block(spec: Any, parent_short: str, child_short: str) -> list[str]:
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
    parent_c2p = _conv2prim_matrix(parent_short)
    child_c2p = _conv2prim_matrix(child_short)
    for tag, role, mat in (
        ("conv2primparentbasis", "parent", parent_c2p),
        ("conv2primchildbasis", "child", child_c2p),
    ):
        lines.append(
            f"#rows_are_primitive_{role}_basis_vectors_relativeto_conventional_{role}_basis "
        )
        lines.append(f"!{tag} ")
        for row in mat:
            lines.append("  " + "  ".join(f"{float(v):8.5f}" for v in row))
        lines.append("")
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


def _atom_blocks(
    spec: Any, parent: Structure | None, child: Structure
) -> tuple[list[str], _AtomLayout]:
    types = _parent_type_list(spec, parent, child)
    stem_to_idx = {stem: i for i, (stem, _el) in enumerate(types, start=1)}
    origin = np.asarray(getattr(spec.subgroup, "origin", None) or [0.0, 0.0, 0.0], dtype=float)

    child_type_idx: list[int] = []
    for site in child:
        stem = _child_site_stem(site, spec, origin)
        tidx = stem_to_idx.get(stem)
        if tidx is None:
            # Fall back: first type with matching element, else type 1
            el = site.species_string
            tidx = next((i for i, (_s, e) in enumerate(types, start=1) if e == el), 1)
        child_type_idx.append(tidx)

    layout = _AtomLayout(types=types, child_type_idx=child_type_idx)
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
    counters: dict[int, int] = {}
    for site, tidx in zip(child, child_type_idx, strict=True):
        counters[tidx] = counters.get(tidx, 0) + 1
        x, y, z = (float(c) for c in site.frac_coords)
        lines.append(
            f"    {tidx:d}    1    {counters[tidx]}  {x:8.5f}  {y:8.5f}  {z:8.5f} "
        )
    lines.append("")
    n_atoms = len(child)
    lines.append("!atomsinunitcell ")
    for _ in range(n_atoms):
        lines.append("  1")
    lines.append("")
    lines.append("!atomocclist ")
    for _ in range(n_atoms):
        lines.append("  1.00000")
    lines.append("")
    for empty in ("!atommaglist ", "!atomrotlist ", "!bondlist "):
        lines.append(empty)
        lines.append("")
    return lines, layout


def _mode_blocks(spec: Any, layout: _AtomLayout) -> list[str]:
    from .distortion_formats import (
        _centering_multiplicity,
        cart_normalized_mode_matrix,
        compact_mode_label,
    )

    sg = spec.subgroup
    irrep = (sg.irrep_label or "IR").strip() or "IR"
    child = spec.structure
    n_c = _centering_multiplicity(sg.space_group_symbol or "")
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
    stem_to_idx = {stem: i for i, (stem, _el) in enumerate(layout.types, start=1)}
    n_child = len(layout.child_type_idx)
    for n, (key, arr) in enumerate(disp.items(), start=1):
        pretty = labels.get(key, key)
        amp = float(amps.get(key, 0.0))
        mat = np.asarray(arr, dtype=float)
        if mat.ndim != 2 or mat.shape[1] != 3:
            mat = np.zeros((n_child, 3), dtype=float)
        if mat.shape[0] != n_child:
            padded = np.zeros((n_child, 3), dtype=float)
            rows = min(n_child, mat.shape[0])
            padded[:rows] = mat[:rows]
            mat = padded
        parentatom = _mode_parentatom(pretty, mat, layout, stem_to_idx)
        type_rows = [i for i, t in enumerate(layout.child_type_idx) if t == parentatom]
        if not type_rows:
            type_rows = list(range(n_child))
            parentatom = 1
        scaled, maxamp_hint = cart_normalized_mode_matrix(
            mat, child.lattice.matrix, centering_mult=n_c
        )
        vecs = scaled[type_rows]
        maxamp = float(maxamp_hint) if maxamp_hint > 1e-12 else 1.0
        short_label = compact_mode_label(_isoviz_mode_label(pretty, irrep))
        lines.append(
            f"    {parentatom}    {n}  {amp:8.5f}  {maxamp:8.5f}    1 {short_label} "
        )
        for vec in vecs:
            lines.append("  " + "  ".join(f"{float(v):8.5f}" for v in vec))
    lines.append("")
    return lines


def _parent_type_list(
    spec: Any, parent: Structure | None, child: Structure
) -> list[tuple[str, str]]:
    """Ordered unique parent Wyckoff stems (Eu1, Al1, Al2)."""
    from .isodistort_cif import _parent_site_names

    types: list[tuple[str, str]] = []
    seen: set[str] = set()
    if parent is not None:
        names = _parent_site_names(spec)
        for i, site in enumerate(parent):
            stem = names.get(i) or _stem(site.label, site.species_string)
            if stem in seen:
                continue
            seen.add(stem)
            types.append((stem, site.species_string))
    if types:
        return types
    for site in child:
        stem = _stem(site.label, site.species_string)
        if stem in seen:
            continue
        seen.add(stem)
        types.append((stem, site.species_string))
    return types or [("X1", "X")]


def _child_site_stem(site: Any, spec: Any, origin: np.ndarray) -> str:
    from .isodistort_cif import _site_stem

    try:
        return _site_stem(site, np.asarray(site.frac_coords, dtype=float), spec, origin)
    except Exception:  # noqa: BLE001
        return _stem(getattr(site, "label", None), site.species_string)


def _mode_parentatom(
    pretty: str,
    mat: np.ndarray,
    layout: _AtomLayout,
    stem_to_idx: dict[str, int],
) -> int:
    m = _MODE_SITE_RE.search(pretty or "")
    if m:
        stem = m.group(1)
        if stem in stem_to_idx:
            return stem_to_idx[stem]
        # Element-only label (e.g. [Eu:a:dsp]) → first matching type stem
        for s, idx in stem_to_idx.items():
            if s == stem or s.startswith(stem):
                return idx
    norms = np.linalg.norm(mat, axis=1)
    active = [i for i, n in enumerate(norms) if n > 1e-8]
    if not active:
        return 1
    votes: dict[int, float] = {}
    for i in active:
        t = layout.child_type_idx[i] if i < len(layout.child_type_idx) else 1
        votes[t] = votes.get(t, 0.0) + float(norms[i])
    return max(votes, key=votes.get)  # type: ignore[arg-type]


def _isoviz_mode_label(pretty: str, irrep: str) -> str:
    """Prefer compact ``GM1+[Al2:e:dsp]A1(a)``-style labels when present."""
    text = (pretty or "").strip()
    if not text:
        return irrep
    # I4/mmm[0,0,0]GM1+(a)[Al2:e:dsp]A1(a) → GM1+[Al2:e:dsp]A1(a)
    m = re.search(
        r"((?:GM|LD|DT|X|M|A|R|Z|N|P|H|V)\d*[+-]?)"
        r"(?:\([^)]*\))?"
        r"(\[[^\]]+\].*)",
        text,
    )
    if m:
        return f"{m.group(1)}{m.group(2)}"
    m2 = re.search(rf"({re.escape(irrep)}\S*)", text)
    if m2:
        return m2.group(1)
    return text


def _conv2prim_matrix(symbol: str) -> np.ndarray:
    """Conventional → primitive basis rows (ISODISTORT / IsoVIZ convention)."""
    letter = (symbol or "P").lstrip("0123456789 ").strip()[:1].upper()
    if letter == "I":
        return np.array(
            [
                [-0.5, 0.5, 0.5],
                [0.5, -0.5, 0.5],
                [0.5, 0.5, -0.5],
            ],
            dtype=float,
        )
    if letter == "F":
        return np.array(
            [
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ],
            dtype=float,
        )
    if letter == "C":
        return np.array(
            [
                [0.5, 0.5, 0.0],
                [-0.5, 0.5, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=float,
        )
    if letter == "A":
        return np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 0.5, -0.5],
                [0.0, 0.5, 0.5],
            ],
            dtype=float,
        )
    if letter == "B":
        return np.array(
            [
                [0.5, 0.0, -0.5],
                [0.0, 1.0, 0.0],
                [0.5, 0.0, 0.5],
            ],
            dtype=float,
        )
    if letter == "R":
        # Hexagonal setting of rhombohedral — identity is safest for IsoVIZ
        return np.eye(3)
    return np.eye(3)


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
