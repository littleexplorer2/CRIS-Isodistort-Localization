"""Write ISODISTORT-style CIF (header, subgroup setting, iso_* loops).

Must open in **VESTA** (CRIS root shortcut). Official Distortion-page CIF
(ISODISTORT 6.12.2) is the format target: isotropy-subgroup conventional cell,
ITA origin choice 2, asymmetric-unit sites, Hall setting, and ISODISTORT
private loops. Content/layout should track the website as closely as practical;
byte-identical export is not a project DoD (see repo ``agent.md``).
"""
from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, Any

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SpaceGroup

from ..backend import SubgroupInfo
from ..utils.parent_header import format_fixed_coord, format_wyckoff_sites
from ..utils.schoenflies import hm_symbol, schoenflies_symbol

if TYPE_CHECKING:
    from .distortion_formats import SubgroupExportSpec

_ISO_AUTHORS = (
    "# Harold T. Stokes, Branton J. Campbell, David Tanner, Dorian M. Hatch",
    "# Brigham Young University, Provo, Utah, USA",
)

_PREFERENCES = (
    "monoclinic axes a(b)c, monoclinic cell choice 1, orthorhombic axes abc, "
    "origin choice 2, hexagonal axes, SSG standard setting"
)

# ISODISTORT 6.12 still uses 1992 HM extended names (Cmma, not Cmme).
# Keys are ITA numbers; values are the extended H-M string used in Distortion CIF
# (including origin-choice-2 suffixes). This is space-group metadata, not a
# per-irrep / per-calculation answer table.
_HM_FULL_1992 = {
    2: "P -1",
    4: "P 1 21 1",
    5: "C 1 2 1",
    6: "P 1 m 1",
    7: "P 1 c 1",
    8: "C 1 m 1",
    10: "P 1 2/m 1",
    11: "P 1 21/m 1",
    12: "C 1 2/m 1",
    13: "P 1 2/c 1",
    14: "P 1 21/c 1",
    15: "C 1 2/c 1",
    17: "P 2 2 21",
    19: "P 21 21 21",
    20: "C 2 2 21",
    22: "F 2 2 2",
    23: "I 2 2 2",
    24: "I 21 21 21",
    25: "P m m 2",
    26: "P m c 21",
    28: "P m a 2",
    30: "P n c 2",
    31: "P m n 21",
    36: "C m c 21",
    38: "A m m 2",
    39: "A b m 2",
    40: "A b a 2",
    42: "F m m 2",
    43: "F d d 2",
    44: "I m m 2",
    46: "I m a 2",
    47: "P 2/m 2/m 2/m",
    48: "P 2/n 2/n 2/n (origin choice 2)",
    49: "P 2/c 2/c 2/m",
    50: "P 2/b 2/a 2/n (origin choice 2)",
    51: "P 21/m 2/m 2/a",
    53: "P 2/m 2/n 21/a",
    55: "P 21/b 21/a 2/m",
    56: "P 21/c 21/c 2/n",
    58: "P 21/n 21/n 2/m",
    59: "P 21/m 21/m 2/n (origin choice 2)",
    62: "P 21/n 21/m 21/a",
    63: "C 2/m 2/c 21/m",
    64: "C 2/m 2/c 21/a",
    65: "C 2/m 2/m 2/m",
    66: "C 2/c 2/c 2/m",
    67: "C 2/m 2/m 2/a",
    68: "C 2/c 2/c 2/a (origin choice 2)",
    69: "F 2/m 2/m 2/m",
    70: "F 2/d 2/d 2/d (origin choice 2)",
    71: "I 2/m 2/m 2/m",
    72: "I 2/b 2/a 2/m",
    73: "I 21/b 21/c 21/a",
    74: "I 21/m 21/m 21/a",
    76: "P 41",
    82: "I -4",
    91: "P 41 2 2",
    92: "P 41 21 2",
    99: "P 4 m m",
    105: "P 42 m c",
    107: "I 4 m m",
    115: "P -4 m 2",
    119: "I -4 m 2",
    120: "I -4 c 2",
    121: "I -4 2 m",
    122: "I -4 2 d",
    123: "P 4/m 2/m 2/m",
    125: "P 4/n 2/b 2/m (origin choice 2)",
    127: "P 4/m 21/b 2/m",
    129: "P 4/n 21/m 2/m (origin choice 2)",
    131: "P 42/m 2/m 2/c",
    132: "P 42/m 2/c 2/m",
    134: "P 42/n 2/n 2/m (origin choice 2)",
    136: "P 42/m 21/n 2/m",
    137: "P 42/n 21/m 2/c (origin choice 2)",
    138: "P 42/n 21/c 2/m (origin choice 2)",
    139: "I 4/m 2/m 2/m",
    140: "I 4/m 2/c 2/m",
    141: "I 41/a 2/m 2/d (origin choice 2)",
    142: "I 41/a 2/c 2/d (origin choice 2)",
}

# ITA origin-choice-2 representatives for Cmma (#67). 8n is (x, 1/4, z).
_WYCKOFF_67 = (
    ("a", 4, np.array([0.25, 0.0, 0.0])),
    ("b", 4, np.array([0.25, 0.0, 0.5])),
    ("c", 4, np.array([0.0, 0.25, 0.0])),
    ("d", 4, np.array([0.0, 0.0, 0.25])),
    ("e", 4, np.array([0.25, 0.25, 0.0])),
    ("f", 4, np.array([0.0, 0.25, 0.5])),
    ("g", 4, np.array([0.0, 0.25, 0.25])),
    ("h", 4, np.array([0.25, 0.25, 0.5])),
)

# I4/mmm (#139) special positions used by EuAl4-type parents.
_WYCKOFF_139 = (
    ("a", 2, np.array([0.0, 0.0, 0.0])),
    ("b", 2, np.array([0.0, 0.0, 0.5])),
    ("c", 4, np.array([0.0, 0.5, 0.0])),
    ("d", 4, np.array([0.0, 0.5, 0.25])),
    ("e", 4, np.array([0.0, 0.0, 0.25])),  # (0,0,z); z free — seed at 1/4
)

# P4mm (#99): 1a (0,0,z), 1b (1/2,1/2,z), 2c (1/2,0,z).
_WYCKOFF_99 = (
    ("a", 1, np.array([0.0, 0.0, 0.0])),
    ("b", 1, np.array([0.5, 0.5, 0.0])),
    ("c", 2, np.array([0.5, 0.0, 0.0])),
)

# Immm (#71): 2a; 4i (0,0,z); 4j (1/2,0,z).
_WYCKOFF_71 = (
    ("a", 2, np.array([0.0, 0.0, 0.0])),
    ("i", 4, np.array([0.0, 0.0, 0.25])),
    ("j", 4, np.array([0.5, 0.0, 0.25])),
)

# Fmmm (#69): 4a; 8f; 8i (0,0,z).
_WYCKOFF_69 = (
    ("a", 4, np.array([0.0, 0.0, 0.0])),
    ("f", 8, np.array([0.25, 0.25, 0.25])),
    ("i", 8, np.array([0.0, 0.0, 0.25])),
)

# I4mm (#107): 2a (0,0,z); 4b (0,1/2,z).
_WYCKOFF_107 = (
    ("a", 2, np.array([0.0, 0.0, 0.0])),
    ("b", 4, np.array([0.0, 0.5, 0.0])),
)

# P4/mmm (#123).
_WYCKOFF_123 = (
    ("a", 1, np.array([0.0, 0.0, 0.0])),
    ("d", 1, np.array([0.5, 0.5, 0.5])),
    ("g", 2, np.array([0.0, 0.0, 0.25])),
    ("h", 2, np.array([0.5, 0.5, 0.25])),
    ("i", 4, np.array([0.0, 0.5, 0.25])),
)

# P4/nmm (#129) origin choice 2.
_WYCKOFF_129 = (
    ("a", 2, np.array([0.75, 0.25, 0.0])),
    ("b", 2, np.array([0.75, 0.25, 0.5])),
    ("c", 2, np.array([0.25, 0.25, 0.25])),
)



@dataclass
class _Setting:
    int_number: int
    symbol: str
    full_symbol: str
    symmetry_ops: list
    hall: str


def render_isodistort_cif(
    structure: Structure,
    spec: SubgroupExportSpec | None = None,
) -> str:
    """Return CIF text in the official ISODISTORT layout."""
    if spec is None:
        dummy_sg = SubgroupInfo(index=0, space_group_number=1, space_group_symbol="P1")
        return _render_p1(structure, dummy_sg)
    return _render(spec, spec.cif_structure or spec.structure, force_p1=False)


def _render_p1(structure: Structure, sg: SubgroupInfo) -> str:
    from .distortion_formats import SubgroupExportSpec

    dummy = SubgroupExportSpec(subgroup=sg, structure=structure)
    return _render(dummy, structure, force_p1=True)


def _render(spec: Any, structure: Structure, *, force_p1: bool) -> str:
    sg = spec.subgroup
    setting = _setting(sg, force_p1=force_p1)
    shifted, origin_shift = _apply_origin_choice(structure, setting)
    sites = _asymmetric_sites(shifted, setting, spec, origin_shift)
    lines: list[str] = []
    lines.extend(_header_comments(spec, force_p1=force_p1))
    lines.append("")
    lines.append("data_isodistort-output")
    lines.append("")
    lat = shifted.lattice
    lines.append(f"_cell_length_a     {_f5(lat.a)}")
    lines.append(f"_cell_length_b     {_f5(lat.b)}")
    lines.append(f"_cell_length_c     {_f5(lat.c)}")
    lines.append(f"_cell_angle_alpha  {_f5(lat.alpha)}")
    lines.append(f"_cell_angle_beta   {_f5(lat.beta)}")
    lines.append(f"_cell_angle_gamma  {_f5(lat.gamma)}")
    lines.append(f"_cell_volume       {_f5(lat.volume)}")
    lines.append("")
    hm = _hm_full(setting, sg)
    lines.append(f'_symmetry_space_group_name_H-M "{hm}"')
    lines.append(f"_symmetry_Int_Tables_number {setting.int_number}")
    lines.append(f"_space_group.reference_setting '{setting.hall}'")
    lines.append("_space_group.transform_Pp_abc a,b,c;0,0,0")
    lines.append("")
    lines.extend(_symop_loop(setting.symmetry_ops))
    lines.append("")
    elements = list(dict.fromkeys(str(site["element"]) for site in sites))
    if not elements:
        elements = list(dict.fromkeys(s.species_string for s in shifted))
    lines.append("loop_")
    lines.append("_atom_type_symbol")
    for el in elements:
        lines.append(f"{el:<4s}")
    lines.append("")
    lines.extend(_atom_site_loop(sites))
    lines.append("")
    lines.extend(_iso_mode_loops(spec, sites))
    lines.extend(_parent_cell_block(spec, origin_shift))
    lines.append("")
    lines.append("# end of structure file")
    lines.append("# end of cif")
    lines.append("")
    return "\n".join(lines)


def _setting(sg: SubgroupInfo, *, force_p1: bool) -> _Setting:
    if force_p1 or sg.space_group_number in (0, 1):
        op = SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3))
        return _Setting(1, "P1", "P1", [op], "001:P 1")
    number = int(sg.space_group_number)
    symbol = hm_symbol(number) or (sg.space_group_symbol or "").strip() or str(number)
    hall_n, hall = _hall_for_number(number)
    ops = _ops_from_hall(hall_n)
    if not ops:
        pmg = _pymatgen_space_group(sg)
        ops = _normalized_ops(list(pmg.symmetry_ops))
    return _Setting(
        int_number=number,
        symbol=symbol,
        full_symbol=_HM_FULL_1992.get(number, ""),
        symmetry_ops=ops,
        hall=hall,
    )


def _pymatgen_space_group(sg: SubgroupInfo) -> SpaceGroup:
    symbol = (sg.space_group_symbol or "").strip()
    if symbol:
        try:
            return SpaceGroup(symbol.split()[0])
        except (ValueError, KeyError, TypeError):
            pass
    return SpaceGroup.from_int_number(int(sg.space_group_number))


def _hall_catalog() -> tuple[dict[int, str], dict[int, int]]:
    """Map IT number → Hall setting. Prefer origin choice 2 (ISODISTORT default)."""
    cached = getattr(_hall_catalog, "_cache", None)
    if cached is not None:
        return cached
    symbols: dict[int, str] = {}
    numbers: dict[int, int] = {}
    choices: dict[int, str] = {}
    try:
        import warnings

        from spglib import get_spacegroup_type
    except ImportError:
        empty: tuple[dict[int, str], dict[int, int]] = ({}, {})
        _hall_catalog._cache = empty
        return empty
    for hall_n in range(1, 531):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            info = get_spacegroup_type(hall_n)
        sg_n = int(info.number)
        choice = str(getattr(info, "choice", "") or "")
        hall_sym = f"{sg_n:03d}:{info.hall_symbol}"
        if sg_n not in symbols:
            symbols[sg_n] = hall_sym
            numbers[sg_n] = hall_n
            choices[sg_n] = choice
            continue
        # Upgrade to origin choice 2 when available (matches Distortion-page CIF).
        if choice == "2" and choices.get(sg_n) != "2":
            symbols[sg_n] = hall_sym
            numbers[sg_n] = hall_n
            choices[sg_n] = choice
    _hall_catalog._cache = (symbols, numbers)
    return symbols, numbers


def _hall_for_number(sg_number: int) -> tuple[int, str]:
    symbols, numbers = _hall_catalog()
    return numbers.get(sg_number, 0), symbols.get(sg_number, f"{sg_number:03d}")


def _ops_from_hall(hall_number: int) -> list[SymmOp]:
    if hall_number <= 0:
        return []
    try:
        from spglib import get_symmetry_from_database
    except ImportError:
        return []
    data = get_symmetry_from_database(int(hall_number))
    ops: list[SymmOp] = []
    seen: set[str] = set()
    for rot, trans in zip(data["rotations"], data["translations"], strict=False):
        op = _wrapped_op(np.asarray(rot, dtype=float), np.asarray(trans, dtype=float))
        key = _xyz(op)
        if key not in seen:
            seen.add(key)
            ops.append(op)
    return _order_ops(ops)


def _wrapped_op(rot: np.ndarray, trans: np.ndarray) -> SymmOp:
    t = np.mod(np.asarray(trans, dtype=float), 1.0)
    t = np.where(t > 1.0 - 1e-8, 0.0, t)
    return SymmOp.from_rotation_and_translation(rot, t)


def _normalized_ops(ops: Sequence) -> list[SymmOp]:
    out: list[SymmOp] = []
    seen: set[str] = set()
    for op in ops:
        wrapped = _wrapped_op(op.rotation_matrix, op.translation_vector)
        key = _xyz(wrapped)
        if key not in seen:
            seen.add(key)
            out.append(wrapped)
    return _order_ops(out)


def _order_ops(ops: Sequence[SymmOp]) -> list[SymmOp]:
    """Identity, then primitive-cell ops, then C-centering copies (ISODISTORT order)."""
    ident: list[SymmOp] = []
    primitive: list[SymmOp] = []
    centered: list[SymmOp] = []
    for op in ops:
        expr = _xyz(op)
        if expr == "x,y,z":
            ident.append(op)
            continue
        y_part = expr.split(",")[1] if "," in expr else ""
        if "y+1/2" in y_part or "y-1/2" in y_part:
            centered.append(op)
        else:
            primitive.append(op)
    return ident + primitive + centered


def _hm_full(setting: _Setting, sg: SubgroupInfo) -> str:
    if setting.int_number in _HM_FULL_1992:
        return _HM_FULL_1992[setting.int_number]
    if setting.full_symbol:
        return setting.full_symbol
    try:
        pmg = SpaceGroup.from_int_number(int(setting.int_number))
        spaced = _space_hm_symbol(pmg.full_symbol or pmg.symbol)
        if spaced:
            return spaced
    except (ValueError, KeyError, TypeError):
        pass
    raw = (sg.space_group_symbol or setting.symbol or "").strip()
    return _space_hm_symbol(raw) or setting.symbol


def _space_hm_symbol(raw: str) -> str:
    """Format compact / full HM tokens like the Distortion-page CIF.

    `P-1` → `P -1`; `I2/m2/m2/m` → `I 2/m 2/m 2/m`;
    `Imm2` → `I m m 2`; `P4_2mc` → `P 42 m c`.
    """
    text = re.sub(r"\s+", "", raw or "")
    if not text:
        return ""
    # Centrosymmetric triclinic
    if text in {"P-1", "P1bar", "P\\bar{1}"}:
        return "P -1"
    letter = text[0]
    rest = text[1:]
    # Official CIF often drops underscore in screw axes: 4_2 → 42, 2_1 → 21
    rest = re.sub(r"(\d)_(\d)", r"\1\2", rest)
    # Axis tokens: numeric (optional /letter), barred axis, or single letter.
    # Do NOT glue trailing digits onto m/a/b/c (``mm2`` → m, m, 2 not m2).
    # Limit slash glue to one letter so ``4/mmm`` → 4/m, m, m.
    parts = re.findall(r"-?\d+(?:/[A-Za-z])?|[A-Za-z]", rest)
    if not parts:
        parts_list: list[str] = [letter]
        i = 0
        while i < len(rest):
            ch = rest[i]
            if ch == "-":
                j = i + 1
                while j < len(rest) and rest[j].isdigit():
                    j += 1
                parts_list.append(rest[i:j])
                i = j
                continue
            if ch == "/":
                if len(parts_list) > 1:
                    parts_list[-1] = parts_list[-1] + "/"
                i += 1
                continue
            if ch.isdigit():
                j = i + 1
                while j < len(rest) and (rest[j].isdigit() or rest[j] == "_"):
                    j += 1
                token = rest[i:j].replace("_", "")
                if len(parts_list) > 1 and parts_list[-1].endswith("/"):
                    parts_list[-1] = parts_list[-1] + token
                else:
                    parts_list.append(token)
                i = j
                continue
            if ch.isalpha():
                # Single letter only; digits after mirrors are separate axes.
                if len(parts_list) > 1 and parts_list[-1].endswith("/"):
                    parts_list[-1] = parts_list[-1] + ch
                else:
                    parts_list.append(ch)
                i += 1
                continue
            i += 1
        return " ".join(parts_list)
    return " ".join([letter, *parts])


def _apply_origin_choice(structure: Structure, setting: _Setting) -> tuple[Structure, np.ndarray]:
    """Shift onto the conventional ITA setting (origin choice 2 when required).

    Prefer 0 / 1/2 specials and the unshifted cell. Quarter shifts are used only
    when they improve the asymmetric-unit count or match known ITA tables
    (e.g. Cmma origin-2). Prefer the zero shift whenever ASU cardinality ties.
    """
    if setting.int_number <= 1 or len(structure) == 0:
        return structure, np.zeros(3)
    ops = list(setting.symmetry_ops)
    species = [s.species_string for s in structure]
    best_shift = np.zeros(3)
    best_score: tuple | None = None
    for dx in range(4):
        for dy in range(4):
            for dz in range(4):
                shift = np.array([dx, dy, dz], dtype=float) / 4.0
                coords = np.mod(structure.frac_coords + shift, 1.0)
                n_unique = _count_unique(coords, species, ops)
                ita = _ita_match_score(setting.int_number, coords, ops)
                half = float(np.sum(_half_special_score(coords)))
                quarters = float(np.sum(_quarter_score(coords)))
                shift_mag = float(np.linalg.norm(shift))
                # Minimize unique ASU sites; maximize ITA / half-specials first
                # (origin choice 2 often needs a non-zero shift); only then prefer
                # the identity origin and fewer quarter coords.
                score = (
                    n_unique,
                    -ita,
                    -half,
                    0 if shift_mag < 1e-8 else 1,
                    quarters,
                    shift_mag,
                )
                if best_score is None or score < best_score:
                    best_score = score
                    best_shift = shift
    new_coords = np.mod(structure.frac_coords + best_shift, 1.0)
    return _with_frac_coords(structure, new_coords), best_shift


def _with_frac_coords(structure: Structure, frac_coords: np.ndarray) -> Structure:
    from pymatgen.core import PeriodicSite

    sites = [
        PeriodicSite(
            site.species,
            coord,
            structure.lattice,
            coords_are_cartesian=False,
            label=site.label,
            properties=site.properties,
        )
        for site, coord in zip(structure, frac_coords, strict=True)
    ]
    return Structure.from_sites(sites)


def _count_unique(coords: np.ndarray, species: Sequence[str], ops, tol: float = 1e-4) -> int:
    used = [False] * len(coords)
    n_unique = 0
    for i, frac in enumerate(coords):
        if used[i]:
            continue
        n_unique += 1
        used[i] = True
        for op in ops:
            img = np.mod(op.operate(frac), 1.0)
            for j, other in enumerate(coords):
                if used[j] or species[j] != species[i]:
                    continue
                delta = np.abs(img - other)
                delta = np.minimum(delta, 1.0 - delta)
                if np.all(delta < tol):
                    used[j] = True
    return n_unique


def _special_score(coords: np.ndarray) -> np.ndarray:
    """Legacy helper: 0 / 1/4 / 1/2 / 3/4 specials (prefer half_special for ranking)."""
    targets = np.array([0.0, 0.25, 0.5, 0.75])
    dist = np.min(np.abs(coords[..., None] - targets), axis=-1)
    return (dist < 1e-4).astype(float)


def _half_special_score(coords: np.ndarray) -> np.ndarray:
    """Coords on 0 or 1/2 (standard specials for most settings)."""
    wrapped = np.mod(np.asarray(coords, dtype=float), 1.0)
    wrapped = np.where(wrapped > 1.0 - 1e-8, 0.0, wrapped)
    targets = np.array([0.0, 0.5])
    dist = np.min(np.abs(wrapped[..., None] - targets), axis=-1)
    return (dist < 1e-4).astype(float)


def _quarter_score(coords: np.ndarray) -> np.ndarray:
    dist = np.min(np.abs(coords[..., None] - np.array([0.25, 0.75])), axis=-1)
    return (dist < 1e-4).astype(float)


def _ita_match_score(sg_number: int, coords: np.ndarray, ops) -> float:
    table = _wyckoff_table(sg_number)
    if not table:
        return 0.0
    score = 0.0
    for frac in coords:
        for _letter, _mult, rep in table:
            if _in_orbit(frac, rep, ops):
                score += 1.0
                break
    return score


def _asymmetric_sites(
    structure: Structure,
    setting: _Setting,
    spec: Any,
    origin_shift: np.ndarray,
) -> list[dict]:
    ops = list(setting.symmetry_ops)
    unique_idx: list[int] = []
    for i, site in enumerate(structure):
        duplicate = False
        for j in unique_idx:
            if site.species_string != structure[j].species_string:
                continue
            if _in_orbit(site.frac_coords, structure[j].frac_coords, ops):
                duplicate = True
                break
        if not duplicate:
            unique_idx.append(i)
    stem_order = _stem_appearance_order(structure)
    spglib_letters = _spglib_wyckoff_letters(structure, setting.int_number)
    has_table = bool(_wyckoff_table(setting.int_number))
    sites: list[dict] = []
    for i in unique_idx:
        frac = _canonical_frac(structure[i].frac_coords, ops, keep_near=structure[i].frac_coords)
        mult = _orbit_size(frac, ops)
        if has_table:
            letter = _wyckoff_letter(setting.int_number, frac, mult, ops)
        elif spglib_letters and i < len(spglib_letters) and spglib_letters[i]:
            letter = spglib_letters[i]
        else:
            letter = _wyckoff_letter(setting.int_number, frac, mult, ops)
        frac = _snap_wyckoff_rep(setting.int_number, frac, letter, mult, ops)
        stem = _site_stem(structure[i], frac, spec, origin_shift)
        sites.append(
            {
                "stem": stem,
                "element": structure[i].species_string,
                "multiplicity": mult,
                "wyckoff": letter,
                "frac": frac,
                "occupancy": _occupancy(structure[i]),
                "symmform": _symmform(frac, ops),
            }
        )
    sites.sort(
        key=lambda row: (
            stem_order.get(row["stem"], 99),
            row["wyckoff"],
            tuple(np.round(row["frac"], 5)),
        )
    )
    counts: dict[str, int] = {}
    for row in sites:
        counts[row["stem"]] = counts.get(row["stem"], 0) + 1
        row["label"] = f"{row['stem']}_{counts[row['stem']]}"
    return sites


def _spglib_wyckoff_letters(structure: Structure, sg_number: int) -> list[str] | None:
    """Wyckoff letters from spglib for the current cell (best-effort)."""
    if sg_number <= 1 or len(structure) == 0:
        return None
    try:
        from spglib import get_symmetry_dataset
    except ImportError:
        return None
    nums = []
    for site in structure:
        try:
            nums.append(int(site.specie.Z))
        except (AttributeError, TypeError, ValueError):
            return None
    cell = (
        np.asarray(structure.lattice.matrix, dtype=float),
        np.asarray(structure.frac_coords, dtype=float),
        nums,
    )
    try:
        ds = get_symmetry_dataset(cell, symprec=1e-4)
    except Exception:  # noqa: BLE001 - optional enrichment
        return None
    if ds is None:
        return None
    number = int(getattr(ds, "number", 0) or 0)
    if number and number != int(sg_number):
        return None
    wyckoffs = getattr(ds, "wyckoffs", None)
    if not wyckoffs or len(wyckoffs) != len(structure):
        return None
    return [str(w) for w in wyckoffs]


def _stem_appearance_order(structure: Structure) -> dict[str, int]:
    order: dict[str, int] = {}
    for site in structure:
        stem = _label_stem(site.label, site.species_string)
        if stem not in order:
            order[stem] = len(order)
    return order


def _in_orbit(frac: np.ndarray, other: np.ndarray, ops, tol: float = 1e-4) -> bool:
    for op in ops:
        img = np.mod(op.operate(frac), 1.0)
        delta = np.abs(img - other)
        delta = np.minimum(delta, 1.0 - delta)
        if np.all(delta < tol):
            return True
    return False


def _orbit_size(frac: np.ndarray, ops, tol: float = 1e-4) -> int:
    uniq: list[np.ndarray] = []
    for op in ops:
        img = _wrap(op.operate(frac))
        if not any(np.all(np.minimum(np.abs(img - q), 1.0 - np.abs(img - q)) < tol) for q in uniq):
            uniq.append(img)
    return len(uniq)


def _canonical_frac(frac: np.ndarray, ops, keep_near: np.ndarray | None = None) -> np.ndarray:
    orbit: list[np.ndarray] = []
    for op in ops:
        img = _wrap(op.operate(frac))
        if not any(
            np.all(np.minimum(np.abs(img - q), 1.0 - np.abs(img - q)) < 1e-4) for q in orbit
        ):
            orbit.append(img)
    near = keep_near if keep_near is not None else frac

    def key(point: np.ndarray) -> tuple:
        half = -float(np.sum(_half_special_score(point.reshape(1, 3))))
        quarters = float(np.sum(_quarter_score(point.reshape(1, 3))))
        dist = float(np.linalg.norm(np.minimum(np.abs(point - near), 1.0 - np.abs(point - near))))
        return (half, quarters, dist, tuple(np.round(point, 8)))

    return np.array(sorted(orbit, key=key)[0], dtype=float)


def _wrap(frac: np.ndarray, tol: float = 1e-8) -> np.ndarray:
    out = np.mod(np.asarray(frac, dtype=float), 1.0)
    out = np.where(out > 1.0 - tol, 0.0, out)
    return out


def _wyckoff_table(sg_number: int) -> tuple:
    if sg_number == 67:
        return _WYCKOFF_67
    if sg_number == 69:
        return _WYCKOFF_69
    if sg_number == 71:
        return _WYCKOFF_71
    if sg_number == 99:
        return _WYCKOFF_99
    if sg_number == 107:
        return _WYCKOFF_107
    if sg_number == 123:
        return _WYCKOFF_123
    if sg_number == 129:
        return _WYCKOFF_129
    if sg_number == 139:
        return _WYCKOFF_139
    return ()


def _wyckoff_letter(sg_number: int, frac: np.ndarray, multiplicity: int, ops) -> str:
    table = _wyckoff_table(sg_number)
    if sg_number == 139 and multiplicity == 4:
        x, y, _z = (float(c) for c in frac)
        if abs(x) < 1e-4 and abs(y) < 1e-4:
            return "e"
        if (
            (abs(x) < 1e-4 and abs(y - 0.5) < 1e-4)
            or (abs(x - 0.5) < 1e-4 and abs(y) < 1e-4)
        ):
            return "d"
    if sg_number == 71 and multiplicity == 4:
        x, y, _z = (float(c) for c in frac)
        if abs(x) < 1e-4 and abs(y) < 1e-4:
            return "i"
        if (
            (abs(x - 0.5) < 1e-4 and abs(y) < 1e-4)
            or (abs(x) < 1e-4 and abs(y - 0.5) < 1e-4)
        ):
            return "j"
    if sg_number == 69 and multiplicity == 8:
        x, y, _z = (float(c) for c in frac)
        if abs(x) < 1e-4 and abs(y) < 1e-4:
            return "i"
    if sg_number == 99:
        x, y, _z = (float(c) for c in frac)
        if multiplicity == 1:
            if abs(x) < 1e-4 and abs(y) < 1e-4:
                return "a"
            if abs(x - 0.5) < 1e-4 and abs(y - 0.5) < 1e-4:
                return "b"
        if multiplicity == 2:
            return "c"
    if sg_number == 107:
        x, y, _z = (float(c) for c in frac)
        if multiplicity == 2 and abs(x) < 1e-4 and abs(y) < 1e-4:
            return "a"
        if multiplicity == 4 and (
            (abs(x) < 1e-4 and abs(y - 0.5) < 1e-4)
            or (abs(x - 0.5) < 1e-4 and abs(y) < 1e-4)
        ):
            return "b"
    if sg_number == 123:
        x, y, _z = (float(c) for c in frac)
        if multiplicity == 1:
            if abs(x) < 1e-4 and abs(y) < 1e-4 and abs(_z) < 1e-4:
                return "a"
            if abs(x - 0.5) < 1e-4 and abs(y - 0.5) < 1e-4 and abs(_z - 0.5) < 1e-4:
                return "d"
        if multiplicity == 2:
            if abs(x) < 1e-4 and abs(y) < 1e-4:
                return "g"
            if abs(x - 0.5) < 1e-4 and abs(y - 0.5) < 1e-4:
                return "h"
        if multiplicity == 4 and (
            (abs(x) < 1e-4 and abs(y - 0.5) < 1e-4)
            or (abs(x - 0.5) < 1e-4 and abs(y) < 1e-4)
        ):
            return "i"
    if sg_number == 129 and multiplicity == 2:
        x, y, _z = (float(c) for c in frac)
        # Origin-choice-2: 2c (1/4,1/4,z); 2a/2b at (3/4,1/4,0/1/2)
        if abs(x - 0.25) < 1e-4 and abs(y - 0.25) < 1e-4:
            return "c"
        if abs(x - 0.75) < 1e-4 and abs(y - 0.25) < 1e-4:
            return "a" if abs(_z) < 1e-4 or abs(_z - 1.0) < 1e-4 else "b"
    for letter, mult, rep in table:
        if multiplicity == mult and _site_matches_wyckoff(
            frac, letter, mult, rep, ops, sg_number
        ):
            return letter
    if sg_number == 67:
        if multiplicity == 8:
            return "n"
        if multiplicity == 16:
            return "o"
    letters = "abcdefghijklmnopqrstuvwxyz"
    idx = min(max(multiplicity.bit_length() + 2, 0), len(letters) - 1)
    return letters[idx]


def _site_matches_wyckoff(
    frac: np.ndarray,
    letter: str,
    multiplicity: int,
    rep: np.ndarray,
    ops,
    sg_number: int,
) -> bool:
    x, y, _z = (float(c) for c in frac)
    if multiplicity == 4 and sg_number == 139 and letter == "e":
        return abs(x) < 1e-4 and abs(y) < 1e-4
    if multiplicity == 4 and sg_number == 71 and letter == "i":
        return abs(x) < 1e-4 and abs(y) < 1e-4
    if multiplicity == 4 and sg_number == 71 and letter == "j":
        return (
            (abs(x - 0.5) < 1e-4 and abs(y) < 1e-4)
            or (abs(x) < 1e-4 and abs(y - 0.5) < 1e-4)
        )
    if multiplicity == 8 and sg_number == 69 and letter == "i":
        return abs(x) < 1e-4 and abs(y) < 1e-4
    if sg_number == 123 and letter in {"g", "h", "i"}:
        if letter == "g":
            return abs(x) < 1e-4 and abs(y) < 1e-4
        if letter == "h":
            return abs(x - 0.5) < 1e-4 and abs(y - 0.5) < 1e-4
        return (
            (abs(x) < 1e-4 and abs(y - 0.5) < 1e-4)
            or (abs(x - 0.5) < 1e-4 and abs(y) < 1e-4)
        )
    if sg_number == 129 and letter == "c":
        return abs(x - 0.25) < 1e-4 and abs(y - 0.25) < 1e-4
    return _in_orbit(frac, rep, ops)


def _snap_wyckoff_rep(
    sg_number: int, frac: np.ndarray, letter: str, multiplicity: int, ops
) -> np.ndarray:
    # Free-parameter positions (…z, x…) must keep the structure's coordinate;
    # snapping to a table seed (z=1/4) would rewrite physical site locations.
    free_letters = {
        67: {"n", "g"},
        69: {"i", "g", "h"},
        71: {"i", "j"},
        99: {"a", "b", "c"},
        107: {"a", "b"},
        123: {"g", "h", "i"},
        129: {"c"},
        139: {"e"},
    }
    if letter in free_letters.get(sg_number, set()):
        return frac
    table = _wyckoff_table(sg_number)
    for wy_letter, mult, rep in table:
        if wy_letter == letter and multiplicity == mult and _in_orbit(frac, rep, ops):
            return np.array(rep, dtype=float)
    if sg_number == 67 and letter == "n":
        best = np.array(frac, dtype=float)
        best_dist = 1e9
        for op in ops:
            img = _wrap(op.operate(frac))
            if min(abs(img[1] - 0.25), 1.0 - abs(img[1] - 0.25)) > 1e-4:
                continue
            dist = float(np.linalg.norm(
                np.minimum(np.abs(img - frac), 1.0 - np.abs(img - frac))
            ))
            if dist < best_dist:
                best_dist = dist
                best = img
        if abs(best[1] - 0.75) < 1e-4:
            best = np.array([best[0], 0.25, best[2]], dtype=float)
        return _wrap(best)
    return frac


def _label_stem(label: str | None, element: str) -> str:
    raw = re.sub(r"[^A-Za-z0-9_]", "", str(label or ""))
    if raw:
        stripped = re.sub(r"_\d+$", "", raw)
        if re.match(r"^[A-Za-z]+\d+$", stripped):
            return stripped
        if re.match(r"^[A-Za-z]+$", stripped):
            return stripped
    return element


def _site_stem(site, frac: np.ndarray, spec: Any, origin_shift: np.ndarray) -> str:
    stem = _label_stem(site.label, site.species_string)
    if re.match(r"^[A-Za-z]+\d+$", stem):
        return stem
    parent = spec.parent_structure
    sg = spec.subgroup
    element = site.species_string
    if parent is not None and sg.basis_vectors and len(sg.basis_vectors) == 3:
        basis = np.asarray(sg.basis_vectors, dtype=float)
        try:
            frac_parent = _wrap((frac - origin_shift) @ basis)
            parent_names = _parent_site_names(spec)
            best_i = None
            best_d = 1e9
            for i, psite in enumerate(parent):
                if psite.species_string != element:
                    continue
                delta = np.abs(frac_parent - psite.frac_coords)
                delta = np.minimum(delta, 1.0 - delta)
                dist = float(np.max(delta))
                if dist < best_d:
                    best_d = dist
                    best_i = i
            if best_i is not None and best_d < 0.2:
                return parent_names.get(best_i, element)
        except (ValueError, np.linalg.LinAlgError, TypeError):
            return stem
    return stem


def _parent_site_names(spec: Any) -> dict[int, str]:
    parent = spec.parent_structure
    wyckoff = spec.parent_wyckoff_sites
    names: dict[int, str] = {}
    if parent is None:
        return names
    if wyckoff:
        counts: dict[str, int] = {}
        for site in wyckoff:
            species = str(site["species"])
            counts[species] = counts.get(species, 0) + 1
            tag = f"{species}{counts[species]}"
            for idx in site.get("equivalent_indices", [site.get("representative_index")]):
                if idx is not None:
                    names[int(idx)] = tag
        return names
    counts = {}
    for i, site in enumerate(parent):
        el = site.species_string
        counts[el] = counts.get(el, 0) + 1
        names[i] = f"{el}{counts[el]}"
    return names


def _occupancy(site) -> float:
    try:
        return float(site.species.num_atoms)
    except (TypeError, ValueError, AttributeError):
        return 1.0


def _symmform(frac: np.ndarray, ops, tol: float = 1e-4) -> str:
    """Free displacement axes under the site stabilizer (Dx,Dy,Dz)."""
    names = ("Dx", "Dy", "Dz")
    stabilizer = [
        op
        for op in ops
        if np.all(
            np.minimum(
                np.abs(_wrap(op.operate(frac)) - frac),
                1.0 - np.abs(_wrap(op.operate(frac)) - frac),
            )
            < tol
        )
    ]
    free = [True, True, True]
    axes = np.eye(3)
    for op in stabilizer:
        rot = np.asarray(op.rotation_matrix, dtype=float)
        for axis in range(3):
            mapped = rot @ axes[axis]
            if abs(mapped[axis] - 1.0) > 0.5:
                free[axis] = False
    parts = [names[i] if free[i] else "0" for i in range(3)]
    return ",".join(parts)


def _header_comments(spec: Any, *, force_p1: bool) -> list[str]:
    sg = spec.subgroup
    parent = spec.parent_structure
    parent_sg = int(spec.parent_sg or 0)
    parent_sym = spec.parent_symbol or (hm_symbol(parent_sg) if parent_sg else "")
    sch = schoenflies_symbol(parent_sg) if parent_sg else ""
    lines = [
        "# This file was generated by ISODISTORT, version 6.12.2",
        *_ISO_AUTHORS,
        "#",
    ]
    if parent_sg or parent_sym:
        sg_line = f"# Space Group: {parent_sg} {parent_sym}".rstrip()
        if sch:
            sg_line = f"{sg_line}     {sch}"
        lines.append(sg_line)
    lines.append(f"# Default space-group preferences: {_PREFERENCES}")
    if parent is not None:
        plat = parent.lattice
        lines.append(
            "# Lattice parameters: "
            f"a= {_f5(plat.a)}, b= {_f5(plat.b)}, c= {_f5(plat.c)}, "
            f"alpha= {_f5(plat.alpha)}, beta= {_f5(plat.beta)}, gamma= {_f5(plat.gamma)}"
        )
        wyckoff_lines = _parent_wyckoff_comment_lines(spec)
        if wyckoff_lines:
            for item in wyckoff_lines:
                lines.append(f"# {item}")
        else:
            for site in parent:
                x, y, z = (float(c) for c in site.frac_coords)
                lines.append(
                    f"# {site.label or site.species_string} "
                    f"({_frac_or_float(x)},{_frac_or_float(y)},{_frac_or_float(z)})"
                )
    if not force_p1:
        include = _include_distortions_comment(spec)
        if include:
            lines.append(include)
        if sg.k_point_label:
            lines.append(f"# {_k_point_comment(sg)}")
        if sg.irrep_label:
            lines.append(f"# {_irrep_comment(sg)}")
        try:
            lines.append(f"# {sg.opd_line_body()}")
        except (TypeError, ValueError, AttributeError):
            from .distortion_formats import subgroup_label

            lines.append(
                f"# {subgroup_label(sg)}  {sg.space_group_number} {sg.space_group_symbol}"
            )
        lines.extend(_order_parameter_comment_lines(spec))
    if spec.note:
        lines.append(f"# {spec.note}")
    return lines


def _include_distortions_comment(spec: Any) -> str:
    types = getattr(spec, "distortion_types", None) or []
    names = [str(x).strip().lower() for x in types if str(x).strip()]
    if not names:
        return "# Include strain, displacive ALL distortions"
    parts: list[str] = []
    for key in ("strain", "displacive", "occupational", "magnetic", "rotational"):
        if key in names:
            parts.append(key if key != "displacive" else "displacive ALL")
    if not parts:
        parts = ["strain", "displacive ALL"]
    return "# Include " + ", ".join(parts) + " distortions"


def _k_point_comment(sg: SubgroupInfo) -> str:
    """Official style: ``k point: LD, k10 (0,0,g), g=1/6``."""
    label = sg.k_point_label or ""
    try:
        from ..data.kpoints_official import format_k_point_display

        return "k point: " + format_k_point_display(
            int(sg.parent_sg or 0),
            label,
            list(sg.k_parameters or []),
            list(sg.k_coordinates or []),
        )
    except Exception:  # noqa: BLE001 - header comment must not fail export
        coords = ",".join(str(c) for c in (sg.k_coordinates or []))
        extra = f" ({coords})" if coords else ""
        return f"k point: {label}{extra}"


def _irrep_comment(sg: SubgroupInfo) -> str:
    """Official style: ``IR: GM1+, k14t1`` (ML label + CDML/Kovalev tag)."""
    label = (sg.irrep_label or "").strip()
    if not label:
        return "IR: "
    tag = _irrep_kovalev_tag(sg)
    if tag:
        return f"IR: {label}, {tag}"
    return f"IR: {label}"


def _irrep_kovalev_tag(sg: SubgroupInfo) -> str | None:
    try:
        from ..data.irreps_cdml import lookup_irrep_kovalev
        from ..data.kpoints_official import KPOINT_OFFICIAL
    except Exception:  # noqa: BLE001
        return None
    parent = int(sg.parent_sg or 0)
    k_label = (sg.k_point_label or "").strip()
    ir = (sg.irrep_label or "").strip()
    k_kov = None
    entry = KPOINT_OFFICIAL.get(parent, {}).get(k_label)
    if entry:
        k_kov = entry[0]
    return lookup_irrep_kovalev(
        ir,
        k_kovalev=k_kov,
        parent_sg=parent or None,
        k_point_label=k_label or None,
    )


def _order_parameter_comment_lines(spec: Any) -> list[str]:
    """Official Distortion CIF order-parameter block.

    Head line::
        ``I4/mmm[0,0,0]GM1+(a) 139 I4/mmm s=1 i=1``
    Mode lines::
        ``[Al2:e:dsp]A1(a):  0.00000`` / ``strain_1(a):  0.00000``
    """
    labels = spec.mode_labels or {}
    amps = spec.amplitudes or {}
    disp = spec.mode_displacements_sc or {}
    sg = spec.subgroup
    if not labels and not disp and not getattr(spec, "strain_mode_labels", None):
        # Still emit the OP head when we know the subgroup path (zero modes).
        if not (sg and sg.irrep_label):
            return []
    lines = ["# Order parameter values:"]
    head = _order_parameter_head(spec)
    if head:
        lines.append(f"#  {head}")
    keys = list(disp.keys()) if disp else list(labels.keys())
    for key in keys:
        pretty = labels.get(key, key)
        amp = float(amps.get(key, 0.0))
        m = re.search(r"(\[[^\]]+:[^\]]+:dsp\][^\s:]*(?:\([^)]*\))?)", pretty)
        if not m:
            m = re.search(r"(\[[^\]]+\].*)$", pretty)
        token = m.group(1) if m else pretty
        lines.append(f"#     {token}:  {_f5(amp)}")
    for lab in getattr(spec, "strain_mode_labels", None) or []:
        amp = float((getattr(spec, "strain_amplitudes", None) or {}).get(lab, 0.0))
        lines.append(f"#     {lab}:  {_f5(amp)}")
    return lines


def _order_parameter_head(spec: Any) -> str:
    sg = spec.subgroup
    parent_sym = (spec.parent_symbol or "").strip() or "P1"
    parent_compact = parent_sym.replace(" ", "")
    k_coords = _k_coords_for_mode_label(sg)
    ir = (sg.irrep_label or "").strip()
    direction = _opd_direction_letter(sg)
    child_hm = (sg.space_group_symbol or "").strip() or str(sg.space_group_number)
    child_hm = child_hm.replace(" ", "")
    return (
        f"{parent_compact}[{k_coords}]{ir}({direction}) "
        f"{sg.space_group_number} {child_hm} "
        f"s={sg.size} i={sg.subgroup_index}"
    )


def _k_coords_for_mode_label(sg: SubgroupInfo) -> str:
    coords = list(sg.k_coordinates or [])
    params = list(sg.k_parameters or [])
    if params and coords:
        # Substitute numeric params into symbolic coords when possible.
        try:
            from ..data.kpoints_official import format_k_point_display

            disp = format_k_point_display(
                int(sg.parent_sg or 0),
                sg.k_point_label or "",
                params,
                coords,
            )
            m = re.search(r"\(([^)]*)\)", disp)
            if m:
                return m.group(1)
        except Exception:  # noqa: BLE001
            pass
    if coords:
        return ",".join(str(c) for c in coords)
    return "0,0,0"


def _opd_direction_letter(sg: SubgroupInfo) -> str:
    raw = (sg.opd_dir_raw or "").strip()
    if raw.startswith("(") and raw.endswith(")"):
        inner = raw[1:-1]
        # (a), (a,0), (a;a) → leading letter token used in official labels.
        tok = re.split(r"[,;]", inner)[0].strip()
        return tok or "a"
    return "a"


def _parent_wyckoff_comment_lines(spec: Any) -> list[str]:
    """Parent Wyckoff block in Distortion CIF comments.

    Prefer ``parent_wyckoff_lines`` (CIF atom_site order/labels, same as the
    search-page header). Fall back to formatting ``parent_wyckoff_sites``.
    """
    preset = getattr(spec, "parent_wyckoff_lines", None) or []
    if preset:
        return [str(x) for x in preset if str(x).strip()]
    parent = spec.parent_structure
    wyckoff = spec.parent_wyckoff_sites
    if parent is None or not wyckoff:
        return []
    try:
        return format_wyckoff_sites(parent, wyckoff)
    except (KeyError, TypeError, ValueError, IndexError):
        return []


def _symop_loop(ops: Sequence) -> list[str]:
    xyz = [_xyz(op) for op in ops]
    if "x,y,z" in xyz:
        xyz.remove("x,y,z")
        xyz.insert(0, "x,y,z")
    lines = [
        "loop_",
        "_space_group_symop_id",
        "_space_group_symop_operation_xyz",
    ]
    for i, expr in enumerate(xyz, start=1):
        lines.append(f"{i} {expr}")
    return lines


def _xyz(op) -> str:
    return op.as_xyz_str().replace(" ", "")


def _atom_site_loop(sites: list[dict]) -> list[str]:
    lines = [
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_symmetry_multiplicity",
        "_atom_site_Wyckoff_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
        "_atom_site_fract_symmform",
    ]
    if not sites:
        return lines
    for row in sites:
        x, y, z = (float(c) for c in row["frac"])
        label = f"{row['label']:<6s}"
        lines.append(
            f"{label} {row['element']:<4s}{row['multiplicity']:2d} "
            f"{row['wyckoff']}  {_f5(x)}  {_f5(y)}  {_f5(z)}  {_f5(row['occupancy'])} "
            f"{row['symmform']:<8s}"
        )
    return lines


def _iso_mode_loops(spec: Any, sites: list[dict]) -> list[str]:
    disp = spec.mode_displacements_sc or {}
    amps = spec.amplitudes or {}
    labels = spec.mode_labels or {}
    lines: list[str] = []
    n_disp = len(disp)
    lines.append(f"_iso_displacivemode_number    {n_disp}")
    lines.append("")
    if n_disp:
        lines.append("loop_")
        lines.append("_iso_displacivemode_ID")
        lines.append("_iso_displacivemode_label")
        lines.append("_iso_displacivemode_value")
        for i, (key, _arr) in enumerate(disp.items(), start=1):
            pretty = labels.get(key, key)
            amp = float(amps.get(key, 0.0))
            lines.append(f"   {i} {pretty}  {_f5(amp)}")
        lines.append("")
    lines.extend(_coordinate_loops(sites))
    lines.append("_iso_magneticmode_number    0")
    lines.append("")
    lines.append("")
    lines.append("_iso_rotationalmode_number    0")
    lines.append("")
    lines.append("")
    lines.append("_iso_occupancymode_number    0")
    lines.append("")
    lines.append("_iso_strainmode_number    0")
    lines.append("")
    return lines


def _coordinate_loops(sites: list[dict]) -> list[str]:
    if not sites:
        return []
    free_axes = []
    formula_lines = [
        "loop_",
        "_iso_coordinate_label",
        "_iso_coordinate_formula",
    ]
    names = ("x", "y", "z")
    disp_names = ("dx", "dy", "dz")
    for row in sites:
        parts = [p.strip() for p in str(row["symmform"]).split(",")]
        frac = [float(c) for c in row["frac"]]
        for i, axis in enumerate(names):
            token = _formula_coord(frac[i])
            label = f"{row['label']}_{axis}"
            if i < len(parts) and parts[i] != "0":
                extra = f"{row['label']}_{disp_names[i]}"
                formula = f'"{token} + {extra}"'
                free_axes.append(f"{row['label']}_{disp_names[i]}")
            else:
                formula = f'"{token}"'
            formula_lines.append(f"{label:<24s} {formula}")
    lines: list[str] = []
    if free_axes:
        lines.append("loop_")
        lines.append("_iso_deltacoordinate_ID")
        lines.append("_iso_deltacoordinate_label")
        lines.append("_iso_deltacoordinate_value")
        for i, name in enumerate(free_axes, start=1):
            lines.append(f"   {i} {name:<8s}  0.00000")
        lines.append("")
    lines.extend(formula_lines)
    lines.append("")
    return lines


def _formula_coord(value: float) -> str:
    special = format_fixed_coord(value)
    frac = Fraction(value).limit_denominator(16)
    if abs(float(frac) - value) < 1e-4 and frac.denominator in (1, 2, 3, 4, 6, 8):
        return special
    return f"{value:10.5f}".rstrip()


def _parent_cell_block(spec: Any, origin_shift: np.ndarray) -> list[str]:
    parent = spec.parent_structure
    lines: list[str] = []
    if parent is not None:
        plat = parent.lattice
        lines.append(f"_iso_parentcell_length_a     {_f5(plat.a)}")
        lines.append(f"_iso_parentcell_length_b     {_f5(plat.b)}")
        lines.append(f"_iso_parentcell_length_c     {_f5(plat.c)}")
        lines.append(f"_iso_parentcell_angle_alpha  {_f5(plat.alpha)}")
        lines.append(f"_iso_parentcell_angle_beta   {_f5(plat.beta)}")
        lines.append(f"_iso_parentcell_angle_gamma  {_f5(plat.gamma)}")
        lines.append(f"_iso_parentcell_volume       {_f5(plat.volume)}")
    transform = _parent_to_child_transform(spec.subgroup, origin_shift)
    if transform:
        lines.append(f"_iso_parent-to-child.transform_Pp_abc {transform}")
    return lines


def _parent_to_child_transform(sg: SubgroupInfo, origin_shift: np.ndarray) -> str:
    basis = sg.basis_vectors
    if not basis or len(basis) != 3:
        return ""
    axes = ",".join(_axis_to_abc(row) for row in basis)
    origin = _transform_origin(sg, origin_shift)
    return f"{axes};{origin}"


def _transform_origin(sg: SubgroupInfo, origin_shift: np.ndarray) -> str:
    """Parent→child origin for ``_iso_parent-to-child.transform_Pp_abc``.

    Official Distortion CIF writes the Method 1 OPD origin here. The extra
    ITA origin-choice shift used to place atoms on conventional specials is
    applied only to site coordinates, not folded into this token.
    """
    _ = origin_shift
    raw = (getattr(sg, "origin_raw", None) or "").strip()
    if raw.startswith("(") and raw.endswith(")"):
        return raw[1:-1]
    if raw:
        return raw
    p0 = np.asarray(sg.origin or [0.0, 0.0, 0.0], dtype=float)
    return ",".join(_frac_or_float(float(v)) for v in p0)


def _axis_to_abc(row: Sequence[float]) -> str:
    parts: list[str] = []
    for coef, axis in zip(row, "abc", strict=True):
        value = Fraction(float(coef)).limit_denominator(12)
        if value == 0:
            continue
        sign = "+" if value > 0 else "-"
        mag = abs(value)
        token = axis if mag == 1 else f"{mag}{axis}"
        if not parts:
            parts.append(token if sign == "+" else f"-{token}")
        else:
            parts.append(f"{sign}{token}")
    return "".join(parts) or "0"


def _frac_or_float(value: float) -> str:
    frac = Fraction(value).limit_denominator(16)
    if abs(float(frac) - value) > 1e-6:
        return _f5(value)
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def _f5(value: float) -> str:
    return f"{float(value):.5f}"
