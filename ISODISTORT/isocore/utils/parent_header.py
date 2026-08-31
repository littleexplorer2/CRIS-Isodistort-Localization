"""Format the parent-structure header to match the official ISODISTORT search page.

Official example (EuAl4, I4/mmm)::

    Eu1 2a (0,0,0),
    Al1 4d (0,1/2,1/4),
    Al2 4e (0,0,z), z= 0.38000

Display labels and site order come from the parent CIF ``_atom_site_*`` loop
when available (e.g. ``ND`` / ``NI`` as written), not from memorized structures.
Wyckoff multiplicity/letter still come from symmetry analysis of the loaded cell.
"""
from __future__ import annotations

import re
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

# Special fractional values treated as fixed Wyckoff coordinates (IT tables).
_SPECIAL_FRACTIONS: tuple[Fraction, ...] = (
    Fraction(0), Fraction(1, 8), Fraction(1, 6), Fraction(1, 4), Fraction(1, 3),
    Fraction(3, 8), Fraction(1, 2), Fraction(5, 8), Fraction(2, 3), Fraction(3, 4),
    Fraction(5, 6), Fraction(7, 8), Fraction(1),
)
_AXIS_LETTERS = ("x", "y", "z")


def _near_special(value: float, tol: float = 1e-4) -> Fraction | None:
    """Return the matching IT special fraction, or None if the coord is free."""
    x = float(value) % 1.0
    if x > 1.0 - tol:
        x = 0.0
    for special in _SPECIAL_FRACTIONS:
        if abs(x - float(special)) <= tol:
            return Fraction(0) if special == 1 else special
        if abs(x - float(special) + 1.0) <= tol:
            return Fraction(0) if special == 1 else special
    return None


def format_fixed_coord(value: float) -> str:
    """Render a fixed coordinate as ``0``, ``1/2``, ``1/4``, …"""
    special = _near_special(value)
    if special is None:
        f = Fraction(float(value)).limit_denominator(24)
    else:
        f = special
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"


def format_wyckoff_site(
    species: str,
    species_index: int,
    multiplicity: int,
    letter: str,
    frac_coords: Sequence[float],
    *,
    label: str | None = None,
) -> str:
    """One official site token, e.g. ``Al2 4e (0,0,z), z= 0.38000``."""
    coords = [float(x) for x in frac_coords]
    parts: list[str] = []
    free: list[tuple[str, float]] = []
    for axis, value in zip(_AXIS_LETTERS, coords, strict=True):
        if _near_special(value) is None:
            parts.append(axis)
            free.append((axis, value))
        else:
            parts.append(format_fixed_coord(value))
    name = (label or "").strip() or f"{species}{species_index}"
    body = f"{name} {int(multiplicity)}{letter} ({','.join(parts)})"
    if free:
        # Official prints ``z= 0.38000`` (five decimals, one space after =).
        extras = ", ".join(f"{ax}= {val:.5f}" for ax, val in free)
        body = f"{body}, {extras}"
    return body


def format_wyckoff_sites(
    structure: Structure,
    wyckoff_sites: Sequence[dict],
) -> list[str]:
    """Format parent Wyckoff sites (fallback when no CIF path is available)."""
    species_count: dict[str, int] = {}
    lines: list[str] = []
    for site in wyckoff_sites:
        species = str(site["species"])
        species_count[species] = species_count.get(species, 0) + 1
        idx = int(site["representative_index"])
        coords = structure[idx].frac_coords
        label = site.get("display_label")
        lines.append(format_wyckoff_site(
            species=species,
            species_index=species_count[species],
            multiplicity=int(site["multiplicity"]),
            letter=str(site["wyckoff_letter"]),
            frac_coords=coords,
            label=str(label) if label else None,
        ))
    return lines


def parse_cif_atom_site_rows(cif_path: str | Path) -> list[dict]:
    """Read asymmetric ``_atom_site_*`` rows from a CIF in file order.

    Preserves ``_atom_site_label`` / ``_atom_site_type_symbol`` text as written
    (e.g. ``ND`` / ``NI``), without pymatgen element normalization.
    """
    text = Path(cif_path).read_text(encoding="utf-8", errors="replace")
    rows: list[dict] = []
    for loop_m in re.finditer(
        r"(?is)loop_\s*((?:_[a-z0-9_.]+\s*)+)([^_]*?)(?=\nloop_|\n#|\ndata_|\Z)",
        text,
    ):
        tags = re.findall(r"_[a-z0-9_.]+", loop_m.group(1), flags=re.I)
        lower = [t.lower() for t in tags]
        if "_atom_site_fract_x" not in lower:
            continue
        body = loop_m.group(2)
        tokens: list[str] = []
        for line in body.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens.extend(line.split())
        n = len(tags)
        if n == 0 or len(tokens) < n:
            continue
        ix = {name: i for i, name in enumerate(lower)}
        for start in range(0, len(tokens) - n + 1, n):
            chunk = tokens[start : start + n]
            try:
                x = float(chunk[ix["_atom_site_fract_x"]])
                y = float(chunk[ix["_atom_site_fract_y"]])
                z = float(chunk[ix["_atom_site_fract_z"]])
            except (KeyError, ValueError, IndexError):
                continue
            label = chunk[ix["_atom_site_label"]] if "_atom_site_label" in ix else ""
            typ = (
                chunk[ix["_atom_site_type_symbol"]]
                if "_atom_site_type_symbol" in ix
                else label
            )
            rows.append({
                "label": label,
                "type_symbol": typ,
                "frac": (x, y, z),
            })
        if rows:
            break
    return rows


def _frac_dist(a: Sequence[float], b: Sequence[float]) -> float:
    delta = np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))
    delta = np.minimum(delta, 1.0 - delta)
    return float(np.max(delta))


def format_wyckoff_sites_from_cif(
    cif_path: str | Path,
    structure: Structure,
    wyckoff_sites: Sequence[dict],
    *,
    tol: float = 0.05,
) -> list[str] | None:
    """Build display lines from CIF atom-site order + symmetry Wyckoff letters.

    Returns ``None`` if the CIF has no usable atom_site loop (caller falls back).
    """
    asu = parse_cif_atom_site_rows(cif_path)
    if not asu or not wyckoff_sites:
        return None
    used: set[int] = set()
    lines: list[str] = []
    species_count: dict[str, int] = {}
    for atom in asu:
        best_i = None
        best_d = 1e9
        for i, site in enumerate(wyckoff_sites):
            if i in used:
                continue
            idx = int(site["representative_index"])
            d = _frac_dist(atom["frac"], structure[idx].frac_coords)
            if d < best_d:
                best_d = d
                best_i = i
        if best_i is None or best_d > tol:
            continue
        used.add(best_i)
        site = wyckoff_sites[best_i]
        species = str(site["species"])
        species_count[species] = species_count.get(species, 0) + 1
        display = (atom.get("label") or atom.get("type_symbol") or "").strip()
        lines.append(format_wyckoff_site(
            species=species,
            species_index=species_count[species],
            multiplicity=int(site["multiplicity"]),
            letter=str(site["wyckoff_letter"]),
            frac_coords=atom["frac"],
            label=display or None,
        ))
    if len(lines) != len(wyckoff_sites):
        return None
    return lines


def parent_wyckoff_display(
    structure: Structure,
    wyckoff_sites: Sequence[dict],
    cif_path: str | Path | None = None,
) -> list[str]:
    """Preferred entry: CIF-ordered labels when ``cif_path`` is set, else structure order."""
    if cif_path is not None and Path(cif_path).is_file():
        from_cif = format_wyckoff_sites_from_cif(cif_path, structure, wyckoff_sites)
        if from_cif:
            return from_cif
    return format_wyckoff_sites(structure, wyckoff_sites)


def format_parent_header(
    *,
    space_group_number: int,
    space_group_symbol: str,
    schoenflies: str,
    lattice: dict[str, float],
    preferences: str,
    wyckoff_lines: Sequence[str],
    html: bool = False,
) -> str:
    """Full parent block under ``Done.`` (Space Group / lattice / prefs / Wyckoff)."""
    br = "<br>" if html else "\n"
    sg = f"{space_group_number} {space_group_symbol}"
    if schoenflies:
        sg = f"{sg} {schoenflies}"
    lat = (
        f"a= {lattice['a']:.5f}, b= {lattice['b']:.5f}, c= {lattice['c']:.5f}, "
        f"alpha= {lattice['alpha']:.5f}, beta= {lattice['beta']:.5f}, "
        f"gamma= {lattice['gamma']:.5f}"
    )
    wy = ("," + br).join(wyckoff_lines)
    return (
        f"Space Group: {sg}{br}"
        f"Lattice parameters: {lat}{br}"
        f"Default space-group preferences: {preferences}{br}"
        f"{wy}"
    )
