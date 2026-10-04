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
import shlex
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path

import numpy as np
from pymatgen.core import Element, Structure

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
    symmform: str | None = None,
    parameters: dict[str, float] | None = None,
) -> str:
    """One official site token, e.g. ``Al2 4e (0,0,z), z= 0.38000``."""
    if symmform:
        parts = [part.strip() for part in str(symmform).split(",")]
        if len(parts) != 3:
            raise ValueError(f"invalid Wyckoff symmetry form: {symmform!r}")
        free = [
            (axis, float((parameters or {})[axis]))
            for axis in _AXIS_LETTERS
            if axis in (parameters or {})
        ]
    else:
        coords = [float(x) for x in frac_coords]
        parts = []
        free = []
        for axis, value in zip(_AXIS_LETTERS, coords, strict=True):
            if _near_special(value) is None:
                parts.append(axis)
                free.append((axis, value))
            else:
                parts.append(format_fixed_coord(value))
    name = (label or "").strip() or f"{species}{species_index}"
    body = f"{name} {int(multiplicity)}{letter} ({','.join(parts)})"
    if free:
        # Official uses a sign column: positive values have one leading blank
        # (``z= 0.38000``), while negatives have none (``z=-0.13885``).
        extras = ", ".join(f"{ax}={val: .5f}" for ax, val in free)
        body = f"{body}, {extras}"
    return body


def format_wyckoff_sites(
    structure: Structure,
    wyckoff_sites: Sequence[dict],
) -> list[str]:
    """Format parent Wyckoff sites (fallback when no CIF path is available)."""
    species_count: dict[str, int] = {}
    lines: list[str] = []
    ordered_sites = sorted(
        enumerate(wyckoff_sites),
        key=lambda item: (int(item[1].get("display_order", item[0])), item[0]),
    )
    for _, site in ordered_sites:
        species = str(site["species"])
        species_count[species] = species_count.get(species, 0) + 1
        idx = int(site["representative_index"])
        coords = site.get("standard_representative_frac_coords")
        if coords is None:
            coords = structure[idx].frac_coords
        label = site.get("display_label")
        lines.append(format_wyckoff_site(
            species=species,
            species_index=species_count[species],
            multiplicity=int(site["multiplicity"]),
            letter=str(site["wyckoff_letter"]),
            frac_coords=coords,
            label=str(label) if label else None,
            symmform=site.get("standard_representative_symmform"),
            parameters=site.get("standard_representative_parameters"),
        ))
    return lines


def parse_cif_atom_site_rows(cif_path: str | Path) -> list[dict]:
    """Read asymmetric ``_atom_site_*`` rows from a CIF in file order.

    Preserves ``_atom_site_label`` / ``_atom_site_type_symbol`` text as written
    (e.g. ``ND`` / ``NI``), without pymatgen element normalization.
    """
    text = Path(cif_path).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    rows: list[dict] = []
    line_index = 0
    while line_index < len(lines):
        if lines[line_index].strip().lower() != "loop_":
            line_index += 1
            continue
        line_index += 1
        tags: list[str] = []
        while line_index < len(lines):
            stripped = lines[line_index].strip()
            if not stripped.startswith("_"):
                break
            tags.append(stripped.split(maxsplit=1)[0])
            line_index += 1
        lower = [t.lower() for t in tags]
        if "_atom_site_fract_x" not in lower:
            # Continue scanning by lines. The former whole-file regex could
            # catastrophically backtrack on FullProf templates containing
            # many non-structural loops and underscore-rich values.
            continue
        tokens: list[str] = []
        while line_index < len(lines):
            stripped = lines[line_index].strip()
            lowered = stripped.lower()
            if (
                lowered == "loop_"
                or lowered.startswith("data_")
                or lowered.startswith("save_")
                or stripped.startswith("_")
            ):
                break
            line_index += 1
            if not stripped or stripped.startswith("#"):
                continue
            try:
                tokens.extend(shlex.split(stripped, comments=True, posix=True))
            except ValueError:
                # Keep permissive behavior for imperfect legacy CIF rows; the
                # numeric atom-site fields are validated below before use.
                tokens.extend(stripped.split())
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


def _element_symbol(value: object) -> str | None:
    """Normalize CIF/site spellings such as ``ND``, ``Fe3+`` or ``O1``."""
    match = re.match(r"[A-Za-z]+", str(value or "").strip())
    if match is None:
        return None
    letters = match.group(0)
    for width in (2, 1):
        if len(letters) < width:
            continue
        candidate = letters[:width].capitalize()
        try:
            return Element(candidate).symbol
        except ValueError:
            continue
    return None


def _orbit_frac_dist(
    frac: Sequence[float],
    structure: Structure,
    site: dict,
) -> float | None:
    """Minimum periodic distance to any atom in one physical Wyckoff orbit."""
    indices = [int(value) for value in (site.get("equivalent_indices") or ())]
    representative = int(site["representative_index"])
    if representative not in indices:
        indices.append(representative)
    distances = [
        _frac_dist(frac, structure[index].frac_coords)
        for index in indices
        if 0 <= index < len(structure)
    ]
    return min(distances) if distances else None


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
    matched_labels: list[tuple[int, str]] = []
    species_count: dict[str, int] = {}
    for atom in asu:
        atom_species = _element_symbol(atom.get("type_symbol") or atom.get("label"))
        if atom_species is None:
            continue
        best_i = None
        best_d = 1e9
        for i, site in enumerate(wyckoff_sites):
            if i in used:
                continue
            if _element_symbol(site.get("species")) != atom_species:
                continue
            d = _orbit_frac_dist(atom["frac"], structure, site)
            if d is None:
                continue
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
        matched_labels.append((best_i, display))
        display_coords = site.get("standard_representative_frac_coords")
        if display_coords is None:
            display_coords = atom["frac"]
        lines.append(format_wyckoff_site(
            species=species,
            species_index=species_count[species],
            multiplicity=int(site["multiplicity"]),
            letter=str(site["wyckoff_letter"]),
            frac_coords=display_coords,
            label=display or None,
            symmform=site.get("standard_representative_symmform"),
            parameters=site.get("standard_representative_parameters"),
        ))
    if len(lines) != len(wyckoff_sites):
        return None
    for order, (site_idx, label) in enumerate(matched_labels):
        # Preserve the source CIF's atom-site names and order for all export
        # formats. Pymatgen can normalize or reorder both (ND -> Nd1).
        wyckoff_sites[site_idx]["display_label"] = label
        wyckoff_sites[site_idx]["display_order"] = order
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
