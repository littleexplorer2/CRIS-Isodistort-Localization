"""Format the parent-structure header to match the official ISODISTORT search page.

Official example (EuAl4, I4/mmm)::

    Eu1 2a (0,0,0),
    Al1 4d (0,1/2,1/4),
    Al2 4e (0,0,z), z= 0.38000
"""
from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

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
    body = (
        f"{species}{species_index} {int(multiplicity)}{letter} "
        f"({','.join(parts)})"
    )
    if free:
        # Official prints ``z= 0.38000`` (five decimals, one space after =).
        extras = ", ".join(f"{ax}= {val:.5f}" for ax, val in free)
        body = f"{body}, {extras}"
    return body


def format_wyckoff_sites(
    structure: Structure,
    wyckoff_sites: Sequence[dict],
) -> list[str]:
    """Format all parent Wyckoff sites in official order (one string per site)."""
    species_count: dict[str, int] = {}
    lines: list[str] = []
    for site in wyckoff_sites:
        species = str(site["species"])
        species_count[species] = species_count.get(species, 0) + 1
        idx = int(site["representative_index"])
        coords = structure[idx].frac_coords
        lines.append(format_wyckoff_site(
            species=species,
            species_index=species_count[species],
            multiplicity=int(site["multiplicity"]),
            letter=str(site["wyckoff_letter"]),
            frac_coords=coords,
        ))
    return lines


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
