"""Audit archived official ISODISTORT Method 4 runs without modifying them.

The archive under ``output_compare`` is protected and read-only.  This tool
hashes and parses it, compares the selected transformation and exported
structure with the frozen daughter input, and writes a machine report only to
``ISODISTORT/output/validation``.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import sys
from collections import Counter
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from isocore.structure.cif_io import read_cif  # noqa: E402

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method4_download_manifest.json"
LOCAL_REPORT_PATH = PACKAGE_ROOT / "output" / "validation" / "method4_local_validation.json"
OFFICIAL_ROOT = WORKSPACE / "output_compare"
DEFAULT_REPORT = PACKAGE_ROOT / "output" / "validation" / "method4_official_audit.json"
SOURCE_PATHS = [
    Path(__file__).resolve(),
    MANIFEST_PATH,
    LOCAL_REPORT_PATH,
    PACKAGE_ROOT / "isocore" / "api" / "core_api.py",
    PACKAGE_ROOT / "isocore" / "distortion" / "search_methods.py",
    PACKAGE_ROOT / "isocore" / "distortion" / "strain.py",
    PACKAGE_ROOT / "isocore" / "structure" / "cif_io.py",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_signature(paths: list[Path]) -> dict[str, Any]:
    records = []
    combined = hashlib.sha256()
    for path in paths:
        digest = _sha256(path)
        try:
            label = path.resolve().relative_to(WORKSPACE.resolve()).as_posix()
        except ValueError:
            label = str(path.resolve())
        records.append({"path": label, "sha256": digest})
        combined.update(label.encode("utf-8"))
        combined.update(b"\0")
        combined.update(digest.encode("ascii"))
        combined.update(b"\n")
    return {"sha256": combined.hexdigest(), "files": records}


def _html_text(path: Path) -> str:
    raw = path.read_text(encoding="utf-8", errors="replace")
    raw = re.sub(r"(?i)<br\s*/?>", "\n", raw)
    raw = re.sub(r"(?i)</(?:p|pre|div|tr|h\d)>", "\n", raw)
    raw = re.sub(r"(?s)<[^>]+>", "", raw)
    return html.unescape(raw).replace("\xa0", " ")


def _expected_failure_matches(expectation: str, page_text: str) -> bool:
    """Recognize an expected scientific rejection, never a generic server error."""
    patterns = {
        "reject_species_or_stoichiometry": (
            r"types? of atoms?.*do not match",
            r"species.*(?:do not match|mismatch)",
            r"stoichiometr(?:y|ic).*(?:do not match|mismatch)",
        ),
        "reject_incompatible_lattice": (
            r"lattice.*(?:incompatible|not compatible|does not match)",
            r"(?:incompatible|not compatible).*lattice",
        ),
        "reject_unmatched_atom": (
            r"failed to find match",
            r"larger value for dmax",
            r"unable to match.*atom",
        ),
    }
    return any(
        re.search(pattern, page_text, flags=re.IGNORECASE | re.DOTALL)
        for pattern in patterns.get(expectation, ())
    )


def _number(token: str) -> float:
    return float(Fraction(token.strip()))


def _parse_vector(text: str) -> list[float]:
    return [_number(token) for token in text.split(",")]


def _parse_basis(text: str) -> list[list[float]]:
    rows = re.findall(r"\(([^()]*)\)", text)
    if len(rows) != 3:
        raise ValueError(f"Expected three basis vectors, found {rows!r}")
    return [_parse_vector(row) for row in rows]


def _parse_subgroup(text: str) -> dict[str, Any]:
    patterns = [
        r"Subgroup:\s*(\d+)\s+([^,\n]+),\s*basis=\{([^}]+)\},\s*"
        r"origin=\(([^)]+)\),\s*s=(\d+),\s*i=(\d+)",
        r"Subgroup details\s*(\d+)\s+([^,\n]+),\s*basis=\{([^}]+)\},\s*"
        r"origin=\(([^)]+)\),\s*s=(\d+),\s*i=(\d+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return {
                "space_group_number": int(match.group(1)),
                "space_group_symbol": match.group(2).strip(),
                "basis": _parse_basis(match.group(3)),
                "origin": _parse_vector(match.group(4)),
                "s": int(match.group(5)),
                "i": int(match.group(6)),
            }
    raise ValueError("Could not parse official Method 4 subgroup identity")


def _amplitude_rows(details_text: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    displacement_block = details_text.split("Displacive mode amplitudes", 1)[1]
    displacement_block, strain_tail = displacement_block.split(
        "Parent-cell strain mode definitions", 1,
    )
    displacement_rows = []
    for line in displacement_block.splitlines():
        match = re.match(
            r"^\s*(\[[^\n]+?\))\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
            line,
        )
        if match and "  all" not in line:
            displacement_rows.append({
                "label": match.group(1).strip(),
                "As": float(match.group(2)),
                "Ap": float(match.group(3)),
                "dmax_angstrom": float(match.group(4)),
            })
    strain_rows = []
    if "Parent-cell strain mode amplitudes" in strain_tail:
        strain_block = strain_tail.split("Parent-cell strain mode amplitudes", 1)[1]
        for line in strain_block.splitlines():
            match = re.match(r"^\s*(\[[^\n]+?\))\s+([-+0-9.Ee]+)\s*$", line)
            if match:
                strain_rows.append({
                    "label": match.group(1).strip(),
                    "amplitude": float(match.group(2)),
                })
    return displacement_rows, strain_rows


def _definition_normfactors(details_text: str) -> dict[str, float]:
    block = details_text.split("Displacive mode definitions", 1)[1]
    block = block.split("Displacive mode amplitudes", 1)[0]
    result = {}
    for line in block.splitlines():
        match = re.match(
            r"^\s*[^\s\[]+(\[[^\n]+?\))\s+normfactor\s*=\s*([-+0-9.Ee]+)\s*$",
            line,
        )
        if match:
            result[match.group(1).strip()] = float(match.group(2))
    return result


def _rounding_half_unit(token: str) -> float:
    mantissa, _, exponent_text = token.lower().partition("e")
    decimals = len(mantissa.rsplit(".", 1)[1]) if "." in mantissa else 0
    exponent = int(exponent_text) if exponent_text else 0
    return 0.5 * 10.0 ** (exponent - decimals)


def _isoviz_applied_strain_with_rounding(
    text: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the six Voigt components actually applied by IsoVIZ.

    The displayed/CIF strain amplitudes are symmetry-mode coordinates.  The
    viewer applies the following six-component mode vector before converting
    engineering shear components to a symmetric 3x3 matrix, so comparing the
    amplitude numbers alone is not a tensor comparison.
    """

    marker = "!strainmodelist"
    if marker not in text:
        raise ValueError("IsoVIZ export lacks a strain-mode list")
    block = text.split(marker, 1)[1]
    block = block.split("!displacivemodelist", 1)[0]
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    applied = np.zeros(6, dtype=float)
    rounding_bound = np.zeros(6, dtype=float)
    index = 0
    parsed = 0
    while index < len(lines):
        header = re.match(
            r"^\d+\s+([-+0-9.Ee]+)\s+[-+0-9.Ee]+\s+\d+\s+\S+\s*$",
            lines[index],
        )
        if not header:
            index += 1
            continue
        if index + 1 >= len(lines):
            raise ValueError("IsoVIZ strain mode is missing its six-component vector")
        vector_tokens = lines[index + 1].split()
        if len(vector_tokens) != 6:
            raise ValueError("IsoVIZ strain mode vector does not have six components")
        amplitude_token = header.group(1)
        amplitude = float(amplitude_token)
        amplitude_error = _rounding_half_unit(amplitude_token)
        vector = np.asarray([float(token) for token in vector_tokens], dtype=float)
        vector_error = np.asarray(
            [_rounding_half_unit(token) for token in vector_tokens], dtype=float
        )
        applied += amplitude * vector
        rounding_bound += (
            np.abs(vector) * amplitude_error
            + abs(amplitude) * vector_error
            + amplitude_error * vector_error
        )
        parsed += 1
        index += 2
    if parsed == 0:
        raise ValueError("IsoVIZ strain-mode list contains no parseable modes")
    return applied, rounding_bound


def _isoviz_applied_strain(text: str) -> np.ndarray:
    return _isoviz_applied_strain_with_rounding(text)[0]


def _canonical_wavevector(text: str) -> str:
    values = [Fraction(token.strip()) for token in text.split(",")]
    if len(values) != 3:
        raise ValueError(f"Expected a three-component wave vector, got {text!r}")
    if all(value == 0 for value in values):
        return ""
    return ",".join(str(value) for value in values)


def _local_mode_identity(label: str) -> tuple[str, str, str] | None:
    parts = label.split("__")
    if len(parts) < 3:
        return None
    ir_match = re.match(r"^(.+)\[([^]]+)\]$", parts[0])
    if ir_match:
        ir = ir_match.group(1)
        wavevector = _canonical_wavevector(ir_match.group(2))
    else:
        ir = parts[0]
        wavevector = ""
    tail_match = re.match(r"^(.+)\(([^()]*)\)$", parts[-1])
    if tail_match:
        site_ir, component = tail_match.groups()
    else:
        site_ir, component = "", parts[-1]
    # At Gamma the local subgroup-restricted basis and the official P1 basis
    # can use different site-irrep names for the same Cartesian column, so the
    # component is the stable identity.  Away from Gamma, retain both k and the
    # site-irrep name: LD1 harmonics and repeated A1_1/A1_2 branches otherwise
    # collide even though they are distinct physical modes.
    ir_key = f"{ir}[{wavevector}]" if wavevector else ir
    component_key = f"{site_ir}({component})" if wavevector else component
    return ir_key, parts[1], component_key


def _official_mode_identity(label: str) -> tuple[str, str, str] | None:
    match = re.match(
        r"^\[([^]]+)\]([^(]+)\([^)]*\)\[[^:\]]+:([^:\]]+):dsp\]"
        r"([^()]+)\(([^)]+)\)$",
        label,
    )
    if not match:
        return None
    wavevector = _canonical_wavevector(match.group(1))
    ir_key = f"{match.group(2)}[{wavevector}]" if wavevector else match.group(2)
    component_key = (
        f"{match.group(4)}({match.group(5)})" if wavevector else match.group(5)
    )
    return ir_key, match.group(3), component_key


def _parameter_mode_subspace_key(
    identity: tuple[str, str, str] | None,
) -> tuple[str, str, str] | None:
    """Return the basis-independent branch identity for a parameter-k mode.

    Equivalent multidimensional irrep bases may rotate components such as
    ``A2u(a)`` and ``A2u(b)``. Their individual amplitudes are basis
    dependent, while the Euclidean norm of a complete, equally normalized
    subspace is invariant. Gamma-mode identities intentionally return None:
    their local identity deliberately omits the official site-irrep name and
    cannot establish a complete official subspace by itself.
    """
    if identity is None or "[" not in identity[0]:
        return None
    component_match = re.fullmatch(r"(.+)\(([^()]*)\)", identity[2])
    if component_match is None:
        return None
    return identity[0], identity[1], component_match.group(1)


def _gamma_site_subspace_key(
    identity: tuple[str, str, str] | None,
) -> tuple[str, str] | None:
    """Return the complete Gamma irrep/site block for basis-invariant checks.

    A repeated site representation need not have the same component names in
    the two programs.  For example, the four local columns ``a,b,c,d`` for the
    NdNiO2 oxygen ``f`` orbit are the same space that official ISODISTORT calls
    ``B3u(a,b) + B2u(a,b)``.  Individual columns are not identifiable from the
    labels alone, but a complete equally normalized block has invariant As/Ap
    norms.  Parameter-k modes deliberately keep their existing, stricter key.
    """
    if identity is None or "[" in identity[0]:
        return None
    return identity[0], identity[1]


def _complete_gamma_site_subspaces(
    local_rows: list[dict[str, Any]],
    official_rows: list[dict[str, Any]],
) -> list[tuple[tuple[str, str], list[dict[str, Any]], list[dict[str, Any]]]]:
    """Return complete ambiguous Gamma site blocks from both mode bases."""
    local_groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    official_groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in local_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _gamma_site_subspace_key(identity)
        if key is not None:
            local_groups.setdefault(key, []).append(row)
    for row in official_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _gamma_site_subspace_key(identity)
        if key is not None:
            official_groups.setdefault(key, []).append(row)

    complete = []
    for key, rows in local_groups.items():
        official = official_groups.get(key, [])
        # Use a block comparison only when one-to-one component labels are
        # genuinely ambiguous and both sides span the entire same-sized block.
        official_identities = [tuple(row["identity"]) for row in official]
        if (
            len(rows) > 1
            and len(rows) == len(official)
            and len(set(official_identities)) < len(official_identities)
        ):
            complete.append((key, rows, official))
    return complete


def _complete_required_amplitude_subspaces(
    rows: list[dict[str, Any]],
    required_labels: set[str] | None,
    official_identities: set[tuple[str, str, str]],
) -> list[list[dict[str, Any]]]:
    """Select only fully covered multidimensional subspaces for norm checks."""
    if not required_labels:
        return []

    local_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    required_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    official_groups: dict[
        tuple[str, str, str], set[tuple[str, str, str]]
    ] = {}
    for row in rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _parameter_mode_subspace_key(identity)
        if key is None:
            continue
        local_groups.setdefault(key, []).append(row)
        if row["local_label"] in required_labels:
            required_groups.setdefault(key, []).append(row)
    for identity in official_identities:
        key = _parameter_mode_subspace_key(identity)
        if key is not None:
            official_groups.setdefault(key, set()).add(identity)

    complete = []
    for key, required_rows in required_groups.items():
        required_identities = {tuple(row["identity"]) for row in required_rows}
        local_identities = {tuple(row["identity"]) for row in local_groups[key]}
        if (
            len(required_rows) > 1
            and required_identities == local_identities
            and required_identities == official_groups.get(key, set())
            and all(row["official_label"] is not None for row in required_rows)
        ):
            complete.append(required_rows)
    return complete


def _parameter_origin_phase_turns(
    subspace_key: tuple[str, str, str],
    origin_shift_parent_fractional: list[float] | tuple[float, ...],
) -> Fraction | None:
    """Return ``k dot origin_shift`` modulo one for a parameter-k branch."""
    match = re.search(r"\[([^]]+)\]$", subspace_key[0])
    if match is None:
        return None
    coordinates = [
        Fraction(value.strip()) for value in match.group(1).split(",")
    ]
    if len(coordinates) != len(origin_shift_parent_fractional):
        return None
    shift = [
        Fraction(str(float(value))).limit_denominator(1_000_000)
        for value in origin_shift_parent_fractional
    ]
    return sum(
        (coordinate * value for coordinate, value in zip(coordinates, shift, strict=True)),
        Fraction(0),
    ) % 1


def _complete_origin_phase_amplitude_subspaces(
    local_rows: list[dict[str, Any]],
    required_labels: set[str] | None,
    official_rows: list[dict[str, Any]],
    origin_shift_parent_fractional: list[float] | tuple[float, ...],
) -> list[
    tuple[
        tuple[str, str, str],
        list[dict[str, Any]],
        list[dict[str, Any]],
        Fraction,
    ]
]:
    """Return complete parameter branches rotated by an official origin shift.

    A parent-lattice origin shift changes a non-Gamma complex mode by the phase
    ``exp(-2*pi*i*k.dot(t))``. The real ``a,b`` components may therefore rotate
    or change sign even when both programs reconstruct the same daughter. The
    Euclidean norm is invariant, but it is safe to use only when every local
    component in the branch is an explicit case target and the full official
    branch is present in Complete modes details.
    """
    if not required_labels:
        return []

    local_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    official_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in local_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _parameter_mode_subspace_key(identity)
        if key is not None:
            local_groups.setdefault(key, []).append(row)
    for row in official_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _parameter_mode_subspace_key(identity)
        if key is not None:
            official_groups.setdefault(key, []).append(row)

    complete = []
    for key, rows in local_groups.items():
        official = official_groups.get(key, [])
        local_identities = {tuple(row["identity"]) for row in rows}
        official_identities = {tuple(row["identity"]) for row in official}
        phase_turns = _parameter_origin_phase_turns(
            key, origin_shift_parent_fractional,
        )
        if (
            official
            and phase_turns not in {None, Fraction(0)}
            and all(row["local_label"] in required_labels for row in rows)
            and len(local_identities) == len(rows)
            and len(official_identities) == len(official)
            and local_identities <= official_identities
            and all(row["official_label"] is not None for row in rows)
        ):
            complete.append((key, rows, official, phase_turns))
    return complete


def _species_counts(structure) -> Counter[str]:
    return Counter(site.species_string for site in structure)


def _compare_structures(expected_path: Path, actual_path: Path) -> dict[str, Any]:
    expected = read_cif(expected_path)
    actual = read_cif(actual_path)
    expected_metric = np.asarray(expected.lattice.metric_tensor, dtype=float)
    actual_metric = np.asarray(actual.lattice.metric_tensor, dtype=float)
    metric_residual = float(
        np.linalg.norm(expected_metric - actual_metric)
        / max(float(np.linalg.norm(expected_metric)), 1.0)
    )
    result: dict[str, Any] = {
        "expected_atom_count": len(expected),
        "actual_atom_count": len(actual),
        "expected_species": dict(_species_counts(expected)),
        "actual_species": dict(_species_counts(actual)),
        "lattice_metric_relative_residual": metric_residual,
        "direct_max_species_matched_distance_angstrom": None,
        "max_species_matched_distance_angstrom": None,
        "origin_aligned_root_sum_squared_distance_angstrom": None,
        "origin_aligned_fractional_shift": None,
    }
    if len(expected) != len(actual) or _species_counts(expected) != _species_counts(actual):
        return result
    assignments = [-1] * len(expected)
    direct_max_distance = 0.0
    expected_species = [site.species_string for site in expected]
    actual_species = [site.species_string for site in actual]
    for species in sorted(set(expected_species)):
        left = [index for index, value in enumerate(expected_species) if value == species]
        right = [index for index, value in enumerate(actual_species) if value == species]
        distances = expected.lattice.get_all_distances(
            np.asarray(expected.frac_coords)[left],
            np.asarray(actual.frac_coords)[right],
        )
        rows, columns = linear_sum_assignment(distances)
        for row, column in zip(rows, columns, strict=True):
            assignments[left[int(row)]] = right[int(column)]
            direct_max_distance = max(
                direct_max_distance, float(distances[row, column]),
            )
    result["direct_max_species_matched_distance_angstrom"] = direct_max_distance

    # Official automatic-origin detection can translate the complete P1 atom
    # list while preserving the same physical structure.  Compare every
    # species-compatible anchor translation and retain the globally optimal
    # species-wise Hungarian assignment.  This is a single global origin
    # shift, not an independent per-atom adjustment.
    expected_frac = np.asarray(expected.frac_coords, dtype=float)
    actual_frac = np.asarray(actual.frac_coords, dtype=float)
    best: tuple[float, float, list[int], np.ndarray] | None = None
    for species in sorted(set(expected_species)):
        expected_anchors = [
            index for index, value in enumerate(expected_species) if value == species
        ]
        actual_anchors = [
            index for index, value in enumerate(actual_species) if value == species
        ]
        for expected_anchor in expected_anchors:
            for actual_anchor in actual_anchors:
                shift = expected_frac[expected_anchor] - actual_frac[actual_anchor]
                # Refine the anchor translation to the least-squares global
                # origin shift for the resulting periodic assignment.  An
                # anchor alone forces one rounded CIF coordinate to match
                # exactly and can leave a larger artificial maximum residual.
                for _ in range(5):
                    shifted_actual = actual_frac + shift
                    paired: list[tuple[int, int]] = []
                    for trial_species in sorted(set(expected_species)):
                        left = [
                            index for index, value in enumerate(expected_species)
                            if value == trial_species
                        ]
                        right = [
                            index for index, value in enumerate(actual_species)
                            if value == trial_species
                        ]
                        distances = expected.lattice.get_all_distances(
                            expected_frac[left], shifted_actual[right],
                        )
                        rows, columns = linear_sum_assignment(distances)
                        paired.extend(
                            (left[int(row)], right[int(column)])
                            for row, column in zip(rows, columns, strict=True)
                        )
                    fractional_deltas = np.asarray([
                        expected_frac[left] - shifted_actual[right]
                        for left, right in paired
                    ])
                    fractional_deltas -= np.round(fractional_deltas)
                    correction = np.mean(fractional_deltas, axis=0)
                    shift = shift + correction
                    if float(np.linalg.norm(correction)) < 1e-14:
                        break
                shifted_actual = actual_frac + shift
                trial_assignments = [-1] * len(expected)
                squared_sum = 0.0
                max_distance = 0.0
                for trial_species in sorted(set(expected_species)):
                    left = [
                        index for index, value in enumerate(expected_species)
                        if value == trial_species
                    ]
                    right = [
                        index for index, value in enumerate(actual_species)
                        if value == trial_species
                    ]
                    distances = expected.lattice.get_all_distances(
                        expected_frac[left], shifted_actual[right],
                    )
                    rows, columns = linear_sum_assignment(distances)
                    for row, column in zip(rows, columns, strict=True):
                        value = float(distances[row, column])
                        trial_assignments[left[int(row)]] = right[int(column)]
                        squared_sum += value * value
                        max_distance = max(max_distance, value)
                candidate = (squared_sum, max_distance, trial_assignments, shift)
                if best is None or candidate[:2] < best[:2]:
                    best = candidate
    if best is not None:
        aligned_squared_sum, aligned_max, aligned_assignments, aligned_shift = best
        aligned_shift = aligned_shift - np.round(aligned_shift)
        result["max_species_matched_distance_angstrom"] = aligned_max
        result["origin_aligned_root_sum_squared_distance_angstrom"] = float(
            np.sqrt(aligned_squared_sum)
        )
        result["origin_aligned_fractional_shift"] = aligned_shift.tolist()
        result["assignment"] = aligned_assignments
    else:
        result["max_species_matched_distance_angstrom"] = direct_max_distance
        result["assignment"] = assignments
    return result


def _compare_exported_cif(
    expected_path: Path,
    actual_path: Path,
) -> tuple[dict[str, Any], float]:
    comparison = _compare_structures(expected_path, actual_path)
    cif_text = actual_path.read_text(encoding="utf-8", errors="replace")
    coordinate_rows = re.findall(
        r"(?m)^\s*\S+\s+\S+\s+\d+\s+\S+\s+"
        r"([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+",
        cif_text,
    )
    coordinate_decimals = (
        min(
            len(token.rsplit(".", 1)[1])
            for row in coordinate_rows
            for token in row
        )
        if coordinate_rows else None
    )
    expected_structure = read_cif(expected_path)
    coordinate_rounding_bound = (
        0.5 * 10.0 ** (-coordinate_decimals)
        * float(np.sum(np.linalg.norm(expected_structure.lattice.matrix, axis=1)))
        if coordinate_decimals is not None else 0.0
    )
    distance_tolerance = max(1e-4, 1.05 * coordinate_rounding_bound)
    comparison["official_fractional_coordinate_decimals"] = coordinate_decimals
    comparison["coordinate_rounding_bound_angstrom"] = coordinate_rounding_bound
    comparison["distance_tolerance_angstrom"] = distance_tolerance
    return comparison, distance_tolerance


def _amplitude_comparison_tolerance(
    case: dict[str, Any],
    structure_comparison: dict[str, Any],
) -> tuple[float, dict[str, Any]]:
    """Return a case-derived bound for local/official shared-mode amplitudes.

    Exact inputs are limited by the website's five-decimal amplitude display.
    For a deliberately noisy input, official P1 decomposition can select a
    slightly different global origin and distribute noise into additional
    modes.  For any unit-normalized shared mode, Cauchy--Schwarz bounds the
    amplitude change by the Cartesian root-sum-squared difference between the
    aligned official export and the frozen input.  This uses the actual
    archived structures instead of a material- or case-specific constant.
    """
    display_rounding = 5.1e-6
    structure_bound = 0.0
    if case.get("noise_sigma_fractional") is not None:
        value = structure_comparison.get(
            "origin_aligned_root_sum_squared_distance_angstrom"
        )
        if value is not None:
            structure_bound = float(value)
        else:
            maximum = structure_comparison.get(
                "max_species_matched_distance_angstrom"
            )
            count = structure_comparison.get("expected_atom_count")
            if maximum is not None and count is not None:
                structure_bound = float(maximum) * float(np.sqrt(int(count)))
    total = display_rounding + structure_bound
    return total, {
        "display_rounding_angstrom": display_rounding,
        "aligned_structure_rss_bound_angstrom": structure_bound,
        "total_angstrom": total,
        "derivation": (
            "display rounding only"
            if structure_bound == 0.0
            else "display rounding + aligned-export Cartesian RSS "
                 "(Cauchy-Schwarz bound for a unit-normalized mode)"
        ),
    }


def _file_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _checked(raw: str, name: str, value: str) -> bool:
    pattern = (
        rf"<input\b(?=[^>]*\bname=[\"']{re.escape(name)}[\"'])"
        rf"(?=[^>]*\bvalue=[\"']{re.escape(value)}[\"'])[^>]*\bchecked(?:=[\"'][^\"']*[\"'])?[^>]*>"
    )
    return re.search(pattern, raw, flags=re.IGNORECASE) is not None


def _audit_g05_auto_origin(
    directory: Path,
    expected_daughter: Path,
    primary_result: dict[str, Any],
) -> dict[str, Any]:
    """Audit the independent G05 automatic-origin run.

    The explicit-origin run remains the primary case record.  This secondary
    record proves that website origin discovery converges to the same physical
    decomposition; missing screenshots are evidence warnings, not numerical
    failures.
    """
    result: dict[str, Any] = {
        "directory": str(directory.resolve()),
        "status": "not_downloaded",
        "issues": [],
        "warnings": [],
    }
    if not directory.is_dir():
        result["issues"].append("auto-origin directory is missing")
        result["status"] = "fail"
        return result
    core = {
        "result_page": directory / "ISODISTORT_ distortion.html",
        "details_page": directory / "ISODISTORT_ complete modes details.html",
        "cif": directory / "subgroup.cif",
        "topas": directory / "topas.str",
        "isoviz": directory / "data.isoviz",
    }
    missing = [name for name, path in core.items() if not path.is_file()]
    if missing:
        result["issues"].append(f"missing auto-origin core file(s): {missing}")
        result["status"] = "fail"
        return result
    basis_page = directory / "ISODISTORT_ distorted structure (basis).html"
    files = {name: _file_record(path) for name, path in core.items()}
    if basis_page.is_file():
        files["basis_page"] = _file_record(basis_page)
    else:
        result["warnings"].append(
            "auto-origin basis/matching page was not archived"
        )
    result["files"] = files
    control_captures = [
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg", ".pdf"}
    ]
    result["dynamic_control_captures"] = [
        _file_record(path) for path in control_captures
    ]
    if not control_captures:
        result["warnings"].append(
            "auto-origin selection and completed controls lack a screenshot/PDF"
        )

    result_text = _html_text(core["result_page"])
    details_text = _html_text(core["details_page"])
    result_subgroup = _parse_subgroup(result_text)
    details_subgroup = _parse_subgroup(details_text)
    result["official_subgroup"] = result_subgroup
    result["details_subgroup"] = details_subgroup
    if result_subgroup != details_subgroup:
        result["issues"].append(
            "auto-origin result and details subgroup identities differ"
        )
    if result_subgroup != primary_result.get("official_subgroup"):
        result["issues"].append(
            "auto-origin and explicit-origin subgroup identities differ"
        )

    displacement_rows, strain_rows = _amplitude_rows(details_text)
    result["official_amplitudes"] = {
        "displacive_mode_count": len(displacement_rows),
        "strain_mode_count": len(strain_rows),
        "displacive": displacement_rows,
        "strain": strain_rows,
    }
    primary_amplitudes = primary_result.get("official_amplitudes", {})
    if displacement_rows != primary_amplitudes.get("displacive"):
        result["issues"].append(
            "auto-origin and explicit-origin displacive amplitudes differ"
        )
    if strain_rows != primary_amplitudes.get("strain"):
        result["issues"].append(
            "auto-origin and explicit-origin strain amplitudes differ"
        )

    structure_comparison, structure_tolerance = _compare_exported_cif(
        expected_daughter, core["cif"],
    )
    result["exported_cif_vs_frozen_daughter"] = structure_comparison
    if structure_comparison["expected_atom_count"] != structure_comparison["actual_atom_count"]:
        result["issues"].append("auto-origin CIF atom count differs from frozen daughter")
    if structure_comparison["expected_species"] != structure_comparison["actual_species"]:
        result["issues"].append("auto-origin CIF species differ from frozen daughter")
    if structure_comparison["lattice_metric_relative_residual"] > 1e-6:
        result["issues"].append("auto-origin CIF lattice differs from frozen daughter")
    distance = structure_comparison["max_species_matched_distance_angstrom"]
    if distance is None or distance > structure_tolerance:
        result["issues"].append(
            "auto-origin CIF does not reproduce the frozen daughter within its "
            f"coordinate-precision tolerance {structure_tolerance:.8g} angstrom: {distance}"
        )

    result["status"] = (
        "pass" if not result["issues"] and not result["warnings"]
        else "pass_with_warnings" if not result["issues"]
        else "fail"
    )
    return result


def _audit_case(
    parent_name: str,
    case_id: str,
    case: dict[str, Any],
    parent: dict[str, Any],
    readable_folder: str,
    local_case: dict[str, Any] | None,
) -> dict[str, Any]:
    directory = OFFICIAL_ROOT / parent_name / "官网" / "Method4" / readable_folder
    result: dict[str, Any] = {
        "parent": parent_name,
        "case_id": case_id,
        "directory": str(directory.resolve()),
        "status": "not_downloaded",
        "issues": [],
        "warnings": [],
    }
    if not directory.is_dir():
        return result
    basis_page = directory / "ISODISTORT_ distorted structure (basis).html"
    required = {
        "result_page": directory / "ISODISTORT_ distortion.html",
        "details_page": directory / "ISODISTORT_ complete modes details.html",
        "cif": directory / "subgroup.cif",
        "topas": directory / "topas.str",
        "isoviz": directory / "data.isoviz",
    }
    expected_rejection = str(case.get("expect", "")).startswith("reject_")
    if expected_rejection:
        result_page = required["result_page"]
        if not result_page.is_file():
            if any(path.is_file() for path in required.values()) or basis_page.is_file():
                result.update(
                    status="fail",
                    issues=["expected-failure archive is missing its official result/error page"],
                )
            return result

        result["files"] = {"result_page": _file_record(result_page)}
        for name, path in required.items():
            if name != "result_page" and path.is_file():
                result["files"][name] = _file_record(path)
        if basis_page.is_file():
            result["files"]["basis_page"] = _file_record(basis_page)

        uploaded_copies = [
            path for path in directory.glob("*.cif")
            if path.name.lower() != "subgroup.cif"
        ]
        result["uploaded_daughter_copies"] = [
            _file_record(path) for path in uploaded_copies
        ]
        expected_daughter = WORKSPACE / case["daughter_cif"]
        actual_hash = _sha256(expected_daughter)
        result["frozen_daughter"] = {
            "path": str(expected_daughter.resolve()),
            "sha256_expected": case["sha256"],
            "sha256_actual": actual_hash,
        }
        if actual_hash != case["sha256"]:
            result["issues"].append("frozen daughter CIF hash mismatch")
        if not uploaded_copies:
            result["warnings"].append(
                "The exact uploaded daughter CIF was not copied into the official case directory"
            )
        elif not any(_sha256(path) == case["sha256"] for path in uploaded_copies):
            result["issues"].append(
                "no archived uploaded daughter CIF matches the frozen SHA-256"
            )

        control_captures = [
            path for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg", ".pdf"}
        ]
        result["dynamic_control_captures"] = [
            _file_record(path) for path in control_captures
        ]
        if basis_page.is_file():
            basis_raw = basis_page.read_text(encoding="utf-8", errors="replace")
            result["input_page_controls"] = {
                "automatic_origin_checked_in_saved_html": _checked(
                    basis_raw, "chooseorigin", "false",
                ),
                "nearest_site_checked_in_saved_html": _checked(
                    basis_raw, "trynearest", "true",
                ),
                "robust_checked_in_saved_html": _checked(
                    basis_raw, "trynearest", "false",
                ),
            }
        if not control_captures:
            result["warnings"].append(
                "expected-failure controls/error message lack a screenshot/PDF"
            )

        page_text = _html_text(result_page)
        expectation = str(case["expect"])
        observed = _expected_failure_matches(expectation, page_text)
        success_artifacts = [
            name for name in ("details_page", "cif", "topas", "isoviz")
            if required[name].is_file()
        ]
        result["expected_rejection"] = {
            "expectation": expectation,
            "recognized": observed,
            "success_artifacts": success_artifacts,
        }
        if success_artifacts:
            result["issues"].append(
                "official Method 4 accepted a case frozen as an expected rejection; "
                f"success artifacts present: {success_artifacts}"
            )
        if not observed:
            result["issues"].append(
                f"official page does not contain the expected rejection semantics {expectation!r}"
            )
        result["local_case"] = local_case
        result["status"] = (
            "pass" if not result["issues"] and not result["warnings"]
            else "pass_with_warnings" if not result["issues"]
            else "fail"
        )
        return result

    # The preparation script intentionally creates all readable case folders
    # before any website work begins.  An empty/pre-created folder is therefore
    # not a failed download.  Once at least one core artifact appears, however,
    # a missing companion artifact is a genuinely incomplete archive.
    if not any(path.is_file() for path in required.values()):
        return result
    missing = [name for name, path in required.items() if not path.is_file()]
    if missing:
        result.update(status="fail", issues=[f"missing required file(s): {missing}"])
        return result
    result["files"] = {
        name: {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": _sha256(path)}
        for name, path in required.items()
    }
    if basis_page.is_file():
        result["files"]["basis_page"] = _file_record(basis_page)
    uploaded_copies = [
        path for path in directory.glob("*.cif")
        if path.name.lower() not in {"subgroup.cif"}
    ]
    result["uploaded_daughter_copies"] = [
        {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in uploaded_copies
    ]
    if not uploaded_copies:
        result["warnings"].append(
            "The exact uploaded daughter CIF was not copied into the official case directory"
        )

    expected_daughter = WORKSPACE / case["daughter_cif"]
    actual_hash = _sha256(expected_daughter)
    result["frozen_daughter"] = {
        "path": str(expected_daughter.resolve()),
        "sha256_expected": case["sha256"],
        "sha256_actual": actual_hash,
    }
    if actual_hash != case["sha256"]:
        result["issues"].append("frozen daughter CIF hash mismatch")
    if uploaded_copies and not any(_sha256(path) == case["sha256"] for path in uploaded_copies):
        result["issues"].append(
            "no archived uploaded daughter CIF matches the frozen SHA-256"
        )

    basis_text = _html_text(basis_page) if basis_page.is_file() else ""
    result_text = _html_text(required["result_page"])
    details_text = _html_text(required["details_page"])
    basis_raw = (
        basis_page.read_text(encoding="utf-8", errors="replace")
        if basis_page.is_file() else ""
    )
    result_subgroup = _parse_subgroup(result_text)
    details_subgroup = _parse_subgroup(details_text)
    result["official_subgroup"] = result_subgroup
    result["details_subgroup"] = details_subgroup
    if result_subgroup != details_subgroup:
        result["issues"].append("result and Complete modes details subgroup identities differ")

    expected_context = parent["resolved_contexts"][case["context"]]
    expected_basis = np.asarray(expected_context["basis_vectors"], dtype=float)
    official_basis = np.asarray(result_subgroup["basis"], dtype=float)
    if not np.allclose(expected_basis, official_basis, atol=1e-12, rtol=0):
        result["issues"].append(
            f"official basis {official_basis.tolist()} differs from frozen context "
            f"{expected_basis.tolist()}"
        )
    result["input_page_controls"] = {
        "automatic_origin_checked_in_saved_html": _checked(
            basis_raw, "chooseorigin", "false",
        ),
        "nearest_site_checked_in_saved_html": _checked(basis_raw, "trynearest", "true"),
        "robust_checked_in_saved_html": _checked(basis_raw, "trynearest", "false"),
    }
    control_captures = [
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg", ".pdf"}
    ]
    result["dynamic_control_captures"] = [
        {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in control_captures
    ]
    if not basis_page.is_file():
        if control_captures:
            result["warnings"].append(
                "basis/matching HTML was not archived; the completed-controls screenshot/PDF "
                "is the surviving input-form evidence"
            )
        else:
            result["issues"].append(
                "neither basis/matching HTML nor a completed-controls screenshot/PDF was archived"
            )
    if "Specify basis as:" in basis_text and not control_captures:
        result["warnings"].append(
            "Saved HTML contains server-default form values, not proof of dynamic manual basis input; "
            "no screenshot/PDF was archived; the accepted basis is instead proven by the "
            "result/details/export artifacts"
        )

    displacement_rows, strain_rows = _amplitude_rows(details_text)
    definition_normfactors = _definition_normfactors(details_text)
    result["official_amplitudes"] = {
        "displacive_mode_count": len(displacement_rows),
        "strain_mode_count": len(strain_rows),
        "displacive": displacement_rows,
        "strain": strain_rows,
    }
    isoviz_text = required["isoviz"].read_text(
        encoding="utf-8", errors="replace"
    )
    try:
        (
            official_applied_strain,
            official_strain_rounding_bound,
        ) = _isoviz_applied_strain_with_rounding(isoviz_text)
    except ValueError as exc:
        official_applied_strain = None
        official_strain_rounding_bound = None
        result["issues"].append(str(exc))
    result["official_applied_strain_voigt_engineering"] = (
        None
        if official_applied_strain is None
        else {
            label: float(value)
            for label, value in zip(
                ("xx", "yy", "zz", "2yz", "2xz", "2xy"),
                official_applied_strain,
                strict=True,
            )
        }
    )
    result["official_applied_strain_rounding_bound"] = (
        None
        if official_strain_rounding_bound is None
        else official_strain_rounding_bound.tolist()
    )
    if not displacement_rows:
        result["issues"].append("no displacive amplitude rows parsed")
    is_zero_case = case_id in {"G01-zero", "P01-zero-supercell"}
    if is_zero_case:
        nonzero = [
            row for row in displacement_rows
            if max(abs(row["As"]), abs(row["Ap"]), abs(row["dmax_angstrom"])) > 5e-6
        ]
        nonzero_strain = [row for row in strain_rows if abs(row["amplitude"]) > 5e-6]
        if nonzero or nonzero_strain:
            result["issues"].append(
                f"zero case has nonzero official amplitudes: displacive={nonzero}, strain={nonzero_strain}"
            )

    cif_text = required["cif"].read_text(encoding="utf-8", errors="replace")
    cif_displacive = re.search(r"_iso_displacivemode_number\s+(\d+)", cif_text)
    cif_strain = re.search(r"_iso_strainmode_number\s+(\d+)", cif_text)
    result["export_consistency"] = {
        "cif_displacive_mode_count": int(cif_displacive.group(1)) if cif_displacive else None,
        "cif_strain_mode_count": int(cif_strain.group(1)) if cif_strain else None,
        "topas_mode_parameter_count": len(re.findall(
            r"(?m)^\s*prm\s+!a\d+\s+[-+0-9.Ee]+", required["topas"].read_text(
                encoding="utf-8", errors="replace",
            ),
        )),
        "isoviz_has_displacive_modes": "!displacivemodelist" in isoviz_text,
    }
    if result["export_consistency"]["cif_displacive_mode_count"] != len(displacement_rows):
        result["issues"].append("CIF and Complete modes details displacive-mode counts differ")
    if result["export_consistency"]["cif_strain_mode_count"] != len(strain_rows):
        result["issues"].append("CIF and Complete modes details strain-mode counts differ")
    if result["export_consistency"]["topas_mode_parameter_count"] != len(displacement_rows):
        result["issues"].append("TOPAS and Complete modes details displacive-mode counts differ")
    if not result["export_consistency"]["isoviz_has_displacive_modes"]:
        result["issues"].append("IsoVIZ export lacks a displacive-mode list")

    structure_comparison, structure_distance_tolerance = _compare_exported_cif(
        expected_daughter, required["cif"],
    )
    result["exported_cif_vs_frozen_daughter"] = structure_comparison
    amplitude_tolerance, amplitude_tolerance_evidence = (
        _amplitude_comparison_tolerance(case, structure_comparison)
    )
    result["shared_mode_amplitude_tolerance"] = amplitude_tolerance_evidence
    if structure_comparison["expected_atom_count"] != structure_comparison["actual_atom_count"]:
        result["issues"].append("exported CIF atom count differs from frozen daughter")
    if structure_comparison["expected_species"] != structure_comparison["actual_species"]:
        result["issues"].append("exported CIF species differ from frozen daughter")
    lattice_metric_tolerance = max(
        1e-6,
        4.0 * float(np.linalg.norm(official_strain_rounding_bound))
        if official_strain_rounding_bound is not None else 1e-6,
    )
    structure_comparison["lattice_metric_relative_tolerance"] = (
        lattice_metric_tolerance
    )
    if (
        structure_comparison["lattice_metric_relative_residual"]
        > lattice_metric_tolerance
    ):
        result["issues"].append(
            "exported CIF lattice differs from frozen daughter beyond the "
            "archived strain-vector rounding bound"
        )
    distance = structure_comparison["max_species_matched_distance_angstrom"]
    if distance is None or distance > structure_distance_tolerance:
        result["issues"].append(
            "exported CIF does not reproduce the frozen daughter within its "
            f"coordinate-precision tolerance {structure_distance_tolerance:.8g} angstrom: "
            f"{distance}"
        )

    result["local_case"] = local_case
    if local_case is None or local_case.get("status") != "pass":
        result["issues"].append("matching local validation case is absent or not passing")
    else:
        local_strain = local_case.get("strain_voigt_engineering")
        if official_applied_strain is None:
            pass
        elif not isinstance(local_strain, dict):
            result["issues"].append(
                "local validation case lacks applied strain components"
            )
        else:
            labels = ("xx", "yy", "zz", "2yz", "2xz", "2xy")
            local_vector = np.asarray(
                [float(local_strain[label]) for label in labels]
            )
            difference = np.abs(local_vector - official_applied_strain)
            strain_tolerance = official_strain_rounding_bound + 5e-10
            result["local_vs_official_applied_strain"] = {
                "local": local_vector.tolist(),
                "official_isoviz": official_applied_strain.tolist(),
                "absolute_difference": difference.tolist(),
                "max_absolute_difference": float(np.max(difference)),
                "componentwise_tolerance": strain_tolerance.tolist(),
            }
            if bool(np.any(difference > strain_tolerance)):
                result["issues"].append(
                    "local applied strain differs from the official IsoVIZ tensor: "
                    f"component differences {difference.tolist()} exceed archived "
                    f"rounding bounds {strain_tolerance.tolist()}"
                )
    if (
        local_case is not None
        and local_case.get("status") == "pass"
        and definition_normfactors
    ):
        official_normfactor_rows: list[dict[str, Any]] = []
        official_by_identity: dict[
            tuple[str, str, str], list[dict[str, Any]]
        ] = {}
        for label, value in definition_normfactors.items():
            identity = _official_mode_identity(label)
            if identity is not None:
                official_row = {
                    "identity": list(identity),
                    "label": label,
                    "value": float(value),
                }
                official_normfactor_rows.append(official_row)
                official_by_identity.setdefault(identity, []).append(official_row)
        normalization_matches = []
        required_normfactor_labels = (
            set(case.get("contributions", {}))
            if case.get("context") == "parameter"
            else None
        )
        for local_label, local_value in local_case.get(
            "mode_normfactors_inverse_angstrom", {},
        ).items():
            identity = _local_mode_identity(local_label)
            candidates = official_by_identity.get(identity, []) if identity else []
            official = candidates[0] if len(candidates) == 1 else None
            row = {
                "local_label": local_label,
                "identity": list(identity) if identity is not None else None,
                "local_normfactor_inverse_angstrom": float(local_value),
                "official_label": official["label"] if official else None,
                "official_normfactor_inverse_angstrom": (
                    official["value"] if official else None
                ),
                "absolute_difference": (
                    abs(float(local_value) - official["value"])
                    if official else None
                ),
                "comparison_required": (
                    required_normfactor_labels is None
                    or local_label in required_normfactor_labels
                ),
                "comparison_method": "componentwise",
            }
            normalization_matches.append(row)

        gamma_normalization_comparisons = []
        gamma_normalization_labels: set[str] = set()
        for key, local_rows, official_rows in _complete_gamma_site_subspaces(
            normalization_matches, official_normfactor_rows,
        ):
            local_values = sorted(
                row["local_normfactor_inverse_angstrom"] for row in local_rows
            )
            official_values = sorted(row["value"] for row in official_rows)
            differences = np.abs(
                np.asarray(local_values, dtype=float)
                - np.asarray(official_values, dtype=float)
            )
            comparison = {
                "subspace_identity": list(key),
                "local_labels": [row["local_label"] for row in local_rows],
                "official_labels": [row["label"] for row in official_rows],
                "sorted_local_normfactors_inverse_angstrom": local_values,
                "sorted_official_normfactors_inverse_angstrom": official_values,
                "max_absolute_difference": float(np.max(differences)),
                "tolerance_inverse_angstrom": 5.1e-6,
            }
            gamma_normalization_comparisons.append(comparison)
            gamma_normalization_labels.update(
                row["local_label"] for row in local_rows
            )
            for row in local_rows:
                row["comparison_method"] = "complete_gamma_site_subspace"
                row["official_subspace_labels"] = comparison["official_labels"]
            if (
                any(row["comparison_required"] for row in local_rows)
                and comparison["max_absolute_difference"]
                > comparison["tolerance_inverse_angstrom"]
            ):
                result["issues"].append(
                    "normfactor multiset mismatch for complete Gamma site subspace "
                    f"{comparison['subspace_identity']!r}"
                )

        for row in normalization_matches:
            if (
                row["official_label"] is None
                and row["local_label"] not in gamma_normalization_labels
            ):
                result["issues"].append(
                    f"no unique official mode matches local mode {row['local_label']!r}"
                )
            elif (
                row["official_label"] is not None
                and row["comparison_required"]
                and row["absolute_difference"] > 5.1e-6
            ):
                result["issues"].append(
                    f"normfactor mismatch for {row['local_label']!r}: local "
                    f"{row['local_normfactor_inverse_angstrom']}, official "
                    f"{row['official_normfactor_inverse_angstrom']}"
                )
        result["overlapping_mode_normalization"] = normalization_matches
        result["complete_gamma_site_normalization_comparisons"] = (
            gamma_normalization_comparisons
        )
        official_amplitude_rows: list[dict[str, Any]] = []
        official_amplitudes_by_identity: dict[
            tuple[str, str, str], list[dict[str, Any]]
        ] = {}
        for official_row in displacement_rows:
            identity = _official_mode_identity(official_row["label"])
            if identity is None:
                continue
            indexed_row = {**official_row, "identity": list(identity)}
            official_amplitude_rows.append(indexed_row)
            official_amplitudes_by_identity.setdefault(identity, []).append(
                indexed_row
            )
        amplitude_matches = []
        required_amplitude_labels = (
            set(case.get("contributions", {}))
            if (
                case.get("context") == "parameter"
                or case.get("noise_sigma_fractional") is not None
            )
            else None
        )
        for local_label, local_as in local_case.get("actual_As_angstrom", {}).items():
            identity = _local_mode_identity(local_label)
            candidates = (
                official_amplitudes_by_identity.get(identity, [])
                if identity is not None else []
            )
            official_row = candidates[0] if len(candidates) == 1 else None
            local_ap = float(local_case.get("actual_Ap_angstrom", {})[local_label])
            row = {
                "local_label": local_label,
                "identity": list(identity) if identity is not None else None,
                "local_As_angstrom": float(local_as),
                "local_Ap_angstrom": local_ap,
                "official_label": official_row["label"] if official_row else None,
                "official_As_angstrom": official_row["As"] if official_row else None,
                "official_Ap_angstrom": official_row["Ap"] if official_row else None,
                "As_absolute_difference": (
                    abs(float(local_as) - official_row["As"])
                    if official_row else None
                ),
                "Ap_absolute_difference": (
                    abs(local_ap - official_row["Ap"])
                    if official_row else None
                ),
                "comparison_required": (
                    required_amplitude_labels is None
                    or local_label in required_amplitude_labels
                ),
                "comparison_method": "componentwise",
            }
            amplitude_matches.append(row)
        complete_subspaces = _complete_required_amplitude_subspaces(
            amplitude_matches,
            required_amplitude_labels,
            set(official_amplitudes_by_identity),
        )
        official_parameter_groups: dict[
            tuple[str, str, str], list[dict[str, Any]]
        ] = {}
        for official_row in official_amplitude_rows:
            identity = tuple(official_row["identity"])
            key = _parameter_mode_subspace_key(identity)
            if key is not None:
                official_parameter_groups.setdefault(key, []).append(official_row)

        comparison_groups: list[
            tuple[
                tuple[str, str, str],
                list[dict[str, Any]],
                list[dict[str, Any]],
                str,
                Fraction | None,
            ]
        ] = []
        for rows in complete_subspaces:
            key = _parameter_mode_subspace_key(tuple(rows[0]["identity"]))
            if key is not None:
                comparison_groups.append((
                    key,
                    rows,
                    official_parameter_groups[key],
                    "complete_subspace_norm",
                    None,
                ))

        expected_origin = np.asarray(
            expected_context.get("origin", [0.0, 0.0, 0.0]), dtype=float,
        )
        official_origin = np.asarray(result_subgroup["origin"], dtype=float)
        origin_shift = (official_origin - expected_origin).tolist()
        already_grouped_labels = {
            row["local_label"]
            for _, rows, _, _, _ in comparison_groups
            for row in rows
        }
        for key, rows, official_rows, phase_turns in (
            _complete_origin_phase_amplitude_subspaces(
                amplitude_matches,
                required_amplitude_labels,
                official_amplitude_rows,
                origin_shift,
            )
        ):
            if any(row["local_label"] in already_grouped_labels for row in rows):
                continue
            comparison_groups.append((
                key,
                rows,
                official_rows,
                "origin_phase_subspace_norm",
                phase_turns,
            ))
            already_grouped_labels.update(row["local_label"] for row in rows)

        subspace_labels = set(already_grouped_labels)
        subspace_comparisons = []
        for key, rows, official_rows, method, phase_turns in comparison_groups:
            for row in rows:
                row["comparison_method"] = method
                row["official_subspace_labels"] = [
                    official_row["label"] for official_row in official_rows
                ]
            local_as_norm = float(np.linalg.norm([row["local_As_angstrom"] for row in rows]))
            official_as_norm = float(
                np.linalg.norm([row["As"] for row in official_rows])
            )
            local_ap_norm = float(np.linalg.norm([row["local_Ap_angstrom"] for row in rows]))
            official_ap_norm = float(
                np.linalg.norm([row["Ap"] for row in official_rows])
            )
            comparison = {
                "subspace_identity": list(key),
                "local_labels": [row["local_label"] for row in rows],
                "official_labels": [row["label"] for row in official_rows],
                "local_dimension": len(rows),
                "official_dimension": len(official_rows),
                "local_As_norm_angstrom": local_as_norm,
                "official_As_norm_angstrom": official_as_norm,
                "As_norm_absolute_difference": abs(local_as_norm - official_as_norm),
                "local_Ap_norm_angstrom": local_ap_norm,
                "official_Ap_norm_angstrom": official_ap_norm,
                "Ap_norm_absolute_difference": abs(local_ap_norm - official_ap_norm),
                "tolerance_angstrom": amplitude_tolerance,
                "comparison_method": method,
                "origin_phase_turns": (
                    None if phase_turns is None else str(phase_turns)
                ),
            }
            subspace_comparisons.append(comparison)
            if max(
                comparison["As_norm_absolute_difference"],
                comparison["Ap_norm_absolute_difference"],
            ) > comparison["tolerance_angstrom"]:
                result["issues"].append(
                    "As/Ap complete-subspace norm mismatch for "
                    f"{comparison['subspace_identity']!r}: local "
                    f"({local_as_norm}, {local_ap_norm}), official "
                    f"({official_as_norm}, {official_ap_norm})"
                )

        for key, local_rows, official_rows in _complete_gamma_site_subspaces(
            amplitude_matches, official_amplitude_rows,
        ):
            for row in local_rows:
                row["comparison_method"] = "complete_gamma_site_subspace_norm"
                row["official_subspace_labels"] = [
                    official_row["label"] for official_row in official_rows
                ]
            local_as_norm = float(np.linalg.norm([
                row["local_As_angstrom"] for row in local_rows
            ]))
            official_as_norm = float(np.linalg.norm([
                row["As"] for row in official_rows
            ]))
            local_ap_norm = float(np.linalg.norm([
                row["local_Ap_angstrom"] for row in local_rows
            ]))
            official_ap_norm = float(np.linalg.norm([
                row["Ap"] for row in official_rows
            ]))
            comparison = {
                "subspace_identity": list(key),
                "local_labels": [row["local_label"] for row in local_rows],
                "official_labels": [row["label"] for row in official_rows],
                "local_As_norm_angstrom": local_as_norm,
                "official_As_norm_angstrom": official_as_norm,
                "As_norm_absolute_difference": abs(local_as_norm - official_as_norm),
                "local_Ap_norm_angstrom": local_ap_norm,
                "official_Ap_norm_angstrom": official_ap_norm,
                "Ap_norm_absolute_difference": abs(local_ap_norm - official_ap_norm),
                "tolerance_angstrom": amplitude_tolerance,
                "comparison_method": "complete_gamma_site_subspace_norm",
            }
            subspace_comparisons.append(comparison)
            subspace_labels.update(row["local_label"] for row in local_rows)
            if (
                required_amplitude_labels is None
                or any(
                    row["local_label"] in required_amplitude_labels
                    for row in local_rows
                )
            ) and max(
                comparison["As_norm_absolute_difference"],
                comparison["Ap_norm_absolute_difference"],
            ) > comparison["tolerance_angstrom"]:
                result["issues"].append(
                    "As/Ap complete Gamma-site subspace norm mismatch for "
                    f"{comparison['subspace_identity']!r}: local "
                    f"({local_as_norm}, {local_ap_norm}), official "
                    f"({official_as_norm}, {official_ap_norm})"
                )

        for row in amplitude_matches:
            if (
                row["official_label"] is None
                and row["local_label"] not in subspace_labels
            ):
                result["issues"].append(
                    "no unique official amplitude row matches local mode "
                    f"{row['local_label']!r}"
                )
            if (
                row["official_label"] is not None
                and row["comparison_required"]
                and row["local_label"] not in subspace_labels
                and max(row["As_absolute_difference"], row["Ap_absolute_difference"])
                > amplitude_tolerance
            ):
                result["issues"].append(
                    f"As/Ap mismatch for {row['local_label']!r}: local "
                    f"({row['local_As_angstrom']}, {row['local_Ap_angstrom']}), official "
                    f"({row['official_As_angstrom']}, {row['official_Ap_angstrom']}), "
                    f"tolerance {amplitude_tolerance}"
                )
        result["overlapping_mode_amplitudes"] = amplitude_matches
        result["complete_subspace_amplitude_comparisons"] = subspace_comparisons
    if is_zero_case:
        result["scientific_scope"] = (
            f"{case_id} proves the selected basis, atom mapping, zero-amplitude decomposition, "
            "normalization factors, and export consistency at the website's displayed precision. "
            "It does not test nonzero As/Ap values. The official P1 basis can span more modes "
            "than the local subgroup-restricted basis."
        )
    else:
        result["scientific_scope"] = (
            f"{case_id} compares signed nonzero As/Ap values when the official and local mode "
            "components share a fixed basis, and compares the invariant Euclidean norm when a "
            "complete multidimensional parameter-k or repeated Gamma site subspace is present. "
            "For a nontrivial parameter-k phase caused by the archived official origin, the "
            "complete targeted branch norm is compared because its real components may rotate "
            "or change sign. "
            "Normalization factors are checked for target modes; extra official P1 modes are "
            "retained but are not treated as local columns."
        )
    if case_id == "G05-origin-shift":
        auto_origin = _audit_g05_auto_origin(
            directory / "auto-origin", expected_daughter, result,
        )
        result["auto_origin_run"] = auto_origin
        result["issues"].extend(
            f"auto-origin: {issue}" for issue in auto_origin["issues"]
        )
        result["warnings"].extend(
            f"auto-origin: {warning}" for warning in auto_origin["warnings"]
        )
    result["status"] = (
        "pass" if not result["issues"] and not result["warnings"]
        else "pass_with_warnings" if not result["issues"]
        else "fail"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--local-report", type=Path, default=LOCAL_REPORT_PATH)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--parent")
    parser.add_argument("--case-id")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    local = json.loads(args.local_report.read_text(encoding="utf-8"))
    cases = []
    for parent_name, parent in manifest["parents"].items():
        if args.parent and parent_name != args.parent:
            continue
        for case_id, case in parent["cases"].items():
            if args.case_id and case_id != args.case_id:
                continue
            local_case = (
                local.get("parents", {}).get(parent_name, {}).get("cases", {}).get(case_id)
            )
            cases.append(_audit_case(
                parent_name,
                case_id,
                case,
                parent,
                manifest["readable_case_folders"][case_id],
                local_case,
            ))
    downloaded = [case for case in cases if case["status"] != "not_downloaded"]
    failures = [case for case in downloaded if case["status"] == "fail"]
    warnings = [case for case in downloaded if case["status"] == "pass_with_warnings"]
    report = {
        "schema": 1,
        "kind": "method4_official_archive_audit",
        "generated_at": datetime.now(UTC).isoformat(),
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": _sha256(args.manifest),
        "local_report": str(args.local_report.resolve()),
        "local_report_sha256": _sha256(args.local_report),
        "source_signature": _source_signature([
            Path(__file__).resolve(),
            args.manifest,
            args.local_report,
            *SOURCE_PATHS[3:],
        ]),
        "scope_note": (
            "Downloaded official runs are audited read-only. Zero cases validate routing and "
            "zero preservation, not nonzero amplitude normalization."
        ),
        "summary": {
            "selected_case_count": len(cases),
            "downloaded_case_count": len(downloaded),
            "passed_case_count": sum(case["status"] == "pass" for case in downloaded),
            "warning_case_count": len(warnings),
            "failed_case_count": len(failures),
            "not_downloaded_case_count": sum(
                case["status"] == "not_downloaded" for case in cases
            ),
            "status": (
                "fail" if failures
                else "pass_with_warnings" if warnings
                else "pass"
            ),
        },
        "cases": cases,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report["summary"], ensure_ascii=False))
    print(args.json_output.resolve())
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
