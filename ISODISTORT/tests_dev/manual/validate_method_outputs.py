"""Read-only scientific audit of saved Method 1/2/3 output pairs.

Candidates are paired by crystallographic header identity rather than directory
order or the printed basis. Point-group-equivalent bases may use different rows
while generating the same lattice.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import math
import platform
import re
import sys
import time
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterator, Mapping
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from pymatgen.io.cif import CifFile, CifParser

WORKSPACE = Path(__file__).resolve().parents[3]
VALIDATOR_ROOT = WORKSPACE / "ISODISTORT_VALIDATE"
if str(VALIDATOR_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATOR_ROOT))

from isodistort_validate.compare_cif import compare_cif  # noqa: E402

DEFAULT_OUTPUT_ROOT = WORKSPACE / "output_compare"
OFFICIAL_DIR = "官网"
LOCAL_DIR = "现有网页版交互"
CORE_LOCAL_FILES = (
    "subgroup.cif",
    "data.isoviz",
    "topas.str",
    "Complete modes details.txt",
)
LIVE12_CHECKPOINT_SCHEMA = 2
ISO_DATA_BASENAMES = (
    "const.dat",
    "data_diperiodic.txt",
    "data_images.txt",
    "data_irreps.txt",
    "data_isotropy.txt",
    "data_little.txt",
    "data_magnetic.txt",
    "data_space.txt",
    "data_ssg.txt",
    "data_ssgmag.txt",
    "data_wyckoff.txt",
)
RUNTIME_DISTRIBUTIONS = ("numpy", "pymatgen", "spglib", "PyYAML")
LIVE12_GROUP_COUNT = 4


@dataclass(frozen=True)
class CandidateHeader:
    ir: str
    identity: str
    basis: tuple[tuple[Fraction, Fraction, Fraction], ...]
    origin: tuple[Fraction, Fraction, Fraction]
    parametric_k: bool


@dataclass(frozen=True)
class CifSettingComparison:
    """Precision-aware comparison of two CIFs in candidate-defined settings."""

    equivalent: bool
    conclusive: bool
    details: dict[str, Any]


class AmbiguousCandidateHeaderError(ValueError):
    """The CIF declares more than one candidate setting transform."""


@dataclass
class AuditIssue:
    severity: str
    code: str
    case: str
    candidate: str
    message: str


@dataclass
class AuditSummary:
    output_root: str
    groups: dict[str, dict[str, int]] = field(default_factory=dict)
    paired_candidates: int = 0
    direct_cif_matches: int = 0
    equivalent_setting_matches: int = 0
    basis_equivalent_pairs: int = 0
    basis_exact_pairs: int = 0
    bases_normalized_to_official: int = 0
    format_files_checked: int = 0
    mode_files_compared: int = 0
    mode_count_matches: int = 0
    issues: list[AuditIssue] = field(default_factory=list)

    @property
    def failures(self) -> int:
        return sum(issue.severity == "error" for issue in self.issues)

    @property
    def warnings(self) -> int:
        return sum(issue.severity == "warning" for issue in self.issues)


_MATRIX_RE = re.compile(r"\(([^()]*)\)")
_IR_RE = re.compile(r"^#\s*IR:\s*([^,\s]+)", re.MULTILINE)
_PARAMETER_RE = re.compile(
    r"^#\s*k point:.*\b([a-z])\s*=\s*[-+]?\d+(?:/\d+)?\s*$",
    re.MULTILINE,
)


def _fraction(value: str) -> Fraction:
    return Fraction(value.strip())


def _canonical_number(value: Any) -> str:
    """Format an integer/float/string as the exact short rational used by the site."""
    number = Fraction(str(value)).limit_denominator(10_000)
    return str(number.numerator) if number.denominator == 1 else f"{number.numerator}/{number.denominator}"


def _canonical_origin(origin: tuple[Fraction, Fraction, Fraction]) -> str:
    """Return an origin modulo integer parent-lattice translations."""
    return ",".join(_canonical_number(value % 1) for value in origin)


def _parse_vectors(value: str, expected: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    vectors: list[tuple[Fraction, Fraction, Fraction]] = []
    for raw in _MATRIX_RE.findall(value):
        parts = tuple(_fraction(item) for item in raw.split(","))
        if len(parts) == 3:
            vectors.append(parts)
    if len(vectors) != expected:
        raise ValueError(f"expected {expected} three-component vectors, got {len(vectors)}: {value}")
    return tuple(vectors)


def _candidate_line(text: str) -> str:
    candidates = [
        line.removeprefix("#").strip()
        for line in text.splitlines()
        if line.startswith("#") and "basis={" in line and ", origin=" in line
    ]
    if not candidates:
        raise ValueError("CIF header has no candidate line containing basis and origin")
    if len(candidates) != 1:
        raise AmbiguousCandidateHeaderError(
            f"CIF header has {len(candidates)} candidate basis/origin lines"
        )
    return candidates[0]


def parse_candidate_header(path: Path) -> CandidateHeader:
    text = path.read_text(encoding="utf-8-sig", errors="strict")
    line = _candidate_line(text)
    ir_matches = _IR_RE.findall(text)
    if not ir_matches:
        raise ValueError("CIF header has no IR line")
    if len(ir_matches) != 1:
        raise AmbiguousCandidateHeaderError(
            f"CIF header has {len(ir_matches)} IR declarations"
        )
    basis_matches = re.findall(r"basis=\{(.*?)\},\s*origin=", line)
    origin_matches = re.findall(r"origin=(\([^)]*\))", line)
    if len(basis_matches) != 1 or len(origin_matches) != 1:
        raise ValueError("candidate line has an invalid basis/origin")
    basis = _parse_vectors(basis_matches[0], 3)
    origin = _parse_vectors(origin_matches[0], 1)[0]
    identity = re.sub(r",?\s*basis=\{.*?\}(?=,\s*origin=)", "", line)
    identity = re.sub(r"origin=\([^)]*\)", f"origin=({_canonical_origin(origin)})", identity)
    identity = " ".join(identity.split())
    return CandidateHeader(
        ir=ir_matches[0],
        identity=f"{ir_matches[0]} | {identity}",
        basis=basis,
        origin=origin,
        parametric_k=_PARAMETER_RE.search(text) is not None,
    )


def _determinant(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _inverse(matrix: tuple[tuple[Fraction, ...], ...]) -> tuple[tuple[Fraction, ...], ...]:
    det = _determinant(matrix)
    if det == 0:
        raise ValueError("candidate basis is singular")
    a, b, c = matrix
    cofactors = (
        (b[1] * c[2] - b[2] * c[1], -(b[0] * c[2] - b[2] * c[0]), b[0] * c[1] - b[1] * c[0]),
        (-(a[1] * c[2] - a[2] * c[1]), a[0] * c[2] - a[2] * c[0], -(a[0] * c[1] - a[1] * c[0])),
        (a[1] * b[2] - a[2] * b[1], -(a[0] * b[2] - a[2] * b[0]), a[0] * b[1] - a[1] * b[0]),
    )
    return tuple(tuple(cofactors[col][row] / det for col in range(3)) for row in range(3))


def _multiply(
    first: tuple[tuple[Fraction, ...], ...],
    second: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(sum(first[row][k] * second[k][col] for k in range(3)) for col in range(3))
        for row in range(3)
    )


def bases_generate_same_lattice(
    first: tuple[tuple[Fraction, ...], ...],
    second: tuple[tuple[Fraction, ...], ...],
) -> bool:
    """Return whether two row bases differ by an integer unimodular matrix."""
    transform = _multiply(first, _inverse(second))
    return all(value.denominator == 1 for row in transform for value in row) and abs(_determinant(transform)) == 1


_CIF_NUMBER_RE = re.compile(
    r"^\s*(?P<mantissa>[+-]?(?:\d+(?:\.\d*)?|\.\d+))"
    r"(?:\(\d+\))?(?:[Ee](?P<exponent>[+-]?\d+))?\s*$"
)
_CELL_TAGS = (
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma",
)
_COORDINATE_TAGS = (
    "_atom_site_fract_x",
    "_atom_site_fract_y",
    "_atom_site_fract_z",
)
_OCCUPANCY_TAG = "_atom_site_occupancy"
_SPECIES_TAG = "_atom_site_type_symbol"
_SPACE_GROUP_OPERATION_TAGS = (
    "_space_group_symop_operation_xyz",
    "_symmetry_equiv_pos_as_xyz",
)
_MAGNETIC_MOMENT_PREFIX = "_atom_site_moment"
_GRAM_LABELS = ("a", "b", "c")
_SPACE_GROUP_NUMBER_TAGS = (
    "_symmetry_Int_Tables_number",
    "_space_group_IT_number",
)
_SPACE_GROUP_SYMBOL_TAGS = (
    "_symmetry_space_group_name_H-M",
    "_space_group_name_H-M_alt",
)


@dataclass(frozen=True)
class _CifPrecision:
    cell_values: tuple[float, ...]
    cell_quantums: tuple[float, ...]
    coordinate_quantums: tuple[float, float, float]
    occupancy_quantum: float


@dataclass(frozen=True)
class _AtomSiteRow:
    """One unexpanded atom-site loop row with its serialized precision."""

    species: str
    coordinates: tuple[float, float, float]
    coordinate_quantums: tuple[float, float, float]
    occupancy: float
    occupancy_quantum: float


@dataclass(frozen=True)
class _ExactSymmetryOperation:
    """A CIF Seitz operation in exact fractional-coordinate arithmetic."""

    rotation: tuple[tuple[Fraction, Fraction, Fraction], ...]
    translation: tuple[Fraction, Fraction, Fraction]


@dataclass(frozen=True)
class _CifBlock(Mapping[str, Any]):
    """Parsed scalar values plus the original CIF loop-column membership."""

    data: dict[str, Any]
    loops: tuple[frozenset[str], ...]

    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.data)

    def __len__(self) -> int:
        return len(self.data)


def _cif_block(path: Path) -> _CifBlock:
    blocks = CifFile.from_file(path).data
    if len(blocks) != 1:
        raise ValueError(f"expected one CIF data block, got {len(blocks)}")
    parsed = next(iter(blocks.values()))
    return _CifBlock(
        data=parsed.data,
        loops=tuple(frozenset(str(tag) for tag in loop) for loop in parsed.loops),
    )


def _cif_number_and_quantum(value: Any) -> tuple[float, float]:
    """Return a CIF number and the resolution of its final printed digit."""
    match = _CIF_NUMBER_RE.fullmatch(str(value))
    if match is None:
        raise ValueError(f"unsupported CIF numeric value: {value!r}")
    mantissa = match.group("mantissa")
    exponent = int(match.group("exponent") or 0)
    decimal_places = len(mantissa.partition(".")[2]) if "." in mantissa else 0
    number = float(f"{mantissa}e{exponent}")
    return number, 10.0 ** (exponent - decimal_places)


def _cif_number_fraction(value: Any) -> Fraction:
    """Return the exact nominal decimal carried by a CIF numeric token."""
    match = _CIF_NUMBER_RE.fullmatch(str(value))
    if match is None:
        raise ValueError(f"unsupported CIF numeric value: {value!r}")
    result = Fraction(match.group("mantissa"))
    exponent = int(match.group("exponent") or 0)
    return result * (Fraction(10) ** exponent)


def _tag_values(block: Mapping[str, Any], tag: str) -> list[Any]:
    if tag not in block:
        raise ValueError(f"CIF is missing required tag {tag}")
    values = block[tag]
    return values if isinstance(values, list) else [values]


def _cif_precision(block: Mapping[str, Any]) -> _CifPrecision:
    cell = [_cif_number_and_quantum(_tag_values(block, tag)[0]) for tag in _CELL_TAGS]
    coordinate_quantums = tuple(
        max(_cif_number_and_quantum(value)[1] for value in _tag_values(block, tag))
        for tag in _COORDINATE_TAGS
    )
    occupancy_values = block.get(_OCCUPANCY_TAG)
    if occupancy_values is None:
        occupancy_quantum = 0.0
    else:
        if not isinstance(occupancy_values, list):
            occupancy_values = [occupancy_values]
        occupancy_quantum = max(
            _cif_number_and_quantum(value)[1] for value in occupancy_values
        )
    return _CifPrecision(
        cell_values=tuple(item[0] for item in cell),
        cell_quantums=tuple(item[1] for item in cell),
        coordinate_quantums=coordinate_quantums,
        occupancy_quantum=occupancy_quantum,
    )


def _column_values(
    block: Mapping[str, Any],
    tag: str,
    row_count: int,
    *,
    default: Any | None = None,
) -> list[Any]:
    """Return a strict atom-loop column without silently recycling values."""
    if tag not in block:
        if default is None:
            raise ValueError(f"CIF is missing required atom-site tag {tag}")
        return [default] * row_count
    values = block[tag]
    values = values if isinstance(values, list) else [values]
    if len(values) != row_count:
        raise ValueError(
            f"CIF atom-site column {tag} has {len(values)} rows; expected {row_count}"
        )
    return values


def _atom_site_rows(block: Mapping[str, Any]) -> tuple[_AtomSiteRow, ...]:
    """Parse raw ASU rows, retaining every row's coordinate/occupancy precision."""
    if isinstance(block, _CifBlock):
        required_columns = {_SPECIES_TAG, *_COORDINATE_TAGS}
        if _OCCUPANCY_TAG in block:
            required_columns.add(_OCCUPANCY_TAG)
        atom_loops = [
            loop for loop in block.loops if _SPECIES_TAG in loop
        ]
        matching_loops = [
            loop for loop in atom_loops if required_columns.issubset(loop)
        ]
        if len(matching_loops) != 1:
            raise ValueError(
                "CIF must contain exactly one atom-site loop holding species, "
                "fractional coordinates, and occupancy"
            )
    species_values = _tag_values(block, _SPECIES_TAG)
    row_count = len(species_values)
    if row_count == 0:
        raise ValueError("CIF atom-site loop is empty")
    coordinate_columns = [
        _column_values(block, tag, row_count) for tag in _COORDINATE_TAGS
    ]
    occupancy_values = _column_values(
        block,
        _OCCUPANCY_TAG,
        row_count,
        default="1",
    )
    rows: list[_AtomSiteRow] = []
    nominal_site_occupancies: defaultdict[
        tuple[Fraction, Fraction, Fraction], Fraction
    ] = defaultdict(Fraction)
    for index, raw_species in enumerate(species_values):
        species = str(raw_species).strip().strip("'\"")
        if not species or species in {".", "?"}:
            raise ValueError(f"atom-site row {index + 1} has no concrete species")
        coordinate_items = tuple(
            _cif_number_and_quantum(column[index]) for column in coordinate_columns
        )
        occupancy, occupancy_quantum = _cif_number_and_quantum(
            occupancy_values[index]
        )
        if not np.isfinite(occupancy) or not 0.0 <= occupancy <= 1.0:
            raise ValueError(
                f"atom-site row {index + 1} has invalid occupancy {occupancy!r}"
            )
        coordinates = tuple(float(item[0] % 1.0) for item in coordinate_items)
        rows.append(
            _AtomSiteRow(
                species=species,
                coordinates=coordinates,
                coordinate_quantums=tuple(float(item[1]) for item in coordinate_items),
                occupancy=float(occupancy),
                occupancy_quantum=float(occupancy_quantum),
            )
        )
        coordinate_key = tuple(
            _cif_number_fraction(column[index]) % 1
            for column in coordinate_columns
        )
        nominal_site_occupancies[coordinate_key] += _cif_number_fraction(
            occupancy_values[index]
        )
    invalid_totals = {
        key: value
        for key, value in nominal_site_occupancies.items()
        if value > 1
    }
    if invalid_totals:
        preview = next(iter(invalid_totals.items()))
        raise ValueError(
            "atom-site rows at the same nominal coordinate exceed total occupancy 1: "
            f"{preview[0]} -> {preview[1]}"
        )
    return tuple(rows)


def _magnetic_moment_tags(block: Mapping[str, Any]) -> list[str]:
    return sorted(
        str(tag)
        for tag in block
        if str(tag).casefold().startswith(_MAGNETIC_MOMENT_PREFIX)
    )


def _parse_symmetry_component(
    expression: str,
) -> tuple[tuple[Fraction, Fraction, Fraction], Fraction]:
    """Parse one CIF x/y/z affine expression without float roundoff."""
    compact = (
        expression.strip().strip("'\"").replace(" ", "").replace("*", "").casefold()
    )
    if not compact:
        raise ValueError("empty symmetry-operation component")
    normalized = compact if compact[0] in "+-" else f"+{compact}"
    terms = re.findall(r"[+-][^+-]+", normalized)
    if "".join(terms) != normalized:
        raise ValueError(f"unsupported symmetry-operation component {expression!r}")
    coefficients = [Fraction(0), Fraction(0), Fraction(0)]
    translation = Fraction(0)
    for term in terms:
        variable_indices = [axis for axis, name in enumerate("xyz") if name in term]
        if variable_indices:
            if len(variable_indices) != 1 or sum(term.count(name) for name in "xyz") != 1:
                raise ValueError(f"nonlinear symmetry-operation term {term!r}")
            axis = variable_indices[0]
            coefficient_text = term.replace("xyz"[axis], "")
            if coefficient_text in {"+", ""}:
                coefficient = Fraction(1)
            elif coefficient_text == "-":
                coefficient = Fraction(-1)
            else:
                coefficient = Fraction(coefficient_text)
            coefficients[axis] += coefficient
        else:
            translation += Fraction(term)
    return tuple(coefficients), translation


def _space_group_operations(
    block: Mapping[str, Any],
) -> tuple[_ExactSymmetryOperation, ...]:
    present_tags = [tag for tag in _SPACE_GROUP_OPERATION_TAGS if tag in block]
    if not present_tags:
        raise ValueError("CIF has no explicit space-group operation loop")
    if len(present_tags) != 1:
        raise ValueError(
            "CIF has multiple alternative space-group operation declarations"
        )
    values = block[present_tags[0]]
    raw_values = values if isinstance(values, list) else [values]
    operations: list[_ExactSymmetryOperation] = []
    for raw_value in raw_values:
        components = [part.strip() for part in str(raw_value).strip("'\"").split(",")]
        if len(components) != 3:
            raise ValueError(f"invalid symmetry operation {raw_value!r}")
        parsed = [_parse_symmetry_component(component) for component in components]
        operations.append(
            _ExactSymmetryOperation(
                rotation=tuple(item[0] for item in parsed),
                translation=tuple(item[1] for item in parsed),
            )
        )
    if not operations:
        raise ValueError("CIF space-group operation loop is empty")
    return tuple(operations)


def _first_tag_value(
    block: Mapping[str, Any], tags: tuple[str, ...]
) -> Any | None:
    for tag in tags:
        if tag in block:
            value = block[tag]
            return value[0] if isinstance(value, list) else value
    return None


def _declared_space_group(block: Mapping[str, Any]) -> dict[str, Any]:
    number = _first_tag_value(block, _SPACE_GROUP_NUMBER_TAGS)
    symbol = _first_tag_value(block, _SPACE_GROUP_SYMBOL_TAGS)
    normalized_symbol = (
        " ".join(str(symbol).strip("'\"").split()) if symbol is not None else None
    )
    numeric_number = None
    if number is not None:
        parsed_number = _cif_number_and_quantum(number)[0]
        if not float(parsed_number).is_integer() or not 1 <= parsed_number <= 230:
            raise ValueError(f"invalid declared space-group IT number {number!r}")
        numeric_number = int(parsed_number)
    return {"number": numeric_number, "symbol": normalized_symbol}


def _normalized_space_group_symbol(symbol: Any) -> str | None:
    if symbol is None:
        return None
    return re.sub(r"[\s_]", "", str(symbol)).casefold()


def _compare_declared_space_groups(
    local: dict[str, Any],
    official: dict[str, Any],
) -> tuple[bool, str, list[str]]:
    """Compare declarations without treating H-M typography as physics."""
    local_number = local.get("number")
    official_number = official.get("number")
    local_symbol = _normalized_space_group_symbol(local.get("symbol"))
    official_symbol = _normalized_space_group_symbol(official.get("symbol"))
    symbol_equal = (
        local_symbol is not None
        and official_symbol is not None
        and local_symbol == official_symbol
    )
    diagnostics: list[str] = []
    if local_number is not None and official_number is not None:
        if local_number != official_number:
            return False, "different-it-number", diagnostics
        if not symbol_equal:
            diagnostics.append("declared-symbol-differs-with-same-it-number")
        return True, "same-it-number", diagnostics
    if symbol_equal:
        diagnostics.append("declared-it-number-missing")
        return True, "same-normalized-symbol", diagnostics
    diagnostics.append("declared-space-group-comparison-inconclusive")
    return True, "inconclusive-without-comparable-it-numbers", diagnostics


def _interval_product(*intervals: tuple[float, float]) -> tuple[float, float]:
    values = [1.0]
    for low, high in intervals:
        values = [value * endpoint for value in values for endpoint in (low, high)]
    return min(values), max(values)


def _gram_intervals(
    parameters: tuple[float, ...],
    quantums: tuple[float, ...],
) -> tuple[tuple[tuple[float, float], ...], ...]:
    ranges = tuple(
        (value - quantum / 2.0, value + quantum / 2.0)
        for value, quantum in zip(parameters, quantums, strict=True)
    )
    lengths = ranges[:3]
    angle_cosines: list[tuple[float, float]] = []
    for low, high in ranges[3:]:
        if not 0.0 < low <= high < 180.0:
            raise ValueError(f"invalid CIF cell-angle interval [{low}, {high}]")
        endpoints = (np.cos(np.deg2rad(low)), np.cos(np.deg2rad(high)))
        angle_cosines.append((float(min(endpoints)), float(max(endpoints))))

    gram: list[list[tuple[float, float] | None]] = [[None] * 3 for _ in range(3)]
    for index in range(3):
        gram[index][index] = _interval_product(lengths[index], lengths[index])
    for first, second, angle_index in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        interval = _interval_product(
            lengths[first],
            lengths[second],
            angle_cosines[angle_index],
        )
        gram[first][second] = interval
        gram[second][first] = interval
    return tuple(
        tuple(item for item in row if item is not None)
        for row in gram
    )


def _transform_gram_intervals(
    gram: tuple[tuple[tuple[float, float], ...], ...],
    transform: np.ndarray,
) -> tuple[tuple[tuple[float, float], ...], ...]:
    result: list[list[tuple[float, float]]] = []
    for row in range(3):
        output_row: list[tuple[float, float]] = []
        for column in range(3):
            low = 0.0
            high = 0.0
            terms: list[tuple[float, tuple[float, float]]] = [
                (
                    transform[row, index] * transform[column, index],
                    gram[index][index],
                )
                for index in range(3)
            ]
            terms.extend(
                (
                    transform[row, first] * transform[column, second]
                    + transform[row, second] * transform[column, first],
                    gram[first][second],
                )
                for first, second in ((0, 1), (0, 2), (1, 2))
            )
            for coefficient, (source_low, source_high) in terms:
                if coefficient >= 0:
                    low += coefficient * source_low
                    high += coefficient * source_high
                else:
                    low += coefficient * source_high
                    high += coefficient * source_low
            output_row.append((float(low), float(high)))
        result.append(output_row)
    return tuple(tuple(row) for row in result)


def _gram_matrix(parameters: tuple[float, ...]) -> np.ndarray:
    a, b, c, alpha, beta, gamma = parameters
    cos_alpha, cos_beta, cos_gamma = np.cos(np.deg2rad((alpha, beta, gamma)))
    return np.array(
        (
            (a * a, a * b * cos_gamma, a * c * cos_beta),
            (a * b * cos_gamma, b * b, b * c * cos_alpha),
            (a * c * cos_beta, b * c * cos_alpha, c * c),
        ),
        dtype=float,
    )


def _cell_parameters_from_gram(gram: np.ndarray) -> np.ndarray:
    """Recover the conventional six cell parameters from a positive Gram matrix."""
    symmetric = 0.5 * (np.asarray(gram, dtype=float) + np.asarray(gram, dtype=float).T)
    eigenvalues = np.linalg.eigvalsh(symmetric)
    if not np.all(np.isfinite(eigenvalues)) or float(np.min(eigenvalues)) <= 0.0:
        raise ValueError("transformed Gram matrix is not positive definite")
    a, b, c = (math.sqrt(float(symmetric[index, index])) for index in range(3))

    def angle(value: float) -> float:
        return math.degrees(math.acos(max(-1.0, min(1.0, value))))

    return np.asarray(
        (
            a,
            b,
            c,
            angle(float(symmetric[1, 2]) / (b * c)),
            angle(float(symmetric[0, 2]) / (a * c)),
            angle(float(symmetric[0, 1]) / (a * b)),
        ),
        dtype=float,
    )


def _joint_metric_witness(
    local: _CifPrecision,
    official: _CifPrecision,
    transform: np.ndarray,
) -> dict[str, Any]:
    """Find a single pair of rounded cells satisfying G_l = U G_o U.T.

    A concrete witness proves compatibility. Failure to find one is deliberately
    not treated as proof of incompatibility; callers report it as inconclusive.
    """
    try:
        from scipy.optimize import least_squares
    except ImportError as exc:
        return {
            "found": False,
            "method": "solver-unavailable",
            "error": f"{type(exc).__module__}.{type(exc).__qualname__}: {exc}",
        }

    local_center = np.asarray(local.cell_values, dtype=float)
    local_half = 0.5 * np.asarray(local.cell_quantums, dtype=float)
    official_center = np.asarray(official.cell_values, dtype=float)
    official_half = 0.5 * np.asarray(official.cell_quantums, dtype=float)
    local_low = local_center - local_half
    local_high = local_center + local_half
    official_low = official_center - official_half
    official_high = official_center + official_half
    arithmetic_scale = np.maximum(1.0, np.abs(local_center))
    arithmetic_slack = 256.0 * np.finfo(float).eps * arithmetic_scale

    def derived_local(parameters: np.ndarray) -> np.ndarray:
        gram = _gram_matrix(tuple(float(value) for value in parameters))
        return _cell_parameters_from_gram(transform @ gram @ transform.T)

    local_gram_intervals = _gram_intervals(
        local.cell_values,
        local.cell_quantums,
    )
    official_gram_intervals = _gram_intervals(
        official.cell_values,
        official.cell_quantums,
    )

    def witness_payload(
        official_parameters: np.ndarray,
        local_parameters: np.ndarray,
    ) -> dict[str, Any]:
        official_gram = _gram_matrix(
            tuple(float(value) for value in official_parameters)
        )
        transformed_gram = transform @ official_gram @ transform.T
        reconstructed_local_gram = _gram_matrix(
            tuple(float(value) for value in local_parameters)
        )
        component_names: list[str] = []
        official_components: list[float] = []
        local_components: list[float] = []
        official_inside: list[bool] = []
        local_inside: list[bool] = []
        for row in range(3):
            for column in range(row, 3):
                component_names.append(
                    f"{_GRAM_LABELS[row]}{_GRAM_LABELS[column]}"
                )
                official_value = float(official_gram[row, column])
                local_value = float(transformed_gram[row, column])
                official_components.append(official_value)
                local_components.append(local_value)
                official_overlap, _ = _intervals_overlap(
                    (official_value, official_value),
                    official_gram_intervals[row][column],
                )
                local_overlap, _ = _intervals_overlap(
                    (local_value, local_value),
                    local_gram_intervals[row][column],
                )
                official_inside.append(official_overlap)
                local_inside.append(local_overlap)
        reconstruction_error = float(
            np.max(np.abs(transformed_gram - reconstructed_local_gram))
        )
        scale = max(1.0, float(np.max(np.abs(transformed_gram))))
        reconstruction_ok = reconstruction_error <= (
            512.0 * np.finfo(float).eps * scale
        )
        return {
            "gram_component_order": component_names,
            "official_gram_components": official_components,
            "local_gram_components": local_components,
            "official_components_inside_printed_box_intervals": official_inside,
            "local_components_inside_printed_box_intervals": local_inside,
            "all_six_gram_components_jointly_validated": bool(
                all(official_inside) and all(local_inside) and reconstruction_ok
            ),
            "local_gram_reconstruction_max_abs_error": reconstruction_error,
        }

    def accepted(
        parameters: np.ndarray,
    ) -> tuple[bool, np.ndarray, dict[str, Any]]:
        derived = derived_local(parameters)
        parameter_boxes_ok = bool(
            np.all(parameters >= official_low - arithmetic_slack)
            and np.all(parameters <= official_high + arithmetic_slack)
            and np.all(derived >= local_low - arithmetic_slack)
            and np.all(derived <= local_high + arithmetic_slack)
        )
        payload = witness_payload(parameters, derived)
        return (
            parameter_boxes_ok
            and payload["all_six_gram_components_jointly_validated"],
            derived,
            payload,
        )

    accepted_nominal, nominal_local, nominal_payload = accepted(official_center)
    if accepted_nominal:
        return {
            "found": True,
            "method": "official-nominal",
            "official_cell_parameters": official_center.tolist(),
            "local_cell_parameters": nominal_local.tolist(),
            "max_local_interval_excess": 0.0,
            **nominal_payload,
        }

    scale = np.maximum(local_half, 64.0 * np.finfo(float).eps * arithmetic_scale)

    def residual(parameters: np.ndarray) -> np.ndarray:
        return (derived_local(parameters) - local_center) / scale

    starts = [
        official_center,
        official_low,
        official_high,
    ]
    best: tuple[float, np.ndarray, np.ndarray, Any] | None = None
    for start in starts:
        try:
            result = least_squares(
                residual,
                np.clip(start, official_low, official_high),
                bounds=(official_low, official_high),
                xtol=1e-14,
                ftol=1e-14,
                gtol=1e-14,
                max_nfev=4000,
            )
            candidate = np.asarray(result.x, dtype=float)
            ok, derived, candidate_payload = accepted(candidate)
            excess = np.maximum(local_low - derived, derived - local_high)
            maximum_excess = float(np.max(np.maximum(excess, 0.0)))
            if best is None or maximum_excess < best[0]:
                best = (maximum_excess, candidate, derived, result)
            if ok:
                return {
                    "found": True,
                    "method": "bounded-least-squares-witness",
                    "official_cell_parameters": candidate.tolist(),
                    "local_cell_parameters": derived.tolist(),
                    "max_local_interval_excess": maximum_excess,
                    **candidate_payload,
                    "optimizer": {
                        "success": bool(result.success),
                        "status": int(result.status),
                        "cost": float(result.cost),
                        "optimality": float(result.optimality),
                        "nfev": int(result.nfev),
                    },
                }
        except (FloatingPointError, ValueError):
            continue
    payload: dict[str, Any] = {
        "found": False,
        "method": "bounded-least-squares-no-witness",
    }
    if best is not None:
        payload.update(
            {
                "best_official_cell_parameters": best[1].tolist(),
                "best_derived_local_cell_parameters": best[2].tolist(),
                "best_max_local_interval_excess": best[0],
                "optimizer": {
                    "success": bool(best[3].success),
                    "status": int(best[3].status),
                    "cost": float(best[3].cost),
                    "optimality": float(best[3].optimality),
                    "nfev": int(best[3].nfev),
                },
            }
        )
    return payload


def _intervals_overlap(
    first: tuple[float, float],
    second: tuple[float, float],
) -> tuple[bool, float]:
    scale = max(1.0, *(abs(value) for value in first + second))
    arithmetic_slack = 128.0 * np.finfo(float).eps * scale
    gap = max(first[0], second[0]) - min(first[1], second[1])
    return gap <= arithmetic_slack, max(0.0, float(gap))


def _fraction_json(value: Fraction) -> int | str:
    return value.numerator if value.denominator == 1 else str(value)


def _fraction_matrix_json(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> list[list[int | str]]:
    return [[_fraction_json(value) for value in row] for row in matrix]


def _float_matrix(matrix: tuple[tuple[Fraction, ...], ...]) -> np.ndarray:
    return np.asarray([[float(value) for value in row] for row in matrix], dtype=float)


def _vector_multiply(
    vector: tuple[Fraction, ...],
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[Fraction, ...]:
    return tuple(
        sum(vector[index] * matrix[index][column] for index in range(3))
        for column in range(3)
    )


def _transpose(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(matrix[row][column] for row in range(3)) for column in range(3))


def _matrix_vector(
    matrix: tuple[tuple[Fraction, ...], ...],
    vector: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    )


def _vector_subtract(
    first: tuple[Fraction, ...],
    second: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    return tuple(left - right for left, right in zip(first, second, strict=True))


def _vector_add(
    first: tuple[Fraction, ...],
    second: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    return tuple(left + right for left, right in zip(first, second, strict=True))


def _operation_key(
    operation: _ExactSymmetryOperation,
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    return (
        tuple(value for row in operation.rotation for value in row),
        tuple(value % 1 for value in operation.translation),
    )


def _transform_operation_to_official(
    operation: _ExactSymmetryOperation,
    local_to_official: tuple[tuple[Fraction, ...], ...],
    origin_shift: tuple[Fraction, ...],
) -> _ExactSymmetryOperation:
    """Conjugate a local-setting Seitz operation through x_o = x_l U + q."""
    column_transform = _transpose(local_to_official)
    inverse_column_transform = _inverse(column_transform)
    rotation = _multiply(
        _multiply(column_transform, operation.rotation),
        inverse_column_transform,
    )
    translated = _matrix_vector(column_transform, operation.translation)
    shifted_origin = _matrix_vector(rotation, origin_shift)
    translation = _vector_add(
        translated,
        _vector_subtract(origin_shift, shifted_origin),
    )
    return _ExactSymmetryOperation(rotation=rotation, translation=translation)


def _operation_json(
    key: tuple[tuple[Fraction, ...], tuple[Fraction, ...]],
) -> dict[str, Any]:
    flat_rotation, translation = key
    rotation = [
        [_fraction_json(flat_rotation[3 * row + column]) for column in range(3)]
        for row in range(3)
    ]
    return {
        "rotation": rotation,
        "translation_modulo_one": [_fraction_json(value) for value in translation],
    }


def _compare_space_group_operations(
    local_block: Mapping[str, Any],
    official_block: Mapping[str, Any],
    local_to_official: tuple[tuple[Fraction, ...], ...],
    origin_shift: tuple[Fraction, ...],
) -> tuple[bool | None, tuple[_ExactSymmetryOperation, ...], dict[str, Any]]:
    """Compare declared Seitz sets after the exact candidate setting change."""
    try:
        local_operations = _space_group_operations(local_block)
        official_operations = _space_group_operations(official_block)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        return (
            None,
            (),
            {
                "status": "inconclusive",
                "error": f"{type(exc).__module__}.{type(exc).__qualname__}: {exc}",
            },
        )
    transformed_keys = {
        _operation_key(
            _transform_operation_to_official(
                operation,
                local_to_official,
                origin_shift,
            )
        )
        for operation in local_operations
    }
    official_keys = {_operation_key(operation) for operation in official_operations}
    missing = sorted(official_keys - transformed_keys)
    extra = sorted(transformed_keys - official_keys)
    equivalent = not missing and not extra
    return (
        equivalent,
        official_operations,
        {
            "status": "equivalent" if equivalent else "different",
            "local_declared_operation_count": len(local_operations),
            "official_declared_operation_count": len(official_operations),
            "local_unique_operation_count": len(
                {_operation_key(op) for op in local_operations}
            ),
            "official_unique_operation_count": len(official_keys),
            "missing_operation_count": len(missing),
            "extra_operation_count": len(extra),
            "missing_operation_samples": [_operation_json(key) for key in missing[:8]],
            "extra_operation_samples": [_operation_json(key) for key in extra[:8]],
            "convention": "R_o=A R_l A^-1; t_o=A t_l+q-R_o q, A=U^T",
        },
    )


def _perfect_matching(
    candidates: list[list[tuple[float, int]]],
) -> list[int] | None:
    """Return a deterministic perfect bipartite matching, if one exists."""
    for values in candidates:
        values.sort(key=lambda item: (item[0], item[1]))
    reference_to_local: dict[int, int] = {}

    def assign(local_index: int, visited: set[int]) -> bool:
        for _, reference_index in candidates[local_index]:
            if reference_index in visited:
                continue
            visited.add(reference_index)
            previous = reference_to_local.get(reference_index)
            if previous is None or assign(previous, visited):
                reference_to_local[reference_index] = local_index
                return True
        return False

    order = sorted(range(len(candidates)), key=lambda index: (len(candidates[index]), index))
    if any(not assign(local_index, set()) for local_index in order):
        return None
    local_to_reference = {
        local_index: reference_index
        for reference_index, local_index in reference_to_local.items()
    }
    return [local_to_reference[index] for index in range(len(candidates))]


def _site_match_diagnostics(
    local: tuple[_AtomSiteRow, ...],
    official: tuple[_AtomSiteRow, ...],
    transform: np.ndarray,
    origin_shift: np.ndarray,
    official_operations: tuple[_ExactSymmetryOperation, ...],
) -> dict[str, Any]:
    local_counts = Counter(row.species for row in local)
    official_counts = Counter(row.species for row in official)
    count_equal = len(local) == len(official)
    species_counts_equal = local_counts == official_counts
    coordinate_slack = 128.0 * np.finfo(float).eps * (
        1.0 + np.sum(np.abs(transform), axis=0)
    )
    occupancy_slack = 128.0 * np.finfo(float).eps
    coordinate_assignment: list[int] | None = None
    joint_assignment: list[int] | None = None
    maximum_by_component = np.zeros(3, dtype=float)
    maximum_bound_by_component = np.zeros(3, dtype=float)
    maximum_occupancy = 0.0
    maximum_occupancy_bound = 0.0
    coordinate_bound_ambiguous = False
    edge_details: dict[tuple[int, int], dict[str, Any]] = {}

    if count_equal and species_counts_equal:
        coordinate_candidates: list[list[tuple[float, int]]] = [[] for _ in local]
        joint_candidates: list[list[tuple[float, int]]] = [[] for _ in local]
        for local_index, local_row in enumerate(local):
            transformed = (
                np.asarray(local_row.coordinates, dtype=float) @ transform
                + origin_shift
            ) % 1.0
            local_half_bound = (
                0.5
                * np.asarray(local_row.coordinate_quantums, dtype=float)
                @ np.abs(transform)
            )
            for official_index, official_row in enumerate(official):
                if local_row.species != official_row.species:
                    continue
                best: dict[str, Any] | None = None
                for operation_index, operation in enumerate(official_operations):
                    rotation = _float_matrix(operation.rotation)
                    translation = np.asarray(
                        [float(value) for value in operation.translation], dtype=float
                    )
                    target = (
                        rotation @ np.asarray(official_row.coordinates, dtype=float)
                        + translation
                    ) % 1.0
                    official_half_bound = (
                        np.abs(rotation)
                        @ (0.5 * np.asarray(official_row.coordinate_quantums, dtype=float))
                    )
                    coordinate_bound = local_half_bound + official_half_bound
                    if np.any(coordinate_bound >= 0.5):
                        coordinate_bound_ambiguous = True
                    delta = target - transformed
                    delta -= np.round(delta)
                    absolute = np.abs(delta)
                    ratio_components = np.divide(
                        absolute,
                        coordinate_bound,
                        out=np.full(3, np.inf),
                        where=coordinate_bound > 0.0,
                    )
                    ratio_components[(coordinate_bound == 0.0) & (absolute == 0.0)] = 0.0
                    ratio = float(np.max(ratio_components))
                    if not np.all(absolute <= coordinate_bound + coordinate_slack):
                        continue
                    candidate = {
                        "ratio": ratio,
                        "operation_index": operation_index,
                        "absolute": absolute,
                        "coordinate_bound": coordinate_bound,
                    }
                    if best is None or (ratio, operation_index) < (
                        best["ratio"],
                        best["operation_index"],
                    ):
                        best = candidate
                if best is None:
                    continue
                edge_details[(local_index, official_index)] = best
                coordinate_candidates[local_index].append(
                    (best["ratio"], official_index)
                )
                occupancy_difference = abs(
                    local_row.occupancy - official_row.occupancy
                )
                occupancy_bound = 0.5 * (
                    local_row.occupancy_quantum + official_row.occupancy_quantum
                )
                if occupancy_difference <= occupancy_bound + occupancy_slack:
                    joint_candidates[local_index].append(
                        (best["ratio"], official_index)
                    )
        coordinate_assignment = _perfect_matching(coordinate_candidates)
        joint_assignment = _perfect_matching(joint_candidates)

        assignment = joint_assignment or coordinate_assignment
        if assignment is not None:
            for local_index, official_index in enumerate(assignment):
                edge = edge_details[(local_index, official_index)]
                maximum_by_component = np.maximum(
                    maximum_by_component, edge["absolute"]
                )
                maximum_bound_by_component = np.maximum(
                    maximum_bound_by_component, edge["coordinate_bound"]
                )
                maximum_occupancy = max(
                    maximum_occupancy,
                    abs(
                        local[local_index].occupancy
                        - official[official_index].occupancy
                    ),
                )
                maximum_occupancy_bound = max(
                    maximum_occupancy_bound,
                    0.5
                    * (
                        local[local_index].occupancy_quantum
                        + official[official_index].occupancy_quantum
                    ),
                )

    def composition(rows: tuple[_AtomSiteRow, ...]) -> dict[str, float]:
        totals: defaultdict[str, float] = defaultdict(float)
        for row in rows:
            totals[row.species] += row.occupancy
        return dict(sorted(totals.items()))

    local_composition = composition(local)
    official_composition = composition(official)

    return {
        "local_atom_count": len(local),
        "official_atom_count": len(official),
        "count_semantics": "raw atom-site loop rows before symmetry expansion",
        "atom_count_equal": count_equal,
        "composition_equal": local_composition == official_composition,
        "local_composition": local_composition,
        "official_composition": official_composition,
        "local_species_counts": dict(sorted(local_counts.items())),
        "official_species_counts": dict(sorted(official_counts.items())),
        "species_counts_equal": species_counts_equal,
        "coordinates_equal_at_output_precision": coordinate_assignment is not None,
        "occupancies_equal_at_output_precision": joint_assignment is not None,
        "coordinate_max_periodic_difference_by_component": maximum_by_component.tolist(),
        "coordinate_max_bound_by_component": maximum_bound_by_component.tolist(),
        "coordinate_bound_ambiguous_on_periodic_cell": coordinate_bound_ambiguous,
        "occupancy_max_abs_difference": maximum_occupancy,
        "occupancy_max_bound": maximum_occupancy_bound,
        "coordinate_arithmetic_slack_by_component": coordinate_slack.tolist(),
        "occupancy_arithmetic_slack": occupancy_slack,
    }


def compare_cif_alternate_settings(
    local_path: Path,
    official_path: Path,
    *,
    direct_details: dict[str, Any] | None = None,
) -> CifSettingComparison:
    """Compare candidate CIFs through their exact local-to-official setting map.

    Printed cell, coordinate, and occupancy values are interpreted as rounded
    values.  Their half-last-digit intervals are propagated through the exact
    integer basis transform instead of being replaced by a global tolerance.
    """
    local_path = Path(local_path)
    official_path = Path(official_path)
    try:
        local_header = parse_candidate_header(local_path)
        official_header = parse_candidate_header(official_path)
    except AmbiguousCandidateHeaderError as exc:
        return CifSettingComparison(
            equivalent=False,
            conclusive=False,
            details={
                "transform": {
                    "unique": False,
                    "error": f"{type(exc).__module__}.{type(exc).__qualname__}: {exc}",
                },
                "issues": [],
                "diagnostics": [],
                "inconclusive_reasons": [
                    "setting-transform-not-uniquely-declared"
                ],
            },
        )
    official_inverse = _inverse(official_header.basis)
    exact_transform = _multiply(local_header.basis, official_inverse)
    exact_origin_shift = _vector_multiply(
        tuple(
            local_header.origin[index] - official_header.origin[index]
            for index in range(3)
        ),
        official_inverse,
    )
    determinant = _determinant(exact_transform)
    integer_transform = all(
        value.denominator == 1 for row in exact_transform for value in row
    )
    unimodular = integer_transform and abs(determinant) == 1
    transform_details: dict[str, Any] = {
        "local_to_official_matrix": _fraction_matrix_json(exact_transform),
        "origin_shift_in_official_cell": [
            _fraction_json(value) for value in exact_origin_shift
        ],
        "determinant": _fraction_json(determinant),
        "integer": integer_transform,
        "unimodular": unimodular,
        "unique": True,
        "coordinate_convention": "x_official = x_local @ U + q (mod 1)",
    }
    if not unimodular:
        return CifSettingComparison(
            equivalent=False,
            conclusive=True,
            details={
                "transform": transform_details,
                "issues": ["non-unimodular-basis-transform"],
                "diagnostics": [],
            },
        )

    local_block = _cif_block(local_path)
    official_block = _cif_block(official_path)
    local_precision = _cif_precision(local_block)
    official_precision = _cif_precision(official_block)
    local_rows = _atom_site_rows(local_block)
    official_rows = _atom_site_rows(official_block)
    transform = _float_matrix(exact_transform)
    origin_shift = np.asarray([float(value) for value in exact_origin_shift])

    local_gram_intervals = _gram_intervals(
        local_precision.cell_values,
        local_precision.cell_quantums,
    )
    official_gram_intervals = _gram_intervals(
        official_precision.cell_values,
        official_precision.cell_quantums,
    )
    transformed_official_intervals = _transform_gram_intervals(
        official_gram_intervals,
        transform,
    )
    metric_failures: list[dict[str, Any]] = []
    for row in range(3):
        for column in range(row, 3):
            overlap, gap = _intervals_overlap(
                local_gram_intervals[row][column],
                transformed_official_intervals[row][column],
            )
            if not overlap:
                metric_failures.append(
                    {
                        "component": f"{_GRAM_LABELS[row]}{_GRAM_LABELS[column]}",
                        "local_interval": list(local_gram_intervals[row][column]),
                        "transformed_official_interval": list(
                            transformed_official_intervals[row][column]
                        ),
                        "gap": gap,
                    }
                )
    local_gram = _gram_matrix(local_precision.cell_values)
    official_gram = _gram_matrix(official_precision.cell_values)
    transformed_official_gram = transform @ official_gram @ transform.T
    metric_witness = (
        {
            "found": False,
            "method": "component-interval-separation-proves-mismatch",
        }
        if metric_failures
        else _joint_metric_witness(local_precision, official_precision, transform)
    )

    operation_equivalent, official_operations, operation_details = (
        _compare_space_group_operations(
            local_block,
            official_block,
            exact_transform,
            exact_origin_shift,
        )
    )
    if not official_operations:
        official_operations = (
            _ExactSymmetryOperation(
                rotation=(
                    (Fraction(1), Fraction(0), Fraction(0)),
                    (Fraction(0), Fraction(1), Fraction(0)),
                    (Fraction(0), Fraction(0), Fraction(1)),
                ),
                translation=(Fraction(0), Fraction(0), Fraction(0)),
            ),
        )
    sites = _site_match_diagnostics(
        local_rows,
        official_rows,
        transform,
        origin_shift,
        official_operations,
    )
    coordinate_bound_ambiguous = bool(
        sites["coordinate_bound_ambiguous_on_periodic_cell"]
    )
    local_magnetic_tags = _magnetic_moment_tags(local_block)
    official_magnetic_tags = _magnetic_moment_tags(official_block)

    direct_error: str | None = None
    if direct_details is None:
        try:
            direct_details = compare_cif(
                local_path,
                official_path,
                ignore_atom_order=True,
            ).details
        except (AttributeError, TypeError, ValueError) as exc:
            direct_error = f"{type(exc).__module__}.{type(exc).__qualname__}: {exc}"
            direct_details = {}
    local_declared = _declared_space_group(local_block)
    official_declared = _declared_space_group(official_block)
    declared_compatible, declared_status, declared_diagnostics = (
        _compare_declared_space_groups(local_declared, official_declared)
    )
    inferred_equal = direct_details.get("inferred_space_group_equal")

    issues: list[str] = []
    if metric_failures:
        issues.append("lattice-metric-outside-output-precision")
    if not sites["atom_count_equal"]:
        issues.append("atom-count-mismatch")
    if not sites["species_counts_equal"]:
        issues.append("species-count-mismatch")
    if not sites["coordinates_equal_at_output_precision"]:
        issues.append("coordinate-mismatch-outside-output-precision")
    elif not sites["occupancies_equal_at_output_precision"]:
        issues.append("occupancy-mismatch-outside-output-precision")
    if not declared_compatible:
        issues.append("declared-space-group-it-number-mismatch")
    if operation_equivalent is False:
        issues.append("declared-space-group-operation-mismatch")

    diagnostics = list(declared_diagnostics)
    if inferred_equal is False:
        diagnostics.append("spglib-inferred-space-group-differs")
    if direct_error is not None:
        diagnostics.append("direct-comparator-diagnostic-unavailable")
    if "declared space group does not match spglib inference" in direct_details.get(
        "issues", []
    ):
        diagnostics.append("declared-space-group-differs-from-spglib-inference")
    inconclusive_reasons: list[str] = []
    if not metric_failures and not metric_witness["found"]:
        inconclusive_reasons.append(
            "lattice-output-precision-joint-feasibility-not-proven"
        )
    if coordinate_bound_ambiguous:
        inconclusive_reasons.append("coordinate-output-precision-spans-periodic-cell")
    if operation_equivalent is None:
        inconclusive_reasons.append("declared-space-group-operation-equivalence-not-proven")
    if local_magnetic_tags or official_magnetic_tags:
        inconclusive_reasons.append("magnetic-moment-setting-transform-not-implemented")

    def quantum_ranges(
        rows: tuple[_AtomSiteRow, ...],
    ) -> tuple[list[list[float]], list[float]]:
        coordinate_ranges = [
            [
                min(row.coordinate_quantums[axis] for row in rows),
                max(row.coordinate_quantums[axis] for row in rows),
            ]
            for axis in range(3)
        ]
        occupancy_range = [
            min(row.occupancy_quantum for row in rows),
            max(row.occupancy_quantum for row in rows),
        ]
        return coordinate_ranges, occupancy_range

    local_coordinate_ranges, local_occupancy_range = quantum_ranges(local_rows)
    official_coordinate_ranges, official_occupancy_range = quantum_ranges(official_rows)

    def interval_payload(
        matrix: tuple[tuple[tuple[float, float], ...], ...],
    ) -> list[list[list[float]]]:
        return [[list(interval) for interval in row] for row in matrix]

    details = {
        "transform": transform_details,
        "precision": {
            "local_cell_quantums": list(local_precision.cell_quantums),
            "official_cell_quantums": list(official_precision.cell_quantums),
            "local_coordinate_quantums": list(
                local_precision.coordinate_quantums
            ),
            "official_coordinate_quantums": list(
                official_precision.coordinate_quantums
            ),
            "local_coordinate_quantum_ranges_by_component": local_coordinate_ranges,
            "official_coordinate_quantum_ranges_by_component": official_coordinate_ranges,
            "coordinate_bound_by_official_component": sites[
                "coordinate_max_bound_by_component"
            ],
            "coordinate_bound_ambiguous_on_periodic_cell": coordinate_bound_ambiguous,
            "local_occupancy_quantum": local_precision.occupancy_quantum,
            "official_occupancy_quantum": official_precision.occupancy_quantum,
            "local_occupancy_quantum_range": local_occupancy_range,
            "official_occupancy_quantum_range": official_occupancy_range,
            "occupancy_bound": sites["occupancy_max_bound"],
            "atom_site_precision_semantics": "per raw CIF atom-site row",
        },
        "lattice": {
            "metric_equal_at_output_precision": (
                not metric_failures and bool(metric_witness["found"])
            ),
            "nominal_max_abs_gram_difference": float(
                np.max(np.abs(local_gram - transformed_official_gram))
            ),
            "local_gram_intervals": interval_payload(local_gram_intervals),
            "transformed_official_gram_intervals": interval_payload(
                transformed_official_intervals
            ),
            "failed_components": metric_failures,
            "joint_feasibility_witness": metric_witness,
        },
        "sites": sites,
        "space_group": {
            "declared_compatible": declared_compatible,
            "declared_comparison_status": declared_status,
            "inferred_equal": inferred_equal,
            "inference_is_diagnostic_only": True,
            "local_declared": local_declared,
            "official_declared": official_declared,
            "local_inferred": direct_details.get("space_group", {}).get(
                "local_inferred"
            ),
            "official_inferred": direct_details.get("space_group", {}).get(
                "reference_inferred"
            ),
            "operation_equivalent_under_setting_transform": operation_equivalent,
            "operation_comparison": operation_details,
            "direct_comparator_error": direct_error,
            "local_magnetic_moment_tags": local_magnetic_tags,
            "official_magnetic_moment_tags": official_magnetic_tags,
        },
        "issues": issues,
        "diagnostics": diagnostics,
        "inconclusive_reasons": inconclusive_reasons,
    }
    conclusive = not inconclusive_reasons
    return CifSettingComparison(
        equivalent=not issues and not inconclusive_reasons,
        conclusive=conclusive,
        details=details,
    )


def _equivalent_structures(local: Path, official: Path) -> bool:
    """Compatibility wrapper for older saved-output audit callers."""
    return compare_cif_alternate_settings(local, official).equivalent


def _index_candidates(root: Path) -> tuple[dict[str, Path], list[tuple[str, list[Path]]]]:
    grouped: dict[str, list[Path]] = defaultdict(list)
    for cif in sorted(root.rglob("subgroup.cif")):
        grouped[parse_candidate_header(cif).identity].append(cif.parent)
    unique: dict[str, Path] = {}
    for key, values in grouped.items():
        if len(values) == 1:
            unique[key] = values[0]
            continue
        # A saved reference set can contain a copied/misnamed directory. Keep
        # the one whose directory label agrees with its internal IR/OPD header,
        # and report the remaining aliases as a reference conflict.
        ir, body = (part.strip() for part in key.split("|", 1))
        opd = body.split(None, 1)[0]
        matching = [path for path in values if path.name.startswith(f"{ir} {opd}")]
        if len(matching) == 1:
            unique[key] = matching[0]
    duplicates = [(key, values) for key, values in grouped.items() if len(values) > 1]
    return unique, duplicates


def _add_issue(
    summary: AuditSummary,
    severity: str,
    code: str,
    case: str,
    candidate: str,
    message: str,
) -> None:
    summary.issues.append(AuditIssue(severity, code, case, candidate, message))


def _check_formats(summary: AuditSummary, case: str, candidate: str, folder: Path, parametric: bool) -> None:
    for name in CORE_LOCAL_FILES:
        path = folder / name
        if not path.is_file():
            _add_issue(summary, "error", "missing-local-format", case, candidate, f"missing {name}")
            continue
        summary.format_files_checked += 1
        if path.stat().st_size == 0:
            _add_issue(summary, "error", "empty-local-format", case, candidate, f"empty {name}")
    markers = {
        "data.isoviz": ("!isoversion", "!atomcoordlist"),
        "topas.str": ("space_group", "'{{{distorted parameters"),
        "Complete modes details.txt": ("Complete modes details", "Subgroup:"),
    }
    for name, required in markers.items():
        path = folder / name
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8-sig", errors="replace")
        missing = [marker for marker in required if marker not in text]
        if missing:
            _add_issue(
                summary,
                "error",
                "invalid-local-format",
                case,
                candidate,
                f"{name} missing markers: {', '.join(missing)}",
            )
    if not parametric:
        return
    details = folder / "Complete modes details.txt"
    if details.is_file():
        text = details.read_text(encoding="utf-8-sig", errors="replace")
        if _count_local_displacive_modes(text) == 0:
            _add_issue(
                summary,
                "error",
                "empty-parametric-modes",
                case,
                candidate,
                "Complete modes details.txt contains no displacive parameter-k modes",
            )
    isoviz = folder / "data.isoviz"
    if isoviz.is_file():
        text = isoviz.read_text(encoding="utf-8-sig", errors="replace")
        marker = text.lower().find("!displacivemodelist")
        tail = text[marker:] if marker >= 0 else ""
        if marker < 0 or re.search(r"(?m)^\s*\d+\s+\d+\s+[-+0-9.]", tail) is None:
            _add_issue(
                summary,
                "error",
                "empty-parametric-isoviz-modes",
                case,
                candidate,
                "data.isoviz contains no displacive parameter-k mode rows",
            )
    topas = folder / "topas.str"
    if topas.is_file():
        text = topas.read_text(encoding="utf-8-sig", errors="replace")
        if re.search(r"(?m)^\s*prm\s+!a\d+\b", text) is None:
            _add_issue(
                summary,
                "error",
                "empty-parametric-topas-modes",
                case,
                candidate,
                "topas.str contains no active displacive mode parameters",
            )


def _count_local_displacive_modes(text: str) -> int:
    """Count actual mode definitions, excluding prose/table headings."""
    return len(
        re.findall(
            r"(?m)^Mode\s+(?!definitions\b|vectors\b|amplitudes\b)\S+",
            text,
        )
    )


def _count_official_displacive_modes(folder: Path) -> int | None:
    html_files = sorted(folder.glob("*complete modes details*.html"))
    if not html_files:
        return None
    text = html_files[0].read_text(encoding="utf-8-sig", errors="replace")
    return len(re.findall(r":dsp\].*?\bnormfactor\s*=", text, flags=re.IGNORECASE))


def _compare_mode_counts(
    summary: AuditSummary,
    case: str,
    candidate: str,
    official_folder: Path,
    local_folder: Path,
) -> None:
    official_count = _count_official_displacive_modes(official_folder)
    details = local_folder / "Complete modes details.txt"
    if official_count is None or not details.is_file():
        return
    local_text = details.read_text(encoding="utf-8-sig", errors="replace")
    local_count = _count_local_displacive_modes(local_text)
    summary.mode_files_compared += 1
    if local_count == official_count:
        summary.mode_count_matches += 1
        return
    _add_issue(
        summary,
        "error",
        "displacive-mode-count-mismatch",
        case,
        candidate,
        f"official={official_count}; local={local_count}",
    )


def audit_output_root(output_root: Path) -> AuditSummary:
    output_root = output_root.expanduser().resolve()
    if not output_root.is_dir():
        raise FileNotFoundError(f"output comparison root not found: {output_root}")
    summary = AuditSummary(output_root=str(output_root))
    for case_root in sorted(path for path in output_root.iterdir() if path.is_dir()):
        # Method 3 rows are affine embeddings and can legitimately have no
        # single-IR header.  They are audited by audit_method3_downloads.py and
        # compare_method3_official_local.py, which bind the saved result page to
        # the manifest and compare exact affine subgroup operators.  Reusing the
        # Method 1/2 IR-path parser here produced false "CIF header has no IR"
        # failures for scientifically valid coupled/route-less embeddings.
        for method in ("Method1", "Method2"):
            official_root = case_root / OFFICIAL_DIR / method
            local_root = case_root / LOCAL_DIR / method
            if not official_root.exists() and not local_root.exists():
                continue
            case = f"{case_root.name}/{method}"
            if not official_root.is_dir() or not local_root.is_dir():
                _add_issue(summary, "error", "missing-side", case, "", "official or local method directory is missing")
                continue
            try:
                official, official_duplicates = _index_candidates(official_root)
                local, local_duplicates = _index_candidates(local_root)
            except (OSError, ValueError) as exc:
                _add_issue(summary, "error", "inventory-parse", case, "", str(exc))
                continue
            summary.groups[case] = {"official": len(official), "local": len(local)}
            official_conflict_names = {
                path.name for _key, paths in official_duplicates for path in paths
            }
            for key, paths in official_duplicates:
                _add_issue(
                    summary,
                    "warning",
                    "reference-conflict",
                    case,
                    key,
                    "duplicate internal identity: " + "; ".join(map(str, paths)),
                )
            for key, paths in local_duplicates:
                _add_issue(summary, "error", "duplicate-identity", case, key, "; ".join(map(str, paths)))
            for key in sorted(official.keys() - local.keys()):
                _add_issue(summary, "error", "missing-local-candidate", case, key, official[key].name)
            for key in sorted(local.keys() - official.keys()):
                if local[key].name in official_conflict_names:
                    _add_issue(
                        summary,
                        "warning",
                        "reference-conflict",
                        case,
                        key,
                        f"official directory {local[key].name!r} contains another candidate's CIF",
                    )
                else:
                    _add_issue(summary, "error", "extra-local-candidate", case, key, local[key].name)
            for key in sorted(official.keys() & local.keys()):
                summary.paired_candidates += 1
                official_folder = official[key]
                local_folder = local[key]
                official_cif = official_folder / "subgroup.cif"
                local_cif = local_folder / "subgroup.cif"
                official_header = parse_candidate_header(official_cif)
                local_header = parse_candidate_header(local_cif)
                try:
                    same_lattice = bases_generate_same_lattice(local_header.basis, official_header.basis)
                except ValueError as exc:
                    _add_issue(summary, "error", "invalid-basis", case, key, str(exc))
                    same_lattice = False
                if same_lattice:
                    summary.basis_equivalent_pairs += 1
                    if local_header.basis == official_header.basis:
                        summary.basis_exact_pairs += 1
                    else:
                        # All subsequent semantic checks use the official file
                        # as the target setting.  This is an explicit
                        # local->official unimodular normalization, not a
                        # tolerance-based waiver of the representation.
                        summary.bases_normalized_to_official += 1
                else:
                    _add_issue(
                        summary,
                        "error",
                        "non-equivalent-basis",
                        case,
                        key,
                        f"official={official_folder.name}; local={local_folder.name}",
                    )
                try:
                    direct = compare_cif(local_cif, official_cif, ignore_atom_order=True)
                    if direct.structure_equal:
                        summary.direct_cif_matches += 1
                    elif same_lattice and compare_cif_alternate_settings(
                        local_cif,
                        official_cif,
                        direct_details=direct.details,
                    ).equivalent:
                        summary.equivalent_setting_matches += 1
                    else:
                        _add_issue(
                            summary,
                            "error",
                            "cif-semantic-mismatch",
                            case,
                            key,
                            "; ".join(direct.details["issues"]),
                        )
                except (OSError, TypeError, ValueError) as exc:
                    _add_issue(summary, "error", "cif-compare-error", case, key, str(exc))
                _check_formats(summary, case, key, local_folder, local_header.parametric_k)
                _compare_mode_counts(summary, case, key, official_folder, local_folder)
    return summary


def _payload(summary: AuditSummary) -> dict[str, Any]:
    result = asdict(summary)
    result.update(failures=summary.failures, warnings=summary.warnings, passed=summary.failures == 0)
    return result


def _print_summary(summary: AuditSummary) -> None:
    print("Method output scientific audit")
    for case, counts in summary.groups.items():
        print(f"  {case}: official={counts['official']} local={counts['local']}")
    print(f"Paired candidates: {summary.paired_candidates}")
    print(f"Direct CIF semantic matches: {summary.direct_cif_matches}")
    print(f"Equivalent-setting CIF matches: {summary.equivalent_setting_matches}")
    print(f"Equivalent basis pairs: {summary.basis_equivalent_pairs}")
    print(f"Exact official basis representations: {summary.basis_exact_pairs}")
    print(f"Bases normalized to official representation: {summary.bases_normalized_to_official}")
    print(f"Local core files checked: {summary.format_files_checked}")
    print(
        "Displacive mode-count matches: "
        f"{summary.mode_count_matches}/{summary.mode_files_compared}"
    )
    for code, count in sorted(Counter(issue.code for issue in summary.issues).items()):
        print(f"  {code}: {count}")
    for issue in summary.issues[:40]:
        candidate = f" [{issue.candidate}]" if issue.candidate else ""
        print(
            f"{issue.severity.upper()} {issue.case} {issue.code}{candidate}: "
            f"{issue.message}"
        )
    if len(summary.issues) > 40:
        print(f"... {len(summary.issues) - 40} additional issues (use --json for full details)")
    print(f"Conclusion: {'PASS' if summary.failures == 0 else 'FAIL'}")


def run_live_method3_smoke() -> list[dict[str, Any]]:
    """Run the two fixed Method 3 cases through the real iso/WSL backend.

    The expected fields are a test manifest captured from the official site,
    never consulted by production code.  ZIPs stay in memory so the check does
    not leave downloads or extracted trees behind.
    """
    project_root = WORKSPACE / "ISODISTORT"
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from isocore.api import IsoDistort

    cases = (
        {
            "name": "EuAl4 Parent.cif",
            "target": 99,
            "basis": [[1, 0, 0], [0, 1, 0], [0, 0, 6]],
            "expected": {("LD1", "C1", "(0,0,1/6)", "0,0,0")},
        },
        {
            "name": "NdNiO2 own.cif",
            "target": 47,
            "basis": [[-3, 0, 0], [0, 0, 1], [0, 2, 0]],
            "expected": {
                ("Y1", "P1", "(1/3,1/2,0)", "0,0,0"),
                ("Y4", "P1", "(1/3,1/2,0)", "0,1/2,0"),
            },
        },
    )
    reports: list[dict[str, Any]] = []
    for case in cases:
        api = IsoDistort()
        api.load_structure(WORKSPACE / "experiment_data" / case["name"])
        items = api.search_method_3(
            distortion_types=["strain", "displacive"],
            space_group_type=case["target"],
            supercell_basis=case["basis"],
            direct_sublattice_centering="d",
            lattice_type="direct",
            generate_if_missing=True,
        )
        actual = {
            (
                item.subgroup.irrep_label,
                item.subgroup.opd_symbol,
                item.subgroup.official_fields()["k_active"],
                ",".join(_canonical_number(value) for value in item.subgroup.origin),
            )
            for item in items
        }
        if actual != case["expected"]:
            raise AssertionError(
                f"{case['name']} Method 3 candidates differ from the official manifest: "
                f"expected={sorted(case['expected'])}, actual={sorted(actual)}"
            )
        requested_basis = tuple(
            tuple(Fraction(value) for value in row) for row in case["basis"]
        )
        for item in items:
            subgroup = item.subgroup
            if subgroup.space_group_number != case["target"]:
                raise AssertionError(
                    f"{case['name']}: expected SG {case['target']}, got {subgroup.space_group_number}"
                )
            candidate_basis = tuple(
                tuple(Fraction(str(value)) for value in row)
                for row in subgroup.basis_vectors
            )
            if not bases_generate_same_lattice(candidate_basis, requested_basis):
                raise AssertionError(f"{case['name']}: candidate basis does not generate the requested lattice")

        body = api.export_subgroups_zip(
            formats=["cif"],
            subgroups=[item.subgroup for item in items],
            compute_missing_modes=False,
        )
        with zipfile.ZipFile(io.BytesIO(body)) as archive:
            cif_names = [name for name in archive.namelist() if name.endswith("/subgroup.cif")]
            if len(cif_names) != len(items):
                raise AssertionError(
                    f"{case['name']}: ZIP contains {len(cif_names)} CIFs for {len(items)} candidates"
                )
            for name in cif_names:
                cif_text = archive.read(name).decode("utf-8-sig")
                match = re.search(r"(?im)^_symmetry_Int_Tables_number\s+(\d+)\s*$", cif_text)
                if match is None or int(match.group(1)) != case["target"]:
                    raise AssertionError(f"{case['name']}: {name} has the wrong declared space group")
                parsed = CifParser.from_str(cif_text, occupancy_tolerance=100).parse_structures(
                    primitive=False
                )
                if len(parsed) != 1:
                    raise AssertionError(f"{case['name']}: {name} is not a single parseable structure")
        reports.append(
            {
                "case": case["name"],
                "target_space_group": case["target"],
                "candidates": len(items),
                "zip_cifs": len(cif_names),
                "identities": sorted(" | ".join(fields) for fields in actual),
            }
        )
    return reports


def _live_subgroup_identity(subgroup: Any) -> str:
    """Build the same basis-independent identity used for saved CIF headers."""
    line = subgroup.opd_line_body()
    line = re.sub(r",?\s*basis=\{.*?\}(?=,\s*origin=)", "", line)
    origin = tuple(
        Fraction(str(value)).limit_denominator(10_000)
        for value in subgroup.origin
    )
    line = re.sub(
        r"origin=\([^)]*\)",
        f"origin=({_canonical_origin(origin)})",
        line,
    )
    return f"{subgroup.irrep_label} | {' '.join(line.split())}"


def _live12_candidate_fingerprint(subgroup: Any) -> str:
    """Sign all subgroup fields that can steer the downstream mode calculation."""
    fields = (
        "index",
        "space_group_number",
        "space_group_symbol",
        "subgroup_index",
        "size",
        "is_maximal",
        "opd_symbol",
        "opd_vector",
        "basis_vectors",
        "origin",
        "k_point_label",
        "irrep_label",
        "k_parameters",
        "k_coordinates",
        "parent_sg",
        "opd_dir_raw",
        "basis_raw",
        "origin_raw",
        "k_active_raw",
    )
    payload = {
        "opd_line_body": subgroup.opd_line_body(),
        "fields": {name: getattr(subgroup, name, None) for name in fields},
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(encoded, usedforsecurity=False).hexdigest()


def _live12_now() -> str:
    return datetime.now(UTC).isoformat()


def _live12_runtime_versions() -> dict[str, Any]:
    """Describe the interpreter and scientific packages that affect results."""
    distributions: dict[str, str] = {}
    for name in RUNTIME_DISTRIBUTIONS:
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "<missing>"
        except Exception as exc:  # noqa: BLE001 - signatures must survive broken package metadata
            error_type = f"{type(exc).__module__}.{type(exc).__qualname__}"
            distributions[name] = f"<unavailable:{error_type}>"
    return {
        "python": {
            "implementation": sys.implementation.name,
            "version": list(sys.version_info[:5]),
            "cache_tag": sys.implementation.cache_tag,
            "build": sys.version,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "distributions": distributions,
    }


def _live12_configured_path(config_path: Path, value: object, fallback: str) -> Path:
    raw = value if isinstance(value, str) and value.strip() else fallback
    path = Path(raw)
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def _live12_backend_dependencies(package_root: Path) -> list[tuple[str, Path]]:
    """Return immutable ISO/smodes inputs; generated WSL ``i*.iso`` is excluded."""
    config_path = package_root / "config" / "settings.yaml"
    isobyu: dict[str, Any] = {}
    try:
        loaded = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict) and isinstance(loaded.get("isobyu"), dict):
            isobyu = loaded["isobyu"]
    except (OSError, UnicodeError, yaml.YAMLError):
        # The configuration file is itself content-signed below.  Deterministic
        # fallbacks preserve an explicit missing-file marker for a broken config.
        pass

    bin_dir = _live12_configured_path(config_path, isobyu.get("bin_dir"), "../isobyu")
    data_dir = _live12_configured_path(config_path, isobyu.get("data_dir"), "../isobyu")
    iso_name = isobyu.get("iso_bin")
    if not isinstance(iso_name, str) or not iso_name.strip():
        iso_name = "iso"
    smodes_name = isobyu.get("smodes_bin")
    if not isinstance(smodes_name, str) or not smodes_name.strip():
        smodes_name = "smodes"

    dependencies = [
        ("iso-binary", bin_dir / iso_name),
        ("smodes-binary", bin_dir / smodes_name),
    ]
    data_names = set(ISO_DATA_BASENAMES)
    try:
        data_names.update(path.name for path in data_dir.glob("data_*.txt"))
    except OSError:
        pass
    dependencies.extend(
        (f"iso-data:{name}", data_dir / name)
        for name in sorted(data_names, key=str.casefold)
    )
    return dependencies


def _live12_checkpoint_policy() -> dict[str, Any]:
    """Document why mutable generated subgroup caches are not batch-signed."""
    return {
        "scientific_inputs": "content_signed",
        "candidate_mode_counts": (
            "reuse_only_after_the_same_identity_and_full_candidate_fingerprint_are_rediscovered_"
            "by_current_live_enumeration"
        ),
        "generated_wsl_i_star_iso_cache": {
            "included_in_global_signature": False,
            "reason": (
                "generate_if_missing may create or update ~/.id/tmp/i*.iso during this audit; "
                "signing that mutable global directory would invalidate the audit's own checkpoint"
            ),
            "safety_guard": (
                "the live candidate set is enumerated on every invocation and a saved count is "
                "used only when identity plus basis/origin/k/OPD/subgroup fields match the current "
                "live object and the identity is also present in the official set"
            ),
        },
    }


def _live12_hash_path(digest: Any, path: Path, label: str) -> None:
    digest.update(label.encode("utf-8"))
    digest.update(b"\0")
    if path.is_file():
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
    else:
        digest.update(b"<missing>")
    digest.update(b"\0")


def _live12_signature_state(
    output_root: Path,
    *,
    workspace: Path | None = None,
    package_root: Path | None = None,
) -> dict[str, Any]:
    """Fingerprint every immutable scientific input needed for safe resume.

    This deliberately does not walk WSL's mutable ``~/.id/tmp/i*.iso`` cache.
    The accompanying policy explains the identity-level guard that makes
    interrupted Method 1/2 count reuse safe without self-invalidating a run.
    """
    workspace = (workspace or WORKSPACE).resolve()
    package_root = (package_root or workspace / "ISODISTORT").resolve()
    output_root = output_root.resolve()
    runtime = _live12_runtime_versions()
    policy = _live12_checkpoint_policy()
    digest = hashlib.sha256(usedforsecurity=False)
    digest.update(
        json.dumps(
            {"runtime_dependencies": runtime, "checkpoint_reuse_policy": policy},
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    digest.update(b"\0")

    inputs: list[tuple[str, Path]] = [
        ("validator", Path(__file__).resolve()),
        ("config", package_root / "config" / "settings.yaml"),
    ]
    inputs.extend(
        (f"isocore:{path.relative_to(package_root).as_posix()}", path)
        for path in sorted(
            (package_root / "isocore").rglob("*.py"),
            key=lambda item: item.as_posix().casefold(),
        )
    )
    for case_name in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        inputs.append((f"parent:{case_name}", workspace / "experiment_data" / case_name))

    backend_dependencies = _live12_backend_dependencies(package_root)
    inputs.extend(backend_dependencies)
    official_reference_count = 0
    for case_name in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        for method in ("Method1", "Method2"):
            official_root = output_root / case_name / OFFICIAL_DIR / method
            prefix = f"official:{case_name}:{method}"
            if not official_root.is_dir():
                inputs.append((f"{prefix}:<missing-directory>", official_root))
                continue
            references = [
                path
                for path in official_root.rglob("*")
                if path.is_file()
                and (
                    path.name.casefold() == "subgroup.cif"
                    or (
                        "complete modes details" in path.name.casefold()
                        and path.suffix.casefold() in {".html", ".htm"}
                    )
                )
            ]
            if not references:
                inputs.append((f"{prefix}:<no-reference-files>", official_root / ".missing"))
            for path in sorted(
                references,
                key=lambda item: item.relative_to(official_root).as_posix().casefold(),
            ):
                relative = path.relative_to(official_root).as_posix()
                inputs.append((f"{prefix}:{relative}", path))
                official_reference_count += 1

    for label, path in inputs:
        _live12_hash_path(digest, path, label)
    return {
        "signature": digest.hexdigest(),
        "runtime_dependencies": runtime,
        "checkpoint_reuse_policy": policy,
        "signature_inputs": {
            "config": str((package_root / "config" / "settings.yaml").resolve()),
            "python_source_file_count": 1 + sum(label.startswith("isocore:") for label, _ in inputs),
            "official_reference_file_count": official_reference_count,
            "backend_dependencies": [
                {
                    "label": label,
                    "path": str(path.resolve()),
                    "exists": path.is_file(),
                    "content_signed": True,
                }
                for label, path in backend_dependencies
            ],
        },
    }


def _live12_signature(
    output_root: Path,
    *,
    workspace: Path | None = None,
    package_root: Path | None = None,
) -> str:
    """Return the content/runtime signature used by the Method 1/2 checkpoint."""
    return str(
        _live12_signature_state(
            output_root,
            workspace=workspace,
            package_root=package_root,
        )["signature"]
    )


def _save_live12_checkpoint(path: Path, payload: dict[str, Any]) -> None:
    """Atomically persist audit progress so terminal interruption is harmless."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(path)


def _live12_report_payload(checkpoint: dict[str, Any]) -> dict[str, Any]:
    reports = list(checkpoint.get("reports", {}).values())
    return {
        "schema": LIVE12_CHECKPOINT_SCHEMA,
        "signature": checkpoint["signature"],
        "status": checkpoint["status"],
        "execution_complete": checkpoint["execution_complete"],
        "validation_complete": checkpoint.get("validation_complete", False),
        "passed": checkpoint.get("passed", False),
        "created_at": checkpoint["created_at"],
        "updated_at": checkpoint["updated_at"],
        "completed_at": checkpoint.get("completed_at"),
        "last_invocation_started_at": checkpoint["last_invocation_started_at"],
        "last_invocation_elapsed_seconds": checkpoint["last_invocation_elapsed_seconds"],
        "cumulative_elapsed_seconds": checkpoint["cumulative_elapsed_seconds"],
        "resumed_from_checkpoint": checkpoint["resumed_from_checkpoint"],
        "checkpoint_reuse_policy": checkpoint["signature_state"]["checkpoint_reuse_policy"],
        "runtime_dependencies": checkpoint["signature_state"]["runtime_dependencies"],
        "signature_inputs": checkpoint["signature_state"]["signature_inputs"],
        "progress": checkpoint["progress"],
        "last_error": checkpoint.get("last_error"),
        "reports": reports,
    }


def _update_live12_timing(
    checkpoint: dict[str, Any],
    invocation_started: float,
    previous_elapsed: float,
) -> None:
    elapsed = max(0.0, time.perf_counter() - invocation_started)
    checkpoint["updated_at"] = _live12_now()
    checkpoint["last_invocation_elapsed_seconds"] = round(elapsed, 6)
    checkpoint["cumulative_elapsed_seconds"] = round(previous_elapsed + elapsed, 6)


def _save_live12_state(
    checkpoint_path: Path,
    report_path: Path | None,
    checkpoint: dict[str, Any],
) -> None:
    _save_live12_checkpoint(checkpoint_path, checkpoint)
    if report_path is not None:
        _save_live12_checkpoint(report_path, _live12_report_payload(checkpoint))


def run_live_method12_audit(
    output_root: Path,
    checkpoint_path: Path,
    report_path: Path | None = None,
    *,
    resume: bool = True,
) -> dict[str, Any]:
    """Recompute every saved Method 1/2 candidate and compare mode counts.

    The official folders are read only for candidate identities and expected
    displacive-mode counts.  All local results stay in memory, which separates
    stale saved downloads from defects in the current computation engine.
    Successful candidate counts are atomically checkpointed after every item.
    Reuse requires an exact signature over source, configuration, interpreter,
    scientific package versions, parent/reference files, ISO/smodes binaries,
    and bundled data tables.  Mutable generated WSL ``i*.iso`` files are not
    globally signed; instead, identities are always enumerated afresh and only
    a current identity with the same complete subgroup fingerprint can reuse
    its saved mode count.
    """
    project_root = WORKSPACE / "ISODISTORT"
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from isocore.api import IsoDistort

    cases = (
        ("EuAl4 Parent.cif", "LD", ["1/6"]),
        ("NdNiO2 own.cif", "Y", ["1/3"]),
    )
    invocation_started = time.perf_counter()
    invocation_started_at = _live12_now()
    signature_state = _live12_signature_state(output_root)
    signature = str(signature_state["signature"])
    checkpoint: dict[str, Any] = {
        "schema": LIVE12_CHECKPOINT_SCHEMA,
        "signature": signature,
        "signature_state": signature_state,
        "status": "incomplete",
        "execution_complete": False,
        "validation_complete": False,
        "passed": False,
        "created_at": invocation_started_at,
        "updated_at": invocation_started_at,
        "last_invocation_started_at": invocation_started_at,
        "last_invocation_elapsed_seconds": 0.0,
        "cumulative_elapsed_seconds": 0.0,
        "resumed_from_checkpoint": False,
        "counts": {},
        "reports": {},
        "progress": {
            "completed_groups": 0,
            "total_groups": LIVE12_GROUP_COUNT,
            "current_group": None,
            "current_candidate": 0,
            "current_group_candidates": 0,
            "checkpointed_candidate_counts": 0,
        },
    }
    if resume and checkpoint_path.is_file():
        try:
            saved = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            saved = None
        if (
            isinstance(saved, dict)
            and saved.get("schema") == LIVE12_CHECKPOINT_SCHEMA
            and saved.get("signature") == signature
            and isinstance(saved.get("counts"), dict)
            and isinstance(saved.get("reports"), dict)
        ):
            checkpoint = saved
            checkpoint["signature_state"] = signature_state
            checkpoint["resumed_from_checkpoint"] = any(
                isinstance(values, dict) and bool(values)
                for values in checkpoint["counts"].values()
            )
            print(f"Resuming live Method 1/2 audit from {checkpoint_path}")
        else:
            print("Ignoring incompatible live Method 1/2 checkpoint")
    checkpoint.setdefault("created_at", invocation_started_at)
    previous_elapsed_raw = checkpoint.get("cumulative_elapsed_seconds", 0.0)
    previous_elapsed = (
        float(previous_elapsed_raw)
        if isinstance(previous_elapsed_raw, (int, float)) and not isinstance(previous_elapsed_raw, bool)
        else 0.0
    )
    checkpoint.update(
        {
            "signature_state": signature_state,
            "status": "incomplete",
            "execution_complete": False,
            "validation_complete": False,
            "passed": False,
            "last_invocation_started_at": invocation_started_at,
            "last_invocation_elapsed_seconds": 0.0,
            "reports": {},
            "progress": {
                "completed_groups": 0,
                "total_groups": LIVE12_GROUP_COUNT,
                "current_group": None,
                "current_candidate": 0,
                "current_group_candidates": 0,
                "checkpointed_candidate_counts": sum(
                    len(values)
                    for values in checkpoint["counts"].values()
                    if isinstance(values, dict)
                ),
            },
        }
    )
    checkpoint.pop("completed_at", None)
    checkpoint.pop("last_error", None)
    _save_live12_state(checkpoint_path, report_path, checkpoint)

    try:
        for case_name, k_label, k_parameters in cases:
            for method in ("Method1", "Method2"):
                group_started = time.perf_counter()
                group_key = f"{case_name}|{method}"
                official_root = output_root / case_name / OFFICIAL_DIR / method
                official, duplicates = _index_candidates(official_root)
                expected = {
                    identity: count
                    for identity, folder in official.items()
                    if (count := _count_official_displacive_modes(folder)) is not None
                }

                api = IsoDistort()
                with contextlib.redirect_stdout(io.StringIO()):
                    api.load_structure(WORKSPACE / "experiment_data" / case_name)
                    if method == "Method1":
                        rows = api.search_method_1(
                            distortion_types=["strain", "displacive"],
                        )
                        candidates = [row.subgroup for row in rows]
                    else:
                        candidates = api.list_subgroups_at_kpoint(
                            k_label,
                            k_parameters=k_parameters,
                            generate_if_missing=True,
                        )

                live: dict[str, Any] = {}
                duplicate_live: set[str] = set()
                for subgroup in candidates:
                    identity = _live_subgroup_identity(subgroup)
                    if identity in live:
                        duplicate_live.add(identity)
                    live[identity] = subgroup

                missing = sorted(official.keys() - live.keys())
                extra = sorted(live.keys() - official.keys())
                mismatches: list[dict[str, Any]] = []
                matched = 0
                identities = sorted(expected.keys() & live.keys())
                saved_counts = checkpoint["counts"].get(group_key, {})
                counts = {
                    identity: record
                    for identity, record in saved_counts.items()
                    if identity in identities
                    and isinstance(record, dict)
                    and isinstance(record.get("count"), int)
                    and not isinstance(record.get("count"), bool)
                    and record.get("candidate_fingerprint")
                    == _live12_candidate_fingerprint(live[identity])
                } if isinstance(saved_counts, dict) else {}
                checkpoint["counts"][group_key] = counts
                computed_count = 0
                reused_count = 0
                checkpoint["progress"].update(
                    {
                        "current_group": group_key,
                        "current_candidate": 0,
                        "current_group_candidates": len(identities),
                    }
                )
                _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
                _save_live12_state(checkpoint_path, report_path, checkpoint)

                for position, identity in enumerate(identities, start=1):
                    subgroup = live[identity]
                    candidate_fingerprint = _live12_candidate_fingerprint(subgroup)
                    saved_record = counts.get(identity)
                    if (
                        isinstance(saved_record, dict)
                        and isinstance(saved_record.get("count"), int)
                        and not isinstance(saved_record.get("count"), bool)
                        and saved_record.get("candidate_fingerprint") == candidate_fingerprint
                    ):
                        actual_count = saved_record["count"]
                        reused_count += 1
                    else:
                        with contextlib.redirect_stdout(io.StringIO()):
                            api.search_method_2(
                                subgroup_idx=subgroup.index,
                                distortion_type=["displacive"],
                                number_of_independent_modulations=0,
                                candidates=candidates,
                            )
                        actual_count = len(api.mode_displacements)
                        counts[identity] = {
                            "count": actual_count,
                            "candidate_fingerprint": candidate_fingerprint,
                        }
                        computed_count += 1
                    checkpoint["progress"].update(
                        {
                            "current_candidate": position,
                            "checkpointed_candidate_counts": sum(
                                len(values)
                                for values in checkpoint["counts"].values()
                                if isinstance(values, dict)
                            ),
                        }
                    )
                    _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
                    _save_live12_state(checkpoint_path, report_path, checkpoint)
                    if position % 10 == 0 or position == len(identities):
                        print(
                            f"Live audit progress {case_name} {method}: "
                            f"{position}/{len(identities)}"
                        )
                    expected_count = expected[identity]
                    if actual_count == expected_count:
                        matched += 1
                    else:
                        mismatches.append(
                            {
                                "candidate": identity,
                                "official": expected_count,
                                "live": actual_count,
                            }
                        )

                duplicate_official = [
                    {
                        "candidate": identity,
                        "directories": [str(path) for path in paths],
                    }
                    for identity, paths in duplicates
                ]
                reference_complete = bool(official) and not duplicate_official and len(expected) == len(official)
                validation_passed = bool(
                    reference_complete
                    and not missing
                    and not extra
                    and not duplicate_live
                    and not mismatches
                )
                report = {
                    "case": case_name,
                    "method": method,
                    "status": (
                        "passed"
                        if validation_passed
                        else "failed"
                        if reference_complete
                        else "incomplete_reference"
                    ),
                    "execution_complete": True,
                    "validation_complete": reference_complete,
                    "validation_passed": validation_passed,
                    "official_candidates": len(official),
                    "official_candidates_with_mode_reference": len(expected),
                    "live_candidates": len(live),
                    "mode_files_compared": len(expected.keys() & live.keys()),
                    "mode_count_matches": matched,
                    "checkpoint_counts_reused": reused_count,
                    "mode_counts_computed": computed_count,
                    "elapsed_seconds": round(time.perf_counter() - group_started, 6),
                    "missing": missing,
                    "extra": extra,
                    "duplicate_official": duplicate_official,
                    "duplicate_live": sorted(duplicate_live),
                    "mismatches": mismatches,
                }
                checkpoint["reports"][group_key] = report
                checkpoint["progress"].update(
                    {
                        "completed_groups": len(checkpoint["reports"]),
                        "current_group": None,
                        "current_candidate": 0,
                        "current_group_candidates": 0,
                    }
                )
                _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
                _save_live12_state(checkpoint_path, report_path, checkpoint)

        reports = list(checkpoint["reports"].values())
        validation_complete = (
            len(reports) == LIVE12_GROUP_COUNT
            and all(report["validation_complete"] for report in reports)
        )
        passed = validation_complete and all(report["validation_passed"] for report in reports)
        checkpoint.update(
            {
                "status": (
                    "complete_passed"
                    if passed
                    else "complete_failed"
                    if validation_complete
                    else "complete_validation_incomplete"
                ),
                "execution_complete": True,
                "validation_complete": validation_complete,
                "passed": passed,
                "completed_at": _live12_now(),
            }
        )
        _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
        _save_live12_state(checkpoint_path, report_path, checkpoint)
        return _live12_report_payload(checkpoint)
    except BaseException as exc:
        checkpoint.update(
            {
                "status": "incomplete_error",
                "execution_complete": False,
                "validation_complete": False,
                "passed": False,
                "last_error": f"{type(exc).__name__}: {exc}",
            }
        )
        _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
        _save_live12_state(checkpoint_path, report_path, checkpoint)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--json", action="store_true", help="write the complete audit result to stdout as JSON")
    parser.add_argument(
        "--live-method3",
        action="store_true",
        help="also run the fixed EuAl4/NdNiO2 Method 3 cases through iso/WSL",
    )
    parser.add_argument(
        "--live-method12",
        action="store_true",
        help="recompute every saved Method 1/2 candidate and compare current mode counts",
    )
    parser.add_argument(
        "--live-method12-checkpoint",
        type=Path,
        default=WORKSPACE / "ISODISTORT" / "output" / "validation" / "live_method12_checkpoint.json",
        help="atomic resume/checkpoint JSON for --live-method12",
    )
    parser.add_argument(
        "--live-method12-report",
        type=Path,
        default=WORKSPACE / "ISODISTORT" / "output" / "validation" / "live_method12_report.json",
        help="atomic complete/incomplete report JSON for --live-method12",
    )
    parser.add_argument(
        "--live-method12-restart",
        action="store_true",
        help="ignore an otherwise compatible Method 1/2 checkpoint and recompute all mode counts",
    )
    args = parser.parse_args(argv)
    try:
        summary = audit_output_root(args.output_root)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    try:
        live_method3 = run_live_method3_smoke() if args.live_method3 else []
        live_method12 = (
            run_live_method12_audit(
                args.output_root.expanduser().resolve(),
                args.live_method12_checkpoint.expanduser().resolve(),
                args.live_method12_report.expanduser().resolve(),
                resume=not args.live_method12_restart,
            )
            if args.live_method12 else {}
        )
    except (AssertionError, OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: live validation failed: {exc}", file=sys.stderr)
        return 1
    if args.json:
        payload = _payload(summary)
        payload["live_method3"] = live_method3
        payload["live_method12"] = live_method12
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        _print_summary(summary)
        for report in live_method3:
            print(
                f"Live Method 3 PASS {report['case']}: "
                f"SG={report['target_space_group']} candidates={report['candidates']} "
                f"ZIP_CIFs={report['zip_cifs']}"
            )
        if live_method12:
            print(
                "Live Method 1/2 audit: "
                f"status={live_method12['status']} "
                f"execution_complete={live_method12['execution_complete']} "
                f"validation_complete={live_method12['validation_complete']} "
                f"elapsed={live_method12['last_invocation_elapsed_seconds']:.3f}s"
            )
            print(
                "  Checkpoint reuse: "
                f"resumed={live_method12['resumed_from_checkpoint']}; "
                "generated WSL i*.iso cache is excluded; live identities and candidate "
                "fingerprints are re-evaluated"
            )
        for report in live_method12.get("reports", []):
            print(
                f"Live {report['method']} {report['case']}: "
                f"status={report['status']} "
                f"candidates={report['live_candidates']}/"
                f"{report['official_candidates']} mode-counts="
                f"{report['mode_count_matches']}/"
                f"{report['mode_files_compared']} missing={len(report['missing'])} "
                f"extra={len(report['extra'])} mismatches={len(report['mismatches'])}"
            )
            for mismatch in report["mismatches"][:10]:
                print(
                    "  MODE MISMATCH "
                    f"{mismatch['candidate']}: official={mismatch['official']} "
                    f"live={mismatch['live']}"
                )
    live12_failed = bool(live_method12) and not bool(live_method12.get("passed"))
    return 0 if summary.failures == 0 and not live12_failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
