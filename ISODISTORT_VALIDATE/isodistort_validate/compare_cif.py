"""Compare a local ISODISTORT CIF with an independent reference CIF."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import spglib
from pymatgen.core import Structure
from pymatgen.io.cif import CifFile, CifParser

from . import compare_paths as cpaths

SPACE_GROUP_NUMBER_TAGS = ("_symmetry_Int_Tables_number", "_space_group_IT_number")
SPACE_GROUP_SYMBOL_TAGS = ("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt")
OCCUPANCY_TAG = "_atom_site_occupancy"
LABEL_TAG = "_atom_site_label"
MAGMOM_TAGS = (
    "_atom_site_moment_crystalaxis_x",
    "_atom_site_moment_crystalaxis_y",
    "_atom_site_moment_crystalaxis_z",
)


@dataclass
class ComparisonResult:
    local_path: str
    reference_path: str
    byte_exact: bool
    text_exact: bool
    structure_equal: bool
    reference_hash_valid: bool | None
    local_sha256: str
    reference_sha256: str
    details: dict[str, Any]

    @property
    def passed(self) -> bool:
        """Semantic equality is the default; formatting is reported separately."""
        return self.structure_equal and self.reference_hash_valid is not False


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8-sig")


def _parse_structure(path: Path) -> Structure:
    parser = CifParser(str(path), occupancy_tolerance=100)
    structures = parser.parse_structures(primitive=False)
    if not structures:
        raise ValueError(f"CIF does not contain a parseable structure: {path}")
    if len(structures) > 1:
        raise ValueError(f"CIF contains {len(structures)} structures; expected one: {path}")
    return structures[0]


def _read_block(path: Path) -> dict[str, Any]:
    cif = CifFile.from_file(path)
    if len(cif.data) != 1:
        raise ValueError(f"CIF must contain exactly one data block: {path}")
    return next(iter(cif.data.values())).data


def _first_value(block: dict[str, Any], tags: tuple[str, ...]) -> Any | None:
    for tag in tags:
        if tag in block:
            value = block[tag]
            return value[0] if isinstance(value, list) else value
    return None


def _numeric(value: Any) -> float:
    text = re.sub(r"\([^)]*\)$", "", str(value).strip().strip("'\""))
    return float(text)


def _numeric_list(block: dict[str, Any], tag: str) -> list[float] | None:
    if tag not in block:
        return None
    values = block[tag]
    if not isinstance(values, list):
        values = [values]
    return [_numeric(value) for value in values]


def _periodic_delta(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    return (second - first) - np.round(second - first)


def _match_atoms(
    local: Structure,
    reference: Structure,
    coordinate_tolerance: float,
    ignore_atom_order: bool,
) -> tuple[list[int] | None, float | None]:
    if len(local) != len(reference):
        return None, None
    if not ignore_atom_order:
        delta = _periodic_delta(local.frac_coords, reference.frac_coords)
        maximum = float(np.max(np.abs(delta))) if len(local) else 0.0
        if maximum > coordinate_tolerance:
            return None, maximum
        return list(range(len(local))), maximum

    candidates: list[tuple[float, int, int]] = []
    for local_index, local_site in enumerate(local):
        for reference_index, reference_site in enumerate(reference):
            if local_site.species_string != reference_site.species_string:
                continue
            delta = _periodic_delta(
                np.asarray(local_site.frac_coords),
                np.asarray(reference_site.frac_coords),
            )
            distance = float(np.max(np.abs(delta)))
            if distance <= coordinate_tolerance:
                candidates.append((distance, local_index, reference_index))

    assignments: dict[int, int] = {}
    used_local: set[int] = set()
    used_reference: set[int] = set()
    for _, local_index, reference_index in sorted(candidates):
        if local_index not in used_local and reference_index not in used_reference:
            assignments[local_index] = reference_index
            used_local.add(local_index)
            used_reference.add(reference_index)
    if len(assignments) != len(local):
        return None, None
    ordered = [assignments[index] for index in range(len(local))]
    distances = [
        distance
        for distance, local_index, reference_index in candidates
        if assignments.get(local_index) == reference_index
    ]
    return ordered, max(distances, default=0.0)


def _spglib_space_group(structure: Structure, symprec: float) -> dict[str, Any]:
    numbers = [site.specie.Z for site in structure]
    dataset = spglib.get_symmetry_dataset(
        (structure.lattice.matrix, structure.frac_coords, numbers), symprec=symprec
    )
    if dataset is None:
        return {"number": None, "symbol": None}
    return {"number": int(dataset.number), "symbol": str(dataset.international)}


def _compare_mapped_scalars(
    local: list[float] | None,
    reference: list[float] | None,
    assignments: list[int] | None,
    tolerance: float,
) -> bool | None:
    if local is None and reference is None:
        return None
    if local is None or reference is None or assignments is None:
        return False
    if len(local) != len(reference) or len(local) != len(assignments):
        return False
    return all(abs(local[i] - reference[assignments[i]]) <= tolerance for i in range(len(local)))


def _compare_mapped_vectors(
    local: list[list[float] | None] | None,
    reference: list[list[float] | None] | None,
    assignments: list[int] | None,
    tolerance: float,
) -> bool | None:
    if local is None and reference is None:
        return None
    if local is None or reference is None or assignments is None:
        return False
    if any(value is None for value in local + reference):
        return False
    local_columns = [value for value in local if value is not None]
    reference_columns = [value for value in reference if value is not None]
    if len(local_columns) != 3 or len(reference_columns) != 3:
        return False
    if any(len(column) != len(assignments) for column in local_columns + reference_columns):
        return False
    return all(
        abs(local_columns[axis][index] - reference_columns[axis][reference_index]) <= tolerance
        for index, reference_index in enumerate(assignments)
        for axis in range(3)
    )


def _compare_mapped_strings(
    local: Any, reference: Any, assignments: list[int] | None
) -> bool | None:
    if local is None and reference is None:
        return None
    if not isinstance(local, list) or not isinstance(reference, list) or assignments is None:
        return False
    if len(local) != len(reference) or len(local) != len(assignments):
        return False
    return all(str(local[i]) == str(reference[assignments[i]]) for i in range(len(local)))


def _declared_space_group(block: dict[str, Any]) -> dict[str, Any]:
    number = _first_value(block, SPACE_GROUP_NUMBER_TAGS)
    symbol = _first_value(block, SPACE_GROUP_SYMBOL_TAGS)
    normalized_symbol = None
    if symbol is not None:
        normalized_symbol = " ".join(str(symbol).strip("'\"").split())
    return {
        "number": int(_numeric(number)) if number is not None else None,
        "symbol": normalized_symbol,
    }


def _compare_optional_metadata(first: Any, second: Any) -> bool | None:
    if first is None and second is None:
        return None
    return first == second


def _validate_declared_space_group(declared: dict[str, Any], inferred: dict[str, Any]) -> bool:
    return declared["number"] is None or declared["number"] == inferred["number"]


def _compare_structures(
    local: Structure,
    reference: Structure,
    local_block: dict[str, Any],
    reference_block: dict[str, Any],
    lattice_tolerance: float,
    coordinate_tolerance: float,
    scalar_tolerance: float,
    ignore_atom_order: bool,
) -> tuple[bool, dict[str, Any]]:
    details: dict[str, Any] = {
        "lattice_tolerance": lattice_tolerance,
        "coordinate_tolerance": coordinate_tolerance,
        "scalar_tolerance": scalar_tolerance,
        "ignore_atom_order": ignore_atom_order,
        "local_atom_count": len(local),
        "reference_atom_count": len(reference),
        "lattice_max_abs_difference": None,
        "coordinate_max_periodic_difference": None,
        "species_order_equal": None,
        "atom_order_equal": None,
        "lattice_equal": None,
        "coordinates_equal": None,
        "occupancies_equal": None,
        "magnetic_moments_equal": None,
        "labels_equal": None,
        "declared_space_group_equal": None,
        "inferred_space_group_equal": None,
        "space_group": {},
        "issues": [],
    }
    if len(local) != len(reference):
        details["issues"].append("atom count differs")
        return False, details

    lattice_max = float(np.max(np.abs(local.lattice.matrix - reference.lattice.matrix)))
    species_local = [site.species_string for site in local]
    species_reference = [site.species_string for site in reference]
    species_order_equal = species_local == species_reference
    assignments, coordinate_max = _match_atoms(
        local, reference, coordinate_tolerance, ignore_atom_order
    )
    coordinates_equal = assignments is not None
    if not coordinates_equal and coordinate_max is None:
        coordinate_max = float("inf")
    lattice_equal = lattice_max <= lattice_tolerance
    atom_order_equal = species_order_equal and not ignore_atom_order

    occupancies_equal = _compare_mapped_scalars(
        _numeric_list(local_block, OCCUPANCY_TAG),
        _numeric_list(reference_block, OCCUPANCY_TAG),
        assignments,
        scalar_tolerance,
    )
    local_moments = [_numeric_list(local_block, tag) for tag in MAGMOM_TAGS]
    reference_moments = [_numeric_list(reference_block, tag) for tag in MAGMOM_TAGS]
    magnetic_present = any(value is not None for value in local_moments + reference_moments)
    magnetic_equal = _compare_mapped_vectors(
        local_moments if magnetic_present else None,
        reference_moments if magnetic_present else None,
        assignments,
        scalar_tolerance,
    )
    labels_equal = _compare_mapped_strings(
        local_block.get(LABEL_TAG), reference_block.get(LABEL_TAG), assignments
    )

    declared_local = _declared_space_group(local_block)
    declared_reference = _declared_space_group(reference_block)
    declared_equal = _compare_optional_metadata(declared_local, declared_reference)
    symprec = max(lattice_tolerance, coordinate_tolerance * max(local.lattice.abc))
    inferred_local = _spglib_space_group(local, symprec)
    inferred_reference = _spglib_space_group(reference, symprec)
    inferred_equal = inferred_local == inferred_reference
    declared_valid = _validate_declared_space_group(declared_local, inferred_local)
    declared_valid = declared_valid and _validate_declared_space_group(
        declared_reference, inferred_reference
    )

    issues: list[str] = []
    if not lattice_equal:
        issues.append("lattice differs beyond lattice tolerance")
    if not species_order_equal and not ignore_atom_order:
        issues.append("species or atom order differs")
    if not coordinates_equal:
        issues.append("fractional coordinates differ beyond tolerance")
    if occupancies_equal is False:
        issues.append("occupancies differ beyond scalar tolerance")
    if magnetic_equal is False:
        issues.append("magnetic moments differ beyond scalar tolerance")
    if labels_equal is False:
        issues.append("atom labels differ")
    if declared_equal is False:
        issues.append("declared space group differs")
    if not inferred_equal:
        issues.append("spglib-inferred space group differs")
    if not declared_valid:
        issues.append("declared space group does not match spglib inference")

    structure_equal = (
        lattice_equal
        and coordinates_equal
        and occupancies_equal is not False
        and magnetic_equal is not False
        and labels_equal is not False
        and declared_equal is not False
        and inferred_equal
    )
    details.update(
        {
            "lattice_max_abs_difference": lattice_max,
            "coordinate_max_periodic_difference": coordinate_max,
            "species_order_equal": species_order_equal,
            "atom_order_equal": atom_order_equal,
            "lattice_equal": lattice_equal,
            "coordinates_equal": coordinates_equal,
            "occupancies_equal": occupancies_equal,
            "magnetic_moments_equal": magnetic_equal,
            "labels_equal": labels_equal,
            "declared_space_group_equal": declared_equal,
            "inferred_space_group_equal": inferred_equal,
            "space_group": {
                "local_declared": declared_local,
                "reference_declared": declared_reference,
                "local_inferred": inferred_local,
                "reference_inferred": inferred_reference,
            },
            "issues": issues,
        }
    )
    return structure_equal, details


def compare_cif(
    local_path: str | Path,
    reference_path: str | Path,
    lattice_tolerance: float = 1e-5,
    coordinate_tolerance: float = 1e-5,
    scalar_tolerance: float = 1e-5,
    ignore_atom_order: bool = False,
    reference_sha256: str | None = None,
) -> ComparisonResult:
    """Compare two CIF files at byte, text, metadata, and structure levels."""
    local = Path(local_path).expanduser().resolve()
    reference = Path(reference_path).expanduser().resolve()
    for path in (local, reference):
        if not path.is_file():
            raise FileNotFoundError(f"CIF file not found: {path}")
    if min(lattice_tolerance, coordinate_tolerance, scalar_tolerance) < 0:
        raise ValueError("tolerances must be non-negative")

    local_bytes = local.read_bytes()
    reference_bytes = reference.read_bytes()
    local_sha = _sha256(local)
    reference_sha = _sha256(reference)
    local_structure = _parse_structure(local)
    reference_structure = _parse_structure(reference)
    local_block = _read_block(local)
    reference_block = _read_block(reference)
    structure_equal, details = _compare_structures(
        local_structure,
        reference_structure,
        local_block,
        reference_block,
        lattice_tolerance,
        coordinate_tolerance,
        scalar_tolerance,
        ignore_atom_order,
    )
    hash_valid = None if reference_sha256 is None else reference_sha.lower() == reference_sha256.lower()
    if hash_valid is False:
        details["issues"].append("reference SHA-256 does not match expected hash")

    return ComparisonResult(
        local_path=str(local),
        reference_path=str(reference),
        byte_exact=local_bytes == reference_bytes,
        text_exact=_read_text(local) == _read_text(reference),
        structure_equal=structure_equal,
        reference_hash_valid=hash_valid,
        local_sha256=local_sha,
        reference_sha256=reference_sha,
        details=details,
    )


def _list_relative_cifs(root: Path, pattern: str = "*.cif") -> list[str]:
    return cpaths.list_relative_cifs(root, pattern)


def _default_relative_path() -> str:
    paired, item_only, true_only = cpaths.pairing_status()
    warning = cpaths.format_unpaired_warning(item_only, true_only)
    if len(paired) == 1:
        return paired[0]
    if not paired:
        raise FileNotFoundError(
            warning
            or (
                "compare/item 与 compare/true 中没有 CIF。"
                "请把官网参考 CIF 放入 compare/true/，把待验证 CIF 放入 compare/item/，"
                "并保持相对路径相同。"
            )
        )
    preview = ", ".join(paired[:8])
    extra = "" if len(paired) <= 8 else f" ... ({len(paired)} files)"
    raise ValueError(
        "compare/item 与 compare/true 中有多对 CIF，请给出相对路径，"
        f"或改用 python main.py batch。可用文件: {preview}{extra}"
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="main.py compare",
        description=(
            "Compare a local ISODISTORT CIF in compare/item with the matching "
            "official CIF in compare/true (same relative path)."
        )
    )
    parser.add_argument(
        "relative_path",
        nargs="?",
        help="CIF path relative to compare/item and compare/true; omit if only one pair exists",
    )
    parser.add_argument("--lattice-tol", type=float, default=1e-5)
    parser.add_argument("--coord-tol", type=float, default=1e-5)
    parser.add_argument("--scalar-tol", type=float, default=1e-5)
    parser.add_argument(
        "--ignore-atom-order",
        action="store_true",
        help=(
            "ignore CIF atom row order and match by element plus fractional coordinates; "
            "default compares row-by-row. Does not ignore real coordinate or species differences"
        ),
    )
    parser.add_argument("--reference-sha256", help="expected SHA-256 for the trusted reference CIF")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="require byte-identical files in addition to semantic equality",
    )
    parser.add_argument(
        "--structure-only",
        action="store_true",
        help="deprecated alias; semantic equality is the default",
    )
    parser.add_argument("--json", action="store_true", help="print the complete result as JSON")
    return parser


def _print_result(result: ComparisonResult, strict: bool) -> None:
    print("=" * 72)
    print("CIF 比对结果")
    print(f"本地文件: {result.local_path}")
    print(f"官网文件: {result.reference_path}")
    print(f"字节完全一致: {'是' if result.byte_exact else '否'}")
    print(f"文本内容一致: {'是' if result.text_exact else '否'}")
    print(f"解析后结构与元数据一致: {'是' if result.structure_equal else '否'}")
    print(f"本地 SHA-256: {result.local_sha256}")
    print(f"官网 SHA-256: {result.reference_sha256}")
    print(f"参考文件哈希校验: {'通过' if result.reference_hash_valid is not False else '失败'}")
    details = result.details
    print("差异摘要:")
    print(f"  原子数: {details['local_atom_count']} / {details['reference_atom_count']}")
    print(f"  晶格最大绝对误差: {details['lattice_max_abs_difference']}")
    print(f"  分数坐标最大周期误差: {details['coordinate_max_periodic_difference']}")
    print(f"  占据率一致: {details['occupancies_equal']}")
    print(f"  磁矩一致: {details['magnetic_moments_equal']}")
    print(f"  空间群声明一致: {details['declared_space_group_equal']}")
    print(f"  spglib 推断一致: {details['inferred_space_group_equal']}")
    if details["issues"]:
        print("差异原因:")
        for issue in details["issues"]:
            print(f"  - {issue}")
    passed = result.byte_exact and result.structure_equal if strict else result.passed
    print("结论: " + ("PASS" if passed else "FAIL"))
    if result.structure_equal and not result.byte_exact:
        print("提示: 结构和关键元数据一致，但 CIF 排版或非关键文本不同。")


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    cpaths.warn_unpaired_filenames()
    try:
        relative = args.relative_path or _default_relative_path()
        local_cif, reference_cif = cpaths.resolve_pair(relative)
        result = compare_cif(
            local_cif,
            reference_cif,
            args.lattice_tol,
            args.coord_tol,
            args.scalar_tol,
            args.ignore_atom_order,
            args.reference_sha256,
        )
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2

    passed = result.byte_exact and result.structure_equal if args.strict else result.passed
    if args.json:
        payload = asdict(result)
        payload["passed"] = passed
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        _print_result(result, args.strict)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
