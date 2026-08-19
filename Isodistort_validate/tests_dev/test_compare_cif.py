from pathlib import Path

import pytest
from pymatgen.core import Lattice, Structure

from compare_cif import compare_cif


def _write_structure(path: Path, shift: float = 0.0) -> None:
    structure = Structure(
        Lattice.cubic(3.5),
        ["Si", "Si"],
        [[0, 0, 0], [0.25 + shift, 0.25, 0.25]],
    )
    structure.to(filename=str(path), fmt="cif")


def test_identical_cif_passes(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    _write_structure(first)
    second.write_bytes(first.read_bytes())

    result = compare_cif(first, second)

    assert result.byte_exact
    assert result.text_exact
    assert result.structure_equal
    assert result.passed


def test_format_difference_can_still_have_same_structure(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    _write_structure(first)
    second.write_text("\n" + first.read_text(encoding="utf-8"), encoding="utf-8")

    result = compare_cif(first, second)

    assert not result.byte_exact
    assert result.structure_equal
    assert result.passed


def test_strict_mode_can_still_detect_format_difference(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    _write_structure(first)
    second.write_text("\n" + first.read_text(encoding="utf-8"), encoding="utf-8")

    result = compare_cif(first, second)

    assert result.structure_equal
    assert not result.byte_exact


def test_reference_hash_is_part_of_the_pass_condition(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    _write_structure(first)

    result = compare_cif(first, first, reference_sha256="0" * 64)

    assert result.structure_equal
    assert result.reference_hash_valid is False
    assert not result.passed


def test_coordinate_difference_fails(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    _write_structure(first)
    _write_structure(second, shift=0.01)

    result = compare_cif(first, second)

    assert not result.structure_equal
    assert "fractional coordinates differ beyond tolerance" in result.details["issues"]


def test_missing_file_is_reported(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        compare_cif(tmp_path / "missing.cif", tmp_path / "other.cif")
