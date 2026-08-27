import json
from pathlib import Path

import pytest
from pymatgen.core import Lattice, Structure

from isodistort_validate.compare_cif import compare_cif, main
from isodistort_validate.compare_paths import (
    ensure_compare_dirs,
    format_unpaired_warning,
    pairing_status,
    resolve_pair,
)


def _write_structure(path: Path, shift: float = 0.0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def test_resolve_pair_rejects_absolute_and_parent_paths() -> None:
    with pytest.raises(ValueError):
        resolve_pair(r"C:\absolute\file.cif")
    with pytest.raises(ValueError):
        resolve_pair("../secret.cif")


def test_ensure_compare_dirs_creates_item_and_true(tmp_path: Path, monkeypatch) -> None:
    compare = tmp_path / "compare"
    monkeypatch.setattr("isodistort_validate.compare_paths.COMPARE_ROOT", compare)
    monkeypatch.setattr("isodistort_validate.compare_paths.ITEM_DIR", compare / "item")
    monkeypatch.setattr("isodistort_validate.compare_paths.TRUE_DIR", compare / "true")
    item, true = ensure_compare_dirs()
    assert item.is_dir()
    assert true.is_dir()
    assert item.name == "item"
    assert true.name == "true"


def test_cli_compares_fixed_compare_dirs(tmp_path: Path, monkeypatch, capsys) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _write_structure(item / "pair.cif")
    _write_structure(true / "pair.cif")
    monkeypatch.setattr("isodistort_validate.compare_paths.ITEM_DIR", item)
    monkeypatch.setattr("isodistort_validate.compare_paths.TRUE_DIR", true)

    assert main(["pair.cif", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["passed"] is True


def test_cli_auto_picks_single_pair(tmp_path: Path, monkeypatch, capsys) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _write_structure(item / "only.cif")
    _write_structure(true / "only.cif")
    monkeypatch.setattr("isodistort_validate.compare_paths.ITEM_DIR", item)
    monkeypatch.setattr("isodistort_validate.compare_paths.TRUE_DIR", true)

    assert main(["--json"]) == 0
    assert "only.cif" in capsys.readouterr().out


def test_pairing_status_lists_unmatched_filenames(tmp_path: Path, monkeypatch) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _write_structure(item / "X4- P3 CIF.cif")
    _write_structure(true / "subgroup.cif")
    monkeypatch.setattr("isodistort_validate.compare_paths.ITEM_DIR", item)
    monkeypatch.setattr("isodistort_validate.compare_paths.TRUE_DIR", true)

    paired, item_only, true_only = pairing_status()
    assert paired == []
    assert item_only == ["X4- P3 CIF.cif"]
    assert true_only == ["subgroup.cif"]
    text = format_unpaired_warning(item_only, true_only)
    assert text is not None
    assert "请先修改文件名" in text
    assert "X4- P3 CIF.cif" in text
    assert "subgroup.cif" in text


def test_cli_warns_immediately_when_filenames_do_not_pair(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _write_structure(item / "X4- P3 CIF.cif")
    _write_structure(true / "subgroup.cif")
    monkeypatch.setattr("isodistort_validate.compare_paths.ITEM_DIR", item)
    monkeypatch.setattr("isodistort_validate.compare_paths.TRUE_DIR", true)

    assert main(["--json"]) == 2
    err = capsys.readouterr().err
    assert "发现文件名对不上" in err
    assert "X4- P3 CIF.cif" in err
    assert "subgroup.cif" in err
