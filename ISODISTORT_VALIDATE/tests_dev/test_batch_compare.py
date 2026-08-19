import json
from pathlib import Path

from pymatgen.core import Lattice, Structure

from batch_compare import main


def _write_structure(path: Path) -> None:
    structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
    structure.to(filename=str(path), fmt="cif")


def test_batch_compare_matches_relative_paths(tmp_path: Path, capsys) -> None:
    local = tmp_path / "local"
    reference = tmp_path / "reference"
    (local / "case").mkdir(parents=True)
    (reference / "case").mkdir(parents=True)
    _write_structure(local / "case" / "sample.cif")
    _write_structure(reference / "case" / "sample.cif")

    assert main([str(local), str(reference), "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["total"] == 1
    assert payload["passed"] == 1
    assert payload["failed"] == 0


def test_batch_compare_reports_missing_pair(tmp_path: Path, capsys) -> None:
    local = tmp_path / "local"
    reference = tmp_path / "reference"
    local.mkdir()
    reference.mkdir()
    _write_structure(local / "only-local.cif")

    assert main([str(local), str(reference), "--json"]) == 1
    payload = json.loads(capsys.readouterr().out)
    assert payload["failed"] == 1
    assert payload["results"][0]["issues"] == ["missing matching CIF"]
