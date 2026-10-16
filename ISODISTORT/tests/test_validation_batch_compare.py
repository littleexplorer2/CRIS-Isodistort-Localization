import argparse
import json
from pathlib import Path

from pymatgen.core import Lattice, Structure

from tests.validation.frontend.cli.batch import main
from tests.validation.application.batch import run_batch
from tests.validation.backend.paths.service import BATCH_PAIRING_HINT


def _write_structure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
    structure.to(filename=str(path), fmt="cif")


def _batch_args(local: Path, reference: Path, json_output: bool = True) -> argparse.Namespace:
    return argparse.Namespace(
        local_dir=str(local),
        reference_dir=str(reference),
        pattern="*.cif",
        lattice_tol=1e-5,
        coord_tol=1e-5,
        scalar_tol=1e-5,
        ignore_atom_order=False,
        hash_manifest=None,
        strict=False,
        json=json_output,
    )


def test_batch_compare_matches_relative_paths(tmp_path: Path) -> None:
    local = tmp_path / "item"
    reference = tmp_path / "true"
    _write_structure(local / "case" / "sample.cif")
    _write_structure(reference / "case" / "sample.cif")

    results, failed = run_batch(_batch_args(local, reference))
    assert failed == 0
    assert len(results) == 1
    assert results[0]["passed"]


def test_batch_compare_reports_missing_pair(tmp_path: Path) -> None:
    local = tmp_path / "item"
    reference = tmp_path / "true"
    local.mkdir()
    reference.mkdir()
    _write_structure(local / "only-local.cif")

    results, failed = run_batch(_batch_args(local, reference))
    assert failed == 1
    assert results[0]["issues"] == ["missing matching CIF"]


def test_batch_cli_uses_fixed_compare_dirs(tmp_path: Path, monkeypatch, capsys) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _write_structure(item / "sample.cif")
    _write_structure(true / "sample.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    assert main(["--json"]) == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["total"] == 1
    assert payload["passed"] == 1
    assert payload["failed"] == 0
    assert "compare/true" in captured.err
    assert BATCH_PAIRING_HINT in captured.err


def test_entry_main_batch_subcommand(tmp_path: Path, monkeypatch, capsys) -> None:
    from tests.validation.frontend import cli as module

    item = tmp_path / "item"
    true = tmp_path / "true"
    _write_structure(item / "sample.cif")
    _write_structure(true / "sample.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    assert module.main(["batch", "--json"]) == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["passed"] == 1
    assert BATCH_PAIRING_HINT in captured.err
