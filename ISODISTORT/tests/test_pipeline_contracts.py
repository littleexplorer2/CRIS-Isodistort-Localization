"""Independent invariants for shared geometry and read-only verification."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from backend.input_cif import build_supercell, map_supercell_sites
from backend.method4 import ModeDecomposer
from backend.utils.lattice import centering_multiplicity, centering_primitive_matrix, determinant
from tests.validation.application.batch import run_batch
from tests.validation.application.single import compare_cif
from tests.validation.frontend.cli import main


@pytest.mark.parametrize("letter,index", [("P", 1), ("A", 2), ("B", 2), ("C", 2), ("I", 2), ("F", 4), ("R", 3)])
def test_centering_index_matches_exact_volume_ratio(letter, index):
    # Independent crystallographic conventional-cell coset counts.
    assert centering_multiplicity(letter) == index
    assert abs(determinant(centering_primitive_matrix(letter))) * index == 1


def test_supercell_mapping_conserves_sites_translations_and_species():
    parent = Structure(Lattice.from_parameters(3, 4, 5, 75, 85, 65), ["Si", "O"], [[0, 0, 0], [.2, .3, .4]])
    basis = np.array([[2, 1, 0], [0, 1, 0], [0, 0, 1]])
    child = build_supercell(parent, basis)
    indices, images = map_supercell_sites(parent, child, basis)
    assert np.bincount(indices).tolist() == [2, 2]
    assert np.allclose(child.frac_coords @ basis, parent.frac_coords[indices] + images)
    assert all(site.species == parent[index].species for site, index in zip(child, indices))
    reordered = Structure.from_sites(list(reversed(child)))
    other, translations = map_supercell_sites(parent, reordered, basis)
    assert other.tolist() == indices[::-1].tolist()
    assert np.array_equal(translations, images[::-1])


@pytest.mark.parametrize("species,coords", [(["O"], [[0, 0, 0]]), (["Si"], [[.2, .2, .2]])])
def test_mapping_rejects_missing_parent_images(species, coords):
    parent = Structure(Lattice.cubic(3), ["Si"], [[0, 0, 0]])
    with pytest.raises(ValueError, match="compatible parent sites"):
        map_supercell_sites(parent, Structure(parent.lattice, species, coords), np.eye(3))


def test_mapping_rejects_two_translations_of_one_parent_site():
    parent = Structure(Lattice.cubic(1), ["Si"], [[0, 0, 0]])
    child = Structure(parent.lattice, ["Si"], [[.5, 0, 0]])
    with pytest.raises(ValueError, match="2 compatible parent sites"):
        map_supercell_sites(parent, child, np.eye(3), tolerance=.51)


def test_method4_facade_and_owner_share_the_same_algorithm():
    from backend.method1 import IsoSearchEngine
    assert IsoSearchEngine.method_4_decompose is ModeDecomposer.method_4_decompose


def _pair(tmp_path):
    item, reference = tmp_path / "item", tmp_path / "true"
    item.mkdir()
    reference.mkdir()
    structure = Structure(Lattice.cubic(3), ["Si"], [[0, 0, 0]])
    structure.to(filename=str(item / "pair.cif"), fmt="cif")
    (reference / "pair.cif").write_bytes((item / "pair.cif").read_bytes())
    return item, reference


def test_strict_cli_cannot_bypass_reference_hash(tmp_path, monkeypatch, capsys):
    item, reference = _pair(tmp_path)
    from tests.validation.backend.paths import service
    monkeypatch.setattr(service, "ITEM_DIR", item)
    monkeypatch.setattr(service, "TRUE_DIR", reference)
    before = {path: path.read_bytes() for root in (item, reference) for path in root.iterdir()}
    assert main(["compare", "pair.cif", "--strict", "--reference-sha256", "0" * 64, "--json"]) == 1
    result = json.loads(capsys.readouterr().out)
    assert result["byte_exact"] and not result["reference_hash_valid"] and not result["passed"]
    assert all(path.read_bytes() == data for path, data in before.items())


def _batch(item, reference, manifest=None):
    return argparse.Namespace(local_dir=item, reference_dir=reference, pattern="*.cif",
                              lattice_tol=1e-5, coord_tol=1e-5, scalar_tol=1e-5,
                              ignore_atom_order=False, strict=True, hash_manifest=manifest)


def test_strict_batch_cannot_bypass_reference_hash(tmp_path):
    item, reference = _pair(tmp_path)
    manifest = tmp_path / "hashes.json"
    manifest.write_text(json.dumps({"pair.cif": "0" * 64}), encoding="utf-8")
    results, failed = run_batch(_batch(item, reference, manifest))
    assert failed == 1 and results[0]["passed"] is False


def test_cli_reads_explicit_paths_without_creating_default_directories(tmp_path, monkeypatch, capsys):
    item, reference = _pair(tmp_path)
    from tests.validation.backend.paths import service
    missing = tmp_path / "must-not-be-created"
    monkeypatch.setattr(service, "ITEM_DIR", missing / "item")
    monkeypatch.setattr(service, "TRUE_DIR", missing / "true")
    assert main(["compare", "--local", str(item / "pair.cif"), "--reference", str(reference / "pair.cif"), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["passed"]
    assert main(["batch", "--local-dir", str(item), "--reference-dir", str(reference), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["passed"] == 1
    assert not missing.exists()


def test_empty_batch_is_not_a_successful_validation(tmp_path):
    item, reference = tmp_path / "item", tmp_path / "true"
    item.mkdir()
    reference.mkdir()
    with pytest.raises(ValueError, match="no CIF comparisons executed"):
        run_batch(_batch(item, reference))


def test_cli_configuration_defaults_and_explicit_override(tmp_path, monkeypatch, capsys):
    from tests.validation.backend.config import get_config
    item, reference = _pair(tmp_path)
    with (item / "pair.cif").open("a", encoding="utf-8") as stream:
        stream.write("\n# formatting-only difference\n")
    monkeypatch.setitem(get_config()._cfg["defaults"], "strict", True)
    arguments = ["compare", "--local", str(item / "pair.cif"), "--reference", str(reference / "pair.cif"), "--json"]
    assert main(arguments) == 1
    assert json.loads(capsys.readouterr().out)["structure_equal"]
    assert main(arguments + ["--no-strict"]) == 0
    assert json.loads(capsys.readouterr().out)["passed"]
    assert main(["batch", "--local-dir", str(item), "--reference-dir", str(reference), "--json"]) == 1
    assert json.loads(capsys.readouterr().out)["failed"] == 1


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1.0])
def test_invalid_tolerance_is_rejected(tmp_path, value):
    item, reference = _pair(tmp_path)
    with pytest.raises(ValueError, match="finite and non-negative"):
        compare_cif(item / "pair.cif", reference / "pair.cif", coordinate_tolerance=value)


def test_agent_entry_is_noninteractive_and_does_not_import_runtime():
    project = Path(__file__).resolve().parents[1]
    process = subprocess.run([sys.executable, "-c",
        "import sys; from tests.validation.frontend.cli import main; "
        "assert main([]) == 0; assert 'application.session' not in sys.modules; "
        "assert not any(n.startswith('backend.wrappers') for n in sys.modules)"],
        cwd=project, capture_output=True, text=True, timeout=30)
    assert process.returncode == 0, process.stderr
