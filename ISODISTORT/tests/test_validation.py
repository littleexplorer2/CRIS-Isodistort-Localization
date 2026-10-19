"""Tests for test_validation; migrated without changing scientific expectations."""
from __future__ import annotations

from pathlib import Path
from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests.support_paths import TESTS_ROOT
from tests.validation.application.batch import run_batch
from tests.validation.application.single import compare_cif
from tests.validation.backend.config import CONFIG_PATH
from tests.validation.backend.config import get_config
from tests.validation.backend.config import service as config_loader
from tests.validation.backend.paths.service import BATCH_PAIRING_HINT
from tests.validation.backend.paths.service import COMPARE_ROOT
from tests.validation.backend.paths.service import ITEM_DIR
from tests.validation.backend.paths.service import TRUE_DIR
from tests.validation.backend.paths.service import VALIDATE_ROOT
from tests.validation.backend.paths.service import format_unpaired_warning
from tests.validation.backend.paths.service import pairing_status
from tests.validation.backend.paths.service import resolve_pair
from tests.validation.backend.provenance import FileSnapshot
from tests.validation.backend.provenance import build_provenance_state
from tests.validation.backend.provenance import changed_provenance_surfaces
from tests.validation.backend.provenance import production_source_paths
from tests.validation.backend.provenance import validation_source_paths
from tests.validation.backend.provenance import write_json_atomic
from tests.validation.frontend.cli import main as _pipeline_contracts_main
from tests.validation.frontend.cli.batch import main as _validation_batch_compare_main
from tests.validation.frontend.cli.compare import main as _validation_compare_cif_main
from types import SimpleNamespace
import argparse
import hashlib
import importlib
import json
import pytest
import subprocess
import sys


# Original coverage: test_manual_signature_paths.py.

PROJECT_ROOT = TESTS_ROOT.parent

MANUAL_ROOT = PROJECT_ROOT / "tests"

def test_manual_signature_targets_exist_in_the_three_layer_layout() -> None:
    expected = (
        PROJECT_ROOT / "application" / "session" / "core_api.py",
        PROJECT_ROOT / "backend" / "wrappers" / "iso_wrapper.py",
        PROJECT_ROOT / "backend" / "wrappers" / "smodes_wrapper.py",
        PROJECT_ROOT / "backend" / "input_cif" / "__init__.py",
        PROJECT_ROOT / "backend" / "input_cif" / "cif_io.py",
        PROJECT_ROOT / "backend" / "utils" / "affine_embeddings.py",
        PROJECT_ROOT / "backend" / "method1" / "search.py",
        PROJECT_ROOT / "backend" / "method2" / "search.py",
        PROJECT_ROOT / "backend" / "method3" / "search.py",
        PROJECT_ROOT / "backend" / "utils" / "k_lattice.py",
        PROJECT_ROOT / "backend" / "models" / "phase_path.py",
        PROJECT_ROOT / "backend" / "tables" / "spacegroup_meta.py",
        PROJECT_ROOT / "backend" / "method2" / "distortion_mapper.py",
        PROJECT_ROOT / "backend" / "method2" / "superspace.py",
        PROJECT_ROOT / "backend" / "method3" / "inverse_landau.py",
        PROJECT_ROOT / "backend" / "method3" / "inverse_landau_adapter.py",
        PROJECT_ROOT / "backend" / "method4" / "distortion_engine.py",
        PROJECT_ROOT / "backend" / "method4" / "decomposition.py",
        PROJECT_ROOT / "backend" / "method4" / "strain.py",
        PROJECT_ROOT / "backend" / "export" / "isodistort_cif.py",
        MANUAL_ROOT / "official_html_resolver.py",
        MANUAL_ROOT / "validate_method_outputs.py",
    )

    missing = [str(path) for path in expected if not path.is_file()]
    assert not missing, f"manual signature targets are missing: {missing}"

def test_manual_audits_do_not_sign_retired_package_paths() -> None:
    retired_fragments = (
        '/ "backend" / "distortion" /',
        '/ "backend" / "backend" /',
        '/ "backend" / "structure" /',
        '/ "backend" / "io" /',
        "ISODISTORT/isocore/",
        "ISODISTORT/tests_dev/",
    )
    offenders: list[str] = []
    tools = sorted(path for path in MANUAL_ROOT.glob("*.py")
                   if path.name.startswith(("audit_", "compare_", "diagnose_", "prepare_", "validate_", "verify_"))
                   or path.name in {"method4_provenance.py", "method3_affine_equivalence.py",
                                    "official_html_resolver.py", "register_official_archive_reuse.py"})
    assert tools and MANUAL_ROOT / "validate_method_outputs.py" in tools
    for path in tools:
        source = path.read_text(encoding="utf-8-sig")
        for fragment in retired_fragments:
            if fragment in source:
                offenders.append(f"{path.name}: {fragment}")

    assert not offenders, "retired manual signature paths remain: " + "; ".join(offenders)

def test_comparison_signatures_cover_the_relocated_scientific_implementation():
    from tests.validation.backend.provenance import validation_source_paths

    paths = validation_source_paths(PROJECT_ROOT)
    assert "cif-validator:application/single/service.py" in paths
    assert "cif-validator:backend/cif/service.py" in paths
    assert "cif-validator:backend/comparison/service.py" in paths
    assert "cif-validator:backend/results/__init__.py" in paths
    assert "cif-validator:backend/provenance/service.py" in paths
    assert all(path.is_file() for path in paths.values())
    assert not any(path.name == "compare_cif.py" for path in paths.values())

def test_numeric_audit_inventory_includes_shared_verification_policy_and_snapshot():
    from tests.audit_method12_numeric_semantics import _file_inventory

    paths = set(_file_inventory(PROJECT_ROOT / "output", []))
    assert PROJECT_ROOT / "tests/validation/backend/results/model.py" in paths
    assert PROJECT_ROOT / "tests/validation/backend/provenance/service.py" in paths

def _preparer_fixture(tmp_path, monkeypatch, name):
    """Sign a disposable copy of actual production sources, without executing it."""
    module = importlib.import_module(f"tests.{name}")
    root = tmp_path / "ISODISTORT"
    for layer in ("backend", "application"):
        for source in (PROJECT_ROOT / layer).rglob("*.py"):
            target = root / source.relative_to(PROJECT_ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())
    script = root / "tests" / f"{name}.py"
    script.parent.mkdir(parents=True)
    script.write_bytes(Path(module.__file__).read_bytes())
    config = root / "resources/config/settings.yaml"
    config.parent.mkdir(parents=True)
    config.write_bytes((PROJECT_ROOT / "resources/config/settings.yaml").read_bytes())
    inputs = tmp_path / "experiment_data"
    inputs.mkdir()
    for parent in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        (inputs / parent).write_text("independent signature test input\n", encoding="utf-8")
    monkeypatch.setattr(module, "__file__", str(script))
    monkeypatch.setattr(module, "PACKAGE_ROOT", root)
    monkeypatch.setattr(module, "WORKSPACE", tmp_path)
    return module, root

@pytest.mark.parametrize("name", ["prepare_method3_manifest", "prepare_method4_inputs"])
@pytest.mark.parametrize("relative", [
    "backend/utils/affine_embeddings.py", "backend/wrappers/base_wrapper.py",
    "application/session/publication.py", "backend/export/child_frame.py",
    "backend/models/phase_path.py", "backend/utils/k_lattice.py",
    "resources/config/settings.yaml",
])
def test_preparer_signatures_change_with_migrated_implementations(tmp_path, monkeypatch, name, relative):
    module, root = _preparer_fixture(tmp_path, monkeypatch, name)
    before = module._signature()
    target = root / relative
    target.write_bytes(target.read_bytes() + b"\n# changed implementation\n")
    assert module._signature() != before

def test_method3_changed_affine_implementation_invalidates_pass_checkpoint(tmp_path, monkeypatch):
    module, root = _preparer_fixture(tmp_path, monkeypatch, "prepare_method3_manifest")
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(json.dumps({
        "signature": module._signature(), "results": {"case": {"status": "pass"}},
    }), encoding="utf-8")
    target = root / "backend/utils/affine_embeddings.py"
    target.write_bytes(target.read_bytes() + b"\n# changed affine calculation\n")
    calls = []
    subgroup = SimpleNamespace(
        irrep_label="GM1", opd_symbol="P1", space_group_number=1,
        space_group_symbol="P1", basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0], official_fields=lambda: {"k_active": "GM"},
    )

    def search(**_kwargs):
        calls.append("searched")
        return [SimpleNamespace(subgroup=subgroup)]

    monkeypatch.setattr(module, "IsoDistort", lambda: SimpleNamespace(
        load_structure=lambda _path: None, search_method_3=search,
    ))
    monkeypatch.setattr(module, "CHECKPOINT_PATH", checkpoint)
    monkeypatch.setattr(module, "MANIFEST_PATH", tmp_path / "manifest.json")
    monkeypatch.setattr(module, "CASES", {"EuAl4 Parent.cif": [{
        "id": "case", "expected_candidate": {"irrep": "GM1", "opd": "P1", "sg": 1},
        "distortion_types": ["displacive"], "space_group_type": 1,
        "supercell_basis": subgroup.basis_vectors, "direct_sublattice_centering": "P",
        "lattice_type": "direct", "category": "identity",
    }]})
    assert module._run_preflight({"signature": module._signature()}) == 0
    assert calls == ["searched"]

def test_method4_changed_affine_implementation_refuses_old_frozen_input_batch(tmp_path, monkeypatch):
    module, root = _preparer_fixture(tmp_path, monkeypatch, "prepare_method4_inputs")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"signature": module._signature(), "parents": {}}), encoding="utf-8")
    before = manifest.read_bytes()
    target = root / "backend/utils/affine_embeddings.py"
    target.write_bytes(target.read_bytes() + b"\n# changed affine calculation\n")
    monkeypatch.setattr(module, "MANIFEST_PATH", manifest)
    monkeypatch.setattr(module, "PARENTS", {})
    monkeypatch.setattr("sys.argv", ["prepare_method4_inputs"])
    with pytest.raises(RuntimeError, match="different source signature"):
        module.main()
    assert manifest.read_bytes() == before

# Original coverage: test_pipeline_contracts.py.

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
    assert _pipeline_contracts_main(["compare", "pair.cif", "--strict", "--reference-sha256", "0" * 64, "--json"]) == 1
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
    assert _pipeline_contracts_main(["compare", "--local", str(item / "pair.cif"), "--reference", str(reference / "pair.cif"), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["passed"]
    assert _pipeline_contracts_main(["batch", "--local-dir", str(item), "--reference-dir", str(reference), "--json"]) == 0
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
    assert _pipeline_contracts_main(arguments) == 1
    assert json.loads(capsys.readouterr().out)["structure_equal"]
    assert _pipeline_contracts_main(arguments + ["--no-strict"]) == 0
    assert json.loads(capsys.readouterr().out)["passed"]
    assert _pipeline_contracts_main(["batch", "--local-dir", str(item), "--reference-dir", str(reference), "--json"]) == 1
    assert json.loads(capsys.readouterr().out)["failed"] == 1

@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1.0])
def test_invalid_tolerance_is_rejected(tmp_path, value):
    item, reference = _pair(tmp_path)
    with pytest.raises(ValueError, match="finite and non-negative"):
        compare_cif(item / "pair.cif", reference / "pair.cif", coordinate_tolerance=value)

def test_agent_entry_is_noninteractive_and_does_not_import_runtime():
    project = TESTS_ROOT.parent
    process = subprocess.run([sys.executable, "-c",
        "import sys; from tests.validation.frontend.cli import main; "
        "assert main([]) == 0; assert 'application.session' not in sys.modules; "
        "assert not any(n.startswith('backend.wrappers') for n in sys.modules)"],
        cwd=project, capture_output=True, text=True, timeout=30)
    assert process.returncode == 0, process.stderr

# Original coverage: test_validation_batch_compare.py.

def _validation_batch_compare_write_structure(path: Path) -> None:
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
    _validation_batch_compare_write_structure(local / "case" / "sample.cif")
    _validation_batch_compare_write_structure(reference / "case" / "sample.cif")

    results, failed = run_batch(_batch_args(local, reference))
    assert failed == 0
    assert len(results) == 1
    assert results[0]["passed"]

def test_batch_compare_reports_missing_pair(tmp_path: Path) -> None:
    local = tmp_path / "item"
    reference = tmp_path / "true"
    local.mkdir()
    reference.mkdir()
    _validation_batch_compare_write_structure(local / "only-local.cif")

    results, failed = run_batch(_batch_args(local, reference))
    assert failed == 1
    assert results[0]["issues"] == ["missing matching CIF"]

def test_strict_batch_rejects_missing_hash_entries(tmp_path: Path) -> None:
    local = tmp_path / "item"
    reference = tmp_path / "true"
    for name in ("signed.cif", "unsigned.cif"):
        _validation_batch_compare_write_structure(local / name)
        _validation_batch_compare_write_structure(reference / name)
    manifest = tmp_path / "hashes.json"
    manifest.write_text(json.dumps({
        "signed.cif": hashlib.sha256((reference / "signed.cif").read_bytes()).hexdigest(),
    }), encoding="utf-8")
    args = _batch_args(local, reference)
    args.strict = True
    args.hash_manifest = manifest
    results, failed = run_batch(args)
    assert failed == 1
    by_name = {row["relative_path"]: row for row in results}
    assert by_name["signed.cif"]["passed"] is True
    assert by_name["unsigned.cif"]["passed"] is False
    assert by_name["unsigned.cif"]["verification"]["status"] == "inconclusive"

def test_batch_cli_uses_fixed_compare_dirs(tmp_path: Path, monkeypatch, capsys) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _validation_batch_compare_write_structure(item / "sample.cif")
    _validation_batch_compare_write_structure(true / "sample.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    assert _validation_batch_compare_main(["--json"]) == 0
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
    _validation_batch_compare_write_structure(item / "sample.cif")
    _validation_batch_compare_write_structure(true / "sample.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    assert module.main(["batch", "--json"]) == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["passed"] == 1
    assert BATCH_PAIRING_HINT in captured.err

# Original coverage: test_validation_compare_cif.py.

def _validation_compare_cif_write_structure(path: Path, shift: float = 0.0) -> None:
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
    _validation_compare_cif_write_structure(first)
    second.write_bytes(first.read_bytes())

    result = compare_cif(first, second)

    assert result.byte_exact
    assert result.text_exact
    assert result.structure_equal
    assert result.passed

def test_format_difference_can_still_have_same_structure(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    _validation_compare_cif_write_structure(first)
    second.write_text("\n" + first.read_text(encoding="utf-8"), encoding="utf-8")

    result = compare_cif(first, second)

    assert not result.byte_exact
    assert result.structure_equal
    assert result.passed

def test_ignore_order_compares_expanded_occupancies_not_asu_label_rows(
    tmp_path: Path,
) -> None:
    """ASU metadata lists must not be indexed by expanded-structure indices."""
    header = """data_ordered
_cell_length_a 4.1
_cell_length_b 5.2
_cell_length_c 6.3
_cell_angle_alpha 81
_cell_angle_beta 93
_cell_angle_gamma 104
_symmetry_space_group_name_H-M 'P -1'
_symmetry_Int_Tables_number 2
loop_
_space_group_symop_operation_xyz
'x,y,z'
'-x,-y,-z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
"""
    local = tmp_path / "local.cif"
    reference = tmp_path / "reference.cif"
    local.write_text(
        header + "Si_local Si 0.10 0.20 0.30 1\nO_local O 0.27 0.36 0.41 1\n",
        encoding="utf-8",
    )
    reference.write_text(
        header + "O_ref O 0.27 0.36 0.41 1\nSi_ref Si 0.10 0.20 0.30 1\n",
        encoding="utf-8",
    )

    result = compare_cif(local, reference, ignore_atom_order=True)

    assert result.structure_equal
    assert result.details["local_atom_count"] == 4  # two P-1 ASU rows expand to four
    assert result.details["occupancies_equal"] is True
    assert result.details["labels_equal"] is None
    assert "occupancies differ beyond scalar tolerance" not in result.details["issues"]
    assert "atom labels differ" not in result.details["issues"]

def test_ignore_order_uses_complete_bipartite_matching(tmp_path: Path) -> None:
    """A flexible site must be reassigned so a constrained site can match."""
    local = tmp_path / "local.cif"
    reference = tmp_path / "reference.cif"
    Structure(
        Lattice.cubic(3.5),
        ["C", "C"],
        [[0.10, 0.0, 0.0], [0.85, 0.0, 0.0]],
    ).to(filename=str(local), fmt="cif")
    Structure(
        Lattice.cubic(3.5),
        ["C", "C"],
        [[0.00, 0.0, 0.0], [0.30, 0.0, 0.0]],
    ).to(filename=str(reference), fmt="cif")

    result = compare_cif(
        local,
        reference,
        coordinate_tolerance=0.21,
        ignore_atom_order=True,
    )

    assert result.details["coordinates_equal"] is True

def test_strict_mode_can_still_detect_format_difference(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    _validation_compare_cif_write_structure(first)
    second.write_text("\n" + first.read_text(encoding="utf-8"), encoding="utf-8")

    result = compare_cif(first, second)

    assert result.structure_equal
    assert not result.byte_exact

def test_reference_hash_is_part_of_the_pass_condition(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    _validation_compare_cif_write_structure(first)

    result = compare_cif(first, first, reference_sha256="0" * 64)

    assert result.structure_equal
    assert result.reference_hash_valid is False
    assert not result.passed

def test_strict_cli_requires_trusted_reference_hash(tmp_path: Path, capsys) -> None:
    local = tmp_path / "local.cif"
    reference = tmp_path / "reference.cif"
    _validation_compare_cif_write_structure(local)
    reference.write_bytes(local.read_bytes())
    arguments = ["--local", str(local), "--reference", str(reference), "--strict", "--json"]

    assert _validation_compare_cif_main(arguments) == 1
    payload = json.loads(capsys.readouterr().out)
    assert payload["structure_equal"] is True
    assert payload["byte_exact"] is True
    assert payload["passed"] is False
    assert payload["verification"]["status"] == "inconclusive"

    expected = hashlib.sha256(reference.read_bytes()).hexdigest()
    assert _validation_compare_cif_main([*arguments, "--reference-sha256", expected]) == 0
    assert json.loads(capsys.readouterr().out)["verification"]["status"] == "pass"

def test_comparison_uses_one_captured_content_per_file(tmp_path: Path, monkeypatch) -> None:
    local = tmp_path / "local.cif"
    reference = tmp_path / "reference.cif"
    _validation_compare_cif_write_structure(local)
    original = local.read_bytes()
    reference.write_bytes(original)
    original_read = Path.read_bytes
    calls = []

    def capture_then_replace(path):
        content = original_read(path)
        if path in (local, reference):
            calls.append(path)
            path.write_text("not a CIF anymore", encoding="utf-8")
        return content

    monkeypatch.setattr(Path, "read_bytes", capture_then_replace)
    expected = hashlib.sha256(original).hexdigest()
    result = compare_cif(local, reference, reference_sha256=expected)
    assert result.passed
    assert result.byte_exact
    assert result.local_sha256 == result.reference_sha256 == expected
    assert calls == [local, reference]

def test_coordinate_difference_fails(tmp_path: Path) -> None:
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    _validation_compare_cif_write_structure(first)
    _validation_compare_cif_write_structure(second, shift=0.01)

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

def test_pairing_missing_directories_does_not_write_inputs(tmp_path: Path, monkeypatch) -> None:
    compare = tmp_path / "compare"
    monkeypatch.setattr("tests.validation.backend.paths.service.COMPARE_ROOT", compare)
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", compare / "item")
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", compare / "true")
    assert pairing_status() == ([], [], [])
    assert resolve_pair("a.cif") == (compare / "item" / "a.cif", compare / "true" / "a.cif")
    assert not compare.exists()

def test_cli_compares_fixed_compare_dirs(tmp_path: Path, monkeypatch, capsys) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _validation_compare_cif_write_structure(item / "pair.cif")
    _validation_compare_cif_write_structure(true / "pair.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    assert _validation_compare_cif_main(["pair.cif", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["passed"] is True

def test_cli_auto_picks_single_pair(tmp_path: Path, monkeypatch, capsys) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _validation_compare_cif_write_structure(item / "only.cif")
    _validation_compare_cif_write_structure(true / "only.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    assert _validation_compare_cif_main(["--json"]) == 0
    assert "only.cif" in capsys.readouterr().out

def test_pairing_status_lists_unmatched_filenames(tmp_path: Path, monkeypatch) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _validation_compare_cif_write_structure(item / "X4- P3 CIF.cif")
    _validation_compare_cif_write_structure(true / "subgroup.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    paired, item_only, true_only = pairing_status()
    assert paired == []
    assert item_only == ["X4- P3 CIF.cif"]
    assert true_only == ["subgroup.cif"]
    text = format_unpaired_warning(item_only, true_only)
    assert text is not None
    assert "不要修改参考原件" in text
    assert "X4- P3 CIF.cif" in text
    assert "subgroup.cif" in text

def test_cli_warns_immediately_when_filenames_do_not_pair(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    item = tmp_path / "item"
    true = tmp_path / "true"
    _validation_compare_cif_write_structure(item / "X4- P3 CIF.cif")
    _validation_compare_cif_write_structure(true / "subgroup.cif")
    monkeypatch.setattr("tests.validation.backend.paths.service.ITEM_DIR", item)
    monkeypatch.setattr("tests.validation.backend.paths.service.TRUE_DIR", true)

    assert _validation_compare_cif_main(["--json"]) == 2
    err = capsys.readouterr().err
    assert "发现文件名对不上" in err
    assert "X4- P3 CIF.cif" in err
    assert "subgroup.cif" in err

# Original coverage: test_validation_config.py.

def test_settings_yaml_defines_compare_dirs() -> None:
    cfg = get_config()
    assert CONFIG_PATH.is_file()
    assert cfg.project_root == VALIDATE_ROOT
    assert cfg.item_dir == ITEM_DIR
    assert cfg.true_dir == TRUE_DIR
    assert ITEM_DIR.parent == COMPARE_ROOT
    assert ITEM_DIR.name == "item"
    assert TRUE_DIR.name == "true"
    assert cfg.pattern == "*.cif"
    assert cfg.lattice_tolerance == 1e-5
    assert cfg.coordinate_tolerance == 1e-5
    assert cfg.scalar_tolerance == 1e-5
    assert cfg.ignore_atom_order is False
    assert cfg.strict is False

def test_config_path_override_keeps_relative_paths_beside_selected_file(tmp_path, monkeypatch) -> None:
    selected = tmp_path / "custom-config" / "settings.yaml"
    monkeypatch.setattr(config_loader, "CONFIG_PATH", selected)

    assert get_config().resolve_path("../custom-compare") == tmp_path / "custom-compare"

# Original coverage: test_validation_provenance.py.

def _write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path

def _fixture(tmp_path: Path):
    root = tmp_path / "ISODISTORT"
    config = _write(root / "resources/config/settings.yaml", json.dumps({
        "isobyu": {"bin_dir": "../../runtime", "data_dir": "../../runtime"},
    }))
    _write(root / "backend/method3/solver.py", "algorithm = 1\n")
    _write(root / "application/session/service.py", "contract = 1\n")
    _write(root / "pyproject.toml", "[project]\nname='fixture'\n")
    _write(root / "requirements.txt", "numpy\n")
    _write(root / "runtime/iso", "independent-program-v1")
    _write(root / "runtime/data_isotropy.txt", "generic-data-v1")
    source = _write(tmp_path / "input/parent.cif", "independent-input-v1")

    def state(**kwargs):
        return build_provenance_state(
            schema=1, kind="fixture-validation", workspace=tmp_path, package_root=root,
            inputs=[("parent", source)], metadata={"selection": "case"}, **kwargs,
        )

    return root, config, source, state

@pytest.mark.parametrize("relative", [
    "backend/method3/solver.py", "application/session/service.py", "requirements.txt",
    "resources/config/settings.yaml", "runtime/iso", "runtime/data_isotropy.txt",
])
def test_common_provenance_detects_all_executable_input_surfaces(tmp_path, relative):
    root, _config, _source, state = _fixture(tmp_path)
    before = state()
    target = root / relative
    target.write_text(target.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    after = state()
    assert before["signature"] != after["signature"]
    assert changed_provenance_surfaces(before, after)

def test_common_provenance_ignores_outputs_but_tracks_new_code_and_runtime(tmp_path, monkeypatch):
    root, _config, _source, state = _fixture(tmp_path)
    original = state()
    _write(root / "output/cache/result.json", "generated")
    _write(root / "backend/method3/__pycache__/solver.pyc", "cache")
    assert state()["signature"] == original["signature"]
    _write(root / "backend/method3/new_solver.py", "algorithm = 2")
    added = state()
    assert added["signature"] != original["signature"]
    monkeypatch.setattr("tests.validation.backend.provenance.service.importlib.metadata.version", lambda _name: "fixture-v2")
    assert "runtime" in changed_provenance_surfaces(added, state())

def test_provenance_signs_the_same_captured_inputs_used_for_parsing(tmp_path):
    _root, _config, source, state = _fixture(tmp_path)
    captured = FileSnapshot.capture(source)
    source.write_text("changed after parsing", encoding="utf-8")
    frozen = state(captured_inputs=[captured])
    record = next(row for row in frozen["inputs"] if row["label"] == "parent")
    assert record["sha256"] == captured.sha256
    assert record["bytes"] == len(captured.content)
    assert "parent" in changed_provenance_surfaces(frozen, state())

def test_registry_tracks_new_comparison_modules_without_a_hardcoded_file_list(tmp_path):
    root = tmp_path / "ISODISTORT"
    first = _write(root / "tests/validation/backend/cif/service.py", "parser = 1")
    assert set(validation_source_paths(root).values()) == {first}
    second = _write(root / "tests/validation/backend/provenance/service.py", "signer = 1")
    _write(root / "tests/validation/backend/provenance/__pycache__/service.pyc", "cache")
    assert set(validation_source_paths(root).values()) == {first, second}

def test_common_provenance_consumes_the_shared_production_registry(tmp_path, monkeypatch):
    from tests.validation.backend.provenance import service

    root, _config, _source, state = _fixture(tmp_path)
    expected = production_source_paths(root)
    first = state()
    actual = {row["label"] for row in first["inputs"] if row["label"].startswith("production:")}
    assert actual == set(expected)
    added = _write(root / "application/session/publication.py", "publication = 1")
    assert added in production_source_paths(root).values()
    after_add = state()
    assert "production:application/session/publication.py" in changed_provenance_surfaces(first, after_add)
    added.unlink()
    assert state()["signature"] == first["signature"]
    monkeypatch.setattr(service, "production_source_paths", lambda _root: {})
    assert not any(row["label"].startswith("production:") for row in state()["inputs"])

def test_atomic_report_failure_preserves_previous_content_and_removes_temporary(tmp_path, monkeypatch):
    path = _write(tmp_path / "report.json", '{"complete": true}\n')
    before = path.read_bytes()
    with pytest.raises(TypeError):
        write_json_atomic(path, {"unsupported": object()})
    assert path.read_bytes() == before

    def fail_publish(_source, _target):
        raise OSError("simulated publication failure")

    monkeypatch.setattr("tests.validation.backend.provenance.service.os.replace", fail_publish)
    with pytest.raises(OSError, match="publication failure"):
        write_json_atomic(path, {"complete": False})
    assert path.read_bytes() == before
    assert sorted(item.name for item in tmp_path.iterdir()) == ["report.json"]

def test_method3_uses_shared_provenance_for_config_and_iso_inputs(tmp_path, monkeypatch):
    from tests import audit_method3_basis_metamorphism as audit

    root, config, _source, _state = _fixture(tmp_path)
    parent = _write(tmp_path / "experiment_data/parent.cif", "parent-input")
    manifest = {"parents": [{"parent_cif": parent.name, "cases": []}]}
    manifest_path = _write(tmp_path / "manifest.json", json.dumps(manifest))
    snapshot = FileSnapshot.capture(manifest_path)
    monkeypatch.setattr(audit, "PACKAGE_ROOT", root)
    monkeypatch.setattr(audit, "WORKSPACE", tmp_path)
    before = audit._provenance_state(manifest, snapshot)
    config.write_text(config.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    after = audit._provenance_state(manifest, snapshot)
    assert "runtime-config" in changed_provenance_surfaces(before, after)
    _write(root / "runtime/iso", "new-independent-program")
    assert "backend:iso-binary" in changed_provenance_surfaces(after, audit._provenance_state(manifest, snapshot))

@pytest.mark.parametrize("manifest", [
    {"status": "unfrozen_pending_independent_inputs", "parents": []},
    {"status": "frozen", "parents": []},
])
def test_method3_rejects_unfrozen_or_empty_scope_before_execution(tmp_path, monkeypatch, manifest):
    from tests import audit_method3_basis_metamorphism as audit

    path = _write(tmp_path / "manifest.json", json.dumps(manifest))
    monkeypatch.setattr("sys.argv", ["audit", "--manifest", str(path)])
    with pytest.raises(SystemExit) as failure:
        audit.main()
    assert failure.value.code == 2

@pytest.mark.parametrize("update_manifest", [False, True])
def test_method4_manifest_replacement_after_parse_cannot_receive_new_signature(tmp_path, monkeypatch, update_manifest):
    from tests import validate_method4_local as audit
    from tests.method4_provenance import build_method4_provenance_state

    root, _config, _source, _state = _fixture(tmp_path)
    manifest_path = _write(tmp_path / "manifest.json", json.dumps({
        "signature": "parsed-input", "parents": {}, "local_validation_tolerances": {},
    }))
    report_path = tmp_path / "report.json"
    original_read_text, original_read_bytes = Path.read_text, Path.read_bytes
    replaced = False

    def replace_after_read(path, value):
        nonlocal replaced
        if path.resolve() == manifest_path.resolve() and not replaced:
            replaced = True
            path.write_text(json.dumps({
                "signature": "replacement-input", "parents": {}, "local_validation_tolerances": {},
            }), encoding="utf-8")
        return value

    monkeypatch.setattr(Path, "read_text", lambda path, *a, **kw: replace_after_read(path, original_read_text(path, *a, **kw)))
    monkeypatch.setattr(Path, "read_bytes", lambda path: replace_after_read(path, original_read_bytes(path)))
    monkeypatch.setattr(audit, "build_method4_provenance_state", lambda source, selection: build_method4_provenance_state(
        source, selection, workspace=tmp_path, package_root=root,
    ))
    monkeypatch.setattr(audit, "WORKSPACE", tmp_path)
    monkeypatch.setattr("sys.argv", ["audit", "--manifest", str(manifest_path), "--json-output", str(report_path)]
                        + (["--update-manifest-status"] if update_manifest else []))
    assert audit.main() == 1
    report = json.loads(original_read_text(report_path, encoding="utf-8"))
    assert report["provenance_failure_count"] == 1
    assert "input-manifest" in report["provenance"]["changed_surfaces_during_run"]
    assert json.loads(original_read_text(manifest_path, encoding="utf-8"))["signature"] == "replacement-input"
    assert report["provenance"]["manifest_self_update"]["performed"] is False

@pytest.mark.parametrize("stage", ["a", "b"])
def test_stage_signatures_use_the_consumed_json_snapshots_after_disk_replacement(tmp_path, stage):
    from tests import audit_method3_stage_a_diagnostics as stage_a
    from tests import audit_method3_stage_b_feasibility as stage_b

    manifest = _write(tmp_path / "manifest.json", '{"parents": []}')
    route = _write(tmp_path / "routes.json", '{"cases": ["consumed"]}')
    previous = _write(tmp_path / "stage-a.json", '{"cases": ["consumed"]}')
    captured = tuple(FileSnapshot.capture(path) for path in (manifest, route, previous))
    for path in (manifest, route, previous):
        path.write_text('{"parents": [], "cases": ["replacement"]}', encoding="utf-8")
    if stage == "a":
        state = stage_a._signature_state(
            manifest, route, tmp_path, manifest=json.loads(captured[0].text),
            runtime_versions={}, tolerances={}, captured_inputs=captured,
        )
        expected_labels = ("manifest", "official-route-audit")
    else:
        state = stage_b._signature_state(
            manifest_path=manifest, route_report_path=route, stage_a_report_path=previous,
            input_root=tmp_path, package_root=stage_b.PACKAGE_ROOT, tool_path=Path(stage_b.__file__),
            runtime_versions={}, tolerances={}, max_quotient_order=1, case_ids=None,
            captured_inputs=captured,
        )
        expected_labels = ("manifest", "official-route-report", "stage-a-all-report")
        assert stage_b._load_json(captured[1], "route")["cases"] == ["consumed"]
    records = {row["label"]: row for row in state["inputs"]}
    for label, snapshot in zip(expected_labels, captured):
        assert records[label]["sha256"] == snapshot.sha256
        assert records[label]["sha256"] != FileSnapshot.capture(snapshot.path).sha256

def test_method4_official_gate_rejects_a_parsed_manifest_signed_from_replacement_bytes(tmp_path, monkeypatch):
    from tests import audit_method4_official as audit
    from tests.method4_provenance import build_method4_provenance_state

    root, _config, _source, _state = _fixture(tmp_path)
    path = _write(tmp_path / "manifest.json", '{"parents": {}, "signature": "consumed"}')
    snapshot = FileSnapshot.capture(path)
    path.write_text('{"parents": {}, "signature": "replacement"}', encoding="utf-8")
    observed = []

    def sign(source, selection):
        observed.append(source)
        return build_method4_provenance_state(source, selection, workspace=tmp_path, package_root=root)

    monkeypatch.setattr(audit, "build_method4_provenance_state", sign)
    result = audit._local_report_provenance_gate({}, json.loads(snapshot.text), path)
    assert "parsed manifest does not match the captured signed content" in result["issues"]
    result = audit._local_report_provenance_gate({}, json.loads(snapshot.text), snapshot)
    assert "parsed manifest does not match the captured signed content" not in result["issues"]
    assert observed[-1] is snapshot
