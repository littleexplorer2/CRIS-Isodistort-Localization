"""Tests for test_method4_provenance; migrated without changing scientific expectations."""
from __future__ import annotations

from pathlib import Path
from tests.method4_provenance import build_method4_provenance_state
import json


# Original coverage: test_method4_official_audit.py.

def test_method4_provenance_covers_code_config_runtime_backend_and_inputs(
    tmp_path: Path,
) -> None:
    package_root = tmp_path / "ISODISTORT"
    files = {
        "tests/validate_method4_local.py": "validator-v1",
        "tests/prepare_method4_inputs.py": "preparer-v1",
        "application/session/core_api.py": "core-v1",
        "backend/method4/strain_modes.py": "strain-modes-v1",
        "pyproject.toml": "[project]\nname='fixture'\n",
        "requirements.txt": "numpy\nscipy\n",
        "runtime/iso": "iso-v1",
        "runtime/findsym": "findsym-v1",
        "runtime/smodes": "smodes-v1",
        "runtime/const.dat": "const-v1",
        "runtime/data_isotropy.txt": "data-v1",
    }
    for relative, content in files.items():
        path = package_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    config_path = package_root / "resources" / "config" / "settings.yaml"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(
        json.dumps({
            "isobyu": {
                "bin_dir": "../../runtime",
                "data_dir": "../../runtime",
                "iso_bin": "iso",
                "findsym_bin": "findsym",
                "smodes_bin": "smodes",
            }
        }),
        encoding="utf-8",
    )
    input_root = tmp_path / "inputs"
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("parent-v1", encoding="utf-8")
    (input_root / "Daughter.cif").write_text("daughter-v1", encoding="utf-8")
    manifest_path = package_root / "docs" / "manifests" / "method4.json"
    manifest_path.parent.mkdir(parents=True)
    manifest = {
        "signature": "manifest-source-v1",
        "parents": {
            "fixture": {
                "parent_cif": "inputs/Parent.cif",
                "cases": {
                    "C01": {"daughter_cif": "inputs/Daughter.cif"},
                },
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    selection = {"parent": None, "context": None, "case_id": None}

    def signature_state() -> dict:
        return build_method4_provenance_state(
            manifest_path,
            selection,
            workspace=tmp_path,
            package_root=package_root,
        )

    state = signature_state()
    labels = {record["label"] for record in state["inputs"]}
    assert {
        "runtime-config",
        "input-manifest",
        "production:application/session/core_api.py",
        "production:backend/method4/strain_modes.py",
        "backend:iso-binary",
        "backend:findsym-binary",
        "backend:smodes-binary",
        "backend:iso-data:const.dat",
        "backend:iso-data:data_isotropy.txt",
        "parent:fixture",
        "daughter:fixture:C01",
    } <= labels
    assert state["runtime"]["python"]["executable"]["sha256"]
    assert set(state["runtime"]["distributions"]) == {
        "numpy", "scipy", "pymatgen", "spglib", "PyYAML",
    }

    signatures = [state["signature"]]
    mutation_path = package_root / "backend" / "method4" / "strain_modes.py"
    mutation_path.write_text("strain-modes-v2", encoding="utf-8")
    signatures.append(signature_state()["signature"])
    config_path.write_text(config_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    signatures.append(signature_state()["signature"])
    manifest["marker"] = "manifest-v2"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    signatures.append(signature_state()["signature"])
    (package_root / "runtime" / "findsym").unlink()
    signatures.append(signature_state()["signature"])

    assert len(set(signatures)) == len(signatures)
