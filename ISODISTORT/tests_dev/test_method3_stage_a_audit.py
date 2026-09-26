from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests_dev.manual import audit_method3_stage_a_diagnostics as stage_a_module


def _tolerances() -> dict[str, object]:
    return {
        "configured": {
            "symmetry_cartesian_tolerance_angstrom": 1.25e-3,
            "symmetry_angle_tolerance_degrees": 4.5,
            "affine_exact_cartesian_tolerance_angstrom": 2.5e-7,
            "fractional_coordinate_tolerance": 3.5e-8,
        },
        "semantics": {"test": "fixed scientific inputs"},
    }


def _manifest(parent_names: tuple[str, ...] = ("Parent.cif",)) -> dict:
    return {
        "parents": [
            {
                "parent_cif": parent_name,
                "parent_space_group_type": 1,
                "cases": [
                    {
                        "id": f"M3-T-{index:02d}",
                        "space_group_type": 1,
                        "supercell_basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                        "direct_sublattice_centering": "d",
                        "lattice_type": "direct",
                        "distortion_types": ["strain", "displacive"],
                    }
                ],
            }
            for index, parent_name in enumerate(parent_names, start=1)
        ]
    }


def _route_audit(case_id: str = "M3-T-01") -> dict:
    return {
        "cases": [
            {
                "parent_cif": "Parent.cif",
                "case_id": case_id,
                "official_source": "result_table_html",
                "official_completeness_proven": True,
                "embeddings": [
                    {
                        "classification": "coupled_ir_required",
                        "official_identity": {
                            "space_group_number": 1,
                            "symbol": "P1",
                            "basis": [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]],
                            "origin": ["0", "0", "0"],
                            "s": 1,
                            "i": 1,
                        },
                    }
                ],
            }
        ]
    }


def test_stage_a_signature_covers_sources_config_parents_runtime_and_tolerances(
    tmp_path: Path,
) -> None:
    package_root = tmp_path / "ISODISTORT"
    input_root = tmp_path / "experiment_data"
    input_root.mkdir()
    manifest_path = tmp_path / "manifest.json"
    route_path = tmp_path / "route.json"
    manifest = _manifest(("Eu parent.cif", "Nd parent.cif"))
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    route_path.write_text("{}", encoding="utf-8")
    for name in ("Eu parent.cif", "Nd parent.cif"):
        (input_root / name).write_text(f"data_{name}", encoding="utf-8")

    tool_path = package_root / "tests_dev" / "manual" / "stage_a.py"
    relevant = {
        tool_path: "stage-a-v1",
        package_root / "config" / "settings.yaml": "config-v1",
        package_root / "isocore" / "utils" / "lattice.py": "lattice-v1",
        package_root / "isocore" / "structure" / "cif_io.py": "cif-reader-v1",
        package_root
        / "tests_dev"
        / "manual"
        / "method3_affine_equivalence.py": "affine-v1",
    }
    for path, content in relevant.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    runtime = {
        "python": {"version": [3, 12, 5]},
        "distributions": {
            "numpy": "2.0",
            "spglib": "2.5",
            "pymatgen": "2026.1",
        },
    }

    def signature(
        *,
        runtime_state: dict | None = None,
        tolerance_state: dict | None = None,
    ) -> dict:
        return stage_a_module._signature_state(
            manifest_path,
            route_path,
            input_root,
            manifest=manifest,
            package_root=package_root,
            tool_path=tool_path,
            runtime_versions=runtime_state or runtime,
            tolerances=tolerance_state or _tolerances(),
        )

    initial = signature()
    labels = {item["label"] for item in initial["inputs"]}
    assert {
        "config",
        "source:lattice",
        "source:cif-reader",
        "source:affine-comparator",
        "parent-cif:Eu parent.cif",
        "parent-cif:Nd parent.cif",
    } <= labels

    signatures = [initial["signature"]]
    mutations = [
        (package_root / "config" / "settings.yaml", "config-v2"),
        (package_root / "isocore" / "utils" / "lattice.py", "lattice-v2"),
        (package_root / "isocore" / "structure" / "cif_io.py", "cif-reader-v2"),
        (input_root / "Eu parent.cif", "eu-v2"),
        (input_root / "Nd parent.cif", "nd-v2"),
    ]
    for path, content in mutations:
        path.write_text(content, encoding="utf-8")
        signatures.append(signature()["signature"])

    changed_runtime = json.loads(json.dumps(runtime))
    changed_runtime["distributions"]["spglib"] = "2.6"
    signatures.append(signature(runtime_state=changed_runtime)["signature"])
    changed_tolerances = _tolerances()
    changed_tolerances["configured"][
        "fractional_coordinate_tolerance"
    ] = 9.5e-8
    signatures.append(signature(tolerance_state=changed_tolerances)["signature"])

    assert len(set(signatures)) == len(signatures)


def test_stage_a_runner_resumes_atomically_and_rejects_signature_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest_path = tmp_path / "manifest.json"
    route_path = tmp_path / "route.json"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    input_root.mkdir()
    manifest = _manifest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    route_path.write_text(json.dumps(_route_audit()), encoding="utf-8")
    (input_root / "Parent.cif").write_text("data_parent\n", encoding="utf-8")

    captured: dict[str, object] = {"search_count": 0}
    parent_structure = SimpleNamespace(name="parent")
    parent_model = SimpleNamespace(name="parent-model")

    class FakeRotation:
        def tolist(self) -> list[list[int]]:
            return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    class FakeAnalyzer:
        def __init__(
            self,
            structure,
            *,
            symprec: float,
            angle_tolerance: float,
        ) -> None:
            assert structure is parent_structure
            captured["analyzer_tolerances"] = (symprec, angle_tolerance)

        def get_symmetry_dataset(self):
            return SimpleNamespace(rotations=[FakeRotation()])

    official_identity = stage_a_module._identity_from_json(
        _route_audit()["cases"][0]["embeddings"][0]["official_identity"]
    )
    stage_item = SimpleNamespace(
        route_resolution="affine_only_unresolved_coupled_route",
        routes=[],
    )

    class FakeSearchEngine:
        def __init__(self, _backend) -> None:
            pass

        def method_3_search(self, parent_sg: int, query):
            assert parent_sg == 1
            captured["search_count"] = int(captured["search_count"]) + 1
            captured["query"] = query
            return [stage_item]

    monkeypatch.setattr(stage_a_module, "_scientific_tolerances", _tolerances)
    monkeypatch.setattr(stage_a_module, "read_cif", lambda _path: parent_structure)
    monkeypatch.setattr(stage_a_module, "SpacegroupAnalyzer", FakeAnalyzer)
    monkeypatch.setattr(
        stage_a_module,
        "parent_operations",
        lambda structure, symprec, angle: (
            parent_model
            if structure is parent_structure and symprec > 0 and angle > 0
            else None
        ),
    )
    monkeypatch.setattr(stage_a_module, "IsoSearchEngine", FakeSearchEngine)
    monkeypatch.setattr(
        stage_a_module,
        "_identity_from_local",
        lambda _item: official_identity,
    )
    monkeypatch.setattr(
        stage_a_module,
        "embedding_operations",
        lambda identity, _parent: identity,
    )
    monkeypatch.setattr(
        stage_a_module,
        "affine_equivalence_witness",
        lambda _left, _right, _parent: {"relationship": "same_affine_subgroup"},
    )

    first = stage_a_module.run_audit(
        manifest_path=manifest_path,
        route_audit_path=route_path,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        input_root=input_root,
    )
    configured = _tolerances()["configured"]
    query = captured["query"]
    assert first["summary"]["complete"] is True
    assert first["summary"]["resumed_case_count"] == 0
    assert captured["search_count"] == 1
    assert captured["analyzer_tolerances"] == (
        configured["symmetry_cartesian_tolerance_angstrom"],
        configured["symmetry_angle_tolerance_degrees"],
    )
    assert query.symmetry_tolerance == configured[
        "symmetry_cartesian_tolerance_angstrom"
    ]
    assert query.symmetry_angle_tolerance_degrees == configured[
        "symmetry_angle_tolerance_degrees"
    ]
    assert query.affine_exact_tolerance == configured[
        "affine_exact_cartesian_tolerance_angstrom"
    ]
    assert query.fractional_coordinate_tolerance == configured[
        "fractional_coordinate_tolerance"
    ]
    assert checkpoint_path.is_file()
    assert report_path.is_file()
    assert not checkpoint_path.with_suffix(".json.tmp").exists()
    assert not report_path.with_suffix(".json.tmp").exists()

    class BackendMustNotRun:
        def __init__(self, _backend) -> None:
            raise AssertionError("a compatible checkpoint must be reused")

    monkeypatch.setattr(stage_a_module, "IsoSearchEngine", BackendMustNotRun)
    second = stage_a_module.run_audit(
        manifest_path=manifest_path,
        route_audit_path=route_path,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        input_root=input_root,
    )
    assert second["summary"]["complete"] is True
    assert second["summary"]["resumed_case_count"] == 1
    assert captured["search_count"] == 1

    manifest["scientific_input_revision"] = 2
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(RuntimeError, match="signature mismatch"):
        stage_a_module.run_audit(
            manifest_path=manifest_path,
            route_audit_path=route_path,
            report_path=report_path,
            checkpoint_path=checkpoint_path,
            input_root=input_root,
        )
