from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from tests_dev.manual import audit_method3_embedding_routes as route_module
from tests_dev.manual.audit_method3_downloads import parse_subgroup_identity
from tests_dev.manual.audit_method3_embedding_routes import (
    audit_embedding,
    build_display_direction_session,
    classify_embedding_routes,
    parse_display_direction_output,
    run_audit,
)


def _identity():
    return parse_subgroup_identity(
        "12 C2/m, basis={(1,-1,0),(0,0,1),(-1,0,0)}, "
        "origin=(0,0,0), s=1, i=4"
    )


def test_build_display_direction_session_preserves_exact_fractions():
    identity = parse_subgroup_identity(
        "12 C2/m, basis={(-1,-1,0),(1,-1,0),(1/2,1/2,1/2)}, "
        "origin=(-1/4,0,1/2), s=1, i=4"
    )

    commands = build_display_direction_session(139, identity)

    assert "VALUE PARENT 139" in commands
    assert "VALUE SUBGROUP 12" in commands
    assert "VALUE BASIS -1,-1,0 1,-1,0 1/2,1/2,1/2" in commands
    assert "VALUE ORIGIN -1/4,0,1/2" in commands
    assert commands[-2:] == ["DISPLAY DIRECTION", "QUIT"]


def test_parse_and_classify_exact_single_ir_route():
    stdout = """Isotropy, Version 9.6.1, Jan 2022
**********Irrep (ML) k vector Dir    Subgroup   Size
GM1+       (0,0,0)  (a)    139 I4/mmm 1
GM5+       (0,0,0)  (a,-a) 12 C2/m    1
*
"""

    routes, errors = parse_display_direction_output(stdout)
    status, exact = classify_embedding_routes(_identity(), routes, errors)

    assert errors == []
    assert status == "single_ir_exact"
    assert [route["irrep"] for route in exact] == ["GM5+"]


def test_same_space_group_with_wrong_size_is_not_exact():
    stdout = """**********Irrep (ML) k vector    Dir       Subgroup Size
GM5+       (0,0,0)     (a,b)     2 P-1    1
N1+        (1/2,0,1/2) (a,0,b,0) 2 P-1    4
*
"""
    identity = parse_subgroup_identity(
        "2 P-1, basis={(0,2,0),(-2,0,0),(1,-1,1)}, "
        "origin=(-1/4,-1/4,1/4), s=8, i=64"
    )

    routes, errors = parse_display_direction_output(stdout)
    status, exact = classify_embedding_routes(identity, routes, errors)

    assert status == "coupled_ir_required"
    assert exact == []


def test_parse_parametric_k_row_with_k_parameter_column():
    stdout = """**********Irrep (ML) k params k vector  Dir   Subgroup   Size
GM1+                (0,0,0)   (a)   139 I4/mmm 1
LD1        1/12     (0,0,1/6) (a,0) 123 P4/mmm 12
*
"""
    identity = parse_subgroup_identity(
        "123 P4/mmm, basis={(1,0,0),(0,1,0),(0,0,6)}, "
        "origin=(0,0,0), s=12, i=12"
    )

    routes, errors = parse_display_direction_output(stdout)
    status, exact = classify_embedding_routes(identity, routes, errors)

    assert status == "single_ir_exact"
    assert exact[0]["irrep"] == "LD1"
    assert exact[0]["k_parameters"] == "1/12"
    assert exact[0]["k_vector"] == "(0,0,1/6)"


def test_missing_direction_table_is_indeterminate_not_coupled():
    routes, errors = parse_display_direction_output(
        "DISPLAY DIRECTION has not been implemented for superspacegroups\n"
    )
    status, exact = classify_embedding_routes(_identity(), routes, errors)

    assert status == "indeterminate"
    assert exact == []
    assert errors


def test_audit_embedding_records_commands_and_raw_provenance():
    stdout = """Isotropy, Version 9.6.1, Jan 2022
**********Irrep (ML) k vector Dir Subgroup Size
GM5+ (0,0,0) (a,-a) 12 C2/m 1
*
"""
    seen: list[str] = []

    def run_iso(commands: str) -> str:
        seen.append(commands)
        return stdout

    record = audit_embedding(139, _identity(), run_iso)

    assert record["classification"] == "single_ir_exact"
    assert record["raw_stdout"] == stdout
    assert record["stdout_sha256"]
    assert record["iso_input"]["target_size"] == 1
    assert record["iso_input"]["origin"] == ["0", "0", "0"]
    assert "DISPLAY DIRECTION\nQUIT\n" in seen[0]
    assert Fraction(record["iso_input"]["basis"][0][0]) == 1


def test_candidate_inventory_route_audit_is_explicitly_provisional(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = {
        "signature": "test",
        "parents": [
            {
                "parent_cif": "Parent.cif",
                "parent_space_group_type": 139,
                "cases": [
                    {
                        "id": "M3-T-01",
                        "readable_case_folder": "case-one",
                        "space_group_type": 12,
                        "supercell_basis": [[1, -1, 0], [0, 0, 1], [-1, 0, 0]],
                        "direct_sublattice_centering": "d",
                        "distortion_types": ["strain", "displacive"],
                    }
                ],
            }
        ],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    output_root = tmp_path / "output_compare"
    candidate_dir = (
        output_root
        / "Parent.cif"
        / "官网"
        / "Method3"
        / "case-one"
        / "candidate"
    )
    candidate_dir.mkdir(parents=True)
    (candidate_dir / "subgroup.cif").write_text(
        "# Subgroup: 12 C2/m, basis={(1,-1,0),(0,0,1),(-1,0,0)}, "
        "origin=(0,0,0), s=1, i=4\n",
        encoding="utf-8",
    )
    runtime_dir = tmp_path / "runtime"
    runtime_dir.mkdir()
    binary = runtime_dir / "iso"
    binary.write_bytes(b"fake-iso")

    class FakeIsoWrapper:
        def __init__(self) -> None:
            self.binary = binary

        def run_stdin(self, _binary, commands, timeout=None):
            if commands == "QUIT\n":
                return "Isotropy, Version test\n"
            return (
                "**********Irrep (ML) k vector Dir Subgroup Size\n"
                "GM5+ (0,0,0) (a,-a) 12 C2/m 1\n*\n"
            )

    monkeypatch.setattr(route_module, "IsoWrapper", FakeIsoWrapper)
    report = run_audit(
        manifest_path,
        output_root,
        include_raw_output=False,
        allow_candidate_inventory=True,
    )

    assert report["summary"]["audited_case_count"] == 0
    assert report["summary"]["provisional_audited_case_count"] == 1
    assert report["summary"]["official_embedding_count"] == 0
    assert report["summary"]["provisional_inventory_embedding_count"] == 1
    assert report["summary"]["audited_embedding_count"] == 1
    case = report["cases"][0]
    assert case["status"] == "audited_provisional"
    assert case["official_source"] == "downloaded_candidate_inventory"
    assert case["official_completeness_proven"] is False
    assert case["embeddings"][0]["classification"] == "single_ir_exact"


def _route_manifest() -> dict:
    return {
        "signature": "test",
        "parents": [
            {
                "parent_cif": "Parent.cif",
                "parent_space_group_type": 139,
                "cases": [
                    {
                        "id": "M3-T-01",
                        "readable_case_folder": "case-one",
                        "space_group_type": 12,
                        "supercell_basis": [[1, -1, 0], [0, 0, 1], [-1, 0, 0]],
                        "direct_sublattice_centering": "d",
                        "distortion_types": ["strain", "displacive"],
                    }
                ],
            }
        ],
    }


def _write_route_result(path: Path, *, target_sg: int = 12) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "<html><body>"
        "Space Group: 139 I4/mmm D4h-17<br>"
        "Include strain, displacive ALL distortions<br>"
        f"Subgroup: {target_sg} C2/m, basis of lattice = "
        "(1,-1,0),(0,0,1),(-1,0,0), C centering<br>"
        '<input type="hidden" name="spacegroup" value="139 I4/mmm D4h-17">'
        '<input type="hidden" name="includestrain" value="true">'
        '<input type="hidden" name="includedisplacive001" value="true">'
        "Finish selecting the distortion mode by selecting one of the following choices:"
        '<input type="radio" name="orderparam" value="x">'
        "12 C2/m, basis={(1,-1,0),(0,0,1),(-1,0,0)}, "
        "origin=(0,0,0), s=1, i=4<br>"
        "</body></html>",
        encoding="utf-8",
    )


def _install_fake_wrapper(tmp_path: Path, monkeypatch, calls: dict[str, int]) -> None:
    runtime_dir = tmp_path / "runtime"
    runtime_dir.mkdir()
    binary = runtime_dir / "iso"
    binary.write_bytes(b"fake-iso")

    class FakeIsoWrapper:
        def __init__(self) -> None:
            self.binary = binary

        def run_stdin(self, _binary, commands, timeout=None):
            if commands == "QUIT\n":
                return "Isotropy, Version test\n"
            calls["embedding"] += 1
            raise AssertionError("invalid official evidence must not reach DISPLAY DIRECTION")

    monkeypatch.setattr(route_module, "IsoWrapper", FakeIsoWrapper)


def test_route_audit_rejects_result_page_bound_to_wrong_manifest_context(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(_route_manifest()), encoding="utf-8")
    output_root = tmp_path / "output_compare"
    case_dir = output_root / "Parent.cif" / "官网" / "Method3" / "case-one"
    _write_route_result(case_dir / "result.html", target_sg=8)
    calls = {"embedding": 0}
    _install_fake_wrapper(tmp_path, monkeypatch, calls)

    report = run_audit(manifest_path, output_root, include_raw_output=False)

    assert report["summary"]["error_case_count"] == 1
    assert report["cases"][0]["status"] == "error"
    assert "manifest query context" in report["cases"][0]["issues"][0]
    assert calls["embedding"] == 0


def test_route_audit_rejects_invalid_candidate_cif_fallback(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(_route_manifest()), encoding="utf-8")
    output_root = tmp_path / "output_compare"
    candidate = (
        output_root
        / "Parent.cif"
        / "官网"
        / "Method3"
        / "case-one"
        / "candidate"
    )
    candidate.mkdir(parents=True)
    (candidate / "subgroup.cif").write_text("data_invalid\n", encoding="utf-8")
    calls = {"embedding": 0}
    _install_fake_wrapper(tmp_path, monkeypatch, calls)

    report = run_audit(
        manifest_path,
        output_root,
        include_raw_output=False,
        allow_candidate_inventory=True,
    )

    case = report["cases"][0]
    assert report["summary"]["skipped_case_count"] == 1
    assert case["status"] == "skipped_no_unique_official_result"
    assert any("invalid subgroup.cif" in issue for issue in case["issues"])
    assert calls["embedding"] == 0


def test_route_audit_rejects_missing_candidate_cif_fallback(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(_route_manifest()), encoding="utf-8")
    output_root = tmp_path / "output_compare"
    case_dir = output_root / "Parent.cif" / "官网" / "Method3" / "case-one"
    case_dir.mkdir(parents=True)
    calls = {"embedding": 0}
    _install_fake_wrapper(tmp_path, monkeypatch, calls)

    report = run_audit(
        manifest_path,
        output_root,
        include_raw_output=False,
        allow_candidate_inventory=True,
    )

    case = report["cases"][0]
    assert report["summary"]["skipped_case_count"] == 1
    assert case["status"] == "skipped_no_unique_official_result"
    assert case["inventory_candidate_count"] == 0
    assert calls["embedding"] == 0


def test_route_audit_exit_code_only_accepts_clean_provisional_evidence() -> None:
    summary = {
        "selected_case_count": 1,
        "error_case_count": 0,
        "skipped_case_count": 0,
        "provisional_audited_case_count": 1,
    }
    report = {"summary": summary}

    assert route_module._audit_exit_code(report) == 2
    assert route_module._audit_exit_code(report, accept_provisional=True) == 0

    summary["skipped_case_count"] = 1
    assert route_module._audit_exit_code(report, accept_provisional=True) == 2
    summary["skipped_case_count"] = 0
    summary["error_case_count"] = 1
    assert route_module._audit_exit_code(report, accept_provisional=True) == 1


def test_route_cli_accept_provisional_flag_is_narrow(
    tmp_path: Path,
    monkeypatch,
) -> None:
    report = {
        "summary": {
            "selected_case_count": 1,
            "error_case_count": 0,
            "skipped_case_count": 0,
            "provisional_audited_case_count": 1,
        },
        "cases": [],
    }
    monkeypatch.setattr(route_module, "run_audit", lambda *_args, **_kwargs: report)
    monkeypatch.setattr(route_module, "_print_summary", lambda _report: None)
    output = tmp_path / "report.json"

    assert route_module.main(["--json-output", str(output)]) == 2
    assert (
        route_module.main(["--json-output", str(output), "--accept-provisional"])
        == 0
    )

    report["summary"]["skipped_case_count"] = 1
    assert (
        route_module.main(["--json-output", str(output), "--accept-provisional"])
        == 2
    )
