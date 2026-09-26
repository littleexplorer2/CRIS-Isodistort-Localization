from __future__ import annotations

from pathlib import Path

from tests_dev.manual import audit_method3_stage_b_feasibility as audit


def test_stage_b_partial_audit_rejects_three_novel_orbits_and_resumes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    report_path = tmp_path / "stage_b_report.json"
    checkpoint_path = tmp_path / "stage_b_checkpoint.json"
    selected = ("M3-ND-12", "M3-ND-13")

    report = audit.run_audit(
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        case_ids=selected,
        restart=True,
    )

    assert report["summary"]["complete"]
    assert report["summary"]["official_embedding_count"] == 5
    assert report["summary"]["official_reachable_count"] == 5
    assert report["summary"]["novel_orbit_count"] == 3
    assert report["summary"]["novel_orbit_infeasible_count"] == 3
    assert report["summary"]["novel_orbits_by_case"] == {
        "M3-ND-12": 1,
        "M3-ND-13": 2,
    }
    assert report_path.is_file()
    assert checkpoint_path.is_file()
    assert not report_path.with_suffix(".json.tmp").exists()
    assert all(
        case["cache_contract"]["quotient_build_count"] == 1
        and case["cache_contract"]["representation_build_count"] == 1
        and case["cache_contract"]["feasibility_analyzer_build_count"] == 1
        for case in report["cases"]
    )
    assert all(
        not candidate["selectable"]
        and candidate["route_status"] == "unresolved_irrep_decomposition"
        for case in report["cases"]
        for candidate in case["candidate_results"]
    )
    signed_labels = {
        item["label"] for item in report["signature_state"]["inputs"]
    }
    assert {
        "manifest",
        "official-route-report",
        "stage-a-all-report",
        "source:stage-b-audit",
        "source:inverse-landau",
        "source:inverse-landau-adapter",
        "parent-cif:NdNiO2 own.cif",
    }.issubset(signed_labels)

    def _must_not_recompute(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("compatible checkpoint should resume both cases")

    monkeypatch.setattr(audit, "_audit_case", _must_not_recompute)
    resumed = audit.run_audit(
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        case_ids=selected,
        restart=False,
    )

    assert resumed["summary"]["complete"]
    assert resumed["checkpoint"]["resumed_case_count"] == 2
