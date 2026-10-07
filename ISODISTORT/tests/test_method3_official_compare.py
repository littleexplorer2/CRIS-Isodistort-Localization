from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import pytest
from pymatgen.core import Lattice, Structure

from backend.utils.lattice import identity_matrix, rational_matrix
from backend.wrappers import SubgroupInfo
from features.method1.affine_embeddings import AffineEmbedding, AffineOperation
from features.method1.search_methods import Method3ResultItem
from features.method3.coupled_routes import CoupledRouteWitness
from tests.manual import compare_method3_official_local as comparison_module
from tests.manual.audit_method3_downloads import parse_subgroup_identity
from tests.manual.compare_method3_official_local import (
    compare_embedding_sets,
    run_comparison,
    select_manifest_cases,
)
from tests.manual.method3_affine_equivalence import (
    compare_affine_embedding_sets,
)


def _local_item(
    *,
    basis: list[list[int]] | None = None,
    origin: list[int] | None = None,
) -> SimpleNamespace:
    subgroup = SimpleNamespace(
        space_group_number=139,
        space_group_symbol="I4/mmm",
        basis_vectors=basis or [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=origin or [0, 0, 0],
        size=1,
        subgroup_index=1,
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        irrep_label="GM1+",
        opd_symbol="P1",
        k_point_label="GM",
        k_parameters=[],
        k_active_raw="(0,0,0)",
    )
    return SimpleNamespace(subgroup=subgroup, basis=subgroup.basis_vectors, routes=[subgroup])


def _exact_fixed_space_item() -> Method3ResultItem:
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        subgroup_index=1,
        size=1,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
        parent_sg=139,
    )
    first_component = SubgroupInfo(
        index=1,
        space_group_number=123,
        irrep_label="GM2+",
        opd_symbol="P1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    second_component = SubgroupInfo(
        index=2,
        space_group_number=123,
        irrep_label="GM4+",
        opd_symbol="P1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    zero = (Fraction(0), Fraction(0), Fraction(0))
    embedding = AffineEmbedding(
        identity_matrix(),
        (AffineOperation(rational_matrix(identity_matrix()), zero),),
    )
    embedding_id = "m3e1-" + "a" * 64
    subgroup._method3_embedding_id = embedding_id
    subgroup._method3_route_resolution = "exact_fixed_space"
    return Method3ResultItem(
        subgroup=subgroup,
        point_group="4/mmm",
        basis=subgroup.basis_vectors,
        routes=[],
        affine_embedding=embedding,
        route_resolution="exact_fixed_space",
        embedding_id=embedding_id,
        coupled_witness=CoupledRouteWitness(
            components=(first_component, second_component),
            target_subgroup_order=1,
            quotient_order=4,
            explored_intersection_count=3,
        ),
        feasibility_status="embedding_feasible",
        fixed_space_dimension=2,
    )


def _manifest() -> dict:
    return {
        "official_download_root": "output_compare/<parent>/官网/Method3/<readable-case-folder>/",
        "parents": [
            {
                "parent_cif": "Parent.cif",
                "parent_space_group_type": 139,
                "cases": [
                    {
                        "id": "M3-T-01",
                        "readable_case_folder": "case-one",
                        "space_group_type": 139,
                        "supercell_basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                        "direct_sublattice_centering": "d",
                        "lattice_type": "direct",
                        "distortion_types": ["strain", "displacive"],
                        "category": "identity",
                    },
                    {
                        "id": "M3-T-02",
                        "readable_case_folder": "case-two",
                        "space_group_type": 71,
                        "supercell_basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                        "direct_sublattice_centering": "d",
                        "lattice_type": "direct",
                        "distortion_types": ["displacive"],
                        "category": "identity",
                    },
                ],
            }
        ],
    }


def _write_result_page(path: Path, *, parent_sg: int = 139) -> None:
    path.parent.mkdir(parents=True)
    path.write_text(
        "<html><body>"
        f"Space Group: {parent_sg} I4/mmm D4h-17<br>"
        "Include strain, displacive ALL distortions<br>"
        "Subgroup: 139 I4/mmm, basis of lattice = (1,0,0),(0,1,0),(0,0,1), I centering<br>"
        f'<input type="hidden" name="spacegroup" value="{parent_sg} I4/mmm D4h-17">'
        '<input type="hidden" name="includestrain" value="true">'
        '<input type="hidden" name="includedisplacive001" value="true">'
        "Finish selecting the distortion mode by selecting one of the following choices:"
        '<input type="radio" name="orderparam" value="x">'
        "139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=(0,0,0), s=1, i=1<br>"
        "</body></html>",
        encoding="utf-8",
    )


def test_exact_comparison_does_not_conflate_alternative_basis() -> None:
    official = [
        parse_subgroup_identity(
            "139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=(0,0,0), s=1, i=1"
        )
    ]
    comparison = compare_embedding_sets(
        official,
        [_local_item(basis=[[0, 1, 0], [1, 0, 0], [0, 0, 1]])],
    )

    assert comparison["count_match"] is True
    assert comparison["space_group_multiset_match"] is True
    assert comparison["exact_representation_match"] is False
    assert len(comparison["missing_exact_official_embeddings"]) == 1
    assert len(comparison["unexpected_exact_local_embeddings"]) == 1
    assert comparison["local_embeddings"][0]["raw_fields"]["basis"] == [
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 1],
    ]


def test_embedding_match_is_inconclusive_when_second_stage_route_is_unresolved() -> None:
    official = [
        parse_subgroup_identity(
            "139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1"
        )
    ]
    item = _local_item()
    item.routes = []
    item.route_resolution = "affine_only_unresolved_coupled_route"
    item.feasibility_status = "embedding_feasible"
    item.fixed_space_dimension = 2

    comparison = compare_embedding_sets(official, [item])

    assert comparison["embedding_status"] == "exact_match"
    assert comparison["status"] == "inconclusive_route_resolution"
    assert comparison["route_ready"] is False
    assert comparison["route_readiness"][0]["issues"] == [
        "empty routes require route_resolution=exact_fixed_space",
        "exact coupled route requires a witness with at least two components",
        "exact coupled route requires an exact affine stabilizer-intersection witness",
    ]


def test_exact_coupled_fixed_space_is_route_ready() -> None:
    official = [
        parse_subgroup_identity(
            "139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1"
        )
    ]
    item = _exact_fixed_space_item()

    comparison = compare_embedding_sets(official, [item])

    assert comparison["status"] == "exact_match"
    assert comparison["route_ready"] is True
    readiness = comparison["route_readiness"][0]
    assert readiness["embedding_id"] == item.embedding_id
    assert readiness["subgroup_embedding_id"] == item.embedding_id
    assert readiness["subgroup_route_resolution"] == "exact_fixed_space"
    assert readiness["affine_operation_count"] == 1
    assert readiness["target_subgroup_order"] == 1
    assert readiness["quotient_order"] == 4
    assert readiness["explored_intersection_count"] == 3


def test_exact_coupled_fixed_space_rejects_inconsistent_embedding_provenance() -> None:
    official = [
        parse_subgroup_identity(
            "139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1"
        )
    ]
    item = _exact_fixed_space_item()
    item.embedding_id = "m3e1-" + "b" * 64
    item.affine_embedding = None
    item.subgroup._method3_route_resolution = "known_single_ir"

    comparison = compare_embedding_sets(official, [item])

    assert comparison["status"] == "inconclusive_route_resolution"
    issues = comparison["route_readiness"][0]["issues"]
    assert "exact coupled route requires an AffineEmbedding instance" in issues
    assert (
        "exact coupled route requires subgroup embedding marker to match embedding_id"
        in issues
    )
    assert (
        "exact coupled route requires subgroup route marker to match route_resolution"
        in issues
    )


@pytest.mark.parametrize(
    ("witness", "expected_issue"),
    [
        (
            CoupledRouteWitness(
                components=tuple(_exact_fixed_space_item().coupled_witness.components),
                target_subgroup_order=2,
                quotient_order=4,
                explored_intersection_count=3,
            ),
            "witness target_subgroup_order must equal affine embedding operation count",
        ),
        (
            replace(
                _exact_fixed_space_item().coupled_witness,
                quotient_order=1,
            ),
            "witness quotient_order must be a larger multiple of target_subgroup_order",
        ),
        (
            replace(
                _exact_fixed_space_item().coupled_witness,
                explored_intersection_count=2,
            ),
            "witness explored_intersection_count is too small for its component path",
        ),
        (
            CoupledRouteWitness(
                components=(
                    _exact_fixed_space_item().coupled_witness.components[0],
                    _exact_fixed_space_item().coupled_witness.components[0],
                ),
                target_subgroup_order=1,
                quotient_order=4,
                explored_intersection_count=3,
            ),
            "exact coupled route witness components must be distinct routes",
        ),
    ],
)
def test_exact_coupled_fixed_space_rejects_inconsistent_witness_counts(
    witness: CoupledRouteWitness,
    expected_issue: str,
) -> None:
    official = [
        parse_subgroup_identity(
            "139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1"
        )
    ]
    item = _exact_fixed_space_item()
    item.coupled_witness = witness

    comparison = compare_embedding_sets(official, [item])

    assert comparison["status"] == "inconclusive_route_resolution"
    assert expected_issue in comparison["route_readiness"][0]["issues"]


def _i4mmm_parent() -> Structure:
    return Structure.from_spacegroup(
        "I4/mmm",
        Lattice.tetragonal(4.0, 10.0),
        ["Eu"],
        [[0, 0, 0]],
    )


def test_affine_comparison_proves_alternative_basis_is_same_subgroup() -> None:
    official = parse_subgroup_identity(
        "71 Immm, basis={(1,0,0),(0,0,-1),(0,1,0)}, "
        "origin=(0,0,0), s=1, i=2"
    )
    local = parse_subgroup_identity(
        "71 Immm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
        "origin=(0,0,0), s=1, i=2"
    )

    result = compare_affine_embedding_sets(
        [official],
        [local],
        _i4mmm_parent(),
        symprec=1e-5,
        angle_tolerance_degrees=3.25,
    )

    assert result["complete_match"] is True
    assert result["matches"][0]["relationship"] == "same_affine_subgroup"
    assert result["parent"]["symmetry_cartesian_tolerance_angstrom"] == 1e-5
    assert result["parent"]["symmetry_angle_tolerance_degrees"] == 3.25


def test_affine_comparison_uses_parent_translation_quotient_for_origins() -> None:
    official = parse_subgroup_identity(
        "140 I4/mcm, basis={(-1,1,0),(-1,-1,0),(0,0,2)}, "
        "origin=(-1/2,-1/2,0), s=4, i=4"
    )
    local = parse_subgroup_identity(
        "140 I4/mcm, basis={(-1,1,0),(-1,-1,0),(0,0,2)}, "
        "origin=(0,0,1/2), s=4, i=4"
    )

    result = compare_affine_embedding_sets(
        [official], [local], _i4mmm_parent(), symprec=1e-5
    )

    assert result["complete_match"] is True
    assert result["matches"][0]["relationship"] == "parent_affine_conjugate"
    assert result["matches"][0]["translation_coset"] != ["0", "0", "0"]


def test_affine_comparison_preserves_missing_embedding_multiplicity() -> None:
    identity = parse_subgroup_identity(
        "71 Immm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
        "origin=(0,0,0), s=1, i=2"
    )

    result = compare_affine_embedding_sets(
        [identity, identity], [identity], _i4mmm_parent(), symprec=1e-5
    )

    assert result["complete_match"] is False
    assert result["match_count"] == 1
    assert len(result["unmatched_official_indices"]) == 1
    assert result["unmatched_local_indices"] == []


def test_manifest_filters_are_exact_and_case_insensitive() -> None:
    manifest = _manifest()
    selected = select_manifest_cases(manifest, parents=["parent.CIF"], case_ids=["m3-t-02"])
    assert [item["case"]["id"] for item in selected] == ["M3-T-02"]


def test_runner_skips_missing_result_page_and_resumes_without_backend(tmp_path: Path) -> None:
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(_manifest()), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    _write_result_page(output_root / "Parent.cif" / "官网" / "Method3" / "case-one" / "result.html")
    invalid_dir = output_root / "Parent.cif" / "官网" / "Method3" / "case-two"
    invalid_dir.mkdir(parents=True)
    (invalid_dir / "search.html").write_text("Method 3: Search over arbitrary k points", encoding="utf-8")

    calls = {"factory": 0, "load": 0, "search": 0}

    class FakeIsoDistort:
        def __init__(self) -> None:
            calls["factory"] += 1

        def load_structure(self, path: Path) -> None:
            assert path == input_root / "Parent.cif"
            calls["load"] += 1

        def search_method_3(self, **query):
            assert query["space_group_type"] == 139
            calls["search"] += 1
            return [_local_item()]

    first = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=FakeIsoDistort,
    )
    assert first["summary"] == {
        "complete": True,
        "validation_status": "complete",
        "source_signature_stable": True,
        "source_drift_detected": False,
        "selected_case_count": 2,
        "comparable_case_count": 1,
        "official_eligible_case_count": 1,
        "authoritative_result_table_case_count": 1,
        "provisional_candidate_inventory_case_count": 0,
        "skipped_no_unique_official_result_count": 1,
        "exact_match_count": 1,
        "affine_equivalent_count": 0,
        "inconclusive_route_resolution_count": 0,
        "difference_count": 0,
        "local_error_count": 0,
        "authoritative_exact_match_count": 1,
        "authoritative_affine_equivalent_count": 0,
        "authoritative_inconclusive_route_resolution_count": 0,
        "authoritative_difference_count": 0,
        "authoritative_local_error_count": 0,
        "provisional_exact_match_count": 0,
        "provisional_affine_equivalent_count": 0,
        "provisional_inconclusive_route_resolution_count": 0,
        "provisional_difference_count": 0,
        "provisional_local_error_count": 0,
        "resumed_case_count": 0,
    }
    assert first["source_signature"] == {
        "start": first["run_signature"],
        "end": first["run_signature"],
        "stable": True,
    }
    assert calls == {"factory": 1, "load": 1, "search": 1}
    assert json.loads(report_path.read_text(encoding="utf-8"))["summary"]["complete"] is True

    class BackendMustNotRun:
        def __init__(self) -> None:
            raise AssertionError("compatible checkpoint should have been reused")

    second = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=BackendMustNotRun,
    )
    assert second["summary"]["resumed_case_count"] == 1
    assert second["summary"]["exact_match_count"] == 1

    # Filtering the same manifest must keep the full-batch signature so a
    # one-case retry reuses the compatible checkpoint rather than replacing it.
    third = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        case_ids=["M3-T-01"],
        iso_factory=BackendMustNotRun,
    )
    assert third["summary"]["selected_case_count"] == 1
    assert third["summary"]["resumed_case_count"] == 1


def test_runner_fails_closed_and_invalidates_checkpoint_when_source_drifts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    _write_result_page(
        output_root / "Parent.cif" / "官网" / "Method3" / "case-one" / "result.html"
    )

    def manifest_signature(path: Path, *_args, **_kwargs) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    monkeypatch.setattr(comparison_module, "_run_signature", manifest_signature)
    calls = {"search": 0}
    mutate_during_search = {"enabled": True}

    class FakeIsoDistort:
        def load_structure(self, path: Path) -> None:
            assert path == input_root / "Parent.cif"

        def search_method_3(self, **query):
            assert query["space_group_type"] == 139
            calls["search"] += 1
            if mutate_during_search["enabled"]:
                manifest_path.write_text(
                    manifest_path.read_text(encoding="utf-8") + "\n",
                    encoding="utf-8",
                )
                mutate_during_search["enabled"] = False
            return [_local_item()]

    drifted = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=FakeIsoDistort,
    )

    assert drifted["source_signature"]["stable"] is False
    assert drifted["source_signature"]["start"] != drifted["source_signature"]["end"]
    assert drifted["summary"]["complete"] is False
    assert drifted["summary"]["validation_status"] == "source_drift"
    assert drifted["summary"]["source_drift_detected"] is True
    assert comparison_module._comparison_exit_code(drifted) == 1
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert checkpoint["signature"] == "invalidated:source_signature_drift"
    assert checkpoint["invalidation_reason"] == "source_signature_drift"
    assert checkpoint["source_signature"] == drifted["source_signature"]
    assert checkpoint["results"] == {}

    stable = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=FakeIsoDistort,
    )

    assert calls["search"] == 2
    assert stable["summary"]["complete"] is True
    assert stable["summary"]["source_signature_stable"] is True
    assert stable["summary"]["resumed_case_count"] == 0


def test_runner_consumes_the_same_manifest_and_html_snapshot_that_it_signs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    result_path = output_root / "Parent.cif" / "官网" / "Method3" / "case-one" / "result.html"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    _write_result_page(result_path)

    original_parser = comparison_module._parse_official_snapshots

    def parse_then_mutate_live_sources(*args, **kwargs):
        parsed = original_parser(*args, **kwargs)
        changed = json.loads(manifest_path.read_text(encoding="utf-8"))
        changed["parents"][0]["cases"][0]["space_group_type"] = 123
        manifest_path.write_text(json.dumps(changed), encoding="utf-8")
        result_path.write_text("not the signed result table", encoding="utf-8")
        return parsed

    monkeypatch.setattr(
        comparison_module,
        "_parse_official_snapshots",
        parse_then_mutate_live_sources,
    )
    observed_queries: list[dict] = []

    class SnapshotIsoDistort:
        def load_structure(self, _path: Path) -> None:
            pass

        def search_method_3(self, **query):
            observed_queries.append(query)
            return [_local_item()]

    report = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=SnapshotIsoDistort,
    )

    assert observed_queries[0]["space_group_type"] == 139
    assert report["cases"][0]["status"] == "exact_match"
    assert report["source_signature"]["stable"] is False
    assert report["summary"]["complete"] is False
    assert comparison_module._comparison_exit_code(report) == 1


def test_keyboard_interrupt_immediately_replaces_stale_complete_report_and_checkpoint(
    tmp_path: Path,
) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    _write_result_page(
        output_root / "Parent.cif" / "官网" / "Method3" / "case-one" / "result.html"
    )
    report_path.write_text(
        json.dumps({"schema": comparison_module.REPORT_SCHEMA, "summary": {"complete": True}}),
        encoding="utf-8",
    )

    class InterruptedIsoDistort:
        def load_structure(self, _path: Path) -> None:
            pass

        def search_method_3(self, **_query):
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        run_comparison(
            manifest_path=manifest_path,
            output_compare_root=output_root,
            input_root=input_root,
            report_path=report_path,
            checkpoint_path=checkpoint_path,
            iso_factory=InterruptedIsoDistort,
        )

    interrupted_report = json.loads(report_path.read_text(encoding="utf-8"))
    assert interrupted_report["schema"] == comparison_module.REPORT_SCHEMA
    assert interrupted_report["summary"]["complete"] is False
    assert interrupted_report["summary"]["validation_status"] == "in_progress"
    assert interrupted_report["summary"]["source_signature_stable"] is None
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert checkpoint["source_signature"]["stable"] is None


def test_unstable_checkpoint_results_are_never_reused() -> None:
    signature = "signed-A"
    record = {"status": "exact_match"}
    for stable in (None, False):
        checkpoint = comparison_module._load_checkpoint(
            signature,
            existing={
                "schema": comparison_module.CHECKPOINT_SCHEMA,
                "signature": signature,
                "source_signature": {
                    "start": signature,
                    "end": signature,
                    "stable": stable,
                },
                "results": {"Parent.cif::M3-T-01": record},
            },
        )
        assert checkpoint["results"] == {}
        assert checkpoint["source_signature"]["stable"] is None


def test_empty_selection_is_final_but_never_complete(tmp_path: Path) -> None:
    manifest_path = tmp_path / "manifest.json"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    input_root = tmp_path / "experiment_data"
    input_root.mkdir()
    manifest_path.write_text(
        json.dumps(
            {
                "official_download_root": (
                    "output_compare/<parent>/官网/Method3/<readable-case-folder>/"
                ),
                "parents": [],
            }
        ),
        encoding="utf-8",
    )

    report = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=tmp_path / "output_compare",
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=lambda: pytest.fail("empty selection must not construct a backend"),
    )

    assert report["summary"]["selected_case_count"] == 0
    assert report["summary"]["comparable_case_count"] == 0
    assert report["summary"]["source_signature_stable"] is True
    assert report["summary"]["complete"] is False
    assert report["summary"]["validation_status"] == "no_comparable_cases"
    assert comparison_module._comparison_exit_code(report) == 2


def test_signature_recheck_exception_invalidates_report_and_checkpoint(
    tmp_path: Path,
    monkeypatch,
) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    _write_result_page(
        output_root / "Parent.cif" / "官网" / "Method3" / "case-one" / "result.html"
    )
    signature_calls = 0

    def signature_then_error(*_args, **_kwargs) -> str:
        nonlocal signature_calls
        signature_calls += 1
        if signature_calls == 1:
            return "signed-snapshot"
        raise RuntimeError("recheck unavailable")

    monkeypatch.setattr(comparison_module, "_run_signature", signature_then_error)

    class FakeIsoDistort:
        def load_structure(self, _path: Path) -> None:
            pass

        def search_method_3(self, **_query):
            return [_local_item()]

    report = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=FakeIsoDistort,
    )

    assert report["source_signature"]["stable"] is False
    assert report["source_signature"]["end"] is None
    assert report["source_signature"]["recheck_error"] == (
        "RuntimeError: recheck unavailable"
    )
    assert report["summary"]["validation_status"] == "source_drift"
    assert comparison_module._comparison_exit_code(report) == 1
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert checkpoint["signature"] == "invalidated:source_signature_drift"
    assert checkpoint["results"] == {}


def test_runner_compares_explicit_empty_official_result(tmp_path: Path) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    case_dir = output_root / "Parent.cif" / "官网" / "Method3" / "case-one"
    case_dir.mkdir(parents=True)
    (case_dir / "empty-result.html").write_text(
        "Space Group: 139 I4/mmm D4h-17 "
        "Include strain, displacive ALL distortions "
        "Subgroup: 139 I4/mmm, basis of lattice = (1,0,0),(0,1,0),(0,0,1), I centering "
        '<input type="hidden" name="spacegroup" value="139 I4/mmm D4h-17">'
        '<input type="hidden" name="includestrain" value="true">'
        '<input type="hidden" name="includedisplacive001" value="true">'
        "Finish selecting the distortion mode by selecting one of the following choices: "
        "There are no subgroups with the space-group symmetry you selected. Try again.",
        encoding="utf-8",
    )

    class EmptyIsoDistort:
        def load_structure(self, path: Path) -> None:
            assert path == input_root / "Parent.cif"

        def search_method_3(self, **query):
            assert query["space_group_type"] == 139
            return []

    report = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=EmptyIsoDistort,
    )

    assert report["summary"]["comparable_case_count"] == 1
    assert report["summary"]["official_eligible_case_count"] == 1
    assert report["summary"]["exact_match_count"] == 1
    assert report["cases"][0]["official_download_audit_status"] == "verified"
    assert report["cases"][0]["comparison"]["official_count"] == 0
    assert report["cases"][0]["comparison"]["local_embedding_count"] == 0


def test_runner_can_compare_candidate_inventory_without_claiming_completeness(
    tmp_path: Path,
) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    case_dir = output_root / "Parent.cif" / "官网" / "Method3" / "case-one"
    candidate_dir = case_dir / "candidate"
    candidate_dir.mkdir(parents=True)
    (case_dir / "details.html").write_text(
        "<h1>ISODISTORT: complete modes details</h1>", encoding="utf-8"
    )
    (candidate_dir / "subgroup.cif").write_text(
        "# Subgroup: 139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
        "origin=(0,0,0), s=1, i=1\n",
        encoding="utf-8",
    )

    class FakeIsoDistort:
        def load_structure(self, path: Path) -> None:
            assert path == input_root / "Parent.cif"

        def search_method_3(self, **query):
            assert query["space_group_type"] == 139
            return [_local_item()]

    report = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        allow_candidate_inventory=True,
        iso_factory=FakeIsoDistort,
    )

    assert report["summary"]["comparable_case_count"] == 1
    assert report["summary"]["official_eligible_case_count"] == 0
    assert report["summary"]["authoritative_result_table_case_count"] == 0
    assert report["summary"]["provisional_candidate_inventory_case_count"] == 1
    assert report["cases"][0]["status"] == "exact_match"
    assert report["cases"][0]["official_source"] == "downloaded_candidate_inventory"
    assert report["cases"][0]["official_completeness_proven"] is False
    assert report["cases"][0]["provisional_reason"] == "missing_result_page"


def test_runner_rejects_official_page_bound_to_wrong_parent(tmp_path: Path) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    _write_result_page(
        output_root / "Parent.cif" / "官网" / "Method3" / "case-one" / "result.html",
        parent_sg=123,
    )

    class BackendMustNotRun:
        def __init__(self) -> None:
            raise AssertionError("manifest-mismatched official HTML must be rejected before local execution")

    report = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=BackendMustNotRun,
    )

    assert report["summary"]["official_eligible_case_count"] == 0
    assert report["summary"]["skipped_no_unique_official_result_count"] == 1
    assert report["cases"][0]["reason"] == "result_context_mismatch"


def test_run_signature_invalidates_for_core_and_iso_runtime_inputs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    package_root = tmp_path / "ISODISTORT"
    config_path = package_root / "resources" / "config" / "settings.yaml"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(
        json.dumps(
            {
                "isobyu": {
                    "bin_dir": "../../runtime",
                    "data_dir": "../../runtime",
                    "iso_bin": "iso",
                }
            }
        ),
        encoding="utf-8",
    )
    source_paths = [
        package_root / "backend" / "api" / "core_api.py",
        package_root / "backend" / "distortion" / "search_methods.py",
        package_root / "backend" / "distortion" / "superspace.py",
        package_root / "backend" / "backend" / "iso_wrapper.py",
        package_root / "backend" / "backend" / "smodes_wrapper.py",
        package_root / "backend" / "utils" / "lattice.py",
        package_root / "backend" / "data" / "kpoints_official.py",
        package_root / "tests" / "manual" / "audit_method3_downloads.py",
    ]
    for path in source_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"initial:{path.name}", encoding="utf-8")

    runtime_dir = package_root / "runtime"
    runtime_dir.mkdir()
    iso_path = runtime_dir / "iso"
    data_path = runtime_dir / "data_isotropy.txt"
    iso_path.write_bytes(b"iso-v1")
    data_path.write_bytes(b"data-v1")

    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text("{}", encoding="utf-8")
    input_root = tmp_path / "input"
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("parent", encoding="utf-8")
    selected = [{"parent_cif": "Parent.cif", "case": {}}]
    official_case = tmp_path / "official-case"
    official_case.mkdir()
    official_result = official_case / "result.html"
    official_result.write_text("official-v1", encoding="utf-8")
    official_subgroup = official_case / "subgroup.cif"
    official_subgroup.write_text("subgroup-v1", encoding="utf-8")
    official_inputs = {
        "Parent.cif::M3-T-01": {"case_directory": str(official_case)}
    }
    monkeypatch.setattr(comparison_module, "PACKAGE_ROOT", package_root)

    def signature() -> str:
        return comparison_module._run_signature(
            manifest_path, input_root, selected, official_inputs
        )

    signatures = [signature()]
    for path, replacement in (
        (package_root / "backend" / "utils" / "lattice.py", "lattice-v2"),
        (package_root / "backend" / "data" / "kpoints_official.py", "kpoints-v2"),
        (package_root / "backend" / "backend" / "smodes_wrapper.py", "smodes-v2"),
        (
            package_root / "tests" / "manual" / "audit_method3_downloads.py",
            "audit-v2",
        ),
    ):
        path.write_text(replacement, encoding="utf-8")
        signatures.append(signature())

    official_subgroup.write_text("subgroup-v2", encoding="utf-8")
    signatures.append(signature())

    iso_path.write_bytes(b"iso-v2")
    signatures.append(signature())
    data_path.write_bytes(b"data-v2")
    signatures.append(signature())
    iso_path.unlink()
    signatures.append(signature())

    assert len(set(signatures)) == len(signatures)


def test_run_signature_tracks_runtime_versions_and_marks_unsafe_parameter_k_cache(
    tmp_path: Path,
    monkeypatch,
) -> None:
    package_root = tmp_path / "ISODISTORT"
    (package_root / "resources" / "config").mkdir(parents=True)
    (package_root / "resources" / "config" / "settings.yaml").write_text("{}", encoding="utf-8")
    (package_root / "tests" / "manual").mkdir(parents=True)
    for name in ("audit_method3_downloads.py", "method3_affine_equivalence.py"):
        (package_root / "tests" / "manual" / name).write_text(name, encoding="utf-8")
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text("{}", encoding="utf-8")
    input_root = tmp_path / "input"
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("parent", encoding="utf-8")
    official_case = tmp_path / "official"
    official_case.mkdir()
    official_inputs = {"case": {"case_directory": str(official_case)}}
    ordinary = [{"parent_cif": "Parent.cif", "case": {"category": "identity"}}]
    parameter_k = [
        {"parent_cif": "Parent.cif", "case": {"category": "parameter-k"}}
    ]
    monkeypatch.setattr(comparison_module, "PACKAGE_ROOT", package_root)

    runtime = {"version": "runtime-v1"}
    monkeypatch.setattr(
        comparison_module,
        "_runtime_dependency_versions",
        lambda: dict(runtime),
    )

    ordinary_v1 = comparison_module._run_signature(
        manifest_path,
        input_root,
        ordinary,
        official_inputs,
    )
    runtime["version"] = "runtime-v2"
    ordinary_v2 = comparison_module._run_signature(
        manifest_path,
        input_root,
        ordinary,
        official_inputs,
    )
    assert ordinary_v1 != ordinary_v2

    policy = comparison_module._isotropy_cache_signature_policy(parameter_k)
    assert policy == {
        "status": "content_fingerprint_unavailable",
        "failure_marker": "public_cache_api_has_no_content_sha256",
        "resume_policy": "disable_cross_process_reuse_for_affected_cases",
        "affected_case_keys": ["Parent.cif::"],
    }
    assert comparison_module._isotropy_cache_signature_policy(ordinary) == {
        "status": "not_required",
        "affected_case_keys": [],
    }


def test_parameter_k_checkpoint_is_never_reused_across_invocations(
    tmp_path: Path,
) -> None:
    manifest = _manifest()
    manifest["parents"][0]["cases"] = manifest["parents"][0]["cases"][:1]
    manifest["parents"][0]["cases"][0]["category"] = "parameter-k"
    manifest["parents"][0]["cases"][0]["generate_if_missing"] = True
    manifest_path = tmp_path / "manifest.json"
    output_root = tmp_path / "output_compare"
    input_root = tmp_path / "experiment_data"
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    input_root.mkdir()
    (input_root / "Parent.cif").write_text("test parent", encoding="utf-8")
    _write_result_page(
        output_root / "Parent.cif" / "官网" / "Method3" / "case-one" / "result.html"
    )
    calls = {"search": 0}

    class FakeIsoDistort:
        def load_structure(self, _path: Path) -> None:
            pass

        def search_method_3(self, **_query):
            calls["search"] += 1
            return [_local_item()]

    first = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=FakeIsoDistort,
    )
    second = run_comparison(
        manifest_path=manifest_path,
        output_compare_root=output_root,
        input_root=input_root,
        report_path=report_path,
        checkpoint_path=checkpoint_path,
        iso_factory=FakeIsoDistort,
    )

    assert calls["search"] == 2
    assert first["cases"][0]["checkpoint_reuse_policy"] == (
        "disabled_mutable_generated_isotropy_cache"
    )
    assert second["summary"]["resumed_case_count"] == 0


def test_comparison_exit_code_only_accepts_clean_provisional_evidence() -> None:
    summary = {
        "complete": True,
        "source_signature_stable": True,
        "comparable_case_count": 1,
        "difference_count": 0,
        "local_error_count": 0,
        "skipped_no_unique_official_result_count": 0,
        "provisional_candidate_inventory_case_count": 1,
    }
    report = {"schema": comparison_module.REPORT_SCHEMA, "summary": summary}

    assert comparison_module._comparison_exit_code(report) == 2
    assert (
        comparison_module._comparison_exit_code(report, accept_provisional=True) == 0
    )

    assert comparison_module._comparison_exit_code({"summary": summary}) == 2
    report["schema"] = comparison_module.REPORT_SCHEMA - 1
    assert (
        comparison_module._comparison_exit_code(report, accept_provisional=True) == 2
    )
    report["schema"] = comparison_module.REPORT_SCHEMA
    summary["source_signature_stable"] = None
    assert (
        comparison_module._comparison_exit_code(report, accept_provisional=True) == 2
    )
    summary["source_signature_stable"] = True
    summary["complete"] = False
    assert (
        comparison_module._comparison_exit_code(report, accept_provisional=True) == 2
    )
    summary["complete"] = True

    summary["skipped_no_unique_official_result_count"] = 1
    assert (
        comparison_module._comparison_exit_code(report, accept_provisional=True) == 2
    )
    summary["skipped_no_unique_official_result_count"] = 0
    summary["difference_count"] = 1
    assert (
        comparison_module._comparison_exit_code(report, accept_provisional=True) == 1
    )


def test_comparison_cli_accept_provisional_flag_is_narrow(monkeypatch) -> None:
    report = {
        "schema": comparison_module.REPORT_SCHEMA,
        "summary": {
            "complete": True,
            "source_signature_stable": True,
            "comparable_case_count": 1,
            "difference_count": 0,
            "local_error_count": 0,
            "skipped_no_unique_official_result_count": 0,
            "provisional_candidate_inventory_case_count": 1,
        }
    }
    monkeypatch.setattr(comparison_module, "run_comparison", lambda **_kwargs: report)
    monkeypatch.setattr(comparison_module, "_print_summary", lambda *_args: None)

    assert comparison_module.main([]) == 2
    assert comparison_module.main(["--accept-provisional"]) == 0

    report["summary"]["skipped_no_unique_official_result_count"] = 1
    assert comparison_module.main(["--accept-provisional"]) == 2
