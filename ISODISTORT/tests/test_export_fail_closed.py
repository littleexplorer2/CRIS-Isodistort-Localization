"""Fail-closed boundaries for batch and IsoVIZ mode export."""
from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import threading
import time
import zipfile
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

import backend.api.core_api as core_api_module
from backend.api.core_api import (
    ExportBatchPlan,
    ExportCandidateFailure,
    ExportPreparationError,
    IsoDistort,
    ModeIdentityFailure,
    UnresolvedModeIdentityError,
)
from backend.wrappers import SubgroupInfo
from features.export.isodistort_isoviz import (
    _AtomLayout,
    _authoritative_mode_parentatom,
)


def _candidate(index: int, irrep: str) -> SubgroupInfo:
    return SubgroupInfo(
        index=index,
        space_group_number=12,
        subgroup_index=2,
        space_group_symbol="C2/m",
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label=irrep,
        opd_symbol="P1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        parent_sg=139,
    )


def _bare_export_session() -> IsoDistort:
    iso = object.__new__(IsoDistort)
    iso._selected_subgroup = None
    iso.phase_path = None
    iso.distortion_modes = []
    iso.mode_displacements = {}
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.distorted_structure = None
    iso.number_of_independent_modulations = 0
    iso.structure = None
    iso.distortion_types = ["strain", "displacive"]
    iso.distortion_scope = {}
    iso._mode_cache_key = None
    iso._generated_structure_cache_key = None
    return iso


def test_export_snapshot_owns_independent_structure_candidates_and_backends(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(core_api_module, "IsoWrapper", object)
    monkeypatch.setattr(core_api_module, "SmodesWrapper", object)
    iso = IsoDistort()
    iso.structure = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    iso.symmetry_info = {
        "space_group_number": 221,
        "space_group_symbol": "Pm-3m",
        "wyckoff_sites": [],
    }
    iso.subgroups = [_candidate(0, "GM1+")]
    iso.distortion_types = ["displacive"]
    iso.distortion_scope = {"displacive": ["Fe"]}

    worker = iso.snapshot_for_export()
    worker.structure.translate_sites([0], [0.25, 0.0, 0.0], frac_coords=True)
    worker.subgroups[0].index = 99
    worker.distortion_scope["displacive"].append("O")

    assert worker is not iso
    assert worker._iso is not iso._iso
    assert np.allclose(iso.structure[0].frac_coords, [0.0, 0.0, 0.0])
    assert iso.subgroups[0].index == 0
    assert iso.distortion_scope == {"displacive": ["Fe"]}


def test_batch_export_aggregates_mode_failures_before_building_specs() -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    iso = _bare_export_session()
    built: list[SubgroupInfo] = []

    def fail_mode_search(subgroup_idx, **_kwargs):
        raise ValueError(f"mode failure {subgroup_idx}")

    iso.search_method_2 = fail_mode_search
    iso._spec_for_subgroup = lambda subgroup, **_kwargs: built.append(subgroup)

    with pytest.raises(ExportPreparationError) as caught:
        iso._collect_export_specs(candidates, ["modes"], True)

    assert built == []
    assert [failure.candidate_index for failure in caught.value.failures] == [0, 1]
    assert [failure.subgroup_index for failure in caught.value.failures] == [2, 2]
    assert all(
        failure.parent_space_group_number == 139
        and failure.k_coordinates == ("0", "0", "0")
        and failure.k_parameters == ()
        for failure in caught.value.failures
    )
    assert caught.value.failures[0].basis_vectors == (
        ("1", "0", "0"),
        ("0", "1", "0"),
        ("0", "0", "1"),
    )
    assert [failure.error_type for failure in caught.value.failures] == [
        "ValueError",
        "ValueError",
    ]
    assert iso._selected_subgroup is None
    assert iso.mode_displacements == {}


def test_successful_empty_mode_table_is_distinct_from_calculation_failure() -> None:
    candidate = _candidate(0, "GM1+")
    iso = _bare_export_session()

    def empty_mode_search(subgroup_idx, *, candidates, **_kwargs):
        iso._selected_subgroup = candidates[subgroup_idx]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}

    iso.search_method_2 = empty_mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    specs = iso._collect_export_specs([candidate], ["modes"], True)

    assert len(specs) == 1
    assert specs[0].use_current_modes is True
    assert specs[0].note == (
        "no displacement modes for this path (empty BUSH/smodes table)"
    )


def test_collection_can_skip_only_unresolved_identity_candidates() -> None:
    ready = _candidate(0, "GM1+")
    unresolved = _candidate(1, "GM2+")
    iso = _bare_export_session()

    def mode_search(subgroup_idx, *, candidates, **_kwargs):
        candidate = candidates[0]
        iso._selected_subgroup = candidate
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        if subgroup_idx == unresolved.index:
            raise UnresolvedModeIdentityError([
                ModeIdentityFailure(
                    amplitude_key="GM2+__a__unresolved(a)",
                    orbit_id="sg139:1a:Fe:test",
                    global_irrep="GM2+",
                    k_coordinates=("0", "0", "0"),
                    reason=(
                        "iso_microscopic_block_missing; parent_SG=139; "
                        "KVALUE=(); direction_selector=(a)"
                    ),
                )
            ])
        return SimpleNamespace(subgroup=candidate)

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    plan = iso._collect_export_specs(
        [ready, unresolved],
        ["modes"],
        True,
        allow_unresolved_identities=True,
    )

    assert isinstance(plan, ExportBatchPlan)
    assert [spec.subgroup for spec in plan.specs] == [ready]
    assert [failure.candidate_index for failure in plan.skipped] == [1]
    assert "iso_microscopic_block_missing" in plan.skipped[0].message
    assert "parent_SG=139" in plan.skipped[0].message


def test_collection_cancels_at_candidate_boundary_and_keeps_completed_spec() -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    iso = _bare_export_session()
    cancelled = False
    events: list[dict] = []

    def mode_search(_subgroup_idx, *, candidates, **_kwargs):
        iso._selected_subgroup = candidates[0]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        return SimpleNamespace(subgroup=candidates[0])

    def progress(event):
        nonlocal cancelled
        events.append(event)
        if event.get("phase") == "candidate_complete":
            cancelled = True

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    plan = iso._collect_export_specs(
        candidates,
        ["modes"],
        True,
        allow_unresolved_identities=True,
        progress_callback=progress,
        cancel_check=lambda: cancelled,
    )

    assert isinstance(plan, ExportBatchPlan)
    assert [spec.subgroup for spec in plan.specs] == [candidates[0]]
    assert [failure.error_type for failure in plan.skipped] == [
        "ExportCancelledError"
    ]
    assert any(event["phase"] == "cancelled" for event in events)
    completed = next(
        event for event in events if event["phase"] == "candidate_complete"
    )
    assert completed["completed"] == 1
    assert completed["estimated_remaining_seconds"] >= 0


def test_collection_uses_singleton_candidate_pools_for_duplicate_local_indices() -> None:
    first = _candidate(0, "GM1+")
    second = _candidate(0, "GM2+")
    iso = _bare_export_session()
    pools: list[tuple[SubgroupInfo, ...]] = []

    def mode_search(subgroup_idx, *, candidates, **_kwargs):
        assert subgroup_idx == 0
        pools.append(tuple(candidates))
        iso._selected_subgroup = candidates[0]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        return SimpleNamespace(subgroup=candidates[0])

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    specs = iso._collect_export_specs([first, second], ["modes"], True)

    assert pools == [(first,), (second,)]
    assert [spec.subgroup for spec in specs] == [first, second]


def test_collection_rejects_search_result_with_different_scientific_identity() -> None:
    requested = _candidate(0, "GM1+")
    wrong = _candidate(0, "GM2+")
    iso = _bare_export_session()

    def wrong_search(_subgroup_idx, **_kwargs):
        iso._selected_subgroup = wrong
        return SimpleNamespace(subgroup=wrong)

    iso.search_method_2 = wrong_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    with pytest.raises(ExportPreparationError) as caught:
        iso._collect_export_specs([requested], ["modes"], True)

    assert caught.value.failures[0].candidate_index == 0
    assert "different scientific subgroup" in caught.value.failures[0].message


def test_disk_export_does_not_create_root_when_preparation_fails(tmp_path) -> None:
    candidate = _candidate(0, "GM1+")
    failure = ExportCandidateFailure(
        position=1,
        candidate_index=0,
        subgroup_index=2,
        parent_space_group_number=139,
        k_point_label="GM",
        k_coordinates=("0", "0", "0"),
        k_parameters=(),
        irrep_label="GM1+",
        opd_symbol="P1",
        space_group_number=12,
        basis_vectors=(("1", "0", "0"), ("0", "1", "0"), ("0", "0", "1")),
        origin=("0.0", "0.0", "0.0"),
        embedding_id="embedding-test",
        error_type="ValueError",
        message="mode failure",
    )
    error = ExportPreparationError(2, [failure])
    iso = object.__new__(IsoDistort)
    iso.structure = object()
    iso.subgroups = [candidate]

    def fail_collection(*_args, **_kwargs):
        raise error

    iso._collect_export_specs = fail_collection
    destination = tmp_path / "batch"

    with pytest.raises(ExportPreparationError):
        iso.export_subgroups(destination, formats=["modes"])

    assert not destination.exists()


def test_disk_batch_prerenders_every_candidate_before_publication(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    specs = [
        SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}")
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(IsoDistort)
    iso.structure = object()
    iso.subgroups = candidates
    iso.symmetry_info = {"space_group_number": 139}
    iso._collect_export_specs = lambda *_args, **_kwargs: specs
    rendered: list[str] = []

    def render(spec, _formats):
        rendered.append(spec.folder_name)
        if spec is specs[1]:
            raise ValueError("later writer validation failed")
        return (("candidate-0.txt", b"ready"),)

    monkeypatch.setattr(core_api_module, "render_subgroup_files", render)
    destination = tmp_path / "batch"

    with pytest.raises(ExportPreparationError, match="later writer validation failed"):
        iso.export_subgroups(destination, formats=["modes"])

    assert rendered == ["candidate-0", "candidate-1"]
    assert not destination.exists()


def test_anonymous_nonempty_modes_fail_during_collection() -> None:
    candidate = _candidate(4, "GM1+")
    iso = _bare_export_session()
    iso._selected_subgroup = candidate
    iso.mode_displacements_sc = {"anonymous": np.ones((1, 3))}
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=iso.mode_displacements_sc,
        displacive_data=None,
        **kwargs,
    )

    with pytest.raises(ExportPreparationError) as caught:
        iso._collect_export_specs([candidate], ["modes"], False)

    failure = caught.value.failures[0]
    assert failure.candidate_index == 4
    assert failure.subgroup_index == 2
    assert failure.parent_space_group_number == 139
    assert failure.k_coordinates == ("0", "0", "0")
    assert "lack verified identity" in failure.message


def test_supercell_modes_prevent_false_empty_mode_note() -> None:
    candidate = _candidate(0, "GM1+")
    iso = _bare_export_session()

    def mode_search(_subgroup_idx, *, candidates, **_kwargs):
        iso._selected_subgroup = candidates[0]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {"mode": np.ones((1, 3))}

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=iso.mode_displacements_sc,
        displacive_data=object(),
        **kwargs,
    )

    specs = iso._collect_export_specs([candidate], ["modes"], True)

    assert len(specs) == 1
    assert specs[0].note == ""


def test_export_recompute_receives_current_types_and_scope_explicitly() -> None:
    candidate = _candidate(0, "GM1+")
    iso = _bare_export_session()
    iso.structure = Structure(
        Lattice.cubic(4.0),
        ["Fe", "O"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    iso.distortion_types = ["displacive"]
    iso.distortion_scope = {"displacive": ["Fe"]}
    iso._selected_subgroup = candidate
    old_scope = {"displacive": ["O"]}
    iso._mode_cache_key = iso._mode_context_key(
        candidate,
        0,
        ["displacive"],
        old_scope,
    )
    calls: list[dict] = []

    def mode_search(_subgroup_idx, **kwargs):
        calls.append(kwargs)
        iso._selected_subgroup = candidate
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        return SimpleNamespace(subgroup=candidate)

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    specs = iso._collect_export_specs([candidate], ["modes"], True)

    assert len(specs) == 1
    assert calls == [
        {
            "distortion_type": ["displacive"],
            "number_of_independent_modulations": 0,
            "candidates": [candidate],
            "distortion_scope": {"displacive": ["Fe"]},
        }
    ]
    assert specs[0].use_current_modes is True
    assert specs[0].use_generated_structure is False


@pytest.mark.parametrize(
    ("types", "scope"),
    (
        (["strain", "displacive"], {"displacive": ["O"]}),
        (["strain"], {"displacive": ["Fe"]}),
    ),
)
def test_cif_only_export_rejects_generated_structure_from_old_context(
    types,
    scope,
) -> None:
    candidate = _candidate(0, "GM1+")
    iso = _bare_export_session()
    iso.structure = Structure(
        Lattice.cubic(4.0),
        ["Fe", "O"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    iso._selected_subgroup = candidate
    iso.distortion_types = types
    iso.distortion_scope = scope
    old_key = iso._mode_context_key(
        candidate,
        0,
        ["strain", "displacive"],
        {"displacive": ["Fe"]},
    )
    iso._mode_cache_key = old_key
    iso._generated_structure_cache_key = old_key
    iso.distorted_structure = object()
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    [spec] = iso._collect_export_specs([candidate], ["cif"], False)

    assert spec.use_current_modes is False
    assert spec.use_generated_structure is False


def test_cif_only_export_reuses_generated_structure_for_exact_context() -> None:
    candidate = _candidate(0, "GM1+")
    iso = _bare_export_session()
    iso.structure = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    iso._selected_subgroup = candidate
    iso.distortion_types = ["strain", "displacive"]
    iso.distortion_scope = {"displacive": ["Fe"]}
    key = iso._mode_context_key(
        candidate,
        0,
        iso.distortion_types,
        iso.distortion_scope,
    )
    iso._mode_cache_key = key
    iso._generated_structure_cache_key = key
    iso.distorted_structure = object()
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    [spec] = iso._collect_export_specs([candidate], ["cif"], False)

    assert spec.use_current_modes is False
    assert spec.use_generated_structure is True


def test_zip_prerenders_and_aggregates_all_writer_failures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    specs = [
        SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}")
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(IsoDistort)
    iso.structure = object()
    iso.subgroups = candidates
    iso._collect_export_specs = lambda *_args, **_kwargs: specs
    archive_called = False

    def render(spec, _formats):
        raise ValueError(f"writer failed for {spec.folder_name}")

    def archive(*_args, **_kwargs):
        nonlocal archive_called
        archive_called = True
        return b""

    monkeypatch.setattr(core_api_module, "render_subgroup_files", render)
    monkeypatch.setattr(core_api_module, "build_export_zip", archive)

    with pytest.raises(ExportPreparationError) as caught:
        iso.export_subgroups_zip(formats=["modes"])

    assert archive_called is False
    assert [failure.candidate_index for failure in caught.value.failures] == [0, 1]
    assert [failure.irrep_label for failure in caught.value.failures] == [
        "GM1+",
        "GM2+",
    ]
    assert all("writer failed" in failure.message for failure in caught.value.failures)


def test_zip_uses_each_prerendered_candidate_payload_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    specs = [
        SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}")
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(IsoDistort)
    iso.structure = object()
    iso.subgroups = candidates
    iso._collect_export_specs = lambda *_args, **_kwargs: specs
    calls: list[str] = []

    def render(spec, _formats):
        calls.append(spec.folder_name)
        return (("result.txt", spec.folder_name.encode("ascii")),)

    monkeypatch.setattr(core_api_module, "render_subgroup_files", render)

    archive = iso.export_subgroups_zip(formats=["modes"])

    assert calls == ["candidate-0", "candidate-1"]
    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        assert zipped.namelist() == [
            "candidate-0/result.txt",
            "candidate-1/result.txt",
        ]
        assert zipped.read("candidate-1/result.txt") == b"candidate-1"


def test_method2_zip_keeps_ready_candidate_and_embeds_skipped_identity_report(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ready = _candidate(0, "GM1+")
    skipped = _candidate(1, "GM2+")
    spec = SimpleNamespace(subgroup=ready, folder_name="candidate-ready")
    failure = IsoDistort._export_candidate_failure(
        _bare_export_session(),
        2,
        skipped,
        UnresolvedModeIdentityError([
            ModeIdentityFailure(
                amplitude_key="GM2+__a__unresolved(a)",
                orbit_id="sg139:1a:Fe:test",
                global_irrep="GM2+",
                k_coordinates=("0", "0", "0"),
                reason="candidate_columns_not_independent; parent_SG=139",
            )
        ]),
    )
    iso = object.__new__(IsoDistort)
    iso.structure = object()
    iso.subgroups = [ready, skipped]
    iso._collect_export_specs = lambda *_args, **_kwargs: ExportBatchPlan(
        (spec,), (failure,)
    )
    monkeypatch.setattr(
        core_api_module,
        "render_subgroup_files",
        lambda _spec, _formats: (("result.txt", b"ready"),),
    )

    archive = iso.export_subgroups_zip(
        formats=["modes"],
        compute_missing_modes=True,
        export_method=2,
    )

    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        assert zipped.read("candidate-ready/result.txt") == b"ready"
        report = json.loads(zipped.read("export_candidate_status.json"))
        assert report["status"] == "partial_success"
        assert report["successful_candidate_count"] == 1
        assert report["ineligible_candidate_count"] == 1
        assert report["failures"][0]["candidate_index"] == 1
        assert "candidate_columns_not_independent" in zipped.read(
            "export_candidate_status.txt"
        ).decode("utf-8")


def test_method2_disk_export_atomically_publishes_ready_candidate_and_report(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ready = _candidate(0, "GM1+")
    skipped = _candidate(1, "GM2+")
    spec = SimpleNamespace(subgroup=ready, folder_name="candidate-ready")
    failure = IsoDistort._export_candidate_failure(
        _bare_export_session(),
        2,
        skipped,
        UnresolvedModeIdentityError([
            ModeIdentityFailure(
                amplitude_key="GM2+__a__unresolved(a)",
                orbit_id="sg139:1a:Fe:test",
                global_irrep="GM2+",
                k_coordinates=("0", "0", "0"),
                reason="rank_mismatch; parent_SG=139; KVALUE=()",
            )
        ]),
    )
    iso = object.__new__(IsoDistort)
    iso.structure = object()
    iso.subgroups = [ready, skipped]
    iso._collect_export_specs = lambda *_args, **_kwargs: ExportBatchPlan(
        (spec,), (failure,)
    )
    monkeypatch.setattr(
        core_api_module,
        "render_subgroup_files",
        lambda _spec, _formats: (("result.txt", b"ready"),),
    )
    destination = tmp_path / "batch"

    paths = iso.export_subgroups(
        destination,
        formats=["modes"],
        compute_missing_modes=True,
        export_method=2,
    )

    assert (destination / "candidate-ready" / "result.txt").read_bytes() == b"ready"
    report = json.loads(
        (destination / "export_candidate_status.json").read_text(encoding="utf-8")
    )
    assert report["status"] == "partial_success"
    assert report["successful_candidate_count"] == 1
    assert report["failures"][0]["candidate_index"] == 1
    assert set(paths) == {
        destination / "candidate-ready" / "result.txt",
        destination / "export_candidate_status.json",
        destination / "export_candidate_status.txt",
    }


def test_disk_batch_staging_write_failure_publishes_nothing(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}"),
            (("result.txt", b"ready"),),
        )
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(IsoDistort)
    original_write = core_api_module.Path.write_bytes

    def fail_second_write(path, payload):
        if path.parent.name == "candidate-1":
            raise OSError("staged write failed")
        return original_write(path, payload)

    monkeypatch.setattr(core_api_module.Path, "write_bytes", fail_second_write)
    destination = tmp_path / "batch"

    with pytest.raises(ExportPreparationError) as caught:
        iso._publish_export_batch(2, destination, rendered)

    assert caught.value.failures[0].candidate_index == 1
    assert "staged write failed" in caught.value.failures[0].message
    assert not destination.exists()
    assert not list(tmp_path.glob(".isodistort-batch-*"))


def test_existing_disk_root_commit_failure_keeps_the_entire_batch_invisible(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}"),
            (("result.txt", b"ready"),),
        )
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    sentinel = destination / "keep.txt"
    sentinel.write_bytes(b"unchanged")
    original_rename = core_api_module.Path.rename

    def fail_batch_commit(path, target):
        if (
            path.name.startswith(".isodistort-batch-")
            and target.parent == destination
            and target.name.endswith(".ready")
        ):
            raise OSError("batch commit failed")
        return original_rename(path, target)

    monkeypatch.setattr(core_api_module.Path, "rename", fail_batch_commit)

    with pytest.raises(ExportPreparationError) as caught:
        iso._publish_export_batch(2, destination, rendered)

    assert [failure.candidate_index for failure in caught.value.failures] == [0, 1]
    assert all(
        "batch commit failed" in failure.message for failure in caught.value.failures
    )
    assert list(destination.iterdir()) == [sentinel]
    assert sentinel.read_bytes() == b"unchanged"
    assert not list(tmp_path.glob(".isodistort-batch-*"))


def test_existing_disk_root_exposes_a_complete_version_in_one_transition(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}"),
            (("result.txt", f"ready-{index}".encode("ascii")),),
        )
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"unchanged")
    original_rename = core_api_module.Path.rename
    commit_entered = threading.Event()
    committed = threading.Event()
    observed_before = threading.Event()
    observed_after = threading.Event()
    stop_reader = threading.Event()
    snapshots: list[tuple[tuple[str, ...], ...]] = []

    def ready_snapshot() -> tuple[tuple[str, ...], ...]:
        batches = sorted(
            path
            for path in destination.iterdir()
            if path.name.startswith(core_api_module._EXPORT_BATCH_PREFIX)
            and path.name.endswith(core_api_module._EXPORT_BATCH_READY_SUFFIX)
        )
        return tuple(
            tuple(
                sorted(
                    entry.name
                    for entry in batch.iterdir()
                    if entry.name != core_api_module._EXPORT_BATCH_MANIFEST
                )
            )
            for batch in batches
        )

    def read_while_publishing() -> None:
        while not stop_reader.is_set():
            if not commit_entered.is_set():
                stop_reader.wait(0.001)
                continue
            snapshot = ready_snapshot()
            snapshots.append(snapshot)
            if not snapshot:
                observed_before.set()
            if committed.is_set() and snapshot:
                observed_after.set()
            stop_reader.wait(0.001)

    def observe_single_commit(path, target):
        if (
            path.name.startswith(".isodistort-batch-")
            and target.parent == destination
            and target.name.endswith(".ready")
        ):
            commit_entered.set()
            assert observed_before.wait(5)
            result = original_rename(path, target)
            committed.set()
            assert observed_after.wait(5)
            return result
        return original_rename(path, target)

    monkeypatch.setattr(core_api_module.Path, "rename", observe_single_commit)
    reader = threading.Thread(target=read_while_publishing)
    reader.start()
    try:
        paths = iso._publish_export_batch(2, destination, rendered)
    finally:
        stop_reader.set()
        reader.join(timeout=5)

    complete = (("candidate-0", "candidate-1"),)
    assert snapshots
    assert () in snapshots
    assert complete in snapshots
    assert all(snapshot in {(), complete} for snapshot in snapshots)
    assert len({path.parents[1] for path in paths}) == 1
    assert all(path.read_bytes().startswith(b"ready-") for path in paths)
    assert core_api_module._published_export_folder_names(destination) >= {
        "candidate-0",
        "candidate-1",
    }
    assert not list(tmp_path.glob(".isodistort-batch-*"))


def test_stale_staging_directory_is_not_discoverable_as_a_committed_batch(
    tmp_path,
) -> None:
    candidate = _candidate(0, "GM1+")
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name="candidate-0"),
            (("result.txt", b"ready"),),
        )
    ]
    iso = object.__new__(IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    stale = destination / ".isodistort-batch-abandoned"
    (stale / "candidate-stale").mkdir(parents=True)
    (stale / "candidate-stale" / "result.txt").write_bytes(b"partial")

    before = core_api_module._published_export_folder_names(destination)
    paths = iso._publish_export_batch(2, destination, rendered)
    after = core_api_module._published_export_folder_names(destination)

    assert "candidate-stale" not in before
    assert "candidate-stale" not in after
    assert "candidate-0" in after
    assert paths[0].read_bytes() == b"ready"
    ready_batches = [
        path
        for path in destination.iterdir()
        if path.name.startswith(core_api_module._EXPORT_BATCH_PREFIX)
        and path.name.endswith(core_api_module._EXPORT_BATCH_READY_SUFFIX)
    ]
    assert len(ready_batches) == 1


@pytest.mark.parametrize(
    ("tamper", "message"),
    (
        ("manifest", "manifest does not match"),
        ("payload", "file hash changed"),
    ),
)
def test_tampered_ready_batch_fails_closed_before_another_publish(
    tmp_path,
    tamper: str,
    message: str,
) -> None:
    candidates = [_candidate(0, "GM1+"), _candidate(1, "GM2+")]
    first = [
        (
            SimpleNamespace(subgroup=candidates[0], folder_name="candidate-0"),
            (("result.txt", b"ready"),),
        )
    ]
    second = [
        (
            SimpleNamespace(subgroup=candidates[1], folder_name="candidate-1"),
            (("result.txt", b"later"),),
        )
    ]
    iso = object.__new__(IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"unchanged")
    [published] = iso._publish_export_batch(2, destination, first)
    ready_root = published.parents[1]
    if tamper == "manifest":
        manifest = ready_root / core_api_module._EXPORT_BATCH_MANIFEST
        manifest.write_bytes(manifest.read_bytes() + b" ")
    else:
        published.write_bytes(b"wrong")

    with pytest.raises(ExportPreparationError) as caught:
        iso._publish_export_batch(2, destination, second)

    assert message in caught.value.failures[0].message
    assert not any(path.name == "candidate-1" for path in destination.rglob("*"))
    assert not any(
        path.name.startswith(".isodistort-batch-")
        and not path.name.endswith(core_api_module._EXPORT_BATCH_READY_SUFFIX)
        for path in destination.iterdir()
    )


def test_concurrent_writers_reject_a_partially_overlapping_candidate_set(
    tmp_path,
) -> None:
    candidate_sets = (
        [
            (_candidate(0, "GM1+"), "candidate-shared", b"writer-a-shared"),
            (_candidate(1, "GM2+"), "candidate-a", b"writer-a-only"),
        ],
        [
            (_candidate(2, "GM3+"), "candidate-shared", b"writer-b-shared"),
            (_candidate(3, "GM4+"), "candidate-b", b"writer-b-only"),
        ],
    )
    destination = tmp_path / "batch"
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"unchanged")
    start = threading.Barrier(3)
    successes: list[list] = []
    failures: list[ExportPreparationError] = []

    def publish(entries) -> None:
        rendered = [
            (
                SimpleNamespace(subgroup=candidate, folder_name=folder),
                (("result.txt", payload),),
            )
            for candidate, folder, payload in entries
        ]
        start.wait()
        try:
            successes.append(
                object.__new__(IsoDistort)._publish_export_batch(
                    2,
                    destination,
                    rendered,
                )
            )
        except ExportPreparationError as exc:
            failures.append(exc)

    writers = [
        threading.Thread(target=publish, args=(entries,))
        for entries in candidate_sets
    ]
    for writer in writers:
        writer.start()
    start.wait()
    for writer in writers:
        writer.join(timeout=10)

    assert all(not writer.is_alive() for writer in writers)
    assert len(successes) == 1
    assert len(failures) == 1
    assert [failure.message for failure in failures[0].failures] == [
        f"target candidate directory exists in a committed batch: "
        f"{destination / 'candidate-shared'}"
    ]
    ready_batches = [
        path
        for path in destination.iterdir()
        if path.name.startswith(core_api_module._EXPORT_BATCH_PREFIX)
        and path.name.endswith(core_api_module._EXPORT_BATCH_READY_SUFFIX)
    ]
    assert len(ready_batches) == 1
    visible = {
        path.name
        for path in ready_batches[0].iterdir()
        if path.name != core_api_module._EXPORT_BATCH_MANIFEST
    }
    assert visible in (
        {"candidate-shared", "candidate-a"},
        {"candidate-shared", "candidate-b"},
    )
    assert not any(
        path.name.startswith(".isodistort-batch-")
        and not path.name.endswith(core_api_module._EXPORT_BATCH_READY_SUFFIX)
        for path in destination.iterdir()
    )


def test_publish_lock_serializes_independent_processes(tmp_path) -> None:
    destination = tmp_path / "batch"
    destination.mkdir()
    first_entered = tmp_path / "first-entered"
    second_entered = tmp_path / "second-entered"
    release_first = tmp_path / "release-first"
    # ``python -c`` does not put the working directory on sys.path, so the
    # child process needs the project root explicitly.
    project_root = core_api_module.Path(__file__).parents[1]
    child_env = dict(os.environ)
    child_env["PYTHONPATH"] = os.pathsep.join(
        [str(project_root), child_env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    script = """
from pathlib import Path
import sys
import time
from backend.api.core_api import _export_publish_lock

root = Path(sys.argv[1])
entered = Path(sys.argv[2])
release = None if sys.argv[3] == "-" else Path(sys.argv[3])
with _export_publish_lock(root):
    entered.write_text("entered", encoding="ascii")
    while release is not None and not release.exists():
        time.sleep(0.01)
"""
    first = subprocess.Popen(  # noqa: S603 - exact current interpreter
        [sys.executable, "-c", script, str(destination), str(first_entered), str(release_first)],
        cwd=project_root,
        env=child_env,
    )
    deadline = time.monotonic() + 5
    while not first_entered.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    assert first_entered.exists()

    second = subprocess.Popen(  # noqa: S603 - exact current interpreter
        [sys.executable, "-c", script, str(destination), str(second_entered), "-"],
        cwd=project_root,
        env=child_env,
    )
    time.sleep(0.2)
    assert not second_entered.exists()

    release_first.write_text("release", encoding="ascii")
    assert first.wait(timeout=5) == 0
    assert second.wait(timeout=5) == 0
    assert second_entered.exists()


def _two_type_layout() -> _AtomLayout:
    return _AtomLayout(
        types=[("A1", "A"), ("B1", "B")],
        child_type_idx=[1, 2],
        image_source_idx=[0, 1],
        orbit_type_idx={"orbit-a": 1, "orbit-b": 2},
    )


def _authoritative_mode(
    rows: np.ndarray,
    *,
    key: str = "mode-1",
    orbit_id: str = "orbit-b",
    type_index: int = 2,
    type_label: str = "B1",
) -> SimpleNamespace:
    return SimpleNamespace(
        key=key,
        parent_orbit_id=orbit_id,
        parent_type_index=type_index,
        parent_type_label=type_label,
        normalized_fractional_per_angstrom=rows,
    )


def test_authoritative_isoviz_mode_rejects_cross_type_column() -> None:
    rows = np.array([[0.25, 0.0, 0.0], [0.0, 0.5, 0.0]])

    with pytest.raises(ValueError, match="spans multiple IsoVIZ parent-atom types"):
        _authoritative_mode_parentatom(
            _authoritative_mode(rows),
            _two_type_layout(),
        )


def test_authoritative_isoviz_mode_rejects_empty_column() -> None:
    with pytest.raises(ValueError, match="has no nonzero child-atom row"):
        _authoritative_mode_parentatom(
            _authoritative_mode(np.zeros((2, 3))),
            _two_type_layout(),
        )


def test_authoritative_isoviz_mode_uses_nonzero_rows_not_display_label() -> None:
    rows = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.0]])

    assert (
        _authoritative_mode_parentatom(
            _authoritative_mode(rows, key="misleading-[A1:a:dsp]"),
            _two_type_layout(),
        )
        == 2
    )


def test_authoritative_isoviz_mode_rejects_rows_on_wrong_orbit_type() -> None:
    rows = np.array([[0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="expected type 2 for orbit 'orbit-b'"):
        _authoritative_mode_parentatom(
            _authoritative_mode(rows),
            _two_type_layout(),
        )


def test_authoritative_isoviz_mode_rejects_unmapped_orbit_identity() -> None:
    rows = np.array([[0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="absent from the IsoVIZ type layout"):
        _authoritative_mode_parentatom(
            _authoritative_mode(
                rows,
                orbit_id="orbit-missing",
                type_index=1,
                type_label="A1",
            ),
            _two_type_layout(),
        )
