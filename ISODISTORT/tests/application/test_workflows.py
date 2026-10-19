"""Tests for application.workflows; migrated without changing scientific expectations."""
from __future__ import annotations

from application import workflows
from application.session import ModeIdentityFailure
from application.workflows import service
from application.workflows import service as workflow_service
from application.workflows.jobs import ExportJobs
from backend.method3 import Method3ResultItem
from backend.wrappers import SubgroupInfo
from pathlib import Path
from tests.support_frontend_web import _SessionIsoStub
from tests.support_paths import TESTS_ROOT
from types import SimpleNamespace
import copy
import os
import pytest
import subprocess
import sys
import threading
import uuid


# Original coverage: test_application_workflows.py.

class _MutableEngine:
    def __init__(self):
        self.structure = None
        self.subgroups = ["old"]
        self.mode_displacements = {"old": [1.0]}
        self.mode_occupancies = {}
        self.distorted_structure = None
        self.closed = False

    def snapshot_for_action(self):
        return copy.deepcopy(self)

    def species(self):
        return []

    def close(self):
        self.closed = True

def test_invalid_scope_cannot_partially_change_the_live_session(monkeypatch):
    session = service.WebSession()
    session.revision = 4
    session._iso = SimpleNamespace()
    session.method1 = ["previous candidate"]
    monkeypatch.setattr(service, "_SESSION", session)

    with pytest.raises(ValueError, match="scope"):
        workflows.execute("set_types", {"types": ["strain"], "scope": []}, revision=4)

    assert session.distortion_types == ["strain", "displacive"]
    assert session.method1 == ["previous candidate"]
    assert session.revision == 4

@pytest.mark.parametrize("nmod", ["invalid", 1.5, True, -1, None])
def test_invalid_nmod_is_rejected_before_engine_access(monkeypatch, nmod):
    session = service.WebSession()
    monkeypatch.setattr(service, "_SESSION", session)
    monkeypatch.setattr(service, "IsoDistort", lambda: pytest.fail("engine constructed before validation"))

    with pytest.raises(ValueError, match="nmod"):
        workflows.execute("method2", {"subgroup_idx": 0, "nmod": nmod})

    assert session._iso is None
    assert session.revision == 0

def test_unsupported_modulation_count_is_explicit_before_worker_creation(monkeypatch):
    live = service.WebSession()
    monkeypatch.setattr(service, "_SESSION", live)
    monkeypatch.setattr(service, "IsoDistort", lambda: pytest.fail("unexpected engine"))
    with pytest.raises(NotImplementedError, match="nmod=0..3"):
        workflows.execute("method2", {"subgroup_idx": 0, "nmod": 4})
    with pytest.raises(NotImplementedError, match="nmod=0..3"):
        workflows.export_archive({"nmod": ["4"]})
    assert live._iso is None

def test_two_committed_actions_preserve_independent_generated_runtime_cache(monkeypatch):
    import application.session.core_api as core

    class RuntimeStub:
        def __init__(self):
            self.generated = {}
            self.closed = False

        def clone_runtime_from(self, source):
            self.generated = copy.deepcopy(source.generated)

        def close(self):
            self.closed = True
            self.generated.clear()

    monkeypatch.setattr(core, "IsoWrapper", RuntimeStub)
    monkeypatch.setattr(core, "SmodesWrapper", RuntimeStub)
    live = service.WebSession()
    live._iso = core.IsoDistort()
    monkeypatch.setattr(service, "_SESSION", live)

    def query(staged, data):
        if data.get("generate"):
            staged.iso._iso.generated["generic-table"] = [1, 2, 3]
        else:
            assert staged.iso._iso.generated["generic-table"] == [1, 2, 3]
        return {}

    monkeypatch.setitem(service._ACTIONS, "method1", query)
    workflows.execute("method1", {"generate": True})
    first_runtime = service._SESSION.iso._iso
    workflows.execute("method1", {"generate": False}, revision=1)

    assert first_runtime.closed
    assert first_runtime.generated == {}
    assert service._SESSION.iso._iso is not first_runtime
    assert service._SESSION.iso._iso.generated == {"generic-table": [1, 2, 3]}
    assert service._SESSION.revision == 2

def test_workflows_import_and_execute_without_frontend(tmp_path):
    project = TESTS_ROOT.parent
    probe = """
import sys
from types import SimpleNamespace
from application import workflows
from application.workflows import service

service._SESSION._iso = SimpleNamespace(list_k_points=lambda: [
    SimpleNamespace(label='GM', coordinates=[0, 0, 0], parameters=[],
                    is_special=True, kovalev='GM')
])
assert workflows.execute('kpoints')['kpoints'][0]['label'] == 'GM'
assert not any(name == 'frontend' or name.startswith('frontend.') for name in sys.modules)
assert 'http.server' not in sys.modules
"""
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["ISODISTORT_RUNTIME_ROOT"] = str(tmp_path)
    completed = subprocess.run(
        [sys.executable, "-B", "-c", probe], cwd=project, env=env,
        capture_output=True, text=True, check=False, timeout=30,
    )
    assert completed.returncode == 0, completed.stderr

def test_application_rejects_stale_revision_before_calling_the_engine(monkeypatch):
    calls = []
    session = service.WebSession()
    session.revision = 4
    session._iso = SimpleNamespace(list_k_points=lambda: calls.append(True))
    monkeypatch.setattr(service, "_SESSION", session)

    with pytest.raises(workflows.StaleContextError):
        workflows.execute("kpoints", revision=3)

    assert calls == []
    assert session.revision == 4

def test_generated_downloads_are_confined_by_application(tmp_path, monkeypatch):
    output = tmp_path / "output"
    output.mkdir()
    generated = output / "result.cif"
    generated.write_bytes(b"data_result\n")
    outside = tmp_path / "input.cif"
    outside.write_bytes(b"read only input")
    monkeypatch.setattr(service, "get_config", lambda: SimpleNamespace(output_dir=output))

    download = workflows.read_download("result.cif")
    assert download.filename == "result.cif"
    assert download.body == b"data_result\n"
    with pytest.raises(PermissionError):
        workflows.read_download("../input.cif")
    assert outside.read_bytes() == b"read only input"

def test_failed_export_releases_the_worker_for_the_next_request(monkeypatch):
    calls = []

    def export_subgroups_zip(**kwargs):
        calls.append(kwargs)
        if len(calls) == 1:
            raise RuntimeError("worker failed")
        return b"archive"

    session = service.WebSession()
    session._iso = SimpleNamespace(
        structure=object(),
        snapshot_for_export=lambda: SimpleNamespace(
            export_subgroups_zip=export_subgroups_zip, close=lambda: None,
        ),
    )
    session.method2_subgroups = [SimpleNamespace(index=0)]
    monkeypatch.setattr(service, "_SESSION", session)
    monkeypatch.setattr(service, "_EXPORT_JOBS", ExportJobs())
    try:
        with pytest.raises(workflows.ExportFailedError, match="worker failed"):
            workflows.export_archive({"method": ["2"], "formats": ["cif"]})
        assert workflows.export_progress_snapshot()["phase"] == "failed"

        download = workflows.export_archive({"method": ["2"], "formats": ["cif"]})
        assert download.body == b"archive"
        assert download.filename == "isodistort_method2.zip"
        assert workflows.export_progress_snapshot()["phase"] == "complete"
    finally:
        workflows.reset_export_progress()

def test_partial_upload_write_is_removed_without_changing_live_state(monkeypatch, tmp_path):
    live = service.WebSession()
    live.revision = 7
    monkeypatch.setattr(service, "_SESSION", live)
    monkeypatch.setattr(service, "get_config", lambda: SimpleNamespace(temp_dir=tmp_path))
    original_write = Path.write_text

    def partial_write(path, content, **kwargs):
        original_write(path, content[:4], **kwargs)
        raise OSError("disk write interrupted")

    monkeypatch.setattr(Path, "write_text", partial_write)
    with pytest.raises(OSError, match="disk write interrupted"):
        workflows.execute("load_cif", {"content": "data_input\n"})
    assert service._SESSION is live and live.revision == 7
    assert live._iso is None and live.parent_upload_path is None
    assert list(tmp_path.iterdir()) == []

@pytest.mark.parametrize("failure", [None, "constructor", "load", "response"])
def test_new_parent_uses_fresh_runtime_after_old_engine_is_disabled(monkeypatch, tmp_path, failure):
    class DisabledEngine(_MutableEngine):
        def snapshot_for_action(self):
            raise RuntimeError("disabled runtime cannot be cloned")

    live = service.WebSession()
    live._iso = DisabledEngine()
    live._iso.closed = True
    live._iso.generated_tables = {"previous-runtime": [1, 2, 3]}
    live.revision = 7
    live.distortion_types = ["occupational"]
    live.distortion_scope = {"occupational": ["X"]}
    live.method1 = ["previous candidates"]
    old_upload = tmp_path / "previous.cif"
    old_upload.write_text("previous input", encoding="utf-8")
    live.parent_upload_path = old_upload
    monkeypatch.setattr(service, "_SESSION", live)
    monkeypatch.setattr(service, "get_config", lambda: SimpleNamespace(temp_dir=tmp_path))

    class FreshEngine(_MutableEngine):
        def __init__(self):
            super().__init__()
            self.generated_tables = {}

        def load_structure(self, path):
            assert service._SESSION is live and old_upload.exists()
            self.input_text = Path(path).read_text(encoding="utf-8")
            if failure == "load":
                raise ValueError("new parent rejected")

        def set_distortion_scope(self, scope):
            self.scope = copy.deepcopy(scope)

        def set_distortion_types(self, types):
            self.types = list(types)

    fresh = FreshEngine()

    def construct():
        if failure == "constructor":
            raise ValueError("new runtime unavailable")
        return fresh

    monkeypatch.setattr(service, "IsoDistort", construct)
    if failure == "response":
        monkeypatch.setattr(service, "_state_summary", lambda _session: {"invalid": float("nan")})
    request = {"content": "data_new_parent\n", "filename": "new.cif"}
    if failure:
        with pytest.raises(ValueError):
            workflows.execute("load_cif", request)
        assert service._SESSION is live and live.revision == 7
        assert live.parent_upload_path == old_upload and old_upload.exists()
        assert list(tmp_path.iterdir()) == [old_upload]
        assert fresh.closed is (failure != "constructor")
    else:
        result = workflows.execute("load_cif", request)
        assert result["ok"] is True and result["state"]["revision"] == 8
        committed = service._SESSION
        assert committed is not live and committed.iso is fresh
        assert fresh.input_text == "data_new_parent\n" and not fresh.closed
        assert fresh.generated_tables == {}
        assert committed.distortion_types == fresh.types == ["occupational"]
        assert committed.distortion_scope == fresh.scope == {"occupational": ["X"]}
        assert committed.method1 == [] and committed.method2 is None and committed.method3 == []
        assert committed.method2_subgroups == []
        assert not old_upload.exists() and committed.parent_upload_path.is_file()
        committed.distortion_scope["occupational"].append("Y")
        assert live.distortion_scope == {"occupational": ["X"]}
        committed.cleanup_uploads()
        committed.close()
    assert live.iso.generated_tables == {"previous-runtime": [1, 2, 3]}
    assert live.method1 == ["previous candidates"]

def test_failed_computation_discards_worker_and_keeps_live_revision(monkeypatch):
    live = service.WebSession()
    live._iso = _MutableEngine()
    live.revision = 5
    live.method1 = ["old"]
    workers = []

    def fail_after_mutation(staged, _data):
        workers.append(staged.iso)
        staged.iso.mode_displacements["old"][0] = 9.0
        staged.method1.clear()
        raise ValueError("calculation failed after changing modes")

    monkeypatch.setattr(service, "_SESSION", live)
    monkeypatch.setitem(service._ACTIONS, "method1", fail_after_mutation)

    with pytest.raises(ValueError, match="calculation failed"):
        workflows.execute("method1", {}, revision=5)

    assert service._SESSION is live
    assert live.revision == 5
    assert live.method1 == ["old"]
    assert live.iso.mode_displacements == {"old": [1.0]}
    assert not live.iso.closed
    assert workers[0].closed

def test_success_commits_once_and_returns_a_detached_dto(monkeypatch):
    live = service.WebSession()
    live._iso = _MutableEngine()
    live.revision = 5

    def compute(staged, _data):
        staged.iso.mode_displacements["old"][0] = 2.0
        staged.method1 = [{"basis": [1, 2, 3]}]
        return {"candidates": staged.method1}

    monkeypatch.setattr(service, "_SESSION", live)
    monkeypatch.setitem(service._ACTIONS, "method1", compute)
    result = workflows.execute("method1", {}, revision=5)

    assert service._SESSION is not live
    assert service._SESSION.revision == result["state"]["revision"] == 6
    assert live.iso.mode_displacements == {"old": [1.0]}
    assert live.iso.closed
    result["candidates"][0]["basis"].append(4)
    assert service._SESSION.method1 == [{"basis": [1, 2, 3]}]

def test_close_session_releases_only_live_resources_and_can_restart(tmp_path, monkeypatch):
    upload = tmp_path / "parent.cif"
    upload.write_text("data_parent", encoding="utf-8")
    live = service.WebSession()
    live._iso = _MutableEngine()
    live.parent_upload_path = upload
    independent_worker = live.iso.snapshot_for_action()
    monkeypatch.setattr(service, "_SESSION", live)

    workflows.close_session()
    workflows.close_session()

    assert live.iso.closed
    assert not independent_worker.closed
    assert not upload.exists()
    assert service._SESSION._iso is None
    assert service._SESSION.revision == 0

def test_failed_response_preparation_does_not_commit_or_delete_previous_upload(tmp_path, monkeypatch):
    previous = tmp_path / "previous.cif"
    previous.write_text("previous input", encoding="utf-8")
    staged_path = tmp_path / "new.cif"
    staged_path.write_text("new input", encoding="utf-8")
    live = service.WebSession()
    live._iso = _MutableEngine()
    live.parent_upload_path = previous
    monkeypatch.setattr(service, "_SESSION", live)

    def load(staged, _data):
        staged.parent_upload_path = staged_path
        return {}

    def fail_summary(_staged):
        raise ValueError("invalid result DTO")

    monkeypatch.setitem(service._ACTIONS, "load_cif", load)
    monkeypatch.setattr(service, "_state_summary", fail_summary)
    with pytest.raises(ValueError, match="invalid result DTO"):
        workflows.execute("load_cif", {"content": "data_test"})

    assert service._SESSION is live
    assert live.revision == 0
    assert previous.read_text(encoding="utf-8") == "previous input"
    assert not staged_path.exists()

@pytest.mark.parametrize("nmod", ["invalid", "1.5", "-1", "", "True"])
def test_export_nmod_validation_precedes_worker_creation(monkeypatch, nmod):
    live = service.WebSession()
    monkeypatch.setattr(service, "_SESSION", live)
    monkeypatch.setattr(service, "IsoDistort", lambda: pytest.fail("unexpected engine"))
    with pytest.raises(ValueError, match="nmod"):
        workflows.export_archive({"nmod": [nmod]})
    assert live._iso is None

def test_late_progress_and_cancellation_cannot_change_the_next_export(monkeypatch):
    registry = ExportJobs()
    monkeypatch.setattr(service, "_EXPORT_JOBS", registry)
    live = service.WebSession()
    live.revision = 7
    live.method2_subgroups = [SimpleNamespace(index=0)]
    started = threading.Event()
    release = threading.Event()
    old_callback = []
    calls = []
    errors = []

    class Worker:
        def export_subgroups_zip(self, **kwargs):
            calls.append(kwargs)
            if len(calls) == 1:
                old_callback.append(kwargs["progress_callback"])
            else:
                started.set()
                assert release.wait(5)
                assert kwargs["cancel_check"]()
                kwargs["progress_callback"]({"phase": "cancelled", "completed": 0})
            return b"zip"

        def close(self):
            pass

    live._iso = SimpleNamespace(structure=object(), snapshot_for_export=Worker)
    monkeypatch.setattr(service, "_SESSION", live)
    first_id, second_id = str(uuid.uuid4()), str(uuid.uuid4())
    workflows.export_archive({"job_id": [first_id], "formats": ["cif"]})
    assert registry.snapshot()["phase"] == "complete"

    def second():
        try:
            workflows.export_archive({"job_id": [second_id], "formats": ["cif"]})
        except BaseException as exc:
            errors.append(exc)

    thread = threading.Thread(target=second)
    thread.start()
    try:
        assert started.wait(5)
        old_callback[0]({"phase": "complete", "completed": 99})
        assert not workflows.cancel_export(first_id)
        assert not workflows.cancel_export()
        progress = workflows.export_progress_snapshot()
        assert progress["job_id"] == second_id
        assert progress["active"] is True
        assert progress["source_revision"] == 7
        assert progress["completed"] == 0
        assert workflows.cancel_export(second_id)
    finally:
        release.set()
        thread.join(5)
    assert not thread.is_alive()
    assert errors == []
    assert registry.snapshot()["phase"] == "cancelled"
    assert not registry.snapshot()["active"]

def test_terminal_state_is_committed_before_the_next_job_can_start(monkeypatch):
    from application.workflows import jobs

    registry = ExportJobs()
    first = registry.start("first", revision=0, method=2, total=1)
    terminal_copy = threading.Event()
    allow_finish = threading.Event()
    second_started = threading.Event()
    errors = []
    original_deepcopy = copy.deepcopy

    def paused_copy(value):
        if isinstance(value, dict) and value.get("job_id") == "first" and value.get("phase") == "complete":
            terminal_copy.set()
            assert allow_finish.wait(5)
        return original_deepcopy(value)

    monkeypatch.setattr(jobs.copy, "deepcopy", paused_copy)

    def finish_first():
        try:
            registry.finish(first)
        except BaseException as exc:
            errors.append(exc)

    def start_second():
        try:
            registry.start("second", revision=1, method=2, total=2)
            second_started.set()
        except BaseException as exc:
            errors.append(exc)

    ending = threading.Thread(target=finish_first)
    beginning = threading.Thread(target=start_second)
    ending.start()
    try:
        assert terminal_copy.wait(5)
        beginning.start()
        assert not second_started.wait(0.05)
    finally:
        allow_finish.set()
        ending.join(5)
        if beginning.ident is not None:
            beginning.join(5)
    assert errors == []
    assert second_started.is_set()
    assert registry.snapshot()["job_id"] == "second"
    assert registry.snapshot()["active"] is True

# Original coverage: test_web.py.

@pytest.mark.usefixtures("_restore_migrated_workflow_session")
def test_subgroup_row_marks_known_unresolved_mode_identity_unselectable():
    subgroup = SubgroupInfo(
        index=4,
        space_group_number=51,
        space_group_symbol="Pmma",
        subgroup_index=36,
        size=6,
        is_maximal=False,
        opd_symbol="P2",
        k_point_label="Y",
        irrep_label="Y1",
        basis_vectors=[[0, 2, 0], [0, 0, 1], [3, 0, 0]],
        origin=["1/2", "1/2", "0"],
        k_parameters=["1/3"],
        mode_identity_status="unresolved",
        mode_identity_failures=(
            ModeIdentityFailure(
                amplitude_key="Y1_mode", orbit_id="orbit-a", global_irrep="Y1",
                k_coordinates=("1/3", "0", "0"),
                reason=(
                    "rank_mismatch; parent_SG=123; KVALUE=(1/3); "
                    "direction_selector=(a)"
                )
            ),
        ),
    )

    [row] = workflow_service._subgroup_rows([subgroup])

    assert row["selectable"] is False
    assert row["mode_identity_status"] == "unresolved"
    assert "rank_mismatch" in row["mode_identity_warning"]
    assert "parent_SG=123" in row["mode_identity_warning"]

@pytest.mark.usefixtures("_restore_migrated_workflow_session")
def test_method3_rows_keep_numbers_and_add_official_fraction_display():
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=12,
        space_group_symbol="C2/m",
        subgroup_index=4,
        size=1,
        basis_vectors=[[0, 0.5, 1], [-1, 0, 0], [0, -0.25, 0]],
        origin=[0.5, -0.25, 0],
        k_point_label="LD",
        irrep_label="GM1+",
        opd_symbol="P1",
        k_parameters=["1/6"],
    )
    item = Method3ResultItem(
        subgroup=subgroup,
        point_group="2/m",
        basis=subgroup.basis_vectors,
        routes=[subgroup],
    )

    row = workflow_service._method3_rows([item])[0]
    assert row["basis_vectors"] == [[0, 0.5, 1], [-1, 0, 0], [0, -0.25, 0]]
    assert row["origin"] == [0.5, -0.25, 0]
    assert row["basis_display"] == "{(0,1/2,1),(-1,0,0),(0,-1/4,0)}"
    assert row["origin_display"] == "(1/2,-1/4,0)"

@pytest.mark.usefixtures("_restore_migrated_workflow_session")
def test_method3_affine_only_row_is_not_selectable_or_fake_route():
    subgroup = SubgroupInfo(
        index=3,
        space_group_number=3,
        space_group_symbol="P2",
        subgroup_index=2,
        size=2,
        basis_vectors=[[0, 0, 1], [0, -1, 0], [2, 0, 0]],
        origin=[0, 0, 0],
        method3_route_resolution="affine_only_unresolved_coupled_route",
    )
    item = Method3ResultItem(
        subgroup=subgroup,
        point_group="2",
        basis=subgroup.basis_vectors,
        routes=[],
    )

    row = workflow_service._method3_rows([item])[0]
    assert row["known_routes"] == []
    assert row["known_route_count"] == 0
    assert row["route_resolution"] == "affine_only_unresolved_coupled_route"
    assert row["selectable"] is False

@pytest.mark.usefixtures("_restore_migrated_workflow_session")
def test_method3_exact_fixed_space_row_is_selectable_without_fake_single_ir_route():
    subgroup = SubgroupInfo(
        index=4,
        space_group_number=12,
        space_group_symbol="C2/m",
        subgroup_index=4,
        size=2,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 2]],
        origin=[0, 0, 0],
        method3_route_resolution="exact_fixed_space",
    )
    item = Method3ResultItem(
        subgroup=subgroup,
        point_group="2/m",
        basis=subgroup.basis_vectors,
        routes=[],
    )

    web_row = workflow_service._method3_rows([item])[0]
    assert web_row["known_routes"] == []
    assert web_row["route_resolution"] == "exact_fixed_space"
    assert web_row["selectable"] is True

@pytest.mark.usefixtures("_restore_migrated_workflow_session")
def test_method4_api_removes_daughter_upload(tmp_path, monkeypatch):
    upload_dir = tmp_path / "web_uploads"
    upload_dir.mkdir()
    uploaded = upload_dir / "daughter.cif"
    uploaded.write_text("data_test\n", encoding="utf-8")

    result = SimpleNamespace(
        amplitudes={"mode": 1.0},
        parent_cell_amplitudes={"mode": 2.0},
        raw_coefficients={"mode": 3.0},
        mode_normfactors={"mode": 4.0},
        strain_mode_amplitudes={"GM1+strain(a)": 0.1},
        strain_modes=[{
            "label": "GM1+strain(a)",
            "irrep_label": "GM1+",
            "q_raw": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "q_unit": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "normfactor": 1.0,
        }],
        strain_raw_coordinate_sum_parent_basis={
            "11": 0.1, "22": 0.2, "33": 0.3,
            "2*23": 0.0, "2*13": 0.0, "2*12": 0.0,
        },
        strain_applied_engineering_q_parent_basis={
            "11": 0.1, "22": 0.2, "33": 0.3,
            "2*23": 0.0, "2*13": 0.0, "2*12": 0.0,
        },
        strain_tensor_parent_basis=[
            [0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]
        ],
        strain_multiplier_parent_basis=[
            [1.1, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 1.3]
        ],
        rms_residual=1e-9,
        max_abs_residual=2e-9,
        metadata={"strain_convention": "test convention"},
    )
    stub = SimpleNamespace(search_method_4=lambda **_kwargs: result)
    workflow_service._SESSION._iso = stub
    monkeypatch.setattr(workflow_service, "_write_upload", lambda *_args: str(uploaded))

    payload = workflow_service._api_method4(workflow_service._SESSION, {"content": "data_test", "filename": "daughter.cif"})

    assert payload["strain_applied_engineering_q_parent_basis"]["11"] == pytest.approx(0.1)
    assert payload["strain_raw_coordinate_sum_parent_basis"]["11"] == pytest.approx(0.1)
    assert payload["strain_tensor_parent_basis"][2][2] == pytest.approx(0.3)
    assert payload["strain_mode_amplitudes"]["GM1+strain(a)"] == pytest.approx(0.1)
    assert payload["strain_voigt_engineering"] == payload[
        "strain_applied_engineering_q_parent_basis"
    ]
    assert payload["strain_tensor"] == payload["strain_tensor_parent_basis"]
    assert not uploaded.exists()

@pytest.mark.usefixtures("_restore_migrated_workflow_session")
def test_method4_api_removes_daughter_upload_after_failure(tmp_path, monkeypatch):
    upload_dir = tmp_path / "web_uploads"
    upload_dir.mkdir()
    uploaded = upload_dir / "bad-daughter.cif"
    uploaded.write_text("data_bad\n", encoding="utf-8")

    def _fail(**_kwargs):
        raise ValueError("matching failed")

    workflow_service._SESSION._iso = SimpleNamespace(search_method_4=_fail)
    monkeypatch.setattr(workflow_service, "_write_upload", lambda *_args: str(uploaded))

    with pytest.raises(ValueError, match="matching failed"):
        workflow_service._api_method4(workflow_service._SESSION, {"content": "data_bad", "filename": "bad.cif"})
    assert not uploaded.exists()

@pytest.mark.usefixtures("_restore_migrated_workflow_session")
def test_web_session_replaces_and_cleans_parent_upload(tmp_path, monkeypatch):
    upload_dir = tmp_path / "web_uploads"
    upload_dir.mkdir()
    first = upload_dir / "first.cif"
    second = upload_dir / "second.cif"
    first.write_text("first", encoding="utf-8")
    second.write_text("second", encoding="utf-8")
    session = workflow_service.WebSession()
    session._iso = _SessionIsoStub()
    session.parent_upload_path = first
    monkeypatch.setattr(workflow_service, "_SESSION", session)
    monkeypatch.setattr(workflow_service, "IsoDistort", _SessionIsoStub)

    def load(staged, _data):
        staged.parent_upload_path = second
        return {}

    monkeypatch.setitem(workflow_service._ACTIONS, "load_cif", load)
    assert first.exists()
    workflow_service.execute("load_cif", {"content": "data_second"})
    assert not first.exists()
    assert second.exists()
    workflow_service.cleanup_uploads()
    assert not second.exists()


@pytest.fixture
def _restore_migrated_workflow_session(monkeypatch):
    """Keep the original function isolation without importing HTTP lifecycle code."""
    previous = workflow_service._SESSION
    workflow_service._SESSION = workflow_service.WebSession()
    workflow_service.reset_export_progress()
    try:
        yield
    finally:
        workflow_service._SESSION.cleanup_uploads()
        # Restore test-local monkeypatches before restoring the previous owner.
        monkeypatch.undo()
        workflow_service._SESSION = previous
        workflow_service.reset_export_progress()
