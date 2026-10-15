"""Application use cases remain callable without starting or importing a UI."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from application import workflows
from application.workflows import service


def test_workflows_import_and_execute_without_frontend(tmp_path):
    project = Path(__file__).resolve().parents[1]
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
        structure=object(), export_subgroups_zip=export_subgroups_zip,
    )
    session.method2_subgroups = [SimpleNamespace(index=0)]
    monkeypatch.setattr(service, "_SESSION", session)
    monkeypatch.setattr(service, "_EXPORT_PROGRESS", service._initial_export_progress())
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
