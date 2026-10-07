from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

import setup_cris


def test_all_projects_have_requirements_and_shared_version() -> None:
    expected = setup_cris.cris_version()
    assert re.fullmatch(r"\d+\.\d+\.\d+", expected)
    for project in setup_cris.PROJECTS.values():
        assert project.requirement_file.is_file()
        assert project.dev_requirement_file.is_file()
        text = (project.directory / "pyproject.toml").read_text(encoding="utf-8")
        match = re.search(r'(?m)^version\s*=\s*"([^"]+)"', text)
        assert match is not None
        assert match.group(1) == expected
        assert 'dynamic = ["dependencies"]' in text
        assert 'dependencies = { file = ["requirements.txt"] }' in text


def test_requirement_file_selection_is_stable_and_unique() -> None:
    projects = setup_cris._selected_projects("all")
    runtime = setup_cris._requirement_files(projects, dev=False)
    development = setup_cris._requirement_files(projects, dev=True)
    assert runtime == [project.requirement_file for project in projects]
    assert len(development) == 2 * len(projects)
    assert len(development) == len(set(development))


def test_recreate_target_is_scoped_to_repository() -> None:
    target = setup_cris.VENV_DIR.resolve()
    assert target.parent == Path(setup_cris.ROOT).resolve()
    assert target.name == ".venv"


def test_wsl_blocked_windows_launcher_is_replaced_after_ab_probe(
    tmp_path,
    monkeypatch,
) -> None:
    venv_root = tmp_path / ".venv"
    scripts = venv_root / "Scripts"
    scripts.mkdir(parents=True)
    python = scripts / "python.exe"
    python.write_bytes(b"venv-launcher")
    base_python = tmp_path / "base-python.exe"
    base_python.write_bytes(b"base-python")
    (venv_root / "pyvenv.cfg").write_text(
        f"executable = {base_python}\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(setup_cris, "_is_windows", lambda: True)
    monkeypatch.setattr(setup_cris.shutil, "which", lambda _name: "wsl.exe")

    def probe(candidate: Path) -> tuple[bool, str]:
        if candidate == python:
            return False, "Error code: Wsl/E_ACCESSDENIED"
        return True, ""

    monkeypatch.setattr(setup_cris, "_probe_wsl_from_python", probe)
    monkeypatch.setattr(setup_cris, "_probe_venv_prefix", lambda _path: venv_root.resolve())

    setup_cris._ensure_windows_wsl_compatible_launcher(python)

    assert python.samefile(base_python)
    assert (scripts / "python-venv-launcher.exe").read_bytes() == b"venv-launcher"


def test_windows_launcher_is_unchanged_without_specific_access_denial(
    tmp_path,
    monkeypatch,
) -> None:
    venv_root = tmp_path / ".venv"
    scripts = venv_root / "Scripts"
    scripts.mkdir(parents=True)
    python = scripts / "python.exe"
    python.write_bytes(b"venv-launcher")
    base_python = tmp_path / "base-python.exe"
    base_python.write_bytes(b"base-python")
    (venv_root / "pyvenv.cfg").write_text(
        f"executable = {base_python}\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(setup_cris, "_is_windows", lambda: True)
    monkeypatch.setattr(setup_cris.shutil, "which", lambda _name: "wsl.exe")
    monkeypatch.setattr(
        setup_cris,
        "_probe_wsl_from_python",
        lambda _path: (False, "The requested WSL distribution is not installed"),
    )

    setup_cris._ensure_windows_wsl_compatible_launcher(python)

    assert python.read_bytes() == b"venv-launcher"
    assert not (scripts / "python-venv-launcher.exe").exists()


def test_cli_exposes_install_download_and_doctor() -> None:
    parser = setup_cris.build_parser()
    assert parser.parse_args(["install"]).command == "install"
    assert parser.parse_args(["download", "--wheelhouse", "wheels"]).command == "download"
    assert parser.parse_args(["doctor"]).command == "doctor"


def test_post_install_doctor_uses_wsl_safe_launcher_on_windows(monkeypatch) -> None:
    monkeypatch.setattr(setup_cris, "_is_windows", lambda: True)
    powershell = r"C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe"
    monkeypatch.setattr(setup_cris.shutil, "which", lambda _name: powershell)
    python = setup_cris.VENV_DIR / "Scripts" / "python.exe"

    command = setup_cris._post_install_doctor_command(
        python,
        project="all",
        dev=True,
    )

    assert command == [
        powershell,
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-File",
        setup_cris.ROOT / "run_cris.ps1",
        str(Path(setup_cris.__file__).resolve()),
        "doctor",
        "--project",
        "all",
        "--dev",
    ]


@pytest.mark.skipif(os.name != "nt", reason="PowerShell launcher is Windows-only")
def test_run_cris_restores_caller_environment() -> None:
    powershell = shutil.which("powershell.exe")
    assert powershell is not None
    launcher = str(setup_cris.ROOT / "run_cris.ps1").replace("'", "''")
    probe = f"""
$env:VIRTUAL_ENV = 'outer-venv'
$env:PYTHONNOUSERSITE = 'outer-no-user-site'
$env:PYTHONPATH = 'outer-python-path'
$beforePath = $env:PATH
& '{launcher}' -c "print('child-ok')"
$childExit = $LASTEXITCODE
[ordered]@{{
    child_exit = $childExit
    virtual_env = $env:VIRTUAL_ENV
    python_no_user_site = $env:PYTHONNOUSERSITE
    python_path = $env:PYTHONPATH
    path_unchanged = ($env:PATH -ceq $beforePath)
}} | ConvertTo-Json -Compress
"""

    completed = subprocess.run(
        [
            powershell,
            "-NoProfile",
            "-ExecutionPolicy",
            "Bypass",
            "-Command",
            probe,
        ],
        cwd=setup_cris.ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    assert "child-ok" in lines
    payload = json.loads(lines[-1])
    assert payload == {
        "child_exit": 0,
        "virtual_env": "outer-venv",
        "python_no_user_site": "outer-no-user-site",
        "python_path": "outer-python-path",
        "path_unchanged": True,
    }


def test_wsl_stream_decoder_handles_utf8_and_utf16le_separately() -> None:
    assert setup_cris._decode_process_output(b"CRIS_WSL_OK\n") == "CRIS_WSL_OK"
    warning = "wsl: proxy warning"
    assert setup_cris._decode_process_output(warning.encode("utf-16-le")) == warning
    assert setup_cris._is_wsl_access_denied("Error code: Wsl/E_ACCESSDENIED")
    assert setup_cris._is_wsl_access_denied("Wsl/Service/E_ACCESSDENIED")
    assert setup_cris._is_wsl_access_denied(
        "Wsl/EnumerateDistros/Service/E_ACCESSDENIED"
    )


def test_native_doctor_distinguishes_direct_staged_and_unready_binaries(
    tmp_path, monkeypatch,
) -> None:
    monkeypatch.setattr(
        setup_cris, "_probe_json",
        lambda *_args, **_kwargs: ({
            "iso": {"ready": True, "staged": False, "executable": "/opt/isobyu/iso"},
            "findsym": {"ready": True, "staged": True, "source": "/opt/isobyu/findsym"},
            "smodes": {"ready": False, "error": "private copy is not executable"},
        }, ""),
    )
    report = setup_cris.DoctorReport(["isodistort"])

    setup_cris._check_native_executables(report, Path("python"), tmp_path)

    assert [(check.name, check.status) for check in report.checks] == [
        ("iso-executable", "pass"),
        ("findsym-executable", "warn"),
        ("smodes-executable", "fail"),
    ]
    assert "source remains unchanged" in report.checks[1].detail
    assert "private copy is not executable" in report.checks[2].detail


def test_native_doctor_probe_failure_never_reports_execute_success(
    tmp_path, monkeypatch,
) -> None:
    monkeypatch.setattr(setup_cris, "_probe_json", lambda *_args, **_kwargs: (None, "permission denied"))
    report = setup_cris.DoctorReport(["isodistort"])

    setup_cris._check_native_executables(report, Path("python"), tmp_path)

    assert len(report.checks) == 3
    assert all(check.status == "fail" for check in report.checks)
    assert all(check.detail == "permission denied" for check in report.checks)
