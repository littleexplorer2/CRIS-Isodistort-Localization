from __future__ import annotations

import re
from pathlib import Path

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


def test_wsl_stream_decoder_handles_utf8_and_utf16le_separately() -> None:
    assert setup_cris._decode_process_output(b"CRIS_WSL_OK\n") == "CRIS_WSL_OK"
    warning = "wsl: proxy warning"
    assert setup_cris._decode_process_output(warning.encode("utf-16-le")) == warning
