from types import SimpleNamespace

import pytest
import yaml

from isoviz_input.backend.config_loader import settings
from isoviz_input.backend.launcher import runtime as launcher


@pytest.mark.parametrize("use_default", [False, True])
def test_relocated_launchers_are_discovered_independently_of_cwd(tmp_path, monkeypatch, use_default):
    payload = yaml.safe_load(settings.CONFIG_PATH.read_text(encoding="utf-8"))
    if use_default:
        payload["isoviz"].pop("launcher_search_dirs")
    config_path = tmp_path / "workspace" / "ISOVIZ_INPUT" / "config" / "settings.yaml"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(yaml.safe_dump(payload), encoding="utf-8")
    monkeypatch.setattr(settings, "CONFIG_PATH", config_path)
    monkeypatch.setattr(settings.Config, "_instance", None)
    monkeypatch.setattr(settings.Config, "_loaded", False)
    for name in settings.get_config().isoviz_env_vars:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.chdir(tmp_path)
    root = config_path.parents[2]
    old_launcher = root / "local_tools" / "launchers" / "ISOViz.lnk"
    old_launcher.parent.mkdir(parents=True)
    old_launcher.write_bytes(b"obsolete-location")
    assert launcher.find_isoviz_launcher() is None
    current_launcher = root / "launchers" / "ISOViz.lnk"
    current_launcher.parent.mkdir()
    current_launcher.write_bytes(b"current-location")
    assert launcher.find_isoviz_launcher() == current_launcher
    override = tmp_path / "explicit.jar"
    override.write_bytes(b"explicit-choice")
    monkeypatch.setenv("ISOVIZ", str(override))
    assert launcher.find_isoviz_launcher() == override


def test_open_isoviz_jar_passes_target_to_java(tmp_path, monkeypatch):
    target = tmp_path / "prepared.isoviz"
    app = tmp_path / "IsoVIZ.jar"
    target.write_text("!isoversion 6.12\n", encoding="utf-8")
    app.write_bytes(b"jar")
    calls: list[tuple[list[str], bool]] = []

    monkeypatch.setattr(launcher, "java_executable", lambda: "javaw.exe")
    monkeypatch.setattr(
        launcher.subprocess,
        "Popen",
        lambda args, close_fds: calls.append((args, close_fds)) or SimpleNamespace(),
    )

    launcher.open_isoviz(target, launcher=app)

    assert calls == [(["javaw.exe", "-jar", str(app), str(target.resolve())], True)]


def test_open_isoviz_windows_shortcut_passes_target_explicitly(tmp_path, monkeypatch):
    target = tmp_path / "prepared.isoviz"
    app = tmp_path / "ISOViz.lnk"
    target.write_text("!isoversion 6.12\n", encoding="utf-8")
    app.write_bytes(b"shortcut")
    calls: list[tuple[list[str], bool]] = []

    monkeypatch.setattr(launcher.sys, "platform", "win32")
    monkeypatch.setattr(
        launcher.subprocess,
        "Popen",
        lambda args, close_fds: calls.append((args, close_fds)) or SimpleNamespace(),
    )

    launcher.open_isoviz(target, launcher=app)

    assert calls == [
        (["cmd", "/c", "start", "", str(app), str(target.resolve())], True)
    ]


def test_open_isoviz_executable_passes_target_explicitly(tmp_path, monkeypatch):
    target = tmp_path / "prepared.isoviz"
    app = tmp_path / "IsoVIZ.exe"
    target.write_text("!isoversion 6.12\n", encoding="utf-8")
    app.write_bytes(b"exe")
    calls: list[tuple[list[str], bool]] = []

    monkeypatch.setattr(
        launcher.subprocess,
        "Popen",
        lambda args, close_fds: calls.append((args, close_fds)) or SimpleNamespace(),
    )

    launcher.open_isoviz(target, launcher=app)

    assert calls == [([str(app), str(target.resolve())], True)]


def test_find_launcher_prefers_configured_environment(tmp_path, monkeypatch):
    app = tmp_path / "custom.jar"
    app.write_bytes(b"jar")
    monkeypatch.setenv("ISOVIZ_TEST", str(app))
    monkeypatch.setattr(
        launcher,
        "get_config",
        lambda: SimpleNamespace(
            isoviz_env_vars=("ISOVIZ_TEST",),
            launcher_search_dirs=(tmp_path,),
            launcher_names=("fallback.jar",),
        ),
    )

    assert launcher.find_isoviz_launcher() == app
