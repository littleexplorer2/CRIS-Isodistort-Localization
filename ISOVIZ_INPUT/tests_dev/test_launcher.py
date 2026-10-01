from types import SimpleNamespace

from isoviz_input import launcher


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
