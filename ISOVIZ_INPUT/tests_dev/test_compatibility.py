"""Legacy imports must share state with their canonical implementations."""
from __future__ import annotations

import importlib
from types import SimpleNamespace

import pytest


@pytest.mark.parametrize("legacy,canonical,public_names", [
    ("amplitudes", "backend.amplitudes.core", (
        "ModeAmplitude", "PatchReport", "apply_amplitudes", "list_mode_headers",
        "parse_amplitude_csv", "patch_isoviz_file", "read_amplitude_csv",
    )),
    ("config_loader", "backend.config_loader.settings", (
        "CONFIG_PATH", "PROJECT_ROOT", "Config", "get_config",
    )),
    ("launcher", "backend.launcher.runtime", (
        "cris_root", "find_isoviz_launcher", "java_executable", "open_isoviz",
    )),
    ("paths", "backend.paths.resolver", (
        "DATA_DIR", "DESKTOP_FOLDER_NAME", "INPUT_ROOT", "ISOVIZ_ROOT", "STRUCTURE_DIR",
        "best_model_root", "ensure_best_model_root", "ensure_input_content",
        "list_csv_files", "list_named_subfolders", "resolve_best_model_csv",
        "resolve_best_model_folder", "resolve_isoviz_path", "strip_user_path", "user_desktop",
    )),
    ("validation", "backend.validation.checks", ("section_scalar", "validate_patched_text")),
])
def test_legacy_modules_share_identity_and_exports(legacy, canonical, public_names):
    original = importlib.import_module(f"isoviz_input.{legacy}")
    current = importlib.import_module(f"isoviz_input.{canonical}")
    assert original is current
    assert all(hasattr(original, name) for name in public_names)


def test_legacy_path_override_controls_directory_creation(tmp_path, monkeypatch):
    legacy = importlib.import_module("isoviz_input.paths")
    inputs = tmp_path / "legacy_inputs"
    monkeypatch.setattr(legacy, "INPUT_ROOT", inputs)
    monkeypatch.setattr(legacy, "DATA_DIR", inputs / "data")
    monkeypatch.setattr(legacy, "STRUCTURE_DIR", inputs / "structure")

    data, structure = legacy.ensure_input_content()

    assert data == inputs / "data"
    assert structure == inputs / "structure"
    assert data.is_dir() and structure.is_dir()


def test_legacy_configuration_override_controls_loading(tmp_path, monkeypatch):
    legacy = importlib.import_module("isoviz_input.config_loader")
    settings = tmp_path / "settings.yaml"
    settings.write_text("paths:\n  best_model_folder_name: Legacy_Model_Directory\n", encoding="utf-8")
    monkeypatch.setattr(legacy, "CONFIG_PATH", settings)
    monkeypatch.setattr(legacy.Config, "_instance", None)
    monkeypatch.setattr(legacy.Config, "_loaded", False)

    configuration = legacy.get_config()
    assert configuration.best_model_folder_name == "Legacy_Model_Directory"
    assert configuration.resolve_path("relative_inputs") == settings.parent / "relative_inputs"


def test_bundled_settings_preserve_logical_source_layout(tmp_path, monkeypatch):
    legacy = importlib.import_module("isoviz_input.config_loader")
    configuration = legacy.get_config()
    installed_root = tmp_path / "installed"
    bundled_settings = installed_root / "isoviz_input" / "_resources" / "settings.yaml"
    monkeypatch.setattr(legacy, "_SOURCE_CONFIG_PATH", installed_root / "config" / "settings.yaml")
    monkeypatch.setattr(legacy, "_PACKAGED_CONFIG_PATH", bundled_settings)
    monkeypatch.setattr(legacy, "CONFIG_PATH", bundled_settings)

    assert configuration.resolve_path("../input_content") == installed_root / "input_content"
    assert configuration.resolve_path(str(tmp_path)) == tmp_path


def test_legacy_java_override_controls_jar_dispatch(tmp_path, monkeypatch):
    legacy = importlib.import_module("isoviz_input.launcher")
    target = tmp_path / "prepared.isoviz"
    launcher = tmp_path / "IsoVIZ.jar"
    target.write_text("!isoversion 6.12\n", encoding="utf-8")
    launcher.write_bytes(b"test launcher placeholder")
    calls = []
    monkeypatch.setattr(legacy, "java_executable", lambda: "legacy-test-java.exe")
    monkeypatch.setattr(
        legacy.subprocess,
        "Popen",
        lambda arguments, close_fds: calls.append((arguments, close_fds)) or SimpleNamespace(),
    )

    legacy.open_isoviz(target, launcher=launcher)

    assert calls == [(["legacy-test-java.exe", "-jar", str(launcher), str(target.resolve())], True)]
