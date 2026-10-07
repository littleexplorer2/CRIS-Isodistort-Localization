"""Writable runtime storage for source trees and read-only wheel installations."""

from __future__ import annotations

from pathlib import Path

import pytest

from backend.utils import config_loader


@pytest.fixture
def layout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    package_root = tmp_path / "arbitrarily-named-installation"
    config_path = package_root / "resources" / "config" / "settings.yaml"
    config_path.parent.mkdir(parents=True)
    monkeypatch.setattr(config_loader, "PROJECT_ROOT", package_root)
    monkeypatch.setattr(config_loader, "CONFIG_PATH", config_path)
    monkeypatch.delenv("ISODISTORT_RUNTIME_ROOT", raising=False)
    monkeypatch.setattr(config_loader.sys, "platform", "win32")
    user_data = tmp_path / "user-data"
    monkeypatch.setenv("LOCALAPPDATA", str(user_data))
    cfg = object.__new__(config_loader.Config)
    cfg._cfg = {
        "runtime": {"output_dir": "../../output", "temp_dir": "../../output/tmp"},
        "isobyu": {"bin_dir": "../isobyu", "iso_bin": "iso"},
    }
    return cfg, package_root, user_data


def _mark_source_tree(package_root: Path) -> None:
    (package_root / "pyproject.toml").touch()
    (package_root / "scripts").mkdir()
    (package_root / "scripts" / "main_web.py").touch()


def test_source_checkout_keeps_project_output_default(layout) -> None:
    cfg, package_root, _ = layout
    _mark_source_tree(package_root)

    assert cfg.runtime_root == package_root
    assert cfg.output_dir == package_root / "output"
    assert cfg.temp_dir == package_root / "output" / "tmp"


def test_installed_read_only_layout_uses_user_state(layout, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg, package_root, user_data = layout
    original_mkdir = Path.mkdir

    def read_only_package_mkdir(path: Path, *args, **kwargs) -> None:
        if path == package_root or package_root in path.parents:
            raise PermissionError("the installed package is read-only")
        original_mkdir(path, *args, **kwargs)

    monkeypatch.setattr(Path, "mkdir", read_only_package_mkdir)

    assert cfg.output_dir == user_data / "ISODISTORT" / "output"
    assert cfg.temp_dir == user_data / "ISODISTORT" / "output" / "tmp"
    assert cfg.output_dir.is_dir()
    assert cfg.temp_dir.is_dir()
    assert cfg.iso_bin == package_root / "resources" / "isobyu" / "iso"
    assert not (package_root / "output").exists()


@pytest.mark.parametrize("source_checkout", [False, True])
def test_runtime_root_override_rebases_relative_paths_only(
    layout, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source_checkout: bool
) -> None:
    cfg, package_root, _ = layout
    if source_checkout:
        _mark_source_tree(package_root)
    override = tmp_path / "external-state"
    monkeypatch.setenv("ISODISTORT_RUNTIME_ROOT", str(override))

    assert cfg.output_dir == override / "output"
    assert cfg.temp_dir == override / "output" / "tmp"
    assert cfg.iso_bin == package_root / "resources" / "isobyu" / "iso"
    explicit_output = tmp_path / "explicit-output"
    cfg._cfg["runtime"]["output_dir"] = str(explicit_output)
    assert cfg.output_dir == explicit_output


def test_relative_override_is_rejected_before_directory_creation(
    layout, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg, package_root, _ = layout
    monkeypatch.setenv("ISODISTORT_RUNTIME_ROOT", "relative-state")
    with pytest.raises(ValueError, match="must be an absolute path"):
        _ = cfg.output_dir
    assert not (package_root / "output").exists()


def test_installed_relative_runtime_path_cannot_escape_root(layout, tmp_path: Path) -> None:
    cfg, _, _ = layout
    with pytest.raises(ValueError, match="use an absolute path"):
        cfg.resolve_runtime_path("../../../outside")
    assert cfg.resolve_runtime_path(str(tmp_path / "outside")) == tmp_path / "outside"


@pytest.mark.parametrize(
    ("platform", "expected_suffix"),
    [
        ("win32", "AppData/Local/ISODISTORT"),
        ("darwin", "Library/Application Support/ISODISTORT"),
        ("linux", ".local/state/isodistort"),
    ],
)
def test_user_state_platform_fallbacks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, platform: str, expected_suffix: str
) -> None:
    monkeypatch.setattr(config_loader.sys, "platform", platform)
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.delenv("LOCALAPPDATA", raising=False)
    monkeypatch.delenv("XDG_STATE_HOME", raising=False)
    assert config_loader._user_runtime_root() == tmp_path / expected_suffix


def test_xdg_state_home_is_used_only_when_absolute(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(config_loader.sys, "platform", "linux")
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "custom-state"))
    assert config_loader._user_runtime_root() == tmp_path / "custom-state" / "isodistort"
    monkeypatch.setenv("XDG_STATE_HOME", "relative-state")
    assert config_loader._user_runtime_root() == tmp_path / ".local" / "state" / "isodistort"
