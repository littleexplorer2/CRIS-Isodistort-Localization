"""Native binary packaging must not mutate the read-only source files."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path, PurePosixPath
from types import SimpleNamespace

import pytest

from backend.utils import WrapperRunError
from backend.wrappers.base_wrapper import BaseWrapper


def _native_wrapper(tmp_path: Path, monkeypatch) -> BaseWrapper:
    temporary_directory = tempfile.TemporaryDirectory
    monkeypatch.setattr(
        tempfile, "TemporaryDirectory",
        lambda *, prefix, dir: temporary_directory(prefix=prefix, dir=tmp_path),
    )
    wrapper = object.__new__(BaseWrapper)
    wrapper._mode = "native"
    wrapper._init_native_environment()
    return wrapper


def _remove_source_execute_permission(source: Path, monkeypatch) -> None:
    source.chmod(0o600)
    actual_access = os.access
    # Windows does not implement POSIX execute bits; emulate only the source
    # permission result while exercising real copying and file metadata.
    monkeypatch.setattr(
        os, "access",
        lambda path, mode: False if Path(path) == source and mode == os.X_OK else actual_access(path, mode),
    )


def test_native_nonexecutable_binary_is_privately_copied_without_source_changes(
    tmp_path, monkeypatch,
) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"binary fixture")
    _remove_source_execute_permission(source, monkeypatch)
    before = (source.read_bytes(), source.stat().st_mode)
    first = _native_wrapper(tmp_path, monkeypatch)
    second = _native_wrapper(tmp_path, monkeypatch)

    executable = Path(first._wsl_bin_path(source))

    assert executable != source
    assert executable.parent == Path(first._stage_dir)
    assert first._stage_dir != second._stage_dir
    assert executable.read_bytes() == source.read_bytes()
    assert (source.read_bytes(), source.stat().st_mode) == before
    assert first._wsl_bin_path(source) == str(executable)
    if os.name == "posix":
        assert executable.stat().st_mode & 0o777 == 0o700
        assert executable.parent.stat().st_mode & 0o777 == 0o700


def test_native_executable_binary_runs_directly(tmp_path, monkeypatch) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"binary fixture")
    monkeypatch.setattr(os, "access", lambda _path, _mode: True)
    wrapper = _native_wrapper(tmp_path, monkeypatch)

    assert wrapper._wsl_bin_path(source) == str(source.resolve())
    assert not list(Path(wrapper._stage_dir).iterdir())


def test_wsl_binary_path_keeps_existing_conversion(monkeypatch) -> None:
    wrapper = object.__new__(BaseWrapper)
    wrapper._mode = "wsl"
    monkeypatch.setattr(wrapper, "_win_to_wsl", lambda _path: "/mnt/c/isobyu/iso")

    assert wrapper._wsl_bin_path(Path("iso")) == "/mnt/c/isobyu/iso"
    assert wrapper._wsl_bin_path(PurePosixPath("/home/user/iso")) == "/home/user/iso"


def test_native_staging_failure_is_reported(monkeypatch) -> None:
    def deny_staging(**_kwargs):
        raise PermissionError("private directory denied")

    monkeypatch.setattr(tempfile, "TemporaryDirectory", deny_staging)
    wrapper = object.__new__(BaseWrapper)

    with pytest.raises(WrapperRunError, match="Cannot create private staging directory"):
        wrapper._init_native_environment()


def test_native_binary_copy_is_invalidated_when_source_is_replaced(
    tmp_path, monkeypatch,
) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"old")
    _remove_source_execute_permission(source, monkeypatch)
    wrapper = _native_wrapper(tmp_path, monkeypatch)
    previous = Path(wrapper._wsl_bin_path(source))

    source.write_bytes(b"new binary")
    current = Path(wrapper._wsl_bin_path(source))

    assert current != previous
    assert current.read_bytes() == b"new binary"
    assert source.read_bytes() == b"new binary"


def test_native_binary_copy_failure_has_actionable_error(tmp_path, monkeypatch) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"binary fixture")
    _remove_source_execute_permission(source, monkeypatch)
    wrapper = _native_wrapper(tmp_path, monkeypatch)
    def deny_chmod(_path, _mode):
        raise PermissionError("blocked")

    monkeypatch.setattr(Path, "chmod", deny_chmod)

    with pytest.raises(WrapperRunError, match=r"Cannot prepare native executable.*blocked"):
        wrapper._wsl_bin_path(source)

    assert source.read_bytes() == b"binary fixture"
    assert not list(Path(wrapper._stage_dir).iterdir())


@pytest.mark.skipif(os.name != "posix", reason="requires native POSIX shell execution")
def test_native_run_executes_private_copy_with_isodata_and_stdin(tmp_path) -> None:
    source = tmp_path / "iso"
    source.write_text('#!/bin/sh\nprintf "%s\\n" "$ISODATA"\ncat\n', encoding="utf-8")
    source.chmod(0o600)
    data = tmp_path / "data"
    data.mkdir()
    temporary_inputs = tmp_path / "inputs"
    temporary_inputs.mkdir()
    wrapper = object.__new__(BaseWrapper)
    wrapper._mode = "native"
    wrapper.cfg = SimpleNamespace(
        timeout=5, temp_dir=temporary_inputs,
        _cfg={"isobyu": {"data_dir": str(data)}},
        resolve_path=Path,
    )
    wrapper._init_native_environment()
    before = source.stat().st_mode

    assert wrapper.run_stdin(source, "input text\n").splitlines() == [str(data) + "/", "input text"]
    assert source.stat().st_mode == before
