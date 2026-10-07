from __future__ import annotations

from types import SimpleNamespace

import pytest

from isoviz_input.paths import (
    best_model_root,
    resolve_best_model_csv,
    resolve_isoviz_path,
    strip_user_path,
)


def test_best_model_root_prefers_environment_override(tmp_path, monkeypatch):
    target = tmp_path / "amplitudes"
    monkeypatch.setenv("BEST_MODEL_PARAMETERS_DIR", str(target))
    monkeypatch.setattr(
        "isoviz_input.paths.get_config",
        lambda: SimpleNamespace(
            best_model_env_var="BEST_MODEL_PARAMETERS_DIR",
            best_model_parent=None,
            best_model_folder_name="Best_Model_Parameters",
            project_root=tmp_path / "CRIS" / "ISOVIZ_INPUT",
        ),
    )
    assert best_model_root() == target.resolve()


def test_best_model_root_uses_sibling_checkout_layout(tmp_path, monkeypatch):
    project_root = tmp_path / "CRIS" / "ISOVIZ_INPUT"
    sibling = tmp_path / "Best_Model_Parameters"
    sibling.mkdir()
    monkeypatch.delenv("BEST_MODEL_PARAMETERS_DIR", raising=False)
    monkeypatch.setattr(
        "isoviz_input.paths.get_config",
        lambda: SimpleNamespace(
            best_model_env_var="BEST_MODEL_PARAMETERS_DIR",
            best_model_parent=None,
            best_model_folder_name="Best_Model_Parameters",
            project_root=project_root,
        ),
    )
    assert best_model_root() == sibling


def test_strip_user_path_ascii_and_curly_quotes():
    assert strip_user_path('"C:\\data\\LD1_C1.isoviz"') == "C:\\data\\LD1_C1.isoviz"
    assert strip_user_path("'C:\\data\\file.csv'") == "C:\\data\\file.csv"
    assert strip_user_path("“C:\\data\\LD1_C1.isoviz”") == "C:\\data\\LD1_C1.isoviz"
    assert strip_user_path("  LD1_C1  ") == "LD1_C1"


def test_resolve_best_model_csv_from_names_and_indices(tmp_path):
    folder = tmp_path / "LD1_C1"
    folder.mkdir()
    csv_path = folder / "LD1_C1_best_model_parameters.csv"
    csv_path.write_text("Mode,Mode Name,Best Model Parameter\na1,x,0.1\n", encoding="utf-8")

    found = resolve_best_model_csv('"LD1_C1"', "LD1_C1_best_model_parameters.csv", root=tmp_path)
    assert found == csv_path
    assert resolve_best_model_csv("LD1_C1", "LD1_C1_best_model_parameters", root=tmp_path) == csv_path
    assert resolve_best_model_csv("1", "1", root=tmp_path) == csv_path


def test_resolve_isoviz_path_strips_quotes(tmp_path):
    target = tmp_path / "crystal.isoviz"
    target.write_text("!isoversion 6.12\n", encoding="utf-8")
    quoted = f'"{target}"'
    assert resolve_isoviz_path(quoted) == target
    other = tmp_path / "notes.txt"
    other.write_text("nope\n", encoding="utf-8")
    with pytest.raises(ValueError):
        resolve_isoviz_path(str(other))
