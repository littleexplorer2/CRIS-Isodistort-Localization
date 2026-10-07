from __future__ import annotations

from pathlib import Path

import save_best_model_parameters as saver


def test_best_model_root_prefers_environment_override(tmp_path, monkeypatch):
    target = tmp_path / "Best_Model_Parameters"
    monkeypatch.setenv(saver.BEST_MODEL_ENV_VAR, str(target))
    assert saver.best_model_root() == target.resolve()
    assert target.is_dir()


def test_best_model_root_uses_existing_sibling(tmp_path, monkeypatch):
    sibling = tmp_path / "Best_Model_Parameters"
    sibling.mkdir()
    monkeypatch.delenv(saver.BEST_MODEL_ENV_VAR, raising=False)
    monkeypatch.setattr(saver, "GD_ROOT", tmp_path / "GD")
    assert saver.best_model_root() == sibling


def test_default_output_path_stays_inside_selected_root(tmp_path, monkeypatch):
    monkeypatch.setenv(saver.BEST_MODEL_ENV_VAR, str(tmp_path))
    output = saver.default_output_path(irrep="LD1", structure_type="C1")
    assert output == Path(tmp_path) / "LD1_C1" / "LD1_C1_best_model_parameters.csv"
