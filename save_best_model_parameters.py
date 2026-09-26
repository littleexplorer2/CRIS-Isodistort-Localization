"""Save gradient-descent mode parameters to Desktop/Best_Model_Parameters.

Layout (created automatically if missing)::

    <Desktop>/Best_Model_Parameters/
      <irrep>_<structure_type>/
        <irrep>_<structure_type>_best_model_parameters.csv

CSV rows are written in mode-index order (a1, a2, …, a10), not lexicographic
string order. Folder and file names are sanitized to ``[A-Za-z0-9._-]``.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

DESKTOP_FOLDER_NAME = "Best_Model_Parameters"


def user_desktop() -> Path:
    """Return the current user's Desktop, preferring OneDrive Desktop on Windows."""
    env_home = Path(os.environ.get("USERPROFILE") or Path.home())
    candidates = [
        env_home / "OneDrive" / "Desktop",
        env_home / "Desktop",
        Path.home() / "OneDrive" / "Desktop",
        Path.home() / "Desktop",
    ]
    for path in candidates:
        if path.is_dir():
            return path
    fallback = env_home / "Desktop"
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def best_model_root() -> Path:
    """``Desktop/Best_Model_Parameters``, created if it does not exist."""
    root = user_desktop() / DESKTOP_FOLDER_NAME
    root.mkdir(parents=True, exist_ok=True)
    return root


def sanitize_group_name(*parts: str) -> str:
    tokens: list[str] = []
    for part in parts:
        text = re.sub(r"[^\w.-]+", "_", str(part).strip(), flags=re.UNICODE)
        text = text.strip("._")
        if text:
            tokens.append(text)
    return "_".join(tokens) or "unnamed"


def group_folder_name(irrep: str = "", structure_type: str = "", group: str = "") -> str:
    if str(group).strip():
        return sanitize_group_name(group)
    return sanitize_group_name(irrep, structure_type)


def csv_file_name(group: str) -> str:
    return f"{group}_best_model_parameters.csv"


def default_output_path(
    irrep: str = "",
    structure_type: str = "",
    group: str = "",
) -> Path:
    folder_name = group_folder_name(irrep, structure_type, group)
    dest_dir = best_model_root() / folder_name
    dest_dir.mkdir(parents=True, exist_ok=True)
    return dest_dir / csv_file_name(folder_name)


def _to_1d(values: Any) -> np.ndarray:
    if hasattr(values, "numpy") and callable(getattr(values, "numpy", None)):
        try:
            values = values.numpy()
        except Exception:
            pass
    return np.asarray(values, dtype=float).reshape(-1)


def _mode_sort_key(alias: str, index: int) -> tuple[int, int]:
    match = re.fullmatch(r"a(\d+)", str(alias).strip().lower())
    if match:
        return (0, int(match.group(1)))
    return (1, index)


def save_best_model_parameters(
    best_parameters: Any,
    mode_names: pd.DataFrame,
    max_mode_amps: Any,
    output_path: str | Path | None = None,
    *,
    irrep: str = "",
    structure_type: str = "",
    group: str = "",
) -> Path:
    """Save the displayed best-model parameters as a CSV file.

    The normalized amplitude matches the value printed by the notebook:
    ``best_parameters[i] / max_mode_amps[i]``.

    If ``output_path`` is omitted, the file is written to
    ``Desktop/Best_Model_Parameters/<irrep>_<structure_type>/``.
    """
    parameters = _to_1d(best_parameters)
    bounds = _to_1d(max_mode_amps)

    if len(parameters) != len(bounds):
        raise ValueError(
            "best_parameters and max_mode_amps must contain the same number of modes"
        )
    if len(mode_names) < len(parameters) or mode_names.shape[1] < 2:
        raise ValueError("mode_names must contain at least two columns for every mode")

    aliases = [f"a{i + 1}" for i in range(len(parameters))]
    names = [str(mode_names.iloc[i, 1]) for i in range(len(parameters))]
    result = pd.DataFrame(
        {
            "Mode": aliases,
            "Mode Name": names,
            "Best Model Parameter": parameters,
            "Maximum Mode Amplitude": bounds,
            "Normalized Amplitude": parameters / bounds,
        }
    )
    order = sorted(
        range(len(result)),
        key=lambda i: _mode_sort_key(result.iloc[i]["Mode"], i),
    )
    result = result.iloc[order].reset_index(drop=True)

    output_file = (
        Path(output_path) if output_path is not None else default_output_path(
            irrep=irrep,
            structure_type=structure_type,
            group=group,
        )
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_file, index=False, encoding="utf-8-sig")
    print(f"Saved best model parameters to: {output_file.resolve()}")
    return output_file


if __name__ == "__main__":
    print(
        "Import save_best_model_parameters and call it from the notebook; "
        "the notebook variables are required."
    )
