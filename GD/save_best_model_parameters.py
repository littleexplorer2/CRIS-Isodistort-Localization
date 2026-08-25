"""Save gradient-descent mode parameters to CSV.

Set OUTPUT_CSV_PATH before running the notebook, or pass output_path explicitly
when calling save_best_model_parameters().
"""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


# Change this path before running the notebook.
OUTPUT_CSV_PATH = Path("best_model_parameters.csv")


def save_best_model_parameters(
    best_parameters: Any,
    mode_names: pd.DataFrame,
    max_mode_amps: Any,
    output_path: str | Path = OUTPUT_CSV_PATH,
) -> Path:
    """Save the displayed best-model parameters as a CSV file.

    The normalized amplitude matches the value printed by the notebook:
    ``best_parameters[i] / max_mode_amps[i]``.
    """
    parameters = np.asarray(best_parameters).reshape(-1)
    bounds = np.asarray(max_mode_amps).reshape(-1)

    if len(parameters) != len(bounds):
        raise ValueError(
            "best_parameters and max_mode_amps must contain the same number of modes"
        )
    if len(mode_names) < len(parameters) or mode_names.shape[1] < 2:
        raise ValueError("mode_names must contain at least two columns for every mode")

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    result = pd.DataFrame(
        {
            "Mode": [f"a{i + 1}" for i in range(len(parameters))],
            "Mode Name": [str(mode_names.iloc[i, 1]) for i in range(len(parameters))],
            "Best Model Parameter": parameters,
            "Maximum Mode Amplitude": bounds,
            "Normalized Amplitude": parameters / bounds,
        }
    )
    result.to_csv(output_file, index=False, encoding="utf-8-sig")
    print(f"Saved best model parameters to: {output_file.resolve()}")
    return output_file


if __name__ == "__main__":
    print(
        "Import save_best_model_parameters and call it from the notebook; "
        "the notebook variables are required."
    )
