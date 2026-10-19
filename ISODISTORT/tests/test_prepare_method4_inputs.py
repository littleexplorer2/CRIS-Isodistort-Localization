"""Tests for test_prepare_method4_inputs; migrated without changing scientific expectations."""
from __future__ import annotations

from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests import prepare_method4_inputs as preparation
import numpy as np


# Original coverage: test_method4_official_audit.py.

def test_dmax_scaling_probes_inside_unwrapped_linear_regime(monkeypatch) -> None:
    lattice = Lattice.cubic(10.0)
    zero = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])

    def wrapped_generator(_iso, contributions):
        coefficient = contributions.get("mode", 0.0)
        return Structure(
            lattice,
            ["H"],
            [[0.0, 0.0, (2.0 * coefficient) % 1.0]],
        )

    monkeypatch.setattr(preparation, "_generated", wrapped_generator)
    contributions, generated = preparation._scaled_contributions_for_dmax(
        object(), zero, {"mode": 1.0}, 0.30,
    )

    assert np.isclose(contributions["mode"], 0.015)
    displacement = np.linalg.norm(generated.cart_coords[0] - zero.cart_coords[0])
    assert np.isclose(displacement, 0.30)
