"""Tests for support_backend_method4; migrated without changing scientific expectations."""
from __future__ import annotations

from pymatgen.core import Lattice
from pymatgen.core import Structure
import numpy as np


# Original coverage: test_real_mode_fallback.py.

def _parent(centered=False, origin=(0.0, 0.0, 0.0)):
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]] if centered else [[0.0, 0.0, 0.0]]
    return Structure(
        Lattice.cubic(4.0), ["Fe"] * len(coords),
        np.asarray(coords) + np.asarray(origin), to_unit_cell=False,
    )
