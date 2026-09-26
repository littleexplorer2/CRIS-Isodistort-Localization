"""Format checks for isodistort_to_gd writers (no WSL / TensorFlow)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from isodistort_to_gd import (
    GdBundle,
    check_bundle,
    render_alris_functions,
    write_gd_bundle,
)


def _toy_bundle() -> GdBundle:
    names = [
        "I4/mmm[0,0,1/6]LD1(a,b)[Eu1:a:dsp]A2u(a)",
        "I4/mmm[0,0,0]GM3-(a)[Eu1:a:dsp]A2u(a)",
    ]
    base = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], dtype=np.float32)
    modes = np.zeros((2, 2, 3), dtype=np.float32)
    modes[0, 0, 2] = 0.01
    modes[1, 1, 2] = 0.02
    return GdBundle(
        names=names,
        norm_factors=np.array([1.0, 1.0], dtype=np.float32),
        max_bounds=np.array([2.4494897, 2.4494897], dtype=np.float32),
        species=["Eu", "Al"],
        base_frac=base,
        mode_frac=modes,
        lattice=np.diag([4.4, 4.4, 67.0]).astype(np.float32),
        hkl_transform=np.diag([1.0, 1.0, 6.0]).astype(np.float32),
        parent_symbol="I4/mmm",
        irrep="LD1",
        opd="C1",
        metadata={"n_modes": 2, "n_atoms": 2},
    )


def test_write_gd_bundle_notebook_contract(tmp_path):
    bundle = _toy_bundle()
    assert check_bundle(bundle) == []
    written = write_gd_bundle(tmp_path, "LD1_C1", bundle)
    names = written["mario_mode_names"].read_text(encoding="utf-8").strip().splitlines()
    assert len(names) == 2
    assert all(names)
    norms = np.loadtxt(written["norm_factors"], delimiter=",")
    bounds = np.loadtxt(written["max_bound_vectors"], delimiter=",")
    assert norms.shape == (2,)
    assert bounds.shape == (2,)
    rows = written["displacive_modes"].read_text(encoding="utf-8").strip().splitlines()
    assert rows[0].split()[0] == "1"
    assert rows[0].split()[1] == names[0]
    alris = written["alris_functions"].read_text(encoding="utf-8")
    for token in (
        "def transform_list_hkl_p63_p65",
        "def atom_position_list",
        "def get_structure_factors",
        "def get_structure_factors_CO",
        "N_MODES = 2",
    ):
        assert token in alris


def test_generated_ld1_c1_matches_notebook_contract_if_present():
    root = Path(__file__).resolve().parents[1] / "generated" / "LD1_C1"
    prefix = root / "LD1_C1"
    names_path = Path(f"{prefix}_mario_mode_names.txt")
    if not names_path.is_file():
        return
    names = [ln.strip() for ln in names_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert names
    norms = np.loadtxt(f"{prefix}_norm_factors.txt", delimiter=",")
    bounds = np.loadtxt(f"{prefix}_max_bound_vectors.txt", delimiter=",")
    rows = Path(f"{prefix}_displacive_modes.txt").read_text(encoding="utf-8").strip().splitlines()
    assert len(norms) == len(names)
    assert len(bounds) == len(names)
    assert len(rows) == len(names)
    assert rows[0].split()[1] == names[0]
    alris = (root / "LD1_C1_alris_functions.py").read_text(encoding="utf-8")
    assert f"N_MODES = {len(names)}" in alris
    assert "def atom_position_list" in alris


def test_render_alris_mentions_child_hkl():
    bundle = _toy_bundle()
    text = render_alris_functions(
        names=bundle.names,
        species=bundle.species,
        base_frac=bundle.base_frac,
        mode_frac=bundle.mode_frac,
        lattice=bundle.lattice,
        hkl_transform=bundle.hkl_transform,
        parent_symbol=bundle.parent_symbol,
        irrep=bundle.irrep,
        opd=bundle.opd,
    )
    assert "hkl_child = hkl_parent @ HKL_TRANSFORM" in text
    compile(text, "LD1_C1_alris_functions.py", "exec")
