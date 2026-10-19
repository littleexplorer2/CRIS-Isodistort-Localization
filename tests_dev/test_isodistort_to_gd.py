"""Format checks for isodistort_to_gd writers (no WSL / TensorFlow)."""
from __future__ import annotations

import json
import sys
import types
from pathlib import Path

import numpy as np
import pytest

from isodistort_to_gd import (
    GdBundle,
    _load_isodistort_components,
    check_bundle,
    default_cris_root,
    render_alris_functions,
    write_gd_bundle,
    _hkl_transform_from_basis,
    _parent_cell_count,
    bundle_from_displacement_snapshot,
)


@pytest.mark.parametrize("basis", [
    [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
    [[2, 0, 1], [0, -1, 0], [0, 0, 3]],
])
def test_hkl_dual_map_preserves_phase_for_nonsymmetric_basis(basis):
    basis = np.asarray(basis, dtype=float)
    parent_h = np.array([[2, -1, 3], [1, 0, -2]], dtype=float)
    child_x = np.array([[0.1, 0.2, 0.3], [0.3, 0.1, 0.4]], dtype=float)
    transform = _hkl_transform_from_basis(basis)
    parent_x = child_x @ basis
    np.testing.assert_allclose(parent_h @ parent_x.T, (parent_h @ transform) @ child_x.T,
                               rtol=0, atol=1e-14)
    assert not np.allclose(parent_h @ parent_x.T, (parent_h @ basis) @ child_x.T)


@pytest.mark.parametrize("basis", [np.eye(2), np.zeros((3, 3)), np.full((3, 3), np.nan)])
def test_invalid_basis_is_rejected_instead_of_replaced_by_identity(basis):
    with pytest.raises(ValueError):
        _hkl_transform_from_basis(basis)


def test_parent_cell_index_is_exact_and_not_rounded_or_clamped():
    assert _parent_cell_count([[2, 1, 0], [0, -3, 0], [0, 0, 1]]) == 6
    with pytest.raises(NotImplementedError, match="integer"):
        _parent_cell_count([["1/2", 0, 0], [0, 1, 0], [0, 0, 1]])
    with pytest.raises(ValueError):
        _parent_cell_count(np.zeros((3, 3)))


def test_default_cris_root_uses_environment_or_sibling(tmp_path, monkeypatch):
    configured = tmp_path / "configured-cris"
    monkeypatch.setenv("CRIS_ROOT", str(configured))
    assert default_cris_root(tmp_path / "GD") == configured.resolve()
    monkeypatch.delenv("CRIS_ROOT")
    assert default_cris_root(tmp_path / "GD") == tmp_path / "CRIS"


def test_runtime_loader_uses_restructured_cris_modules(tmp_path, monkeypatch):
    (tmp_path / "CRIS" / "ISODISTORT").mkdir(parents=True)
    iso_distort = object()

    def normalize_modes(*_args, **_kwargs):
        return None

    modules = {
        "backend": types.ModuleType("backend"),
        "application.session": types.ModuleType("application.session"),
        "application": types.ModuleType("application"),
        "backend.export": types.ModuleType("backend.export"),
        "backend.export.distortion_formats": types.ModuleType(
            "backend.export.distortion_formats"
        ),
    }
    modules["application.session"].IsoDistort = iso_distort
    modules[
        "backend.export.distortion_formats"
    ].cart_normalized_mode_matrix = normalize_modes
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)

    loaded_iso_distort, loaded_normalizer = _load_isodistort_components(
        tmp_path / "CRIS"
    )

    assert loaded_iso_distort is iso_distort
    assert loaded_normalizer is normalize_modes


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


@pytest.mark.parametrize("lattice_scale", [1e40, 1e-40])
def test_snapshot_conversion_rejects_float32_overflow(lattice_scale):
    snapshot = types.SimpleNamespace(
        schema="isodistort.displacement-snapshot.v1",
        coordinate_convention="row_fractional_basis_supercell_before_origin_shift",
        normalization="raw_calculated_supercell_columns",
        candidate_identity_json=json.dumps({
            "k_point_label": "GM", "irrep_label": "GM1", "opd_symbol": "P1", "k_parameters": [],
            "basis_vectors": np.eye(3, dtype=int).tolist(), "origin": [0, 0, 0], "space_group_number": 1,
        }),
        nmod=0, species=("H",), fractional_coordinates=((0., 0., 0.),),
        lattice_angstrom=lattice_scale * np.eye(3),
        parent_space_group_number=1, parent_space_group_symbol="P1", child_space_group_symbol="P1",
        columns=(types.SimpleNamespace(key="one", label="one", fractional_vectors=((1., 0., 0.),)),),
    )

    # Independent closed-form normalization for one vector on a cubic lattice.
    def normalize(raw, lattice, *, centering_mult):
        assert centering_mult == 1
        return raw / lattice[0, 0], 1.0

    with pytest.raises(ValueError, match="float32"):
        bundle_from_displacement_snapshot(snapshot, normalize, cif=Path("independent.cif"))


@pytest.mark.parametrize("field", ["lattice", "base_frac", "hkl_transform", "norm_factors", "max_bounds"])
def test_check_bundle_rejects_nonfinite_output_arrays(field):
    bundle = _toy_bundle()
    values = getattr(bundle, field).copy()
    values.flat[0] = np.nan
    setattr(bundle, field, values)
    assert check_bundle(bundle)


@pytest.mark.parametrize("field", ["mode_frac", "base_frac", "lattice", "hkl_transform", "norm_factors", "max_bounds"])
def test_check_bundle_reports_malformed_shapes_without_indexing_them(field):
    bundle = _toy_bundle()
    setattr(bundle, field, np.array(1.0))
    assert any(field + " shape" in error for error in check_bundle(bundle))


@pytest.mark.parametrize("field", ["norm_factors", "max_bounds"])
@pytest.mark.parametrize("value", [0.0, -1.0, 1e-50])
def test_check_bundle_requires_positive_representable_amplitude_parameters(field, value):
    bundle = _toy_bundle()
    values = getattr(bundle, field).astype(float)
    values[-1] = value
    setattr(bundle, field, values)
    assert any(field + " must be positive" in error for error in check_bundle(bundle))


def test_check_bundle_rejects_a_nonfirst_mode_that_underflows_to_zero():
    bundle = _toy_bundle()
    bundle.mode_frac = bundle.mode_frac.astype(float)
    bundle.mode_frac[-1] *= 1e-50
    assert "mode vector 1 is identically zero in float32" in check_bundle(bundle)


@pytest.mark.parametrize("field", ["lattice", "hkl_transform"])
def test_check_bundle_rejects_singular_stored_coordinate_maps(field):
    bundle = _toy_bundle()
    setattr(bundle, field, np.diag([1., 1., 0.]))
    assert any("nonsingular" in error for error in check_bundle(bundle))


def test_check_bundle_rejects_a_finite_lattice_with_unrepresentable_reciprocal():
    bundle = _toy_bundle()
    bundle.lattice = np.eye(3, dtype=np.float32) * np.float32(1e-40)
    assert "reciprocal lattice must be finite in float32" in check_bundle(bundle)


def test_invalid_bundle_is_rejected_before_any_file_is_written(tmp_path):
    bundle = _toy_bundle()
    bundle.lattice[0, 0] = np.inf
    output = tmp_path / "output"
    output.mkdir()
    previous = output / "metadata.json"
    previous.write_text("previous output", encoding="utf-8")
    with pytest.raises(ValueError):
        write_gd_bundle(output, "invalid", bundle, install_alris=tmp_path / "installed.py")
    assert list(output.iterdir()) == [previous]
    assert previous.read_text(encoding="utf-8") == "previous output"
    assert not (tmp_path / "installed.py").exists()


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
    metadata = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["maxamp_convention"] == (
        "inverse maximum Cartesian displacement at As=1; "
        "bound times dmax equals 1 Angstrom"
    )
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
