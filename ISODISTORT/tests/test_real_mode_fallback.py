"""Scientific boundary of the single-real-parent-column fallback."""

from __future__ import annotations

from fractions import Fraction
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from backend.api import IsoDistort
from backend.utils.lattice import (
    centering_primitive_matrix,
    rational_vector,
    real_translation_character_sign,
    wavevector_translation_phases,
)
from features.input_cif import build_supercell
from features.method4 import DistortionEngine


def _parent(centered=False, origin=(0.0, 0.0, 0.0)):
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]] if centered else [[0.0, 0.0, 0.0]]
    return Structure(
        Lattice.cubic(4.0), ["Fe"] * len(coords),
        np.asarray(coords) + np.asarray(origin), to_unit_cell=False,
    )


def _run_column_operation(operation, parent, values, k, basis=None):
    engine = DistortionEngine()
    basis = [1, 1, 1] if basis is None else basis
    if operation == "lift":
        reference, columns = engine.lift_mode_displacements(
            parent, basis, {"mode": values}, k_vector=k,
        )
        return reference, columns["mode"]
    if operation == "single":
        result = engine.generate_single_mode(parent, values, 1.0, basis, k_vector=k)
    else:
        result = engine.generate_modes(parent, basis, parent_displacements=values, k_vector=k)
    reference = build_supercell(parent, basis)
    observed = result.frac_coords - reference.frac_coords
    observed -= np.round(observed)
    return reference, observed


@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("k", [[0, 0, 0], [2, 0, 0]])
def test_explicit_primitive_gamma_requires_centering_covariance(operation, k):
    parent = _parent(centered=True)
    opposite = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    with pytest.raises(ValueError, match=r"does not prove.*translation character"):
        _run_column_operation(operation, parent, opposite, k)

    same = np.abs(opposite)
    _reference, observed = _run_column_operation(operation, parent, same, k)
    assert np.allclose(observed, same, atol=1e-14, rtol=0.0)


@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("basis", [[1, 1, 1], [2, 1, 1]])
def test_unspecified_wavevector_keeps_the_undeclared_parent_column_contract(operation, basis):
    parent = _parent(centered=True)
    values = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    reference, observed = _run_column_operation(operation, parent, values, None, basis)
    parent_positions = reference.frac_coords @ np.diag(basis)
    expected = np.zeros_like(observed)
    expected[:, 0] = np.where(np.isclose(parent_positions[:, 0] % 1.0, 0.5), -0.01, 0.01)
    assert np.allclose(observed @ np.diag(basis), expected, atol=1e-14, rtol=0.0)


@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
def test_unspecified_column_cannot_lose_a_distinct_vector_on_primitive_reduction(operation):
    parent = _parent(centered=True)
    basis = np.asarray(centering_primitive_matrix("I"), dtype=float)
    opposite = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    with pytest.raises(ValueError, match="displacement column is not periodic"):
        _run_column_operation(operation, parent, opposite, None, basis)


@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("k", [None, [0, 0, Fraction(1, 2)]])
def test_child_translation_lattice_must_be_a_parent_sublattice(operation, k):
    # k·B is integral and det(B)=1, but a/2 is not a translation of this P parent.
    basis = np.diag([0.5, 1, 2])
    with pytest.raises(ValueError, match="not a sublattice"):
        _run_column_operation(operation, _parent(), np.asarray([[0.01, 0, 0]]), k, basis)


@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("k", [None, [0, 0, 0]])
def test_valid_centered_primitive_child_preserves_atom_density_and_vectors(operation, k):
    parent = _parent(centered=True)
    basis = np.asarray(centering_primitive_matrix("I"), dtype=float)
    values = np.asarray([[0.01, 0, 0], [0.01, 0, 0]])
    reference, observed = _run_column_operation(operation, parent, values, k, basis)
    assert len(reference) == 1
    assert reference.volume == pytest.approx(parent.volume / 2)
    assert np.allclose(
        observed @ reference.lattice.matrix,
        values[:1] @ parent.lattice.matrix,
        atol=1e-14, rtol=0.0,
    )


@pytest.mark.parametrize("component", [
    1e-13, -1e-13, 1.0 + 1e-13, 1.0 - 1e-13,
    0.5 + 1e-13, 0.5 - 1e-13,
])
def test_float_wavevector_is_not_rounded_across_a_real_phase_boundary(component):
    k = rational_vector([0, 0, component])
    assert k[2] == Fraction(str(component))
    assert (2 * k[2]).denominator != 1
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0.0, 0.0]]), 1.0, [1, 1, 2],
            k_vector=[0, 0, component],
        )


def test_centered_reciprocal_character_is_checked_on_primitive_translations():
    phases = wavevector_translation_phases(
        [Fraction(1, 2), 0, 0], centering_primitive_matrix("I"),
    )
    assert set(phases) == {Fraction(-1, 4), Fraction(1, 4)}
    with pytest.raises(ValueError, match="paired real/complex"):
        real_translation_character_sign([Fraction(1, 2), 0, 0], [Fraction(1, 2)] * 3)


@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("centered,k,basis", [
    (False, [0, 0, Fraction(1, 4)], [1, 1, 4]),
    (True, [Fraction(1, 2), 0, 0], [2, 1, 1]),
])
def test_missing_independent_complex_phase_is_rejected(operation, centered, k, basis):
    parent = _parent(centered)
    values = np.tile([0.01, 0.0, 0.0], (len(parent), 1))
    engine = DistortionEngine()
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        if operation == "single":
            engine.generate_single_mode(parent, values, 1.0, basis, k_vector=k)
        elif operation == "mixed":
            engine.generate_modes(parent, basis, parent_displacements=values, k_vector=k)
        else:
            engine.lift_mode_displacements(parent, basis, {"mode": values}, k_vector=k)


@pytest.mark.parametrize("opd", [[1, 0], [0, 1]])
def test_opd_components_cannot_be_silently_ignored(opd):
    with pytest.raises(ValueError, match="opd_direction cannot be interpreted"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0.0, 0.0]]), 1.0, [1, 1, 4],
            k_vector=[0, 0, Fraction(1, 4)], opd_direction=opd,
        )


@pytest.mark.parametrize("origin", [(0, 0, 0), (0.23, -0.17, 0.4)])
def test_self_conjugate_real_mode_has_exact_signs_and_origin_covariance(origin):
    parent = _parent(origin=origin)
    engine = DistortionEngine()
    values = np.asarray([[0.01, 0.0, 0.0]])
    reference, lifted = engine.lift_mode_displacements(
        parent, [1, 1, 2], {"mode": values}, k_vector=[0, 0, Fraction(1, 2)],
    )
    parent_positions = reference.frac_coords @ np.diag([1, 1, 2])
    images = np.round(parent_positions - parent.frac_coords[0]).astype(int)
    expected = np.asarray([(-1) ** int(image[2]) for image in images]) * 0.01
    assert np.array_equal(lifted["mode"][:, 0], expected)
    generated = engine.generate_single_mode(
        parent, values, 1.0, [1, 1, 2], k_vector=[0, 0, Fraction(1, 2)],
    )
    observed = generated.frac_coords - reference.frac_coords
    observed -= np.round(observed)
    assert np.allclose(observed, lifted["mode"], atol=1e-14)


def test_centered_self_conjugate_real_mode_requires_internal_site_character():
    parent = _parent(centered=True)
    engine = DistortionEngine()
    k = [1, 0, 0]  # Real -1 on I centering, +1 on conventional translations.
    good = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    reference, lifted = engine.lift_mode_displacements(
        parent, [1, 1, 1], {"mode": good}, k_vector=k,
    )
    assert len(reference) == 2
    assert np.array_equal(lifted["mode"], good)
    with pytest.raises(ValueError, match=r"does not prove.*translation character"):
        engine.generate_modes(
            parent, parent_displacements=np.abs(good), k_vector=k,
        )


def test_sign_only_wavevector_requires_a_periodic_child_cell():
    with pytest.raises(ValueError, match="not periodic"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0, 0]]), k_vector=[0, 0, Fraction(1, 2)],
        )


def test_nonzero_small_wavevector_is_not_fuzzed_into_gamma():
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0, 0]]), k_vector=[0, 0, "1/10000000000"],
        )


def test_api_equal_atom_count_cannot_bypass_non_gamma_guard():
    api = object.__new__(IsoDistort)
    api.structure = _parent()
    api.mode_displacements = {"mode": {"displacements": np.asarray([[0.01, 0, 0]])}}
    api._selected_subgroup = None
    api.phase_path = SimpleNamespace(k_vector=[0, 0, Fraction(1, 4)])
    api._dist_engine = DistortionEngine()
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        api._method4_reference_modes(api.structure.copy())
