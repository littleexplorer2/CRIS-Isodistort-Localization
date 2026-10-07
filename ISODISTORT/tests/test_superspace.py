"""Unit tests for smodes parsing, folding k, nmod, and IsoVIZ parametric labels."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.operations import SymmOp
from test_wsl import _mode_core_token

from backend.models.iso_mode_models import ModeIdentity
from backend.wrappers.iso_wrapper import DistortionMode
from backend.wrappers.smodes_wrapper import (
    SmodesAtom,
    SmodesModeBlock,
    parse_smodes_mode_blocks,
)
from features.export.isodistort_isoviz import _isoviz_mode_label
from features.method1.search_methods import _k_compatible_with_supercell
from features.method2.superspace import (
    ParametricModeResult,
    _child_origin_candidates,
    _collect_smodes_blocks,
    _displacement_actions,
    _harmonic_order,
    _is_self_conjugate_parent_k,
    _KCandidate,
    _linearly_independent_mode_items,
    _map_smodes_block_to_supercell,
    _parent_little_group_orbit,
    _project_invariant,
    _quadrature_supercell_shift,
    _with_unresolved_identity_reason,
    folding_k_candidates,
    format_k_token,
    is_harmonic_of,
    opd_free_components,
    smodes_k_spec,
)

_SAMPLE = (
    Path(__file__).resolve().parent / "fixtures" / "smodes_sample.out"
)


def _kp(label, coords, *, special=True, parameters=None):
    return SimpleNamespace(
        label=label,
        coordinates=list(coords),
        is_special=special,
        parameters=list(parameters or []),
    )


def test_parse_smodes_sample_blocks():
    text = _SAMPLE.read_text(encoding="utf-8", errors="replace")
    blocks = parse_smodes_mode_blocks(text)
    assert blocks
    irreps = {b.irrep for b in blocks}
    assert "GM3-" in irreps
    assert "GM4-" in irreps
    gm3 = next(b for b in blocks if b.irrep == "GM3-")
    assert gm3.k_label == "GM"
    dense = gm3.dense_displacements()
    assert dense.shape[1] == 3
    assert float(np.max(np.abs(dense))) > 0
    z_blocks = [b for b in blocks if b.k_label == "Z"]
    assert z_blocks
    assert z_blocks[0].superlattice_cart.shape == (3, 3)


def test_is_harmonic_of_line_and_gamma():
    q = np.array([0.0, 0.0, 1.0 / 6.0])
    assert is_harmonic_of(np.array([0.0, 0.0, 1.0 / 6.0]), [q])
    assert is_harmonic_of(np.array([0.0, 0.0, 1.0 / 3.0]), [q])
    assert is_harmonic_of(np.array([0.0, 0.0, 0.5]), [q])
    assert is_harmonic_of(np.zeros(3), [q])  # 3D average structure
    # X of I4/mmm is not a harmonic of LD along c*
    assert not is_harmonic_of(np.array([0.5, 0.5, 0.0]), [q])
    assert _harmonic_order(np.array([0.0, 0.0, 5.0 / 6.0]), q) == 5
    assert _harmonic_order(np.array([0.0, 0.0, 1.0 / 6.0]), q) == 1


def test_harmonic_order_uses_the_exact_finite_period_not_a_fixed_cap():
    q = np.array([0.0, 0.0, 1.0 / 30.0])
    thirteenth = np.array([0.0, 0.0, 13.0 / 30.0])

    assert is_harmonic_of(thirteenth, [q])
    assert _harmonic_order(thirteenth, q) == 13
    # Reciprocal equivalence uses the primitive translations of the centered
    # parent.  In an I lattice (5/4,0,0) is the fifth, not the first, harmonic
    # of (1/4,0,0), because (1,0,0) is not an allowed reciprocal vector.
    assert _harmonic_order(
        np.array([5.0 / 4.0, 0.0, 0.0]),
        np.array([1.0 / 4.0, 0.0, 0.0]),
        parent_sg=229,
    ) == 5


def test_format_k_token_preserves_long_exact_commensurate_periods():
    assert format_k_token(1.0 / 49.0) == "1/49"
    assert format_k_token(Fraction(1, 2_000_001)) == "1/2000001"
    assert format_k_token(1e-13) == "1/10000000000000"
    assert format_k_token(-1e-13) == "-1/10000000000000"


def test_opd_free_components_c1():
    assert opd_free_components("(a,b)", "C1") == ["a", "b"]
    assert opd_free_components("(a)", "P1") == ["a"]
    assert opd_free_components("", "C1") == ["a", "b"]


def test_parametric_rejection_reason_is_bound_to_every_unresolved_identity():
    mode = DistortionMode(
        irrep_label="Y1",
        amplitude_key="Y1__f__unresolved(a)",
        wyckoff_site="f",
        wyckoff_orbit_id="sg123:2f:O:test",
        mode_identity=ModeIdentity.unresolved(
            parent_sg=123,
            global_irrep="Y1",
            k_coordinates=("1/3", "1/2", "0"),
            wyckoff_letter="f",
            orbit_id="sg123:2f:O:test",
            reason="smodes_basis_has_no_iso_microscopic_identity",
        ),
    )
    result = ParametricModeResult(
        modes=[mode],
        supercell_displacements={mode.amplitude_key: np.ones((1, 3))},
        labels={},
        nmod=0,
    )
    reason = (
        "parameter_k_microscopic_canonicalization_rejected; "
        "criterion=ValueError: candidate_columns_not_independent; "
        "parent_SG=123; k=(1/3,1/2,0); KVALUE=(2/3); "
        "IR=Y1; direction_selector=(a,0;0,0)"
    )

    updated = _with_unresolved_identity_reason(result, reason)

    assert updated.note == reason
    assert updated.modes[0].mode_identity.reason == reason
    assert updated.modes[0].mode_identity.status == "unresolved"


def test_cached_displacement_actions_equal_reynolds_projection():
    structure = Structure(
        Lattice.cubic(4.0),
        ["Fe", "Fe"],
        [[0.25, 0, 0], [0.75, 0, 0]],
    )
    identity = SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3))
    inversion = SymmOp.from_rotation_and_translation(-np.eye(3), np.zeros(3))
    ops = [identity, inversion]
    disp = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    actions = _displacement_actions(structure, ops)
    assert actions is not None
    projected = _project_invariant(disp, structure, ops, actions=actions)
    expected_inverted = np.asarray([[-4.0, -5.0, -6.0], [-1.0, -2.0, -3.0]])
    assert np.allclose(projected, 0.5 * (disp + expected_inverted))


def test_displacement_action_site_tolerance_uses_cartesian_angstrom_metric():
    """A fractional residual is judged in the actual direct-cell metric."""
    structure = Structure(
        Lattice.orthorhombic(10.0, 2.0, 3.0), ["Fe"], [[0.0, 0.0, 0.0]],
    )
    within_one_millangstrom = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.00005, 0.0, 0.0],
    )
    outside_one_millangstrom = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.00015, 0.0, 0.0],
    )

    assert _displacement_actions(
        structure, [within_one_millangstrom], atol=0.001,
    ) is not None
    assert _displacement_actions(
        structure, [outside_one_millangstrom], atol=0.001,
    ) is None


def test_displacement_action_uses_true_minimum_image_in_skew_cell():
    """Fractional component wrapping is not a nearest-image algorithm."""
    structure = Structure(
        Lattice.from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 10.0),
        ["Fe"],
        [[0.0, 0.0, 0.0]],
    )
    operation = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.49, 0.49, 0.0],
    )

    # The component-wrapped vector is 9.76 Å, while the true periodic image
    # is 0.894 Å because the a and b axes are nearly parallel.
    assert _displacement_actions(structure, [operation], atol=1.0) is not None


def test_displacement_action_solves_global_species_bijection():
    """A perfect matching must not fail because of greedy target ownership."""
    structure = Structure(
        Lattice.cubic(10.0),
        ["Fe", "Fe"],
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]],
    )
    operation = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.11, 0.0, 0.0],
    )

    actions = _displacement_actions(structure, [operation], atol=1.2)
    assert actions is not None
    targets, _rotation = actions[0]
    assert targets.tolist() == [0, 1]


def test_child_origin_candidates_preserve_supercell_integer_lifts():
    basis = np.diag([1.0, 1.0, 6.0])
    translations = np.asarray([[0.0, 0.0, z] for z in range(6)])
    candidates = _child_origin_candidates(
        [0.75, 0.75, 0.25], basis, translations,
    )
    # The website's equivalent parent origin (-1/4,-1/4,9/4) maps to
    # child z=3/8.  Reducing 9/4 modulo the parent cell must not lose it.
    assert any(np.allclose(item, [0.75, 0.75, 0.375]) for item in candidates)
    assert len(candidates) == 6


def test_mode_rank_keeps_near_parallel_independent_and_drops_dependency():
    mode = SimpleNamespace()
    first = np.asarray([[1.0, 0.0, 0.0]])
    near_parallel = np.asarray([[1.0, 1e-3, 0.0]])
    dependent = first + near_parallel
    items = [
        ("a", mode, first, "a"),
        ("b", mode, near_parallel, "b"),
        ("c", mode, dependent, "c"),
    ]
    kept = _linearly_independent_mode_items(items, np.eye(3))
    assert [item[0] for item in kept] == ["a", "b"]


def test_parent_little_group_orbit_completes_vector_irrep():
    structure = Structure(Lattice.tetragonal(4.0, 5.0), ["Fe"], [[0, 0, 0]])
    seed = np.asarray([[1.0, 0.0, 0.0]])
    orbit = _parent_little_group_orbit(
        seed,
        structure,
        np.eye(3),
        123,
        np.zeros(3),
        2,
    )
    assert len(orbit) == 2
    cart = [item @ structure.lattice.matrix for item in orbit]
    assert abs(float(np.vdot(cart[0].ravel(), cart[1].ravel()))) < 1e-8


def test_self_conjugate_k_respects_parent_centering():
    assert _is_self_conjugate_parent_k(np.array([0.5, 0.5, 0.0]), 123)
    assert _is_self_conjugate_parent_k(np.array([0.5, 0.0, 0.5]), 139)
    assert not _is_self_conjugate_parent_k(np.array([0.0, 0.0, 1 / 6]), 139)
    # (0,0,1) is not a reciprocal vector of an I-centred conventional cell.
    assert not _is_self_conjugate_parent_k(np.array([0.0, 0.0, 0.5]), 139)


def test_folding_k_candidates_i4mmm_1x1x6():
    basis = np.diag([1.0, 1.0, 6.0])
    kpoints = [
        _kp("GM", ["0", "0", "0"], special=True),
        _kp("M", ["1", "1", "1"], special=True),
        _kp("X", ["1/2", "1/2", "0"], special=True),
        _kp("LD", ["0", "0", "g"], special=False, parameters=["g"]),
    ]
    primary = SimpleNamespace(
        label="LD",
        coords=np.array([0.0, 0.0, 1.0 / 6.0]),
        params=["1/6"],
        coord_tokens=("0", "0", "1/6"),
    )
    prim = _KCandidate("LD", primary.coords, ["1/6"], ("0", "0", "1/6"))
    lockin = folding_k_candidates(139, kpoints, basis, prim, nmod=0)
    labels = {(c.label, tuple(format_k_token(x) for x in c.coords)) for c in lockin}
    assert ("GM", ("0", "0", "0")) in labels
    assert ("M", ("1", "1", "1")) in labels
    assert ("LD", ("0", "0", "1/6")) in labels
    assert ("LD", ("0", "0", "1/3")) in labels
    assert ("LD", ("0", "0", "5/6")) in labels
    assert ("LD", ("0", "0", "1/2")) in labels
    assert ("X", ("1/2", "1/2", "0")) not in labels

    harm = folding_k_candidates(139, kpoints, basis, prim, nmod=1)
    harm_labels = {(c.label, tuple(format_k_token(x) for x in c.coords)) for c in harm}
    assert ("LD", ("0", "0", "1/6")) in harm_labels
    assert ("LD", ("0", "0", "1/3")) in harm_labels
    assert ("GM", ("0", "0", "0")) in harm_labels
    assert ("X", ("1/2", "1/2", "0")) not in harm_labels


def test_folding_k_candidates_uses_parent_primitive_lattice_for_fractional_basis():
    """An I-parent primitive basis is nonsingular although det(B_conv)=1/2."""
    half = 0.5
    basis = np.array([
        [-half, half, half],
        [half, -half, half],
        [half, half, -half],
    ])
    kpoints = [
        _kp("GM", ["0", "0", "0"], special=True),
        _kp("M", ["1", "1", "1"], special=True),
    ]
    folded = folding_k_candidates(139, kpoints, basis, None, nmod=0)
    assert [(item.label, item.coord_tokens) for item in folded] == [
        ("GM", ("0", "0", "0")),
    ]


def test_folding_candidates_canonicalize_i_boundary_star_before_early_stop():
    """Equivalent P-star boundary arms must not displace folded GM/M classes."""
    basis = np.array([
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [0.0, 0.0, 2.0],
    ])
    kpoints = [
        _kp("GM", ["0", "0", "0"], special=True),
        _kp("M", ["1", "1", "1"], special=True),
        _kp("P", ["1/2", "1/2", "1/2"], special=True),
    ]
    primary = _KCandidate(
        "P", np.array([0.5, 0.5, 0.5]), None, ("1/2", "1/2", "1/2"),
    )

    folded = folding_k_candidates(
        139, kpoints, basis, primary, nmod=0, child_centering="I",
    )

    assert [item.label for item in folded] == ["P", "GM", "M"]
    assert len(folded) == 3


def test_k_compatibility_uses_row_basis_for_axis_permutation():
    """A row-basis cell must test each child translation against k.

    NdNiO2 Method 2 uses this representative for Y=(1/3,1/2,0).  The
    old transpose convention rejected Y even though every child translation
    has an integral Bloch phase.
    """
    basis = np.array([
        [-3.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 2.0, 0.0],
    ])
    y = np.array([1.0 / 3.0, 1.0 / 2.0, 0.0])
    assert np.allclose(basis @ y, [-1.0, 0.0, 1.0])
    assert _k_compatible_with_supercell(y, basis)

    incompatible = np.array([1.0 / 4.0, 1.0 / 2.0, 0.0])
    assert not _k_compatible_with_supercell(incompatible, basis)


def test_folding_candidates_keep_compatible_rotated_star_arm():
    """A displayed line representative may differ from the active star arm."""
    basis = np.array([
        [-3.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 2.0, 0.0],
    ])
    kpoints = [
        _kp("DT", ["0", "b", "0"], special=False, parameters=["b"]),
    ]
    primary = _KCandidate(
        "Y", np.array([1.0 / 3.0, 1.0 / 2.0, 0.0]),
        ["1/3"], ("1/3", "1/2", "0"),
    )
    candidates = folding_k_candidates(123, kpoints, basis, primary, nmod=0)
    dt_thirds = [
        item for item in candidates
        if item.label == "DT" and item.coord_tokens == ("0", "1/3", "0")
    ]
    assert dt_thirds
    assert any(
        _k_compatible_with_supercell(item.active_coords, basis)
        for item in dt_thirds
    )


def test_rotated_star_arm_keeps_physical_k_for_phase_generation():
    """The website representative must not replace a rotated arm's phase k."""
    parent = Structure(Lattice.cubic(4.0), ["O"], [[0, 0, 0]])
    raw = SmodesModeBlock(
        irrep="DT1",
        degeneracy=1,
        k_label="DT",
        k_coords=("0", "1/3", "0"),
        superlattice_cart=np.eye(3) * 4.0,
        atoms=[SmodesAtom(1, "O", np.zeros(3))],
        atom_indices=[1],
        displacements_cart=np.array([[1.0, 0.0, 0.0]]),
    )

    class _FakeSmodes:
        @staticmethod
        def compute_mode_blocks(*_args, **_kwargs):
            return [raw]

    candidate = _KCandidate(
        "DT",
        np.array([0.0, 1.0 / 3.0, 0.0]),
        ["1/3"],
        ("0", "1/3", "0"),
        active_coords=np.array([1.0 / 3.0, 0.0, 0.0]),
        arm_rotation=np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
    )
    blocks = _collect_smodes_blocks(
        _FakeSmodes(), parent, 123, [], [(candidate, ("DT", ["1/3"]))],
    )
    assert len(blocks) == 1
    assert blocks[0].k_coords == ("0", "1/3", "0")
    assert blocks[0].active_k_coords == ("1/3", "0", "0")


def test_smodes_mapping_rejects_an_uncovered_child_atom():
    """An unmatched row must not be emitted as a genuine zero displacement."""
    child = Structure(
        Lattice.cubic(10.0),
        ["Fe", "Fe"],
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )
    incomplete = SmodesModeBlock(
        irrep="X1",
        degeneracy=1,
        k_label="X",
        k_coords=("1/2", "0", "0"),
        superlattice_cart=np.eye(3) * 10.0,
        atoms=[SmodesAtom(1, "Fe", np.zeros(3))],
        atom_indices=[1],
        displacements_cart=np.array([[1.0, 0.0, 0.0]]),
    )

    assert _map_smodes_block_to_supercell(incomplete, child) is None


def test_smodes_mapping_preserves_a_covered_zero_displacement_atom():
    """Coverage is positional; a symmetry-required zero vector remains valid."""
    child = Structure(
        Lattice.cubic(10.0),
        ["Fe", "Fe"],
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )
    complete = SmodesModeBlock(
        irrep="X1",
        degeneracy=1,
        k_label="X",
        k_coords=("1/2", "0", "0"),
        superlattice_cart=np.eye(3) * 10.0,
        atoms=[
            SmodesAtom(1, "Fe", np.zeros(3)),
            SmodesAtom(2, "Fe", np.array([5.0, 0.0, 0.0])),
        ],
        atom_indices=[1],
        displacements_cart=np.array([[1.0, 0.0, 0.0]]),
    )

    mapped = _map_smodes_block_to_supercell(complete, child)
    assert mapped is not None
    assert np.allclose(mapped[0], [0.1, 0.0, 0.0])
    assert np.allclose(mapped[1], 0.0)


def test_smodes_k_spec_converts_official_ld_g():
    iso_kp = _kp("LD", ["0", "0", "2a"], special=False, parameters=["2a"])
    cand = _KCandidate(
        "LD", np.array([0.0, 0.0, 1.0 / 6.0]), ["1/6"], ("0", "0", "1/6"),
    )
    label, params = smodes_k_spec(139, cand, [iso_kp])
    assert label == "LD"
    assert params == ["1/12"]


def test_isoviz_mode_label_keeps_parametric_k_prefix():
    pretty = "I4/mmm[0,0,1/6]LD1(a,b)[Eu0:a:dsp]A2u(a)"
    assert _isoviz_mode_label(pretty, "LD1") == "[0,0,1/6]LD1[Eu0:a:dsp]A2u(a)"
    gamma = "I4/mmm[0,0,0]GM1+(a)[Al2:e:dsp]A1(a)"
    assert _isoviz_mode_label(gamma, "GM1+") == "GM1+[Al2:e:dsp]A1(a)"
    mpoint = "I4/mmm[1,1,1]M3-(a)[Eu0:a:dsp]A2u(a)"
    assert _isoviz_mode_label(mpoint, "M3-") == "[1,1,1]M3-[Eu0:a:dsp]A2u(a)"


def test_quadrature_shift_covers_long_primitive_supercell():
    """A fixed 0..5 parent-cell window misses k·T≈1/4 on a long P lattice."""
    basis = np.diag([1.0, 1.0, 24.0])
    k = np.array([0.0, 0.0, 1.0 / 24.0])
    shift = _quadrature_supercell_shift(k, basis, parent_sg=1)
    assert shift is not None
    phase = float(np.dot(k, np.asarray(shift, dtype=float) @ basis))
    phase -= round(phase)
    assert abs(abs(phase) - 0.25) < 0.02


def test_quadrature_shift_solves_long_period_exactly_without_a_span_cap():
    basis = np.diag([1.0, 1.0, 400.0])
    k = np.array([0.0, 0.0, 1.0 / 400.0])

    shift = _quadrature_supercell_shift(k, basis, parent_sg=1)
    assert shift is not None
    parent_translation = np.asarray(shift) @ basis
    phase = float(np.dot(k, parent_translation)) % 1.0
    assert phase == 0.25


def test_quadrature_shift_fails_closed_when_the_character_has_only_sign_phase():
    basis = np.diag([1.0, 1.0, 2.0])
    k = np.array([0.0, 0.0, 1.0 / 2.0])

    assert _quadrature_supercell_shift(k, basis, parent_sg=1) is None
    # I centering adds quarter-phase primitive translations for the same
    # conventional k point, so the centered-parent case remains resolvable.
    assert _quadrature_supercell_shift(k, basis, parent_sg=139) is not None


def test_gold_compact_m_labels_are_not_gamma():
    """Official IsoVIZ omits ``[1,1,1]`` on I4/mmm M compact labels."""
    compact = "M1+[Al2:d:dsp]B2(a)"
    full = "I4/mmm[1,1,1]M1+(a)[Al2:d:dsp]B2(a)"
    assert _mode_core_token(compact) == "1,1,1|M1+|d|B2(a)"
    assert _mode_core_token(full) == _mode_core_token(compact)
    assert _mode_core_token("M3-[Eu0:a:dsp]A2u(a)") == "1,1,1|M3-|a|A2u(a)"
    assert _mode_core_token("GM1+[Al2:e:dsp]A1(a)") == "0,0,0|GM1+|e|A1(a)"
