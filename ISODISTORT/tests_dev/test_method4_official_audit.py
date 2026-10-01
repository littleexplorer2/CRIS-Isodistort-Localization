from fractions import Fraction

import numpy as np
from pymatgen.core import Lattice, Structure

from tests_dev.manual import prepare_method4_inputs as preparation
from tests_dev.manual.audit_method4_official import (
    _amplitude_comparison_tolerance,
    _complete_gamma_site_subspaces,
    _complete_origin_phase_amplitude_subspaces,
    _complete_required_amplitude_subspaces,
    _expected_failure_matches,
    _isoviz_applied_strain,
    _isoviz_applied_strain_with_rounding,
    _local_mode_identity,
    _official_mode_identity,
    _parameter_mode_subspace_key,
    _parameter_origin_phase_turns,
)


def test_isoviz_strain_parser_applies_archived_mode_vectors() -> None:
    text = """
!strainmodelist
  1    0.06016   0.10000    1 GM1+strain_1(a)
   0.70711   0.70711   0.00000   0.00000   0.00000   0.00000
  2   -0.06016   0.10000    2 GM2+strain(a)
   0.70711  -0.70711   0.00000   0.00000   0.00000   0.00000
  3    0.09683   0.10000    4 GM5+strain(a)
   0.00000   0.00000   0.00000   1.41421   0.00000   0.00000
!displacivemodelist
"""

    np.testing.assert_allclose(
        _isoviz_applied_strain(text),
        [0.0, 0.06016 * 2 * 0.70711, 0.0, 0.09683 * 1.41421, 0.0, 0.0],
        atol=1e-12,
    )
    applied, bound = _isoviz_applied_strain_with_rounding(text)
    exact = np.array(
        [0.0, 0.06016 * np.sqrt(2), 0.0, 0.09683 * np.sqrt(2), 0.0, 0.0]
    )
    assert np.all(np.abs(applied - exact) <= bound + 1e-12)


def test_expected_failure_matching_rejects_generic_server_errors() -> None:
    assert _expected_failure_matches(
        "reject_species_or_stoichiometry",
        "Types of atoms in subgroup do not match types of atoms in parent.",
    )
    assert _expected_failure_matches(
        "reject_unmatched_atom",
        "Failed to find match. Try using a larger value for dmax.",
    )
    assert not _expected_failure_matches(
        "reject_incompatible_lattice",
        "Internal server error. Please try again later.",
    )


def test_gamma_mode_identity_uses_stable_component() -> None:
    assert _local_mode_identity("GM5+__d__P1(1)__a") == _official_mode_identity(
        "[0,0,0]GM5+(a,b)[Al1:d:dsp]E(a)"
    )


def test_complete_gamma_site_subspace_accepts_repeated_official_components() -> None:
    local_rows = [
        {
            "local_label": f"GM5-__f__C1(1)__{component}",
            "identity": list(_local_mode_identity(
                f"GM5-__f__C1(1)__{component}"
            )),
            "official_label": None,
        }
        for component in "abcd"
    ]
    official_labels = [
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(a)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(b)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B2u(a)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B2u(b)",
    ]
    official_rows = [
        {
            "label": label,
            "identity": list(_official_mode_identity(label)),
        }
        for label in official_labels
    ]

    groups = _complete_gamma_site_subspaces(local_rows, official_rows)

    assert len(groups) == 1
    key, grouped_local, grouped_official = groups[0]
    assert key == ("GM5-", "f")
    assert grouped_local == local_rows
    assert grouped_official == official_rows


def test_incomplete_gamma_site_subspace_is_not_compared_by_norm() -> None:
    local_rows = [
        {
            "local_label": f"GM5-__f__C1(1)__{component}",
            "identity": list(_local_mode_identity(
                f"GM5-__f__C1(1)__{component}"
            )),
            "official_label": None,
        }
        for component in "abcd"
    ]
    official_labels = [
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(a)",
        "[0,0,0]GM5-(a,b)[O:f:dsp]B3u(b)",
    ]
    official_rows = [
        {
            "label": label,
            "identity": list(_official_mode_identity(label)),
        }
        for label in official_labels
    ]

    assert _complete_gamma_site_subspaces(local_rows, official_rows) == []


def test_exact_amplitude_tolerance_does_not_absorb_structure_mismatch() -> None:
    tolerance, evidence = _amplitude_comparison_tolerance(
        {},
        {"origin_aligned_root_sum_squared_distance_angstrom": 0.25},
    )

    assert np.isclose(tolerance, 5.1e-6)
    assert evidence["aligned_structure_rss_bound_angstrom"] == 0.0


def test_noisy_amplitude_tolerance_uses_aligned_export_rss_bound() -> None:
    tolerance, evidence = _amplitude_comparison_tolerance(
        {"noise_sigma_fractional": 2e-5},
        {"origin_aligned_root_sum_squared_distance_angstrom": 4.2e-5},
    )

    assert np.isclose(tolerance, 4.71e-5)
    assert np.isclose(evidence["aligned_structure_rss_bound_angstrom"], 4.2e-5)
    assert "Cauchy-Schwarz" in evidence["derivation"]


def test_parameter_mode_identity_retains_wavevector_and_site_irrep() -> None:
    assert _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_2(b)"
    ) == _official_mode_identity(
        "[0,0,1/6]LD1(a,b)[Al2:e:dsp]A1_2(b)"
    )
    assert _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_2(b)"
    ) != _local_mode_identity(
        "LD1[0,0,1/3]__e__A1_2(b)"
    )
    assert _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_1(b)"
    ) != _local_mode_identity(
        "LD1[0,0,1/6]__e__A1_2(b)"
    )


def test_complete_parameter_subspace_uses_basis_independent_group() -> None:
    identity_a = _local_mode_identity("LD1[0,0,1/6]__a__A2u(a)")
    identity_b = _local_mode_identity("LD1[0,0,1/6]__a__A2u(b)")
    assert _parameter_mode_subspace_key(identity_a) == (
        "LD1[0,0,1/6]", "a", "A2u",
    )
    rows = [
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(a)",
            "identity": list(identity_a),
            "official_label": "official-a",
        },
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(b)",
            "identity": list(identity_b),
            "official_label": "official-b",
        },
    ]
    required = {row["local_label"] for row in rows}
    groups = _complete_required_amplitude_subspaces(
        rows, required, {identity_a, identity_b},
    )
    assert groups == [rows]


def test_partial_parameter_subspace_keeps_componentwise_comparison() -> None:
    identity_a = _local_mode_identity("LD1[0,0,1/6]__a__A2u(a)")
    identity_b = _local_mode_identity("LD1[0,0,1/6]__a__A2u(b)")
    rows = [
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(a)",
            "identity": list(identity_a),
            "official_label": "official-a",
        },
        {
            "local_label": "LD1[0,0,1/6]__a__A2u(b)",
            "identity": list(identity_b),
            "official_label": "official-b",
        },
    ]
    groups = _complete_required_amplitude_subspaces(
        rows,
        {"LD1[0,0,1/6]__a__A2u(a)"},
        {identity_a, identity_b},
    )
    assert groups == []


def test_parameter_origin_phase_is_exact_and_periodic() -> None:
    key = ("Y1[1/3,1/2,0]", "d", "Eu")

    assert _parameter_origin_phase_turns(key, [1.0, 0.0, 0.0]) == Fraction(1, 3)
    assert _parameter_origin_phase_turns(key, [0.0, -1.0, 0.0]) == Fraction(1, 2)
    assert _parameter_origin_phase_turns(key, [3.0, 0.0, 0.0]) == 0


def test_origin_phase_subspace_accepts_full_local_branch_in_official_expansion() -> None:
    identity_a = _local_mode_identity("Y1[1/3,1/2,0]__d__Eu(a)")
    identity_b = _official_mode_identity(
        "[1/3,1/2,0]Y1(a,b;0,0)[ND:d:dsp]Eu(b)"
    )
    local_rows = [{
        "local_label": "Y1[1/3,1/2,0]__d__Eu(a)",
        "identity": list(identity_a),
        "official_label": "official-a",
    }]
    official_rows = [
        {"label": "official-a", "identity": list(identity_a)},
        {"label": "official-b", "identity": list(identity_b)},
    ]

    groups = _complete_origin_phase_amplitude_subspaces(
        local_rows,
        {local_rows[0]["local_label"]},
        official_rows,
        [1.0, 0.0, 0.0],
    )

    assert len(groups) == 1
    key, grouped_local, grouped_official, phase_turns = groups[0]
    assert key == ("Y1[1/3,1/2,0]", "d", "Eu")
    assert grouped_local == local_rows
    assert grouped_official == official_rows
    assert phase_turns == Fraction(1, 3)


def test_origin_phase_subspace_rejects_trivial_phase_or_partial_local_branch() -> None:
    identity_a = _local_mode_identity("Y1[1/3,1/2,0]__d__Eu(a)")
    local_rows = [{
        "local_label": "Y1[1/3,1/2,0]__d__Eu(a)",
        "identity": list(identity_a),
        "official_label": "official-a",
    }]
    official_rows = [{"label": "official-a", "identity": list(identity_a)}]

    assert _complete_origin_phase_amplitude_subspaces(
        local_rows,
        {local_rows[0]["local_label"]},
        official_rows,
        [3.0, 0.0, 0.0],
    ) == []
    assert _complete_origin_phase_amplitude_subspaces(
        local_rows,
        set(),
        official_rows,
        [1.0, 0.0, 0.0],
    ) == []


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
