from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Lattice

from isocore.distortion.strain import decompose_homogeneous_strain
from isocore.distortion.strain_modes import engineering_voigt_to_tensor

_REPO = Path(__file__).resolve().parents[2]
_F02 = _REPO / "output_compare" / "EuAl4 Parent.cif" / "官网" / "Method4" / "F02"


def test_homogeneous_strain_identity_is_zero_and_contract_is_explicit() -> None:
    parent = np.diag([4.0, 5.0, 6.0])
    basis = np.array(
        [[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
    )

    result = decompose_homogeneous_strain(parent, basis, basis @ parent)

    assert result.applied_engineering_q_parent_basis == pytest.approx(
        np.zeros(6), abs=1e-14
    )
    assert result.coordinate_frame == "parent_lattice_basis"
    assert result.voigt_order == ("11", "22", "33", "2*23", "2*13", "2*12")
    assert result.column_lattice_action == "C_deformed = C_parent @ M"
    assert result.row_lattice_action == "P_deformed = M @ P_parent"
    assert result.child_lattice_action == "L_child = B @ M @ P_parent"
    assert result.metric_equation == "M @ G_parent @ M = H_parent_basis"
    assert result.relative_metric_residual < 1e-14


def _read_f02_official_strain() -> tuple[np.ndarray, np.ndarray]:
    """Return official raw CIF values and normalized IsoVIZ vector sum."""

    cif = (_F02 / "subgroup.cif").read_text(encoding="utf-8")
    raw = np.array(
        [
            float(
                re.search(
                    rf"^\s*{index}\s+E_{index}\s+([-+0-9.]+)\s*$",
                    cif,
                    re.MULTILINE,
                ).group(1)
            )
            for index in range(1, 7)
        ]
    )
    isoviz = (_F02 / "data.isoviz").read_text(encoding="utf-8")
    block = isoviz.split("!strainmodelist", 1)[1].split(
        "!displacivemodelist", 1
    )[0]
    applied = np.zeros(6)
    pairs = re.findall(
        r"^\s*\d+\s+([-+0-9.]+)\s+[-+0-9.]+\s+\d+\s+\S+\s*$"
        r"\n\s*((?:[-+0-9.]+\s+){5}[-+0-9.]+)\s*$",
        block,
        re.MULTILINE,
    )
    assert len(pairs) == 6
    for amplitude_text, vector_text in pairs:
        amplitude = float(amplitude_text)
        vector = np.array([float(value) for value in vector_text.split()])
        applied += amplitude * vector
    return raw, applied


def test_homogeneous_strain_matches_authoritative_eual4_f02_metric_contract() -> None:
    """F02 separates raw CIF coordinates from the applied normalized q."""

    parent = np.diag([4.402, 4.402, 11.163])
    basis = np.array(
        [[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
    )
    distorted = basis @ parent
    distorted[0] *= 1.08
    raw_cif, applied_isoviz = _read_f02_official_strain()

    result = decompose_homogeneous_strain(parent, basis, distorted)

    assert raw_cif == pytest.approx(
        [0.0, 0.12032, -0.00036, 0.09683, 0.0, 0.0], abs=5e-8
    )
    assert applied_isoviz == pytest.approx(
        [0.0, 0.0850794752, -0.00036, 0.1369379543, 0.0, 0.0],
        abs=6e-6,
    )
    assert not np.allclose(raw_cif, applied_isoviz, atol=1e-3)
    assert result.applied_engineering_q_parent_basis == pytest.approx(
        applied_isoviz, abs=6e-6
    )
    assert result.strain_tensor_parent_basis[1, 2] == pytest.approx(
        result.applied_engineering_q_parent_basis[3] / 2.0, abs=1e-14
    )
    assert result.relative_metric_residual < 1e-14


def test_oblique_nonidentity_basis_recovers_parent_basis_multiplier() -> None:
    parent = Lattice.from_parameters(4.2, 5.1, 6.4, 77.0, 84.0, 68.0).matrix
    basis = np.array(
        [[2.0, 0.0, 0.0], [0.0, 1.0, 1.0], [-1.0, 0.0, 1.0]]
    )
    q = np.array([0.02, -0.02, 0.03, 0.04, -0.012, 0.02])
    multiplier = np.eye(3) + engineering_voigt_to_tensor(q)
    angle = 0.37
    rigid_rotation = np.array(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    distorted = basis @ multiplier @ parent @ rigid_rotation

    result = decompose_homogeneous_strain(parent, basis, distorted)

    assert result.multiplier_parent_basis == pytest.approx(multiplier, abs=2e-9)
    assert result.strain_tensor_parent_basis == pytest.approx(
        multiplier - np.eye(3), abs=2e-9
    )
    assert result.applied_engineering_q_parent_basis == pytest.approx(
        q, abs=2e-9
    )
    assert result.reconstructed_child_metric == pytest.approx(
        distorted @ distorted.T, abs=2e-7
    )
    parent_metric = parent @ parent.T
    child_parent_metric = (
        np.linalg.inv(basis)
        @ (distorted @ distorted.T)
        @ np.linalg.inv(basis).T
    )
    assert (
        result.multiplier_parent_basis
        @ parent_metric
        @ result.multiplier_parent_basis
    ) == pytest.approx(child_parent_metric, abs=2e-7)
    assert result.relative_metric_residual < 2e-9


def test_deprecated_generic_aliases_share_parent_basis_objects() -> None:
    result = decompose_homogeneous_strain(
        np.eye(3), np.eye(3), 1.02 * np.eye(3)
    )

    with pytest.deprecated_call(match="multiplier_parent_basis"):
        assert result.strain_plus_identity is result.multiplier_parent_basis
    with pytest.deprecated_call(match="strain_tensor_parent_basis"):
        assert result.tensor is result.strain_tensor_parent_basis
    with pytest.deprecated_call(match="applied_engineering_q_parent_basis"):
        assert (
            result.voigt_engineering
            is result.applied_engineering_q_parent_basis
        )


def test_homogeneous_strain_rejects_singular_basis() -> None:
    with pytest.raises(ValueError, match="basis must be nonsingular"):
        decompose_homogeneous_strain(np.eye(3), np.zeros((3, 3)), np.eye(3))
