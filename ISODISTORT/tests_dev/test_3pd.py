"""Acceptance tests for the (3+d) superspace kernel and CLI."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml
from pymatgen.symmetry.groups import SpaceGroup

from isocore.superspace import (
    EPS,
    PHYSICAL_DIM,
    KsVector,
    SuperspaceGroup,
    SuperspaceOperation,
    apply_operation_to_mode,
    assert_ks_compatible,
    build_irreps,
    little_group,
    load_superspace_result,
    run_superspace_workflow,
    save_superspace_result,
    validate_nmod,
)
from isocore.superspace.group import wrap_translation
from isocore.superspace.kvector import SuperspaceLattice
from isocore.superspace.representation import generate_opds, transform_opd
from isocore.superspace.validate import require_invertible
from isocore.utils import (
    DimensionMismatchError,
    InputError,
    NumericalSingularError,
    SymmetryIncompatibleError,
    get_config,
)

ROOT = Path(__file__).resolve().parent.parent
CLI = ROOT / "main_terminal.py"


def _check_workflow(result, nmod: int) -> None:
    assert result.nmod == nmod
    assert result.group.dim == 3 + nmod
    assert result.ks.nmod == nmod
    assert result.lattice.metric.shape == (3 + nmod, 3 + nmod)
    assert result.irreps
    assert result.opds
    assert result.modes
    for mode in result.modes:
        for row in mode.project_to_3d():
            assert len(row) == 3
            assert np.all(np.isfinite(row))


def _run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(  # noqa: S603
        [sys.executable, str(CLI), *args],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        check=False,
    )


def test_eps_matches_config_lattice_tolerance():
    cfg = get_config()
    assert abs(EPS - cfg.eps) < 1e-15
    assert abs(cfg.eps - cfg.lattice_tolerance) < 1e-15


def test_nmod_d0_matches_3d_space_group():
    ssg = SuperspaceGroup.from_space_group(139, nmod=0)
    parent = SpaceGroup.from_int_number(139)
    assert ssg.dim == PHYSICAL_DIM
    assert ssg.order == len(parent.symmetry_ops)
    assert ssg.nmod == 0
    for op in ssg.operations:
        assert op.rotation.shape == (3, 3)


def test_nmod1_matches_dedicated_3p1():
    q = [0.0, 0.0, 1.0 / 6.0]
    generic = SuperspaceGroup.from_space_group(139, nmod=1, q_vectors=[q])
    dedicated = SuperspaceGroup.from_space_group_3p1(139, q_vector=q)
    assert generic.order == dedicated.order
    assert generic.dim == 4
    assert dedicated.dim == 4
    for a, b in zip(generic.operations, dedicated.operations, strict=True):
        assert a.equivalent(b, EPS)


def test_nmod2_group_dim_and_closure():
    q_vectors = [[1.0 / 6.0, 0.0, 0.0], [0.0, 1.0 / 6.0, 0.0]]
    ssg = SuperspaceGroup.from_space_group(139, nmod=2, q_vectors=q_vectors)
    assert ssg.dim == 5
    assert ssg.order > 0
    ssg.assert_closed(EPS)
    inv = ssg.find_inverse(ssg.operations[3])
    assert inv.compose(ssg.operations[3]).is_identity(EPS)


def test_nmod3_metric_and_reciprocal():
    result = run_superspace_workflow(
        139,
        3,
        q_vectors=[[1 / 6, 0, 0], [0, 1 / 6, 0], [0, 0, 1 / 6]],
        check_ks=True,
    )
    g = result.lattice.metric
    ginv = result.lattice.reciprocal_metric
    assert g.shape == (6, 6)
    assert np.allclose(g @ ginv, np.eye(6), atol=EPS)


def test_ks_reduce_project_embed():
    ks = KsVector.from_3d([0.0, 0.0, 1.0 / 6.0], nmod=1, space_group_number=139, k_point_label="LD")
    reduced = ks.reduce()
    assert reduced.external.shape == (3,)
    assert np.allclose(reduced.project_to_3d(), ks.external, atol=EPS)
    lat = run_superspace_workflow(139, 1, q_vectors=[[0, 0, 1 / 6]]).lattice
    embedded = lat.embed_3d([0.1, 0.2, 0.3], internal=[0.4])
    assert embedded.shape == (4,)
    assert np.allclose(lat.project_3d(embedded), [0.1, 0.2, 0.3], atol=EPS)


def test_ks_equivalence_and_illegal_wavevector():
    ks1 = KsVector.from_3d([0.0, 0.0, 1 / 6], 1, space_group_number=139)
    ks2 = KsVector.from_3d([0.0, 0.0, 1 / 6], 1, space_group_number=139)
    assert ks1.equivalent(ks2)
    bad = KsVector.from_3d([0.2, 0.3, 0.4], 1, space_group_number=139)
    with pytest.raises(SymmetryIncompatibleError):
        assert_ks_compatible(bad, [[0.0, 0.0, 1 / 6]])


def test_irreps_nmod1_vs_3p1():
    q = [[0.0, 0.0, 1 / 6]]
    result = run_superspace_workflow(139, 1, q_vectors=q, k_point_label="LD")
    reference = SuperspaceGroup.from_space_group_3p1(139, q_vector=q[0])
    irreps = build_irreps(result.group, result.ks, k_point_label="LD")
    reference_irreps = build_irreps(reference, result.ks, k_point_label="LD")
    assert len(irreps) == len(reference_irreps)
    assert irreps[0].dimension >= 1
    assert len(little_group(result.group, result.ks)) >= 1
    for irrep in irreps:
        assert irrep.label.startswith("LD")


def test_opd_and_modes_project_3d():
    result = run_superspace_workflow(139, 1, q_vectors=[[0, 0, 1 / 6]], k_point_label="LD")
    assert result.opds
    assert result.modes
    mode = result.modes[0]
    projected = mode.project_to_3d()
    assert projected
    assert len(projected[0]) == 3
    rotated = apply_operation_to_mode(mode, result.group.operations[0])
    assert rotated.nmod == 1
    dm = mode.to_distortion_mode(wyckoff_letters=["a"])
    assert dm.irrep_label == mode.irrep_label
    assert dm.opd_symbol == mode.opd_symbol


def test_json_yaml_roundtrip_nmod2_nmod3(tmp_path):
    r2 = run_superspace_workflow(139, 2, q_vectors=[[1 / 6, 0, 0], [0, 1 / 6, 0]], k_point_label="SM")
    r3 = run_superspace_workflow(139, 3, q_vectors=[[1 / 6, 0, 0], [0, 1 / 6, 0], [0, 0, 1 / 6]])
    p2 = tmp_path / "ss2.json"
    p3 = tmp_path / "ss3.yaml"
    save_superspace_result(r2, p2)
    save_superspace_result(r3, p3)
    b2 = load_superspace_result(p2)
    b3 = load_superspace_result(p3)
    assert r2.equivalent(b2, EPS)
    assert r3.equivalent(b3, EPS)
    raw = yaml.safe_load(p3.read_text(encoding="utf-8"))
    assert raw["nmod"] == 3
    assert raw["kind"] == "superspace"


def test_load_legacy_3p1_json(tmp_path):
    payload = {
        "kind": "superspace_3p1",
        "nmod": 1,
        "d": 1,
        "space_group_number": 139,
        "space_group_symbol": "I4/mmm",
        "schoenflies": "D4h-17",
        "group": SuperspaceGroup.from_space_group_3p1(139, q_vector=[0, 0, 1 / 6]).to_dict(),
        "ks": KsVector.from_3d([0, 0, 1 / 6], 1, space_group_number=139).to_dict(),
        "irreps": [],
        "opds": [],
        "modes": [],
    }
    path = tmp_path / "old_3p1.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    loaded = load_superspace_result(path)
    assert loaded.nmod == 1
    assert loaded.group.dim == 4
    assert loaded.group.space_group_symbol == "I4/mmm"


def test_wrap_translation_unit_interval():
    values = wrap_translation(np.array([1.0, -0.2, 2.3, 0.0]))
    assert np.all(values >= 0.0 - EPS)
    assert np.all(values < 1.0 + EPS)


def test_workflow_nmod2():
    result = run_superspace_workflow(139, 2, q_vectors=[[1 / 6, 0, 0], [0, 1 / 6, 0]], k_point_label="SM")
    _check_workflow(result, 2)
    assert result.group.space_group_symbol == "I4/mmm"
    assert result.group.schoenflies == "D4h-17"
    assert result.group.setting == "standard (IT-C)"


def test_workflow_nmod3():
    result = run_superspace_workflow(
        139, 3, q_vectors=[[1 / 6, 0, 0], [0, 1 / 6, 0], [0, 0, 1 / 6]], k_point_label="GP"
    )
    _check_workflow(result, 3)


def test_workflow_nmod1_matches_3p1():
    q = [0.0, 0.0, 1 / 6]
    result = run_superspace_workflow(139, 1, q_vectors=[q], k_point_label="LD")
    reference = SuperspaceGroup.from_space_group_3p1(139, q_vector=q)
    _check_workflow(result, 1)
    assert result.group.order == reference.order
    for actual, expected in zip(result.group.operations, reference.operations, strict=True):
        assert np.allclose(actual.rotation, expected.rotation, atol=EPS)
        assert np.allclose(actual.translation, expected.translation, atol=EPS)


@pytest.mark.parametrize("value", [-1, -3, 1.5, None, "abc", "1.25", 100])
def test_invalid_nmod_raises_input_error(value):
    with pytest.raises(InputError):
        validate_nmod(value)


def test_huge_nmod_raises():
    with pytest.raises(InputError):
        validate_nmod(99)


def test_incompatible_ks_raises():
    with pytest.raises(SymmetryIncompatibleError):
        run_superspace_workflow(139, 1, q_vectors=[[0.0, 0.0, 1 / 6]], ks_coords=[0.1, 0.2, 0.3, 1.0], check_ks=True)


def test_custom_operations_not_closed():
    ssg = SuperspaceGroup.from_space_group(2, nmod=0)
    lone = SuperspaceOperation(
        rotation=[[-1.0, 0, 0], [0, -1.0, 0], [0, 0, 1.0]],
        translation=[0.1, 0.0, 0.0],
        nmod=0,
    )
    with pytest.raises(SymmetryIncompatibleError):
        ssg.check_custom_operations_closed([lone])


def test_singular_lattice_raises():
    with pytest.raises(NumericalSingularError):
        SuperspaceLattice(lattice_3d=np.zeros((3, 3)), nmod=1)
    with pytest.raises(NumericalSingularError):
        require_invertible(np.array([[1.0, 2.0], [2.0, 4.0]]), name="metric")


def test_dimension_mismatch_on_ks_and_opd():
    with pytest.raises(DimensionMismatchError):
        KsVector(coords=[0.0, 0.0], nmod=1)
    result = run_superspace_workflow(139, 0)
    irreps = result.irreps or build_irreps(result.group, result.ks)
    opds = generate_opds(irreps, nmod=0)
    if irreps and irreps[0].matrices and opds:
        with pytest.raises(DimensionMismatchError):
            transform_opd(opds[0], irreps[0], op_index=10_000)


def test_q_vector_count_mismatch():
    with pytest.raises(DimensionMismatchError):
        SuperspaceGroup.from_space_group(139, nmod=2, q_vectors=[[0, 0, 1 / 6]])


def test_cli_superspace_d2_export(tmp_path):
    dest = tmp_path / "ss2.json"
    result = _run_cli(
        "--superspace-d", "2", "--space-group", "139", "--q-vectors", "1/6,0,0;0,1/6,0",
        "--k-label", "SM", "--export", str(dest)
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "nmod=2" in result.stdout
    assert "(3+2)" in result.stdout
    assert dest.is_file()
    assert '"nmod": 2' in dest.read_text(encoding="utf-8")


def test_cli_invalid_d_nonzero_exit():
    result = _run_cli("--superspace-d", "-1")
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "nmod" in combined.lower() or "d" in combined.lower()


def test_cli_nmod1_matches_3p1_language():
    result = _run_cli("--superspace-d", "1", "--space-group", "139", "--q-vectors", "0,0,1/6", "--k-label", "LD")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "nmod=1" in result.stdout
    assert "I4/mmm" in result.stdout
    assert "D4h-17" in result.stdout


def test_cli_load_saved_task(tmp_path):
    dest = tmp_path / "task.json"
    first = _run_cli("--superspace-d", "1", "--space-group", "139", "--q-vectors", "0,0,1/6", "--export", str(dest))
    assert first.returncode == 0, first.stdout + first.stderr
    second = _run_cli("--superspace-load", str(dest))
    assert second.returncode == 0, second.stdout + second.stderr
    assert "nmod=1" in second.stdout
    assert "I4/mmm" in second.stdout
