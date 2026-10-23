"""Tests for application.session; migrated without changing scientific expectations."""
from __future__ import annotations

from application.session import IsoDistort as _basics_IsoDistort
from application.session import IsoDistort as _distortion_IsoDistort
from application.session import IsoDistort as _parent_setting_guard_IsoDistort
from application.session import IsoDistort as _partial_microscopic_mapping_IsoDistort
from application.session import IsoDistort as _production_displacive_core_IsoDistort
from application.session import IsoDistort as _real_mode_fallback_IsoDistort
from application.session import IsoDistort as _robustness_IsoDistort
from application.session import IsoDistort as _search_methods_IsoDistort
from application.session import IsoDistort as _session_result_contract_IsoDistort
from application.session import IsoDistort as _structure_IsoDistort
from application.session import IsoDistort as _theory_contracts_IsoDistort
from application.session import IsoDistort as _wsl_IsoDistort
from application.session import UnresolvedModeIdentityError as _production_displacive_core_UnresolvedModeIdentityError
from application.session import core_api
from application.session import core_api as core_api_module
from application.session.core_api import ExportBatchPlan
from application.session.core_api import ExportCandidateFailure
from application.session.core_api import ExportPreparationError
from application.session.core_api import IsoDistort as _export_fail_closed_IsoDistort
from application.session.core_api import IsoDistort as _formats_IsoDistort
from application.session.core_api import ModeIdentityFailure
from application.session.core_api import UnresolvedModeIdentityError as _export_fail_closed_UnresolvedModeIdentityError
from application.session.publication import publish_export_batch
from backend.export.distortion_formats import FORMAT_CIF
from backend.export.distortion_formats import SubgroupExportSpec
from backend.export.distortion_formats import cart_normalized_mode_matrix
from backend.export.distortion_formats import method3_case_folder
from backend.export.distortion_formats import render_cif
from backend.export.distortion_formats import render_complete_modes
from backend.export.distortion_formats import render_isoviz
from backend.export.distortion_formats import render_topas
from backend.input_cif import SymmetryValidator
from backend.input_cif import read_cif
from backend.method1 import Method1ResultItem
from backend.method1.search import _basis_is_same_lattice
from backend.method2 import DistortionMapper as _production_displacive_core_DistortionMapper
from backend.method2 import Method2Query
from backend.method2 import Method2Solver
from backend.method2 import superspace
from backend.method2 import superspace as superspace_module
from backend.method2.distortion_mapper import DistortionMapper as _partial_microscopic_mapping_DistortionMapper
from backend.method3 import Method3ResultItem
from backend.method4 import DistortionEngine
from backend.method4.strain_modes import compute_homogeneous_strain_modes
from backend.models.iso_mode_models import MacroscopicTensorBasis
from backend.models.iso_mode_models import MicroscopicColumnProvenance
from backend.models.iso_mode_models import MicroscopicVectorBlock
from backend.models.iso_mode_models import MicroscopicVectorRow
from backend.models.iso_mode_models import ModeIdentity
from backend.utils import IsodistortError
from backend.utils import get_config
from backend.utils.affine_embeddings import embedding_from_identity
from backend.utils.affine_embeddings import parent_affine_group
from backend.utils.lattice import is_sublattice
from backend.utils.lattice import rational_matrix
from backend.wrappers import BushMode
from backend.wrappers import BushMode as _search_methods_BushMode
from backend.wrappers import DistortionMode
from backend.wrappers import DistortionMode as _search_methods_DistortionMode
from backend.wrappers import SubgroupInfo
from backend.wrappers import SubgroupInfo as _formats_SubgroupInfo
from backend.wrappers import SubgroupInfo as _search_methods_SubgroupInfo
from copy import deepcopy
from data_dir import experiment_data_dir
from dataclasses import FrozenInstanceError
from dataclasses import asdict
from dataclasses import replace
from fractions import Fraction
from html.parser import HTMLParser
from pathlib import Path
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from tests.support_application_session import _mode_core_token
from tests.support_application_session import _robustness_CIFS_DIR
from tests.support_application_session import _wsl_required
from tests.support_application_session import _wsl_wsl_available
from tests.support_backend_export import _TetragonalMacroIso
from tests.support_backend_export import _bare_api
from tests.support_backend_export import _cubic
from tests.support_backend_export import _minimal_export_api
from tests.support_backend_export import _sg
from tests.support_backend_export import _strain_export_subgroup
from tests.support_backend_export import _tetragonal_macro_basis
from tests.support_backend_method2 import _active_arm_mapping_fixture
from tests.support_backend_method3 import _i4mmm_parent
from tests.support_backend_method3 import _method3_route
from tests.support_backend_method4 import _parent
from tests.support_paths import TESTS_ROOT
from types import MappingProxyType, SimpleNamespace
from unittest.mock import Mock
import application.session.publication as publication_module
import contextlib
import copy
import hashlib
import io
import json
import numpy as np
import numpy.linalg as la
import os
import pytest
import re
import shutil
import subprocess
import sys
import threading
import time
import zipfile


# Original coverage: test_basics.py.

@pytest.mark.parametrize("generate_if_missing", [False, True])
def test_method2_all_ir_inventory_rejects_any_partial_backend_failure(
    generate_if_missing: bool,
) -> None:
    """A successful IR must not hide a failed IR in the advertised full table."""

    old_candidate = SimpleNamespace(index=7)
    first_candidate = SimpleNamespace(index=0)

    class PartialBackend:
        @staticmethod
        def list_irreps(*_args):
            return [SimpleNamespace(label="A1"), SimpleNamespace(label="A2")]

        @staticmethod
        def list_subgroups(*_args, **_kwargs):
            if _args[2] == "A1":
                return [first_candidate]
            raise IsodistortError("second IR parse failed")

    api = object.__new__(_basics_IsoDistort)
    api.structure = object()
    api.symmetry_info = {"space_group_number": 123}
    api.subgroups = [old_candidate]
    api._iso = PartialBackend()
    api._resolve_iso_kparams = lambda *_args: []
    api._filter_subgroups_for_search = lambda rows, *_args: rows
    api._tag_official_kparams = lambda *_args: None

    with pytest.raises(IsodistortError, match="second IR parse failed"):
        api.list_subgroups_at_kpoint(
            "X",
            generate_if_missing=generate_if_missing,
        )

    assert api.subgroups == [old_candidate]

def test_method2_all_ir_inventory_allows_a_genuine_empty_ir() -> None:
    candidate = SimpleNamespace(index=99)

    class CompleteBackend:
        @staticmethod
        def list_irreps(*_args):
            return [SimpleNamespace(label="A1"), SimpleNamespace(label="A2")]

        @staticmethod
        def list_subgroups(*_args, **_kwargs):
            return [] if _args[2] == "A1" else [candidate]

    api = object.__new__(_basics_IsoDistort)
    api.structure = object()
    api.symmetry_info = {"space_group_number": 123}
    api.subgroups = []
    api._iso = CompleteBackend()
    api._resolve_iso_kparams = lambda *_args: []
    api._filter_subgroups_for_search = lambda rows, *_args: rows
    api._tag_official_kparams = lambda *_args: None

    result = api.list_subgroups_at_kpoint("X")

    assert result == [candidate]
    assert candidate.index == 0
    assert api.subgroups == [candidate]

# Original coverage: test_distortion.py.

_distortion_CIFS_DIR = TESTS_ROOT / "cifs_30"

_distortion_DATA_DIR = experiment_data_dir()

@pytest.mark.parametrize(
    "values",
    [
        np.zeros((1, 3)),
        np.zeros((3, 3)),
        np.zeros((2, 2)),
        np.array([[0.0, 0.0, np.nan], [0.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, np.inf], [0.0, 0.0, 0.0]]),
    ],
)
def test_supercell_modes_fail_fast_instead_of_padding_or_truncating(values) -> None:
    with pytest.raises(ValueError, match="Mode 'bad'"):
        _distortion_IsoDistort._validated_supercell_displacement("bad", values, 2)

def test_supercell_mode_validator_accepts_exact_finite_shape() -> None:
    values = np.arange(6, dtype=float).reshape((2, 3))
    checked = _distortion_IsoDistort._validated_supercell_displacement("ok", values, 2)
    assert np.array_equal(checked, values)

# --- from test_multi_mode.py ---

def _load_iso(path: Path) -> _distortion_IsoDistort:
    iso = _distortion_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(path)
    return iso

def _subgroup_with_modes(iso: _distortion_IsoDistort, min_modes: int = 2) -> tuple:
    """找首个可产生 >= min_modes 个位移模式的子群。"""
    candidates = iso.search_method_1(
        distortion_types=["displacive", "strain"],
    )
    for cand in candidates:
        try:
            iso.search_method_2(subgroup_idx=cand.subgroup.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
            continue
        if len(iso.mode_displacements) >= min_modes:
            return cand.subgroup, list(iso.mode_displacements.keys())
    raise AssertionError("未找到含 >= 2 个位移模式的子群")

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
def test_mixed_mode_linear_superposition(tmp_path):
    """混合畸变 = 各单模式位移的线性叠加（StructureMatcher 原子对应后比对）。

    位移向量以母相原胞为单位定义（mode_displacements）；混合生成的超胞
    原子位移（相对母相，周期约化）应等于各模式位移的相位调制线性叠加。
    子群基矢可能是旋转幺模胞，原子顺序会重排，故用 StructureMatcher
    建立“混合畸变原子 -> 母相原子”的对应后再比对位移。
    """
    cif = _distortion_DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = _distortion_CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")
    iso = _load_iso(cif)
    subgroup, labels = _subgroup_with_modes(iso, min_modes=2)
    label_a, label_b = labels[:2]

    base_disp = {
        lab: np.asarray(iso.mode_displacements[lab]["displacements"], dtype=float)
        for lab in (label_a, label_b)
    }
    # 若子群超胞因子 > 1，则改用超胞因子 1 的 Γ 点子群验证线性叠加
    if abs(round(float(la.det(subgroup.basis_vectors)), 3)) != 1.0:
        iso2 = _load_iso(cif)
        subs = iso2.list_subgroups_at("GM", "GM1+")
        for sg in subs:
            try:
                iso2.search_method_2(subgroup_idx=sg.index,
                                     distortion_type=["displacive"])
            except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
                continue
            if len(iso2.mode_displacements) >= 2:
                labels2 = list(iso2.mode_displacements)[:2]
                base_disp = {
                    lab: np.asarray(iso2.mode_displacements[lab]["displacements"],
                                    dtype=float)
                    for lab in labels2
                }
                amp_a, amp_b = 0.05, -0.03
                mixed = iso2.generate_mixed_distortion(
                    {labels2[0]: amp_a, labels2[1]: amp_b})
                _assert_linear_superposition(iso2, mixed, base_disp,
                                             labels2[0], labels2[1],
                                             amp_a, amp_b)
                return
        pytest.skip("无超胞因子 1 的多模式子群，线性叠加改由对称性校验覆盖")

    amp_a, amp_b = 0.05, -0.03
    mixed = iso.generate_mixed_distortion(
        {label_a: amp_a, label_b: amp_b}
    )
    _assert_linear_superposition(iso, mixed, base_disp,
                                 label_a, label_b, amp_a, amp_b)

def _assert_linear_superposition(iso, mixed, base_disp, label_a, label_b,
                                 amp_a, amp_b) -> None:
    """按“同物种 + 周期最小镜像”原子对应比对混合畸变位移与线性叠加。

    混合结构以子群超胞坐标表达（超胞 = basis @ 母相格子）；位移以母相
    分数坐标定义，故先把超胞坐标变换回母相坐标（parent = sc @ basis）
    再与期望线性叠加比对。
    """
    expected = amp_a * base_disp[label_a] + amp_b * base_disp[label_b]
    n_parent = len(iso.structure)
    assert len(mixed) == n_parent, \
        f"该子群超胞因子应为 1（{n_parent} -> {len(mixed)}）"

    # 子群基矢（母相格单位）；无超胞时为单位阵
    basis = (np.asarray(iso.phase_path.supercell_basis(), dtype=float)
             if iso.phase_path is not None else np.eye(3))
    if basis.shape != (3, 3):
        basis = np.eye(3)

    pc = np.asarray(iso.structure.frac_coords, dtype=float)
    mc_parent = np.asarray(mixed.frac_coords, dtype=float) @ basis
    max_diff = 0.0
    for i in range(n_parent):
        # 同物种中取最小镜像距离对应的混合原子
        best_j, best_dist = -1, float("inf")
        for j in range(n_parent):
            if mixed[j].species_string != iso.structure[i].species_string:
                continue
            delta = mc_parent[j] - pc[i]
            delta -= np.round(delta)
            dist = float(np.linalg.norm(delta))
            if dist < best_dist:
                best_dist, best_j = dist, j
        delta = mc_parent[best_j] - pc[i]
        delta -= np.round(delta)
        diff = delta - expected[i]
        diff -= np.round(diff)
        max_diff = max(max_diff, float(np.max(np.abs(diff))))
    assert max_diff < 1e-6, \
        f"混合位移与线性叠加偏差 {max_diff:.2e}"

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
def test_mixed_mode_symmetry_conserved(tmp_path):
    """多模式混合畸变：spglib 对称性 == 目标子群（叠加不破坏对称性）。"""
    cif = _distortion_DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = _distortion_CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")
    iso = _load_iso(cif)
    subgroup, labels = _subgroup_with_modes(iso, min_modes=2)

    mixed = iso.generate_mixed_distortion({
        labels[0]: 0.05,
        labels[1]: -0.03,
    })
    sg = SpacegroupAnalyzer(mixed, symprec=1e-3).get_space_group_number()
    assert sg == subgroup.space_group_number, \
        f"混合畸变对称性 #{sg}，应为 #{subgroup.space_group_number}"
    # 原子数 = 母相 × 超胞倍数（基矢行列式）
    det = abs(round(float(la.det(subgroup.basis_vectors)), 3))
    if det >= 1:
        assert len(mixed) == round(det * len(iso.structure)), \
            f"混合畸变原子数 {len(mixed)}，期望 {det * len(iso.structure)}"

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
def test_mixed_mode_invalid_label_raises(tmp_path):
    """混合生成传入未知模式标号：明确报错，不静默。"""
    cif = _distortion_CIFS_DIR / "sg001.cif"
    if not cif.exists():
        pytest.skip("测试 CIF 不存在")
    iso = _load_iso(cif)
    subs = iso.list_subgroups_at("GM", "GM1")
    for sg in subs:
        try:
            iso.search_method_2(subgroup_idx=sg.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
            continue
        if iso.mode_displacements:
            break
    else:
        pytest.skip("该母相无可计算位移模式")
    with pytest.raises(ValueError):
        iso.generate_mixed_distortion({"NOT_A_MODE": 0.1})

# Original coverage: test_export_fail_closed.py.

def _export_fail_closed_candidate(index: int, irrep: str) -> SubgroupInfo:
    return SubgroupInfo(
        index=index,
        space_group_number=12,
        subgroup_index=2,
        space_group_symbol="C2/m",
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label=irrep,
        opd_symbol="P1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        parent_sg=139,
    )

def _bare_export_session() -> _export_fail_closed_IsoDistort:
    iso = object.__new__(_export_fail_closed_IsoDistort)
    iso._selected_subgroup = None
    iso.phase_path = None
    iso.distortion_modes = []
    iso.mode_displacements = {}
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.distorted_structure = None
    iso.number_of_independent_modulations = 0
    iso.structure = None
    iso.distortion_types = ["strain", "displacive"]
    iso.distortion_scope = {}
    iso._mode_cache_key = None
    iso._generated_structure_cache_key = None
    return iso

def test_mode_identity_qualification_survives_declared_candidate_copy_and_serialization():
    from backend.models.iso_mode_models import ModeIdentityFailure as DomainFailure

    assert ModeIdentityFailure is DomainFailure
    failure = DomainFailure("mode-1", "orbit-a", "GM1+", ("0", "0", "0"), "unresolved source")
    candidate = _export_fail_closed_candidate(0, "GM1+")
    candidate.mode_identity_status = "unresolved"
    candidate.mode_identity_failures = (failure,)

    copied = replace(candidate, index=1)
    serialized = asdict(copied)
    assert copied.mode_identity_failures == (failure,)
    assert serialized["mode_identity_status"] == "unresolved"
    assert serialized["mode_identity_failures"][0]["reason"] == "unresolved source"

def test_export_snapshot_owns_independent_structure_candidates_and_backends(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class RuntimeStub:
        def clone_runtime_from(self, source):
            pass

        def close(self):
            pass

    monkeypatch.setattr(core_api_module, "IsoWrapper", RuntimeStub)
    monkeypatch.setattr(core_api_module, "SmodesWrapper", RuntimeStub)
    iso = _export_fail_closed_IsoDistort()
    iso.structure = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    iso.symmetry_info = {
        "space_group_number": 221,
        "space_group_symbol": "Pm-3m",
        "wyckoff_sites": [],
    }
    iso.subgroups = [_export_fail_closed_candidate(0, "GM1+")]
    iso.distortion_types = ["displacive"]
    iso.distortion_scope = {"displacive": ["Fe"]}

    worker = iso.snapshot_for_export()
    worker.structure.translate_sites([0], [0.25, 0.0, 0.0], frac_coords=True)
    worker.subgroups[0].index = 99
    worker.distortion_scope["displacive"].append("O")

    assert worker is not iso
    assert worker._iso is not iso._iso
    assert np.allclose(iso.structure[0].frac_coords, [0.0, 0.0, 0.0])
    assert iso.subgroups[0].index == 0
    assert iso.distortion_scope == {"displacive": ["Fe"]}

@pytest.mark.parametrize("cleanup_fails", [False, True])
def test_failed_runtime_clone_closes_only_new_workers_and_preserves_original_error(
    monkeypatch, cleanup_fails,
):
    created = []

    class RuntimeStub:
        def __init__(self):
            self.index = len(created)
            self.closed = False
            self.generated = {"table": [1, 2, 3]}
            created.append(self)

        def clone_runtime_from(self, source):
            self.generated = copy.deepcopy(source.generated)
            raise ValueError("runtime clone failed")

        def close(self):
            self.closed = True
            if cleanup_fails and self.index == 2:
                raise OSError("new worker cleanup failed")

    monkeypatch.setattr(core_api_module, "IsoWrapper", RuntimeStub)
    monkeypatch.setattr(core_api_module, "SmodesWrapper", RuntimeStub)
    source = _export_fail_closed_IsoDistort()
    with pytest.raises(ValueError, match="runtime clone failed") as caught:
        source.snapshot_for_action()
    assert len(created) == 4
    assert not created[0].closed and not created[1].closed
    assert created[2].closed and created[3].closed
    assert source._iso.generated == {"table": [1, 2, 3]}
    if cleanup_fails:
        assert "new worker cleanup failed" in str(caught.value.__cause__)

@pytest.mark.parametrize("cleanup_fails", [False, True])
def test_session_constructor_failure_closes_prior_wrappers_and_preserves_cause(
    monkeypatch, cleanup_fails,
):
    closed = []

    class FirstWrapper:
        def close(self):
            closed.append("iso")
            if cleanup_fails:
                raise OSError("cleanup failed")

    def fail_second_wrapper():
        raise ValueError("smodes initialization failed")

    monkeypatch.setattr(core_api_module, "IsoWrapper", FirstWrapper)
    monkeypatch.setattr(core_api_module, "SmodesWrapper", fail_second_wrapper)
    with pytest.raises(ValueError, match="smodes initialization failed") as caught:
        _export_fail_closed_IsoDistort()
    assert closed == ["iso"]
    if cleanup_fails:
        assert isinstance(caught.value.__cause__, RuntimeError)
        assert "cleanup failed" in str(caught.value.__cause__)

def test_session_close_attempts_every_wrapper_when_one_cleanup_fails():
    closed = []

    class Wrapper:
        def __init__(self, name, fail=False):
            self.name, self.fail = name, fail

        def close(self):
            closed.append(self.name)
            if self.fail:
                raise OSError("first cleanup failed")

    session = object.__new__(_export_fail_closed_IsoDistort)
    session._iso = Wrapper("iso", fail=True)
    session._findsym = Wrapper("findsym")
    session._smodes = Wrapper("smodes")
    with pytest.raises(RuntimeError, match="first cleanup failed"):
        session.close()
    assert closed == ["iso", "findsym", "smodes"]

def test_batch_export_aggregates_mode_failures_before_building_specs() -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    iso = _bare_export_session()
    built: list[SubgroupInfo] = []

    def fail_mode_search(subgroup_idx, **_kwargs):
        raise ValueError(f"mode failure {subgroup_idx}")

    iso.search_method_2 = fail_mode_search
    iso._spec_for_subgroup = lambda subgroup, **_kwargs: built.append(subgroup)

    with pytest.raises(ExportPreparationError) as caught:
        iso._collect_export_specs(candidates, ["modes"], True)

    assert built == []
    assert [failure.candidate_index for failure in caught.value.failures] == [0, 1]
    assert [failure.subgroup_index for failure in caught.value.failures] == [2, 2]
    assert all(
        failure.parent_space_group_number == 139
        and failure.k_coordinates == ("0", "0", "0")
        and failure.k_parameters == ()
        for failure in caught.value.failures
    )
    assert caught.value.failures[0].basis_vectors == (
        ("1", "0", "0"),
        ("0", "1", "0"),
        ("0", "0", "1"),
    )
    assert [failure.error_type for failure in caught.value.failures] == [
        "ValueError",
        "ValueError",
    ]
    assert iso._selected_subgroup is None
    assert iso.mode_displacements == {}

def test_successful_empty_mode_table_is_distinct_from_calculation_failure() -> None:
    candidate = _export_fail_closed_candidate(0, "GM1+")
    iso = _bare_export_session()

    def empty_mode_search(subgroup_idx, *, candidates, **_kwargs):
        iso._selected_subgroup = candidates[subgroup_idx]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}

    iso.search_method_2 = empty_mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    specs = iso._collect_export_specs([candidate], ["modes"], True)

    assert len(specs) == 1
    assert specs[0].use_current_modes is True
    assert specs[0].note == (
        "no displacement modes for this path (empty BUSH/smodes table)"
    )

def test_collection_can_skip_only_unresolved_identity_candidates() -> None:
    ready = _export_fail_closed_candidate(0, "GM1+")
    unresolved = _export_fail_closed_candidate(1, "GM2+")
    iso = _bare_export_session()
    events: list[dict] = []

    def mode_search(subgroup_idx, *, candidates, **_kwargs):
        candidate = candidates[0]
        iso._selected_subgroup = candidate
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        if subgroup_idx == unresolved.index:
            raise _export_fail_closed_UnresolvedModeIdentityError([
                ModeIdentityFailure(
                    amplitude_key="GM2+__a__unresolved(a)",
                    orbit_id="sg139:1a:Fe:test",
                    global_irrep="GM2+",
                    k_coordinates=("0", "0", "0"),
                    reason=(
                        "iso_microscopic_block_missing; parent_SG=139; "
                        "KVALUE=(); direction_selector=(a)"
                    ),
                )
            ])
        return SimpleNamespace(subgroup=candidate)

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    plan = iso._collect_export_specs(
        [ready, unresolved],
        ["modes"],
        True,
        allow_unresolved_identities=True,
        progress_callback=events.append,
    )

    assert isinstance(plan, ExportBatchPlan)
    assert [spec.subgroup for spec in plan.specs] == [ready]
    assert [failure.candidate_index for failure in plan.skipped] == [1]
    assert "iso_microscopic_block_missing" in plan.skipped[0].message
    assert "parent_SG=139" in plan.skipped[0].message
    completed = [
        event for event in events if event.get("phase") == "candidate_complete"
    ]
    assert completed[-1]["candidate_status"] == "ineligible"
    assert completed[-1]["candidate_error_type"] == "UnresolvedModeIdentityError"
    assert "iso_microscopic_block_missing" in completed[-1][
        "candidate_error_message"
    ]

def test_collection_cancels_at_candidate_boundary_and_keeps_completed_spec() -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    iso = _bare_export_session()
    cancelled = False
    events: list[dict] = []

    def mode_search(_subgroup_idx, *, candidates, **_kwargs):
        iso._selected_subgroup = candidates[0]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        return SimpleNamespace(subgroup=candidates[0])

    def progress(event):
        nonlocal cancelled
        events.append(event)
        if event.get("phase") == "candidate_complete":
            cancelled = True

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    plan = iso._collect_export_specs(
        candidates,
        ["modes"],
        True,
        allow_unresolved_identities=True,
        progress_callback=progress,
        cancel_check=lambda: cancelled,
    )

    assert isinstance(plan, ExportBatchPlan)
    assert [spec.subgroup for spec in plan.specs] == [candidates[0]]
    assert [failure.error_type for failure in plan.skipped] == [
        "ExportCancelledError"
    ]
    assert any(event["phase"] == "cancelled" for event in events)
    completed = next(
        event for event in events if event["phase"] == "candidate_complete"
    )
    assert completed["completed"] == 1
    assert completed["estimated_remaining_seconds"] >= 0
    assert completed["candidate_status"] == "ready"
    assert completed["candidate_error_type"] is None
    assert completed["candidate_error_message"] is None

def test_collection_uses_singleton_candidate_pools_for_duplicate_local_indices() -> None:
    first = _export_fail_closed_candidate(0, "GM1+")
    second = _export_fail_closed_candidate(0, "GM2+")
    iso = _bare_export_session()
    pools: list[tuple[SubgroupInfo, ...]] = []

    def mode_search(subgroup_idx, *, candidates, **_kwargs):
        assert subgroup_idx == 0
        pools.append(tuple(candidates))
        iso._selected_subgroup = candidates[0]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        return SimpleNamespace(subgroup=candidates[0])

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    specs = iso._collect_export_specs([first, second], ["modes"], True)

    assert pools == [(first,), (second,)]
    assert [spec.subgroup for spec in specs] == [first, second]

def test_collection_rejects_search_result_with_different_scientific_identity() -> None:
    requested = _export_fail_closed_candidate(0, "GM1+")
    wrong = _export_fail_closed_candidate(0, "GM2+")
    iso = _bare_export_session()

    def wrong_search(_subgroup_idx, **_kwargs):
        iso._selected_subgroup = wrong
        return SimpleNamespace(subgroup=wrong)

    iso.search_method_2 = wrong_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    with pytest.raises(ExportPreparationError) as caught:
        iso._collect_export_specs([requested], ["modes"], True)

    assert caught.value.failures[0].candidate_index == 0
    assert "different scientific subgroup" in caught.value.failures[0].message

def test_disk_export_does_not_create_root_when_preparation_fails(tmp_path) -> None:
    candidate = _export_fail_closed_candidate(0, "GM1+")
    failure = ExportCandidateFailure(
        position=1,
        candidate_index=0,
        subgroup_index=2,
        parent_space_group_number=139,
        k_point_label="GM",
        k_coordinates=("0", "0", "0"),
        k_parameters=(),
        irrep_label="GM1+",
        opd_symbol="P1",
        space_group_number=12,
        basis_vectors=(("1", "0", "0"), ("0", "1", "0"), ("0", "0", "1")),
        origin=("0.0", "0.0", "0.0"),
        embedding_id="embedding-test",
        error_type="ValueError",
        message="mode failure",
    )
    error = ExportPreparationError(2, [failure])
    iso = object.__new__(_export_fail_closed_IsoDistort)
    iso.structure = object()
    iso.subgroups = [candidate]

    def fail_collection(*_args, **_kwargs):
        raise error

    iso._collect_export_specs = fail_collection
    destination = tmp_path / "batch"

    with pytest.raises(ExportPreparationError):
        iso.export_subgroups(destination, formats=["modes"])

    assert not destination.exists()

def test_disk_batch_prerenders_every_candidate_before_publication(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    specs = [
        SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}")
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    iso.structure = object()
    iso.subgroups = candidates
    iso.symmetry_info = {"space_group_number": 139}
    iso._collect_export_specs = lambda *_args, **_kwargs: specs
    rendered: list[str] = []

    def render(spec, _formats):
        rendered.append(spec.folder_name)
        if spec is specs[1]:
            raise ValueError("later writer validation failed")
        return (("candidate-0.txt", b"ready"),)

    monkeypatch.setattr(publication_module, "render_subgroup_files", render)
    destination = tmp_path / "batch"

    with pytest.raises(ExportPreparationError, match="later writer validation failed"):
        iso.export_subgroups(destination, formats=["modes"])

    assert rendered == ["candidate-0", "candidate-1"]
    assert not destination.exists()

def test_anonymous_nonempty_modes_fail_during_collection() -> None:
    candidate = _export_fail_closed_candidate(4, "GM1+")
    iso = _bare_export_session()
    iso._selected_subgroup = candidate
    iso.mode_displacements_sc = {"anonymous": np.ones((1, 3))}
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=iso.mode_displacements_sc,
        displacive_data=None,
        **kwargs,
    )

    with pytest.raises(ExportPreparationError) as caught:
        iso._collect_export_specs([candidate], ["modes"], False)

    failure = caught.value.failures[0]
    assert failure.candidate_index == 4
    assert failure.subgroup_index == 2
    assert failure.parent_space_group_number == 139
    assert failure.k_coordinates == ("0", "0", "0")
    assert "lack verified identity" in failure.message

def test_supercell_modes_prevent_false_empty_mode_note() -> None:
    candidate = _export_fail_closed_candidate(0, "GM1+")
    iso = _bare_export_session()

    def mode_search(_subgroup_idx, *, candidates, **_kwargs):
        iso._selected_subgroup = candidates[0]
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {"mode": np.ones((1, 3))}

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=iso.mode_displacements_sc,
        displacive_data=object(),
        **kwargs,
    )

    specs = iso._collect_export_specs([candidate], ["modes"], True)

    assert len(specs) == 1
    assert specs[0].note == ""

def test_export_recompute_receives_current_types_and_scope_explicitly() -> None:
    candidate = _export_fail_closed_candidate(0, "GM1+")
    iso = _bare_export_session()
    iso.structure = Structure(
        Lattice.cubic(4.0),
        ["Fe", "O"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    iso.distortion_types = ["displacive"]
    iso.distortion_scope = {"displacive": ["Fe"]}
    iso._selected_subgroup = candidate
    old_scope = {"displacive": ["O"]}
    iso._mode_cache_key = iso._mode_context_key(
        candidate,
        0,
        ["displacive"],
        old_scope,
    )
    calls: list[dict] = []

    def mode_search(_subgroup_idx, **kwargs):
        calls.append(kwargs)
        iso._selected_subgroup = candidate
        iso.mode_displacements = {}
        iso.mode_occupancies = {}
        iso.mode_displacements_sc = {}
        return SimpleNamespace(subgroup=candidate)

    iso.search_method_2 = mode_search
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    specs = iso._collect_export_specs([candidate], ["modes"], True)

    assert len(specs) == 1
    assert calls == [
        {
            "distortion_type": ["displacive"],
            "number_of_independent_modulations": 0,
            "candidates": [candidate],
            "distortion_scope": {"displacive": ["Fe"]},
        }
    ]
    assert specs[0].use_current_modes is True
    assert specs[0].use_generated_structure is False

@pytest.mark.parametrize(
    ("types", "scope"),
    (
        (["strain", "displacive"], {"displacive": ["O"]}),
        (["strain"], {"displacive": ["Fe"]}),
    ),
)
def test_cif_only_export_rejects_generated_structure_from_old_context(
    types,
    scope,
) -> None:
    candidate = _export_fail_closed_candidate(0, "GM1+")
    iso = _bare_export_session()
    iso.structure = Structure(
        Lattice.cubic(4.0),
        ["Fe", "O"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    iso._selected_subgroup = candidate
    iso.distortion_types = types
    iso.distortion_scope = scope
    old_key = iso._mode_context_key(
        candidate,
        0,
        ["strain", "displacive"],
        {"displacive": ["Fe"]},
    )
    iso._mode_cache_key = old_key
    iso._generated_structure_cache_key = old_key
    iso.distorted_structure = object()
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    [spec] = iso._collect_export_specs([candidate], ["cif"], False)

    assert spec.use_current_modes is False
    assert spec.use_generated_structure is False

def test_cif_only_export_reuses_generated_structure_for_exact_context() -> None:
    candidate = _export_fail_closed_candidate(0, "GM1+")
    iso = _bare_export_session()
    iso.structure = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    iso._selected_subgroup = candidate
    iso.distortion_types = ["strain", "displacive"]
    iso.distortion_scope = {"displacive": ["Fe"]}
    key = iso._mode_context_key(
        candidate,
        0,
        iso.distortion_types,
        iso.distortion_scope,
    )
    iso._mode_cache_key = key
    iso._generated_structure_cache_key = key
    iso.distorted_structure = object()
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup,
        mode_displacements_sc=None,
        displacive_data=None,
        **kwargs,
    )

    [spec] = iso._collect_export_specs([candidate], ["cif"], False)

    assert spec.use_current_modes is False
    assert spec.use_generated_structure is True

def test_zip_prerenders_and_aggregates_all_writer_failures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    specs = [
        SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}")
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    iso.structure = object()
    iso.subgroups = candidates
    iso._collect_export_specs = lambda *_args, **_kwargs: specs
    archive_called = False

    def render(spec, _formats):
        raise ValueError(f"writer failed for {spec.folder_name}")

    def archive(*_args, **_kwargs):
        nonlocal archive_called
        archive_called = True
        return b""

    monkeypatch.setattr(publication_module, "render_subgroup_files", render)
    monkeypatch.setattr(core_api_module, "build_export_zip", archive)

    with pytest.raises(ExportPreparationError) as caught:
        iso.export_subgroups_zip(formats=["modes"])

    assert archive_called is False
    assert [failure.candidate_index for failure in caught.value.failures] == [0, 1]
    assert [failure.irrep_label for failure in caught.value.failures] == [
        "GM1+",
        "GM2+",
    ]
    assert all("writer failed" in failure.message for failure in caught.value.failures)

def test_zip_uses_each_prerendered_candidate_payload_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    specs = [
        SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}")
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    iso.structure = object()
    iso.subgroups = candidates
    iso._collect_export_specs = lambda *_args, **_kwargs: specs
    calls: list[str] = []

    def render(spec, _formats):
        calls.append(spec.folder_name)
        return (("result.txt", spec.folder_name.encode("ascii")),)

    monkeypatch.setattr(publication_module, "render_subgroup_files", render)

    archive = iso.export_subgroups_zip(formats=["modes"])

    assert calls == ["candidate-0", "candidate-1"]
    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        assert zipped.namelist() == [
            "candidate-0/result.txt",
            "candidate-1/result.txt",
            "export_candidate_status.json",
            "export_candidate_status.txt",
        ]
        assert zipped.read("candidate-1/result.txt") == b"candidate-1"
        report = json.loads(zipped.read("export_candidate_status.json"))
        assert report == {
            "schema": "isodistort.export-candidate-status.v1",
            "method": 2,
            "requested_candidate_count": 2,
            "successful_candidate_count": 2,
            "ineligible_candidate_count": 0,
            "unresolved_identity_candidate_count": 0,
            "cancelled_candidate_count": 0,
            "status": "success",
            "failures": [],
        }
        assert "requested=2 successful=2 ineligible=0" in zipped.read(
            "export_candidate_status.txt"
        ).decode("utf-8")

def test_method2_zip_keeps_ready_candidate_and_embeds_skipped_identity_report(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ready = _export_fail_closed_candidate(0, "GM1+")
    skipped = _export_fail_closed_candidate(1, "GM2+")
    spec = SimpleNamespace(subgroup=ready, folder_name="candidate-ready")
    failure = _export_fail_closed_IsoDistort._export_candidate_failure(
        _bare_export_session(),
        2,
        skipped,
        _export_fail_closed_UnresolvedModeIdentityError([
            ModeIdentityFailure(
                amplitude_key="GM2+__a__unresolved(a)",
                orbit_id="sg139:1a:Fe:test",
                global_irrep="GM2+",
                k_coordinates=("0", "0", "0"),
                reason="candidate_columns_not_independent; parent_SG=139",
            )
        ]),
    )
    iso = object.__new__(_export_fail_closed_IsoDistort)
    iso.structure = object()
    iso.subgroups = [ready, skipped]
    iso._collect_export_specs = lambda *_args, **_kwargs: ExportBatchPlan(
        (spec,), (failure,)
    )
    monkeypatch.setattr(
        publication_module,
        "render_subgroup_files",
        lambda _spec, _formats: (("result.txt", b"ready"),),
    )

    archive = iso.export_subgroups_zip(
        formats=["modes"],
        compute_missing_modes=True,
        export_method=2,
    )

    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        assert zipped.read("candidate-ready/result.txt") == b"ready"
        report = json.loads(zipped.read("export_candidate_status.json"))
        assert report["status"] == "partial_success"
        assert report["successful_candidate_count"] == 1
        assert report["ineligible_candidate_count"] == 1
        assert report["failures"][0]["candidate_index"] == 1
        assert "candidate_columns_not_independent" in zipped.read(
            "export_candidate_status.txt"
        ).decode("utf-8")

def test_method2_disk_export_atomically_publishes_ready_candidate_and_report(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ready = _export_fail_closed_candidate(0, "GM1+")
    skipped = _export_fail_closed_candidate(1, "GM2+")
    spec = SimpleNamespace(subgroup=ready, folder_name="candidate-ready")
    failure = _export_fail_closed_IsoDistort._export_candidate_failure(
        _bare_export_session(),
        2,
        skipped,
        _export_fail_closed_UnresolvedModeIdentityError([
            ModeIdentityFailure(
                amplitude_key="GM2+__a__unresolved(a)",
                orbit_id="sg139:1a:Fe:test",
                global_irrep="GM2+",
                k_coordinates=("0", "0", "0"),
                reason="rank_mismatch; parent_SG=139; KVALUE=()",
            )
        ]),
    )
    iso = object.__new__(_export_fail_closed_IsoDistort)
    iso.structure = object()
    iso.subgroups = [ready, skipped]
    iso._collect_export_specs = lambda *_args, **_kwargs: ExportBatchPlan(
        (spec,), (failure,)
    )
    monkeypatch.setattr(
        publication_module,
        "render_subgroup_files",
        lambda _spec, _formats: (("result.txt", b"ready"),),
    )
    destination = tmp_path / "batch"

    paths = iso.export_subgroups(
        destination,
        formats=["modes"],
        compute_missing_modes=True,
        export_method=2,
    )

    assert (destination / "candidate-ready" / "result.txt").read_bytes() == b"ready"
    report = json.loads(
        (destination / "export_candidate_status.json").read_text(encoding="utf-8")
    )
    assert report["status"] == "partial_success"
    assert report["successful_candidate_count"] == 1
    assert report["failures"][0]["candidate_index"] == 1
    assert set(paths) == {
        destination / "candidate-ready" / "result.txt",
        destination / "export_candidate_status.json",
        destination / "export_candidate_status.txt",
    }

def test_disk_batch_staging_write_failure_publishes_nothing(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}"),
            (("result.txt", b"ready"),),
        )
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    original_write = core_api_module.Path.write_bytes

    def fail_second_write(path, payload):
        if path.parent.name == "candidate-1":
            raise OSError("staged write failed")
        return original_write(path, payload)

    monkeypatch.setattr(core_api_module.Path, "write_bytes", fail_second_write)
    destination = tmp_path / "batch"

    with pytest.raises(ExportPreparationError) as caught:
        publish_export_batch(2, destination, rendered, iso._export_candidate_failure)

    assert caught.value.failures[0].candidate_index == 1
    assert "staged write failed" in caught.value.failures[0].message
    assert not destination.exists()
    assert not list(tmp_path.glob(".isodistort-batch-*"))

def test_existing_disk_root_commit_failure_keeps_the_entire_batch_invisible(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}"),
            (("result.txt", b"ready"),),
        )
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    sentinel = destination / "keep.txt"
    sentinel.write_bytes(b"unchanged")
    original_rename = core_api_module.Path.rename

    def fail_batch_commit(path, target):
        if (
            path.name.startswith(".isodistort-batch-")
            and target.parent == destination
            and target.name.endswith(".ready")
        ):
            raise OSError("batch commit failed")
        return original_rename(path, target)

    monkeypatch.setattr(core_api_module.Path, "rename", fail_batch_commit)

    with pytest.raises(ExportPreparationError) as caught:
        publish_export_batch(2, destination, rendered, iso._export_candidate_failure)

    assert [failure.candidate_index for failure in caught.value.failures] == [0, 1]
    assert all(
        "batch commit failed" in failure.message for failure in caught.value.failures
    )
    assert list(destination.iterdir()) == [sentinel]
    assert sentinel.read_bytes() == b"unchanged"
    assert not list(tmp_path.glob(".isodistort-batch-*"))

def test_existing_disk_root_exposes_a_complete_version_in_one_transition(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name=f"candidate-{index}"),
            (("result.txt", f"ready-{index}".encode("ascii")),),
        )
        for index, candidate in enumerate(candidates)
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"unchanged")
    original_rename = core_api_module.Path.rename
    commit_entered = threading.Event()
    committed = threading.Event()
    observed_before = threading.Event()
    observed_after = threading.Event()
    stop_reader = threading.Event()
    snapshots: list[tuple[tuple[str, ...], ...]] = []

    def ready_snapshot() -> tuple[tuple[str, ...], ...]:
        batches = sorted(
            path
            for path in destination.iterdir()
            if path.name.startswith(publication_module._EXPORT_BATCH_PREFIX)
            and path.name.endswith(publication_module._EXPORT_BATCH_READY_SUFFIX)
        )
        return tuple(
            tuple(
                sorted(
                    entry.name
                    for entry in batch.iterdir()
                    if entry.name != publication_module._EXPORT_BATCH_MANIFEST
                )
            )
            for batch in batches
        )

    def read_while_publishing() -> None:
        while not stop_reader.is_set():
            if not commit_entered.is_set():
                stop_reader.wait(0.001)
                continue
            snapshot = ready_snapshot()
            snapshots.append(snapshot)
            if not snapshot:
                observed_before.set()
            if committed.is_set() and snapshot:
                observed_after.set()
            stop_reader.wait(0.001)

    def observe_single_commit(path, target):
        if (
            path.name.startswith(".isodistort-batch-")
            and target.parent == destination
            and target.name.endswith(".ready")
        ):
            commit_entered.set()
            assert observed_before.wait(5)
            result = original_rename(path, target)
            committed.set()
            assert observed_after.wait(5)
            return result
        return original_rename(path, target)

    monkeypatch.setattr(core_api_module.Path, "rename", observe_single_commit)
    reader = threading.Thread(target=read_while_publishing)
    reader.start()
    try:
        paths = publish_export_batch(2, destination, rendered, iso._export_candidate_failure)
    finally:
        stop_reader.set()
        reader.join(timeout=5)

    complete = (("candidate-0", "candidate-1"),)
    assert snapshots
    assert () in snapshots
    assert complete in snapshots
    assert all(snapshot in {(), complete} for snapshot in snapshots)
    assert len({path.parents[1] for path in paths}) == 1
    assert all(path.read_bytes().startswith(b"ready-") for path in paths)
    assert publication_module._published_export_folder_names(destination) >= {
        "candidate-0",
        "candidate-1",
    }
    assert not list(tmp_path.glob(".isodistort-batch-*"))

def test_stale_staging_directory_is_not_discoverable_as_a_committed_batch(
    tmp_path,
) -> None:
    candidate = _export_fail_closed_candidate(0, "GM1+")
    rendered = [
        (
            SimpleNamespace(subgroup=candidate, folder_name="candidate-0"),
            (("result.txt", b"ready"),),
        )
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    stale = destination / ".isodistort-batch-abandoned"
    (stale / "candidate-stale").mkdir(parents=True)
    (stale / "candidate-stale" / "result.txt").write_bytes(b"partial")

    before = publication_module._published_export_folder_names(destination)
    paths = publish_export_batch(2, destination, rendered, iso._export_candidate_failure)
    after = publication_module._published_export_folder_names(destination)

    assert "candidate-stale" not in before
    assert "candidate-stale" not in after
    assert "candidate-0" in after
    assert paths[0].read_bytes() == b"ready"
    ready_batches = [
        path
        for path in destination.iterdir()
        if path.name.startswith(publication_module._EXPORT_BATCH_PREFIX)
        and path.name.endswith(publication_module._EXPORT_BATCH_READY_SUFFIX)
    ]
    assert len(ready_batches) == 1

@pytest.mark.parametrize(
    ("tamper", "message"),
    (
        ("manifest", "manifest does not match"),
        ("payload", "file hash changed"),
    ),
)
def test_tampered_ready_batch_fails_closed_before_another_publish(
    tmp_path,
    tamper: str,
    message: str,
) -> None:
    candidates = [_export_fail_closed_candidate(0, "GM1+"), _export_fail_closed_candidate(1, "GM2+")]
    first = [
        (
            SimpleNamespace(subgroup=candidates[0], folder_name="candidate-0"),
            (("result.txt", b"ready"),),
        )
    ]
    second = [
        (
            SimpleNamespace(subgroup=candidates[1], folder_name="candidate-1"),
            (("result.txt", b"later"),),
        )
    ]
    iso = object.__new__(_export_fail_closed_IsoDistort)
    destination = tmp_path / "batch"
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"unchanged")
    [published] = publish_export_batch(2, destination, first, iso._export_candidate_failure)
    ready_root = published.parents[1]
    if tamper == "manifest":
        manifest = ready_root / publication_module._EXPORT_BATCH_MANIFEST
        manifest.write_bytes(manifest.read_bytes() + b" ")
    else:
        published.write_bytes(b"wrong")

    with pytest.raises(ExportPreparationError) as caught:
        publish_export_batch(2, destination, second, iso._export_candidate_failure)

    assert message in caught.value.failures[0].message
    assert not any(path.name == "candidate-1" for path in destination.rglob("*"))
    assert not any(
        path.name.startswith(".isodistort-batch-")
        and not path.name.endswith(publication_module._EXPORT_BATCH_READY_SUFFIX)
        for path in destination.iterdir()
    )

def test_concurrent_writers_reject_a_partially_overlapping_candidate_set(
    tmp_path,
) -> None:
    candidate_sets = (
        [
            (_export_fail_closed_candidate(0, "GM1+"), "candidate-shared", b"writer-a-shared"),
            (_export_fail_closed_candidate(1, "GM2+"), "candidate-a", b"writer-a-only"),
        ],
        [
            (_export_fail_closed_candidate(2, "GM3+"), "candidate-shared", b"writer-b-shared"),
            (_export_fail_closed_candidate(3, "GM4+"), "candidate-b", b"writer-b-only"),
        ],
    )
    destination = tmp_path / "batch"
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"unchanged")
    start = threading.Barrier(3)
    successes: list[list] = []
    failures: list[ExportPreparationError] = []

    def publish(entries) -> None:
        rendered = [
            (
                SimpleNamespace(subgroup=candidate, folder_name=folder),
                (("result.txt", payload),),
            )
            for candidate, folder, payload in entries
        ]
        start.wait()
        try:
            successes.append(
                publish_export_batch(
                    2,
                    destination,
                    rendered,
                    object.__new__(_export_fail_closed_IsoDistort)._export_candidate_failure,
                )
            )
        except ExportPreparationError as exc:
            failures.append(exc)

    writers = [
        threading.Thread(target=publish, args=(entries,))
        for entries in candidate_sets
    ]
    for writer in writers:
        writer.start()
    start.wait()
    for writer in writers:
        writer.join(timeout=10)

    assert all(not writer.is_alive() for writer in writers)
    assert len(successes) == 1
    assert len(failures) == 1
    assert [failure.message for failure in failures[0].failures] == [
        f"target candidate directory exists in a committed batch: "
        f"{destination / 'candidate-shared'}"
    ]
    ready_batches = [
        path
        for path in destination.iterdir()
        if path.name.startswith(publication_module._EXPORT_BATCH_PREFIX)
        and path.name.endswith(publication_module._EXPORT_BATCH_READY_SUFFIX)
    ]
    assert len(ready_batches) == 1
    visible = {
        path.name
        for path in ready_batches[0].iterdir()
        if path.name != publication_module._EXPORT_BATCH_MANIFEST
    }
    assert visible in (
        {"candidate-shared", "candidate-a"},
        {"candidate-shared", "candidate-b"},
    )
    assert not any(
        path.name.startswith(".isodistort-batch-")
        and not path.name.endswith(publication_module._EXPORT_BATCH_READY_SUFFIX)
        for path in destination.iterdir()
    )

def test_publish_lock_serializes_independent_processes(tmp_path) -> None:
    destination = tmp_path / "batch"
    destination.mkdir()
    first_entered = tmp_path / "first-entered"
    second_entered = tmp_path / "second-entered"
    release_first = tmp_path / "release-first"
    # ``python -c`` does not put the working directory on sys.path, so the
    # child process needs the project root explicitly.
    project_root = TESTS_ROOT.parent
    child_env = dict(os.environ)
    child_env["PYTHONPATH"] = os.pathsep.join(
        [str(project_root), child_env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    script = """
from pathlib import Path
import sys
import time
from application.session.publication import _export_publish_lock

root = Path(sys.argv[1])
entered = Path(sys.argv[2])
release = None if sys.argv[3] == "-" else Path(sys.argv[3])
with _export_publish_lock(root):
    entered.write_text("entered", encoding="ascii")
    while release is not None and not release.exists():
        time.sleep(0.01)
"""
    first = subprocess.Popen(  # noqa: S603 - exact current interpreter
        [sys.executable, "-c", script, str(destination), str(first_entered), str(release_first)],
        cwd=project_root,
        env=child_env,
    )
    deadline = time.monotonic() + 5
    while not first_entered.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    assert first_entered.exists()

    second = subprocess.Popen(  # noqa: S603 - exact current interpreter
        [sys.executable, "-c", script, str(destination), str(second_entered), "-"],
        cwd=project_root,
        env=child_env,
    )
    time.sleep(0.2)
    assert not second_entered.exists()

    release_first.write_text("release", encoding="ascii")
    assert first.wait(timeout=5) == 0
    assert second.wait(timeout=5) == 0
    assert second_entered.exists()

# Original coverage: test_formats.py.

def test_strain_only_gamma_keeps_gm4_not_gm3():
    # This exercises the exact strain representation only.  Bypass the public
    # constructor so an unavailable WSL wrapper cannot turn a pure unit test
    # into an environment-dependent integration test.
    api = object.__new__(_formats_IsoDistort)
    api._strain_representation_cache = None
    api.structure = Structure.from_spacegroup(
        "I4/mmm", Lattice.tetragonal(4.0, 10.0), ["Eu"], [[0, 0, 0]],
    )
    gm4 = SimpleNamespace(subgroup=_formats_SubgroupInfo(
        index=0,
        parent_sg=139,
        k_point_label="GM",
        irrep_label="GM4+",
        size=1,
        subgroup_index=2,
        space_group_number=69,
        space_group_symbol="Fmmm",
        basis_vectors=[[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
    ))
    gm3 = SimpleNamespace(subgroup=_formats_SubgroupInfo(
        index=1,
        parent_sg=139,
        k_point_label="GM",
        irrep_label="GM3+",
        size=1,
        subgroup_index=2,
        space_group_number=87,
        space_group_symbol="I4/m",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
    ))
    mirror_only = SimpleNamespace(subgroup=_formats_SubgroupInfo(
        index=2,
        parent_sg=139,
        k_point_label="GM",
        # Synthetic even label isolates the representation test: Cs fixes
        # more tensors than D4h, but inversion acts trivially on every strain
        # tensor, so its exact pointwise stabilizer is strictly larger than Cs.
        irrep_label="GM_FAKE+",
        size=1,
        subgroup_index=8,
        space_group_number=8,
        space_group_symbol="Cm",
        basis_vectors=[[-1, 1, 0], [0, 0, 1], [1, 0, 0]],
        origin=[0, 0, 0],
    ))
    assert api._keep_strain_only_irrep(gm4, 139)
    assert not api._keep_strain_only_irrep(gm3, 139)
    assert not api._keep_strain_only_irrep(mirror_only, 139)

def test_method3_api_zip_uses_generated_case_and_short_candidate_folders():
    parent = _cubic()
    first = _sg("GM1+", "P1", index=0, number=139)
    second = _sg("N1+", "4D1", index=1, number=2)
    iso = _minimal_export_api(parent)

    raw = iso.export_subgroups_zip(
        formats=[FORMAT_CIF],
        subgroups=[first, second],
        export_method=3,
    )
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()
        case = method3_case_folder([first, second])
        status = json.loads(zf.read(f"{case}/export_candidate_status.json"))

    assert names == [
        f"{case}/C01_SG139/subgroup.cif",
        f"{case}/C02_SG2/subgroup.cif",
        f"{case}/export_candidate_status.json",
        f"{case}/export_candidate_status.txt",
    ]
    assert status["status"] == "success"
    assert status["successful_candidate_count"] == status["requested_candidate_count"] == 2

# Original coverage: test_parent_setting_guard.py.

def _diagnostic_api(parent=None):
    if parent is None:
        parent = Structure(
            Lattice.orthorhombic(3, 4, 5), ["Fe", "Fe"],
            [[0, 0, 0], [0, 0.5, 0.5]],
        )
    api = object.__new__(_parent_setting_guard_IsoDistort)
    api.cfg = get_config()
    api.structure = None
    api.structure_path = None
    api.symmetry_info = None
    api._sym_val = SymmetryValidator()
    api._findsym = SimpleNamespace(standardize_wyckoff_orbits=Mock(
        side_effect=OSError("intentional FINDSYM environment failure"),
    ))
    api.set_structure(parent)
    api.distortion_scope = {}
    api.distortion_types = ["displacive"]
    api.number_of_independent_modulations = 7
    return api

def test_failed_setting_verification_keeps_the_upload_available_for_diagnostics():
    api = _diagnostic_api()
    assert api.symmetry_info["space_group_number"] == 65
    assert api.symmetry_info["parent_orbit_standardization"]["status"] == "unavailable"
    assert np.array_equal(api.structure.frac_coords, [[0, 0, 0], [0, 0.5, 0.5]])

def _invoke_mode_entry(api, entry):
    if entry == "scoped_displacive":
        return api._compute_scoped_modes(["displacive"], raw_modes=[])
    if entry == "scoped_occupational":
        return api._compute_scoped_modes(["occupational"], raw_modes=[])
    if entry == "scoped_supplied_raw":
        return api._compute_scoped_modes(["displacive"], raw_modes=[object()])
    if entry == "select":
        return api.select_path(0)
    return api.search_method_2(0, number_of_independent_modulations=1)

@pytest.mark.parametrize("entry", [
    "scoped_displacive", "scoped_occupational", "scoped_supplied_raw",
    "select", "search",
])
def test_unavailable_setting_blocks_all_mode_entries_before_calls_or_mutation(monkeypatch, entry):
    api = _diagnostic_api()
    api._iso = Mock(name="ISO")
    api._smodes = Mock(name="SMODES")
    api._method1 = Mock(name="Method1")
    api._method2 = Mock(name="Method2")
    api._method3 = Mock(name="Method3")
    api._dist_mapper = Mock(name="Mapper")
    api.phase_path = object()
    api._selected_subgroup = object()
    api.distortion_modes = [object()]
    api.mode_displacements = {"existing": object()}
    api.mode_displacements_sc = {"existing": object()}
    api.mode_occupancies = {"existing": object()}
    api._mode_label_overrides = {"existing": "existing label"}
    api._mode_cache_key = ("existing",)
    api.distorted_structure = api.structure.copy()
    api._generated_structure_cache_key = ("existing",)
    parametric = Mock(name="Parametric", side_effect=AssertionError("must not calculate modes"))
    special = Mock(name="Special", side_effect=AssertionError("must not calculate modes"))
    monkeypatch.setattr(superspace, "compute_parametric_modes", parametric)
    monkeypatch.setattr(superspace, "compute_special_modes_with_rootless_supplement", special)
    before = dict(vars(api))
    before_metadata = copy.deepcopy(api.symmetry_info)

    with pytest.raises(IsodistortError, match=r"Parent setting has not been verified.*reload"):
        _invoke_mode_entry(api, entry)

    assert set(vars(api)) == set(before)
    for key, value in before.items():
        assert vars(api)[key] is value
    assert api.symmetry_info == before_metadata
    for backend in (api._iso, api._smodes, api._method1, api._method2, api._method3, api._dist_mapper,
                    parametric, special):
        assert backend.mock_calls == []

@pytest.mark.parametrize("status", ["available", "canonicalized"])
@pytest.mark.parametrize("entry", ["scoped_displacive", "scoped_occupational"])
def test_verified_scope_filter_consumes_modes_without_calculating(monkeypatch, status, entry):
    # Empty downstream stubs prove only that the gate allows these statuses.
    # They do not purport to certify the deliberately nonstandard test parent.
    api = _diagnostic_api()
    api.symmetry_info["parent_orbit_standardization"] = {"status": status}
    api._iso = SimpleNamespace(list_k_points=Mock(return_value=[]))
    api._smodes = object()
    parametric = Mock(return_value=superspace.ParametricModeResult(
        modes=[], supercell_displacements={}, labels={}, nmod=0,
    ))
    monkeypatch.setattr(superspace, "compute_parametric_modes", parametric)
    if entry == "scoped_displacive":
        assert _invoke_mode_entry(api, entry) == []
        assert parametric.call_count == 0
        assert api._iso.list_k_points.call_count == 0
    else:
        # A verified setting cannot supply a missing scalar representation.
        before = dict(vars(api))
        with pytest.raises(NotImplementedError, match="occupational"):
            _invoke_mode_entry(api, entry)
        assert set(vars(api)) == set(before)
        assert all(vars(api)[key] is value for key, value in before.items())
        assert parametric.mock_calls == []

@pytest.mark.parametrize("status", ["available", "canonicalized"])
@pytest.mark.parametrize("entry", ["select", "search"])
def test_verified_status_reaches_the_normal_public_entry_after_guard(status, entry):
    api = _diagnostic_api()
    api.symmetry_info["parent_orbit_standardization"] = {"status": status}
    # This case isolates the setting guard; the public selection contract also
    # rejects unsupported nmod values before attempting candidate discovery.
    api.number_of_independent_modulations = 1
    downstream = RuntimeError("normal candidate listing reached")
    api.list_subgroups = Mock(side_effect=downstream)
    with pytest.raises(RuntimeError, match="normal candidate listing reached"):
        _invoke_mode_entry(api, entry)
    assert api.list_subgroups.call_count == 1

def test_legacy_missing_setting_metadata_is_not_reclassified_as_verified():
    api = _diagnostic_api()
    del api.symmetry_info["parent_orbit_standardization"]
    before = copy.deepcopy(api.symmetry_info)
    api._require_verified_parent_setting()
    assert api.symmetry_info == before
    assert "parent_orbit_standardization" not in api.symmetry_info

def _setting_boundary_candidate(parent_sg=65, size=2, group_index=8):
    return SubgroupInfo(
        index=0, space_group_number=2, space_group_symbol="P-1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0], parent_sg=parent_sg, size=size,
        subgroup_index=group_index,
    )

_SETTING_BOUNDARIES = ["search1", "search3", "zip", "disk", "spec", "collect", "supercell"]

def _invoke_setting_boundary(api, entry, dest, progress, cancel):
    candidate = _setting_boundary_candidate()
    if entry == "search1":
        return api.search_method_1(distortion_types=["occupational"])
    if entry == "search3":
        return api.search_method_3(
            distortion_types=["occupational"], space_group_type=2,
            supercell_basis=candidate.basis_vectors,
        )
    if entry == "zip":
        return api.export_subgroups_zip(
            formats=["cif"], subgroups=[candidate], compute_missing_modes=False,
            export_method=3, progress_callback=progress, cancel_check=cancel,
        )
    if entry == "disk":
        return api.export_subgroups(
            dest, formats=["cif"], subgroups=[candidate], compute_missing_modes=False,
            export_method=3, progress_callback=progress, cancel_check=cancel,
        )
    if entry == "spec":
        return api._spec_for_subgroup(
            candidate, use_current_modes=False, use_generated_structure=False,
        )
    if entry == "collect":
        return api._collect_export_specs(
            [candidate], ["cif"], False, export_method=3,
            progress_callback=progress, cancel_check=cancel,
        )
    return api._supercell_for_subgroup(candidate)

def _intercept_first_downstream(monkeypatch, api, entry, downstream):
    if entry in {"search1", "search3"}:
        monkeypatch.setattr(api, "_parent_rotations", downstream)
    elif entry == "zip":
        monkeypatch.setattr(api, "_collect_export_specs", downstream)
    elif entry == "disk":
        # The public disk entry must stop even before inspecting an export
        # destination, not only before writing or before preparing candidates.
        monkeypatch.setattr(core_api, "_published_export_folder_names", downstream)
    elif entry == "spec":
        monkeypatch.setattr(api, "_supercell_for_subgroup", downstream)
    elif entry == "collect":
        monkeypatch.setattr(api, "_snapshot_distortion_state", downstream)
    else:
        monkeypatch.setattr(core_api, "build_supercell", downstream)

@pytest.mark.parametrize("entry", _SETTING_BOUNDARIES)
def test_unavailable_setting_blocks_search_and_cif_exports_before_any_effect(
    monkeypatch, tmp_path, entry,
):
    api = _diagnostic_api()
    api._iso = Mock(name="ISO")
    api._smodes = Mock(name="SMODES")
    api._method1 = Mock(name="Method1")
    api._method2 = Mock(name="Method2")
    api._method3 = Mock(name="Method3")
    api._dist_mapper = Mock(name="Mapper")
    downstream = Mock(side_effect=AssertionError("unverified setting reached downstream"))
    _intercept_first_downstream(monkeypatch, api, entry, downstream)
    progress = Mock(side_effect=AssertionError("must not call progress"))
    cancel = Mock(side_effect=AssertionError("must not call cancellation observer"))
    dest = tmp_path / "never-created"
    before = dict(vars(api))
    before_metadata = copy.deepcopy(api.symmetry_info)
    before_structure = api.structure.copy()

    with pytest.raises(IsodistortError, match=r"Parent setting has not been verified.*reload"):
        _invoke_setting_boundary(api, entry, dest, progress, cancel)

    assert not dest.exists()
    assert list(tmp_path.iterdir()) == []
    assert set(vars(api)) == set(before)
    assert all(vars(api)[key] is value for key, value in before.items())
    assert api.symmetry_info == before_metadata
    assert api.structure == before_structure
    for observer in (downstream, progress, cancel, api._iso, api._smodes, api._method1, api._method2, api._method3, api._dist_mapper):
        assert observer.mock_calls == []

@pytest.mark.parametrize("status", ["available", "canonicalized"])
@pytest.mark.parametrize("entry", _SETTING_BOUNDARIES)
def test_verified_setting_reaches_the_normal_search_or_export_boundary(
    monkeypatch, tmp_path, status, entry,
):
    # Deliberate downstream stops prove gate reachability only; assigning a
    # status here is not a scientific certification of the diagnostic parent.
    api = _diagnostic_api()
    api.symmetry_info["parent_orbit_standardization"] = {"status": status}
    downstream = Mock(side_effect=RuntimeError("normal downstream reached"))
    _intercept_first_downstream(monkeypatch, api, entry, downstream)
    progress = Mock()
    cancel = Mock()

    with pytest.raises(RuntimeError, match="normal downstream reached"):
        _invoke_setting_boundary(api, entry, tmp_path / "never-created", progress, cancel)

    assert downstream.call_count == 1
    assert progress.mock_calls == []
    assert cancel.mock_calls == []
    assert list(tmp_path.iterdir()) == []

def test_unverified_shifted_parent_cannot_publish_a_cif_that_adds_inversion_images():
    """A valid origin-shifted parent is not in the standard ISO subgroup frame.

    Pmmm has one Ni at this shifted inversion center. If the unchecked writer
    applies standard P-1 at O=0, its CIF creates a second Ni at -x. Only the
    FINDSYM runtime failure is simulated; structure detection and any reached
    ZIP/CIF rendering and parsing are real, with no filesystem publication.
    """
    parent = Structure(Lattice.orthorhombic(3, 4, 5), ["Ni"], [[0.1, 0.2, 0.3]])
    api = _diagnostic_api(parent)
    assert api.symmetry_info["space_group_number"] == 47
    assert api.symmetry_info["parent_orbit_standardization"]["status"] == "unavailable"
    candidate = _setting_boundary_candidate(parent_sg=47, size=1, group_index=4)
    try:
        archive = api.export_subgroups_zip(
            formats=["cif"], subgroups=[candidate], compute_missing_modes=False,
            export_method=3,
        )
    except IsodistortError as exc:
        assert "Parent setting has not been verified" in str(exc)
        assert "reload" in str(exc)
        assert api.structure is parent
        return

    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        cif_names = [name for name in zipped.namelist() if name.endswith("subgroup.cif")]
        assert len(cif_names) == 1
        text = zipped.read(cif_names[0]).decode("utf-8")
    daughter = CifParser.from_str(text).parse_structures(primitive=False)[0]
    pytest.fail(
        "Unverified parent frame published a CIF: "
        f"{parent.composition} ({len(parent)} atoms) -> "
        f"{daughter.composition} ({len(daughter)} atoms)"
    )

# Original coverage: test_partial_microscopic_mapping.py.

WORKSPACE = TESTS_ROOT.parent.parent

def test_method1_candidate_method2_flow_forwards_active_arm_context():
    """The selected Method-1 candidate remains the mapping proof context."""

    parent, sites, subgroup, mode, _basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mapper = _partial_microscopic_mapping_DistortionMapper()
    mapper.map_modes_to_atoms = lambda *_args, **_kwargs: {}

    api = object.__new__(_partial_microscopic_mapping_IsoDistort)
    api.structure = parent
    api.symmetry_info = {
        "space_group_number": 123,
        "space_group_symbol": "P4/mmm",
        "wyckoff_sites": sites,
    }
    api.cfg = get_config()
    api.subgroups = []
    api._unresolved_method3_embedding_keys = set()
    api.distortion_types = ["displacive"]
    api.distortion_scope = {"displacive": ["*"]}
    api._selected_subgroup = None
    api._dist_mapper = mapper
    api._iso = SimpleNamespace(list_k_points=lambda *_args, **_kwargs: [])
    api._smodes = None
    api._method2 = SimpleNamespace(
        search=lambda *_args, **_kwargs: SimpleNamespace(
            subgroup=subgroup,
            modes=[mode],
            metadata={},
        )
    )
    api._compute_scoped_modes = lambda *_args, **_kwargs: [mode]
    api._resolve_k_vector = lambda *_args, **_kwargs: [0.0, 0.5, 0.0]
    api._mode_context_key = lambda *_args, **_kwargs: ("test-context",)
    api.mode_occupancies = {}
    api.mode_displacements_sc = {}
    api._mode_label_overrides = {}

    result = api.search_method_2(
        subgroup.index,
        distortion_type=["displacive"],
        candidates=[subgroup],
    )

    assert result.subgroup is subgroup
    assert api._selected_subgroup is subgroup
    assert set(api.mode_displacements_sc) == {mode.amplitude_key}
    assert np.allclose(
        sorted(api.mode_displacements_sc[mode.amplitude_key][:, 0]),
        [-1 / 3, -1 / 3, -1 / 3, 1 / 3, 1 / 3, 1 / 3],
    )

def _partial_microscopic_mapping_wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        result = subprocess.run(
            ["wsl.exe", "--status"],  # noqa: S607 - platform executable
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return result.returncode == 0

@pytest.mark.skipif(not _partial_microscopic_mapping_wsl_available(), reason="WSL/ISO binary unavailable")
def test_live_4310_secondary_e_modes_export_for_p_and_x_primary_routes():
    parent_path = WORKSPACE / "experiment_data" / "4310_tetra.cif"
    if not parent_path.is_file():
        pytest.skip("4310 parent input unavailable")

    api = _partial_microscopic_mapping_IsoDistort()
    api.set_distortion_scope({
        "displacive": ["*"],
        "occupational": [],
        "strain": ["*"],
        "magnetic": [],
        "rotational": [],
    })
    with contextlib.redirect_stdout(io.StringIO()):
        api.load_structure(parent_path)
        rows = api.search_method_1(distortion_types=["strain", "displacive"])
    candidates = [row.subgroup for row in rows]
    target_specs = (
        ("P2", "C1", 120, "(1/2,1/2,0)", 18),
        ("X1-", "C1", 50, "(0,1/2,0)", 16),
    )
    for irrep, opd, space_group, origin_raw, expected_count in target_specs:
        target = next(
            subgroup
            for subgroup in candidates
            if (
                subgroup.irrep_label == irrep
                and subgroup.opd_symbol == opd
                and subgroup.space_group_number == space_group
                and subgroup.origin_raw == origin_raw
            )
        )
        with contextlib.redirect_stdout(io.StringIO()):
            api.search_method_2(
                subgroup_idx=target.index,
                distortion_type=["displacive"],
                number_of_independent_modulations=0,
                candidates=candidates,
            )
        assert len(api.distortion_modes) == expected_count
        assert all(
            mode.mode_identity is not None
            and mode.mode_identity.status == "verified"
            and mode.mode_identity.source == "iso_microscopic"
            for mode in api.distortion_modes
        )
        spec = api._spec_for_subgroup(
            target,
            use_current_modes=True,
            use_generated_structure=False,
            folder_name=f"test_{irrep}_{opd}",
        )
        cif_text = render_cif(spec.structure, spec)
        assert "_iso_displacivemode_ID" in cif_text

# Original coverage: test_production_displacive_core.py.

def _p1_api_with_canonical_modes() -> tuple[_production_displacive_core_IsoDistort, SubgroupInfo]:
    parent = Structure(
        Lattice.from_parameters(4.1, 5.2, 6.3, 77.0, 83.0, 71.0),
        ["H", "He"],
        [[0.1, 0.2, 0.3], [0.37, 0.41, 0.53]],
    )
    sites = [
        {
            "orbit_id": "orbit-h",
            "wyckoff_letter": "a",
            "species": "H",
            "representative_index": 0,
            "equivalent_indices": [0],
            "multiplicity": 1,
            "display_label": "H1",
            "display_order": 0,
        },
        {
            "orbit_id": "orbit-he",
            "wyckoff_letter": "a",
            "species": "He",
            "representative_index": 1,
            "equivalent_indices": [1],
            "multiplicity": 1,
            "display_label": "He1",
            "display_order": 1,
        },
    ]
    subgroup = SubgroupInfo(
        index=0,
        parent_sg=1,
        space_group_number=1,
        space_group_symbol="P1",
        subgroup_index=1,
        size=1,
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a,b,c)",
        opd_vector=[1],
        basis_vectors=[[1, 1, 0], [0, 1, 0], [0, 0, 1]],
        origin=[Fraction(1, 3), Fraction(1, 5), Fraction(1, 7)],
    )
    exact_vectors = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    block = MicroscopicVectorBlock(
        "GM1+",
        "a",
        "A",
        (
            MicroscopicVectorRow(
                point_raw=("x", "y", "z"),
                displacements=exact_vectors,
            ),
        ),
        source_order=0,
        direction_symbol="VECTOR,A,B,C",
    )
    digest = hashlib.sha256(b"synthetic-p1-production-contract").hexdigest()
    modes = []
    for orbit_id in ("orbit-h", "orbit-he"):
        for component_index, (component, vector) in enumerate(
            zip(("a", "b", "c"), exact_vectors, strict=True)
        ):
            provenance = MicroscopicColumnProvenance.from_exact_column(
                query_digest=digest,
                query_order=0,
                query_irrep_label="GM1+",
                source_parent_sg=1,
                source_k_coordinates=("0", "0", "0"),
                source_subgroup_space_group_number=1,
                source_subgroup_basis=subgroup.basis_vectors,
                source_subgroup_origin=subgroup.origin,
                source_subgroup_irrep_label="GM1+",
                source_subgroup_opd_symbol="P1",
                source_subgroup_primary_direction_selector="VECTOR,A,B,C",
                block=block,
                column_index=component_index,
            )
            identity = ModeIdentity(
                parent_sg=1,
                global_irrep="GM1+",
                k_coordinates=("0", "0", "0"),
                wyckoff_letter="a",
                orbit_id=orbit_id,
                site_irrep="A",
                component_index=component_index,
                component_label=component,
                source="iso_microscopic",
                status="verified",
            )
            key = f"GM1+[0,0,0]__a@{orbit_id}__A({component})"
            modes.append(
                DistortionMode(
                    irrep_label="GM1+",
                    dimension=3,
                    mode_type="displacive",
                    wyckoff_site="a",
                    k_point_label="GM",
                    opd_symbol="P1",
                    opd_dir_raw="(a,b,c)",
                    bush_modes=[
                        BushMode(
                            irrep_label="GM1+",
                            opd_symbol="P1",
                            wyckoff_letter="a",
                            point=[0.0, 0.0, 0.0],
                            point_raw=["x", "y", "z"],
                            displacements=[
                                [float(value) for value in vector]
                            ],
                        )
                    ],
                    amplitude_key=key,
                    site_irrep="A",
                    k_coords_label="0,0,0",
                    opd_component=component,
                    wyckoff_orbit_id=orbit_id,
                    mode_identity=identity,
                    microscopic_provenance=provenance,
                )
            )

    api = object.__new__(_production_displacive_core_IsoDistort)
    api.structure = parent
    api.structure_path = None
    api.symmetry_info = {
        "space_group_number": 1,
        "space_group_symbol": "P1",
        "wyckoff_sites": sites,
    }
    api.cfg = get_config()
    api._dist_mapper = _production_displacive_core_DistortionMapper()
    api.distortion_modes = modes
    api.mode_displacements = {}
    api.mode_occupancies = {}
    api.mode_displacements_sc = {}
    api._mode_label_overrides = {}
    api.phase_path = None
    api._selected_subgroup = subgroup
    api.distortion_types = ["displacive"]
    return api, subgroup

def test_writer_identity_gate_checks_only_columns_emitted_for_candidate():
    api, subgroup = _p1_api_with_canonical_modes()
    api.mode_displacements = api._dist_mapper.map_modes_to_atoms(
        api.structure,
        api.symmetry_info["wyckoff_sites"],
        api.distortion_modes,
        subgroup_context=subgroup,
    )
    api._dist_engine = DistortionEngine(api._dist_mapper)
    lifted_by_key = api._lifted_mode_displacements(subgroup)
    source = api.distortion_modes[0]
    unresolved_key = "unused__unresolved(a)"
    unresolved = replace(
        source,
        amplitude_key=unresolved_key,
        mode_identity=ModeIdentity.unresolved(
            parent_sg=1,
            global_irrep="GM1+",
            k_coordinates=("0", "0", "0"),
            wyckoff_letter="a",
            orbit_id="orbit-h",
            reason="iso_microscopic_block_missing; parent_SG=1",
        ),
        microscopic_provenance=None,
    )
    api.distortion_modes = [*api.distortion_modes, unresolved]

    data = api._displacive_export_data_for_subgroup(subgroup, lifted_by_key)

    assert data is not None
    assert len(data.model.canonical_modes) == len(lifted_by_key)
    assert all(
        mode.mode_identity is not None and mode.mode_identity.status == "verified"
        for mode in data.model.canonical_modes
    )

    with pytest.raises(
        _production_displacive_core_UnresolvedModeIdentityError,
        match=r"iso_microscopic_block_missing.*parent_SG=1",
    ):
        api._displacive_export_data_for_subgroup(
            subgroup,
            {**lifted_by_key, unresolved_key: next(iter(lifted_by_key.values()))},
        )

def test_production_contract_preserves_nonzero_amplitude_across_origin_frame():
    api, subgroup = _p1_api_with_canonical_modes()
    embedding = embedding_from_identity(subgroup, parent_affine_group(api.structure))
    base = api._supercell_for_subgroup(subgroup)
    atom_ids = tuple(f"test-{index}" for index in range(len(base)))
    mapped = api._dist_mapper.map_microscopic_columns_to_supercell(
        api.structure,
        api.symmetry_info["wyckoff_sites"],
        api.distortion_modes,
        subgroup.basis_vectors,
        subgroup_context=subgroup,
        frame_id="test-frame",
        atom_ids=atom_ids,
        subgroup_operations=embedding.operations,
        subgroup_translation_lattice=[
            [float(value) for value in row] for row in embedding.lattice
        ],
    )
    lifted = {column.amplitude_key: column.displacements for column in mapped}
    selected_key = api.distortion_modes[0].amplitude_key
    final = Structure(
        base.lattice,
        base.species,
        (
            np.asarray(base.frac_coords, dtype=float)
            + 0.125 * lifted[selected_key]
        )
        % 1.0,
        coords_are_cartesian=False,
    )

    data = api._displacive_export_data_for_subgroup(
        subgroup,
        lifted,
        final_structure=final,
    )

    assert data is not None
    assert data.subgroup_identity.target_hall_number == embedding.hall_number
    assert len(data.model.canonical_modes) == 6
    assert np.count_nonzero(np.abs(np.asarray(data.model.amplitudes)) > 1.0e-10) == 1
    assert data.final_structure is not None
    assert np.allclose(
        np.asarray(data.final_structure.frac_coords, dtype=float),
        (
            np.asarray(final.frac_coords, dtype=float)
            - np.asarray([float(value) for value in subgroup.origin])
            @ np.linalg.inv(np.asarray(subgroup.basis_vectors, dtype=float))
        )
        % 1.0,
        atol=1.0e-12,
    )

def test_exact_fixed_space_production_contract_reaches_all_four_writers():
    api, subgroup = _p1_api_with_canonical_modes()
    subgroup.irrep_label = ""
    subgroup.opd_symbol = ""
    subgroup.opd_dir_raw = ""
    subgroup.k_point_label = ""
    subgroup.k_coordinates = []
    subgroup.opd_vector = []
    subgroup.method3_route_resolution = "exact_fixed_space"
    exact_modes = []
    for mode in api.distortion_modes:
        provenance = replace(
            mode.microscopic_provenance,
            source_subgroup_irrep_label=None,
            source_subgroup_opd_symbol=None,
            source_subgroup_primary_direction_selector=None,
            source_subgroup_context_kind="exact_fixed_space",
        )
        exact_modes.append(
            replace(
                mode,
                opd_symbol="",
                bush_modes=[
                    replace(bush, opd_symbol="") for bush in mode.bush_modes
                ],
                microscopic_provenance=provenance,
            )
        )
    api.distortion_modes = exact_modes

    embedding = embedding_from_identity(subgroup, parent_affine_group(api.structure))
    base = api._supercell_for_subgroup(subgroup)
    atom_ids = tuple(f"test-{index}" for index in range(len(base)))
    mapped = api._dist_mapper.map_microscopic_columns_to_supercell(
        api.structure,
        api.symmetry_info["wyckoff_sites"],
        api.distortion_modes,
        subgroup.basis_vectors,
        subgroup_context=subgroup,
        frame_id="test-exact-fixed-space-frame",
        atom_ids=atom_ids,
        subgroup_operations=embedding.operations,
        subgroup_translation_lattice=[
            [float(value) for value in row] for row in embedding.lattice
        ],
    )
    lifted = {column.amplitude_key: column.displacements for column in mapped}
    data = api._displacive_export_data_for_subgroup(subgroup, lifted)

    assert data is not None
    assert data.subgroup_identity.context_kind == "exact_fixed_space"
    assert data.subgroup_identity.irrep_label is None
    assert data.subgroup_identity.opd_symbol is None
    assert data.subgroup_identity.primary_direction_selector is None
    assert {
        mode.microscopic_provenance.source_global_irrep
        for mode in data.model.canonical_modes
    } == {"GM1+"}

    spec = SubgroupExportSpec(
        subgroup=subgroup,
        displacive_data=data,
        require_verified_displacive_data=True,
    )
    rendered = (
        render_cif(spec.cif_structure, spec),
        render_isoviz(spec),
        render_complete_modes(spec),
        render_topas(spec),
    )
    assert all(rendered)
    assert "Exact fixed-space coupled embedding" in rendered[0]
    assert "!irreplist \n  1 GM1+" in rendered[1]
    assert "exact fixed-space coupled embedding" in rendered[2]

    tampered_provenance = replace(
        exact_modes[0].microscopic_provenance,
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A,B,C",
        source_subgroup_context_kind="single_irrep",
    )
    tampered_modes = [
        replace(
            exact_modes[0],
            opd_symbol="P1",
            microscopic_provenance=tampered_provenance,
        ),
        *exact_modes[1:],
    ]
    with pytest.raises(
        ValueError,
        match=r"primary direction|query context|invalid exact OPD direction",
    ):
        api._dist_mapper.map_microscopic_columns_to_supercell(
            api.structure,
            api.symmetry_info["wyckoff_sites"],
            tampered_modes,
            subgroup.basis_vectors,
            subgroup_context=subgroup,
            frame_id="test-tampered-exact-fixed-space-frame",
            atom_ids=atom_ids,
            subgroup_operations=embedding.operations,
            subgroup_translation_lattice=[
                [float(value) for value in row] for row in embedding.lattice
            ],
        )

# Original coverage: test_real_mode_fallback.py.

def test_api_equal_atom_count_cannot_bypass_non_gamma_guard():
    api = object.__new__(_real_mode_fallback_IsoDistort)
    api.structure = _parent()
    api.mode_displacements = {"mode": {"displacements": np.asarray([[0.01, 0, 0]])}}
    api._selected_subgroup = None
    api.phase_path = SimpleNamespace(k_vector=[0, 0, Fraction(1, 4)])
    api._dist_engine = DistortionEngine()
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        api._method4_reference_modes(api.structure.copy())

# Original coverage: test_robustness.py.

_robustness_DATA_DIR = experiment_data_dir()

def test_negative_amplitude_at_api_level(tmp_path):
    """API 层负振幅：数学上合法（位移方向翻转），不崩溃。"""
    cif = _robustness_CIFS_DIR / "sg001.cif"
    if not cif.exists():
        pytest.skip("测试 CIF 不存在")
    iso = _robustness_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    subs = iso.list_subgroups_at("GM", "GM1")
    found = False
    for sg in subs:
        try:
            iso.search_method_2(subgroup_idx=sg.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
            continue
        if iso.mode_displacements:
            label = next(iter(iso.mode_displacements))
            d = iso.generate_distortion(irrep_label=label, amplitude=-0.1)
            assert len(d) == len(iso.structure)
            found = True
            break
    if not found:
        pytest.skip("该母相无可计算位移模式")

def test_method1_options_repeatable(tmp_path):
    """API 层 method1_options 两次调用：选项列表一致（可重复性）。"""
    cif = _robustness_CIFS_DIR / "sg139.cif"
    if not cif.exists():
        cif = _robustness_DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")
    iso = _robustness_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    o1 = iso.method1_options()
    o2 = iso.method1_options()
    assert [g["number"] for g in o1["space_groups"]] == \
        [g["number"] for g in o2["space_groups"]]
    assert [g["label"] for g in o1["conventional_lattices"]] == \
        [g["label"] for g in o2["conventional_lattices"]]

def _robustness_candidate(*, index: int, irrep: str) -> SubgroupInfo:
    return SubgroupInfo(
        index=index,
        space_group_number=12,
        subgroup_index=2,
        space_group_symbol="C2/m",
        k_point_label="GM",
        irrep_label=irrep,
        opd_symbol="P1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )

def test_method2_explicit_candidates_do_not_replace_other_result_pool():
    """Method tables may reuse indices; the explicit pool must select the target."""
    wrong = _robustness_candidate(index=0, irrep="GM1+")
    target = _robustness_candidate(index=0, irrep="GM2+")
    seen: dict = {}

    class Search:
        def search(self, parent_sg, candidates, query, wyckoff_letters,
                            **kwargs):
            seen["candidates"] = candidates
            return SimpleNamespace(
                subgroup=candidates[0], modes=[], metadata={},
            )

    class Mapper:
        @staticmethod
        def map_modes_to_atoms(structure, wyckoff_sites, modes):
            return {}

    iso = object.__new__(_robustness_IsoDistort)
    iso.structure = Structure(Lattice.cubic(1.0), ["H"], [[0.0, 0.0, 0.0]])
    iso.symmetry_info = {"space_group_number": 139, "wyckoff_sites": []}
    iso.subgroups = [wrong]
    iso._selected_subgroup = None
    iso._method2 = Search()
    iso._dist_mapper = Mapper()
    iso._iso = SimpleNamespace(list_k_points=lambda *_a, **_k: [])
    iso._smodes = None
    iso._letters_for_species = lambda species: []
    iso._union_scope_species = lambda types, **_kwargs: set()
    iso._resolve_k_vector = lambda *a, **k: [0.0, 0.0, 0.0]
    iso._compute_scoped_modes = lambda *args, **kwargs: []
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.number_of_independent_modulations = 0

    result = iso.search_method_2(0, candidates=[target])

    assert seen["candidates"] == [target]
    assert result.subgroup is target
    assert iso._selected_subgroup is target
    assert iso.subgroups == [wrong]

def test_export_current_modes_require_full_subgroup_identity():
    """Equal display indices from different Method tables must not share modes."""
    selected = _robustness_candidate(index=0, irrep="GM1+")
    other = _robustness_candidate(index=0, irrep="GM2+")
    iso = object.__new__(_robustness_IsoDistort)
    iso._selected_subgroup = selected
    iso.phase_path = SimpleNamespace(subgroup_index=0)
    iso.distortion_modes = []
    iso.mode_displacements = {"GM1+": {}}
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.distorted_structure = None
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    specs = iso._collect_export_specs([other], ["modes"], False)

    assert specs[0].use_current_modes is False
    assert specs[0].use_generated_structure is False

def test_export_current_cache_uses_exact_embedding_identity():
    """Embeddings that differ below display precision must not share cached modes."""

    selected = _robustness_candidate(index=0, irrep="GM1+")
    other = _robustness_candidate(index=0, irrep="GM1+")
    other.basis_vectors[0][0] = "1.00000000001"
    iso = object.__new__(_robustness_IsoDistort)
    iso._selected_subgroup = selected
    iso.phase_path = SimpleNamespace(subgroup_index=0)
    iso.distortion_modes = []
    iso.mode_displacements = {"cached": {}}
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.distorted_structure = None
    iso.number_of_independent_modulations = 0
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    specs = iso._collect_export_specs([other], ["modes"], False)

    assert specs[0].use_current_modes is False
    assert specs[0].use_generated_structure is False

def test_same_affine_embedding_different_routes_never_share_mode_or_final_cache():
    selected = _robustness_candidate(index=0, irrep="GM1+")
    selected.method3_embedding_id = "same-affine-subgroup"
    selected.k_coordinates = ["0", "0", "0"]
    selected.opd_dir_raw = "(a)"
    selected.opd_vector = [1]

    changed_irrep = deepcopy(selected)
    changed_irrep.irrep_label = "GM2+"
    changed_opd = deepcopy(selected)
    changed_opd.opd_symbol = "P2"
    changed_opd.opd_dir_raw = "(a,b)"
    changed_opd.opd_vector = [1, 1]
    changed_k = deepcopy(selected)
    changed_k.k_point_label = "X"
    changed_k.k_coordinates = ["1/2", "0", "0"]

    selected_key = _robustness_IsoDistort._method3_embedding_guard_key(selected)
    assert all(
        _robustness_IsoDistort._method3_embedding_guard_key(candidate) != selected_key
        for candidate in (changed_irrep, changed_opd, changed_k)
    )
    exact_fixed_space = deepcopy(selected)
    exact_fixed_space.method3_route_resolution = "exact_fixed_space"
    assert (
        _robustness_IsoDistort._method3_embedding_guard_key(exact_fixed_space)
        != selected_key
    )

    selected.basis_vectors = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    selected.origin = [0, 0, 0]
    selected_key = _robustness_IsoDistort._method3_embedding_guard_key(selected)
    parent_cell_shift = deepcopy(selected)
    parent_cell_shift.origin = [1, 0, 0]
    assert (
        _robustness_IsoDistort._method3_embedding_guard_key(parent_cell_shift)
        != selected_key
    )
    child_lattice_shift = deepcopy(selected)
    child_lattice_shift.origin = [2, 0, 0]
    assert (
        _robustness_IsoDistort._method3_embedding_guard_key(child_lattice_shift)
        == selected_key
    )

    iso = object.__new__(_robustness_IsoDistort)
    iso._selected_subgroup = selected
    iso.phase_path = SimpleNamespace(subgroup_index=0)
    iso.distortion_modes = []
    iso.mode_displacements = {"cached": {}}
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.distorted_structure = object()
    iso.number_of_independent_modulations = 0
    iso.distortion_types = ["displacive"]
    iso.distortion_scope = {}
    iso._scope_species = lambda _kind, **_kwargs: {"H"}
    selected_context = iso._mode_context_key(
        selected,
        0,
        ["displacive"],
        {},
    )
    iso._mode_cache_key = selected_context
    iso._generated_structure_cache_key = selected_context
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )

    specs = iso._collect_export_specs([changed_irrep], ["modes"], False)

    assert specs[0].use_current_modes is False
    assert specs[0].use_generated_structure is False

def test_method2_rejects_duplicate_local_candidate_indices():
    iso = object.__new__(_robustness_IsoDistort)
    iso.structure = object()
    iso.subgroups = []
    first = _robustness_candidate(index=0, irrep="GM1+")
    second = _robustness_candidate(index=0, irrep="GM2+")

    with pytest.raises(ValueError, match="ambiguous"):
        iso.search_method_2(0, candidates=[first, second])

def test_compute_scoped_modes_does_not_mutate_source_mode():
    raw = DistortionMode(
        irrep_label="GM1+",
        mode_type="magnetic",
        wyckoff_site="a",
        wyckoff_orbit_id="orbit-a",
        bush_modes=[
            BushMode(
                irrep_label="GM1+",
                opd_symbol="P1",
                wyckoff_letter="a",
                point=[0.0, 0.0, 0.0],
            )
        ],
    )
    iso = object.__new__(_robustness_IsoDistort)
    iso.mode_occupancies = {}
    iso._scope_species = lambda _kind, **_kwargs: {"Fe"}
    iso._letters_for_species = lambda _species: ["a"]
    iso._orbit_ids_for_species = lambda _species: ["orbit-a"]

    scoped = iso._compute_scoped_modes(
        ["displacive"],
        raw_modes=[raw],
    )

    assert raw.mode_type == "magnetic"
    assert scoped[0] is not raw
    assert scoped[0].mode_type == "displacive"

@pytest.mark.parametrize(
    ("export_method", "legacy"),
    ((2, True), (3, False)),
)
def test_export_method_rejects_conflicting_legacy_folder_switch(
    export_method,
    legacy,
):
    with pytest.raises(ValueError, match="conflicts"):
        _robustness_IsoDistort._export_folder_method(export_method, legacy)

def test_export_method_legacy_switch_remains_compatible_when_new_argument_omitted():
    assert _robustness_IsoDistort._export_folder_method(None, None) == 2
    assert _robustness_IsoDistort._export_folder_method(None, True) == 1
    assert _robustness_IsoDistort._export_folder_method(None, False) == 2
    assert _robustness_IsoDistort._export_folder_method(1, True) == 1

def test_export_recomputes_every_candidate_for_requested_nmod_and_restores_state():
    """A ZIP must not mix cached modes from one nmod with nmod=0 fallbacks."""
    selected = _robustness_candidate(index=0, irrep="LD1")
    other = _robustness_candidate(index=1, irrep="LD2")
    iso = object.__new__(_robustness_IsoDistort)
    iso._selected_subgroup = selected
    iso.phase_path = SimpleNamespace(subgroup_index=0)
    original_mode = DistortionMode(
        irrep_label="LD1",
        mode_type="magnetic",
        wyckoff_site="a",
        wyckoff_orbit_id="orbit-a",
        bush_modes=[
            BushMode(
                irrep_label="LD1",
                opd_symbol="P1",
                wyckoff_letter="a",
                point=[0.0, 0.0, 0.0],
            )
        ],
    )
    iso.distortion_modes = [original_mode]
    iso.mode_displacements = {"cached-nmod-1": {}}
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.distorted_structure = None
    iso.number_of_independent_modulations = 1
    iso._scope_species = lambda _kind, **_kwargs: {"Fe"}
    iso._letters_for_species = lambda _species: ["a"]
    iso._orbit_ids_for_species = lambda _species: ["orbit-a"]
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SimpleNamespace(
        subgroup=subgroup, **kwargs
    )
    calls: list[tuple[int, int]] = []

    def search_method_2(
        subgroup_idx,
        *,
        distortion_type,
        number_of_independent_modulations,
        candidates,
        distortion_scope,
    ):
        calls.append((subgroup_idx, number_of_independent_modulations))
        assert distortion_type == ["displacive", "strain"]
        assert distortion_scope == {}
        iso._selected_subgroup = candidates[0]
        iso.number_of_independent_modulations = number_of_independent_modulations
        iso.mode_displacements = {f"fresh-{subgroup_idx}": {}}
        iso.distortion_modes = iso._compute_scoped_modes(
            ["displacive"],
            raw_modes=[original_mode],
        )

    iso.search_method_2 = search_method_2

    specs = iso._collect_export_specs(
        [selected, other],
        ["modes"],
        True,
        number_of_independent_modulations=2,
    )

    assert calls == [(0, 2), (1, 2)]
    assert [spec.use_current_modes for spec in specs] == [True, True]
    assert [spec.use_generated_structure for spec in specs] == [False, False]
    assert iso._selected_subgroup is selected
    assert iso.distortion_modes == [original_mode]
    assert iso.distortion_modes[0].mode_type == "magnetic"
    assert iso.mode_displacements == {"cached-nmod-1": {}}
    assert iso.number_of_independent_modulations == 1

# Original coverage: test_search_methods.py.

def _bare_mode_label_api(mode: _search_methods_DistortionMode) -> _search_methods_IsoDistort:
    api = object.__new__(_search_methods_IsoDistort)
    api.structure = Structure(Lattice.cubic(4.0), ["Na"], [[0.0, 0.0, 0.0]])
    api.symmetry_info = {
        "space_group_number": 221,
        "wyckoff_sites": [
            {
                "wyckoff_letter": "a",
                "orbit_id": "orbit-a",
                "species": "Na",
                "display_label": "Na1",
                "representative_index": 0,
            }
        ],
    }
    api.mode_displacements = {
        mode.amplitude_key: {
            "mode": mode,
            "wyckoff_letter": "a",
            "displacements": np.array([[1.0, 0.0, 0.0]]),
        }
    }
    api.distortion_modes = [mode]
    api.mode_occupancies = {}
    api._mode_label_overrides = {}
    api.phase_path = None
    api._selected_subgroup = None
    return api

def test_mode_labels_use_only_verified_iso_microscopic_identity() -> None:
    identity = ModeIdentity(
        parent_sg=221,
        global_irrep="GM4-",
        k_coordinates=("0", "0", "0"),
        wyckoff_letter="a",
        orbit_id="orbit-a",
        site_irrep="T1u",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="verified",
    )
    mode = _search_methods_DistortionMode(
        irrep_label="GM4-",
        wyckoff_site="a",
        wyckoff_orbit_id="orbit-a",
        amplitude_key="verified-mode",
        mode_identity=identity,
    )

    labels = _bare_mode_label_api(mode)._mode_labels_now()

    assert labels == {
        "verified-mode": "Pm-3m[0,0,0]GM4-(a)[Na1:a:dsp]T1u(a)"
    }

def test_mode_labels_use_each_microscopic_irreps_complete_direction() -> None:
    identity = ModeIdentity(
        parent_sg=139,
        global_irrep="X1+",
        k_coordinates=("1/2", "1/2", "0"),
        wyckoff_letter="a",
        orbit_id="orbit-a",
        site_irrep="A1",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="verified",
    )
    mode = _search_methods_DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        wyckoff_orbit_id="orbit-a",
        k_point_label="X",
        opd_symbol="C1",
        opd_dir_raw="(a;a)",
        amplitude_key="secondary-x-mode",
        mode_identity=identity,
    )
    api = _bare_mode_label_api(mode)
    api.symmetry_info["space_group_number"] = 139
    api.phase_path = SimpleNamespace(opd_dir_raw="(a;b;c;d)")

    labels = api._mode_labels_now()

    assert labels == {
        "secondary-x-mode": (
            "I4/mmm[1/2,1/2,0]X1+(a;a)[Na1:a:dsp]A1(a)"
        )
    }

def test_mode_labels_fail_closed_for_unresolved_site_irrep() -> None:
    mode = _search_methods_DistortionMode(
        irrep_label="GM4-",
        wyckoff_site="a",
        wyckoff_orbit_id="orbit-a",
        amplitude_key="unresolved-mode",
        mode_identity=ModeIdentity.unresolved(
            parent_sg=221,
            global_irrep="GM4-",
            k_coordinates=("0", "0", "0"),
            wyckoff_letter="a",
            orbit_id="orbit-a",
            reason="test evidence is absent",
        ),
    )

    api = _bare_mode_label_api(mode)
    api._mode_label_overrides = {
        "unresolved-mode": "Pm-3m[0,0,0]GM4-(a)[Na1:a:dsp]A1(a)"
    }
    with pytest.raises(ValueError, match="test evidence is absent"):
        api._mode_labels_now()

def test_method3_api_propagates_physical_tolerance_model(monkeypatch) -> None:
    # This test exercises only the API/configuration boundary; no external
    # ISO executable or WSL session is part of the tolerance model.
    monkeypatch.setattr(core_api_module, "IsoWrapper", SimpleNamespace)
    monkeypatch.setattr(core_api_module, "SmodesWrapper", SimpleNamespace)
    api = _search_methods_IsoDistort()
    api.structure = _i4mmm_parent()
    api.symmetry_info = {"space_group_number": 139}
    captured = {}

    def _capture(parent_sg, query):
        captured["parent_sg"] = parent_sg
        captured["query"] = query
        return []

    monkeypatch.setattr(api._method3, "search", _capture)
    assert api.search_method_3(space_group_type=139) == []

    cfg = api.cfg
    query = captured["query"]
    assert captured["parent_sg"] == 139
    assert query.symmetry_tolerance == cfg.symmetry_cartesian_tolerance_angstrom
    assert (
        query.symmetry_angle_tolerance_degrees
        == cfg.symmetry_angle_tolerance_degrees
    )
    assert (
        query.affine_exact_tolerance
        == cfg.affine_exact_cartesian_tolerance_angstrom
    )
    assert (
        query.fractional_coordinate_tolerance
        == cfg.fractional_coordinate_tolerance
    )

def test_public_method_2_api_rejects_superposed_ir_request_before_search():
    api = object.__new__(_search_methods_IsoDistort)
    api.structure = object()

    with pytest.raises(NotImplementedError, match="cannot be combined as a union"):
        api.search_method_2(
            subgroup_idx=0,
            number_of_superposed_irreps=2,
            candidates=[],
        )

def _method3_filter_api(smodes) -> _search_methods_IsoDistort:
    api = object.__new__(_search_methods_IsoDistort)
    api.structure = [SimpleNamespace(species_string="Eu")]
    api.symmetry_info = {"space_group_number": 139, "wyckoff_sites": []}
    api.distortion_types = ["strain", "displacive"]
    api.distortion_scope = {
        "displacive": ["Eu"],
        "rotational": [],
        "occupational": ["Eu"],
        "magnetic": [],
        "strain": [],
    }
    api._smodes = smodes
    return api

def test_method_3_filters_routes_with_parameter_specific_activity() -> None:
    class _RecordingSmodes:
        def __init__(self) -> None:
            self.calls = []

        def active_irreps(
            self, structure, parent_sg, wyckoff_sites, k_label,
            k_parameters=None, species_filter=None,
        ):
            self.calls.append((k_label, k_parameters, species_filter))
            return {"LD1"}

    smodes = _RecordingSmodes()
    api = _method3_filter_api(smodes)
    # Route identity stores the website's physical g=1/6, while the bundled
    # backend represents this LD parameter on its internal half scale.
    api._resolve_iso_kparams = lambda label, params: ["1/12"]
    inactive = _method3_route("LD2", index=7)
    active = _method3_route("LD1", index=8)
    missing = _method3_route("LD3", index=9, origin=[0.5, 0, 0])
    items = [
        Method3ResultItem(
            subgroup=inactive,
            point_group="4mm",
            basis=inactive.basis_vectors,
            routes=[inactive, active],
        ),
        Method3ResultItem(
            subgroup=missing,
            point_group="4mm",
            basis=missing.basis_vectors,
            routes=[missing],
        ),
        Method3ResultItem(
            subgroup=_method3_route("LD1", index=10, origin=[0, 0.5, 0]),
            point_group="4mm",
            basis=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
            routes=[],
        ),
    ]

    filtered = api._filter_method3_routes_by_types(items, ["displacive"])

    assert len(filtered) == 1
    assert filtered[0].routes == [active]
    assert filtered[0].subgroup is active
    assert filtered[0].subgroup.index == 0
    assert smodes.calls == [("LD", ["1/12"], {"Eu"})]

def test_method_2_rejects_affine_only_method3_diagnostic_row() -> None:
    unresolved = _method3_route("", index=4)
    api = _method3_filter_api(SimpleNamespace())
    api.subgroups = [unresolved]
    api._unresolved_method3_embedding_keys = {
        api._method3_embedding_guard_key(unresolved)
    }

    # The guard is a scientific path identity, not a Python object address:
    # an equivalent result reconstructed by another API/UI layer remains
    # non-selectable even without the transient marker attribute.
    equivalent = deepcopy(unresolved)
    assert api._method3_embedding_guard_key(equivalent) == (
        api._method3_embedding_guard_key(unresolved)
    )

    close_but_distinct = deepcopy(unresolved)
    close_but_distinct.origin = [1e-11, 0.0, 0.0]
    assert api._method3_embedding_guard_key(close_but_distinct) != (
        api._method3_embedding_guard_key(unresolved)
    )

    with pytest.raises(RuntimeError, match="no resolved single-IR or coupled-IR"):
        api.search_method_2(
            subgroup_idx=4,
            candidates=[equivalent],
        )

@pytest.mark.parametrize("kind", ["magnetic", "occupational", "rotational"])
def test_method_3_route_filter_rejects_unproved_physical_models(kind) -> None:
    class _NoSmodes:
        def active_irreps(self, *args, **kwargs):
            raise AssertionError("displacement activity cannot certify this model")

    api = _method3_filter_api(_NoSmodes())
    route = _method3_route("mLD1" if kind == "magnetic" else "LD1")
    item = Method3ResultItem(route, "4mm", route.basis_vectors, [route])
    with pytest.raises(NotImplementedError, match=kind):
        api._filter_method3_routes_by_types([item], [kind])

def test_method_1_strain_filter_uses_symmetric_tensor_fixed_space() -> None:
    api = _search_methods_IsoDistort()
    api.structure = _i4mmm_parent()
    api.symmetry_info = {
        "space_group_number": 139,
        "wyckoff_sites": [],
    }

    def _item(
        irrep: str,
        child_sg: int,
        child_symbol: str,
        basis: list[list[int]],
    ) -> Method1ResultItem:
        route = _search_methods_SubgroupInfo(
            index=0,
            parent_sg=139,
            k_point_label="GM",
            irrep_label=irrep,
            size=1,
            subgroup_index=1,
            space_group_number=child_sg,
            space_group_symbol=child_symbol,
            basis_vectors=basis,
            origin=[0, 0, 0],
        )
        return Method1ResultItem(
            subgroup=route,
            crystal_system="tetragonal",
            is_maximal=False,
        )

    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    gm1 = _item("GM1+", 139, "I4/mmm", identity)
    gm4 = _item(
        "GM4+",
        69,
        "Fmmm",
        [[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
    )
    gm3 = _item("GM3+", 87, "I4/m", identity)
    polar = _item("GM5-", 25, "Pmm2", identity)

    kept = api._filter_method1_by_types(
        [gm1, gm4, gm3, polar], ["strain"]
    )

    assert [item.subgroup.irrep_label for item in kept] == ["GM1+", "GM4+"]

def test_displacive_species_scope_keeps_only_selected_shared_letter_orbits(
    monkeypatch,
) -> None:
    target = _method3_route("P2")
    target.k_parameters = []
    rooted = _search_methods_DistortionMode(
        irrep_label="P2",
        wyckoff_site="e",
        amplitude_key="P2__e__C1__a",
        bush_modes=[_search_methods_BushMode(
            irrep_label="P2",
            opd_symbol="C1",
            wyckoff_letter="e",
            point=[0, 0, 0],
            displacements=[[1, 0, 0]],
        )],
    )
    secondary = _search_methods_DistortionMode(
        irrep_label="GM1+",
        wyckoff_site="g",
        wyckoff_orbit_id="orbit-h-g",
        amplitude_key="GM1+__orbit-h-g",
        bush_modes=[_search_methods_BushMode(
            irrep_label="GM1+",
            opd_symbol="C1",
            wyckoff_letter="g",
            point=[0, 0, 0],
            displacements=[[0, 1, 0]],
        )],
    )
    sites = [
        {"orbit_id": "orbit-h-e", "wyckoff_letter": "e", "species": "H"},
        {"orbit_id": "orbit-he-e", "wyckoff_letter": "e", "species": "He"},
        {"orbit_id": "orbit-h-g", "wyckoff_letter": "g", "species": "H"},
        {"orbit_id": "orbit-he-g", "wyckoff_letter": "g", "species": "He"},
    ]

    class _Iso:
        @staticmethod
        def list_k_points(parent_sg):
            return []

        @staticmethod
        def calc_distortion_modes(parent_sg, selected, wyckoff_letters):
            assert wyckoff_letters == ["e", "g"]
            return [rooted]

    def _fixed_space(parent, symmetry_info, selected, letters, smodes, **kwargs):
        assert letters == ["g"]
        assert kwargs["wyckoff_orbit_ids"] == ["orbit-h-g"]
        return superspace_module.ParametricModeResult(
            modes=[secondary],
            supercell_displacements={secondary.amplitude_key: np.zeros((2, 3))},
            labels={secondary.amplitude_key: "secondary"},
            nmod=0,
        )

    monkeypatch.setattr(superspace_module, "compute_parametric_modes", _fixed_space)
    api = object.__new__(_search_methods_IsoDistort)
    api.structure = [
        SimpleNamespace(species_string="H"),
        SimpleNamespace(species_string="He"),
    ]
    api.symmetry_info = {"space_group_number": 139, "wyckoff_sites": sites}
    api.distortion_scope = {"displacive": ["H"]}
    api.number_of_independent_modulations = 0
    api._iso = _Iso()
    api._smodes = object()

    selected_species = api._scope_species("displacive")
    result = Method2Solver(api._iso).search(
        139, [target], Method2Query(subgroup_idx=target.index),
        wyckoff_letters=api._letters_for_species(selected_species),
        wyckoff_orbit_ids=api._orbit_ids_for_species(selected_species),
        structure=api.structure,
        wyckoff_sites=sites,
        smodes=api._smodes,
        symmetry_info=api.symmetry_info,
    )
    modes = api._compute_scoped_modes(["displacive"], raw_modes=result.modes)

    assert [mode.wyckoff_orbit_id for mode in modes] == [
        "orbit-h-e", "orbit-h-g",
    ]
    assert "orbit-h-e" in modes[0].amplitude_key
    assert all("orbit-he" not in mode.amplitude_key for mode in modes)

def test_method_4_api_uses_selected_child_cell_and_precomputed_modes():
    api = object.__new__(_search_methods_IsoDistort)
    api.structure = Structure(
        Lattice.orthorhombic(4.0, 5.0, 6.0),
        ["Na"],
        [[0.25, 0.0, 0.0]],
    )
    api._selected_subgroup = SimpleNamespace(
        basis_vectors=[[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    )
    api.mode_displacements = {
        "parent-mode": {"displacements": np.array([[1.0, 0.0, 0.0]])},
    }
    child_mode = np.array([[0.0, 0.25, 0.0]])
    api.mode_displacements_sc = {"child-mode": child_mode}
    api._dist_engine = SimpleNamespace()

    expected_child = api._supercell_for_subgroup(api._selected_subgroup)
    reference, modes = api._method4_reference_modes(expected_child.copy())

    assert np.allclose(reference.lattice.matrix, expected_child.lattice.matrix)
    assert np.allclose(reference.frac_coords, expected_child.frac_coords)
    assert set(modes) == {"child-mode"}
    assert np.array_equal(modes["child-mode"], child_mode)

# Original coverage: test_session_result_contract.py.

def test_action_snapshot_preserves_internal_aliases_and_new_state(monkeypatch):
    class RuntimeStub:
        def clone_runtime_from(self, source):
            self.source = source
        def close(self):
            pass
    monkeypatch.setattr(core_api, "IsoWrapper", RuntimeStub)
    monkeypatch.setattr(core_api, "SmodesWrapper", RuntimeStub)
    session = _session_result_contract_IsoDistort()
    selected = SubgroupInfo(index=0, space_group_number=1)
    session.subgroups = [selected]
    session._selected_subgroup = selected
    session.additional_scientific_state = {"candidate": selected, "values": [1]}
    snapshot = session.snapshot_for_action()
    assert snapshot.structure is None
    assert snapshot._selected_subgroup is snapshot.subgroups[0]
    assert snapshot.additional_scientific_state["candidate"] is snapshot.subgroups[0]
    snapshot.additional_scientific_state["values"].append(2)
    snapshot._selected_subgroup.index = 3
    assert session.subgroups[0].index == 0
    assert session.additional_scientific_state["values"] == [1]
    assert snapshot._iso is not session._iso

def _mode_context_review_session():
    """Pure state fixture; no scientific runtime or structure solver."""
    session = object.__new__(core_api.IsoDistort)

    class SyntheticStructure(list):
        lattice = SimpleNamespace(matrix=np.eye(3), volume=1)

    session.structure = SyntheticStructure([
        SimpleNamespace(species_string="H"), SimpleNamespace(species_string="He"),
    ])
    session.symmetry_info = {"space_group_number": 1}
    session.distortion_types = ["displacive"]
    session.distortion_scope = {"displacive": ["H"]}
    session.number_of_independent_modulations = 0
    session.mode_displacements = {}
    session.mode_displacements_sc = {"H-only": np.array([[1., 0., 0.], [0., 0., 0.]])}
    session._selected_subgroup = SubgroupInfo(index=0, space_group_number=1, parent_sg=1)
    session._supercell_for_subgroup = lambda _sg: session.structure
    session._mode_cache_key = session._mode_context_key(
        session._selected_subgroup, 0, session.distortion_types,
    )
    return session


@pytest.mark.parametrize("scope,expected", [
    (None, {"H", "He"}), ({}, {"H", "He"}),
    ({"occupational": []}, {"H", "He"}),
    ({"displacive": "*"}, {"H", "He"}),
    ({"displacive": "all"}, {"H", "He"}),
    ({"displacive": ["*"]}, {"H", "He"}),
    ({"displacive": ["He"]}, {"He"}),
    ({"displacive": ["He", "absent-species"]}, {"He"}),
    ({"displacive": None}, {"H", "He"}),
    ({"displacive": ("He",)}, {"He"}),
    (MappingProxyType({"displacive": ["He"]}), {"He"}),
    ({"displacive": " ALL "}, {"H", "He"}),
])
def test_scope_missing_all_and_selected_species_remain_distinct(scope, expected):
    session = _mode_context_review_session()
    session.set_distortion_scope(scope)
    session.set_distortion_types(["displacive"])
    assert session._scope_species("displacive") == expected


@pytest.mark.parametrize("scope", [False, 0, 0.0, "", [], (), set(), "all", True, 1, ["He"]])
def test_scope_rejects_nonmapping_without_resetting_existing_selection(scope):
    session = _mode_context_review_session()
    previous = session.distortion_scope
    with pytest.raises(ValueError, match="scope"):
        session.set_distortion_scope(scope)
    assert session.distortion_scope is previous
    assert session.distortion_scope == {"displacive": ["H"]}


@pytest.mark.parametrize("invalid", [
    "He", True, False, 0, 1.2, {}, {"H": 1}, "", set(),
    ["H", 1], ["H", ""], ["H", None], ["*", False],
])
def test_scope_rejects_invalid_explicit_value_without_partial_commit(invalid):
    session = _mode_context_review_session()
    previous = session.distortion_scope
    # The valid first field must not become visible when the later field fails.
    with pytest.raises(ValueError, match="scope"):
        session.set_distortion_scope({"displacive": ["He"], "magnetic": invalid})
    assert session.distortion_scope is previous
    assert session.distortion_scope == {"displacive": ["H"]}


def test_scope_rejects_unknown_type_without_replacing_selection():
    session = _mode_context_review_session()
    previous = session.distortion_scope
    with pytest.raises(ValueError, match="scope"):
        session.set_distortion_scope({"displacive": ["He"], "not-a-distortion-type": ["H"]})
    assert session.distortion_scope is previous


@pytest.mark.parametrize("selected", [["He"], ["absent-species"]])
def test_scope_configuration_defers_parent_species_intersection_until_loaded(selected):
    session = _mode_context_review_session()
    parent = session.structure
    session.structure = None
    session.set_distortion_scope({"displacive": selected})
    session.set_distortion_types(["displacive"])
    assert session.distortion_scope == {"displacive": selected}
    assert session.distortion_types == ["displacive"]
    # No loader or scientific identification is run; only parent availability changes.
    session.structure = parent
    if selected == ["He"]:
        assert session._scope_species("displacive") == {"He"}
    else:
        with pytest.raises(ValueError, match="scope.*empty"):
            session._scope_species("displacive")


def test_scope_configuration_rejects_explicit_empty_even_before_parent_load():
    session = _mode_context_review_session()
    session.structure = None
    session.set_distortion_scope({"displacive": []})
    with pytest.raises(ValueError, match="scope.*empty"):
        session.set_distortion_types(["displacive"])


@pytest.mark.parametrize("selected_scope", [[], ["absent-species"]])
def test_explicit_empty_scope_never_expands_to_all_species(selected_scope):
    session = _mode_context_review_session()
    session.set_distortion_scope({"displacive": selected_scope, "occupational": [], "strain": []})
    assert session.distortion_scope["displacive"] == selected_scope
    with pytest.raises(ValueError, match="empty.*scope|scope.*empty"):
        session._scope_species("displacive")
    session.set_distortion_types(["strain"])
    with pytest.raises(ValueError, match="empty.*scope|scope.*empty"):
        session.set_distortion_types(["displacive"])
    assert session.distortion_types == ["strain"]
    assert session.distortion_scope["occupational"] == []
    assert session._scope_species("strain") == {"H", "He"}


@pytest.mark.parametrize("changed", ["scope", "nmod", "types", "candidate", "missing-key", "missing-selection"])
def test_method4_rejects_stale_mode_context_before_reading_daughter(monkeypatch, changed):
    session = _mode_context_review_session()
    if changed == "scope":
        session.set_distortion_scope({"displacive": ["He"]})
    elif changed == "nmod":
        session.number_of_independent_modulations = 1
    elif changed == "types":
        session.set_distortion_types(["strain"])
    elif changed == "candidate":
        session._selected_subgroup = replace(session._selected_subgroup, irrep_label="different")
    elif changed == "missing-key":
        session._mode_cache_key = None
    else:
        session._selected_subgroup = None
    monkeypatch.setattr(core_api, "read_cif", lambda _path: pytest.fail("stale basis read daughter"))
    with pytest.raises(RuntimeError, match="context changed|No committed"):
        session.search_method_4(Path("not-read.cif"))


@pytest.mark.parametrize("expected", [True, 1.5, "1", -1, 4, 1])
def test_method4_checks_expected_nmod_without_changing_basis(monkeypatch, expected):
    session = _mode_context_review_session()
    previous_key = session._mode_cache_key
    monkeypatch.setattr(core_api, "read_cif", lambda _path: pytest.fail("invalid nmod read daughter"))
    error_type = NotImplementedError if expected == 4 else (ValueError, TypeError)
    with pytest.raises(error_type, match="modulation|nmod"):
        session.search_method_4(Path("not-read.cif"), number_of_independent_modulations=expected)
    assert session.number_of_independent_modulations == 0
    assert session._mode_cache_key == previous_key


@pytest.mark.parametrize("committed_nmod", [0, 1, 3])
@pytest.mark.parametrize("with_expectation", [False, True])
def test_method4_accepts_repeated_daughters_with_same_committed_basis(monkeypatch, committed_nmod, with_expectation):
    session = _mode_context_review_session()
    session.number_of_independent_modulations = committed_nmod
    session._mode_cache_key = session._mode_context_key(
        session._selected_subgroup, committed_nmod, session.distortion_types,
    )
    calls = []

    class DecomposerStub:
        def method_4_decompose(self, parent, daughter, modes, query):
            calls.append((parent, daughter, modes))
            return SimpleNamespace(amplitudes={}, rms_residual=0)

    monkeypatch.setattr(core_api, "ModeDecomposer", DecomposerStub)
    monkeypatch.setattr(core_api, "read_cif", lambda _path: session.structure)
    monkeypatch.setattr(core_api, "read_cif_space_group_number", lambda _path: 1)
    previous_key = session._mode_cache_key
    kwargs = {"number_of_independent_modulations": committed_nmod} if with_expectation else {}
    for filename in ("first.cif", "second.cif"):
        session.search_method_4(Path(filename), **kwargs)
    assert len(calls) == 2
    assert session._mode_cache_key == previous_key
    assert session.number_of_independent_modulations == committed_nmod
    for _basis_parent, _daughter, modes in calls:
        np.testing.assert_array_equal(modes["H-only"], [[1, 0, 0], [0, 0, 0]])


def _calculated_session():
    # Independent P1 full microscopic basis on an oblique lattice. These are
    # analytically constructed unit displacement columns, not website outputs.
    session, selected = _p1_api_with_canonical_modes()
    session.mode_displacements_sc = {}
    for index, mode in enumerate(session.distortion_modes):
        column = np.zeros((len(session.structure), 3))
        column[index // 3, index % 3] = -1 if index % 2 else 1
        session.mode_displacements_sc[mode.amplitude_key] = column
    session.distortion_scope = {}
    session.number_of_independent_modulations = 0
    session._mode_cache_key = session._mode_context_key(selected, 0, session.distortion_types)
    return session

def _selection_session(monkeypatch):
    """Exercise both public adapters with an independent analytic P1 basis.

    Only the backend mode calculation is replaced. Candidate guards, species
    selection, session commit, and the immutable snapshot are real code.
    """
    session = _calculated_session()
    session.subgroups = [session._selected_subgroup]
    session.distortion_types = ["displacive", "strain"]
    session._iso = SimpleNamespace(list_k_points=lambda _parent: [])
    session._smodes = object()
    session._method2 = Method2Solver(session._iso)
    source_modes = session.distortion_modes
    source_columns = session.mode_displacements_sc
    calls = []

    def calculate(_parent, _info, target, letters, *_runtimes, **kwargs):
        calls.append((target.index, tuple(letters), kwargs.get("nmod", 0),
                      tuple(kwargs.get("wyckoff_orbit_ids") or ())))
        allowed = set(kwargs.get("wyckoff_orbit_ids") or ())
        modes = [mode for mode in source_modes
                 if not allowed or mode.wyckoff_orbit_id in allowed]
        return superspace.ParametricModeResult(
            modes=modes,
            supercell_displacements={mode.amplitude_key: source_columns[mode.amplitude_key].copy()
                                     for mode in modes},
            labels={mode.amplitude_key: mode.amplitude_key for mode in modes},
            nmod=kwargs.get("nmod", 0),
        )

    monkeypatch.setattr(superspace, "compute_parametric_modes", calculate)
    monkeypatch.setattr(superspace, "compute_special_modes_with_rootless_supplement", calculate)
    return session, calls

@pytest.mark.parametrize("entry", ["select_path", "search_method_2"])
@pytest.mark.parametrize("invalid", ["ambiguous", "diagnostic", "registered_diagnostic"])
def test_public_selection_rejects_unsafe_candidates_before_calculation(monkeypatch, entry, invalid):
    session, calls = _selection_session(monkeypatch)
    if invalid == "ambiguous":
        session.subgroups.append(replace(session.subgroups[0], irrep_label="GM2+"))
        expected, message = ValueError, "ambiguous"
    elif invalid == "diagnostic":
        session.subgroups[0].method3_route_resolution = "affine_only_unresolved_coupled_route"
        expected, message = RuntimeError, "no resolved single-IR or coupled-IR"
    else:
        session._unresolved_method3_embedding_keys = {
            session._method3_embedding_guard_key(session.subgroups[0])}
        expected, message = RuntimeError, "no resolved single-IR or coupled-IR"
    old_phase = session.phase_path
    old_modes = session.distortion_modes
    with pytest.raises(expected, match=message):
        getattr(session, entry)(0)
    assert calls == []
    assert session.phase_path is old_phase
    assert session.distortion_modes is old_modes

def test_select_path_commits_explicit_types_before_snapshot(monkeypatch):
    session, calls = _selection_session(monkeypatch)
    path = session.select_path(0, distortion_type=["displacive"])
    snapshot = session.displacement_snapshot()
    assert path is session.phase_path
    assert session.distortion_types == ["displacive"]
    assert len(snapshot.columns) == 6
    assert len(calls) == 1

@pytest.mark.parametrize("nmod", [0, 1, 3])
@pytest.mark.parametrize("species", [["*"], ["H"], ["He"]])
@pytest.mark.parametrize("types", [None, ["displacive"]])
def test_public_selection_entries_preserve_current_nmod_and_species_scope(monkeypatch, nmod, species, types):
    snapshots = []
    paths = []
    for entry in ("select_path", "search_method_2"):
        session, calls = _selection_session(monkeypatch)
        session.number_of_independent_modulations = nmod
        session.set_distortion_scope({"displacive": species})
        kwargs = {} if entry == "select_path" else {"number_of_independent_modulations": nmod}
        returned = getattr(session, entry)(0, distortion_type=types, **kwargs)
        snapshots.append(session.displacement_snapshot())
        paths.append(asdict(session.phase_path))
        assert session.number_of_independent_modulations == nmod
        assert session.distortion_scope["displacive"] == species
        assert session.distortion_types == (types or ["displacive", "strain"])
        assert len(calls) == 1
        if entry == "select_path":
            assert returned is session.phase_path
        else:
            assert returned.metadata["number_of_independent_modulations"] == nmod
    assert snapshots[0] == snapshots[1]
    assert paths[0] == paths[1]
    expected_orbits = {"orbit-h", "orbit-he"} if species == ["*"] else {
        "orbit-h" if species == ["H"] else "orbit-he"}
    assert len(snapshots[0].columns) == 3 * len(expected_orbits)
    assert all(any(orbit in column.key for orbit in expected_orbits) for column in snapshots[0].columns)

def test_public_snapshot_preserves_columns_frame_identity_and_is_immutable():
    session = _calculated_session()
    snapshot = session.displacement_snapshot()
    assert snapshot.schema == "isodistort.displacement-snapshot.v1"
    assert [column.key for column in snapshot.columns] == list(session.mode_displacements_sc)
    basis = np.array(session._selected_subgroup.basis_vectors, dtype=float)
    expected_frame = np.mod(session.structure.frac_coords @ np.linalg.inv(basis), 1)
    np.testing.assert_allclose(snapshot.fractional_coordinates, expected_frame, atol=1e-14, rtol=0)
    for column in snapshot.columns:
        np.testing.assert_array_equal(column.fractional_vectors, session.mode_displacements_sc[column.key])
    with pytest.raises(FrozenInstanceError):
        snapshot.nmod = 1
    first = snapshot.columns[0]
    session.mode_displacements_sc[first.key][:] = 25
    assert first.fractional_vectors[0] == (1., 0., 0.)
    assert json.loads(snapshot.candidate_identity_json)["parent_sg"] == 1
    assert asdict(snapshot)["schema"] == snapshot.schema

@pytest.mark.parametrize("corruption", ["nmod", "candidate", "orphan", "nan", "scope", "types", "zero"])
def test_snapshot_rejects_changed_or_incomplete_result(corruption):
    session = _calculated_session()
    if corruption == "nmod":
        session.number_of_independent_modulations = 1
    elif corruption == "candidate":
        session._selected_subgroup.k_parameters = ["1/7"]
    elif corruption == "scope":
        session.set_distortion_scope({"displacive": ["H"]})
    elif corruption == "types":
        session.set_distortion_types(["strain"])
    elif corruption == "zero":
        session.mode_displacements_sc[list(session.mode_displacements_sc)[-1]][:] = 0
    elif corruption == "orphan":
        session.mode_displacements_sc["orphan"] = np.ones((2, 3))
    else:
        session.mode_displacements_sc[next(iter(session.mode_displacements_sc))][0, 0] = np.nan
    with pytest.raises((RuntimeError, ValueError)):
        session.displacement_snapshot()

def test_cris_to_gd_preserves_displacement_and_reciprocal_pairing(monkeypatch):
    gd = TESTS_ROOT.parent.parent.parent / "GD"
    if not (gd / "isodistort_to_gd.py").is_file():
        pytest.skip("Sibling GD converter is not available")
    monkeypatch.syspath_prepend(str(gd))
    from isodistort_to_gd import bundle_from_displacement_snapshot
    session = _calculated_session()
    snapshot = session.displacement_snapshot()
    bundle = bundle_from_displacement_snapshot(
        snapshot, cart_normalized_mode_matrix, cif=Path("independent-P1.cif"),
    )
    assert bundle.irrep == "GM1+" and bundle.opd == "P1"
    assert "unit GD coefficient" in bundle.metadata["maxamp_convention"]
    tiny = replace(snapshot, columns=tuple(replace(c, fractional_vectors=tuple(
        tuple(1e-20 * v for v in row) for row in c.fractional_vectors)) for c in snapshot.columns))
    tiny_bundle = bundle_from_displacement_snapshot(
        tiny, cart_normalized_mode_matrix, cif=Path("independent-P1.cif"))
    np.testing.assert_allclose(tiny_bundle.mode_frac, bundle.mode_frac, rtol=2e-7, atol=2e-8)
    with pytest.raises(TypeError):
        bundle_from_displacement_snapshot(snapshot, cart_normalized_mode_matrix,
                                          cif=Path("input.cif"), irrep="LD1")
    corrupted = replace(snapshot, columns=(*snapshot.columns[:-1], replace(
        snapshot.columns[-1], fractional_vectors=((0., 0., 0.), (0., 0., 0.)))))
    with pytest.raises(ValueError, match="positive finite"):
        bundle_from_displacement_snapshot(corrupted, cart_normalized_mode_matrix,
                                          cif=Path("input.cif"))
    lattice = np.asarray(snapshot.lattice_angstrom)
    # An independent metric calculation on the oblique lattice verifies scale,
    # sign and ordering. GD stores float32, so use its rounding error budget.
    for index, column in enumerate(snapshot.columns):
        raw = np.asarray(column.fractional_vectors)
        metric_norm = np.sqrt(np.sum((raw @ lattice) ** 2))
        expected = raw / metric_norm
        np.testing.assert_allclose(bundle.mode_frac[index], expected, rtol=2e-7, atol=2e-8)
        dmax = np.max(np.linalg.norm(bundle.mode_frac[index] @ lattice, axis=1))
        assert bundle.max_bounds[index] * dmax == pytest.approx(1, rel=2e-7)
    amplitudes = np.linspace(-0.2, 0.3, len(snapshot.columns))
    expected_delta = sum(a * np.asarray(c.fractional_vectors) /
                         np.sqrt(np.sum((np.asarray(c.fractional_vectors) @ lattice) ** 2))
                         for a, c in zip(amplitudes, snapshot.columns, strict=True))
    np.testing.assert_allclose(np.einsum("i,ijk->jk", amplitudes, bundle.mode_frac),
                               expected_delta, rtol=2e-7, atol=2e-8)
    h = np.array([2., -1., 3.])
    x = np.asarray(snapshot.fractional_coordinates[0])
    basis = np.array(session._selected_subgroup.basis_vectors, dtype=float)
    x_parent = x @ basis
    assert h @ x_parent == pytest.approx((h @ bundle.hkl_transform) @ x, abs=1e-12)
    # This non-symmetric basis distinguishes the correct dual map from h @ B.
    assert not np.isclose(h @ x_parent, (h @ basis) @ x)

# Original coverage: test_strain_export.py.

class _FlakyMacroIso(_TetragonalMacroIso):
    def __init__(self) -> None:
        self.macroscopic_calls = 0

    def get_macroscopic_tensor_basis(self, _parent_sg: int) -> MacroscopicTensorBasis:
        self.macroscopic_calls += 1
        if self.macroscopic_calls == 1:
            return MacroscopicTensorBasis.unresolved(
                rank_signature="[12]",
                coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
                reason="transient_test_backend_failure",
            )
        return _tetragonal_macro_basis()

@pytest.mark.parametrize(
    ("symbol", "lattice", "expected"),
    (
        ("P-1", Lattice.from_parameters(4.1, 5.2, 6.3, 76.0, 83.0, 69.0), 6),
        ("P2/m", Lattice.monoclinic(4.1, 5.2, 6.3, 107.0), 4),
        ("Pmmm", Lattice.orthorhombic(4.1, 5.2, 6.3), 3),
        ("P4/mmm", Lattice.tetragonal(4.1, 6.3), 2),
        ("P6/mmm", Lattice.hexagonal(4.1, 6.3), 2),
        ("Pm-3m", Lattice.cubic(4.1), 1),
    ),
)
def test_metric_fixed_space_gives_general_crystal_system_counts(
    symbol: str, lattice: Lattice, expected: int
) -> None:
    structure = Structure.from_spacegroup(symbol, lattice, ["H"], [[0.123, 0.234, 0.345]])
    rotations = [
        operation.rotation_matrix
        for operation in SpacegroupAnalyzer(
            structure, symprec=1.0e-3
        ).get_symmetry_operations(cartesian=False)
    ]
    result = compute_homogeneous_strain_modes(
        rotations, structure.lattice.matrix
    )

    assert len(result.modes) == expected
    assert result.diagnostics.point_group_order >= 2

def test_unresolved_macro_basis_is_not_cached_and_next_attempt_recovers() -> None:
    structure = Structure(
        Lattice.from_parameters(4.0, 5.0, 6.0, 78.0, 83.0, 71.0),
        ["H"],
        [[0.13, 0.27, 0.39]],
    )
    api = _bare_api(structure, ["strain"])
    backend = _FlakyMacroIso()
    api._iso = backend
    subgroup = _strain_export_subgroup(1, "P1")

    with pytest.raises(ValueError, match="macroscopic strain basis is unresolved"):
        api._strain_export_data_for_subgroup(subgroup)

    recovered = api._strain_export_data_for_subgroup(subgroup)

    assert recovered is not None
    assert recovered.result.export_ready
    assert backend.macroscopic_calls == 2

# Original coverage: test_structure.py.

def test_api_load_structure_poscar(tmp_path):
    """API load_structure 直接接受 POSCAR（格式兼容 → 与 CIF 识别一致）。"""
    lattice = [[3.6, 0, 0], [0, 3.6, 0], [0, 0, 3.6]]
    base = Structure(lattice, ["Fe", "Fe"],
                     [[0, 0, 0], [0.5, 0.5, 0.5]])
    poscar_path = tmp_path / "POSCAR"
    Poscar(base).write_file(str(poscar_path))
    cif_path = tmp_path / "fe.cif"
    CifWriter(base).write_file(str(cif_path))

    iso = _structure_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(poscar_path)
    sg_poscar = iso.symmetry_info["space_group_number"]

    iso2 = _structure_IsoDistort()
    iso2.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso2.load_structure(cif_path)
    sg_cif = iso2.symmetry_info["space_group_number"]

    assert sg_poscar == sg_cif == 229
    assert len(iso.structure) == len(iso2.structure)

@pytest.mark.parametrize("entrypoint", ["load", "set"])
def test_failed_parent_install_rolls_back_complete_method_state(
    tmp_path: Path,
    monkeypatch,
    entrypoint: str,
) -> None:
    """A failed new parent must not mix with the prior Method 1--4 caches."""

    # Exercise the state transaction without requiring any external wrapper.
    iso = object.__new__(_structure_IsoDistort)
    iso._sym_val = SimpleNamespace(validate=lambda _structure: {})
    old_structure = Structure(Lattice.cubic(3.0), ["Fe"], [[0, 0, 0]])
    new_structure = Structure(Lattice.cubic(4.0), ["Ni"], [[0, 0, 0]])
    old_path = (tmp_path / "old.cif").resolve()
    old_symmetry = {"space_group_number": 221, "space_group_symbol": "Pm-3m"}
    old_findsym = object()
    old_candidate = object()
    old_mode = object()
    old_distorted = old_structure.copy()
    old_displacements = {"old-mode": np.array([[0.1, 0.0, 0.0]])}
    iso.structure = old_structure
    iso.structure_path = old_path
    iso.symmetry_info = old_symmetry
    iso._findsym = old_findsym
    iso.subgroups = [old_candidate]
    iso._selected_subgroup = old_candidate
    iso.distortion_modes = [old_mode]
    iso.mode_displacements = old_displacements
    iso.distorted_structure = old_distorted
    iso._mode_cache_key = ("old",)
    iso._generated_structure_cache_key = ("old-generated",)

    new_symmetry = {"space_group_number": 221, "space_group_symbol": "Pm-3m"}
    monkeypatch.setattr(iso._sym_val, "validate", lambda _structure: new_symmetry)
    monkeypatch.setattr(iso, "parent_wyckoff_display", lambda: [])

    def reject_new_parent() -> None:
        assert iso.structure is new_structure
        assert iso.symmetry_info is new_symmetry
        # Model a failed canonicalization after it has temporarily replaced
        # both primary fields, not merely a failure before the first write.
        iso.structure = Structure(Lattice.cubic(5.0), ["Co"], [[0, 0, 0]])
        iso.symmetry_info = {"space_group_number": 1, "space_group_symbol": "P1"}
        iso._findsym = object()
        raise ValueError("inconsistent parent setting")

    monkeypatch.setattr(iso, "_attach_standard_parent_orbits", reject_new_parent)
    if entrypoint == "load":
        new_path = tmp_path / "new.cif"
        new_path.write_text("fixture intercepted by read_cif", encoding="utf-8")
        monkeypatch.setattr(core_api_module, "read_cif", lambda _path: new_structure)
        with pytest.raises(ValueError, match="inconsistent parent setting"):
            iso.load_structure(new_path)
    else:
        with pytest.raises(ValueError, match="inconsistent parent setting"):
            iso.set_structure(new_structure)

    assert iso.structure is old_structure
    assert iso.structure_path == old_path
    assert iso.symmetry_info is old_symmetry
    assert iso._findsym is old_findsym
    assert iso.subgroups == [old_candidate]
    assert iso._selected_subgroup is old_candidate
    assert iso.distortion_modes == [old_mode]
    assert iso.mode_displacements is old_displacements
    assert iso.distorted_structure is old_distorted
    assert iso._mode_cache_key == ("old",)
    assert iso._generated_structure_cache_key == ("old-generated",)

# Original coverage: test_theory_contracts.py.

def _session():
    iso = _theory_contracts_IsoDistort.__new__(_theory_contracts_IsoDistort)
    rotation = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
    iso._parent_rotations = lambda: [np.eye(3, dtype=int), rotation, rotation @ rotation]
    return iso, rotation

def test_row_lattice_orbit_uses_transposed_fractional_rotation():
    iso, rotation = _session()
    basis = np.diag([3, 1, 1])
    assert iso._same_lattice_orbit(basis @ rotation.T, basis)
    assert not iso._same_lattice_orbit(basis @ rotation, basis)
    rotations = iso._parent_rotations()
    assert _basis_is_same_lattice(basis @ rotation.T, basis, rotations)
    assert not _basis_is_same_lattice(basis @ rotation, basis, rotations)
    oriented = rational_matrix(basis @ rotation.T)
    assert is_sublattice(oriented, oriented)
    assert not is_sublattice(rational_matrix(basis @ rotation), oriented)

def test_lattice_label_does_not_change_high_denominator_basis():
    label = _theory_contracts_IsoDistort._format_lattice(((Fraction(1, 17), 0, 0), (0, 1, 0), (0, 0, 1)))
    assert label == "(1/17,0,0),(0,1,0),(0,0,1)"

def test_method1_options_do_not_read_reference_archives(monkeypatch):
    iso, _ = _session()
    iso.structure = Structure(Lattice.cubic(4), ["Si"], [[0.13, 0.27, 0.39]])
    iso.symmetry_info = {"space_group_number": 1}
    subgroup = SubgroupInfo(index=0, space_group_number=1, space_group_symbol="P1",
                            basis_vectors=np.eye(3, dtype=int).tolist(),
                            k_point_label="GM", irrep_label="GM1")
    iso._method1_filtered_subgroups = lambda: [subgroup]
    original = Path.iterdir

    def guarded_iterdir(path):
        assert path.name not in {"webpage_info", "experiment_data", "output_compare"}
        return original(path)

    monkeypatch.setattr(Path, "iterdir", guarded_iterdir)
    assert iso.method1_options()["conventional_lattices"]

@pytest.mark.parametrize("kind", ["rotational", "occupational", "magnetic"])
def test_unproven_physical_types_are_not_inferred_from_displacement(kind):
    iso, _ = _session()
    iso.distortion_types = [kind]
    with pytest.raises(NotImplementedError, match=kind):
        iso._filter_method1_by_types([], [kind])

# Original coverage: test_wsl.py.

MATCHER = StructureMatcher(ltol=1e-5, stol=1e-3, angle_tol=0.001)

COORD_TOL = 1e-5

AMP_REL_TOL = 1e-4

_wsl_DATA_DIR = experiment_data_dir()

# --- from test_golden_standard.py ---

def _srtio3_parent(tmp_path: Path) -> tuple[_wsl_IsoDistort, Structure]:
    """构造立方钙钛矿 SrTiO₃（Pm-3m #221，5 原子，a=3.905 Å）。

    位点设置对齐官网 ISODISTORT 示例（isodistortexample.php）：
    Sr 1b(1/2,1/2,1/2)、Ti 1a(0,0,0)、O 3d(1/2,0,0),(0,1/2,0),(0,0,1/2)。
    （注意：若把 O 放在 3c，iso 的 R4+ 根模式落在 3d 位点，
    本地映射器会找不到对应原子——即“空间群原点设置”陷阱。）
    """
    lattice = [[3.905, 0, 0], [0, 3.905, 0], [0, 0, 3.905]]
    parent = Structure(
        lattice,
        ["Sr", "Ti", "O", "O", "O"],
        [[0.5, 0.5, 0.5], [0, 0, 0],
         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
    )
    sg = SpacegroupAnalyzer(parent, symprec=1e-3).get_space_group_number()
    assert sg == 221, f"SrTiO₃ 母相识别为 #{sg}，应为 Pm-3m #221"

    cif = tmp_path / "srtio3.cif"
    CifWriter(parent).write_file(str(cif))

    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    assert iso.symmetry_info["space_group_number"] == 221
    return iso, parent

def _first_mode_at(iso: _wsl_IsoDistort, k: str, ir: str) -> tuple:
    """在 (k, IR) 下取首个可计算位移模式的子群，返回 (子群, 模式标号)。"""
    subs = iso.list_subgroups_at(k, ir)
    for sg in subs:
        try:
            iso.search_method_2(subgroup_idx=sg.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式/计算失败的子群跳过
            continue
        if iso.mode_displacements:
            return sg, next(iter(iso.mode_displacements))
    raise AssertionError(f"{k} {ir}: 无可计算位移模式的子群")

def _distorted_sg(iso: _wsl_IsoDistort, label: str, amplitude: float) -> int:
    d = iso.generate_distortion(irrep_label=label, amplitude=amplitude)
    return SpacegroupAnalyzer(d, symprec=1e-3).get_space_group_number()

@_wsl_required
# ------------------------------------------------------------------
# 金标准相变路径（文献公认子群）
# ------------------------------------------------------------------

def test_srtio3_r4_plus_tilt_path(tmp_path):
    """R₄⁺ 倾转（a⁰a⁰c⁻）→ I4/mcm #140（文献公认）。

    校验链：子群枚举 → 模式计算 → 生成 → spglib 验证 = 子群 → StructureMatcher
    零振幅回退 = 母相。
    """
    iso, parent = _srtio3_parent(tmp_path)
    subgroup, label = _first_mode_at(iso, "R", "R4+")
    assert subgroup.space_group_number == 140, \
        f"R4+ 首个子群应为 I4/mcm #140，实际 {subgroup.space_group_number} {subgroup.space_group_symbol}"

    # 非零振幅：畸变结构 spglib 对称性 == 目标子群
    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 140, f"R4+ 畸变结构对称性 #{d_sg}，应为 #140"

    # 零振幅：StructureMatcher 结构等价 + 空间群回母相
    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert SpacegroupAnalyzer(d0, symprec=1e-3).get_space_group_number() == 221
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"

@_wsl_required
def test_srtio3_m3_plus_inphase_tilt_path(tmp_path):
    """M₃⁺ 同相倾转 → P4/mbm #127（文献公认）。"""
    iso, parent = _srtio3_parent(tmp_path)
    subgroup, label = _first_mode_at(iso, "M", "M3+")
    assert subgroup.space_group_number == 127, \
        f"M3+ 首个子群应为 P4/mbm #127，实际 {subgroup.space_group_number} {subgroup.space_group_symbol}"

    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 127, f"M3+ 畸变结构对称性 #{d_sg}，应为 #127"

    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"

@_wsl_required
def test_srtio3_gm4_minus_polar_path(tmp_path):
    """Γ₄⁻ 极性位移 → P4mm #99（文献公认的铁电软模路径）。"""
    iso, parent = _srtio3_parent(tmp_path)
    subgroup, label = _first_mode_at(iso, "GM", "GM4-")
    assert subgroup.space_group_number == 99, \
        f"GM4- 首个子群应为 P4mm #99，实际 {subgroup.space_group_number} {subgroup.space_group_symbol}"

    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 99, f"GM4- 畸变结构对称性 #{d_sg}，应为 #99"

    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"

@_wsl_required
# ------------------------------------------------------------------
# 数值精度（容差符合科研标准）
# ------------------------------------------------------------------

def test_amplitude_linearity_cartesian(tmp_path):
    """振幅加倍 -> 笛卡尔位移长度线性加倍（相对误差 ≤ 1e-4）。

    用 Γ 点模式（GM4-，超胞因子 1）保证母相/畸变原子一一对应。
    """
    iso, parent = _srtio3_parent(tmp_path)
    _, label = _first_mode_at(iso, "GM", "GM4-")

    d1 = iso.generate_distortion(irrep_label=label, amplitude=0.1)
    d2 = iso.generate_distortion(irrep_label=label, amplitude=0.2)
    assert len(d1) == len(d2) == len(parent)

    cart1 = np.asarray([s.coords for s in d1])
    cart2 = np.asarray([s.coords for s in d2])
    cart0 = np.asarray([s.coords for s in parent])
    disp1, disp2 = [], []

    def _displacement_length(coords: np.ndarray,
                            base: np.ndarray,
                            lattice_matrix: np.ndarray) -> float:
        best = np.inf
        for j in range(len(coords)):
            delta = coords[j] - base
            frac = np.linalg.solve(lattice_matrix.T, delta)
            frac -= np.round(frac)
            d_cart = frac @ lattice_matrix
            best = min(best, float(np.linalg.norm(d_cart)))
        return best

    for i in range(len(parent)):
        base = cart0[i]
        disp1.append(_displacement_length(cart1, base, parent.lattice.matrix))
        disp2.append(_displacement_length(cart2, base, parent.lattice.matrix))
    for v1, v2 in zip(disp1, disp2, strict=True):
        if v1 > 1e-9:
            ratio = v2 / v1
            assert abs(ratio - 2.0) <= AMP_REL_TOL, \
                f"振幅线性失败：位移比 {ratio:.6f}（期望 2.0）"

@_wsl_required
def test_zero_amplitude_coordinate_tolerance(tmp_path):
    """零振幅：畸变结构原子必须与母相一一对应（周期最小镜像 ≤ 1e-5）。

    子群基矢可能是旋转幺模胞（如 GM4- → P4mm），分数坐标随格子基
    变化而重排，不能按行号比对（“原子顺序问题”陷阱）；按同物种 +
    周期最小镜像距离建立原子对应后比较。
    """
    iso, parent = _srtio3_parent(tmp_path)
    _, label = _first_mode_at(iso, "GM", "GM4-")
    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert len(d0) == len(parent)
    # 结构语义等价（含原子重排/周期平移）
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"

    # 原子级对应：同物种 + 最小镜像距离（周期约化）
    pc = np.asarray(parent.frac_coords, dtype=float)
    dc = np.asarray(d0.frac_coords, dtype=float)
    max_disp = 0.0
    for i in range(len(parent)):
        best = float("inf")
        for j in range(len(d0)):
            if parent[i].species_string != d0[j].species_string:
                continue
            delta = dc[j] - pc[i]
            delta -= np.round(delta)
            best = min(best, float(np.linalg.norm(delta)))
        max_disp = max(max_disp, best)
    assert max_disp <= COORD_TOL, \
        f"零振幅原子级最大位移 {max_disp:.2e} 超过 {COORD_TOL}"

@_wsl_required
# ------------------------------------------------------------------
# 区边界 k 点相位（带心母相的副本反号）回归
# ------------------------------------------------------------------

def test_im3m_h4_minus_zone_boundary_phase(tmp_path):
    """Im-3m #229 的 H4-（k=(1,1,1) 惯用坐标）区边界模式。

    回归保护：bcc 惯用胞 2 原子（2a 轨道）相差原始格点平移
    (1/2,1/2,1/2)。BUSH 对该模式输出两个**相位反号**的代表点
    `(0,0,0)→(1,0,1)` 与 `(-1/2,1/2,1/2)→(-1,0,-1)`，二者互为
    “模格点等价”。若映射器按“首个代表胜出”分配，两个原子会得到
    相同位移（刚性平移，畸变退化为母相 #229，对称性校验失败）。
    修复后按“周期等价（mod 1）”分配：原子 0→(1,0,1)、原子 1→(-1,0,-1)，
    畸变结构正确降为子群对称性（H4- 首个子群 P4_2/nmc #129）。
    """
    lattice = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]]
    parent = Structure(
        lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    sg = SpacegroupAnalyzer(parent, symprec=1e-3).get_space_group_number()
    assert sg == 229, f"Im-3m 母相识别为 #{sg}，应为 #229"

    cif = tmp_path / "im3m.cif"
    CifWriter(parent).write_file(str(cif))
    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)

    subgroup, label = _first_mode_at(iso, "H", "H4-")
    assert subgroup.space_group_number == 129, \
        f"H4- 首个子群应为 P4_2/nmc #129，实际 " \
        f"{subgroup.space_group_number} {subgroup.space_group_symbol}"

    # 关键回归：两个原子的位移必须相位反号（非均匀），否则退化为刚性平移
    disp = iso.mode_displacements[label]["displacements"]
    assert not np.allclose(disp, disp[0], atol=1e-6), \
        "H4- 带心副本位移被映射为均匀（刚性平移）——相位丢失回归"
    assert np.allclose(disp[0], -disp[1], atol=1e-6), \
        "H4- 两个 bcc 原子的位移应为相位反号（±(1,1,0)/√2 方向）"

    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 129, f"H4- 畸变结构对称性 #{d_sg}，应为 #129"

    # 零振幅回退母相
    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert SpacegroupAnalyzer(d0, symprec=1e-3).get_space_group_number() == 229

@_wsl_required
# ------------------------------------------------------------------
# 官网参考 CIF（LD1 零振幅框架文件）金标准比对
# ------------------------------------------------------------------

def test_official_ld1_reference_structure(tmp_path):
    """官网 LD1_C1_subgroup.cif：零振幅框架文件，读回 I4/mmm #139。

    与本地 API 对 EuAl4 的识别结果做结构语义比对（StructureMatcher 允许
    超胞倍数差异，此处校验原子数/元素/空间群信息一致）。
    """
    official = _wsl_DATA_DIR / "LD1_C1_subgroup.cif"
    if not official.exists():
        pytest.skip("官网参考 CIF 不存在")

    ref = Structure.from_file(str(official))
    ref_sg = SpacegroupAnalyzer(ref, symprec=1e-3).get_space_group_number()
    # 零振幅框架文件：序参量全 0，读回应为母相 I4/mmm #139
    assert ref_sg == 139, f"官网 LD1 参考读回 #{ref_sg}，应为 I4/mmm #139"
    # 框架文件为 6 倍 c 超胞（s=12 -> 60 原子）；原子数验证
    assert len(ref) == 60, f"官网 LD1 参考应为 60 原子（6×10），实际 {len(ref)}"
    assert {s.species_string for s in ref} <= {"Eu", "Al"}

    # 本地 EuAl4 母相识别
    parent_cif = _wsl_DATA_DIR / "EuAl4 Parent.cif"
    if not parent_cif.exists():
        pytest.skip("EuAl4 母相 CIF 不存在")
    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(parent_cif)
    assert iso.symmetry_info["space_group_number"] == 139
    # 官网框架文件的 6×c 超胞与母相 10 原子：原子总数比例一致
    assert len(ref) == 6 * len(iso.structure), \
        "官网 LD1 参考原子数应为母相 6 倍（c 轴 6 倍超胞）"

@_wsl_required
def test_method1_eual4_matches_official_opd_html():
    """Method 1 OPD identities and exact affine subgroups match official evidence."""
    from tests import validate_method_outputs as validator
    from tests.official_html_resolver import scan_official_html

    from backend.utils.lattice import determinant, inverse, multiply
    from backend.utils.affine_embeddings import affine_equivalence, embedding_from_identity, parent_affine_group

    official_root = (
        TESTS_ROOT.parent.parent
        / "webpage_info"
        / "EuAl4 Parent.cif"
    )
    assert official_root.is_dir(), "official parent evidence directory not present"
    cif = _wsl_DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("EuAl4 Parent.cif not present")

    class _OpdForm(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=True)
            self.fields = {}
            self.lines = []
            self._post = False
            self._line = None

        def handle_starttag(self, tag, attrs):
            attributes = dict(attrs)
            if tag == "form":
                self._post = attributes.get("method", "").casefold() == "post"
            elif tag == "input" and self._post:
                name = attributes.get("name")
                if name:
                    self.fields.setdefault(name, []).append(attributes.get("value", ""))
                if name == "orderparam":
                    assert attributes.get("type", "").casefold() == "radio"
                    assert self._line is None, "unterminated official OPD line"
                    self._line = []
            elif tag == "br" and self._line is not None:
                self.lines.append("".join(self._line).strip())
                self._line = None

        def handle_endtag(self, tag):
            if tag == "form":
                self._post = False

        def handle_data(self, data):
            if self._line is not None:
                self._line.append(data)

    # Capture once through the shared evidence owner. Its generic roles do
    # not yet include OPD pages, so distinguish Method 1 by content, not names.
    inventory = scan_official_html(official_root, recursive=False)
    matches = []
    for page in inventory.pages:
        if page.title != "isodistort: order parameter direction":
            continue
        parsed = _OpdForm()
        parsed.feed(page.raw_html)
        parsed.close()
        if parsed.fields.get("input") == ["distort"] \
                and parsed.fields.get("origintype") == ["method1"]:
            matches.append((page, parsed))
    assert len(matches) == 1, (
        "expected one content-identified Method 1 OPD page; "
        f"matches={[str(page.path) for page, _ in matches]}"
    )
    page, parsed = matches[0]
    official = parsed.lines
    assert official and len(official) == len(parsed.fields.get("orderparam", []))
    assert len(set(official)) == len(official), "duplicate official OPD identities"

    def field(name):
        values = parsed.fields.get(name, [])
        assert len(values) == 1, f"missing or ambiguous official context field {name!r}"
        return values[0].strip()

    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)

    # Verify the complete parent and active distortion scope before comparing
    # candidates. No scientific condition is inferred from the archive path.
    assert int(field("spacegroup").split()[0]) == iso.symmetry_info["space_group_number"]
    assert np.allclose(
        [float(value) for value in field("dlattparam").split()],
        iso.structure.lattice.parameters, atol=1e-8, rtol=0.0,
    )
    def normalize(text):
        return " ".join(text.split())

    count = int(field("wycount"))
    official_sites = [
        f"{field(f'wyatom{index:03d}')} {field(f'wyckoff{index:03d}')}"
        for index in range(1, count + 1)
    ]
    assert [normalize(line) for line in official_sites] == [
        normalize(line) for line in iso.parent_wyckoff_display()
    ]
    assert normalize("Default space-group preferences: " + iso.space_group_preferences()) \
        in normalize(page.raw_html)
    type_species = {}
    for index in range(1, count + 1):
        atom_type = int(field(f"wytype{index:03d}"))
        species = field(f"wyatomtype{index:03d}")
        assert type_species.setdefault(atom_type, species) == species
    assert set(type_species.values()) == set(iso.species())
    enabled = {
        name for name, values in parsed.fields.items()
        if name.startswith("include") and any(value.casefold() == "true" for value in values)
    }
    assert enabled == {"includestrain", *(
        f"includedisplacive{atom_type:03d}" for atom_type in type_species
    )}, "official query is not strain plus displacive on every parent species"

    cands = iso.search_method_1(
        distortion_types=["displacive", "strain"],
        crystal_system=None, subgroup_space_group=None,
        maximal_subgroup_only=False,
    )
    local = [it.subgroup.opd_line().strip() for it in cands]
    assert len(local) == len(official), (
        f"Method 1 count {len(local)} != parsed official {len(official)}"
    )
    assert len(set(local)) == len(local), "duplicate local OPD identities"

    def parse_identity(line):
        details = validator._parse_subgroup_details_identity(line)
        # The shared parser owns rational SG/B/O/s/i/k semantics. This thin
        # text adapter retains the three Method 1 IR/OPD/direction fields.
        head = line.split(", basis=", 1)[0].split()
        assert len(head) == 5 and re.fullmatch(r"\([^)]*\)", head[2])
        assert int(head[3]) == details.space_group_number
        assert head[4].casefold() == details.space_group_symbol
        key = (
            *head[:3], details.space_group_number, details.space_group_symbol,
            details.sublattice_index, details.group_index, details.k_active,
        )
        origin_tokens = re.findall(r"origin=(\([^)]*\))", line)
        assert len(origin_tokens) == 1, "missing or ambiguous candidate origin"
        # The shared identity parser reduces origin modulo one for matching.
        # Restore its exact raw value for the actual Seitz embedding: an
        # integer parent-origin shift need not preserve a child subgroup.
        origin = validator._parse_vectors(origin_tokens[0], 1)[0]
        return key, details, origin

    def index_lines(lines):
        indexed = {}
        for line in lines:
            key, details, origin = parse_identity(line)
            assert key not in indexed, f"duplicate complete OPD identity: {key}"
            indexed[key] = (line, details, origin)
        return indexed

    official_by_identity = index_lines(official)
    local_by_identity = index_lines(local)
    missing = set(official_by_identity) - set(local_by_identity)
    extra = set(local_by_identity) - set(official_by_identity)
    assert not missing, f"missing complete identities: {missing}"
    assert not extra, f"extra complete identities: {extra}"
    parent = parent_affine_group(
        iso.structure,
        symprec=iso.cfg.symmetry_cartesian_tolerance_angstrom,
        angle_tolerance_degrees=iso.cfg.symmetry_angle_tolerance_degrees,
    )
    validation_records = []
    for key, (official_line, official_details, official_origin) in official_by_identity.items():
        local_line, local_details, local_origin = local_by_identity[key]
        assert validator.bases_generate_same_lattice(
            local_details.basis, official_details.basis,
        ), f"non-unimodular conventional basis change: {key}"
        embeddings = [
            embedding_from_identity({
                "space_group_number": details.space_group_number,
                "basis": details.basis,
                "origin": origin,
            }, parent)
            for details, origin in (
                (official_details, official_origin), (local_details, local_origin)
            )
        ]
        official_embedding, local_embedding = embeddings
        assert official_embedding.hall_number == local_embedding.hall_number
        relationship = affine_equivalence(
            local_embedding, official_embedding, parent, allow_parent_conjugacy=False,
        )
        assert relationship == "equal", f"different exact affine subgroup: {key}"
        for embedding in embeddings:
            size = abs(determinant(multiply(embedding.lattice, inverse(parent.lattice))))
            assert size == official_details.sublattice_index
            assert size * len(parent.operations) / len(embedding.operations) == \
                official_details.group_index
        validation_records.append({
            "identity": key,
            "official_line": official_line,
            "local_line": local_line,
            "official_basis": official_details.basis,
            "local_basis": local_details.basis,
            "official_origin": official_origin,
            "local_origin": local_origin,
            "U_local_official": multiply(local_details.basis, inverse(official_details.basis)),
            "hall_number": official_embedding.hall_number,
            "seitz_counts": [len(official_embedding.operations), len(local_embedding.operations)],
            "relationship": relationship,
        })
    assert len(validation_records) == len(official_by_identity)

@_wsl_required
def test_eual4_rootless_wyckoff_secondary_modes_are_completed():
    """BUSH root modes and child-fixed rootless-orbit modes form one basis."""
    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(_wsl_DATA_DIR / "EuAl4 Parent.cif")
    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    for irrep, opd, expected in (
        ("GM4+", "P1", 1),
        ("P2", "C1", 5),
        ("X1-", "P1", 4),
    ):
        target = next(
            subgroup for subgroup in candidates
            if subgroup.irrep_label == irrep and subgroup.opd_symbol == opd
        )
        iso.search_method_2(
            target.index,
            distortion_type=["displacive"],
            candidates=candidates,
        )
        assert len(iso.mode_displacements) == expected
        assert any(mode.wyckoff_site == "e" for mode in iso.distortion_modes)

@_wsl_required
def test_ndnio2_rootless_oxygen_secondary_modes_are_completed():
    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(_wsl_DATA_DIR / "NdNiO2 own.cif")
    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    target = next(
        subgroup for subgroup in candidates
        if subgroup.irrep_label == "R1-" and subgroup.opd_symbol == "C1"
    )
    iso.search_method_2(
        target.index,
        distortion_type=["displacive"],
        candidates=candidates,
    )
    assert len(iso.mode_displacements) == 4
    oxygen = [mode.irrep_label for mode in iso.distortion_modes if mode.wyckoff_site == "f"]
    assert oxygen == ["M1+", "M2+"]

@_wsl_required
def test_4310_n1plus_4d1_complete_modes_cover_the_child_cell():
    """A multi-arm special-k route must not lose part of a rooted orbit."""
    cif = _wsl_DATA_DIR / "4310_tetra.cif"
    official = (
        TESTS_ROOT.parent.parent
        / "output_compare"
        / "4310_tetra.cif"
        / "官网"
        / "Method1"
        / "N1+_4D1_SG2"
        / "data.isoviz"
    )
    if not cif.is_file() or not official.is_file():
        pytest.skip("validated 4310 parent/official Method 1 reference unavailable")

    official_text = official.read_text(encoding="utf-8-sig")
    mode_block = official_text.split("!displacivemodelist", 1)[1]
    expected = len(re.findall(r"(?m)^\s*\d+\s+\d+\s+[-+0-9.Ee]+\s+", mode_block))
    assert expected > 0

    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    target = next(
        subgroup for subgroup in candidates
        if subgroup.irrep_label == "N1+" and subgroup.opd_symbol == "4D1"
    )

    iso.search_method_2(
        target.index,
        distortion_type=["displacive"],
        candidates=candidates,
    )

    assert len(iso.mode_displacements_sc) == expected

@_wsl_required
def test_shifted_4310_parent_is_canonicalized_before_method1_and_modes():
    """An equivalent global origin shift must not change Method 1 completeness."""
    cif = _wsl_DATA_DIR / "4310_tetra.cif"
    official = (
        TESTS_ROOT.parent.parent
        / "output_compare"
        / "4310_tetra.cif"
        / "官网"
        / "Method1"
        / "N1+_4D1_SG2"
        / "data.isoviz"
    )
    if not cif.is_file() or not official.is_file():
        pytest.skip("validated 4310 parent/official Method 1 reference unavailable")

    shifted = read_cif(cif)
    shifted.translate_sites(
        range(len(shifted)),
        [0.0, 0.0, 0.1],
        frac_coords=True,
        to_unit_cell=True,
    )
    iso = _wsl_IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.set_structure(shifted)

    assert iso.symmetry_info["parent_orbit_standardization"]["status"] == "canonicalized"
    assert iso.symmetry_info["space_group_number"] == 139
    assert "Ni1 2a (0,0,0)" in iso.parent_wyckoff_display()

    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    assert len(candidates) == 125
    target = next(
        subgroup for subgroup in candidates
        if subgroup.irrep_label == "N1+" and subgroup.opd_symbol == "4D1"
    )
    iso.search_method_2(
        target.index,
        distortion_type=["displacive"],
        candidates=candidates,
    )

    official_text = official.read_text(encoding="utf-8-sig")
    mode_block = official_text.split("!displacivemodelist", 1)[1]
    expected = len(re.findall(r"(?m)^\s*\d+\s+\d+\s+[-+0-9.Ee]+\s+", mode_block))
    assert len(iso.mode_displacements_sc) == expected

@_wsl_required
def test_primitive_centered_parent_is_canonicalized_as_a_complete_cell():
    conventional = Structure.from_spacegroup(
        225,
        [[5.64, 0, 0], [0, 5.64, 0], [0, 0, 5.64]],
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    primitive = SpacegroupAnalyzer(
        conventional, symprec=1e-3
    ).get_primitive_standard_structure()
    assert len(primitive) == 2

    iso = _wsl_IsoDistort()
    standardized = iso.set_structure(primitive)

    assert iso.symmetry_info["space_group_number"] == 225
    assert iso.symmetry_info["parent_orbit_standardization"]["status"] == "canonicalized"
    assert len(standardized) == 8
    assert {
        (site["multiplicity"], site["wyckoff_letter"], site["species"])
        for site in iso.symmetry_info["wyckoff_sites"]
    } == {(4, "a", "Na"), (4, "b", "Cl")}

def _gold_ld1_c1_labels() -> list[str]:
    root = TESTS_ROOT.parent.parent / "output_compare"
    if not root.is_dir():
        return []
    hits = list(root.rglob("data.isoviz"))
    for path in hits:
        if "ld1" in str(path).lower() and "c1" in str(path).lower():
            text = path.read_text(encoding="utf-8", errors="replace")
            if "!displacivemodelist" not in text:
                continue
            block = text.split("!displacivemodelist", 1)[1]
            labels = []
            for line in block.splitlines():
                if "dsp" not in line:
                    continue
                bits = line.split()
                for bit in bits:
                    if ":dsp]" in bit:
                        labels.append(bit)
                        break
            if labels:
                return labels
    return []

@_wsl_required
def test_eual4_ld1_c1_parametric_complete_modes():
    """smodes/(3+d) lock-in modes for EuAl4 LD g=1/6, LD1 C1.

    Does not hardcode 48.  Compares core tokens to official gold isoviz
    when ``output_compare`` is present; otherwise checks label families.
    """
    cif = _wsl_DATA_DIR / "EuAl4 Parent.cif"
    if not cif.is_file():
        pytest.skip("EuAl4 Parent.cif is not available")
    iso = _wsl_IsoDistort()
    iso.set_distortion_types(["strain", "displacive"])
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    try:
        subs = iso.list_subgroups_at(
            "LD", "LD1", k_parameters=["1/6"], generate_if_missing=False,
        )
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"LD1 g=1/6 subgroup list unavailable: {exc}")
    target = next(
        (s for s in subs if s.irrep_label == "LD1" and s.opd_symbol == "C1"),
        None,
    )
    if target is None:
        pytest.skip("LD1 C1 is not in the cached subgroup list")
    result = iso.search_method_2(
        target.index,
        candidates=subs,
        number_of_independent_modulations=0,
    )
    assert result.modes
    assert all(
        mode.mode_identity is not None
        and mode.mode_identity.status == "verified"
        and mode.mode_identity.source == "iso_microscopic"
        and mode.microscopic_provenance is not None
        for mode in result.modes
    )
    labels = list((iso._mode_label_overrides or {}).values())
    assert labels, result.metadata.get("parametric_note") or "no complete modes"
    joined = " ".join(labels)
    assert "LD1" in joined
    assert any("0,0,1/6" in lab or "1/6" in lab for lab in labels)
    assert any("GM" in lab for lab in labels)
    assert any("M1+" in lab or "M3-" in lab for lab in labels)
    gold = _gold_ld1_c1_labels()
    if gold:
        local_core = {_mode_core_token(lab) for lab in labels}
        gold_core = {_mode_core_token(lab) for lab in gold}
        missing = gold_core - local_core
        extra = local_core - gold_core
        assert not missing, f"missing vs official: {sorted(missing)[:12]}"
        assert not extra, f"extra vs official: {sorted(extra)[:12]}"

    payload = iso.export_subgroups_zip(
        formats=["cif", "isoviz", "modes", "topas"],
        subgroups=[target],
        compute_missing_modes=True,
        export_method=2,
        number_of_independent_modulations=0,
    )
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        names = archive.namelist()
        assert names == [
            "LD1_C1/subgroup.cif",
            "LD1_C1/data.isoviz",
            "LD1_C1/Complete modes details.txt",
            "LD1_C1/topas.str",
            "export_candidate_status.json",
            "export_candidate_status.txt",
        ]
        assert all(archive.getinfo(name).file_size > 0 for name in names)
        assert b"_iso_displacivemode_ID" in archive.read(names[0])
        assert b"!displacivemodelist" in archive.read(names[1])
        assert b"Displacive mode definitions" in archive.read(names[2])
        topas = archive.read(names[3]).lower()
        assert b"site " in topas and b"prm " in topas
        status = json.loads(archive.read("export_candidate_status.json"))
        assert status["schema"] == "isodistort.export-candidate-status.v1"
        assert status["method"] == 2
        assert status["status"] == "success"
        assert status["requested_candidate_count"] == status["successful_candidate_count"] == 1
        assert status["ineligible_candidate_count"] == 0
        assert status["unresolved_identity_candidate_count"] == 0
        assert status["cancelled_candidate_count"] == 0
        assert status["failures"] == []
        assert b"requested=1 successful=1 ineligible=0" in archive.read(
            "export_candidate_status.txt"
        )
