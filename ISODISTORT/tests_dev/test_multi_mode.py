"""
多模式耦合测试（第一层 §1 多模式用例 + 第二层接口一致性）。

验证：
1. 同时激活 2~3 个不可约表示的混合畸变：位移 = 各模式线性叠加；
2. 混合畸变空间群仍为所选子群（对称性不被破坏）；
3. API / 引擎 两路径结果一致（混合生成走同一 generate_modes 内核）。
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import numpy.linalg as la
import pytest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.api import IsoDistort

CIFS_DIR = Path(__file__).resolve().parent / "cifs_30"
DATA_DIR = Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码")


def _wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        result = subprocess.run(  # noqa: PLW1510
            ["wsl.exe", "--status"],  # noqa: S607
            capture_output=True, text=True, timeout=30,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, OSError):
        return False


pytestmark = pytest.mark.skipif(
    not _wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试"
)


def _load_iso(path: Path) -> IsoDistort:
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(path)
    return iso


def _subgroup_with_modes(iso: IsoDistort, min_modes: int = 2) -> tuple:
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


def test_mixed_mode_linear_superposition(tmp_path):
    """混合畸变 = 各单模式位移的线性叠加（StructureMatcher 原子对应后比对）。

    位移向量以母相原胞为单位定义（mode_displacements）；混合生成的超胞
    原子位移（相对母相，周期约化）应等于各模式位移的相位调制线性叠加。
    子群基矢可能是旋转幺模胞，原子顺序会重排，故用 StructureMatcher
    建立“混合畸变原子 -> 母相原子”的对应后再比对位移。
    """
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
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


def test_mixed_mode_symmetry_conserved(tmp_path):
    """多模式混合畸变：spglib 对称性 == 目标子群（叠加不破坏对称性）。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
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


def test_mixed_mode_invalid_label_raises(tmp_path):
    """混合生成传入未知模式标号：明确报错，不静默。"""
    cif = CIFS_DIR / "sg001.cif"
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
