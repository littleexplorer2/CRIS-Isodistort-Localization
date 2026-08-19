"""
物理自洽性自检（Scientific Self-Consistency Checks）

无需与官网比对即可发现计算错误的内置校验（对应“科研级测试方案”第一层第 3 节）：

1. 零振幅回退：振幅=0 时畸变结构必须与母相完全一致（机器精度内）
2. 子群规则：畸变后结构的空间群必须是母相空间群的子群
   （其对称操作必须属于母相对称操作的子集）
3. 模式正交性：不同不可约表示的位移模式向量应正交
4. 线性关系：振幅加倍 -> 位移向量同步线性加倍
5. 对称性守恒：畸变结构的所有对称操作必须是母相对称操作的子集
   （spglib 旋转矩阵包含性判定）

任意一项失败即说明计算逻辑存在 bug（不是“误差”）。
"""
from __future__ import annotations

import numpy as np
import spglib
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

SYMPREC = 1e-3


def _symmetry_rotations(structure: Structure) -> set[tuple]:
    """结构全部对称操作的旋转矩阵（元组集合，用于子集判定）。"""
    dataset = spglib.get_symmetry_dataset(
        (structure.lattice.matrix, structure.frac_coords,
         structure.atomic_numbers),
        symprec=SYMPREC, angle_tolerance=1.0,
    )
    rots = set()
    for rot in dataset["rotations"]:
        rots.add(tuple(map(tuple, rot.tolist())))
    return rots


def check_symmetry_conservation(parent: Structure,
                                distorted: Structure) -> tuple[bool, str]:
    """畸变结构的旋转操作 ⊆ 母相旋转操作。"""
    p_rots = _symmetry_rotations(parent)
    d_rots = _symmetry_rotations(distorted)
    missing = d_rots - p_rots
    ok = not missing
    detail = "" if ok else f"畸变结构有 {len(missing)} 个操作不在母相中"
    return ok, detail


def check_subgroup_rule(parent_sg: int, distorted_sg: int,
                        parent: Structure, distorted: Structure) -> tuple[bool, str]:
    """畸变空间群必须与目标子群一致，且其操作 ⊆ 母相操作。"""
    ok_sym, detail = check_symmetry_conservation(parent, distorted)
    if not ok_sym:
        return False, detail
    # 畸变结构空间群号不允许高于母相（降对称）
    if distorted_sg > 0 and parent_sg > 0 and distorted_sg > parent_sg:
        return False, f"畸变空间群 #{distorted_sg} 高于母相 #{parent_sg}"
    return True, ""


def check_zero_amplitude(iso, irrep_label: str,
                         tolerance: float = 1e-6) -> tuple[bool, str]:
    """振幅=0 时输出结构必须与母相完全一致。"""
    parent = iso.structure
    distorted = iso.generate_distortion(irrep_label=irrep_label,
                                        amplitude=0.0)
    if len(distorted) != len(parent):
        return False, "原子数不一致"
    d = distorted.frac_coords - parent.frac_coords
    d -= np.round(d)
    max_disp = float(np.max(np.abs(d)))
    ok = max_disp <= tolerance
    return ok, f"最大坐标差 {max_disp:.2e}"


def check_linearity(iso, irrep_label: str) -> tuple[bool, str]:
    """振幅加倍 -> 位移向量长度线性加倍（畸变模式是线性分解）。"""
    d1 = iso.generate_distortion(irrep_label=irrep_label, amplitude=0.1)
    d2 = iso.generate_distortion(irrep_label=irrep_label, amplitude=0.2)
    if len(d1) != len(d2):
        return False, "原子数不一致"
    delta1 = d1.frac_coords - iso.structure.frac_coords
    delta2 = d2.frac_coords - iso.structure.frac_coords
    delta1 -= np.round(delta1)
    delta2 -= np.round(delta2)
    if np.max(np.abs(delta1)) < 1e-9:
        return False, "零位移模式（无物理畸变）"
    ratio = np.linalg.norm(delta2) / np.linalg.norm(delta1)
    ok = abs(ratio - 2.0) < 1e-3
    return ok, f"位移比 = {ratio:.6f}（期望 2.0）"


def check_mode_orthogonality(mode_displacements: dict) -> tuple[bool, str]:
    """不同 IR 的位移模式向量正交（内积为 0）。"""
    labels = list(mode_displacements.keys())
    if len(labels) < 2:
        return True, "仅 1 个模式，无需正交性校验"
    max_abs = 0.0
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a = np.asarray(mode_displacements[labels[i]]["displacements"]).reshape(-1)
            b = np.asarray(mode_displacements[labels[j]]["displacements"]).reshape(-1)
            if a.shape != b.shape:
                return False, "模式向量维度不一致"
            max_abs = max(max_abs, abs(float(np.dot(a, b))))
    ok = max_abs < 1e-6
    return ok, f"最大内积绝对值 = {max_abs:.2e}"


def run_self_checks(iso, irrep_label: str,
                    distorted: Structure | None = None) -> dict:
    """对指定模式运行全部自检，返回 {检查名: (通过, 说明)}。"""
    parent = iso.structure
    if distorted is None:
        distorted = iso.generate_distortion(irrep_label=irrep_label,
                                            amplitude=0.1)
    parent_sg = SpacegroupAnalyzer(parent, symprec=SYMPREC).get_space_group_number()
    dist_sg = SpacegroupAnalyzer(distorted, symprec=SYMPREC).get_space_group_number()
    results = {
        "zero_amplitude": check_zero_amplitude(iso, irrep_label),
        "subgroup_rule": check_subgroup_rule(parent_sg, dist_sg,
                                             parent, distorted),
        "symmetry_conservation": check_symmetry_conservation(parent, distorted),
        "linearity": check_linearity(iso, irrep_label),
        "mode_orthogonality": check_mode_orthogonality(iso.mode_displacements),
    }
    return {
        # bool(ok)：各检查函数可能返回 numpy 布尔（如 np.True_），
        # 统一转为原生 bool，保证 JSON 序列化等下游使用
        name: {"ok": bool(ok), "detail": detail}
        for name, (ok, detail) in results.items()
    }
