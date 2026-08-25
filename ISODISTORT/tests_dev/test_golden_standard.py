"""
金标准回归测试（第一层：核心算法正确性）——pymatgen StructureMatcher 语义比对。

科研比对原则（见 README「30-CIF 科学验证」）：
- 结构比对用 ``StructureMatcher(ltol=1e-5, stol=1e-3, angle_tol=0.001)``，
  自动处理原子顺序 / 对称等价 / 数值容差，不逐行比对坐标；
- 非数值项（空间群编号、IR 标记）完全精确匹配；
- 数值项（分数坐标 ≤1e-5、振幅线性 ≤1e-4、零振幅回退 机器精度）。

本文件以经典钙钛矿 SrTiO₃（Pm-3m #221）为金标准母相，验证：
  R₄⁺（八面体倾转 a⁰a⁰c⁻，子群 I4/mcm #140）
  M₃⁺（面内同相倾转，子群 P4/mbm #127）
  Γ₄⁻（极性位移，子群 P4mm #99）
三条文献公认的相变路径；并做零振幅回退 / 振幅线性 / StructureMatcher
结构等价断言。依赖 WSL（iso/findsym 二进制），不可用时自动跳过。
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.api import IsoDistort

# 科研通用容差（对应测试方案第一层 §2）
MATCHER = StructureMatcher(ltol=1e-5, stol=1e-3, angle_tol=0.001)
COORD_TOL = 1e-5
AMP_REL_TOL = 1e-4

from data_dir import experiment_data_dir

DATA_DIR = experiment_data_dir()


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
    not _wsl_available(), reason="WSL 不可用，跳过金标准真实计算测试"
)


def _srtio3_parent(tmp_path: Path) -> tuple[IsoDistort, Structure]:
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

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    assert iso.symmetry_info["space_group_number"] == 221
    return iso, parent


def _first_mode_at(iso: IsoDistort, k: str, ir: str) -> tuple:
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


def _distorted_sg(iso: IsoDistort, label: str, amplitude: float) -> int:
    d = iso.generate_distortion(irrep_label=label, amplitude=amplitude)
    return SpacegroupAnalyzer(d, symprec=1e-3).get_space_group_number()


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
    iso = IsoDistort(language="en")
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


# ------------------------------------------------------------------
# 官网参考 CIF（LD1 零振幅框架文件）金标准比对
# ------------------------------------------------------------------

def test_official_ld1_reference_structure(tmp_path):
    """官网 LD1_C1_subgroup.cif：零振幅框架文件，读回 I4/mmm #139。

    与本地 API 对 EuAl4 的识别结果做结构语义比对（StructureMatcher 允许
    超胞倍数差异，此处校验原子数/元素/空间群信息一致）。
    """
    official = DATA_DIR / "LD1_C1_subgroup.cif"
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
    parent_cif = DATA_DIR / "EuAl4 Parent.cif"
    if not parent_cif.exists():
        pytest.skip("EuAl4 母相 CIF 不存在")
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(parent_cif)
    assert iso.symmetry_info["space_group_number"] == 139
    # 官网框架文件的 6×c 超胞与母相 10 原子：原子总数比例一致
    assert len(ref) == 6 * len(iso.structure), \
        "官网 LD1 参考原子数应为母相 6 倍（c 轴 6 倍超胞）"
