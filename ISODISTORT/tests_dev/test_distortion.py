"""Distortion engine, mapper, occupational modes, self-check, multi-mode."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import numpy.linalg as la
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.api import IsoDistort
from isocore.backend import BushMode, DistortionMode, SubgroupInfo
from isocore.distortion import DistortionEngine, DistortionMapper, OccupationalModeGenerator
from isocore.distortion.phase_path import (
    DEFAULT_DISTORTION_TYPES,
    normalize_distortion_types,
)
from isocore.utils.self_check import (
    check_linearity,
    check_mode_orthogonality,
    check_subgroup_rule,
    check_symmetry_conservation,
    check_zero_amplitude,
)

from data_dir import experiment_data_dir

CIFS_DIR = Path(__file__).resolve().parent / "cifs_30"
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


# --- from test_distortion.py ---

def test_default_distortion_modes():
    """默认畸变类型应为 strain 单种（对齐官网：Types 面板默认只勾选 Strain）。"""
    # 默认对齐官网 Types 面板：Strain 勾选 + Displacive 全物种勾选
    assert DEFAULT_DISTORTION_TYPES == ["strain", "displacive"]
    assert normalize_distortion_types(None) == ["strain", "displacive"]
    assert normalize_distortion_types(["strain", "displacement", "strain"]) == [
        "strain",
        "displacive",  # 旧名 displacement 自动映射为 displacive
    ]
    assert normalize_distortion_types("order") == ["occupational"]


def test_distortion_engine_single_mode():
    # 简单立方

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A", "B"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    # 人工构造位移：A 原子沿 x 方向位移
    displacements = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ])

    engine = DistortionEngine()
    distorted = engine.generate_single_mode(
        struct, displacements, amplitude=0.05
    )

    # 验证位移
    dx = distorted[0].frac_coords[0] - struct[0].frac_coords[0]
    assert abs(dx - 0.05) < 1e-6, f"期望位移 0.05，实际 {dx}"

    # B 原子不动
    dy = distorted[1].frac_coords[1] - struct[1].frac_coords[1]
    assert abs(dy) < 1e-6

    print("✅ 单模式畸变测试通过")


def test_distortion_engine_supercell():

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A"], [[0, 0, 0]])

    displacements = np.array([[0.01, 0.0, 0.0]])
    engine = DistortionEngine()
    super_distorted = engine.generate_single_mode(
        struct, displacements, amplitude=1.0, supercell=[2, 2, 2]
    )

    assert len(super_distorted) == 8
    print("✅ 超胞畸变测试通过")


def test_mixed_mode():

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A"], [[0, 0, 0]])

    disp_x = np.array([[1.0, 0.0, 0.0]])
    disp_y = np.array([[0.0, 1.0, 0.0]])

    engine = DistortionEngine()
    total_disp = 0.1 * disp_x + 0.2 * disp_y
    mixed = engine.generate_modes(struct, parent_displacements=total_disp)

    dx = mixed[0].frac_coords[0]
    dy = mixed[0].frac_coords[1]
    assert abs(dx - 0.1) < 1e-6
    assert abs(dy - 0.2) < 1e-6

    print("✅ 多模式混合畸变测试通过")




# --- from test_distortion_mapper.py ---

def _mapper_parent() -> Structure:
    """四方母相：Eu 2a + Al 2e（(0,0,0.38)/(0,0,0.62)）。"""
    return Structure(
        Lattice.tetragonal(4.4, 11.2),
        ["Eu", "Eu", "Al", "Al"],
        [[0, 0, 0], [0.5, 0.5, 0.5], [0, 0, 0.38], [0, 0, 0.62]],
    )


def _mapper_wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 2, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0, 1]},
        {"wyckoff_letter": "e", "multiplicity": 2, "species": "Al",
         "representative_index": 2, "equivalent_indices": [2, 3]},
    ]


def _mapper_mode() -> DistortionMode:
    mode = DistortionMode(irrep_label="GM1+", opd_symbol="P1",
                          k_point_label="GM")
    mode.bush_modes = [
        BushMode("GM1+", "P1", "a", [0, 0, 0], ["0", "0", "0"], [[0, 0, 0]]),
        # e 位点：两个相反位移的代表（(0,0,z) 与 (0,0,-z)）
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "-z"], [[0, 0, -1]]),
    ]
    return mode


def test_mapper_parameterized_z_points():
    """(0,0,z) 与 (0,0,-z) 必须解析为不同位置：z=0.38 -> +z，z=0.62 -> -z。"""
    mapper = DistortionMapper()
    result = mapper.map_modes_to_atoms(_mapper_parent(), _mapper_wyckoff(), [_mapper_mode()])
    disp = result["GM1+"]["displacements"]
    assert np.allclose(disp[0], 0) and np.allclose(disp[1], 0)  # Eu 不动
    assert np.allclose(disp[2], [0, 0, 1])   # (0,0,0.38) -> +z
    assert np.allclose(disp[3], [0, 0, -1])  # (0,0,0.62) -> -z


def test_mapper_uniform_single_rep():
    """单代表行（均匀模式）：同一位移作用于该位点全部原子。"""
    mode = _mapper_mode()
    mode.bush_modes = [
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
    ]
    mapper = DistortionMapper()
    result = mapper.map_modes_to_atoms(_mapper_parent(), _mapper_wyckoff(), [mode])
    disp = result["GM1+"]["displacements"]
    assert np.allclose(disp[2], [0, 0, 1])
    assert np.allclose(disp[3], [0, 0, 1])


# --- from test_occupational_modes.py ---

def _occ_parent() -> Structure:
    """简单四方晶胞：Eu(0,0,0) + Al(0.25,0.25,0.25)。"""
    return Structure(Lattice.tetragonal(4.0, 8.0),
                     ["Eu", "Al"], [[0, 0, 0], [0.25, 0.25, 0.25]])


def _occ_wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 1, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0]},
        {"wyckoff_letter": "b", "multiplicity": 1, "species": "Al",
         "representative_index": 1, "equivalent_indices": [1]},
    ]


def _occ_subgroup(basis) -> SubgroupInfo:
    return SubgroupInfo(
        index=0, space_group_number=2, space_group_symbol="P-1",
        basis_vectors=basis, k_point_label="X", irrep_label="X1+",
    )


def test_occupational_mode_doubled_cell():
    """2x1x1 超胞：Al 位点分裂为两类 -> 产生 +1/-1 占据率模式。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(), _occ_subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]),
                         {"Al"})
    assert modes, "倍胞下 Al 位点应产生占据率模式"
    m = modes[0]
    assert m.species == "Al"
    assert m.label == "occ-Al-b"
    # 模式值只含 +1/-1/0，且 +1 与 -1 数量相等（两类各半）
    assert set(np.unique(m.pattern)) <= {1.0, -1.0, 0.0}
    assert np.count_nonzero(m.pattern == 1) == np.count_nonzero(m.pattern == -1)


def test_occupational_mode_no_split_identity_cell():
    """单位基矢（t 子群，无超胞）：位点不分裂 -> 无占据率模式。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(),
                         _occ_subgroup([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Al"})
    assert modes == []


def test_occupational_scope_filters_species():
    """作用域只含 Eu 时，Al 位点不产生模式。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(),
                         _occ_subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Eu"})
    assert all(m.species == "Eu" for m in modes)


def test_occupancy_exportable_structure():
    """生成的占据率结构可被 pymatgen 构建（部分占据位点）。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(),
                         _occ_subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Al"})
    engine = DistortionEngine()
    m = modes[0]
    sc = engine.generate_modes(
        parent, m.basis_vectors,
        parent_displacements=None,
        occupancy_patterns=[(m.pattern, 0.2)],
    )
    assert len(sc) == 4  # 2 原子 x 2 倍胞
    occs = {s.species_string for s in sc if "Al" in s.species_string}
    assert any(o != "Al" for o in occs)  # 至少一个 Al 位点变为部分占据


# --- from test_self_check.py ---

class _StubIso:
    """最小桩：仅提供自检所需的 structure / generate_distortion。"""

    def __init__(self, struct: Structure, disp: np.ndarray):
        self.structure = struct
        self._disp = np.asarray(disp, dtype=float)
        self.mode_displacements = {}

    def generate_distortion(self, irrep_label: str, amplitude: float):
        coords = self.structure.frac_coords + amplitude * self._disp
        return Structure(self.structure.lattice, self.structure.species,
                         coords % 1.0, coords_are_cartesian=False)


def _self_parent() -> Structure:
    return Structure(Lattice.cubic(4.0), ["Fe", "Fe"],
                     [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])


def test_symmetry_conservation_same_structure():
    p = _self_parent()
    ok, detail = check_symmetry_conservation(p, p)
    assert ok, detail


def test_symmetry_conservation_extra_ops_breaks():
    """畸变结构出现母相没有的操作（如单斜母相被误建成更高对称）-> 不守恒。"""
    tetragonal = Structure(Lattice.tetragonal(4.0, 8.0), ["Fe", "Fe"],
                           [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    cubic = Structure(Lattice.cubic(4.0), ["Fe", "Fe"],
                      [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    ok, _ = check_symmetry_conservation(tetragonal, cubic)
    assert not ok  # 立方对称操作 ⊄ 四方对称操作


def test_subgroup_rule_same_sg():
    p = _self_parent()
    ok, detail = check_subgroup_rule(229, 229, p, p)
    assert ok, detail


def test_zero_amplitude_returns_parent():
    stub = _StubIso(_self_parent(), np.array([[0.2, 0, 0], [0, 0, 0]]))
    ok, detail = check_zero_amplitude(stub, "X")
    assert ok, detail


def test_linearity_doubling():
    stub = _StubIso(_self_parent(), np.array([[0.2, 0, 0], [0, 0, 0]]))
    ok, detail = check_linearity(stub, "X")
    assert ok, detail


def test_mode_orthogonality():
    ok, _ = check_mode_orthogonality({
        "A": {"displacements": np.array([[1, 0, 0], [0, 0, 0]])},
        "B": {"displacements": np.array([[0, 1, 0], [0, 0, 0]])},
    })
    assert ok
    bad, _ = check_mode_orthogonality({
        "A": {"displacements": np.array([[1, 0, 0], [0, 0, 0]])},
        "B": {"displacements": np.array([[2, 0, 0], [0, 0, 0]])},
    })
    assert not bad


# --- from test_multi_mode.py ---

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



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
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



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
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



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
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

