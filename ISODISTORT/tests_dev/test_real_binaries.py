"""
测试：真实二进制（iso / findsym）端到端冒烟。

依赖 WSL（isobyu 为 Linux ELF 二进制）。WSL 不可用时自动跳过。
"""

import shutil
import subprocess

import pytest

from isocore.backend import FindsymWrapper, IsoWrapper, SubgroupInfo


def _wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        # 探测 WSL 可用性：返回码非零即视为不可用，不抛异常
        result = subprocess.run(  # noqa: PLW1510
            ["wsl.exe", "--status"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, OSError):
        return False


pytestmark = pytest.mark.skipif(
    not _wsl_available(), reason="WSL 不可用，跳过真实二进制测试"
)


def test_findsym_identifies_nacl():
    """findsym 识别 NaCl（F 心）应为 Fm-3m #225。"""
    fs = FindsymWrapper()
    result = fs.identify(
        lattice_params=[5.63, 5.63, 5.63, 90, 90, 90],
        atom_types=["Na", "Cl"],
        atom_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        centering="F",
    )
    assert result.space_group_number == 225
    assert result.space_group_symbol == "Fm-3m"
    assert {s["wyckoff_letter"] for s in result.wyckoff_sites} == {"a", "b"}


def test_iso_kpoints_and_subgroups():
    """iso 枚举 SG 225 的 k 点与 GM5- 子群。"""
    iso = IsoWrapper()
    kpoints = iso.list_k_points(225)
    labels = {kp.label for kp in kpoints}
    assert {"GM", "L", "X", "W"}.issubset(labels)

    subgroups = iso.list_subgroups(225, "GM", "GM5-")
    assert len(subgroups) >= 1
    # GM5- P1 对应子群 I-42m (#121)，指数 6
    p1 = next(sg for sg in subgroups if sg.opd_symbol == "P1")
    assert p1.space_group_number == 121
    assert p1.subgroup_index == 6


def test_iso_modes_and_domains():
    """BUSH 模式基矢与畴列表。"""
    iso = IsoWrapper()

    target = SubgroupInfo(
        index=0, space_group_number=107, space_group_symbol="I4mm",
        subgroup_index=6, size=1, is_maximal=True,
        opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
        basis_vectors=[[0, 0.5, -0.5], [0, 0.5, 0.5], [1, 0, 0]],
        origin=[0, 0, 0], k_point_label="GM", irrep_label="GM4-",
    )
    modes = iso.calc_distortion_modes(225, target, wyckoff_letters=["a"])
    assert len(modes) >= 1
    assert modes[0].bush_modes, "GM4- P1 在 4a 位点应存在位移模式"

    domains = iso.get_domains(225, target)
    assert len(domains) == 6
    assert domains[0].domain_number == 1
