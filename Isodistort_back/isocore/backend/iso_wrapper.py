"""
iso 封装：子群枚举、畸变模式计算、位点分裂、畴变体

对应：
- 阶段二，步骤4：生成母群的所有各向同性子群列表
- 阶段三，步骤6：计算目标相变的畸变模式（序参量基矢）
- 阶段三，步骤7（理论部分）：Wyckoff 位点分裂分析
- 阶段五，步骤10（对称部分）：畴变体对称操作矩阵
"""

import re
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Sequence

from .base_wrapper import BaseWrapper
from ..utils import parse_subgroup_list


def normalize_distortion_types(distortion_types=None) -> List[str]:
    """
    标准化畸变类型输入，默认返回 displacement + strain。

    Args:
        distortion_types: 可为 None、单个字符串或字符串列表
    
    Returns:
        List[str]: 标准化后的畸变类型列表
    
    相对路径： isocore/backend/iso_wrapper.py
    """
    from ..distortion.phase_path import DISTORTION_TYPES, DEFAULT_DISTORTION_TYPES

    # 生成result列表，包含所有传入的畸变类型
    if distortion_types is None:
        result = DEFAULT_DISTORTION_TYPES.copy()
    elif isinstance(distortion_types, str):
        result = [distortion_types]  # 传入单个字符串：封装成单元素列表
    else:
        result = list(distortion_types)  # 传入列表 / 可迭代对象：强制转为 list 实例

    # 标准化：去除空字符串、忽略不合法畸变类型、去重
    normalized = []
    seen = set()  # 收集所有遇到过的畸变类型，用于去重
    for raw_type in result:
        key = str(raw_type).strip().lower()
        if not key or key not in DISTORTION_TYPES:  # 忽略空字符串或不合法畸变类型
            continue
        if key not in seen:  # 去重：只保留第一次出现的畸变类型
            normalized.append(key)
            seen.add(key)

    if not normalized:
        return DEFAULT_DISTORTION_TYPES.copy()
    return normalized


@dataclass
class SubgroupInfo:
    """子群信息"""

    index: int                        # iso 输出表里的子群序号，也是后续给 iso 交互命令时要输入的编号
    space_group_number: int           # 空间群号
    space_group_symbol: str = ""      # 空间群符号
    subgroup_index: int = 0           # 路径对象里保存的“已选子群编号”，与 index 值相同，只是容器不同
    origin_choice: int = 1            # 原点选择
    transformation_matrix: List[List[float]] = field(default_factory=list)  # 坐标变换矩阵


@dataclass
class DistortionMode:
    """畸变模式（序参量）"""

    irrep_label: str                  # 不可约表示标号，如 GM4-
    dimension: int = 1                # 模式维度
    mode_type: str = "displacement"   # 类型: displacement/order/strain/magnetic
    basis_vectors: List[List[float]] = field(default_factory=list)  # 基矢
    wyckoff_site: str = ""            # 对应 Wyckoff 位点


@dataclass
class IsoResult:
    """iso 计算结果基类"""

    raw_output: str = ""


class IsoWrapper(BaseWrapper):
    """
    ISOTROPY (iso) 命令行程序封装

    iso 是交互式程序，通过输入命令序列来控制计算流程。
    本封装将常见功能封装为独立方法，自动生成对应的命令序列。

    核心能力：
    1. list_subgroups: 枚举母群的所有各向同性子群
    2. calc_distortion_modes: 计算指定相变路径的畸变模式基矢
    3. get_site_splitting: 计算 Wyckoff 位点在子群下的分裂
    4. get_domain_operations: 获取畴变体的对称操作
    """

    def __init__(self):
        """Relative path: isocore/backend/iso_wrapper.py"""

        super().__init__()
        self.binary = self.cfg.iso_bin

    # ================================================================
    # 功能 1：枚举所有各向同性子群
    # ================================================================

    def list_subgroups(self, parent_sg: int,
                    distortion_type: str | Sequence[str] | None = None) -> List[SubgroupInfo]:
        """
        列出给定空间群的所有各向同性子群。
        默认采用 displacement + strain 两种模式；也兼容单一字符串输入。
        
        Args:
            parent_sg: 母相空间群号
            distortion_type: 畸变类型，可为 None、单个字符串或字符串列表

        计算过程：
        1. 进入母群
        2. 设置畸变类型
        3. 列出子群

        Returns:
            List[SubgroupInfo]: 子群信息列表

        相对路径： isocore/backend/iso_wrapper.py
        """

        types = normalize_distortion_types(distortion_type)
        all_subgroups: List[SubgroupInfo] = []  # 输出的结果列表
        seen: set[tuple[int,int]] = set()  # 用于收录所有已见过的子群 (space_group_number, index) 元组，避免重复添加

        for mode in types:
            cmd_seq = self._cmd_enter_sg(parent_sg)
            cmd_seq += self._cmd_set_distortion_type(mode)
            cmd_seq += "s\n"  # 列出子群
            cmd_seq += "q\n"  # 退出

            stdout = self._run_binary(self.binary, cmd_seq)
            parsed = self._parse_subgroup_list(stdout)
            # 过滤重复子群：iso 输出中可能出现同一子群在不同畸变类型下的重复条目
            for sg in parsed:
                key = (sg.space_group_number, sg.index)
                if key not in seen:
                    all_subgroups.append(sg)
                    seen.add(key)

        return all_subgroups

    # ================================================================
    # 功能 2：计算畸变模式（序参量基矢）
    # ================================================================

    def calc_distortion_modes(self, parent_sg: int, subgroup_idx: int,
                            distortion_type: str = "displacement",
                            wyckoff_letter: str = None) -> List[DistortionMode]:
        """
        计算指定相变路径下的畸变模式基矢

        Args:
            parent_sg: 母相空间群号
            subgroup_idx: 子群序号（由 list_subgroups 返回的 index）
            distortion_type: 畸变类型
            wyckoff_letter: 指定 Wyckoff 位点字母（可选，不指定则全部）

        计算过程：
        1. 进入母群
        2. 设置畸变类型
        3. 选择子群
        4. 可选：指定 Wyckoff 位点
        5. 输出畸变模式

        Returns:
            List[DistortionMode]: 畸变模式列表
        
        Relative path: isocore/backend/iso_wrapper.py"""

        cmd_seq = self._cmd_enter_sg(parent_sg)
        cmd_seq += self._cmd_set_distortion_type(distortion_type)
        cmd_seq += f"{subgroup_idx}\n"  # 选择子群
        if wyckoff_letter:
            cmd_seq += f"w {wyckoff_letter}\n"
        cmd_seq += "o\n"  # 显示序参量/畸变模式
        cmd_seq += "q\n"  # 退出

        stdout = self._run_binary(self.binary, cmd_seq)
        return self._parse_distortion_modes(stdout)

    # ================================================================
    # 功能 3：Wyckoff 位点分裂分析
    # ================================================================

    def get_site_splitting(self, parent_sg: int, subgroup_idx: int,
                        wyckoff_letter: str) -> List[Dict]:
        """
        计算母相 Wyckoff 位点在子群下的分裂情况

        Args:
            parent_sg: 母相空间群号
            subgroup_idx: 子群序号
            wyckoff_letter: 母相位点 Wyckoff 字母

        计算过程：
        1. 进入母群
        2. 选择子群
        3. 指定位点进行分裂分析

        Returns:
            list of dict: 每个子位点包含 wyckoff_letter, multiplicity, site_symmetry
        
        Relative path: isocore/backend/iso_wrapper.py"""

        cmd_seq = self._cmd_enter_sg(parent_sg)
        cmd_seq += f"{subgroup_idx}\n"
        cmd_seq += f"p {wyckoff_letter}\n"  # 位点分裂
        cmd_seq += "q\n"

        stdout = self._run_binary(self.binary, cmd_seq)
        return self._parse_site_splitting(stdout)

    # ================================================================
    # 功能 4：畴变体对称操作
    # ================================================================

    def get_domain_operations(self, parent_sg: int,
                            subgroup_idx: int) -> List[np.ndarray]:
        """
        获取各畴变体对应的对称操作矩阵

        Args:
            parent_sg: 母相空间群号
            subgroup_idx: 子群序号

        计算过程：
        1. 进入母群
        2. 选择子群

        Returns:
            list of 4x4 仿射变换矩阵
        
        Relative path: isocore/backend/iso_wrapper.py"""


        cmd_seq = self._cmd_enter_sg(parent_sg)
        cmd_seq += f"{subgroup_idx}\n"
        cmd_seq += "d\n"  # 畴信息
        cmd_seq += "q\n"

        stdout = self._run_binary(self.binary, cmd_seq)
        return self._parse_domain_ops(stdout)

    # ================================================================
    # 命令序列生成辅助
    # ================================================================

    @staticmethod
    def _cmd_enter_sg(sg_number: int) -> str:
        """进入指定空间群的命令序列

        Relative path: isocore/backend/iso_wrapper.py"""

        return f"{sg_number}\n"

    @staticmethod
    def _cmd_set_distortion_type(dtype: str) -> str:
        """设置畸变类型命令

        Relative path: isocore/backend/iso_wrapper.py"""

        mapping = {
            "displacement": "1",  # 原子位移
            "order": "2",         # 原子有序
            "strain": "3",        # 晶格应变
            "magnetic": "4",      # 磁矩
        }
        return f"t {mapping.get(dtype, '1')}\n"

    # ================================================================
    # 输出解析
    # ================================================================

    @staticmethod
    def _parse_subgroup_list(text: str) -> List[SubgroupInfo]:
        """解析子群列表输出

        Relative path: isocore/backend/iso_wrapper.py"""

        subgroups_raw = parse_subgroup_list(text)
        result = []
        for sg in subgroups_raw:
            result.append(SubgroupInfo(
                index=sg["index"],
                space_group_number=sg["space_group_number"],
                space_group_symbol=sg["space_group_symbol"],
            ))
        return result

    @staticmethod
    def _parse_distortion_modes(text: str) -> List[DistortionMode]:
        """解析畸变模式输出
        
        Relative path: isocore/backend/iso_wrapper.py"""

        modes = []
        lines = text.splitlines()
        current_mode = None

        for line in lines:
            stripped = line.strip()
            # 匹配不可约表示行，如 "GM1+  1-dim"
            ir_match = re.match(r"^([A-Z][A-Z0-9+-]+)\s+(\d+)-dim", stripped)
            # 填充modes列表
            if ir_match:
                # 如果当前没有模式，则创建新的 DistortionMode
                if current_mode:
                    modes.append(current_mode)
                # 如果当前已有模式，则先保存
                current_mode = DistortionMode(
                    irrep_label=ir_match.group(1),
                    dimension=int(ir_match.group(2)),
                )
                continue

            # 匹配 Wyckoff 位点行，如 "Wyckoff site: 4a"
            if current_mode and stripped and re.match(r"^[-0-9.\s/]+$", stripped):
                # 基矢行
                vals = []
                for tok in stripped.split():
                    if "/" in tok:
                        num, den = tok.split("/")
                        vals.append(float(num) / float(den))
                    else:
                        vals.append(float(tok))
                if vals:
                    current_mode.basis_vectors.append(vals)

        if current_mode:
            modes.append(current_mode)

        return modes

    @staticmethod
    def _parse_site_splitting(text: str) -> List[Dict]:
        """解析位点分裂输出

        Relative path: isocore/backend/iso_wrapper.py"""

        sites = []
        lines = text.splitlines()
        in_section = False
        for line in lines:
            stripped = line.strip()
            if "Wyckoff" in stripped and "split" in stripped.lower():
                in_section = True
                continue
            if in_section and re.match(r"^\s*[a-z]\s+\d+", stripped):
                parts = stripped.split()
                if len(parts) >= 2:
                    sites.append({
                        "wyckoff_letter": parts[0],
                        "multiplicity": int(parts[1]),
                    })
            if in_section and stripped.startswith("---"):
                break
        return sites

    @staticmethod
    def _parse_domain_ops(text: str) -> list:
        """解析畴对称操作（占位实现，待完善）

        Relative path: isocore/backend/iso_wrapper.py"""

        import numpy as np
        # 简化实现：返回单位矩阵列表，实际需解析对称操作
        return [np.eye(4)]
