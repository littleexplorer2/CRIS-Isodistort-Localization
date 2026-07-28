"""
FINDSYM 封装 - 空间群识别与 Wyckoff 位点分析

对应阶段一，步骤2：母相空间群与 Wyckoff 位点识别
封装方式：✅ 直接封装 findsym 二进制
"""
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple

from .base_wrapper import BaseWrapper
from ..utils import (
    parse_space_group_number,
    parse_wyckoff_sites,
    OutputParseError,
)


@dataclass
class FindsymResult:
    """findsym 计算结果"""
    space_group_number: int          # 空间群号
    space_group_symbol: str = ""     # 空间群 Hermann-Mauguin 符号
    wyckoff_sites: List[dict] = field(default_factory=list)  # Wyckoff 位点信息
    lattice_params: Tuple[float, ...] = ()  # 晶格参数 a,b,c,alpha,beta,gamma
    raw_output: str = ""             # 原始输出（调试用）


class FindsymWrapper(BaseWrapper):
    """
    FINDSYM 程序封装

    功能：给定晶格参数与原子坐标，识别空间群、分配 Wyckoff 位点
    """

    def __init__(self):
        super().__init__()
        self.binary = self.cfg.findsym_bin

    def identify(self, lattice_params: List[float],
                 atom_types: List[str],
                 atom_positions: List[List[float]],
                 centering: str = "P",
                 title: str = "findsym input") -> FindsymResult:
        """
        调用 findsym 识别空间群

        Args:
            lattice_params: [a, b, c, alpha, beta, gamma] 晶格参数
            atom_types: 原子种类列表，如 ["Na", "Cl"]
            atom_positions: 原子分数坐标列表，与 atom_types 一一对应
            centering: 点阵中心类型 P/I/F/A/B/C/R
            title: 输入文件标题

        Returns:
            FindsymResult: 识别结果
        """
        assert len(atom_types) == len(atom_positions), \
            "atom_types 与 atom_positions 长度必须一致"

        input_text = self._build_input(
            lattice_params, atom_types, atom_positions, centering, title
        )

        stdout = self._run_binary(self.binary, input_text)

        return self._parse_output(stdout)

    # ---- 输入生成 ----

    @staticmethod
    def _build_input(lattice_params, atom_types, atom_positions,
                     centering, title) -> str:
        """生成 findsym 输入文件内容（关键字格式）"""
        a, b, c, alpha, beta, gamma = lattice_params
        n_atoms = len(atom_types)

        lines = [
            "!useKeyWords",
            "!title",
            title,
            "!latticeParameters",
            f"{a} {b} {c} {alpha} {beta} {gamma}",
            "!unitCellCentering",
            centering,
            "!atomCount",
            str(n_atoms),
            "!atomType",
            " ".join(atom_types),
            "!atomPosition",
        ]
        for pos in atom_positions:
            lines.append(f"{pos[0]} {pos[1]} {pos[2]}")

        return "\n".join(lines) + "\n"

    # ---- 输出解析 ----

    @staticmethod
    def _parse_output(text: str) -> FindsymResult:
        """解析 findsym 输出"""
        sg_num = parse_space_group_number(text)
        if sg_num is None:
            raise OutputParseError("findsym", "未找到空间群号")

        # 空间群符号
        import re
        sym_match = re.search(r"Space Group\s+\d+\s+(\S+)\s+(\S+)", text)
        sg_symbol = sym_match.group(2) if sym_match else ""

        # Wyckoff 位点
        sites = parse_wyckoff_sites(text)

        # 晶格参数
        lat_match = re.search(
            r"Lattice parameters.*?:\s*\n\s*([\d.\s]+)", text
        )
        lattice_params = ()
        if lat_match:
            vals = [float(x) for x in lat_match.group(1).split()]
            if len(vals) == 6:
                lattice_params = tuple(vals)

        return FindsymResult(
            space_group_number=sg_num,
            space_group_symbol=sg_symbol,
            wyckoff_sites=sites,
            lattice_params=lattice_params,
            raw_output=text,
        )
