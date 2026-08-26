"""smodes 封装：参数 k 点上的位移模式 IR 活性（对齐官网 Search 阶段过滤）。

官网在 Method 2 枚举子群前，会按 Distortion Types + 物种作用域过滤
仅对结构有位移模式的 IR；对参数 k 点（LD/DT 等）使用 (3+d) 超空间 /
smodes 机制，而 iso 的 DISPLAY BUSH 无法完成该判断。

本模块调用 ``isobyu/smodes``，解析各 IR 下哪些元素有非零对称位移，
供 ``IsoDistort.list_subgroups_at_kpoint`` 在搜索阶段剔除无相关模式的 IR
（例如 EuAl4 LD g=1/6 + displacive(Eu,Al) 仅保留 LD1/LD2/LD5）。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

from pymatgen.core import Structure

from .base_wrapper import BaseWrapper


@dataclass
class SmodesIrrepActivity:
    """单个 IR 在 smodes 输出中的活性摘要。"""

    label: str
    species: set[str] = field(default_factory=set)


class SmodesWrapper(BaseWrapper):
    """调用 smodes 判断参数 k 点各 IR 的位移模式活性。"""

    def __init__(self) -> None:
        super().__init__()
        self.binary = str(Path(self.cfg.iso_bin).parent / "smodes")

    def active_irreps(self,
                      structure: Structure,
                      space_group_number: int,
                      wyckoff_sites: list[dict],
                      k_label: str,
                      k_parameters: list[str] | None = None,
                      species_filter: set[str] | None = None) -> set[str] | None:
        """
        返回在指定 k 点上对结构（可选：限定物种）有位移模式的 IR 标签集合。

        Args:
            structure: 母相结构（取晶格参数）
            space_group_number: 空间群号
            wyckoff_sites: symmetry_validator 输出的 Wyckoff 列表
            k_label: Miller-Love k 点（如 ``LD``）
            k_parameters: 官网约定的参数值（如 ``["1/6"]`` 表示 g=1/6）
            species_filter: 若给定，仅当 IR 在该物种集合内至少一处有非零
                位移时才视为活性

        Returns:
            IR 标签集合；成功但无活性 IR 时返回空集。smodes 失败时返回
            ``None``（调用方应跳过该 k 点的过滤，避免把“无模式”当成“未知”）。
        """
        try:
            stdout = self.run_stdin(self.binary, self._build_input(
                structure, space_group_number, wyckoff_sites,
                k_label, k_parameters,
            ))
        except Exception:  # noqa: BLE001 - 过滤为尽力而为，失败时不阻断枚举
            return None

        activities = self._parse_irrep_activities(stdout)
        if not species_filter:
            return {a.label for a in activities}

        active: set[str] = set()
        for act in activities:
            if act.species & species_filter:
                active.add(act.label)
        return active

    @staticmethod
    def _build_input(structure: Structure,
                     space_group_number: int,
                     wyckoff_sites: list[dict],
                     k_label: str,
                     k_parameters: list[str] | None) -> str:
        lat = structure.lattice
        lines = [
            structure.composition.reduced_formula,
            str(space_group_number),
            f"{lat.a:.6f} {lat.b:.6f} {lat.c:.6f} "
            f"{lat.alpha:.6f} {lat.beta:.6f} {lat.gamma:.6f}",
            str(len(wyckoff_sites)),
        ]
        for site in wyckoff_sites:
            idx = site["representative_index"]
            frac = structure[idx].frac_coords
            letter = site["wyckoff_letter"]
            species = site["species"]
            # smodes：Wyckoff 字母 + 分数坐标（无自由度坐标可省略）
            coord_str = " ".join(f"{x:.6f}".rstrip("0").rstrip(".") for x in frac)
            if coord_str.replace("0", "").replace(".", "").replace(" ", "") == "":
                lines.append(f"{species} {letter}")
            else:
                lines.append(f"{species} {letter} {coord_str}")
        lines.append("1")
        k_line = k_label.strip()
        if k_parameters:
            k_line += " " + " ".join(str(p).strip() for p in k_parameters)
        lines.append(k_line)
        return "\n".join(lines) + "\n"

    @staticmethod
    def _parse_irrep_activities(stdout: str) -> list[SmodesIrrepActivity]:
        """解析 smodes 输出，提取每个 IR 涉及的非零位移物种。"""
        activities: list[SmodesIrrepActivity] = []
        current: SmodesIrrepActivity | None = None
        in_modes = False
        irrep_re = re.compile(r"^Irrep\s+(\S+)")
        atom_re = re.compile(
            r"^\s*\d+\s+(\S+)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
            r"\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
            r"\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$"
        )

        for line in stdout.splitlines():
            m_ir = irrep_re.match(line.strip())
            if m_ir:
                current = SmodesIrrepActivity(label=m_ir.group(1))
                activities.append(current)
                in_modes = False
                continue
            if current is None:
                continue
            if "Symmetry modes:" in line:
                in_modes = True
                continue
            if not in_modes:
                continue
            stripped = line.strip()
            if stripped.startswith("*") or stripped.startswith("---"):
                continue
            if "displacement" in stripped and "atom" in stripped:
                continue
            m_atom = atom_re.match(line)
            if not m_atom:
                if stripped.startswith("Irrep") or stripped.startswith("*"):
                    in_modes = False
                continue
            species = m_atom.group(1)
            disp = [float(m_atom.group(i)) for i in range(2, 5)]
            if any(abs(d) > 1e-8 for d in disp):
                # 物种符号可能带电荷后缀，取元素前缀
                elem = re.match(r"([A-Z][a-z]?)", species)
                if elem:
                    current.species.add(elem.group(1))

        return activities

    @staticmethod
    def smodes_kparams(k_parameters: list[str] | None) -> list[str] | None:
        """smodes 输入使用官网物理 k 参数（分数），与 UI 一致。"""
        if not k_parameters:
            return None
        return [str(Fraction(str(p).strip())) for p in k_parameters]
