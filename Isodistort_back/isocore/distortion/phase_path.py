"""
相变路径管理 - 子群选择、畸变类型、参数组装与校验

对应阶段二，步骤5：相变路径参数组装
实现方式：❌ 自研
"""
from dataclasses import dataclass, field
from typing import List, Optional

from ..utils import PhasePathError
from ..backend import SubgroupInfo


DISTORTION_TYPES = {
    "displacement": "原子位移畸变",
    "order": "原子有序化畸变",
    "strain": "晶格应变畸变",
    "magnetic": "磁矩畸变",
}


@dataclass
class PhasePath:
    """
    相变路径参数对象

    描述从母相到子相的完整相变路径，包含：
    - 母相空间群
    - 目标子群（序号 + 空间群信息）
    - 畸变类型
    - 原点选择
    - 超胞参数
    """
    parent_sg_number: int
    subgroup_index: int
    subgroup_sg_number: int = 0
    subgroup_symbol: str = ""
    distortion_type: str = "displacement"
    origin_choice: int = 1
    supercell: List[int] = field(default_factory=lambda: [1, 1, 1])
    selected_irreps: List[str] = field(default_factory=list)  # 指定不可约表示

    def validate(self) -> bool:
        """参数合法性校验

        Relative path: isocore/distortion/phase_path.py"""

        if self.parent_sg_number < 1 or self.parent_sg_number > 230:
            raise PhasePathError(f"母相空间群号 {self.parent_sg_number} 不合法")

        if self.distortion_type not in DISTORTION_TYPES:
            raise PhasePathError(
                f"不支持的畸变类型: {self.distortion_type}, "
                f"可选: {list(DISTORTION_TYPES.keys())}"
            )

        if len(self.supercell) != 3 or any(s < 1 for s in self.supercell):
            raise PhasePathError(f"超胞参数不合法: {self.supercell}")

        return True

    @classmethod
    def from_subgroup(cls, parent_sg: int, subgroup: SubgroupInfo,
                    distortion_type: str = "displacement") -> "PhasePath":
        """从子群信息快速构建相变路径

        Relative path: isocore/distortion/phase_path.py"""

        return cls(
            parent_sg_number=parent_sg,
            subgroup_index=subgroup.index,
            subgroup_sg_number=subgroup.space_group_number,
            subgroup_symbol=subgroup.space_group_symbol,
            distortion_type=distortion_type,
        )

    def describe(self) -> str:
        """人类可读的路径描述

        Relative path: isocore/distortion/phase_path.py"""

        return (
            f"空间群 #{self.parent_sg_number} → #{self.subgroup_sg_number} "
            f"({self.subgroup_symbol}), "
            f"畸变类型: {DISTORTION_TYPES.get(self.distortion_type, self.distortion_type)}, "
            f"超胞: {self.supercell}"
        )
