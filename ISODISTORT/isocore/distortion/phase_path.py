"""
相变路径管理 - 子群选择、畸变类型、参数组装与校验

对应阶段二，步骤5：相变路径参数组装（纯 Python 实现，不调用 iso 二进制）。
"""
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from ..utils import PhasePathError

if TYPE_CHECKING:
    from ..backend import SubgroupInfo


# 畸变类型（对齐官网 Types of distortions 复选框：Displacive/Occupational/Strain/
# Magnetic/Rotational）。旧名 displacement/order 作为别名保留（API/脚本兼容）。
DISTORTION_TYPES = {
    "displacive": "原子位移畸变",
    "occupational": "占据率（原子有序）畸变",
    "strain": "晶格应变畸变",
    "magnetic": "磁矩畸变",
    "rotational": "转动畸变",
}

TYPE_ALIASES = {
    "displacement": "displacive",   # 旧名（官网用 Displacive）
    "order": "occupational",        # 旧名（官网用 Occupational）
}

# 默认畸变类型（对齐官网：Types 面板默认只勾选 Strain，其余默认不勾选）
DEFAULT_DISTORTION_TYPES = ["strain"]


def normalize_distortion_types(distortion_types=None) -> list[str]:
    """
    标准化畸变类型输入。允许传入单个字符串、列表或 None。

    Args:
        distortion_types: 可为 None、单个字符串或字符串列表
            （displacement/order 旧名会自动映射为 displacive/occupational）

    Returns:
        List[str]: 标准化后的畸变类型列表，默认返回 ["displacive", "strain"]
    """
    if distortion_types is None:
        result = DEFAULT_DISTORTION_TYPES.copy()
    elif isinstance(distortion_types, str):
        result = [distortion_types]
    else:
        result = list(distortion_types)

    normalized = []
    seen = set()
    for raw_type in result:
        key = str(raw_type).strip().lower()
        key = TYPE_ALIASES.get(key, key)
        if not key or key not in DISTORTION_TYPES:
            continue
        if key not in seen:
            normalized.append(key)
            seen.add(key)

    if not normalized:
        return DEFAULT_DISTORTION_TYPES.copy()
    return normalized


@dataclass
class PhasePath:
    """
    相变路径参数对象

    描述从母相到子相的完整相变路径，包含：
    - 母相空间群
    - 目标子群（序号 + 空间群信息）
    - 畸变类型
    - 原点选择
    - 超胞参数（子群基矢）
    """
    parent_sg_number: int
    subgroup_index: int
    subgroup_sg_number: int = 0
    subgroup_symbol: str = ""
    distortion_type: str = "displacive"
    origin_choice: int = 1
    supercell: list[int] = field(default_factory=lambda: [1, 1, 1])
    basis_vectors: list[list[float]] = field(default_factory=list)  # 子群超胞基矢
    selected_irreps: list[str] = field(default_factory=list)  # 指定不可约表示
    k_vector: list[float] = field(default_factory=list)  # k 点坐标（母相倒格分数单位）

    def validate(self) -> bool:
        """参数合法性校验"""
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

    def supercell_basis(self) -> list | list[list[float]]:
        """
        返回用于生成畸变结构的超胞规格。

        优先返回子群基矢（3x3 矩阵，母相格单位）；无基矢时返回 [1,1,1]。
        """
        if self.basis_vectors and len(self.basis_vectors) == 3:
            return [list(row) for row in self.basis_vectors]
        return [1, 1, 1]

    @classmethod
    def from_subgroup(cls, parent_sg: int, subgroup: "SubgroupInfo",
                      distortion_type: str = "displacive") -> "PhasePath":
        """从子群信息快速构建相变路径（distortion_type 自动标准化）。"""
        return cls(
            parent_sg_number=parent_sg,
            subgroup_index=subgroup.index,
            subgroup_sg_number=subgroup.space_group_number,
            subgroup_symbol=subgroup.space_group_symbol,
            distortion_type=normalize_distortion_types(distortion_type)[0]
            if normalize_distortion_types(distortion_type) else "displacive",
            basis_vectors=[list(row) for row in subgroup.basis_vectors]
            if subgroup.basis_vectors else [],
        )

    def describe(self) -> str:
        """人类可读的路径描述"""
        return (
            f"空间群 #{self.parent_sg_number} → #{self.subgroup_sg_number} "
            f"({self.subgroup_symbol}), "
            f"畸变类型: {DISTORTION_TYPES.get(self.distortion_type, self.distortion_type)}, "
            f"超胞: {self.supercell}"
        )
