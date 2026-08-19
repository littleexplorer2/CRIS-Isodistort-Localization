"""
畴变体生成器 - 基于 iso（SHOW DOMAIN）获取各畴的对称信息

对应阶段五，步骤10：畴变体生成与切换

官网 “Domains” 输出为畴列表：畴总数 = 子群在母相中的指数，
每行包含畴编号、生成元、子群空间群、基矢与原点。
本模块直接复用 iso 的真实计算结果。
"""

from ..backend import DomainInfo, IsoWrapper, SubgroupInfo
from .phase_path import PhasePath


class DomainGenerator:
    """
    畴变体生成器

    功能：
    1. 获取相变对应的全部畴变体信息（编号、生成元、基矢、原点）
    2. 畴总数 = 子群指数（与官网一致）
    """

    def __init__(self, iso_wrapper: IsoWrapper | None = None) -> None:
        self.iso = iso_wrapper or IsoWrapper()

    def get_domains(self, path: PhasePath, subgroup: SubgroupInfo) -> list[DomainInfo]:
        """
        获取相变路径的全部畴变体信息。

        Args:
            path: 相变路径（提供母相空间群号）
            subgroup: 目标子群（提供 k 点/IR/OPD）

        Returns:
            List[DomainInfo]
        """
        return self.iso.get_domains(path.parent_sg_number, subgroup)

    def get_domain_count(self, path: PhasePath, subgroup: SubgroupInfo) -> int:
        """获取畴变体总数（= 子群指数，与官网一致）。"""
        return len(self.get_domains(path, subgroup))

    def generate_domains(self, path: PhasePath,
                         subgroup: SubgroupInfo) -> list[DomainInfo]:
        """
        生成全部畴变体描述（与官网 Domains 输出一致）。

        Args:
            path: 相变路径
            subgroup: 目标子群

        Returns:
            List[DomainInfo]
        """
        return self.get_domains(path, subgroup)
