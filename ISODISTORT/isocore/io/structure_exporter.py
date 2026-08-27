"""
结构文件导出 - CIF / POSCAR / xyz 格式，自动识别格式并导出

对应阶段六，步骤11：标准结构文件导出
"""
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.xyz import XYZ

from ..utils import get_config


class StructureExporter:
    """晶体结构导出器"""

    def __init__(self, output_dir: str | Path | None = None):
        """初始化结构导出器

        Args:
            output_dir: 输出目录；None 时使用配置中的 output_dir
        """
        cfg = get_config()
        self.output_dir = Path(output_dir) if output_dir else cfg.output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def to_cif(self, structure: Structure, filename: str,
               symprec: float | None = None) -> Path:
        """
        导出为 CIF 格式

        Args:
            structure: 晶体结构
            filename: 文件名（不含后缀）
            symprec: 对称性精度，为 None 则不做对称化

        Returns:
            Path: 输出文件路径

        """
        from .isodistort_cif import render_isodistort_cif

        text = render_isodistort_cif(structure)
        path = self.output_dir / f"{filename}.cif"
        path.write_text(text, encoding="utf-8", newline="\n")
        return path

    def to_poscar(self, structure: Structure, filename: str,
                comment: str = "") -> Path:
        """导出为 VASP POSCAR 格式

        """
        poscar = Poscar(structure, comment=comment)
        path = self.output_dir / f"{filename}.vasp"
        poscar.write_file(str(path))
        return path

    def to_xyz(self, structure: Structure, filename: str) -> Path:
        """导出为 xyz 格式

        """
        xyz = XYZ(structure)
        path = self.output_dir / f"{filename}.xyz"
        xyz.write_file(str(path))
        return path

    def auto_export(self, structure: Structure, filename: str,
                    formats: list | None = None) -> list:
        """
        批量导出多种格式

        Args:
            formats: 格式列表，如 ["cif", "poscar", "xyz"]

        Returns:
            list of Path: 所有输出文件路径

        """
        formats = formats or ["cif"]
        paths = []
        for fmt in formats:
            fmt_lower = fmt.lower()
            if fmt_lower == "cif":
                paths.append(self.to_cif(structure, filename))
            elif fmt_lower in ("poscar", "vasp"):
                paths.append(self.to_poscar(structure, filename))
            elif fmt_lower == "xyz":
                paths.append(self.to_xyz(structure, filename))
        return paths
