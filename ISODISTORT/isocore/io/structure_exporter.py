"""
结构文件导出 - CIF / POSCAR / xyz 格式，自动识别格式并导出

对应阶段六，步骤11：标准结构文件导出
"""
import re
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

    @staticmethod
    def _safe_filename(filename: str) -> str:
        """Return one cross-platform-safe basename without changing file data."""
        # Match the website's Windows-download convention: illegal characters
        # are deleted (for example I4/mmm -> I4mmm and 1/2 -> 12).
        safe = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "", str(filename)).rstrip(" .")
        if not safe:
            safe = "structure"
        reserved = {
            "CON", "PRN", "AUX", "NUL",
            *(f"COM{i}" for i in range(1, 10)),
            *(f"LPT{i}" for i in range(1, 10)),
        }
        if safe.upper() in reserved:
            safe = f"_{safe}"
        return safe[:180].rstrip(" .") or "structure"

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
        from .isodistort_cif import render_isodistort_cif  # noqa: PLC0415

        text = render_isodistort_cif(structure)
        path = self.output_dir / f"{self._safe_filename(filename)}.cif"
        path.write_text(text, encoding="utf-8", newline="\n")
        return path

    def to_poscar(self, structure: Structure, filename: str,
                comment: str = "") -> Path:
        """导出为 VASP POSCAR 格式

        """
        poscar = Poscar(structure, comment=comment)
        path = self.output_dir / f"{self._safe_filename(filename)}.vasp"
        poscar.write_file(str(path))
        return path

    def to_xyz(self, structure: Structure, filename: str) -> Path:
        """导出为 xyz 格式

        """
        xyz = XYZ(structure)
        path = self.output_dir / f"{self._safe_filename(filename)}.xyz"
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
