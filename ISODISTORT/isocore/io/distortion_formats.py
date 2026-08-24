"""
Distortion Page 批量导出：CIF / Save interactive distortion / Complete modes details / TOPAS.STR

对齐官网第 6 页（webpage_info/6. ISODISTORT_ distortion.html）与手册
https://landau3.byu.edu/isodistorthelp.php#modeparams 中的导出选项。

命名约定（子群、文件格式）：
    文件夹  ``LD1 C1``
    文件    ``LD1 C1 CIF.cif``
            ``LD1 C1 Save interactive distortion.isoviz``
            ``LD1 C1 Complete modes details.txt``
            ``LD1 C1 TOPAS.STR``
"""
from __future__ import annotations

import io
import re
import zipfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter

from ..backend import SubgroupInfo

# 官网第 6 页 origintype 与本地键的对应
FORMAT_CIF = "cif"
FORMAT_ISOVIZ = "isoviz"
FORMAT_MODES = "modes"
FORMAT_TOPAS = "topas"

SUPPORTED_FORMATS = (FORMAT_CIF, FORMAT_ISOVIZ, FORMAT_MODES, FORMAT_TOPAS)

_FORMAT_ALIASES = {
    "cif": FORMAT_CIF,
    "ciffile": FORMAT_CIF,
    "structurefile": FORMAT_CIF,
    "cif file": FORMAT_CIF,
    "isoviz": FORMAT_ISOVIZ,
    "isovizdistortion": FORMAT_ISOVIZ,
    "interactive": FORMAT_ISOVIZ,
    "save interactive distortion": FORMAT_ISOVIZ,
    "modes": FORMAT_MODES,
    "completemodes": FORMAT_MODES,
    "completemodesdetails": FORMAT_MODES,
    "modesdetails": FORMAT_MODES,
    "complete modes details": FORMAT_MODES,
    "topas": FORMAT_TOPAS,
    "topas.str": FORMAT_TOPAS,
    "str": FORMAT_TOPAS,
}

# 官网第 6 页 Visualization 默认参数
_DEFAULT_ATOMIC_RADIUS = 0.4
_DEFAULT_BOND_MIN = 0.0
_DEFAULT_BOND_MAX = 2.5
_DEFAULT_APPLET_WIDTH = 1024

_WINDOWS_BAD = re.compile(r'[<>:"/\\|?*]')


def parse_export_formats(raw: str | Sequence[str] | None) -> list[str]:
    """把查询参数 / 列表规范化为 SUPPORTED_FORMATS 中的键（去重、保序）。"""
    if raw is None:
        return [FORMAT_CIF]
    if isinstance(raw, str):
        tokens = [p.strip() for p in raw.replace(";", ",").split(",") if p.strip()]
    else:
        tokens = [str(p).strip() for p in raw if str(p).strip()]
    out: list[str] = []
    seen: set[str] = set()
    for tok in tokens:
        key = _FORMAT_ALIASES.get(tok.lower())
        if key is None:
            raise ValueError(
                f"未知导出格式 {tok!r}；可选: CIF file / Save interactive distortion / "
                "Complete modes details / TOPAS.STR"
            )
        if key not in seen:
            out.append(key)
            seen.add(key)
    if not out:
        raise ValueError("未选择任何导出格式")
    return out


def parse_export_method(raw: str | int | None) -> int:
    """解析批量导出的 Method 来源；只允许单独选择 1 / 2 / 3。

    Method 4 是畸变结构分解，不产生可按子群打包的候选列表。
    """
    if raw is None or raw == "":
        return 2
    if isinstance(raw, int):
        text = str(raw)
    else:
        text = str(raw).strip().lower()
    if any(sep in text for sep in (",", "+", ";", "|", " ")):
        raise ValueError(
            "只能选择一个 Method 导出，不能多选 / select exactly one Method"
        )
    aliases = {
        "1": 1, "method1": 1, "m1": 1,
        "2": 2, "method2": 2, "m2": 2,
        "3": 3, "method3": 3, "m3": 3,
    }
    if text not in aliases:
        raise ValueError(
            f"未知 Method {raw!r}；请选择 Method 1、2 或 3 "
            "(Method 4 无子群列表可批量导出)"
        )
    return aliases[text]


def safe_name(text: str, fallback: str = "subgroup") -> str:
    """去掉 Windows 非法文件名字符，压缩连续空白。"""
    cleaned = _WINDOWS_BAD.sub(" ", str(text or ""))
    cleaned = re.sub(r"\s+", " ", cleaned).strip(" .")
    return cleaned or fallback


def subgroup_label(subgroup: SubgroupInfo) -> str:
    """子群显示名：``IR OPD``，例如 ``LD1 C1``。"""
    ir = safe_name(subgroup.irrep_label or "", "IR")
    opd = safe_name(subgroup.opd_symbol or "", "OPD")
    return f"{ir} {opd}"


def unique_folder_name(subgroup: SubgroupInfo, used: set[str]) -> str:
    """保证文件夹名在一次导出内唯一（同 IR+OPD 时追加空间群 / 序号）。"""
    base = subgroup_label(subgroup)
    name = base
    if name in used:
        sg = safe_name(subgroup.space_group_symbol or str(subgroup.space_group_number))
        name = f"{base} {sg}".strip()
    if name in used:
        name = f"{base} #{subgroup.index}"
    used.add(name)
    return name


def format_filename(folder_label: str, fmt: str) -> str:
    """按「子群、文件格式」生成文件名。"""
    if fmt == FORMAT_CIF:
        return f"{folder_label} CIF.cif"
    if fmt == FORMAT_ISOVIZ:
        return f"{folder_label} Save interactive distortion.isoviz"
    if fmt == FORMAT_MODES:
        return f"{folder_label} Complete modes details.txt"
    if fmt == FORMAT_TOPAS:
        return f"{folder_label} TOPAS.STR"
    raise ValueError(f"未知格式: {fmt}")


def _fmt_vec(values: Sequence[float], ndigits: int = 5) -> str:
    return "(" + ",".join(f"{float(v):.{ndigits}g}" for v in values) + ")"


def _fmt_basis(basis: Sequence[Sequence[float]]) -> str:
    rows = [_fmt_vec(row) for row in basis]
    return "{" + ",".join(rows) + "}"


@dataclass
class SubgroupExportSpec:
    """单个子群一次导出所需的全部数据。"""

    subgroup: SubgroupInfo
    structure: Structure
    parent_structure: Structure | None = None
    parent_sg: int = 0
    parent_symbol: str = ""
    mode_displacements_sc: dict[str, np.ndarray] | None = None
    mode_labels: dict[str, str] | None = None
    amplitudes: dict[str, float] | None = None
    note: str = ""
    folder_name: str = ""
    cif_structure: Structure | None = None  # Generate 后 CIF 用畸变结构；其它格式用未畸变超胞


def render_cif(structure: Structure) -> str:
    """CIF 文本（不做对称化，与 StructureExporter.to_cif(symprec=None) 一致）。"""
    return str(CifWriter(structure, symprec=None))


def render_isoviz(spec: SubgroupExportSpec) -> str:
    """Save interactive distortion：官网 ISOVIZ ascii 数据文件的本地可读实现。

    内容覆盖手册所述：母相/子群、超胞原子、各对称模式位移、可视化参数
    （atomic radius / bond length / window / viewing range）。
    """
    sg = spec.subgroup
    parent = spec.parent_structure
    lines = [
        "# Save interactive distortion (ISOVIZ ascii data file)",
        "# Official option: origintype=isovizdistortion",
        "# See https://landau3.byu.edu/isodistorthelp.php#savedist",
        "origintype=isovizdistortion",
        f"spacegroup={spec.parent_sg} {spec.parent_symbol}".rstrip(),
        f"subgroupsym={sg.space_group_number}",
        f"subgroupsymbol={sg.space_group_symbol}",
        f"irrep={sg.irrep_label}",
        f"orderparam={sg.opd_symbol}",
        f"kpoint={sg.k_point_label}",
        f"kparameters={' '.join(str(p) for p in sg.k_parameters)}",
        f"basis={_fmt_basis(sg.basis_vectors) if sg.basis_vectors else ''}",
        f"origin={_fmt_vec(sg.origin) if sg.origin else '(0,0,0)'}",
        f"s={sg.size}",
        f"i={sg.subgroup_index}",
        f"atomicradius={_DEFAULT_ATOMIC_RADIUS}",
        f"bondlengthmin={_DEFAULT_BOND_MIN:.2f}",
        f"bondlength={_DEFAULT_BOND_MAX:.2f}",
        f"appletwidth={_DEFAULT_APPLET_WIDTH}",
        "supercellxmin=0.000",
        "supercellxmax=1.000",
        "supercellymin=0.000",
        "supercellymax=1.000",
        "supercellzmin=0.000",
        "supercellzmax=1.000",
    ]
    if spec.note:
        lines.append(f"note={spec.note}")
    lat = spec.structure.lattice
    lines.append(
        "lattparam="
        f"a={lat.a:.5f}, b={lat.b:.5f}, c={lat.c:.5f}, "
        f"alpha={lat.alpha:.5f}, beta={lat.beta:.5f}, gamma={lat.gamma:.5f}"
    )
    if parent is not None:
        plat = parent.lattice
        lines.append(
            "parentlattparam="
            f"a={plat.a:.5f}, b={plat.b:.5f}, c={plat.c:.5f}, "
            f"alpha={plat.alpha:.5f}, beta={plat.beta:.5f}, gamma={plat.gamma:.5f}"
        )
        lines.append("BEGIN PARENT_ATOMS")
        for i, site in enumerate(parent, start=1):
            x, y, z = (float(c) for c in site.frac_coords)
            lines.append(
                f"{i:4d}  {site.species_string:<4s}  {x:10.6f} {y:10.6f} {z:10.6f}"
            )
        lines.append("END PARENT_ATOMS")
    lines.append("BEGIN ATOMS")
    for i, site in enumerate(spec.structure, start=1):
        x, y, z = (float(c) for c in site.frac_coords)
        try:
            occ = float(site.species.num_atoms)
        except (TypeError, ValueError, AttributeError):
            occ = 1.0
        lines.append(
            f"{i:4d}  {site.species_string:<6s}  {x:10.6f} {y:10.6f} {z:10.6f}  occ={occ:.4f}"
        )
    lines.append("END ATOMS")
    lines.append("BEGIN MODES")
    if spec.mode_displacements_sc:
        for label, disp in spec.mode_displacements_sc.items():
            pretty = (spec.mode_labels or {}).get(label, label)
            amp = float((spec.amplitudes or {}).get(label, 0.0))
            lines.append(f"MODE {label}  {pretty}  amplitude={amp:.6f}")
            arr = np.asarray(disp, dtype=float)
            for j, vec in enumerate(arr, start=1):
                lines.append(
                    f"  atom {j:4d}  {float(vec[0]):10.6f} {float(vec[1]):10.6f} {float(vec[2]):10.6f}"
                )
    else:
        lines.append("# (no displacive modes available for this subgroup)")
    lines.append("END MODES")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _parent_primitive_volume(parent: Structure | None) -> float:
    if parent is None:
        return 1.0
    try:
        return float(parent.get_primitive_structure().lattice.volume)
    except (ValueError, TypeError, np.linalg.LinAlgError):
        return float(parent.lattice.volume)


def render_complete_modes(spec: SubgroupExportSpec) -> str:
    """Complete modes details：列出超胞内每个原子（不仅对称独立原子）。

    手册 #modesdetails：超胞 xyz 表、模式定义（全原子）、As / Ap / dmax。
    """
    sg = spec.subgroup
    sc = spec.structure
    lines = [
        "Complete modes details",
        "Official option: origintype=completemodesdetails",
        "See https://landau3.byu.edu/isodistorthelp.php#modesdetails",
        "",
        f"Parent space group: {spec.parent_sg} {spec.parent_symbol}".rstrip(),
        f"Subgroup: {subgroup_label(sg)}   "
        f"{sg.space_group_number} {sg.space_group_symbol}   "
        f"k={sg.k_point_label}  IR={sg.irrep_label}  OPD={sg.opd_symbol}",
        f"Basis = {_fmt_basis(sg.basis_vectors) if sg.basis_vectors else '(identity)'}  "
        f"origin = {_fmt_vec(sg.origin) if sg.origin else '(0,0,0)'}  "
        f"s={sg.size}  i={sg.subgroup_index}",
        "",
    ]
    if spec.note:
        lines += [f"Note: {spec.note}", ""]

    lat = sc.lattice
    lines.append("Lattice parameters of the supercell:")
    lines.append(
        f"  a={lat.a:.5f}  b={lat.b:.5f}  c={lat.c:.5f}  "
        f"alpha={lat.alpha:.5f}  beta={lat.beta:.5f}  gamma={lat.gamma:.5f}"
    )
    lines.append("")
    lines.append("Superstructure in the traditional atomic-xyz-coordinate basis")
    lines.append(f"{'atom':>6s} {'el':<4s} {'x':>10s} {'y':>10s} {'z':>10s}")
    for i, site in enumerate(sc, start=1):
        x, y, z = (float(c) for c in site.frac_coords)
        lines.append(
            f"{i:6d} {site.species_string:<4s} {x:10.6f} {y:10.6f} {z:10.6f}"
        )

    vp = _parent_primitive_volume(spec.parent_structure)
    vs = float(lat.volume)
    # 原胞体积比：超胞可能含多个原胞；用体积比近似 Vp/Vs
    scale_ap = float(np.sqrt(max(vp, 1e-30) / max(vs, 1e-30)))

    lines.append("")
    lines.append("Mode definitions (every atom in the unit cell)")
    lines.append(
        "Mode vectors are given in unitless superlattice coordinates; "
        "the largest component of each mode is scaled to 1.0. "
        "normfactor makes the sum of squares of Cartesian changes over all "
        "atoms equal to 1.0 (Angstrom for displacive modes)."
    )
    if not spec.mode_displacements_sc:
        lines.append("(no displacive modes available for this subgroup)")
        return "\n".join(lines) + "\n"

    bmat = np.asarray(lat.matrix, dtype=float)
    amp_rows: list[str] = []
    for label, disp in spec.mode_displacements_sc.items():
        pretty = (spec.mode_labels or {}).get(label, label)
        arr = np.asarray(disp, dtype=float)
        if arr.size == 0:
            continue
        max_comp = float(np.max(np.abs(arr)))
        unit = arr / max_comp if max_comp > 1e-16 else arr
        cart = unit @ bmat
        ssq = float(np.sum(cart * cart))
        norm = (1.0 / np.sqrt(ssq)) if ssq > 1e-30 else 0.0
        as_amp = float((spec.amplitudes or {}).get(label, 0.0))
        ap_amp = as_amp * scale_ap
        dmax = float(np.max(np.linalg.norm(cart * as_amp, axis=1))) if as_amp else \
            float(np.max(np.linalg.norm(cart, axis=1)))
        lines.append("")
        lines.append(f"Mode {label}  {pretty}")
        lines.append(f"  normfactor = {norm:.6g} Angstrom^-1")
        for j, vec in enumerate(unit, start=1):
            if np.max(np.abs(vec)) < 1e-10:
                continue
            lines.append(
                f"  atom {j:4d}  {sc[j - 1].species_string:<4s}  "
                f"({float(vec[0]): .6f}, {float(vec[1]): .6f}, {float(vec[2]): .6f})"
            )
        amp_rows.append(
            f"  {label:<16s}  As={as_amp:10.6f}  Ap={ap_amp:10.6f}  "
            f"dmax={dmax:10.6f} Angstrom"
        )

    lines.append("")
    lines.append("Mode amplitudes (As = supercell-normalized, Ap = parent-cell-normalized)")
    lines.append("Ap = As * sqrt(Vp/Vs); dmax is the largest Cartesian atomic displacement.")
    lines.extend(amp_rows or ["  (none)"])
    return "\n".join(lines) + "\n"


def _unique_site_indices(structure: Structure) -> list[int]:
    """对称独立位点在超胞中的代表下标；失败时退回全部原子（P1）。"""
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: PLC0415
        sga = SpacegroupAnalyzer(structure, symprec=1e-3)
        eq = sga.get_symmetrized_structure().equivalent_indices
        return [group[0] for group in eq if group]
    except (ValueError, TypeError, np.linalg.LinAlgError, AttributeError):
        return list(range(len(structure)))


def _site_tag(structure: Structure, idx: int) -> str:
    """TOPAS 位点名：Eu_1 / Al_2。"""
    el = structure[idx].species_string
    el = re.sub(r"[^A-Za-z0-9]", "", el) or "X"
    return f"{el}_{idx + 1}"


def render_topas(spec: SubgroupExportSpec) -> str:
    """TOPAS.STR：官网 distortion-mode 精修输入（见手册 #topas）。"""
    sg = spec.subgroup
    sc = spec.structure
    lat = sc.lattice
    lines = [
        "'Topas .str file generated by ISODISTORT",
        "'Remember to add the appropriate peak shape line when passing this into an input file",
        "",
        " str",
        f" '{sg.space_group_symbol or sg.space_group_number}",
        f" space_group {sg.space_group_number}",
        f" a {lat.a:.5f}",
        f" b {lat.b:.5f}",
        f" c {lat.c:.5f}",
        f" al {lat.alpha:.5f}",
        f" be {lat.beta:.5f}",
        f" ga {lat.gamma:.5f}",
        " scale @ 0.00001",
        "",
        "'{{{mode definitions",
    ]
    mode_items = list((spec.mode_displacements_sc or {}).items())
    unique = _unique_site_indices(sc)
    tags = [_site_tag(sc, i) for i in unique]

    if not mode_items:
        lines.append("' (no displacive modes available for this subgroup)")
        if spec.note:
            lines.append(f"' note: {spec.note}")
    else:
        for n, (label, _disp) in enumerate(mode_items, start=1):
            pretty = (spec.mode_labels or {}).get(label, label)
            amp = float((spec.amplitudes or {}).get(label, 0.0))
            lines.append(
                f" prm !a{n} {amp:.5f} min -2.00 max 2.00 "
                f"'{sg.k_point_label}{sg.irrep_label}({sg.opd_symbol}) {pretty}"
            )
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{mode-amplitude to delta transformation")

    # 每个独立位点、每个模式 -> 超胞分数位移，写 dx/dy/dz 线性组合
    deltas: dict[tuple[str, str], list[str]] = {}  # (tag, axis) -> ["+ 0.12*a1", ...]
    if mode_items:
        for n, (_label, disp) in enumerate(mode_items, start=1):
            arr = np.asarray(disp, dtype=float)
            for tag, idx in zip(tags, unique, strict=True):
                if idx >= arr.shape[0]:
                    continue
                vec = arr[idx]
                for axis, comp in zip(("x", "y", "z"), vec, strict=True):
                    if abs(float(comp)) < 1e-8:
                        continue
                    sign = "+" if float(comp) >= 0 else "-"
                    term = f"{sign} {abs(float(comp)):.5f}*a{n}"
                    deltas.setdefault((tag, axis), []).append(term)
        for (tag, axis), terms in deltas.items():
            expr = " ".join(terms)
            lines.append(f" prm {tag}_d{axis} = {expr};: 0.00000")
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{distorted parameters")
    for tag, idx in zip(tags, unique, strict=True):
        x, y, z = (float(c) for c in sc[idx].frac_coords)
        for axis, val in zip(("x", "y", "z"), (x, y, z), strict=True):
            dprm = f"{tag}_d{axis}"
            has_d = (tag, axis) in deltas
            if has_d:
                lines.append(
                    f" prm {tag}_{axis} = {val:.6f} + {dprm};: {val:.5f}"
                )
            else:
                lines.append(
                    f" prm !{tag}_{axis} = {val:.6f};: {val:.5f}"
                )
        lines.append(f" prm !{tag}_occ = 1;: 1.00000")
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{mode-dependent sites")
    for tag, idx in zip(tags, unique, strict=True):
        el = sc[idx].species_string
        el = re.sub(r"[^A-Za-z]", "", el) or "X"
        lines.append(
            f" site {tag} num_posns 0 "
            f"x = {tag}_x;:0 y = {tag}_y;:0 z = {tag}_z;:0 "
            f"occ {el} = {tag}_occ;:0 beq 0.0"
        )
    lines.append(" site origin num_posns 0 x 0.00000 y 0.00000 z 0.00000 occ D 0")
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{difference restraints for interconnected rigid bodies")
    lines.append("'}}}")
    lines.append("")
    return "\n".join(lines)


def render_format(fmt: str, spec: SubgroupExportSpec) -> bytes:
    if fmt == FORMAT_CIF:
        text = render_cif(spec.cif_structure or spec.structure)
    elif fmt == FORMAT_ISOVIZ:
        text = render_isoviz(spec)
    elif fmt == FORMAT_MODES:
        text = render_complete_modes(spec)
    elif fmt == FORMAT_TOPAS:
        text = render_topas(spec)
    else:
        raise ValueError(f"未知格式: {fmt}")
    return text.encode("utf-8")


def write_subgroup_files(
    dest_dir: Path,
    spec: SubgroupExportSpec,
    formats: Sequence[str],
) -> list[Path]:
    """把选定格式写入 ``dest_dir``（应为该子群的文件夹）。"""
    dest_dir.mkdir(parents=True, exist_ok=True)
    folder = spec.folder_name or subgroup_label(spec.subgroup)
    paths: list[Path] = []
    for fmt in formats:
        name = format_filename(folder, fmt)
        path = dest_dir / name
        path.write_bytes(render_format(fmt, spec))
        paths.append(path)
    return paths


def build_export_zip(
    specs: Iterable[SubgroupExportSpec],
    formats: Sequence[str],
    wrapping: str = "isodistort_outputs",
) -> bytes:
    """打包为 ZIP：解压得到 wrapping/ 下各子群文件夹。

    只写入 ``specs`` 给出的子群文件，不读取、不混入 output_dir 中的其它文件。
    """
    wrap = safe_name(wrapping, "isodistort_outputs")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        used: set[str] = set()
        for spec in specs:
            folder = spec.folder_name or unique_folder_name(spec.subgroup, used)
            spec.folder_name = folder
            for fmt in formats:
                fname = format_filename(folder, fmt)
                arcname = f"{wrap}/{folder}/{fname}"
                zf.writestr(arcname, render_format(fmt, spec))
    return buf.getvalue()
