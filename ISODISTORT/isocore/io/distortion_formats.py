"""
Distortion Page 批量导出：CIF / Save interactive distortion / Complete modes details / TOPAS.STR

对齐官网第 6 页（webpage_info/6. ISODISTORT_ distortion.html）与手册
https://landau3.byu.edu/isodistorthelp.php#modeparams 中的导出选项。

命名约定（对齐官网下载）：
    文件夹  Method 1: 完整 OPD 行（Windows 删除 ``/``：``I4mmm``、``1/2``→``12``）；
            Method 2/3: ``LD1 C1``
    文件    ``subgroup.cif`` / ``data.isoviz`` /
            ``Complete modes details.txt``（官网为 HTML；本地用 .txt）/
            ``topas.str``
"""
from __future__ import annotations

import io
import re
import zipfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

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
    """去掉 Windows 非法文件名字符，压缩连续空白。

    官网 Windows 下载会直接删除 ``/``（``I4/mmm``→``I4mmm``，``1/2``→``12``），
    而不是替换成空格；其余非法字符仍替换为空格。
    """
    cleaned = str(text or "").replace("/", "")
    cleaned = _WINDOWS_BAD.sub(" ", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned).strip(" .")
    return cleaned or fallback


def subgroup_label(subgroup: SubgroupInfo) -> str:
    """子群显示名：``IR OPD``，例如 ``LD1 C1``。"""
    ir = safe_name(subgroup.irrep_label or "", "IR")
    opd = safe_name(subgroup.opd_symbol or "", "OPD")
    return f"{ir} {opd}"


def folder_label_for_subgroup(subgroup: SubgroupInfo, *, use_opd_line: bool = False) -> str:
    """ZIP / 导出文件夹名：Method 1 用完整 OPD 行，Method 2/3 用 ``IR OPD``。"""
    if use_opd_line:
        try:
            line = subgroup.opd_line()
        except (TypeError, ValueError, AttributeError):
            line = ""
        if line:
            return safe_name(line, subgroup_label(subgroup))
    return subgroup_label(subgroup)


def unique_folder_name(
    subgroup: SubgroupInfo,
    used: set[str],
    *,
    use_opd_line: bool = False,
) -> str:
    """保证文件夹名在一次导出内唯一（同名时追加空间群 / 序号）。"""
    base = folder_label_for_subgroup(subgroup, use_opd_line=use_opd_line)
    name = base
    if name in used:
        sg = safe_name(subgroup.space_group_symbol or str(subgroup.space_group_number))
        name = f"{base} {sg}".strip()
    if name in used:
        name = f"{base} #{subgroup.index}"
    used.add(name)
    return name


def format_filename(folder_label: str, fmt: str) -> str:
    """官网同款文件名（``folder_label`` 仅保留调用兼容，不写入文件名）。"""
    _ = folder_label
    if fmt == FORMAT_CIF:
        return "subgroup.cif"
    if fmt == FORMAT_ISOVIZ:
        return "data.isoviz"
    if fmt == FORMAT_MODES:
        return "Complete modes details.txt"
    if fmt == FORMAT_TOPAS:
        return "topas.str"
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
    parent_wyckoff_sites: list | None = None
    # Pre-formatted parent Wyckoff comment lines (CIF order/labels when available).
    parent_wyckoff_lines: list[str] | None = None
    distortion_types: list[str] = field(default_factory=list)


def render_cif(structure: Structure, spec: SubgroupExportSpec | None = None) -> str:
    """ISODISTORT-style CIF (subgroup setting when ``spec`` is given)."""
    from .isodistort_cif import render_isodistort_cif

    return render_isodistort_cif(structure, spec)


def render_isoviz(spec: SubgroupExportSpec) -> str:
    """Save interactive distortion：官网 ``!tag`` ISOVIZ 数据文件布局。"""
    from .isodistort_isoviz import render_isodistort_isoviz

    return render_isodistort_isoviz(spec)


def _parent_primitive_volume(parent: Structure | None) -> float:
    if parent is None:
        return 1.0
    try:
        return float(parent.get_primitive_structure().lattice.volume)
    except (ValueError, TypeError, np.linalg.LinAlgError):
        return float(parent.lattice.volume)


def render_complete_modes(spec: SubgroupExportSpec) -> str:
    """Write Complete modes details as a self-contained ``.txt``.

    Acceptance (see repo ``agent.md``): put **all locally computed** supercell /
    mode information into this text file inside each subgroup folder of the
    Distortion ZIP. Byte-level match to the official HTML modes page is **not**
    required.

    Content outline (handbook #modesdetails inspiration): supercell xyz table,
    mode definitions (every atom), As / Ap / dmax when displacements exist.
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


def _centering_multiplicity(symbol: str | None) -> int:
    """Conventional-cell centering order (P=1, I/C/A/B=2, F=4, R=3)."""
    letter = (symbol or "P").lstrip("0123456789 ").strip()[:1].upper()
    if letter == "F":
        return 4
    if letter == "R":
        return 3
    if letter in {"I", "C", "A", "B"}:
        return 2
    return 1


def compact_mode_label(pretty: str, *, topas_space_before_ir: bool = False) -> str:
    """Align mode-label cosmetics with official exports (``[0,0,0]``, optional space)."""
    text = (pretty or "").strip()
    if not text:
        return text

    def _tri(match: re.Match[str]) -> str:
        parts = []
        for raw in match.group(1).split(","):
            raw = raw.strip()
            try:
                val = float(raw)
            except ValueError:
                parts.append(raw)
                continue
            if abs(val - round(val)) < 1e-8:
                parts.append(str(int(round(val))))
            else:
                parts.append(f"{val:g}")
        return "[" + ",".join(parts) + "]"

    text = re.sub(r"\[([^\[\]]+)\]", _tri, text, count=1)
    if topas_space_before_ir:
        # Official TOPAS: ...[Al2:e:dsp] A1(a)
        text = re.sub(r"(\[[^\]]*:(?:dsp|rot|mag|occ)\])([A-Za-z])", r"\1 \2", text)
    return text


def cart_normalized_mode_matrix(
    arr: np.ndarray,
    lattice_matrix: np.ndarray,
    *,
    centering_mult: int = 1,
) -> tuple[np.ndarray, float]:
    """Scale unitless mode vectors so Σ‖Δr‖² over the *primitive* cell ≈ 1.

    Local iso often emits max-component=1 vectors. Official TOPAS / IsoVIZ use a
    Cartesian primitive-cell normalization (I-centering → factor √2 vs full cell).
    Returns ``(scaled_matrix, maxamp_hint)``.
    """
    mat = np.asarray(arr, dtype=float)
    if mat.ndim != 2 or mat.shape[1] != 3 or mat.size == 0:
        return mat, 1.0
    max_comp = float(np.max(np.abs(mat)))
    unit = mat / max_comp if max_comp > 1e-16 else mat
    bmat = np.asarray(lattice_matrix, dtype=float)
    cart = unit @ bmat
    ssq = float(np.sum(cart * cart))
    n_c = max(int(centering_mult), 1)
    # Attribute equal share of conventional-cell images to the primitive cell.
    ssq_prim = ssq / n_c
    scale = (1.0 / np.sqrt(ssq_prim)) if ssq_prim > 1e-30 else 1.0
    scaled = unit * scale
    maxamp = float(np.sqrt(n_c)) if n_c > 1 else 1.0
    return scaled, maxamp


def render_topas(spec: SubgroupExportSpec) -> str:
    """TOPAS.STR：官网 distortion-mode 精修输入（见手册 #topas）。"""
    from .isodistort_cif import _parent_to_child_transform
    from ..utils.parent_header import format_fixed_coord

    sg = spec.subgroup
    sc = spec.structure
    lat = sc.lattice
    transform = ""
    try:
        transform = _parent_to_child_transform(sg, np.zeros(3))
    except Exception:  # noqa: BLE001
        transform = "a,b,c;0,0,0"
    hm = (sg.space_group_symbol or str(sg.space_group_number or "")).strip()
    n_c = _centering_multiplicity(hm)
    lines = [
        "'Topas .str file generated by ISODISTORT",
        "'Remember to add the appropriate peak shape line when passing this into an input file",
        "",
        "\tstr",
        f"\t\t'{hm} ",
        f"\t\tspace_group {sg.space_group_number} "
        f"'transformPp {transform or 'a,b,c;0,0,0'}",
        f"\t\ta  {lat.a:10.5f}",
        f"\t\tb  {lat.b:10.5f}",
        f"\t\tc  {lat.c:10.5f}",
        f"\t\tal {lat.alpha:10.5f}",
        f"\t\tbe {lat.beta:10.5f}",
        f"\t\tga {lat.gamma:10.5f}",
        "\t\tscale @ 0.00001",
        "",
        "'{{{mode definitions",
    ]
    mode_items = list((spec.mode_displacements_sc or {}).items())
    unique = _unique_site_indices(sc)
    tags = [_site_tag_official(sc, i, spec) for i in unique]
    scaled_modes: list[np.ndarray] = []
    amp_bound = 1.41 if n_c > 1 else 2.00

    if not mode_items:
        lines.append("\t\t' (no displacive modes available for this subgroup)")
        if spec.note:
            lines.append(f"\t\t' note: {spec.note}")
    else:
        for n, (label, disp) in enumerate(mode_items, start=1):
            pretty = compact_mode_label(
                (spec.mode_labels or {}).get(label, label),
                topas_space_before_ir=True,
            )
            amp = float((spec.amplitudes or {}).get(label, 0.0))
            scaled, _hint = cart_normalized_mode_matrix(
                np.asarray(disp, dtype=float), lat.matrix, centering_mult=n_c
            )
            scaled_modes.append(scaled)
            lines.append(
                f"\t\tprm  !a{n:<4d} {amp:10.5f} min  -{amp_bound:.2f} max  {amp_bound:.2f} '{pretty}"
            )
    lines.append("")
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{mode-amplitude to delta transformation")

    deltas: dict[tuple[str, str], list[str]] = {}
    if mode_items and scaled_modes:
        for n, scaled in enumerate(scaled_modes, start=1):
            for tag, idx in zip(tags, unique, strict=True):
                if idx >= scaled.shape[0]:
                    continue
                vec = scaled[idx]
                for axis, comp in zip(("x", "y", "z"), vec, strict=True):
                    if abs(float(comp)) < 1e-8:
                        continue
                    sign = "+" if float(comp) >= 0 else "-"
                    term = f"{sign}  {abs(float(comp)):.5f}*a{n}"
                    deltas.setdefault((tag, axis), []).append(term)
        for (tag, axis), terms in deltas.items():
            expr = " ".join(terms)
            lines.append(f"\t\tprm  {tag}_d{axis}   = {expr};:  0.00000")
    lines.append("")
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{distorted parameters")
    for tag, idx in zip(tags, unique, strict=True):
        x, y, z = (float(c) for c in sc[idx].frac_coords)
        for axis, val in zip(("x", "y", "z"), (x, y, z), strict=True):
            dprm = f"{tag}_d{axis}"
            has_d = (tag, axis) in deltas
            fixed = format_fixed_coord(val)
            # Prefer IT specials (0, 1/2, …) like official TOPAS when no free delta.
            if has_d:
                lines.append(
                    f"\t\tprm  {tag}_{axis}    =    {val:.5f} + {dprm};:  {val:.5f}"
                )
            else:
                try:
                    float(fixed)
                    rhs = f"{float(fixed):.6f}" if "." in fixed else fixed
                except ValueError:
                    rhs = fixed
                if re.fullmatch(r"-?\d+", fixed) or "/" in fixed:
                    rhs = fixed
                else:
                    rhs = f"{val:.6f}"
                lines.append(
                    f"\t\tprm !{tag}_{axis}    = {rhs};:  {val:.5f}"
                )
    for tag, _idx in zip(tags, unique, strict=True):
        lines.append(f"\t\tprm !{tag}_occ  = 1;:  1.00000")
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{mode-dependent sites")
    for tag, idx in zip(tags, unique, strict=True):
        el = sc[idx].species_string
        el = re.sub(r"[^A-Za-z]", "", el) or "X"
        lines.append(
            f"\t\tsite {tag}    num_posns  0  "
            f"x = {tag}_x;:0    y = {tag}_y;:0    z = {tag}_z;:0    "
            f"occ {el:<5s} = {tag}_occ;:0 beq 0.0"
        )
    lines.append(
        "\t\t'site origin   num_posns  0  x  0.00000      y  0.00000      z  0.00000      occ D  0"
    )
    lines.append("'}}}")
    lines.append("")
    lines.append("'{{{difference restraints for interconnected rigid bodies")
    lines.append("'}}}")
    lines.append("")
    lines.append("")
    return "\n".join(lines)


def _site_tag_official(structure: Structure, idx: int, spec: SubgroupExportSpec) -> str:
    """Prefer parent Wyckoff stems (``Eu1_1``) when available."""
    site = structure[idx]
    label = (site.label or "").strip()
    if label and re.match(r"^[A-Za-z]+\d*", label):
        stem = re.split(r"[_\s]", label)[0]
        return f"{stem}_1" if "_" not in label else label
    el = re.sub(r"[^A-Za-z0-9]", "", site.species_string) or "X"
    return f"{el}_{idx + 1}"


def render_format(fmt: str, spec: SubgroupExportSpec) -> bytes:
    if fmt == FORMAT_CIF:
        text = render_cif(spec.cif_structure or spec.structure, spec)
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
    wrapping: str | None = None,
) -> bytes:
    """打包为 ZIP：各子群文件夹在 ZIP 根下（官网同款）；可选 wrapping 前缀。

    只写入 ``specs`` 给出的子群文件，不读取、不混入 output_dir 中的其它文件。
    """
    wrap = safe_name(wrapping) if wrapping else ""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        used: set[str] = set()
        for spec in specs:
            folder = spec.folder_name or unique_folder_name(spec.subgroup, used)
            spec.folder_name = folder
            for fmt in formats:
                fname = format_filename(folder, fmt)
                arcname = f"{wrap}/{folder}/{fname}" if wrap else f"{folder}/{fname}"
                zf.writestr(arcname, render_format(fmt, spec))
    return buf.getvalue()
