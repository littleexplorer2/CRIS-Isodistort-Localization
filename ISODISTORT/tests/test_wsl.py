"""WSL-dependent golden-standard and real-binary smoke tests."""
from __future__ import annotations

import io
import re
import shutil
import subprocess
import zipfile
from html.parser import HTMLParser
from pathlib import Path

import numpy as np
import pytest
from data_dir import experiment_data_dir
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from backend.api import IsoDistort
from backend.wrappers import FindsymWrapper, IsoWrapper, SubgroupInfo
from features.input_cif import read_cif

MATCHER = StructureMatcher(ltol=1e-5, stol=1e-3, angle_tol=0.001)
COORD_TOL = 1e-5
AMP_REL_TOL = 1e-4
DATA_DIR = experiment_data_dir()


def _wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        result = subprocess.run(  # noqa: PLW1510
            ["wsl.exe", "--status"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, OSError):
        return False


pytestmark = pytest.mark.skipif(
    not _wsl_available(), reason="WSL 不可用，跳过真实二进制 / 金标准测试"
)


# --- from test_golden_standard.py ---

def _srtio3_parent(tmp_path: Path) -> tuple[IsoDistort, Structure]:
    """构造立方钙钛矿 SrTiO₃（Pm-3m #221，5 原子，a=3.905 Å）。

    位点设置对齐官网 ISODISTORT 示例（isodistortexample.php）：
    Sr 1b(1/2,1/2,1/2)、Ti 1a(0,0,0)、O 3d(1/2,0,0),(0,1/2,0),(0,0,1/2)。
    （注意：若把 O 放在 3c，iso 的 R4+ 根模式落在 3d 位点，
    本地映射器会找不到对应原子——即“空间群原点设置”陷阱。）
    """
    lattice = [[3.905, 0, 0], [0, 3.905, 0], [0, 0, 3.905]]
    parent = Structure(
        lattice,
        ["Sr", "Ti", "O", "O", "O"],
        [[0.5, 0.5, 0.5], [0, 0, 0],
         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
    )
    sg = SpacegroupAnalyzer(parent, symprec=1e-3).get_space_group_number()
    assert sg == 221, f"SrTiO₃ 母相识别为 #{sg}，应为 Pm-3m #221"

    cif = tmp_path / "srtio3.cif"
    CifWriter(parent).write_file(str(cif))

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    assert iso.symmetry_info["space_group_number"] == 221
    return iso, parent


def _first_mode_at(iso: IsoDistort, k: str, ir: str) -> tuple:
    """在 (k, IR) 下取首个可计算位移模式的子群，返回 (子群, 模式标号)。"""
    subs = iso.list_subgroups_at(k, ir)
    for sg in subs:
        try:
            iso.search_method_2(subgroup_idx=sg.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式/计算失败的子群跳过
            continue
        if iso.mode_displacements:
            return sg, next(iter(iso.mode_displacements))
    raise AssertionError(f"{k} {ir}: 无可计算位移模式的子群")


def _distorted_sg(iso: IsoDistort, label: str, amplitude: float) -> int:
    d = iso.generate_distortion(irrep_label=label, amplitude=amplitude)
    return SpacegroupAnalyzer(d, symprec=1e-3).get_space_group_number()


# ------------------------------------------------------------------
# 金标准相变路径（文献公认子群）
# ------------------------------------------------------------------

def test_srtio3_r4_plus_tilt_path(tmp_path):
    """R₄⁺ 倾转（a⁰a⁰c⁻）→ I4/mcm #140（文献公认）。

    校验链：子群枚举 → 模式计算 → 生成 → spglib 验证 = 子群 → StructureMatcher
    零振幅回退 = 母相。
    """
    iso, parent = _srtio3_parent(tmp_path)
    subgroup, label = _first_mode_at(iso, "R", "R4+")
    assert subgroup.space_group_number == 140, \
        f"R4+ 首个子群应为 I4/mcm #140，实际 {subgroup.space_group_number} {subgroup.space_group_symbol}"

    # 非零振幅：畸变结构 spglib 对称性 == 目标子群
    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 140, f"R4+ 畸变结构对称性 #{d_sg}，应为 #140"

    # 零振幅：StructureMatcher 结构等价 + 空间群回母相
    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert SpacegroupAnalyzer(d0, symprec=1e-3).get_space_group_number() == 221
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"


def test_srtio3_m3_plus_inphase_tilt_path(tmp_path):
    """M₃⁺ 同相倾转 → P4/mbm #127（文献公认）。"""
    iso, parent = _srtio3_parent(tmp_path)
    subgroup, label = _first_mode_at(iso, "M", "M3+")
    assert subgroup.space_group_number == 127, \
        f"M3+ 首个子群应为 P4/mbm #127，实际 {subgroup.space_group_number} {subgroup.space_group_symbol}"

    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 127, f"M3+ 畸变结构对称性 #{d_sg}，应为 #127"

    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"


def test_srtio3_gm4_minus_polar_path(tmp_path):
    """Γ₄⁻ 极性位移 → P4mm #99（文献公认的铁电软模路径）。"""
    iso, parent = _srtio3_parent(tmp_path)
    subgroup, label = _first_mode_at(iso, "GM", "GM4-")
    assert subgroup.space_group_number == 99, \
        f"GM4- 首个子群应为 P4mm #99，实际 {subgroup.space_group_number} {subgroup.space_group_symbol}"

    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 99, f"GM4- 畸变结构对称性 #{d_sg}，应为 #99"

    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"


# ------------------------------------------------------------------
# 数值精度（容差符合科研标准）
# ------------------------------------------------------------------

def test_amplitude_linearity_cartesian(tmp_path):
    """振幅加倍 -> 笛卡尔位移长度线性加倍（相对误差 ≤ 1e-4）。

    用 Γ 点模式（GM4-，超胞因子 1）保证母相/畸变原子一一对应。
    """
    iso, parent = _srtio3_parent(tmp_path)
    _, label = _first_mode_at(iso, "GM", "GM4-")

    d1 = iso.generate_distortion(irrep_label=label, amplitude=0.1)
    d2 = iso.generate_distortion(irrep_label=label, amplitude=0.2)
    assert len(d1) == len(d2) == len(parent)

    cart1 = np.asarray([s.coords for s in d1])
    cart2 = np.asarray([s.coords for s in d2])
    cart0 = np.asarray([s.coords for s in parent])
    disp1, disp2 = [], []

    def _displacement_length(coords: np.ndarray,
                            base: np.ndarray,
                            lattice_matrix: np.ndarray) -> float:
        best = np.inf
        for j in range(len(coords)):
            delta = coords[j] - base
            frac = np.linalg.solve(lattice_matrix.T, delta)
            frac -= np.round(frac)
            d_cart = frac @ lattice_matrix
            best = min(best, float(np.linalg.norm(d_cart)))
        return best

    for i in range(len(parent)):
        base = cart0[i]
        disp1.append(_displacement_length(cart1, base, parent.lattice.matrix))
        disp2.append(_displacement_length(cart2, base, parent.lattice.matrix))
    for v1, v2 in zip(disp1, disp2, strict=True):
        if v1 > 1e-9:
            ratio = v2 / v1
            assert abs(ratio - 2.0) <= AMP_REL_TOL, \
                f"振幅线性失败：位移比 {ratio:.6f}（期望 2.0）"


def test_zero_amplitude_coordinate_tolerance(tmp_path):
    """零振幅：畸变结构原子必须与母相一一对应（周期最小镜像 ≤ 1e-5）。

    子群基矢可能是旋转幺模胞（如 GM4- → P4mm），分数坐标随格子基
    变化而重排，不能按行号比对（“原子顺序问题”陷阱）；按同物种 +
    周期最小镜像距离建立原子对应后比较。
    """
    iso, parent = _srtio3_parent(tmp_path)
    _, label = _first_mode_at(iso, "GM", "GM4-")
    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert len(d0) == len(parent)
    # 结构语义等价（含原子重排/周期平移）
    assert MATCHER.fit(parent, d0), "零振幅畸变结构必须与母相结构等价"

    # 原子级对应：同物种 + 最小镜像距离（周期约化）
    pc = np.asarray(parent.frac_coords, dtype=float)
    dc = np.asarray(d0.frac_coords, dtype=float)
    max_disp = 0.0
    for i in range(len(parent)):
        best = float("inf")
        for j in range(len(d0)):
            if parent[i].species_string != d0[j].species_string:
                continue
            delta = dc[j] - pc[i]
            delta -= np.round(delta)
            best = min(best, float(np.linalg.norm(delta)))
        max_disp = max(max_disp, best)
    assert max_disp <= COORD_TOL, \
        f"零振幅原子级最大位移 {max_disp:.2e} 超过 {COORD_TOL}"


# ------------------------------------------------------------------
# 区边界 k 点相位（带心母相的副本反号）回归
# ------------------------------------------------------------------

def test_im3m_h4_minus_zone_boundary_phase(tmp_path):
    """Im-3m #229 的 H4-（k=(1,1,1) 惯用坐标）区边界模式。

    回归保护：bcc 惯用胞 2 原子（2a 轨道）相差原始格点平移
    (1/2,1/2,1/2)。BUSH 对该模式输出两个**相位反号**的代表点
    `(0,0,0)→(1,0,1)` 与 `(-1/2,1/2,1/2)→(-1,0,-1)`，二者互为
    “模格点等价”。若映射器按“首个代表胜出”分配，两个原子会得到
    相同位移（刚性平移，畸变退化为母相 #229，对称性校验失败）。
    修复后按“周期等价（mod 1）”分配：原子 0→(1,0,1)、原子 1→(-1,0,-1)，
    畸变结构正确降为子群对称性（H4- 首个子群 P4_2/nmc #129）。
    """
    lattice = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]]
    parent = Structure(
        lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    sg = SpacegroupAnalyzer(parent, symprec=1e-3).get_space_group_number()
    assert sg == 229, f"Im-3m 母相识别为 #{sg}，应为 #229"

    cif = tmp_path / "im3m.cif"
    CifWriter(parent).write_file(str(cif))
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)

    subgroup, label = _first_mode_at(iso, "H", "H4-")
    assert subgroup.space_group_number == 129, \
        f"H4- 首个子群应为 P4_2/nmc #129，实际 " \
        f"{subgroup.space_group_number} {subgroup.space_group_symbol}"

    # 关键回归：两个原子的位移必须相位反号（非均匀），否则退化为刚性平移
    disp = iso.mode_displacements[label]["displacements"]
    assert not np.allclose(disp, disp[0], atol=1e-6), \
        "H4- 带心副本位移被映射为均匀（刚性平移）——相位丢失回归"
    assert np.allclose(disp[0], -disp[1], atol=1e-6), \
        "H4- 两个 bcc 原子的位移应为相位反号（±(1,1,0)/√2 方向）"

    d_sg = _distorted_sg(iso, label, 0.1)
    assert d_sg == 129, f"H4- 畸变结构对称性 #{d_sg}，应为 #129"

    # 零振幅回退母相
    d0 = iso.generate_distortion(irrep_label=label, amplitude=0.0)
    assert SpacegroupAnalyzer(d0, symprec=1e-3).get_space_group_number() == 229


# ------------------------------------------------------------------
# 官网参考 CIF（LD1 零振幅框架文件）金标准比对
# ------------------------------------------------------------------

def test_official_ld1_reference_structure(tmp_path):
    """官网 LD1_C1_subgroup.cif：零振幅框架文件，读回 I4/mmm #139。

    与本地 API 对 EuAl4 的识别结果做结构语义比对（StructureMatcher 允许
    超胞倍数差异，此处校验原子数/元素/空间群信息一致）。
    """
    official = DATA_DIR / "LD1_C1_subgroup.cif"
    if not official.exists():
        pytest.skip("官网参考 CIF 不存在")

    ref = Structure.from_file(str(official))
    ref_sg = SpacegroupAnalyzer(ref, symprec=1e-3).get_space_group_number()
    # 零振幅框架文件：序参量全 0，读回应为母相 I4/mmm #139
    assert ref_sg == 139, f"官网 LD1 参考读回 #{ref_sg}，应为 I4/mmm #139"
    # 框架文件为 6 倍 c 超胞（s=12 -> 60 原子）；原子数验证
    assert len(ref) == 60, f"官网 LD1 参考应为 60 原子（6×10），实际 {len(ref)}"
    assert {s.species_string for s in ref} <= {"Eu", "Al"}

    # 本地 EuAl4 母相识别
    parent_cif = DATA_DIR / "EuAl4 Parent.cif"
    if not parent_cif.exists():
        pytest.skip("EuAl4 母相 CIF 不存在")
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(parent_cif)
    assert iso.symmetry_info["space_group_number"] == 139
    # 官网框架文件的 6×c 超胞与母相 10 原子：原子总数比例一致
    assert len(ref) == 6 * len(iso.structure), \
        "官网 LD1 参考原子数应为母相 6 倍（c 轴 6 倍超胞）"


# --- from test_real_binaries.py ---

def test_findsym_identifies_nacl():
    """findsym 识别 NaCl（F 心）应为 Fm-3m #225。"""
    fs = FindsymWrapper()
    result = fs.identify(
        lattice_params=[5.63, 5.63, 5.63, 90, 90, 90],
        atom_types=["Na", "Cl"],
        atom_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        centering="F",
    )
    assert result.space_group_number == 225
    assert result.space_group_symbol == "Fm-3m"
    assert {s["wyckoff_letter"] for s in result.wyckoff_sites} == {"a", "b"}


def test_method1_eual4_matches_official_opd_html():
    """Method 1 OPD identities and exact affine subgroups match official evidence."""
    from manual import validate_method_outputs as validator
    from manual.official_html_resolver import scan_official_html

    from backend.utils.lattice import determinant, inverse, multiply
    from features.method1.affine_embeddings import (
        affine_equivalence,
        embedding_from_identity,
        parent_affine_group,
    )

    official_root = (
        Path(__file__).resolve().parents[2]
        / "webpage_info"
        / "EuAl4 Parent.cif"
    )
    assert official_root.is_dir(), "official parent evidence directory not present"
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("EuAl4 Parent.cif not present")

    class _OpdForm(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=True)
            self.fields = {}
            self.lines = []
            self._post = False
            self._line = None

        def handle_starttag(self, tag, attrs):
            attributes = dict(attrs)
            if tag == "form":
                self._post = attributes.get("method", "").casefold() == "post"
            elif tag == "input" and self._post:
                name = attributes.get("name")
                if name:
                    self.fields.setdefault(name, []).append(attributes.get("value", ""))
                if name == "orderparam":
                    assert attributes.get("type", "").casefold() == "radio"
                    assert self._line is None, "unterminated official OPD line"
                    self._line = []
            elif tag == "br" and self._line is not None:
                self.lines.append("".join(self._line).strip())
                self._line = None

        def handle_endtag(self, tag):
            if tag == "form":
                self._post = False

        def handle_data(self, data):
            if self._line is not None:
                self._line.append(data)

    # Capture once through the shared evidence owner. Its generic roles do
    # not yet include OPD pages, so distinguish Method 1 by content, not names.
    inventory = scan_official_html(official_root, recursive=False)
    matches = []
    for page in inventory.pages:
        if page.title != "isodistort: order parameter direction":
            continue
        parsed = _OpdForm()
        parsed.feed(page.raw_html)
        parsed.close()
        if parsed.fields.get("input") == ["distort"] \
                and parsed.fields.get("origintype") == ["method1"]:
            matches.append((page, parsed))
    assert len(matches) == 1, (
        "expected one content-identified Method 1 OPD page; "
        f"matches={[str(page.path) for page, _ in matches]}"
    )
    page, parsed = matches[0]
    official = parsed.lines
    assert official and len(official) == len(parsed.fields.get("orderparam", []))
    assert len(set(official)) == len(official), "duplicate official OPD identities"

    def field(name):
        values = parsed.fields.get(name, [])
        assert len(values) == 1, f"missing or ambiguous official context field {name!r}"
        return values[0].strip()

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)

    # Verify the complete parent and active distortion scope before comparing
    # candidates. No scientific condition is inferred from the archive path.
    assert int(field("spacegroup").split()[0]) == iso.symmetry_info["space_group_number"]
    assert np.allclose(
        [float(value) for value in field("dlattparam").split()],
        iso.structure.lattice.parameters, atol=1e-8, rtol=0.0,
    )
    def normalize(text):
        return " ".join(text.split())

    count = int(field("wycount"))
    official_sites = [
        f"{field(f'wyatom{index:03d}')} {field(f'wyckoff{index:03d}')}"
        for index in range(1, count + 1)
    ]
    assert [normalize(line) for line in official_sites] == [
        normalize(line) for line in iso.parent_wyckoff_display()
    ]
    assert normalize("Default space-group preferences: " + iso.space_group_preferences()) \
        in normalize(page.raw_html)
    type_species = {}
    for index in range(1, count + 1):
        atom_type = int(field(f"wytype{index:03d}"))
        species = field(f"wyatomtype{index:03d}")
        assert type_species.setdefault(atom_type, species) == species
    assert set(type_species.values()) == set(iso.species())
    enabled = {
        name for name, values in parsed.fields.items()
        if name.startswith("include") and any(value.casefold() == "true" for value in values)
    }
    assert enabled == {"includestrain", *(
        f"includedisplacive{atom_type:03d}" for atom_type in type_species
    )}, "official query is not strain plus displacive on every parent species"

    cands = iso.search_method_1(
        distortion_types=["displacive", "strain"],
        crystal_system=None, subgroup_space_group=None,
        maximal_subgroup_only=False,
    )
    local = [it.subgroup.opd_line().strip() for it in cands]
    assert len(local) == len(official), (
        f"Method 1 count {len(local)} != parsed official {len(official)}"
    )
    assert len(set(local)) == len(local), "duplicate local OPD identities"

    def parse_identity(line):
        details = validator._parse_subgroup_details_identity(line)
        # The shared parser owns rational SG/B/O/s/i/k semantics. This thin
        # text adapter retains the three Method 1 IR/OPD/direction fields.
        head = line.split(", basis=", 1)[0].split()
        assert len(head) == 5 and re.fullmatch(r"\([^)]*\)", head[2])
        assert int(head[3]) == details.space_group_number
        assert head[4].casefold() == details.space_group_symbol
        key = (
            *head[:3], details.space_group_number, details.space_group_symbol,
            details.sublattice_index, details.group_index, details.k_active,
        )
        origin_tokens = re.findall(r"origin=(\([^)]*\))", line)
        assert len(origin_tokens) == 1, "missing or ambiguous candidate origin"
        # The shared identity parser reduces origin modulo one for matching.
        # Restore its exact raw value for the actual Seitz embedding: an
        # integer parent-origin shift need not preserve a child subgroup.
        origin = validator._parse_vectors(origin_tokens[0], 1)[0]
        return key, details, origin

    def index_lines(lines):
        indexed = {}
        for line in lines:
            key, details, origin = parse_identity(line)
            assert key not in indexed, f"duplicate complete OPD identity: {key}"
            indexed[key] = (line, details, origin)
        return indexed

    official_by_identity = index_lines(official)
    local_by_identity = index_lines(local)
    missing = set(official_by_identity) - set(local_by_identity)
    extra = set(local_by_identity) - set(official_by_identity)
    assert not missing, f"missing complete identities: {missing}"
    assert not extra, f"extra complete identities: {extra}"
    parent = parent_affine_group(
        iso.structure,
        symprec=iso.cfg.symmetry_cartesian_tolerance_angstrom,
        angle_tolerance_degrees=iso.cfg.symmetry_angle_tolerance_degrees,
    )
    validation_records = []
    for key, (official_line, official_details, official_origin) in official_by_identity.items():
        local_line, local_details, local_origin = local_by_identity[key]
        assert validator.bases_generate_same_lattice(
            local_details.basis, official_details.basis,
        ), f"non-unimodular conventional basis change: {key}"
        embeddings = [
            embedding_from_identity({
                "space_group_number": details.space_group_number,
                "basis": details.basis,
                "origin": origin,
            }, parent)
            for details, origin in (
                (official_details, official_origin), (local_details, local_origin)
            )
        ]
        official_embedding, local_embedding = embeddings
        assert official_embedding.hall_number == local_embedding.hall_number
        relationship = affine_equivalence(
            local_embedding, official_embedding, parent, allow_parent_conjugacy=False,
        )
        assert relationship == "equal", f"different exact affine subgroup: {key}"
        for embedding in embeddings:
            size = abs(determinant(multiply(embedding.lattice, inverse(parent.lattice))))
            assert size == official_details.sublattice_index
            assert size * len(parent.operations) / len(embedding.operations) == \
                official_details.group_index
        validation_records.append({
            "identity": key,
            "official_line": official_line,
            "local_line": local_line,
            "official_basis": official_details.basis,
            "local_basis": local_details.basis,
            "official_origin": official_origin,
            "local_origin": local_origin,
            "U_local_official": multiply(local_details.basis, inverse(official_details.basis)),
            "hall_number": official_embedding.hall_number,
            "seitz_counts": [len(official_embedding.operations), len(local_embedding.operations)],
            "relationship": relationship,
        })
    assert len(validation_records) == len(official_by_identity)


def test_eual4_rootless_wyckoff_secondary_modes_are_completed():
    """BUSH root modes and child-fixed rootless-orbit modes form one basis."""
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(DATA_DIR / "EuAl4 Parent.cif")
    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    for irrep, opd, expected in (
        ("GM4+", "P1", 1),
        ("P2", "C1", 5),
        ("X1-", "P1", 4),
    ):
        target = next(
            subgroup for subgroup in candidates
            if subgroup.irrep_label == irrep and subgroup.opd_symbol == opd
        )
        iso.search_method_2(
            target.index,
            distortion_type=["displacive"],
            candidates=candidates,
        )
        assert len(iso.mode_displacements) == expected
        assert any(mode.wyckoff_site == "e" for mode in iso.distortion_modes)


def test_ndnio2_rootless_oxygen_secondary_modes_are_completed():
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(DATA_DIR / "NdNiO2 own.cif")
    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    target = next(
        subgroup for subgroup in candidates
        if subgroup.irrep_label == "R1-" and subgroup.opd_symbol == "C1"
    )
    iso.search_method_2(
        target.index,
        distortion_type=["displacive"],
        candidates=candidates,
    )
    assert len(iso.mode_displacements) == 4
    oxygen = [mode.irrep_label for mode in iso.distortion_modes if mode.wyckoff_site == "f"]
    assert oxygen == ["M1+", "M2+"]


def test_4310_n1plus_4d1_complete_modes_cover_the_child_cell():
    """A multi-arm special-k route must not lose part of a rooted orbit."""
    cif = DATA_DIR / "4310_tetra.cif"
    official = (
        Path(__file__).resolve().parents[2]
        / "output_compare"
        / "4310_tetra.cif"
        / "官网"
        / "Method1"
        / "N1+_4D1_SG2"
        / "data.isoviz"
    )
    if not cif.is_file() or not official.is_file():
        pytest.skip("validated 4310 parent/official Method 1 reference unavailable")

    official_text = official.read_text(encoding="utf-8-sig")
    mode_block = official_text.split("!displacivemodelist", 1)[1]
    expected = len(re.findall(r"(?m)^\s*\d+\s+\d+\s+[-+0-9.Ee]+\s+", mode_block))
    assert expected > 0

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    target = next(
        subgroup for subgroup in candidates
        if subgroup.irrep_label == "N1+" and subgroup.opd_symbol == "4D1"
    )

    iso.search_method_2(
        target.index,
        distortion_type=["displacive"],
        candidates=candidates,
    )

    assert len(iso.mode_displacements_sc) == expected


def test_shifted_4310_parent_is_canonicalized_before_method1_and_modes():
    """An equivalent global origin shift must not change Method 1 completeness."""
    cif = DATA_DIR / "4310_tetra.cif"
    official = (
        Path(__file__).resolve().parents[2]
        / "output_compare"
        / "4310_tetra.cif"
        / "官网"
        / "Method1"
        / "N1+_4D1_SG2"
        / "data.isoviz"
    )
    if not cif.is_file() or not official.is_file():
        pytest.skip("validated 4310 parent/official Method 1 reference unavailable")

    shifted = read_cif(cif)
    shifted.translate_sites(
        range(len(shifted)),
        [0.0, 0.0, 0.1],
        frac_coords=True,
        to_unit_cell=True,
    )
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.set_structure(shifted)

    assert iso.symmetry_info["parent_orbit_standardization"]["status"] == "canonicalized"
    assert iso.symmetry_info["space_group_number"] == 139
    assert "Ni1 2a (0,0,0)" in iso.parent_wyckoff_display()

    candidates = [
        item.subgroup
        for item in iso.search_method_1(distortion_types=["strain", "displacive"])
    ]
    assert len(candidates) == 125
    target = next(
        subgroup for subgroup in candidates
        if subgroup.irrep_label == "N1+" and subgroup.opd_symbol == "4D1"
    )
    iso.search_method_2(
        target.index,
        distortion_type=["displacive"],
        candidates=candidates,
    )

    official_text = official.read_text(encoding="utf-8-sig")
    mode_block = official_text.split("!displacivemodelist", 1)[1]
    expected = len(re.findall(r"(?m)^\s*\d+\s+\d+\s+[-+0-9.Ee]+\s+", mode_block))
    assert len(iso.mode_displacements_sc) == expected


def test_primitive_centered_parent_is_canonicalized_as_a_complete_cell():
    conventional = Structure.from_spacegroup(
        225,
        [[5.64, 0, 0], [0, 5.64, 0], [0, 0, 5.64]],
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    primitive = SpacegroupAnalyzer(
        conventional, symprec=1e-3
    ).get_primitive_standard_structure()
    assert len(primitive) == 2

    iso = IsoDistort(language="en")
    standardized = iso.set_structure(primitive)

    assert iso.symmetry_info["space_group_number"] == 225
    assert iso.symmetry_info["parent_orbit_standardization"]["status"] == "canonicalized"
    assert len(standardized) == 8
    assert {
        (site["multiplicity"], site["wyckoff_letter"], site["species"])
        for site in iso.symmetry_info["wyckoff_sites"]
    } == {(4, "a", "Na"), (4, "b", "Cl")}


def test_iso_kpoints_and_subgroups():
    """iso 枚举 SG 225 的 k 点与 GM5- 子群。"""
    iso = IsoWrapper()
    kpoints = iso.list_k_points(225)
    labels = {kp.label for kp in kpoints}
    assert {"GM", "L", "X", "W"}.issubset(labels)

    subgroups = iso.list_subgroups(225, "GM", "GM5-")
    assert len(subgroups) >= 1
    # GM5- P1 对应子群 I-42m (#121)，指数 6
    p1 = next(sg for sg in subgroups if sg.opd_symbol == "P1")
    assert p1.space_group_number == 121
    assert p1.subgroup_index == 6


def test_iso_modes_and_domains():
    """BUSH 模式基矢与畴列表。"""
    iso = IsoWrapper()

    target = SubgroupInfo(
        index=0, space_group_number=107, space_group_symbol="I4mm",
        subgroup_index=6, size=1, is_maximal=True,
        opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
        basis_vectors=[[0, 0.5, -0.5], [0, 0.5, 0.5], [1, 0, 0]],
        origin=[0, 0, 0], k_point_label="GM", irrep_label="GM4-",
    )
    modes = iso.calc_distortion_modes(225, target, wyckoff_letters=["a"])
    assert len(modes) >= 1
    assert modes[0].bush_modes, "GM4- P1 在 4a 位点应存在位移模式"

    domains = iso.get_domains(225, target)
    assert len(domains) == 6
    assert domains[0].domain_number == 1


def _mode_core_token(label: str) -> str:
    """Identity of a complete-mode label ignoring parent-site display names.

    Official IsoVIZ compact labels drop the k prefix for Gamma *and* for
    I4/mmm M=(1,1,1) (written ``M1+[Al:d:dsp]…``).  Local pretty labels keep
    ``[1,1,1]``.  Both must hash to the same core token.
    """
    text = str(label)

    def _looks_like_k(token: str) -> bool:
        return "," in token or "/" in token

    m = re.search(
        r"\[([^]]+)\]([A-Za-z0-9+-]+)(?:\([^)]*\))?"
        r"\[[^:]+:([a-z]):dsp\]([A-Za-z0-9_]+(?:\([^)]*\))?)",
        text,
    )
    if m and _looks_like_k(m.group(1)):
        return f"{m.group(1)}|{m.group(2)}|{m.group(3)}|{m.group(4)}"
    m = re.search(
        r"([A-Za-z0-9+-]+)\[[^:]+:([a-z]):dsp\]([A-Za-z0-9_]+(?:\([^)]*\))?)",
        text,
    )
    if not m:
        return text.strip()
    irrep, letter, site = m.group(1), m.group(2), m.group(3)
    if irrep.startswith("GM"):
        k = "0,0,0"
    elif irrep.startswith("M"):
        k = "1,1,1"
    else:
        k = "0,0,0"
    return f"{k}|{irrep}|{letter}|{site}"


def _gold_ld1_c1_labels() -> list[str]:
    root = Path(__file__).resolve().parents[2] / "output_compare"
    if not root.is_dir():
        return []
    hits = list(root.rglob("data.isoviz"))
    for path in hits:
        if "ld1" in str(path).lower() and "c1" in str(path).lower():
            text = path.read_text(encoding="utf-8", errors="replace")
            if "!displacivemodelist" not in text:
                continue
            block = text.split("!displacivemodelist", 1)[1]
            labels = []
            for line in block.splitlines():
                if "dsp" not in line:
                    continue
                bits = line.split()
                for bit in bits:
                    if ":dsp]" in bit:
                        labels.append(bit)
                        break
            if labels:
                return labels
    return []


def test_eual4_ld1_c1_parametric_complete_modes():
    """smodes/(3+d) lock-in modes for EuAl4 LD g=1/6, LD1 C1.

    Does not hardcode 48.  Compares core tokens to official gold isoviz
    when ``output_compare`` is present; otherwise checks label families.
    """
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.is_file():
        pytest.skip("EuAl4 Parent.cif is not available")
    iso = IsoDistort(language="en")
    iso.set_distortion_types(["strain", "displacive"])
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    try:
        subs = iso.list_subgroups_at(
            "LD", "LD1", k_parameters=["1/6"], generate_if_missing=False,
        )
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"LD1 g=1/6 subgroup list unavailable: {exc}")
    target = next(
        (s for s in subs if s.irrep_label == "LD1" and s.opd_symbol == "C1"),
        None,
    )
    if target is None:
        pytest.skip("LD1 C1 is not in the cached subgroup list")
    result = iso.search_method_2(
        target.index,
        candidates=subs,
        number_of_independent_modulations=0,
    )
    assert result.modes
    assert all(
        mode.mode_identity is not None
        and mode.mode_identity.status == "verified"
        and mode.mode_identity.source == "iso_microscopic"
        and mode.microscopic_provenance is not None
        for mode in result.modes
    )
    labels = list((iso._mode_label_overrides or {}).values())
    assert labels, result.metadata.get("parametric_note") or "no complete modes"
    joined = " ".join(labels)
    assert "LD1" in joined
    assert any("0,0,1/6" in lab or "1/6" in lab for lab in labels)
    assert any("GM" in lab for lab in labels)
    assert any("M1+" in lab or "M3-" in lab for lab in labels)
    gold = _gold_ld1_c1_labels()
    if gold:
        local_core = {_mode_core_token(lab) for lab in labels}
        gold_core = {_mode_core_token(lab) for lab in gold}
        missing = gold_core - local_core
        extra = local_core - gold_core
        assert not missing, f"missing vs official: {sorted(missing)[:12]}"
        assert not extra, f"extra vs official: {sorted(extra)[:12]}"

    payload = iso.export_subgroups_zip(
        formats=["cif", "isoviz", "modes", "topas"],
        subgroups=[target],
        compute_missing_modes=True,
        export_method=2,
        number_of_independent_modulations=0,
    )
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        names = archive.namelist()
        assert names == [
            "LD1_C1/subgroup.cif",
            "LD1_C1/data.isoviz",
            "LD1_C1/Complete modes details.txt",
            "LD1_C1/topas.str",
        ]
        assert all(archive.getinfo(name).file_size > 0 for name in names)
        assert b"_iso_displacivemode_ID" in archive.read(names[0])
        assert b"!displacivemodelist" in archive.read(names[1])
        assert b"Displacive mode definitions" in archive.read(names[2])
        topas = archive.read(names[3]).lower()
        assert b"site " in topas and b"prm " in topas
