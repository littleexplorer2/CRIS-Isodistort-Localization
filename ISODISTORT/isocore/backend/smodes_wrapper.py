"""smodes 封装：参数 k 点上的位移模式（IR 活性 + 完整笛卡尔对称模）。

官网在 Method 2 枚举子群前，会按 Distortion Types + 物种作用域过滤
仅对结构有位移模式的 IR；对参数 k 点（LD/DT 等）使用 (3+d) 超空间 /
smodes 机制，而 iso 的 DISPLAY BUSH 无法完成该判断。

同一套 ``smodes`` 输出也给出每个 IR 的笛卡尔对称模（``smodes.txt``），
供 (3+d) / 公度锁定完整模式计算使用。
"""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

from .base_wrapper import BaseWrapper

_KVEC_RE = re.compile(
    r"^k vector\s+(\S+)\s*=\s*\(([^)]*)\)",
    re.IGNORECASE,
)
_IRREP_RE = re.compile(r"^Irrep\s+(\S+)")
_DEG_RE = re.compile(r"^Degeneracy:\s+(\d+)", re.IGNORECASE)
_NMODES_RE = re.compile(r"^Total number of modes:\s+(\d+)", re.IGNORECASE)
_ATOM_POS_RE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$"
)
_VEC_RE = re.compile(
    r"^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$"
)


@dataclass
class SmodesIrrepActivity:
    """单个 IR 在 smodes 输出中的活性摘要。"""

    label: str
    species: set[str] = field(default_factory=set)


@dataclass
class SmodesAtom:
    """smodes 超胞中的一个原子（笛卡尔 Å）。"""

    index: int
    species: str
    cart: np.ndarray


@dataclass
class SmodesModeBlock:
    """一个 ``---`` 分隔的对称模（笛卡尔位移，Å）。"""

    irrep: str
    degeneracy: int
    k_label: str
    k_coords: tuple[str, ...]
    superlattice_cart: np.ndarray
    atoms: list[SmodesAtom]
    atom_indices: list[int]
    displacements_cart: np.ndarray
    # Website labels keep the conventional representative arm, while the
    # actual block can be a rotated star arm.  Phase/quadrature operations
    # must use this physical arm rather than the display coordinates.
    active_k_coords: tuple[str, ...] | None = None

    def dense_displacements(self) -> np.ndarray:
        """Length-``len(atoms)`` cartesian displacements (zeros for omitted atoms)."""
        out = np.zeros((len(self.atoms), 3), dtype=float)
        index_to_row = {atom.index: i for i, atom in enumerate(self.atoms)}
        for idx, vec in zip(self.atom_indices, self.displacements_cart, strict=False):
            row = index_to_row.get(idx)
            if row is not None:
                out[row] = vec
        return out


class SmodesWrapper(BaseWrapper):
    """调用 smodes：IR 活性探测与完整对称模解析。"""

    def __init__(self) -> None:
        super().__init__()
        self.binary = str(Path(self.cfg.iso_bin).parent / "smodes")
        # A Method 2 ZIP computes many OPDs of the same commensurate cell.
        # Their complete-mode decomposition requests the same finite k grid;
        # keep parsed smodes output in this session instead of relaunching WSL
        # for every subgroup.  The cache is deliberately process-local and
        # bounded; no scientific result is persisted as an answer database.
        self._mode_block_cache: dict[tuple, list[SmodesModeBlock]] = {}

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
            k_parameters: smodes/iso（Kovalev）参数尺度；调用方负责从官网
                显示尺度转换（如 I4/mmm LD 官网 g=1/6 对应此处 ``["1/12"]``）
            species_filter: 若给定，仅当 IR 在该物种集合内至少一处有非零
                位移时才视为活性

        Returns:
            IR 标签集合；成功但无活性 IR 时返回空集。smodes 失败时返回
            ``None``（调用方应跳过该 k 点的过滤，避免把“无模式”当成“未知”）。
        """
        try:
            stdout = self.run_stdin(self.binary, self._build_input(
                structure, space_group_number, wyckoff_sites,
                [(k_label, k_parameters)],
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

    def compute_mode_blocks(
        self,
        structure: Structure,
        space_group_number: int,
        wyckoff_sites: list[dict],
        k_specs: list[tuple[str, list[str] | None]],
    ) -> list[SmodesModeBlock] | None:
        """Run smodes for one or more k points and return every symmetry-mode block."""
        if not k_specs:
            return []
        structure_key = (
            tuple(round(float(value), 8) for value in structure.lattice.matrix.ravel()),
            tuple(site.species_string for site in structure),
            tuple(
                round(float(value), 8)
                for value in np.asarray(structure.frac_coords, dtype=float).ravel()
            ),
        )
        cache_key = (
            int(space_group_number),
            structure_key,
            tuple(
                (
                    str(spec[0]),
                    tuple(str(value) for value in (spec[1] or ())),
                )
                for spec in k_specs
            ),
        )
        cached = self._mode_block_cache.get(cache_key)
        if cached is not None:
            return copy.deepcopy(cached)
        try:
            stdout = self.run_stdin(
                self.binary,
                self._build_input(
                    structure, space_group_number, wyckoff_sites, k_specs,
                ),
            )
        except Exception:  # noqa: BLE001 - caller may fall back
            return None
        parsed = parse_smodes_mode_blocks(stdout)
        if len(self._mode_block_cache) >= 32:
            self._mode_block_cache.pop(next(iter(self._mode_block_cache)))
        self._mode_block_cache[cache_key] = copy.deepcopy(parsed)
        return parsed

    @staticmethod
    def _build_input(structure: Structure,
                     space_group_number: int,
                     wyckoff_sites: list[dict],
                     k_specs: list[tuple[str, list[str] | None]] | str,
                     k_parameters: list[str] | None = None) -> str:
        if isinstance(k_specs, str):
            k_specs = [(k_specs, k_parameters)]
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
        lines.append(str(len(k_specs)))
        for k_label, params in k_specs:
            k_line = k_label.strip()
            if params:
                k_line += " " + " ".join(str(p).strip() for p in params)
            lines.append(k_line)
        return "\n".join(lines) + "\n"

    @staticmethod
    def _parse_irrep_activities(stdout: str) -> list[SmodesIrrepActivity]:
        """解析 smodes 输出，提取每个 IR 涉及的非零位移物种。"""
        by_label: dict[str, SmodesIrrepActivity] = {}
        for block in parse_smodes_mode_blocks(stdout):
            act = by_label.get(block.irrep)
            if act is None:
                act = SmodesIrrepActivity(label=block.irrep)
                by_label[block.irrep] = act
            dense = block.dense_displacements()
            for atom, vec in zip(block.atoms, dense, strict=False):
                if np.max(np.abs(vec)) <= 1e-8:
                    continue
                elem = re.match(r"([A-Z][a-z]?)", atom.species)
                if elem:
                    act.species.add(elem.group(1))
        return list(by_label.values())

    @staticmethod
    def smodes_kparams(k_parameters: list[str] | None) -> list[str] | None:
        """Normalize already-converted smodes/iso k parameters as fractions."""
        if not k_parameters:
            return None
        return [str(Fraction(str(p).strip())) for p in k_parameters]


def parse_smodes_mode_blocks(stdout: str) -> list[SmodesModeBlock]:
    """Parse cartesian symmetry-mode blocks from smodes stdout."""
    blocks: list[SmodesModeBlock] = []
    k_label = ""
    k_coords: tuple[str, ...] = ()
    irrep = ""
    degeneracy = 1
    superlattice: list[list[float]] = []
    atoms: list[SmodesAtom] = []
    reading_lattice = False
    reading_atoms = False
    reading_modes = False
    current_indices: list[int] = []
    current_disps: list[list[float]] = []

    def _flush_mode() -> None:
        nonlocal current_indices, current_disps
        if not current_indices or not irrep or not atoms:
            current_indices, current_disps = [], []
            return
        lat = (
            np.asarray(superlattice, dtype=float)
            if len(superlattice) == 3
            else np.eye(3)
        )
        blocks.append(
            SmodesModeBlock(
                irrep=irrep,
                degeneracy=degeneracy,
                k_label=k_label,
                k_coords=k_coords,
                superlattice_cart=lat,
                atoms=list(atoms),
                atom_indices=list(current_indices),
                displacements_cart=np.asarray(current_disps, dtype=float),
            )
        )
        current_indices, current_disps = [], []

    for raw in stdout.splitlines():
        line = raw.rstrip()
        stripped = line.strip()
        m_k = _KVEC_RE.match(stripped)
        if m_k:
            _flush_mode()
            k_label = m_k.group(1)
            k_coords = tuple(p.strip() for p in m_k.group(2).split(","))
            irrep = ""
            atoms = []
            superlattice = []
            reading_lattice = reading_atoms = reading_modes = False
            continue
        m_ir = _IRREP_RE.match(stripped)
        if m_ir:
            _flush_mode()
            irrep = m_ir.group(1)
            degeneracy = 1
            superlattice = []
            atoms = []
            reading_lattice = reading_atoms = reading_modes = False
            continue
        m_deg = _DEG_RE.match(stripped)
        if m_deg:
            degeneracy = int(m_deg.group(1))
            continue
        m_n = _NMODES_RE.match(stripped)
        if m_n:
            continue
        if stripped.lower().startswith("vectors defining superlattice"):
            reading_lattice = True
            superlattice = []
            reading_atoms = reading_modes = False
            continue
        if reading_lattice:
            m_vec = _VEC_RE.match(stripped)
            if m_vec:
                superlattice.append([float(m_vec.group(i)) for i in range(1, 4)])
                if len(superlattice) >= 3:
                    reading_lattice = False
                continue
            if stripped:
                reading_lattice = False
        if stripped.lower().startswith("atom") and "position" in stripped.lower():
            reading_atoms = True
            atoms = []
            reading_modes = False
            continue
        if "Symmetry modes:" in stripped:
            reading_atoms = False
            reading_modes = True
            continue
        if reading_atoms:
            m_atom = _ATOM_POS_RE.match(line)
            if m_atom:
                atoms.append(
                    SmodesAtom(
                        index=int(m_atom.group(1)),
                        species=m_atom.group(2),
                        cart=np.array(
                            [float(m_atom.group(i)) for i in range(3, 6)],
                            dtype=float,
                        ),
                    )
                )
                continue
            if stripped.startswith("*") or stripped.startswith("Irrep"):
                reading_atoms = False
            continue
        if not reading_modes:
            continue
        if stripped.startswith("---"):
            _flush_mode()
            continue
        if "displacement" in stripped and "atom" in stripped:
            continue
        if stripped.lower().startswith("the following include"):
            continue
        m_disp = _ATOM_POS_RE.match(line)
        if m_disp:
            current_indices.append(int(m_disp.group(1)))
            current_disps.append([float(m_disp.group(i)) for i in range(3, 6)])
            continue
        if stripped.startswith("Irrep") or stripped.startswith("*"):
            _flush_mode()
            reading_modes = False
    _flush_mode()
    return blocks
