"""
Distortion Page 批量导出：CIF / Save interactive distortion / Complete modes details / TOPAS.STR

对齐官网第 6 页（webpage_info/6. ISODISTORT_ distortion.html）与手册
https://landau3.byu.edu/isodistorthelp.php#modeparams 中的导出选项。

命名约定：
    文件夹  Method 1: ``<IR>_<OPD>_SG<number>``；
            Method 2: ``<IR>_<OPD>``；
            Method 3: ``<case-id>/C<sequence>_SG<number>``
    文件    ``subgroup.cif`` / ``data.isoviz`` /
            ``Complete modes details.txt``（官网为 HTML；本地用 .txt）/
            ``topas.str``
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import tempfile
import zipfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

from ..backend import SubgroupInfo
from ..distortion.strain_modes import (
    ENGINEERING_VOIGT_ORDER,
    HomogeneousStrainMode,
    HomogeneousStrainModeResult,
    engineering_voigt_to_tensor,
)
from ..utils.config_loader import get_config
from ..utils.schoenflies import hm_symbol
from .displacive_export import DisplaciveExportData, ParentOrbitType

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

_WINDOWS_BAD = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
_WINDOWS_RESERVED = {
    "CON", "PRN", "AUX", "NUL",
    *(f"COM{number}" for number in range(1, 10)),
    *(f"LPT{number}" for number in range(1, 10)),
}
_MAX_FOLDER_COMPONENT_LENGTH = 120


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
    """Return one portable Windows-safe folder component.

    ``+``/``-`` and crystallographic tokens such as ``4D1`` are preserved.
    Path separators, control characters and Windows device names are never
    emitted.  Long components keep a short content digest so truncation cannot
    silently merge distinct names.
    """
    raw = str(text or "").strip()
    cleaned = _WINDOWS_BAD.sub("_", raw)
    cleaned = re.sub(r"\s+", "_", cleaned)
    cleaned = re.sub(r"_+", "_", cleaned).strip(" ._")
    if not cleaned:
        cleaned = fallback
    stem = cleaned.split(".", 1)[0].upper()
    if stem in _WINDOWS_RESERVED:
        cleaned = f"_{cleaned}"
    if len(cleaned) > _MAX_FOLDER_COMPONENT_LENGTH:
        digest = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:10]
        keep = _MAX_FOLDER_COMPONENT_LENGTH - len(digest) - 1
        cleaned = f"{cleaned[:keep].rstrip(' ._')}_{digest}"
    return cleaned


def subgroup_label(subgroup: SubgroupInfo) -> str:
    """子群显示名：``IR OPD``，例如 ``LD1 C1``。"""
    ir = safe_name(subgroup.irrep_label or "", "IR")
    opd = safe_name(subgroup.opd_symbol or "", "OPD")
    return f"{ir} {opd}"


def _identity_payload(subgroup: SubgroupInfo) -> dict:
    """Serializable scientific identity used only for deterministic names."""

    def values(items) -> list[str]:
        if items is None:
            return []
        return [str(value) for value in items]

    def matrix(rows) -> list[list[str]]:
        if rows is None:
            return []
        return [values(row) for row in rows]

    return {
        "parent_sg": int(subgroup.parent_sg or 0),
        "space_group_number": int(subgroup.space_group_number or 0),
        "subgroup_index": int(subgroup.subgroup_index or 0),
        "size": int(subgroup.size or 0),
        "k_point_label": str(subgroup.k_point_label or ""),
        "k_coordinates": values(subgroup.k_coordinates),
        "k_parameters": values(subgroup.k_parameters),
        "irrep_label": str(subgroup.irrep_label or ""),
        "opd_symbol": str(subgroup.opd_symbol or ""),
        "opd_dir_raw": str(subgroup.opd_dir_raw or ""),
        "opd_vector": values(subgroup.opd_vector),
        "basis_raw": str(subgroup.basis_raw or ""),
        "basis_vectors": matrix(subgroup.basis_vectors),
        "origin_raw": str(subgroup.origin_raw or ""),
        "origin": values(subgroup.origin),
        "k_active_raw": str(subgroup.k_active_raw or ""),
        "method3_embedding_id": str(
            getattr(subgroup, "_method3_embedding_id", "") or ""
        ),
    }


def subgroup_identity_digest(subgroup: SubgroupInfo, length: int = 10) -> str:
    """Short deterministic digest of the full subgroup/path identity."""
    payload = json.dumps(
        _identity_payload(subgroup),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def method3_case_folder(
    subgroups: Iterable[SubgroupInfo],
    stable_case_id: str | None = None,
) -> str:
    """Return a portable Method-3 case folder.

    A caller-supplied manifest ID is preserved after sanitizing.  Interactive
    exports have no such ID, so they receive a content-addressed generic ID
    derived from the complete candidate set; no crystal/sample names are
    embedded in production logic.
    """
    if stable_case_id and str(stable_case_id).strip():
        return safe_name(stable_case_id, "M3")
    identities = sorted(
        json.dumps(
            _identity_payload(subgroup),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        for subgroup in subgroups
    )
    digest = hashlib.sha256("\n".join(identities).encode("utf-8")).hexdigest()[:10]
    return f"M3-{digest}"


def folder_label_for_subgroup(
    subgroup: SubgroupInfo,
    *,
    export_method: int = 2,
    sequence: int | None = None,
) -> str:
    """Return the short candidate-folder base name for one export Method."""
    method = parse_export_method(export_method)
    ir = safe_name(subgroup.irrep_label or "", "IR")
    opd = safe_name(subgroup.opd_symbol or "", "OPD")
    sg_number = int(subgroup.space_group_number or 0)
    if method == 1:
        return f"{ir}_{opd}_SG{sg_number}"
    if method == 2:
        return f"{ir}_{opd}"
    ordinal = sequence if sequence is not None else int(subgroup.index) + 1
    ordinal = max(1, int(ordinal))
    return f"C{ordinal:02d}_SG{sg_number}"


def _deduplicate_folder_name(
    base: str,
    subgroup: SubgroupInfo,
    used: set[str],
) -> str:
    """Reserve a case-insensitively unique folder name without overwriting."""
    occupied = {name.casefold() for name in used}
    name = safe_name(base)
    if name.casefold() not in occupied:
        used.add(name)
        return name

    digest = subgroup_identity_digest(subgroup)
    candidate = safe_name(f"{name}_{digest}")
    if candidate.casefold() not in occupied:
        used.add(candidate)
        return candidate

    suffix = 2
    while True:
        candidate = safe_name(f"{name}_{digest}_{suffix}")
        if candidate.casefold() not in occupied:
            used.add(candidate)
            return candidate
        suffix += 1


def unique_folder_name(
    subgroup: SubgroupInfo,
    used: set[str],
    *,
    export_method: int = 2,
    sequence: int | None = None,
) -> str:
    """Reserve one short collision-safe name in a batch or destination."""
    base = folder_label_for_subgroup(
        subgroup,
        export_method=export_method,
        sequence=sequence,
    )
    return _deduplicate_folder_name(base, subgroup, used)


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
class StrainExportData:
    """One authoritative homogeneous-strain contract for strain-aware outputs.

    The basis uses official parent-lattice-basis engineering order
    ``(E11,E22,E33,2E23,2E13,2E12)``. CIF's matrix and ``_iso_strain_value``
    use the website's raw-coordinate sum ``sum(a_i*q_raw_i)``. IsoVIZ and the
    actual multiplier ``M`` use normalized vectors, so applied engineering
    strain is ``sum(a_i*q_unit_i)`` with
    ``q_unit_i=normfactor_i*q_raw_i``. These two six-vectors are deliberately
    separate because the authoritative EuAl4 F02 files show they differ.

    Production export requires the ISO rank-[12] canonical macroscopic basis
    and a successful whole-space check against the independent metric fixed
    space. The metric-only deterministic basis remains diagnostic and cannot
    silently masquerade as website-compatible output.
    """

    result: HomogeneousStrainModeResult
    amplitudes: tuple[float, ...] = ()
    max_amplitude: float = 0.1

    def __post_init__(self) -> None:
        count = len(self.result.modes)
        values = tuple(float(value) for value in self.amplitudes)
        if not values:
            values = (0.0,) * count
        if len(values) != count:
            raise ValueError("strain amplitude count must equal the strain-mode count")
        if not all(np.isfinite(value) for value in values):
            raise ValueError("strain amplitudes must be finite")
        if not np.isfinite(self.max_amplitude) or self.max_amplitude <= 0.0:
            raise ValueError("strain max amplitude must be positive and finite")
        if tuple(self.result.voigt_order) != tuple(ENGINEERING_VOIGT_ORDER):
            raise ValueError("strain result uses an unsupported Voigt ordering")
        if not self.result.export_ready:
            raise ValueError(
                "strain basis is unresolved: canonical ISO rank-[12] modes "
                "were not validated against the metric fixed space"
            )
        for expected, mode in enumerate(self.result.modes, start=1):
            if int(mode.index) != expected:
                raise ValueError("strain mode indices must be contiguous and one-based")
            if not (
                mode.label_status == "canonical_iso_macro"
                and mode.canonical_label
                and mode.irrep_label
                and mode.irrep_direction
            ):
                raise ValueError(
                    "export-ready strain modes require verified ISO labels and directions"
                )
            q_raw = np.asarray(mode.q_raw, dtype=float)
            q_unit = np.asarray(mode.q_unit, dtype=float)
            if q_raw.shape != (6,) or q_unit.shape != (6,):
                raise ValueError("each strain mode must contain six Voigt components")
            if not np.all(np.isfinite(q_raw)) or not np.all(np.isfinite(q_unit)):
                raise ValueError("strain mode components must be finite")
            if not np.allclose(
                q_unit,
                float(mode.normfactor) * q_raw,
                rtol=2.0e-12,
                atol=2.0e-12,
            ):
                raise ValueError("strain q_unit must equal normfactor * q_raw")
        self.amplitudes = values

    @staticmethod
    def label(mode: HomogeneousStrainMode) -> str:
        """Return the validated ISO label, never a guessed Gamma label."""

        if mode.canonical_label and mode.label_status == "canonical_iso_macro":
            return str(mode.canonical_label)
        raise ValueError("strain mode has no verified canonical ISO label")

    @classmethod
    def full_label(
        cls, mode: HomogeneousStrainMode, parent_symbol: str
    ) -> str:
        """Return the CIF/Complete label including the parent prefix."""

        compact = cls.label(mode)
        irrep = str(mode.irrep_label or "")
        direction = str(mode.irrep_direction or "")
        if irrep and direction and compact.startswith(irrep):
            compact = f"{irrep}{direction}{compact[len(irrep):]}"
        symbol = str(parent_symbol).strip()
        return compact if not symbol else f"{symbol}[0,0,0]{compact}"

    def items(self) -> tuple[tuple[HomogeneousStrainMode, str, float], ...]:
        """Return modes, shared labels, and amplitudes in canonical order."""

        return tuple(
            (mode, self.label(mode), self.amplitudes[position])
            for position, mode in enumerate(self.result.modes)
        )

    def cif_raw_component_values(self) -> np.ndarray:
        """Return official CIF raw coordinates ``sum(a_i*q_raw_i)``."""

        if not self.result.modes:
            return np.zeros(6, dtype=float)
        matrix = np.column_stack(
            [np.asarray(mode.q_raw, dtype=float) for mode in self.result.modes]
        )
        return matrix @ np.asarray(self.amplitudes, dtype=float)

    def applied_engineering_q(self) -> np.ndarray:
        """Return physical applied q used in ``M``: ``sum(a_i*q_unit_i)``."""

        if not self.result.modes:
            return np.zeros(6, dtype=float)
        matrix = np.column_stack(
            [np.asarray(mode.q_unit, dtype=float) for mode in self.result.modes]
        )
        return matrix @ np.asarray(self.amplitudes, dtype=float)


@dataclass
class SubgroupExportSpec:
    """单个子群一次导出所需的全部数据。"""

    subgroup: SubgroupInfo
    structure: Structure | None = None
    parent_structure: Structure | None = None
    parent_sg: int = 0
    parent_symbol: str = ""
    mode_displacements_sc: dict[str, np.ndarray] | None = None
    mode_labels: dict[str, str] | None = None
    amplitudes: dict[str, float] | None = None
    note: str = ""
    folder_name: str = ""
    # Optional already-final structure.  With strain_data its metric must match
    # B @ M @ P; CIF/Complete may use its displaced fractional coordinates.
    cif_structure: Structure | None = None
    parent_wyckoff_sites: list | None = None
    # Pre-formatted parent Wyckoff comment lines (CIF order/labels when available).
    parent_wyckoff_lines: list[str] | None = None
    distortion_types: list[str] = field(default_factory=list)
    strain_data: StrainExportData | None = None
    # Validated single source for every mode-aware writer.  Production API
    # specs set ``require_verified_displacive_data`` whenever they carry
    # displacive arrays; direct legacy writer fixtures can remain opt-in while
    # the upstream exact frame/atom mapper is migrated.
    displacive_data: DisplaciveExportData | None = None
    require_verified_displacive_data: bool = False
    mode_keys: tuple[str, ...] = field(default=(), init=False)
    mode_global_irreps: dict[str, str] | None = field(default=None, init=False)
    atom_ids: tuple[str, ...] = field(default=(), init=False)
    parent_orbit_types: Mapping[str, ParentOrbitType] | None = field(
        default=None,
        init=False,
    )
    _displacive_contract_required: bool = field(default=False, init=False, repr=False)

    def __post_init__(self) -> None:
        self._displacive_contract_required = bool(
            self.displacive_data is not None
            or self.require_verified_displacive_data
        )
        if self.strain_data is not None and "strain" not in {
            str(value).strip().lower() for value in self.distortion_types
        }:
            raise ValueError(
                "strain_data is valid only when distortion_types includes 'strain'"
            )
        if self.displacive_data is None:
            if self.structure is None:
                raise ValueError("export structure is required without displacive_data")
        else:
            parallel = {
                "structure": self.structure,
                "parent_structure": self.parent_structure,
                "parent_sg": self.parent_sg or None,
                "parent_symbol": self.parent_symbol or None,
                "mode_displacements_sc": self.mode_displacements_sc,
                "mode_labels": self.mode_labels,
                "amplitudes": self.amplitudes,
                "cif_structure": self.cif_structure,
            }
            supplied = [name for name, value in parallel.items() if value is not None]
            if supplied:
                raise ValueError(
                    "displacive_data is the only source for " + ", ".join(supplied)
                )
            self._bind_displacive_data()

    def _bind_displacive_data(self) -> None:
        """Revalidate and restore every mutable writer-facing projection."""

        data = self.displacive_data
        if data is None:
            return
        data.validate_subgroup(self.subgroup)
        rows = data.modes()
        self.parent_structure = data.snapshot_parent()
        self.parent_sg = data.parent_space_group_number
        self.parent_symbol = hm_symbol(self.parent_sg) or str(self.parent_sg)
        self.structure = data.snapshot_reference()
        self.cif_structure = data.snapshot_final()
        self.mode_displacements_sc = {
            row.key: np.array(row.raw_fractional, dtype=float, copy=True)
            for row in rows
        }
        self.mode_keys = tuple(row.key for row in rows)
        self.mode_labels = {row.key: row.label for row in rows}
        self.mode_global_irreps = {
            row.key: row.global_irrep_label for row in rows
        }
        self.atom_ids = data.atom_ids
        self.parent_orbit_types = data.parent_orbit_types
        self.amplitudes = {
            row.key: row.amplitude_as_angstrom for row in rows
        }

    def assert_displacive_export_ready(self) -> None:
        """Reject production mode export when exact source evidence was lost."""

        if self._displacive_contract_required and self.displacive_data is None:
            raise ValueError(
                "displacive export contract was removed or missing after preparation; "
                "exact emitted-setting frame, child atom order and unmixed ISO "
                "source-column provenance are required"
            )
        if self.displacive_data is not None:
            self._bind_displacive_data()
        if (
            self.require_verified_displacive_data
            and self.mode_displacements_sc
            and self.displacive_data is None
        ):
            raise ValueError(
                "displacive export is unresolved: exact emitted-setting frame, "
                "child atom order and unmixed ISO source-column provenance are required"
            )


def resolve_strained_export_structure(spec: SubgroupExportSpec) -> Structure:
    """Return the structure whose lattice closes the shared strain contract.

    ``spec.structure`` is the undistorted child reference and must have metric
    ``(B @ P) @ (B @ P).T``.  When strain amplitudes are present, the final
    lattice is constructed once as ``B @ M @ P``.  A verified
    ``displacive_data`` final is defined in the reference frame, so its
    displaced fractional coordinates are carried onto that reconstructed
    lattice.  A legacy explicit ``spec.cif_structure`` is instead treated as
    an already-final structure and is accepted only when its metric matches
    the same reconstructed lattice.
    """

    reference = spec.structure
    final = spec.cif_structure
    strain = spec.strain_data
    if strain is None:
        return final if final is not None else reference
    if spec.parent_structure is None:
        raise ValueError(
            "strain export requires the actual parent structure"
        )

    parent_lattice = np.asarray(
        spec.parent_structure.lattice.matrix, dtype=float
    )
    basis_values = spec.subgroup.basis_vectors
    parent_to_child_basis = (
        np.asarray(basis_values, dtype=float)
        if basis_values and len(basis_values) == 3
        else np.eye(3)
    )
    if parent_to_child_basis.shape != (3, 3) or not np.all(
        np.isfinite(parent_to_child_basis)
    ):
        raise ValueError("strain export requires a finite 3x3 child basis")
    if abs(float(np.linalg.det(parent_to_child_basis))) <= np.finfo(float).eps:
        raise ValueError("strain export requires a nonsingular child basis")

    tolerance = float(get_config().lattice_tolerance)

    def metric(matrix: np.ndarray) -> np.ndarray:
        return matrix @ matrix.T

    def relative_metric_residual(
        actual: np.ndarray, expected: np.ndarray
    ) -> float:
        expected_metric = metric(expected)
        return float(
            np.linalg.norm(metric(actual) - expected_metric)
            / max(float(np.linalg.norm(expected_metric)), 1.0)
        )

    expected_reference_lattice = parent_to_child_basis @ parent_lattice
    reference_residual = relative_metric_residual(
        np.asarray(reference.lattice.matrix, dtype=float),
        expected_reference_lattice,
    )
    if reference_residual > tolerance:
        raise ValueError(
            "strain export reference child lattice is incompatible with B @ P "
            f"(relative metric residual {reference_residual:.6g} > "
            f"{tolerance:.6g})"
        )

    applied_q = strain.applied_engineering_q()
    multiplier = np.eye(3) + engineering_voigt_to_tensor(applied_q)
    expected_final_lattice = (
        parent_to_child_basis @ multiplier @ parent_lattice
    )
    if not np.all(np.isfinite(expected_final_lattice)) or abs(
        float(np.linalg.det(expected_final_lattice))
    ) <= np.finfo(float).eps:
        raise ValueError("strain export reconstructed a singular final lattice")

    if spec.displacive_data is not None:
        if final is None:
            raise ValueError(
                "verified displacive export data has no final child structure"
            )
        displacive_final_lattice = np.asarray(final.lattice.matrix, dtype=float)
        displacive_final_residual = float(
            np.linalg.norm(displacive_final_lattice - expected_reference_lattice)
            / max(float(np.linalg.norm(expected_reference_lattice)), 1.0)
        )
        if not np.allclose(
            displacive_final_lattice,
            expected_reference_lattice,
            rtol=2.0e-12,
            atol=tolerance,
        ):
            raise ValueError(
                "verified displacive final lattice is incompatible with its "
                "B @ P reference frame "
                f"(relative lattice residual {displacive_final_residual:.6g} > "
                f"{tolerance:.6g})"
            )
        return Structure(
            Lattice(expected_final_lattice),
            [site.species for site in final],
            final.frac_coords,
            charge=final.charge,
            coords_are_cartesian=False,
            site_properties=final.site_properties,
            labels=[site.label for site in final],
            properties=final.properties,
        )

    if final is not None:
        final_residual = relative_metric_residual(
            np.asarray(final.lattice.matrix, dtype=float),
            expected_final_lattice,
        )
        if final_residual > tolerance:
            raise ValueError(
                "explicit final CIF structure is incompatible with B @ M @ P "
                f"(relative metric residual {final_residual:.6g} > "
                f"{tolerance:.6g})"
            )
        return final

    if np.allclose(applied_q, 0.0, rtol=0.0, atol=1.0e-15):
        return reference
    return Structure(
        Lattice(expected_final_lattice),
        [site.species for site in reference],
        reference.frac_coords,
        charge=reference.charge,
        coords_are_cartesian=False,
        site_properties=reference.site_properties,
        labels=[site.label for site in reference],
        properties=reference.properties,
    )


def render_cif(structure: Structure, spec: SubgroupExportSpec | None = None) -> str:
    """ISODISTORT-style CIF (subgroup setting when ``spec`` is given)."""
    from .isodistort_cif import render_isodistort_cif  # noqa: PLC0415

    if spec is not None:
        spec.assert_displacive_export_ready()
    return render_isodistort_cif(structure, spec)


def render_isoviz(spec: SubgroupExportSpec) -> str:
    """Save interactive distortion：官网 ``!tag`` ISOVIZ 数据文件布局。"""
    from .isodistort_isoviz import render_isodistort_isoviz  # noqa: PLC0415

    spec.assert_displacive_export_ready()
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

    Acceptance (see ``ISODISTORT/agent.md``): put **all locally computed** supercell /
    mode information into this text file inside each subgroup folder of the
    Distortion ZIP. Byte-level match to the official HTML modes page is **not**
    required.

    Content outline (handbook #modesdetails inspiration): supercell xyz table,
    mode definitions (every atom), As / Ap / dmax when displacements exist.
    """
    spec.assert_displacive_export_ready()
    if spec.displacive_data is not None:
        from .isodistort_cif import _validated_subgroup_sites  # noqa: PLC0415

        _validated_subgroup_sites(spec, spec.structure)
    sg = spec.subgroup
    sc = spec.structure
    strain = spec.strain_data
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

    def append_superstructure(title: str, structure: Structure) -> None:
        lattice = structure.lattice
        lines.append(title)
        lines.append(
            f"  a={lattice.a:.5f}  b={lattice.b:.5f}  c={lattice.c:.5f}  "
            f"alpha={lattice.alpha:.5f}  beta={lattice.beta:.5f}  "
            f"gamma={lattice.gamma:.5f}"
        )
        lines.append(f"{'atom':>6s} {'el':<4s} {'x':>10s} {'y':>10s} {'z':>10s}")
        for index, site in enumerate(structure, start=1):
            x, y, z = (float(value) for value in site.frac_coords)
            lines.append(
                f"{index:6d} {site.species_string:<4s} "
                f"{x:10.6f} {y:10.6f} {z:10.6f}"
            )

    if strain is None and spec.displacive_data is None:
        append_superstructure("Lattice parameters of the supercell:", sc)
    else:
        # The official Complete-modes contract shows both states: mode
        # definitions use the undistorted reference child, while the distorted
        # structure uses the same authoritative B @ M @ P lattice as CIF/TOPAS.
        append_superstructure("Undistorted superstructure:", sc)
        lines.append("")
        append_superstructure("Distorted superstructure:", resolve_strained_export_structure(spec))

    if strain is not None:
        lines.append("")
        lines.append("Parent-cell homogeneous strain mode definitions")
        lines.append(
            "Parent-lattice-basis engineering order q = "
            "(E11, E22, E33, 2E23, 2E13, 2E12)."
        )
        lines.append(
            "M=I+E; column basis C -> C@M, row parent lattice P -> M@P, "
            "and child lattice L=B@M@P."
        )
        lines.append(
            "q_unit = normfactor * q_raw and has unit tensor Frobenius norm. "
            "CIF raw components use sum(amplitude*q_raw); applied M uses "
            "sum(amplitude*q_unit)."
        )
        lines.append(
            f"Basis source: {strain.result.basis_source}; fixed-space validation "
            f"max error={strain.result.fixed_space_validation_max_error:.6g}."
        )
        for mode, _label, amplitude in strain.items():
            label = strain.full_label(mode, spec.parent_symbol)
            raw = " ".join(f"{float(value): .8f}" for value in mode.q_raw)
            unit = " ".join(f"{float(value): .8f}" for value in mode.q_unit)
            lines.append("")
            lines.append(f"Mode {label}  label_status={mode.label_status}")
            lines.append(f"  normfactor = {float(mode.normfactor):.10g}")
            lines.append(f"  q_raw  = ({raw})")
            lines.append(f"  q_unit = ({unit})")
            lines.append(f"  amplitude = {amplitude:.10g}")
        raw_components = strain.cif_raw_component_values()
        applied_components = strain.applied_engineering_q()
        lines.append("")
        lines.append("Official CIF raw-coordinate components:")
        lines.append(
            "  "
            + "  ".join(
                f"{name}={float(value):.10g}"
                for name, value in zip(
                    ENGINEERING_VOIGT_ORDER, raw_components, strict=True
                )
            )
        )
        lines.append("Applied engineering q used in M:")
        lines.append(
            "  "
            + "  ".join(
                f"{name}={float(value):.10g}"
                for name, value in zip(
                    ENGINEERING_VOIGT_ORDER, applied_components, strict=True
                )
            )
        )

    lines.append("")
    lines.append("Displacive mode definitions (every atom in the unit cell)")
    lines.append(
        "Mode vectors are given in unitless superlattice coordinates; "
        "the largest component of each mode is scaled to 1.0. "
        "normfactor makes the sum of squares of Cartesian changes over all "
        "atoms equal to 1.0 (Angstrom for displacive modes)."
    )
    if not spec.mode_displacements_sc:
        lines.append("(no displacive modes available for this subgroup)")
        return "\n".join(lines) + "\n"

    primitive_size_ratio = int(sg.size or 0)
    if primitive_size_ratio <= 0:
        raise ValueError(
            "subgroup primitive-cell size ratio s must be positive for mode amplitudes"
        )
    # ISODISTORT defines Ap = As*sqrt(Vp/Vs).  The OPD field ``s`` is exactly
    # Vs/Vp for the primitive subgroup and parent cells, so this remains
    # correct for centred conventional cells and non-standard settings.
    scale_ap = 1.0 / float(np.sqrt(primitive_size_ratio))
    centering_mult = _centering_multiplicity(sg.space_group_symbol)

    # Displacive modes are defined in the undistorted child reference basis,
    # matching the official Complete-modes normfactor convention.
    bmat = np.asarray(sc.lattice.matrix, dtype=float)
    amp_rows: list[str] = []
    authoritative_modes = (
        {row.key: row for row in spec.displacive_data.modes()}
        if spec.displacive_data is not None
        else {}
    )
    for label, disp in spec.mode_displacements_sc.items():
        pretty = (spec.mode_labels or {}).get(label, label)
        arr = np.asarray(disp, dtype=float)
        if arr.size == 0:
            continue
        authoritative = authoritative_modes.get(label)
        if authoritative is None:
            max_comp = float(np.max(np.abs(arr)))
            unit = arr / max_comp if max_comp > 1e-16 else arr
            cart = unit @ bmat
            ssq_conventional = float(np.sum(cart * cart))
            ssq_primitive = ssq_conventional / centering_mult
            norm = (
                (1.0 / np.sqrt(ssq_primitive))
                if ssq_primitive > 1e-30
                else 0.0
            )
            as_amp = float((spec.amplitudes or {}).get(label, 0.0))
            dmax = (
                abs(as_amp) * norm * float(np.max(np.linalg.norm(cart, axis=1)))
                if norm > 0.0
                else 0.0
            )
        else:
            unit = authoritative.raw_fractional
            norm = authoritative.normfactor_per_angstrom
            as_amp = authoritative.amplitude_as_angstrom
            dmax = abs(as_amp) / authoritative.max_amplitude_angstrom
        ap_amp = as_amp * scale_ap
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
        cfg = get_config()
        sga = SpacegroupAnalyzer(
            structure,
            symprec=cfg.symmetry_cartesian_tolerance_angstrom,
            angle_tolerance=cfg.symmetry_angle_tolerance_degrees,
        )
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
        for raw_value in match.group(1).split(","):
            raw = raw_value.strip()
            try:
                val = float(raw)
            except ValueError:
                parts.append(raw)
                continue
            if abs(val - round(val)) < 1e-8:
                parts.append(str(round(val)))
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
    ``maxamp_hint`` is the amplitude at which the largest Cartesian atomic
    displacement reaches 1 Angstrom, i.e. ``1 / dmax(As=1)``.  This follows
    the ISODISTORT definitions of the supercell-normalized amplitude ``As``
    and ``dmax`` rather than a space-group-centring heuristic.

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
    cartesian_at_as_one = scaled @ bmat
    dmax_at_as_one = float(
        np.max(np.linalg.norm(cartesian_at_as_one, axis=1))
    )
    maxamp = (1.0 / dmax_at_as_one) if dmax_at_as_one > 1e-30 else 1.0
    return scaled, maxamp


def render_topas(spec: SubgroupExportSpec) -> str:
    """TOPAS.STR：官网 distortion-mode 精修输入（见手册 #topas）。"""
    from ..utils.parent_header import format_fixed_coord  # noqa: PLC0415
    from .isodistort_cif import (  # noqa: PLC0415
        _parent_to_child_transform,
        _validated_subgroup_sites,
    )

    spec.assert_displacive_export_ready()
    sg = spec.subgroup
    _setting, sc, subgroup_sites, _origin_shift = _validated_subgroup_sites(
        spec, spec.structure
    )
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
    ]
    strain = spec.strain_data
    if strain is not None:
        lat = resolve_strained_export_structure(spec).lattice
    lines.extend(
        [
            f"\t\ta  {lat.a:10.5f}",
            f"\t\tb  {lat.b:10.5f}",
            f"\t\tc  {lat.c:10.5f}",
            f"\t\tal {lat.alpha:10.5f}",
            f"\t\tbe {lat.beta:10.5f}",
            f"\t\tga {lat.gamma:10.5f}",
        ]
    )
    lines.extend(
        [
            "\t\tscale @ 0.00001",
            "",
            "'{{{mode definitions",
        ]
    )
    mode_items = list((spec.mode_displacements_sc or {}).items())
    authoritative_modes = (
        {row.key: row for row in spec.displacive_data.modes()}
        if spec.displacive_data is not None
        else {}
    )
    unique = [int(site["index"]) for site in subgroup_sites]
    tags = [str(site["label"]) for site in subgroup_sites]
    scaled_modes: list[np.ndarray] = []
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
            authoritative = authoritative_modes.get(label)
            if authoritative is None:
                amp = float((spec.amplitudes or {}).get(label, 0.0))
                scaled, amp_bound = cart_normalized_mode_matrix(
                    np.asarray(disp, dtype=float), lat.matrix, centering_mult=n_c
                )
            else:
                amp = authoritative.amplitude_as_angstrom
                scaled = authoritative.normalized_fractional_per_angstrom
                amp_bound = authoritative.max_amplitude_angstrom
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


def render_subgroup_files(
    spec: SubgroupExportSpec,
    formats: Sequence[str],
) -> tuple[tuple[str, bytes], ...]:
    """Render and validate every requested file without touching the filesystem."""

    folder = spec.folder_name or subgroup_label(spec.subgroup)
    names = [format_filename(folder, fmt) for fmt in formats]
    casefolded_names = [name.casefold() for name in names]
    if len(casefolded_names) != len(set(casefolded_names)):
        raise ValueError("导出格式产生了重复文件名")
    return tuple(
        (name, render_format(fmt, spec))
        for name, fmt in zip(names, formats, strict=True)
    )


def write_rendered_subgroup_files(
    dest_dir: Path,
    rendered: Sequence[tuple[str, bytes]],
) -> list[Path]:
    """Atomically publish one already-rendered candidate directory."""

    entries = tuple((str(name), bytes(payload)) for name, payload in rendered)
    names = [name for name, _payload in entries]
    if any(Path(name).name != name for name in names):
        raise ValueError("导出文件名必须是单个安全路径组件")
    casefolded_names = [name.casefold() for name in names]
    if len(casefolded_names) != len(set(casefolded_names)):
        raise ValueError("导出格式产生了重复文件名")
    if dest_dir.exists() or dest_dir.is_symlink():
        raise FileExistsError(f"目标目录已存在: {dest_dir}")

    dest_dir.parent.mkdir(parents=True, exist_ok=True)
    if dest_dir.exists() or dest_dir.is_symlink():
        raise FileExistsError(f"目标目录已存在: {dest_dir}")

    with tempfile.TemporaryDirectory(
        prefix=".isodistort-export-",
        dir=dest_dir.parent,
    ) as temp_name:
        temp_dir = Path(temp_name)
        for name, payload in entries:
            (temp_dir / name).write_bytes(payload)

        # A sibling rename is an atomic directory publication on Windows.  It
        # also refuses an existing destination there; the explicit recheck
        # preserves that contract on other supported platforms.
        if dest_dir.exists() or dest_dir.is_symlink():
            raise FileExistsError(f"目标目录已存在: {dest_dir}")
        temp_dir.rename(dest_dir)

    return [dest_dir / name for name in names]


def write_subgroup_files(
    dest_dir: Path,
    spec: SubgroupExportSpec,
    formats: Sequence[str],
) -> list[Path]:
    """把所有格式原子发布到全新的候选目录，拒绝覆盖已有目录。"""
    # Rendering can fail for a later format.  Keep it entirely in memory so a
    # failure cannot expose a directory containing only the earlier formats.
    rendered = render_subgroup_files(spec, formats)
    return write_rendered_subgroup_files(dest_dir, rendered)


def build_export_zip(
    specs: Iterable[SubgroupExportSpec],
    formats: Sequence[str],
    wrapping: str | None = None,
    *,
    rendered_batch: Sequence[
        tuple[SubgroupExportSpec, Sequence[tuple[str, bytes]]]
    ] | None = None,
) -> bytes:
    """打包为 ZIP：各子群文件夹在 ZIP 根下（官网同款）；可选 wrapping 前缀。

    只写入 ``specs`` 给出的子群文件，不读取、不混入 output_dir 中的其它文件。
    """
    wrap = safe_name(wrapping) if wrapping else ""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        if rendered_batch is not None:
            expected_specs = list(specs)
            entries = list(rendered_batch)
            if len(entries) != len(expected_specs) or any(
                rendered_spec is not expected
                for expected, (rendered_spec, _payloads) in zip(
                    expected_specs,
                    entries,
                    strict=True,
                )
            ):
                raise ValueError(
                    "pre-rendered ZIP payloads do not match the export specs"
                )
            used_folders: set[str] = set()
            used_paths: set[str] = set()
            for spec, payloads in entries:
                folder = str(
                    spec.folder_name or folder_label_for_subgroup(spec.subgroup)
                )
                if not folder or Path(folder).name != folder:
                    raise ValueError(
                        f"ZIP candidate folder must be one safe component: {folder!r}"
                    )
                folded_folder = folder.casefold()
                if folded_folder in used_folders:
                    raise ValueError(f"duplicate ZIP candidate folder: {folder!r}")
                used_folders.add(folded_folder)
                for filename, payload in payloads:
                    name = str(filename)
                    if not name or Path(name).name != name:
                        raise ValueError(
                            f"ZIP filename must be one safe component: {name!r}"
                        )
                    arcname = (
                        f"{wrap}/{folder}/{name}"
                        if wrap
                        else f"{folder}/{name}"
                    )
                    folded_path = arcname.casefold()
                    if folded_path in used_paths:
                        raise ValueError(f"duplicate ZIP member: {arcname!r}")
                    used_paths.add(folded_path)
                    zf.writestr(arcname, bytes(payload))
        else:
            used: set[str] = set()
            for spec in specs:
                requested = spec.folder_name or folder_label_for_subgroup(spec.subgroup)
                folder = _deduplicate_folder_name(requested, spec.subgroup, used)
                spec.folder_name = folder
                for fmt in formats:
                    fname = format_filename(folder, fmt)
                    arcname = (
                        f"{wrap}/{folder}/{fname}"
                        if wrap
                        else f"{folder}/{fname}"
                    )
                    zf.writestr(arcname, render_format(fmt, spec))
    return buf.getvalue()
