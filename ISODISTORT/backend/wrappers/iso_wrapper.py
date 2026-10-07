"""
iso（Isotropy 9.6.1）封装：k 点枚举、不可约表示、各向同性子群、畴、模式基矢

实测命令接口（见 resources/isobyu/iso 与 ISOTROPY 用户手册 isotropy_doc）：
    PAGE NOBREAK            关闭分页（否则大表会停在 “Enter RETURN to continue”）
    SCREEN 200              加宽屏幕列宽（默认 80 列会导致子群表折行）
    VALUE PARENT <sg>       选择母相空间群（数字或符号）
    VALUE KPOINT <label>    选择 k 点（Miller-Love 记号）
    VALUE KVALUE <n>,<v1>,.. 设置 k 点参数（n 为参数个数，如 "1,1/4"）
    VALUE IRREP <label>     选择不可约表示（Miller-Love 记号）
    VALUE DIRECTION <sym>   选择序参量方向（如 P1）
    VALUE WYCKOFF <letter>  选择 Wyckoff 位置（可多次）
    SHOW <flag>             控制 DISPLAY 输出内容
    DISPLAY KPOINT/IRREP/ISOTROPY/BUSH/PARENT
    QUIT

对应 ISODISTORT 官网阶段：
- 阶段二步骤4：枚举各向同性子群（DISPLAY ISOTROPY）
- 阶段三步步骤6/7：计算畸变模式基矢（DISPLAY BUSH + SHOW MODES）
- 阶段五步骤10：畴变体（SHOW DOMAIN）
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from fractions import Fraction
from math import gcd, isfinite, lcm

import numpy as np

from backend.models.iso_mode_models import (
    InvariantDirection,
    MacroscopicTensorBasis,
    MacroscopicTensorBlock,
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
    OrientedSiteCharacterTable,
    SiteIrrepCharacter,
    SiteOperation,
    SymbolicTensorComponent,
    microscopic_float_column_digest,
    validate_microscopic_mode_source,
)
from backend.tables.kpoints_official import official_special_k_coords
from backend.utils import (
    OutputParseError,
    WrapperRunError,
    WrapperTimeoutError,
    detect_blocked_generation,
    detect_missing_subgroup_db,
    get_config,
    parse_bush_table,
    parse_coords_token,
    parse_domain_table,
    parse_irrep_table,
    parse_kpoint_table,
    parse_subgroup_table,
)
from backend.utils.lattice import as_fraction, centering_primitive_matrix, multiply
from backend.utils.opd_format import (
    direction_coefficient_matrix,
    format_k_active,
    format_opd_line,
    format_opd_line_body,
    official_method1_fields,
    resolve_active_k_star_coordinates,
)
from backend.utils.schoenflies import hm_symbol
from backend.utils.text_parser import (
    parse_invariant_direction_table,
    parse_macroscopic_tensor_blocks,
    parse_microscopic_vector_blocks,
    parse_wyckoff_character_tables,
)

from .base_wrapper import BaseWrapper


@dataclass(frozen=True)
class _MicroscopicResolvedQuery:
    """One exact DISPLAY DISTORTION request and its returned source blocks."""

    invariant: InvariantDirection
    direction_selector: str
    blocks: tuple[MicroscopicVectorBlock, ...]

# ================================================================
# 数据模型
# ================================================================

@dataclass
class KPointInfo:
    """k 点信息（DISPLAY KPOINT）"""

    label: str                       # Miller-Love 记号，如 GM / DT
    coordinates: list[str]           # 坐标分量字符串，如 ["0","2a","0"]
    parameters: list[str]            # 自由参数字母，如 ["a"]
    is_special: bool                 # 无自由参数的特殊 k 点
    kovalev: str | None = None       # 官网 Kovalev 编号（如 "k14"），本地无数据时为 None

    def __post_init__(self) -> None:
        self.parameters = sorted(
            {c for c in self.coordinates if re.search(r"[a-zA-Z]", c)}
        )
        self.is_special = not self.parameters


@dataclass
class IrrepInfo:
    """不可约表示信息（DISPLAY IRREP）"""

    label: str                       # Miller-Love 记号，如 GM1+
    dimension: int                   # 维度
    active: bool = True              # Landau active 标记（仅用于参考）


@dataclass
class SubgroupInfo:
    """各向同性子群信息（DISPLAY ISOTROPY）"""

    index: int                       # 本地枚举序号（用户选择句柄，0 起）
    space_group_number: int          # 子群空间群号
    space_group_symbol: str = ""     # 子群空间群短符号
    subgroup_index: int = 0          # iso 输出 Index（相对母相的子群指数 = 畴数）
    size: int = 1                    # 子群原胞相对母相的大小 s
    is_maximal: bool = False         # 是否为 maximal 子群
    opd_symbol: str = ""             # 序参量方向符号（如 P1）
    opd_vector: list[float] = field(default_factory=list)   # 序参量方向向量
    basis_vectors: list[list[float]] = field(default_factory=list)  # 超胞基矢（母相格单位）
    origin: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])  # 超胞原点
    k_point_label: str = ""          # 产生该子群的 k 点
    irrep_label: str = ""            # 产生该子群的不可约表示
    k_parameters: list[str] = field(default_factory=list)  # k 点参数值（带参数 k 点）
    k_coordinates: list[str] = field(default_factory=list)  # DISPLAY KPOINT 坐标分量
    parent_sg: int = 0               # 母相空间群号（用于 k-star / 官网行格式）
    opd_dir_raw: str = ""            # iso 字母序参量方向，如 (a) / (a;a)
    basis_raw: str = ""              # iso 基矢原文，如 (1,0,0),(0,1,0),(0,0,1)
    origin_raw: str = ""             # iso 原点原文，如 (0,0,0)
    k_active_raw: str = ""           # 官网 k-active 列表（含前导空格）

    def describe(self) -> str:
        return (
            f"#{self.index}: SG {self.space_group_number} {self.space_group_symbol} "
            f"k={self.k_point_label} IR={self.irrep_label} OPD={self.opd_symbol}"
        )

    def official_fields(self) -> dict[str, str | int]:
        """Official Method 1 radio tokens as table columns (same values as ``opd_line``)."""
        return official_method1_fields(
            irrep_label=self.irrep_label,
            opd_symbol=self.opd_symbol,
            opd_dir_raw=self.opd_dir_raw,
            space_group_number=self.space_group_number,
            space_group_symbol=self.space_group_symbol,
            basis_raw=self.basis_raw,
            origin_raw=self.origin_raw,
            size=self.size,
            subgroup_index=self.subgroup_index,
            k_coordinates=self.k_coordinates,
            parent_sg=self.parent_sg or None,
            k_active_raw=self.k_active_raw or None,
            basis_vectors=self.basis_vectors,
            origin=self.origin,
        )

    def opd_line(self) -> str:
        """Official Method 1 radio-button line (visible text, no maximal asterisk)."""
        return format_opd_line(
            irrep_label=self.irrep_label,
            opd_symbol=self.opd_symbol,
            opd_dir_raw=self.opd_dir_raw,
            space_group_number=self.space_group_number,
            space_group_symbol=self.space_group_symbol,
            basis_raw=self.basis_raw,
            origin_raw=self.origin_raw,
            size=self.size,
            subgroup_index=self.subgroup_index,
            k_coordinates=self.k_coordinates,
            parent_sg=self.parent_sg or None,
            k_active_raw=self.k_active_raw or None,
            basis_vectors=self.basis_vectors,
            origin=self.origin,
        )

    def opd_line_body(self, *, pad_opd: bool | None = None) -> str:
        """CIF Distortion-page OPD comment (no irrep prefix).

        Method 1 pads the OPD token; Method 2 (parametric / listed IR OPD)
        uses a single trailing space after the OPD symbol.
        """
        if pad_opd is None:
            # Method 2 folders / parametric paths: short ``C1 (a,b)`` style.
            pad_opd = not bool(self.k_parameters)
        return format_opd_line_body(
            irrep_label=self.irrep_label,
            opd_symbol=self.opd_symbol,
            opd_dir_raw=self.opd_dir_raw,
            space_group_number=self.space_group_number,
            space_group_symbol=self.space_group_symbol,
            basis_raw=self.basis_raw,
            origin_raw=self.origin_raw,
            size=self.size,
            subgroup_index=self.subgroup_index,
            k_coordinates=self.k_coordinates,
            parent_sg=self.parent_sg or None,
            k_active_raw=self.k_active_raw or None,
            basis_vectors=self.basis_vectors,
            origin=self.origin,
            pad_opd=pad_opd,
        )


@dataclass
class DomainInfo:
    """畴变体信息（SHOW DOMAIN）"""

    domain_number: int               # 畴编号（1 起）
    generator: str = ""              # 域生成元，如 (C2y|0,0,0)
    space_group_number: int = 0
    space_group_symbol: str = ""
    subgroup_index: int = 0          # 母相中的子群指数（= 畴总数）
    opd_symbol: str = ""
    opd_vector: list[float] = field(default_factory=list)
    basis_vectors: list[list[float]] = field(default_factory=list)
    origin: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])


@dataclass
class BushMode:
    """模式基矢行（DISPLAY BUSH + SHOW MODES）"""

    irrep_label: str
    opd_symbol: str
    wyckoff_letter: str
    point: list[float]                       # 代表位点坐标（自由参数按 0）
    point_raw: list[str] = field(default_factory=list)  # 原始坐标 token（含 x/y/z 等参数）
    displacements: list[list[float]] = field(default_factory=list)  # 位移向量（可多个，对应模式分量）


@dataclass(frozen=True)
class MicroscopicDomainExtensionEvidence:
    """Immutable proof that sparse ISO rows were extended through a BUSH basis.

    ``MicroscopicColumnProvenance`` continues to describe only the direct
    ``DISPLAY DISTORTION`` source column.  This separate record describes the
    verified linear extension onto the complete ``DISPLAY BUSH`` point domain;
    it must never be interpreted as additional rows printed by ISO.
    """

    schema_version: str
    relation_status: str
    comparison_domain_kind: str
    comparison_point_count: int
    full_point_count: int
    mode_count: int
    reference_rank: int
    canonical_output_column_index: int
    comparison_point_digest: str
    full_point_digest: str
    transport_digest: str
    reference_basis_digest: str
    canonical_source_basis_digest: str
    canonical_source_column_token: str
    canonical_source_order_key: tuple[int, int, int]
    change_of_basis_digest: str
    extended_basis_digest: str
    extended_column_digest: str
    coefficient_column: tuple[float, ...]
    reference_condition_number: float
    canonical_condition_number: float
    change_of_basis_condition_number: float
    rank_relative_tolerance: float
    cartesian_absolute_tolerance_angstrom: float
    maximum_condition_number: float
    maximum_cartesian_residual_angstrom: float
    residual_tolerance_angstrom: float
    maximum_normalized_residual: float
    residual_norm: str
    matrix_convention: str

    def __post_init__(self) -> None:
        if self.schema_version != "bush-domain-extension-v1":
            raise ValueError("unknown microscopic domain-extension schema")
        if self.relation_status != "verified_bush_domain_extension":
            raise ValueError("unknown microscopic domain-extension relation")
        if self.comparison_domain_kind not in {
            "exact_common_points",
            "phase_proven_transported_points",
        }:
            raise ValueError("unknown microscopic domain-extension comparison domain")
        if (
            self.comparison_point_count <= 0
            or self.full_point_count < self.comparison_point_count
            or self.mode_count <= 0
            or self.reference_rank != self.mode_count
            or self.canonical_output_column_index < 0
            or self.canonical_output_column_index >= self.mode_count
            or len(self.coefficient_column) != self.mode_count
        ):
            raise ValueError("microscopic domain-extension dimensions are inconsistent")
        if (
            len(self.canonical_source_order_key) != 3
            or any(value < 0 for value in self.canonical_source_order_key)
        ):
            raise ValueError("microscopic domain extension has an invalid source index")
        digests = (
            self.comparison_point_digest,
            self.full_point_digest,
            self.transport_digest,
            self.reference_basis_digest,
            self.canonical_source_basis_digest,
            self.change_of_basis_digest,
            self.extended_basis_digest,
            self.extended_column_digest,
        )
        if any(re.fullmatch(r"[0-9a-f]{64}", value) is None for value in digests):
            raise ValueError("microscopic domain-extension digest is invalid")
        if not self.canonical_source_column_token.strip():
            raise ValueError("microscopic domain extension lacks a source column token")
        numeric = (
            *self.coefficient_column,
            self.reference_condition_number,
            self.canonical_condition_number,
            self.change_of_basis_condition_number,
            self.rank_relative_tolerance,
            self.cartesian_absolute_tolerance_angstrom,
            self.maximum_condition_number,
            self.maximum_cartesian_residual_angstrom,
            self.residual_tolerance_angstrom,
            self.maximum_normalized_residual,
        )
        if any(not isfinite(value) for value in numeric):
            raise ValueError("microscopic domain-extension evidence must be finite")
        if (
            self.reference_condition_number < 1.0
            or self.canonical_condition_number < 1.0
            or self.change_of_basis_condition_number < 1.0
            or self.rank_relative_tolerance <= 0.0
            or self.cartesian_absolute_tolerance_angstrom <= 0.0
            or self.maximum_condition_number < 1.0
            or self.reference_condition_number > self.maximum_condition_number
            or self.canonical_condition_number > self.maximum_condition_number
            or self.change_of_basis_condition_number > self.maximum_condition_number
            or self.maximum_cartesian_residual_angstrom < 0.0
            or self.residual_tolerance_angstrom <= 0.0
            or self.maximum_normalized_residual < 0.0
            or self.maximum_normalized_residual > 1.0
            or self.maximum_cartesian_residual_angstrom
            > self.residual_tolerance_angstrom
        ):
            raise ValueError("microscopic domain-extension numerical proof is invalid")
        if self.residual_norm != "cartesian_frobenius_angstrom":
            raise ValueError("unknown microscopic domain-extension residual norm")
        if self.matrix_convention != "R_common@T=C_common;C_full=R_full@T":
            raise ValueError("unknown microscopic domain-extension matrix convention")


@dataclass
class DistortionMode:
    """畸变模式（序参量 + 原子位移基矢）"""

    irrep_label: str                 # 不可约表示标号，如 GM4-
    dimension: int = 1               # 模式维度
    mode_type: str = "displacive"    # 类型: displacive/occupational/strain/magnetic/rotational
    basis_vectors: list[list[float]] = field(default_factory=list)  # 兼容旧接口
    wyckoff_site: str = ""           # 对应 Wyckoff 位置
    k_point_label: str = ""          # k 点
    opd_symbol: str = ""             # 序参量方向
    opd_dir_raw: str = ""            # exact invariant vector, e.g. ``(a;a)``
    bush_modes: list[BushMode] = field(default_factory=list)  # 原子级位移模式
    amplitude_key: str = ""          # unique export/GD key when one IR splits
    site_irrep: str = ""             # parent-site irrep, e.g. A2u / A1_1
    k_coords_label: str = ""         # ``0,0,1/6``
    opd_component: str = "a"         # ``a`` / ``b`` / …
    # Stable identity of one physical parent Wyckoff orbit.  The letter alone
    # is insufficient when a structure contains several independent sites on
    # the same Wyckoff position (for example several chemically or
    # parametrically distinct ``e`` orbits).
    wyckoff_orbit_id: str = ""
    # Audited scientific identity.  ``None`` is retained for legacy callers;
    # new construction paths use an explicit unresolved identity instead of
    # inventing a site irrep from vector direction or mode count.
    mode_identity: ModeIdentity | None = None
    # Exact, immutable evidence for an ISO microscopic source column.  Legacy
    # BUSH/SMODES columns have no such evidence and must retain ``None``.
    microscopic_provenance: MicroscopicColumnProvenance | None = None
    # Separate proof for rows reconstructed on a complete BUSH point domain.
    # ``None`` means ``bush_modes`` are the rows directly printed by the
    # microscopic query; a value never changes the direct source provenance.
    microscopic_domain_extension: MicroscopicDomainExtensionEvidence | None = None


def microscopic_subgroup_embedding_tokens(
    subgroup: SubgroupInfo,
) -> tuple[tuple[tuple[str, ...], ...], tuple[str, ...]]:
    """Return the exact child basis/origin used by the ISO query.

    The command renderer and immutable microscopic provenance share this one
    canonicalization path. Consequently the recorded context describes the
    actual query rather than a separately rounded reconstruction.
    """

    raw_vectors = re.findall(r"\(([^()]*)\)", str(subgroup.basis_raw or ""))
    if len(raw_vectors) == 3:
        basis_values: Sequence[Sequence[object]] = tuple(
            tuple(token.strip() for token in vector.split(","))
            for vector in raw_vectors
        )
    else:
        basis_values = tuple(tuple(row) for row in (subgroup.basis_vectors or ()))
    if len(basis_values) != 3 or any(len(row) != 3 for row in basis_values):
        raise ValueError("child embedding requires a 3x3 basis")
    basis = tuple(
        tuple(str(as_fraction(value)) for value in row) for row in basis_values
    )

    raw_origin = str(subgroup.origin_raw or "").strip()
    if raw_origin.startswith("(") and raw_origin.endswith(")"):
        origin_values: Sequence[object] = tuple(
            token.strip() for token in raw_origin[1:-1].split(",")
        )
    else:
        origin_values = tuple(subgroup.origin or ())
    if len(origin_values) != 3:
        raise ValueError("child embedding requires a three-component origin")
    origin = tuple(str(as_fraction(value)) for value in origin_values)
    return basis, origin


def microscopic_subgroup_translation_lattice(
    subgroup: SubgroupInfo,
) -> tuple[tuple[str, ...], ...]:
    """Return the exact child primitive translations in the parent frame.

    ISO prints the subgroup basis in its conventional setting.  Centered
    child groups therefore require the conventional-to-primitive centering
    matrix before a Bloch phase can be tested.  This is the same row-vector
    convention used by the search and superspace lattice code.
    """

    basis, _origin = microscopic_subgroup_embedding_tokens(subgroup)
    symbol = str(subgroup.space_group_symbol or "").replace(" ", "")
    if not symbol:
        symbol = hm_symbol(int(subgroup.space_group_number or 1)).replace(" ", "")
    centering = next(
        (char.upper() for char in symbol if char.upper() in "PABCIFR"),
        "P",
    )
    primitive = multiply(
        centering_primitive_matrix(centering),
        tuple(tuple(as_fraction(value) for value in row) for row in basis),
    )
    return tuple(tuple(str(value) for value in row) for row in primitive)


def microscopic_extension_digest(payload: object) -> str:
    """Hash extension evidence with exact rationals and binary floats preserved."""

    def _stable(value: object) -> object:
        if isinstance(value, Fraction):
            return {"fraction": [value.numerator, value.denominator]}
        if isinstance(value, float):
            if not isfinite(value):
                raise ValueError("extension digest cannot contain a non-finite float")
            return {"float_hex": value.hex()}
        if isinstance(value, Mapping):
            return {
                str(key): _stable(item)
                for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            }
        if isinstance(value, (list, tuple)):
            return [_stable(item) for item in value]
        if isinstance(value, (str, int, bool)) or value is None:
            return value
        raise TypeError(f"unsupported extension digest value {type(value).__name__}")

    encoded = json.dumps(
        _stable(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def microscopic_bush_column_digest(bush_modes: Sequence[BushMode]) -> str:
    """Digest one ordered full-domain column without conflating ISO provenance."""

    return microscopic_extension_digest([
        {
            "irrep": bush.irrep_label,
            "opd": bush.opd_symbol,
            "wyckoff": bush.wyckoff_letter,
            "point_raw": tuple(str(value) for value in bush.point_raw),
            "displacements": tuple(
                tuple(float(value) for value in vector)
                for vector in bush.displacements
            ),
        }
        for bush in bush_modes
    ])


def microscopic_bush_row_digest(bush_modes: Sequence[BushMode]) -> str:
    """Digest the ordered symbolic row domain exactly as ISO provenance does."""

    return microscopic_extension_digest([
        [str(value) for value in bush.point_raw] for bush in bush_modes
    ])


def microscopic_bush_source_float_column_digest(
    bush_modes: Sequence[BushMode],
) -> str:
    """Bind current direct-query float rows to their exact source column."""

    if not bush_modes or any(len(bush.displacements) != 1 for bush in bush_modes):
        raise ValueError("direct microscopic rows must contain one split column")
    return microscopic_float_column_digest([
        (bush.point_raw, bush.displacements[0]) for bush in bush_modes
    ])


def microscopic_provenance_basis_digest(
    provenances: Sequence[MicroscopicColumnProvenance],
) -> str:
    """Digest an ordered set of direct ISO source-column identities."""

    ordered: list[tuple[tuple[int, int, int], MicroscopicColumnProvenance]] = []
    for provenance in provenances:
        ordered.append((provenance.source_order_key, provenance))
    ordered.sort(key=lambda item: item[0])
    if len({item[0] for item in ordered}) != len(ordered):
        raise ValueError("microscopic source basis contains duplicate source columns")
    return microscopic_extension_digest([
        {
            "source_order_key": source_order_key,
            "exact_source_token": provenance.exact_source_token,
            "exact_row_digest": provenance.exact_row_digest,
            "exact_vector_digest": provenance.exact_vector_digest,
        }
        for source_order_key, provenance in ordered
    ])


def microscopic_source_basis_digest(modes: Sequence[DistortionMode]) -> str:
    provenances: list[MicroscopicColumnProvenance] = []
    for mode in modes:
        provenance = mode.microscopic_provenance
        if provenance is None:
            raise ValueError("microscopic source basis lacks exact provenance")
        provenances.append(provenance)
    return microscopic_provenance_basis_digest(provenances)


def validate_microscopic_domain_extension_bindings(
    bindings: Sequence[
        tuple[
            ModeIdentity,
            MicroscopicColumnProvenance,
            MicroscopicDomainExtensionEvidence | None,
        ]
    ],
) -> None:
    """Revalidate every complete per-orbit extension evidence group.

    This check is shared by the mapper and the public raw-mode I/O contract so
    neither boundary treats a frozen dataclass as an implicit trust marker.
    Direct ISO provenance is intentionally reusable across physical copies of
    one Wyckoff letter; the physical orbit is therefore part of the group key,
    but never part of the exact ISO source token.
    """

    groups: dict[
        tuple[str, str],
        list[
            tuple[
                ModeIdentity,
                MicroscopicColumnProvenance,
                MicroscopicDomainExtensionEvidence,
            ]
        ],
    ] = {}
    for identity, provenance, extension in bindings:
        if extension is None:
            continue
        orbit_id = str(identity.orbit_id).strip()
        if not orbit_id:
            raise ValueError("domain extension lacks a physical orbit identity")
        if extension.canonical_source_column_token != provenance.exact_source_token:
            raise ValueError("domain extension is bound to another source column")
        if extension.canonical_source_order_key != provenance.source_order_key:
            raise ValueError("domain extension is bound to another source index")
        groups.setdefault(
            (extension.canonical_source_basis_digest, orbit_id), [],
        ).append((identity, provenance, extension))

    for (expected_digest, _orbit_id), group in groups.items():
        resolved = [item[2] for item in group]
        mode_counts = {extension.mode_count for extension in resolved}
        if len(mode_counts) != 1 or len(group) != next(iter(mode_counts)):
            raise ValueError("domain-extension source basis is incomplete")
        common_facts = {
            (
                extension.schema_version,
                extension.relation_status,
                extension.comparison_domain_kind,
                extension.comparison_point_count,
                extension.full_point_count,
                extension.reference_rank,
                extension.comparison_point_digest,
                extension.full_point_digest,
                extension.transport_digest,
                extension.reference_basis_digest,
                extension.change_of_basis_digest,
                extension.extended_basis_digest,
                extension.reference_condition_number,
                extension.canonical_condition_number,
                extension.change_of_basis_condition_number,
                extension.rank_relative_tolerance,
                extension.cartesian_absolute_tolerance_angstrom,
                extension.maximum_condition_number,
                extension.maximum_cartesian_residual_angstrom,
                extension.residual_tolerance_angstrom,
                extension.maximum_normalized_residual,
                extension.residual_norm,
                extension.matrix_convention,
            )
            for extension in resolved
        }
        if len(common_facts) != 1:
            raise ValueError("domain-extension group evidence is inconsistent")
        if (
            microscopic_provenance_basis_digest([item[1] for item in group])
            != expected_digest
        ):
            raise ValueError("domain extension does not match its direct source basis")
        if len({extension.canonical_source_column_token for extension in resolved}) != len(group):
            raise ValueError("domain extension reuses a source column")
        by_output = {
            extension.canonical_output_column_index: extension
            for extension in resolved
        }
        if len(by_output) != len(group) or set(by_output) != set(range(len(group))):
            raise ValueError("domain extension loses canonical column order")
        ordered = [by_output[index] for index in range(len(group))]
        source_order_keys = [
            extension.canonical_source_order_key for extension in ordered
        ]
        if source_order_keys != sorted(source_order_keys):
            raise ValueError(
                "domain extension canonical output order differs from exact source order"
            )
        change_of_basis = [
            [
                ordered[column].coefficient_column[row]
                for column in range(len(group))
            ]
            for row in range(len(group))
        ]
        if (
            microscopic_extension_digest(change_of_basis)
            != ordered[0].change_of_basis_digest
        ):
            raise ValueError("domain extension change-of-basis digest mismatch")
        change_array = np.asarray(change_of_basis, dtype=float)
        singular_values = np.linalg.svd(change_array, compute_uv=False)
        largest = float(singular_values[0]) if singular_values.size else 0.0
        machine_epsilon = np.finfo(float).eps
        machine_threshold = (
            machine_epsilon * max(change_array.shape, default=1) * largest
        )
        numerical_rank = int(np.count_nonzero(
            singular_values > machine_threshold
        ))
        if numerical_rank != len(group) or singular_values.size < len(group):
            raise ValueError("domain extension change-of-basis rank loss")
        smallest = float(singular_values[len(group) - 1])
        rank_relative_tolerance = ordered[0].rank_relative_tolerance
        conditioning_threshold = max(
            machine_epsilon,
            rank_relative_tolerance * largest,
        )
        recalculated_condition = largest / smallest
        expected_maximum_condition = 1.0 / max(
            rank_relative_tolerance,
            machine_epsilon,
        )
        if (
            not np.isfinite(recalculated_condition)
            or smallest <= conditioning_threshold
            or recalculated_condition > ordered[0].maximum_condition_number
        ):
            raise ValueError("domain extension change-of-basis is ill-conditioned")
        comparison_tolerance = (
            64.0 * machine_epsilon * max(1.0, recalculated_condition)
        )
        if not np.isclose(
            ordered[0].change_of_basis_condition_number,
            recalculated_condition,
            rtol=64.0 * machine_epsilon,
            atol=comparison_tolerance,
        ):
            raise ValueError(
                "domain extension change-of-basis condition mismatch"
            )
        if not np.isclose(
            ordered[0].maximum_condition_number,
            expected_maximum_condition,
            rtol=64.0 * machine_epsilon,
            atol=64.0 * machine_epsilon * expected_maximum_condition,
        ):
            raise ValueError("domain extension condition limit mismatch")
        if (
            microscopic_extension_digest({
                "full_point_digest": ordered[0].full_point_digest,
                "ordered_column_digests": tuple(
                    extension.extended_column_digest for extension in ordered
                ),
            })
            != ordered[0].extended_basis_digest
        ):
            raise ValueError("domain extension full-basis digest mismatch")


def _component_label(index: int) -> str:
    """Return a stable spreadsheet-style component label (a..z, aa..).

    A DISPLAY BUSH displacement *column* is one independent vector in the
    subgroup-fixed displacement space.  The label is metadata only; the
    column itself, rather than this name, defines the scientific identity.
    """
    if index < 0:
        raise ValueError("component index must be non-negative")
    chars: list[str] = []
    value = index
    while True:
        chars.append(chr(ord("a") + value % 26))
        value = value // 26 - 1
        if value < 0:
            return "".join(reversed(chars))


# ================================================================
# iso 封装
# ================================================================

class IsoWrapper(BaseWrapper):
    """
    ISOTROPY (iso) 命令行程序封装

    iso 是交互式程序，通过关键字命令序列驱动（见模块 docstring）。
    本封装将 ISODISTORT 所需能力映射为真实的 iso 命令序列：

    1. list_k_points    枚举母相的全部 k 点（Method 2 的 k 点下拉列表）
    2. list_irreps      枚举指定 k 点下的不可约表示（Method 2 的 IR 下拉列表）
    3. list_subgroups   枚举（k 点, IR）对应的各向同性子群（Method 1/2/3 的核心）
    4. calc_distortion_modes  计算指定路径的畸变模式基矢（DISPLAY BUSH）
    5. get_domains      获取畴变体列表（Domains 输出）
    """

    def __init__(self) -> None:
        super().__init__()
        self.binary = self.cfg.iso_bin
        self._direction_symbol_cache: dict[tuple[int, str, str], str] = {}

    # ================================================================
    # 通用命令流构造
    # ================================================================

    @staticmethod
    def _session(text: str) -> str:
        """包装一个 iso 会话：关闭分页、加宽屏幕，最后退出。"""
        return f"PAGE NOBREAK\nSCREEN 200\n{text}QUIT\n"

    def _run_session(self, commands: str, timeout: float | None = None) -> str:
        """运行一次 iso 会话并返回标准输出。

        Args:
            commands: 会话命令流（不含 QUIT）
            timeout: 子进程超时秒数；None 使用配置默认值
        """
        stdout = self.run_stdin(self.binary, self._session(commands), timeout=timeout)
        # 非零返回码下 run_stdin 已抛 WrapperRunError
        return stdout

    @staticmethod
    def _embedding_commands(subgroup: SubgroupInfo) -> list[str]:
        """Render one exact ISO child embedding without decimal re-fitting."""

        if int(subgroup.space_group_number or 0) <= 0:
            raise ValueError("child embedding requires a positive space-group number")
        basis, origin = microscopic_subgroup_embedding_tokens(subgroup)
        basis_value = " ".join(",".join(row) for row in basis)
        origin_value = ",".join(origin)
        return [
            f"VALUE SUBGROUP {int(subgroup.space_group_number)}",
            f"VALUE BASIS {basis_value}",
            f"VALUE ORIGIN {origin_value}",
        ]

    @classmethod
    def _microscopic_query_commands(
        cls,
        parent_sg: int,
        subgroup: SubgroupInfo,
        irrep_labels: Sequence[str],
        wyckoff_letters: Sequence[str],
        *,
        direction_symbols: Mapping[str, str],
    ) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
        """Return normalized request order and the exact ISO command stream."""

        irreps = tuple(dict.fromkeys(
            str(value).strip() for value in irrep_labels if str(value).strip()
        ))
        letters = tuple(dict.fromkeys(
            str(value).strip() for value in wyckoff_letters if str(value).strip()
        ))
        if not irreps or not letters:
            return irreps, letters, ()
        commands = [
            f"VALUE PARENT {int(parent_sg)}",
            "SHOW IRREP",
            "SHOW WYCKOFF",
            "SHOW WYCKOFF IRREP",
            "SHOW MICROSCOPIC VECTOR",
        ]
        for irrep in irreps:
            direction = str(direction_symbols.get(irrep) or "").strip()
            if not direction:
                raise ValueError(
                    f"exact invariant direction selector is missing for {irrep}"
                )
            stem_match = re.match(r"^(?:m)?([A-Z]+)", irrep)
            if stem_match is None:
                raise ValueError(f"cannot resolve the k-point label for {irrep}")
            commands.append(f"VALUE KPOINT {stem_match.group(1)}")
            commands.append(f"VALUE IRREP {irrep}")
            if direction:
                commands.append(f"VALUE DIRECTION {direction}")
            commands.extend(cls._embedding_commands(subgroup))
            commands.append("VALUE WYCKOFF " + " ".join(letters))
            commands.append("DISPLAY DISTORTION")
        return irreps, letters, tuple(commands)

    @staticmethod
    def _normalized_direction_vector(value: str) -> str:
        """Keep ISO star-arm separators while ignoring presentation spacing."""

        return re.sub(r"\s+", "", str(value))

    @classmethod
    def _direction_components(cls, value: str) -> str:
        """Compare ISO outputs that inconsistently print `,` versus `;`."""

        return cls._normalized_direction_vector(value).replace(";", ",")

    @staticmethod
    def _direction_coefficient_matrix(
        value: str,
    ) -> tuple[tuple[Fraction, ...], ...]:
        """Parse one homogeneous symbolic OPD as an exact rational matrix.

        Rows are printed order-parameter components and columns are the free
        symbols (``a``, ``b``, ...).  A parameter rename or an invertible
        change of parameter basis leaves the column space unchanged, so the
        matrix is the appropriate object for comparing ISO direction forms.
        Unsupported nonlinear or inhomogeneous syntax fails closed.
        """

        return direction_coefficient_matrix(value)

    @staticmethod
    def _direction_subspace_projector(
        matrix: tuple[tuple[Fraction, ...], ...],
    ) -> tuple[tuple[Fraction, ...], ...]:
        """Return the exact orthogonal projector onto a matrix column space."""

        if not matrix or not matrix[0]:
            raise ValueError("direction coefficient matrix is empty")
        width = len(matrix[0])
        if any(len(row) != width for row in matrix):
            raise ValueError("direction coefficient matrix is ragged")

        # Row reduction identifies independent columns of the original matrix.
        reduced = [list(row) for row in matrix]
        pivot_columns: list[int] = []
        pivot_row = 0
        for column in range(width):
            pivot = next(
                (
                    row
                    for row in range(pivot_row, len(reduced))
                    if reduced[row][column]
                ),
                None,
            )
            if pivot is None:
                continue
            reduced[pivot_row], reduced[pivot] = reduced[pivot], reduced[pivot_row]
            scale = reduced[pivot_row][column]
            reduced[pivot_row] = [value / scale for value in reduced[pivot_row]]
            for row in range(len(reduced)):
                if row == pivot_row or not reduced[row][column]:
                    continue
                factor = reduced[row][column]
                reduced[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(
                        reduced[row], reduced[pivot_row], strict=True
                    )
                ]
            pivot_columns.append(column)
            pivot_row += 1
            if pivot_row == len(reduced):
                break
        if not pivot_columns:
            raise ValueError("direction coefficient matrix has rank zero")
        basis = [
            [row[column] for column in pivot_columns]
            for row in matrix
        ]
        rank = len(pivot_columns)
        gram = [
            [
                sum(
                    (basis[row][left] * basis[row][right] for row in range(len(basis))),
                    Fraction(0),
                )
                for right in range(rank)
            ]
            for left in range(rank)
        ]
        augmented = [
            [*row, *(Fraction(int(i == j)) for j in range(rank))]
            for i, row in enumerate(gram)
        ]
        for column in range(rank):
            pivot = next(
                (row for row in range(column, rank) if augmented[row][column]),
                None,
            )
            if pivot is None:
                raise ValueError("direction Gram matrix is singular")
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
            scale = augmented[column][column]
            augmented[column] = [value / scale for value in augmented[column]]
            for row in range(rank):
                if row == column or not augmented[row][column]:
                    continue
                factor = augmented[row][column]
                augmented[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(
                        augmented[row], augmented[column], strict=True
                    )
                ]
        inverse = [row[rank:] for row in augmented]
        return tuple(
            tuple(
                sum(
                    (
                        basis[left_row][left]
                        * inverse[left][right]
                        * basis[right_row][right]
                        for left in range(rank)
                        for right in range(rank)
                    ),
                    Fraction(0),
                )
                for right_row in range(len(basis))
            )
            for left_row in range(len(basis))
        )

    @classmethod
    def _direction_subspaces_same_component_order(
        cls, left: str, right: str
    ) -> bool:
        """Compare direction spaces while preserving the canonical IR axes."""

        try:
            return cls._direction_subspace_projector(
                cls._direction_coefficient_matrix(left)
            ) == cls._direction_subspace_projector(
                cls._direction_coefficient_matrix(right)
            )
        except (ValueError, ZeroDivisionError):
            return False

    @staticmethod
    def _direction_resolution_key(
        parent_sg: int, invariant: InvariantDirection
    ) -> tuple[int, str, str, int, int, tuple[str, ...]]:
        return (
            int(parent_sg),
            str(invariant.irrep_label).strip(),
            re.sub(r"\s+", "", str(invariant.direction_raw)),
            int(invariant.subgroup_number),
            int(invariant.size),
            tuple(invariant.k_parameters),
        )

    def _direction_symbol_for_invariant(
        self,
        parent_sg: int,
        invariant: InvariantDirection,
    ) -> str:
        """Resolve an exact-embedding direction vector to ISO's OPD token.

        This is a primary-candidate consistency gate. ``DISPLAY DIRECTION``
        supplies the exact embedding-fixed vector and the isotropy table names
        its selected OPD.  Their spaces must agree in the canonical component
        order; symmetry-equivalent component permutations are distinct domains
        and are deliberately rejected here.  Microscopic queries themselves
        use the exact ``VECTOR,...`` selector and need no secondary OPD token.
        """

        label = str(invariant.irrep_label).strip()
        cache = getattr(self, "_direction_symbol_cache", None)
        if cache is None:
            cache = {}
            self._direction_symbol_cache = cache
        key = self._direction_resolution_key(parent_sg, invariant)
        if key in cache:
            return cache[key]
        stem_match = re.match(r"^(?:m)?([A-Z]+)", label)
        if stem_match is None:
            raise ValueError(f"cannot resolve the k-point label for {label}")
        subgroup_kwargs = (
            {"k_parameters": tuple(invariant.k_parameters)}
            if invariant.k_parameters
            else {}
        )
        rows = self.list_subgroups(
            int(parent_sg),
            stem_match.group(1),
            label,
            **subgroup_kwargs,
        )
        matching_rows = [
            row
            for row in rows
            if int(getattr(row, "space_group_number", 0))
            == int(invariant.subgroup_number)
            and int(getattr(row, "size", 0)) == int(invariant.size)
        ]
        symbols = {
            str(row.opd_symbol).strip()
            for row in matching_rows
            if str(row.opd_symbol).strip()
            and self._direction_subspaces_same_component_order(
                invariant.direction_raw, row.opd_dir_raw
            )
        }
        if len(symbols) != 1:
            raise ValueError(
                "exact invariant direction does not resolve to one ISO OPD token: "
                f"irrep={label}, direction={invariant.direction_raw}, "
                f"tokens={sorted(symbols)}"
            )
        symbol = next(iter(symbols))
        cache[key] = symbol
        return symbol

    @classmethod
    def _direction_vector_selector(cls, value: str) -> str:
        """Serialize an exact homogeneous ISO direction as ``VECTOR,...``.

        ``DISPLAY DIRECTION`` uses parentheses and either commas or semicolons;
        ``VALUE DIRECTION`` accepts the same component expressions after the
        ``VECTOR`` keyword.  Parsing first rejects constants, nonlinear terms,
        functions, decimals, and other syntax that cannot be round-tripped as
        an exact rational homogeneous linear form.
        """

        raw_variables = re.findall(r"[A-Za-z]", str(value))
        by_casefold: dict[str, set[str]] = {}
        for variable in raw_variables:
            by_casefold.setdefault(variable.casefold(), set()).add(variable)
        if any(len(spellings) > 1 for spellings in by_casefold.values()):
            raise ValueError("direction vector mixes case-distinct parameter names")
        matrix = cls._direction_coefficient_matrix(value)
        parameter_count = len(matrix[0])
        if parameter_count > 26:
            raise ValueError("ISO VECTOR supports at most 26 named parameters")
        projector = cls._direction_subspace_projector(matrix)
        if sum((projector[index][index] for index in range(len(projector))), Fraction(0)) != parameter_count:
            raise ValueError("direction vector contains redundant free parameters")

        # Canonicalize every parameter column to primitive integers.  Clearing
        # rational denominators and dividing the integer gcd is an invertible
        # rescaling of that free parameter, so the exact direction space is
        # unchanged.  It also avoids decimals and ISO-version-dependent
        # rational-expression parsing.
        integer_columns: list[list[int]] = []
        for column in range(parameter_count):
            values = [row[column] for row in matrix]
            denominator = lcm(*(value.denominator for value in values))
            integers = [int(value * denominator) for value in values]
            common = gcd(*(abs(value) for value in integers if value))
            integers = [value // common for value in integers]
            first_nonzero = next(value for value in integers if value)
            if first_nonzero < 0:
                integers = [-value for value in integers]
            integer_columns.append(integers)
        normalized = tuple(
            tuple(Fraction(integer_columns[column][row]) for column in range(parameter_count))
            for row in range(len(matrix))
        )

        def term(coefficient: Fraction, symbol: str, *, first: bool) -> str:
            if not coefficient:
                return ""
            sign = "-" if coefficient < 0 else ("" if first else "+")
            magnitude = abs(coefficient)
            factor = "" if magnitude == 1 else str(magnitude.numerator)
            return f"{sign}{factor}{symbol}"

        components: list[str] = []
        for row in normalized:
            expression = ""
            for column, coefficient in enumerate(row):
                piece = term(
                    coefficient,
                    chr(ord("A") + column),
                    first=not expression,
                )
                expression += piece
            components.append(expression or "0")
        return "VECTOR," + ",".join(components)

    @classmethod
    def _direction_selector_rank(cls, selector: str) -> int | None:
        """Return the exact fixed-space rank encoded by a VECTOR selector."""

        prefix = "VECTOR,"
        compact = re.sub(r"\s+", "", str(selector))
        if not compact.upper().startswith(prefix):
            return None
        body = compact[len(prefix):]
        matrix = cls._direction_coefficient_matrix(f"({body})")
        projector = cls._direction_subspace_projector(matrix)
        rank = sum(
            (projector[index][index] for index in range(len(projector))),
            Fraction(0),
        )
        if rank.denominator != 1 or rank != len(matrix[0]):
            raise ValueError("VECTOR selector contains redundant parameters")
        return int(rank)

    def list_invariant_directions(
        self,
        parent_sg: int,
        subgroup: SubgroupInfo,
        *,
        primary_k_parameters: Sequence[str] | None = None,
    ) -> list[InvariantDirection]:
        """List every parent irrep with a fixed vector for the exact embedding."""

        commands = [f"VALUE PARENT {int(parent_sg)}"]
        if primary_k_parameters is not None:
            if not subgroup.k_point_label or not subgroup.irrep_label:
                raise ValueError(
                    "parametric invariant directions require a primary k point and irrep"
                )
            commands.extend((
                f"VALUE KPOINT {subgroup.k_point_label}",
                f"VALUE IRREP {subgroup.irrep_label}",
            ))
            if primary_k_parameters:
                commands.append(self._kvalue_command(primary_k_parameters))
        commands.extend(self._embedding_commands(subgroup))
        commands.extend((
            "SHOW IRREP",
            "SHOW KPOINT",
            "SHOW DIRECTION VEC",
            "DISPLAY DIRECTION",
        ))
        stdout = self._run_session("\n".join(commands) + "\n")
        try:
            rows = parse_invariant_direction_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError(
                "iso", f"解析精确子群嵌入的不变方向失败: {exc}"
            ) from exc
        return [InvariantDirection(**row) for row in rows]

    def get_macroscopic_tensor_basis(
        self,
        parent_sg: int,
        *,
        rank_signature: str = "[12]",
        subgroup: SubgroupInfo | None = None,
    ) -> MacroscopicTensorBasis:
        """Return ISO's canonical macroscopic tensor blocks in printed order.

        For symmetric rank-2 polar tensors the coefficient basis is engineering
        Voigt order ``(xx, yy, zz, yz, xz, xy)``.  If ``subgroup`` is supplied,
        the exact-embedding ``DISPLAY DIRECTION`` rows are attached as separate
        evidence.  They are not silently converted into component selections:
        a direction such as ``(a,a)`` is a linear combination, not two free
        columns, and must be applied and span-validated by the consumer.

        Rational coefficients are exact.  ISO's finite decimal coefficients
        remain :class:`Decimal` evidence and make the result ``partial`` until
        a tolerance-based whole-space validation succeeds.  Missing or
        malformed output returns an explicit ``unresolved`` result.
        """

        coefficient_basis = ("xx", "yy", "zz", "yz", "xz", "xy")
        directions: tuple[InvariantDirection, ...] = ()
        if subgroup is not None:
            try:
                directions = tuple(self.list_invariant_directions(parent_sg, subgroup))
            except (
                OutputParseError,
                WrapperRunError,
                WrapperTimeoutError,
                ValueError,
                IndexError,
            ):
                return MacroscopicTensorBasis.unresolved(
                    rank_signature=rank_signature,
                    coefficient_basis=coefficient_basis,
                    reason="iso_invariant_direction_output_unavailable",
                )
            if not directions:
                return MacroscopicTensorBasis.unresolved(
                    rank_signature=rank_signature,
                    coefficient_basis=coefficient_basis,
                    reason="iso_invariant_direction_output_empty",
                )

        commands = [
            f"VALUE PARENT {int(parent_sg)}",
            f"VALUE RANK {rank_signature}",
            "SHOW MACROSCOPIC",
            "SHOW IRREP",
            "DISPLAY DISTORTION",
        ]
        try:
            stdout = self._run_session("\n".join(commands) + "\n")
        except (WrapperRunError, WrapperTimeoutError):
            return MacroscopicTensorBasis.unresolved(
                rank_signature=rank_signature,
                coefficient_basis=coefficient_basis,
                reason="iso_macroscopic_tensor_output_unavailable",
                invariant_directions=directions,
            )
        try:
            rows = parse_macroscopic_tensor_blocks(stdout, coefficient_basis)
        except (ValueError, IndexError, ZeroDivisionError):
            return MacroscopicTensorBasis.unresolved(
                rank_signature=rank_signature,
                coefficient_basis=coefficient_basis,
                reason="iso_macroscopic_tensor_parse_error",
                invariant_directions=directions,
            )
        if not rows:
            return MacroscopicTensorBasis.unresolved(
                rank_signature=rank_signature,
                coefficient_basis=coefficient_basis,
                reason="iso_macroscopic_tensor_output_empty",
                invariant_directions=directions,
            )

        multiplicities: dict[str, int] = {}
        for row in rows:
            label = str(row["global_irrep"])
            multiplicities[label] = multiplicities.get(label, 0) + 1
        ordinals: dict[str, int] = {}
        blocks: list[MacroscopicTensorBlock] = []
        has_decimal = False
        for row in rows:
            label = str(row["global_irrep"])
            ordinals[label] = ordinals.get(label, 0) + 1
            copy_index = ordinals[label] if multiplicities[label] > 1 else None
            components = tuple(
                SymbolicTensorComponent(**item) for item in row["components"]
            )
            if any(component.quality == "unresolved" for component in components):
                return MacroscopicTensorBasis.unresolved(
                    rank_signature=rank_signature,
                    coefficient_basis=coefficient_basis,
                    reason="unsupported_iso_macroscopic_tensor_expression",
                    invariant_directions=directions,
                )
            has_decimal = has_decimal or any(
                component.quality == "printed_decimal" for component in components
            )
            blocks.append(MacroscopicTensorBlock(
                global_irrep=label,
                components=components,
                source_order=int(row["source_order"]),
                copy_index=copy_index,
            ))
        return MacroscopicTensorBasis(
            rank_signature=rank_signature,
            coefficient_basis=coefficient_basis,
            blocks=tuple(blocks),
            status="partial" if has_decimal else "verified",
            reason="iso_printed_decimal_coefficients" if has_decimal else None,
            invariant_directions=directions,
        )

    def get_microscopic_vector_blocks(
        self,
        parent_sg: int,
        subgroup: SubgroupInfo,
        irrep_labels: Sequence[str],
        wyckoff_letters: Sequence[str],
        *,
        direction_symbols: Mapping[str, str],
    ) -> list[MicroscopicVectorBlock]:
        """Return ISO's ordered, exact local-irrep displacement blocks.

        Every query is conditioned on the exact ``VECTOR,...`` form of the
        embedding-fixed direction and the child SG/basis/origin embedding.
        """

        requested_irreps, requested_letters, commands = self._microscopic_query_commands(
            parent_sg,
            subgroup,
            irrep_labels,
            wyckoff_letters,
            direction_symbols=direction_symbols,
        )
        if not commands:
            return []
        stdout = self._run_session("\n".join(commands) + "\n")
        return self._microscopic_vector_blocks_from_stdout(
            stdout,
            requested_irreps,
            requested_letters,
            direction_symbols=direction_symbols,
        )

    def _microscopic_vector_blocks_from_stdout(
        self,
        stdout: str,
        requested_irreps: Sequence[str],
        requested_letters: Sequence[str],
        *,
        direction_symbols: Mapping[str, str],
    ) -> list[MicroscopicVectorBlock]:
        """Parse and validate one or more complete microscopic query tables."""

        iso_error = re.search(
            r"(?im)^\s*\*+\s*(?:syntax\s+)?error\b[^\r\n]*",
            stdout,
        )
        if iso_error is not None:
            raise OutputParseError(
                "iso",
                "微观位移 VECTOR 查询被 ISO 拒绝: " + iso_error.group(0).strip(),
            )
        try:
            rows = parse_microscopic_vector_blocks(stdout)
        except (ValueError, IndexError, ZeroDivisionError) as exc:
            raise OutputParseError("iso", f"解析微观位移模式失败: {exc}") from exc
        allowed_irreps = set(requested_irreps)
        allowed_letters = set(requested_letters)
        unexpected_blocks = sorted({
            (str(row["global_irrep"]), str(row["wyckoff_letter"]))
            for row in rows
            if str(row["global_irrep"]) not in allowed_irreps
            or str(row["wyckoff_letter"]) not in allowed_letters
        })
        if unexpected_blocks:
            raise OutputParseError(
                "iso",
                "微观位移 VECTOR 查询返回未请求的块: "
                + ", ".join(
                    f"{irrep}/{letter}" for irrep, letter in unexpected_blocks
                ),
            )
        blocks: list[MicroscopicVectorBlock] = []
        for row in rows:
            expected_direction = str(
                direction_symbols.get(row["global_irrep"]) or ""
            ).strip()
            reported_direction = str(row["direction_symbol"] or "").strip()
            if (
                expected_direction
                and reported_direction
                and reported_direction != expected_direction
            ):
                raise OutputParseError(
                    "iso",
                    "微观位移模式的方向与查询不一致: "
                    f"{row['global_irrep']} {reported_direction} != "
                    f"{expected_direction}",
                )
            block_rows = tuple(
                MicroscopicVectorRow(
                    point_raw=tuple(item["point_raw"]),
                    displacements=tuple(item["displacements"]),
                    displacement_half_steps=tuple(
                        item.get("displacement_half_steps", ())
                    ),
                )
                for item in row["rows"]
            )
            blocks.append(MicroscopicVectorBlock(
                global_irrep=row["global_irrep"],
                direction_symbol=reported_direction or expected_direction,
                wyckoff_letter=row["wyckoff_letter"],
                site_irrep=row["site_irrep"],
                source_order=int(row["source_order"]),
                rows=block_rows,
            ))
        expected_ranks = {
            label: self._direction_selector_rank(selector)
            for label, selector in direction_symbols.items()
        }
        for block in blocks:
            expected_rank = expected_ranks.get(block.global_irrep)
            if expected_rank is not None and block.column_count != expected_rank:
                raise OutputParseError(
                    "iso",
                    "微观位移 VECTOR 查询返回错误列数: "
                    f"{block.global_irrep}/{block.wyckoff_letter}/"
                    f"{block.site_irrep} expected={expected_rank} "
                    f"actual={block.column_count}",
                )
        return blocks

    def _parametric_microscopic_query_results(
        self,
        parent_sg: int,
        subgroup: SubgroupInfo,
        invariants: Sequence[InvariantDirection],
        wyckoff_letters: Sequence[str],
    ) -> tuple[list[_MicroscopicResolvedQuery], tuple[str, ...]]:
        """Run repeated-k microscopic requests in one auditable ISO session."""

        letters = tuple(dict.fromkeys(
            str(value).strip() for value in wyckoff_letters if str(value).strip()
        ))
        if not invariants or not letters:
            return [], ()
        commands = [
            f"VALUE PARENT {int(parent_sg)}",
            "SHOW IRREP",
            "SHOW WYCKOFF",
            "SHOW WYCKOFF IRREP",
            "SHOW MICROSCOPIC VECTOR",
        ]
        selectors: list[str] = []
        for invariant in invariants:
            irrep = str(invariant.irrep_label).strip()
            stem_match = re.match(r"^(?:m)?([A-Z]+)", irrep)
            if stem_match is None:
                raise ValueError(f"cannot resolve the k-point label for {irrep}")
            selector = self._direction_vector_selector(invariant.direction_raw)
            selectors.append(selector)
            commands.extend((
                f"VALUE KPOINT {stem_match.group(1)}",
                f"VALUE IRREP {irrep}",
            ))
            if invariant.k_parameters:
                commands.append(self._kvalue_command(invariant.k_parameters))
            commands.append(f"VALUE DIRECTION {selector}")
            commands.extend(self._embedding_commands(subgroup))
            commands.append("VALUE WYCKOFF " + " ".join(letters))
            commands.append("DISPLAY DISTORTION")

        stdout = self._run_session("\n".join(commands) + "\n")
        all_blocks = self._microscopic_vector_blocks_from_stdout(
            stdout,
            tuple(dict.fromkeys(item.irrep_label for item in invariants)),
            letters,
            direction_symbols={},
        )
        # Every non-empty DISPLAY DISTORTION result starts with a starred
        # table header.  Continuation headers inserted by ISO pagination do
        # not contain stars.  Preserve those physical query boundaries: a
        # global-irrep run is not a query boundary because the same irrep can
        # recur at several exact KVALUEs, and intervening queries may have no
        # displacement columns for the requested Wyckoff sites.
        lines = stdout.splitlines()
        starts = [
            index for index, line in enumerate(lines)
            if re.match(r"^\s*\*+\s*Irrep\s+\(ML\)(?:\s|$)", line)
        ]
        table_blocks: list[list[MicroscopicVectorBlock]] = []
        if starts:
            for offset, start in enumerate(starts):
                end = starts[offset + 1] if offset + 1 < len(starts) else len(lines)
                payload = "\n".join(lines[start:end])
                parsed = self._microscopic_vector_blocks_from_stdout(
                    payload,
                    tuple(dict.fromkeys(item.irrep_label for item in invariants)),
                    letters,
                    direction_symbols={},
                )
                if parsed:
                    labels = {block.global_irrep for block in parsed}
                    if len(labels) != 1:
                        raise OutputParseError(
                            "iso",
                            "参数 k 的单次微观位移查询返回了多个 irrep: "
                            + ", ".join(sorted(labels)),
                        )
                    table_blocks.append(parsed)
        elif all_blocks:
            # Retain compatibility with captured/minimal ISO fixtures that
            # omit the decorative stars.  A single irrep is still unambiguous.
            labels = {block.global_irrep for block in all_blocks}
            if len(labels) != 1:
                raise OutputParseError(
                    "iso",
                    "参数 k 微观位移输出缺少查询边界，无法唯一归属多个 irrep",
                )
            table_blocks.append(all_blocks)

        def _block_signature(block: MicroscopicVectorBlock) -> tuple:
            return (
                block.global_irrep,
                block.direction_symbol,
                block.wyckoff_letter,
                block.site_irrep,
                block.rows,
            )

        segmented = [block for table in table_blocks for block in table]
        if (
            [_block_signature(block) for block in segmented]
            != [_block_signature(block) for block in all_blocks]
        ):
            raise OutputParseError(
                "iso",
                "参数 k 微观位移表边界不能无损还原完整 ISO 输出",
            )

        query_positions: dict[str, list[int]] = {}
        for index, invariant in enumerate(invariants):
            query_positions.setdefault(invariant.irrep_label, []).append(index)
        tables_by_label: dict[
            str, list[tuple[int, list[MicroscopicVectorBlock]]]
        ] = {}
        for table_order, table in enumerate(table_blocks):
            tables_by_label.setdefault(table[0].global_irrep, []).append(
                (table_order, table)
            )

        assigned: list[list[MicroscopicVectorBlock]] = [
            [] for _invariant in invariants
        ]
        mapped_table_order: list[tuple[int, int]] = []
        for label, positions in query_positions.items():
            tables = tables_by_label.get(label, [])
            if tables and len(tables) != len(positions):
                raise OutputParseError(
                    "iso",
                    "参数 k 的重复 irrep 查询含有无法唯一归属的空结果: "
                    f"{label} tables={len(tables)} queries={len(positions)}",
                )
            for position, (table_order, table) in zip(
                positions, tables, strict=False,
            ):
                assigned[position] = table
                mapped_table_order.append((table_order, position))
        mapped_positions = [
            position for _table_order, position in sorted(mapped_table_order)
        ]
        if mapped_positions != sorted(mapped_positions):
            raise OutputParseError(
                "iso",
                "参数 k 微观位移查询返回的 irrep 顺序与请求不一致",
            )

        resolved: list[_MicroscopicResolvedQuery] = []
        for invariant, selector, query_blocks in zip(
            invariants, selectors, assigned, strict=True,
        ):
            expected_rank = self._direction_selector_rank(selector)
            if any(
                block.global_irrep != invariant.irrep_label
                or (
                    expected_rank is not None
                    and block.column_count != expected_rank
                )
                for block in query_blocks
            ):
                raise OutputParseError(
                    "iso",
                    "参数 k 微观位移查询返回了错误的 irrep 或列数: "
                    f"expected={invariant.irrep_label}",
                )
            resolved.append(_MicroscopicResolvedQuery(
                invariant=invariant,
                direction_selector=selector,
                blocks=tuple(
                    replace(
                        block,
                        direction_symbol=block.direction_symbol or selector,
                    )
                    for block in query_blocks
                ),
            ))
        return resolved, tuple(commands)

    def get_wyckoff_character_tables(
        self,
        parent_sg: int,
        wyckoff_letters: Sequence[str],
    ) -> list[OrientedSiteCharacterTable]:
        """Return ISO's oriented site-stabilizer operations and characters."""

        letters = list(dict.fromkeys(str(value).strip() for value in wyckoff_letters))
        letters = [value for value in letters if value]
        if not letters:
            return []
        commands = [
            f"VALUE PARENT {int(parent_sg)}",
            "VALUE WYCKOFF " + " ".join(letters),
            "SHOW WYCKOFF POINTGROUP",
            "SHOW WYCKOFF ELEMENTS",
            "SHOW WYCKOFF CHARACTER",
            "DISPLAY PARENT",
        ]
        stdout = self._run_session("\n".join(commands) + "\n")
        try:
            rows = parse_wyckoff_character_tables(stdout)
        except (ValueError, IndexError, ZeroDivisionError) as exc:
            raise OutputParseError("iso", f"解析定向位点字符表失败: {exc}") from exc
        tables: list[OrientedSiteCharacterTable] = []
        for row in rows:
            tables.append(OrientedSiteCharacterTable(
                wyckoff_letter=row["wyckoff_letter"],
                point_group=row["point_group"],
                operations=tuple(SiteOperation(**item) for item in row["operations"]),
                irreps=tuple(SiteIrrepCharacter(**item) for item in row["irreps"]),
            ))
        return tables

    def _build_microscopic_modes_from_queries(
        self,
        parent_sg: int,
        subgroup: SubgroupInfo,
        queries: Sequence[_MicroscopicResolvedQuery],
        *,
        source_subgroup_context_kind: str,
        primary_direction_selector: str | None,
        query_commands: Sequence[str],
    ) -> list[DistortionMode]:
        """Materialize canonical columns while retaining each query's exact k."""

        direction_evidence = "".join(
            f"DIRECTION EVIDENCE {index} {query.invariant.irrep_label} "
            f"kparams={','.join(query.invariant.k_parameters)} "
            f"k={','.join(query.invariant.k_coordinates)} "
            f"raw={query.invariant.direction_raw} "
            f"selector={query.direction_selector} "
            f"normalized={self._normalized_direction_vector(query.invariant.direction_raw)}\n"
            for index, query in enumerate(queries)
        )
        query_digest = hashlib.sha256(
            (
                "iso-microscopic-query-v9\n"
                + f"SUBGROUP CONTEXT {source_subgroup_context_kind}\n"
                + f"PRIMARY DIRECTION SELECTOR {primary_direction_selector}\n"
                + direction_evidence
                + "\n".join(query_commands)
                + "\n"
            ).encode("utf-8")
        ).hexdigest()
        source_subgroup_basis, source_subgroup_origin = (
            microscopic_subgroup_embedding_tokens(subgroup)
        )
        source_subgroup_translation_lattice = (
            microscopic_subgroup_translation_lattice(subgroup)
        )
        listed_kpoints: list[KPointInfo] | None = None
        kpoint_listing_failed = False

        def _listed_special_k(label: str) -> tuple[str, ...]:
            nonlocal listed_kpoints, kpoint_listing_failed
            if listed_kpoints is None and not kpoint_listing_failed:
                try:
                    listed_kpoints = self.list_k_points(parent_sg)
                except (
                    OutputParseError,
                    WrapperRunError,
                    WrapperTimeoutError,
                    ValueError,
                    IndexError,
                ):
                    kpoint_listing_failed = True
                    listed_kpoints = []
            matches: list[tuple[str, ...]] = []
            for item in listed_kpoints or []:
                if item.label != label or item.parameters or not item.is_special:
                    continue
                coordinates = tuple(str(value).strip() for value in item.coordinates)
                if len(coordinates) != 3:
                    continue
                try:
                    exact = tuple(str(Fraction(value)) for value in coordinates)
                except (ValueError, ZeroDivisionError):
                    continue
                matches.append(exact)
            return matches[0] if len(matches) == 1 else ()

        modes: list[DistortionMode] = []
        for query_order, query in enumerate(queries):
            invariant = query.invariant
            blocks = list(query.blocks)
            multiplicities: dict[tuple[str, str, str], int] = {}
            for block in blocks:
                key = (block.global_irrep, block.wyckoff_letter, block.site_irrep)
                multiplicities[key] = multiplicities.get(key, 0) + 1
            ordinals: dict[tuple[str, str, str], int] = {}
            for block in blocks:
                copy_key = (block.global_irrep, block.wyckoff_letter, block.site_irrep)
                ordinals[copy_key] = ordinals.get(copy_key, 0) + 1
                copy_index = (
                    ordinals[copy_key] if multiplicities[copy_key] > 1 else None
                )
                display_irrep = (
                    f"{block.site_irrep}_{copy_index}"
                    if copy_index is not None
                    else block.site_irrep
                )
                stem_match = re.match(r"^(?:m)?([A-Z]+)", block.global_irrep)
                k_label = stem_match.group(1) if stem_match else subgroup.k_point_label
                k_tokens = tuple(str(value) for value in invariant.k_coordinates)
                if not k_tokens:
                    k_tokens = tuple(official_special_k_coords(parent_sg, k_label, []))
                if not k_tokens and block.global_irrep == subgroup.irrep_label:
                    k_tokens = tuple(str(value) for value in subgroup.k_coordinates)
                if not k_tokens:
                    k_tokens = _listed_special_k(k_label)
                try:
                    k_tokens = tuple(str(Fraction(value)) for value in k_tokens)
                except (ValueError, ZeroDivisionError):
                    k_tokens = ()
                k_joined = ",".join(k_tokens)
                for column in range(block.column_count):
                    component = _component_label(column)
                    active_k_resolution = (
                        resolve_active_k_star_coordinates(
                            k_tokens,
                            parent_sg,
                            block.direction_symbol,
                            column,
                            subgroup_translation_lattice=(
                                source_subgroup_translation_lattice
                            ),
                        )
                        if k_tokens
                        else None
                    )
                    provenance = MicroscopicColumnProvenance.from_exact_column(
                        query_digest=query_digest,
                        query_order=query_order,
                        query_irrep_label=block.global_irrep,
                        source_parent_sg=parent_sg,
                        source_k_coordinates=k_tokens,
                        source_active_k_coordinates=(
                            active_k_resolution.coordinates
                            if active_k_resolution is not None
                            else ()
                        ),
                        source_active_k_resolution_kind=(
                            active_k_resolution.resolution_kind
                            if active_k_resolution is not None
                            else None
                        ),
                        source_subgroup_space_group_number=subgroup.space_group_number,
                        source_subgroup_basis=source_subgroup_basis,
                        source_subgroup_origin=source_subgroup_origin,
                        source_subgroup_irrep_label=(
                            subgroup.irrep_label
                            if source_subgroup_context_kind == "single_irrep"
                            else None
                        ),
                        source_subgroup_opd_symbol=(
                            subgroup.opd_symbol
                            if source_subgroup_context_kind == "single_irrep"
                            else None
                        ),
                        source_subgroup_primary_direction_selector=(
                            primary_direction_selector
                        ),
                        source_subgroup_context_kind=(
                            source_subgroup_context_kind
                        ),
                        block=block,
                        column_index=column,
                    )
                    bushes: list[BushMode] = []
                    for row in block.rows:
                        point_token = "(" + ",".join(row.point_raw) + ")"
                        bushes.append(BushMode(
                            irrep_label=block.global_irrep,
                            opd_symbol=subgroup.opd_symbol,
                            wyckoff_letter=block.wyckoff_letter,
                            point=parse_coords_token(point_token),
                            point_raw=list(row.point_raw),
                            displacements=[
                                [float(value) for value in row.displacements[column]]
                            ],
                        ))
                    identity = ModeIdentity(
                        parent_sg=int(parent_sg),
                        global_irrep=block.global_irrep,
                        k_coordinates=k_tokens,
                        wyckoff_letter=block.wyckoff_letter,
                        site_irrep=block.site_irrep,
                        copy_index=copy_index,
                        component_index=column,
                        component_label=component,
                        source="iso_microscopic",
                        status="partial" if k_tokens else "unresolved",
                        reason=(
                            "awaiting_bush_or_smodes_subspace_validation"
                            if k_tokens
                            else "exact_special_k_identity_unresolved"
                        ),
                    )
                    modes.append(DistortionMode(
                        irrep_label=block.global_irrep,
                        dimension=block.column_count,
                        mode_type="displacive",
                        basis_vectors=[
                            list(bush.displacements[0]) for bush in bushes
                        ],
                        wyckoff_site=block.wyckoff_letter,
                        k_point_label=k_label,
                        opd_symbol=subgroup.opd_symbol,
                        opd_dir_raw=invariant.direction_raw,
                        bush_modes=bushes,
                        amplitude_key=(
                            f"{block.global_irrep}[{k_joined}]__"
                            f"{block.wyckoff_letter}__{display_irrep}({component})"
                        ),
                        site_irrep=display_irrep,
                        k_coords_label=k_joined,
                        opd_component=component,
                        mode_identity=identity,
                        microscopic_provenance=provenance,
                    ))
        return modes

    def calc_microscopic_distortion_modes(
        self,
        parent_sg: int,
        subgroup: SubgroupInfo,
        wyckoff_letters: Sequence[str],
        *,
        irrep_labels: Sequence[str] | None = None,
    ) -> list[DistortionMode]:
        """Build canonical mode columns directly from ISO microscopic blocks."""

        invariant_directions = self.list_invariant_directions(parent_sg, subgroup)
        invariants_by_irrep: dict[str, list[InvariantDirection]] = {}
        for invariant in invariant_directions:
            invariants_by_irrep.setdefault(invariant.irrep_label, []).append(invariant)

        primary_label = str(subgroup.irrep_label or "").strip()
        primary_invariants = invariants_by_irrep.get(primary_label, [])
        if len(primary_invariants) != 1:
            raise ValueError(
                "primary irrep does not have one exact invariant direction: "
                f"{primary_label or '<missing>'}"
            )
        primary_direction = self._direction_symbol_for_invariant(
            parent_sg, primary_invariants[0]
        )
        primary_direction_selector = self._direction_vector_selector(
            primary_invariants[0].direction_raw
        )
        if primary_direction != str(subgroup.opd_symbol or "").strip():
            raise ValueError(
                "primary invariant direction does not match the selected subgroup OPD: "
                f"{primary_label} {primary_direction} != {subgroup.opd_symbol}"
            )

        if irrep_labels is None:
            requested_labels = list(dict.fromkeys(
                invariant.irrep_label for invariant in invariant_directions
            ))
        else:
            requested_labels = list(dict.fromkeys(
                str(label).strip() for label in irrep_labels if str(label).strip()
            ))
        selected_invariants: list[InvariantDirection] = []
        for label in requested_labels:
            matches = invariants_by_irrep.get(label, [])
            if len(matches) != 1:
                raise ValueError(
                    "requested irrep does not have one exact invariant direction: "
                    f"{label}"
                )
            selected_invariants.append(matches[0])
        direction_vectors = {
            item.irrep_label: item.direction_raw for item in selected_invariants
        }
        direction_symbols = {
            label: self._direction_vector_selector(direction_vectors[label])
            for label in direction_vectors
        }
        requested_irreps, requested_letters, query_commands = (
            self._microscopic_query_commands(
                parent_sg,
                subgroup,
                requested_labels,
                wyckoff_letters,
                direction_symbols=direction_symbols,
            )
        )
        if not query_commands:
            return []
        blocks = self.get_microscopic_vector_blocks(
            parent_sg,
            subgroup,
            requested_irreps,
            requested_letters,
            direction_symbols=direction_symbols,
        )
        queries = [
            _MicroscopicResolvedQuery(
                invariant=invariant,
                direction_selector=direction_symbols[invariant.irrep_label],
                blocks=tuple(
                    block for block in blocks
                    if block.global_irrep == invariant.irrep_label
                ),
            )
            for invariant in selected_invariants
        ]
        return self._build_microscopic_modes_from_queries(
            parent_sg,
            subgroup,
            queries,
            source_subgroup_context_kind="single_irrep",
            primary_direction_selector=primary_direction_selector,
            query_commands=query_commands,
        )

    def calc_parametric_microscopic_distortion_modes(
        self,
        parent_sg: int,
        subgroup: SubgroupInfo,
        wyckoff_letters: Sequence[str],
        *,
        primary_k_parameters: Sequence[str],
        irrep_labels: Sequence[str] | None = None,
    ) -> list[DistortionMode]:
        """Build canonical columns for every commensurate parameter-k harmonic."""

        invariant_directions = self.list_invariant_directions(
            parent_sg,
            subgroup,
            primary_k_parameters=primary_k_parameters,
        )
        primary_k = tuple(
            str(Fraction(value)) for value in subgroup.k_coordinates
        )
        primary_matches = [
            invariant for invariant in invariant_directions
            if invariant.irrep_label == str(subgroup.irrep_label or "").strip()
            and tuple(invariant.k_coordinates) == primary_k
        ]
        if len(primary_matches) != 1:
            raise ValueError(
                "primary parameter-k irrep does not have one exact invariant direction: "
                f"{subgroup.irrep_label} k={','.join(primary_k)}"
            )
        primary = primary_matches[0]
        if not self._direction_subspaces_same_component_order(
            primary.direction_raw,
            subgroup.opd_dir_raw,
        ):
            raise ValueError(
                "primary parameter-k invariant direction disagrees with the "
                "selected subgroup direction"
            )
        primary_selector = self._direction_vector_selector(primary.direction_raw)

        requested = None if irrep_labels is None else {
            str(label).strip() for label in irrep_labels if str(label).strip()
        }
        selected = [
            invariant for invariant in invariant_directions
            if requested is None or invariant.irrep_label in requested
        ]
        if requested is not None:
            missing = requested.difference(
                invariant.irrep_label for invariant in selected
            )
            if missing:
                raise ValueError(
                    "requested parameter-k irreps have no exact invariant direction: "
                    + ", ".join(sorted(missing))
                )
        queries, query_commands = self._parametric_microscopic_query_results(
            parent_sg,
            subgroup,
            selected,
            wyckoff_letters,
        )
        return self._build_microscopic_modes_from_queries(
            parent_sg,
            subgroup,
            queries,
            source_subgroup_context_kind="single_irrep",
            primary_direction_selector=primary_selector,
            query_commands=query_commands,
        )

    def calc_embedding_microscopic_distortion_modes(
        self,
        parent_sg: int,
        subgroup: SubgroupInfo,
        wyckoff_letters: Sequence[str],
    ) -> list[DistortionMode]:
        """Build canonical columns for a route-less exact child fixed space.

        The selected Method 3 embedding has no distinguished primary order
        parameter.  ``DISPLAY DIRECTION`` supplies every parent irrep/k fixed
        by that exact child SG/basis/origin, and each row is queried through
        its own exact ``VECTOR`` selector.  No synthetic primary IR, OPD, or
        direction is introduced into the resulting provenance.
        """

        if str(getattr(subgroup, "_method3_route_resolution", "")) != (
            "exact_fixed_space"
        ):
            raise ValueError(
                "embedding microscopic modes require an exact-fixed-space route"
            )
        if any((
            str(subgroup.k_point_label or "").strip(),
            tuple(subgroup.k_coordinates or ()),
            tuple(subgroup.k_parameters or ()),
            str(subgroup.k_active_raw or "").strip(),
            str(subgroup.irrep_label or "").strip(),
            str(subgroup.opd_symbol or "").strip(),
            str(subgroup.opd_dir_raw or "").strip(),
            tuple(subgroup.opd_vector or ()),
        )):
            raise ValueError(
                "exact-fixed-space subgroup cannot claim one primary k, IR, OPD, "
                "or direction"
            )
        invariants = self.list_invariant_directions(parent_sg, subgroup)
        if not invariants:
            raise ValueError(
                "exact-fixed-space embedding has no ISO invariant directions"
            )
        queries, query_commands = self._parametric_microscopic_query_results(
            parent_sg,
            subgroup,
            invariants,
            wyckoff_letters,
        )
        return self._build_microscopic_modes_from_queries(
            parent_sg,
            subgroup,
            queries,
            source_subgroup_context_kind="exact_fixed_space",
            primary_direction_selector=None,
            query_commands=query_commands,
        )

    # ================================================================
    # k 点 / IR 枚举
    # ================================================================

    def list_k_points(self, parent_sg: int) -> list[KPointInfo]:
        """
        枚举母相空间群的全部 k 点（Method 2 的下拉列表数据源）。

        Args:
            parent_sg: 母相空间群号 (1-230)

        Returns:
            List[KPointInfo]
        """
        stdout = self._run_session(
            f"VALUE PARENT {parent_sg}\nSHOW KPOINT\nDISPLAY KPOINT\n"
        )
        try:
            rows = parse_kpoint_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析 k 点列表失败: {exc}") from exc
        return [KPointInfo(**row) for row in rows]

    def list_irreps(self, parent_sg: int, k_point: str,
                    k_parameters: Sequence[str] | None = None) -> list[IrrepInfo]:
        """
        枚举指定 k 点下的不可约表示（Method 2 的 IR 下拉列表数据源）。

        Args:
            parent_sg: 母相空间群号
            k_point: k 点标签（Miller-Love 记号）
            k_parameters: k 点参数值序列（如 ["1/4"]）；带参数 k 点必须提供

        Returns:
            List[IrrepInfo]
        """
        commands = [f"VALUE PARENT {parent_sg}", f"VALUE KPOINT {k_point}"]
        if k_parameters:
            commands.append(self._kvalue_command(k_parameters))
        commands += ["SHOW IRREP", "SHOW DIMENSION", "SHOW ACTIVE", "DISPLAY IRREP"]
        stdout = self._run_session("\n".join(commands) + "\n")
        try:
            rows = parse_irrep_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析不可约表示列表失败: {exc}") from exc
        return [IrrepInfo(**row) for row in rows]

    @staticmethod
    def _kvalue_command(k_parameters: Sequence[str]) -> str:
        """
        构造 VALUE KVALUE 命令。

        iso 语法：``VALUE KVALUE <参数个数>,<v1>,<v2>...``
        例如单个参数 1/4 -> "VALUE KVALUE 1,1/4"。
        """
        values = ",".join(str(v).strip() for v in k_parameters)
        return f"VALUE KVALUE {len(list(k_parameters))},{values}"

    # ================================================================
    # 子群枚举（Method 1/2/3 的核心）
    # ================================================================

    def list_subgroups(self, parent_sg: int,
                       k_point: str,
                       irrep_label: str,
                       k_parameters: Sequence[str] | None = None,
                       opd_symbol: str | None = None,
                       start_index: int = 0,
                       generate_if_missing: bool = False) -> list[SubgroupInfo]:
        """
        枚举指定 (k 点, IR) 下的各向同性子群。

        对应官网 Method 1（遍历全部特殊 k 点与 IR）与 Method 2
        （指定 k 点/IR）中“子群+序参量方向”列表。

        参数 k 点（如 LD、DT 等带 a/b/g 的点）说明：
        - DISPLAY ISOTROPY 流程要求 **先选择 IR 再设置 KVALUE**（实测 iso 9.6.1：
          若在 IR 之前设置 KVALUE，DISPLAY ISOTROPY 会报
          “parameters not selected for k vector”；注意与 list_irreps 的
          DISPLAY IRREP 流程顺序相反）；
        - 参数 k 点的子群数据库默认不存在，iso 会询问是否在线生成
          （对应官网 “Generate isotropy subgroups”，可能耗时数分钟到数小时）。
          当 generate_if_missing=True 时自动应答并等待生成完成（生成的数据库
          会保存到暂存目录，后续查询立即返回）；否则抛出 WrapperRunError。

        Args:
            parent_sg: 母相空间群号
            k_point: k 点标签
            irrep_label: 不可约表示标签
            k_parameters: k 点参数（带参数 k 点必须提供）
            opd_symbol: 若指定则只返回该序参量方向对应的子群
            start_index: 本地枚举序号的起始值（用于多 IR 合并时连续编号）
            generate_if_missing: 子群数据库缺失时是否自动在线生成
                （默认 False；生成可能耗时很长，请谨慎开启）

        Returns:
            List[SubgroupInfo]

        Raises:
            OutputParseError: 解析失败
            WrapperRunError: 需要在线生成子群数据库等无法自动完成的情况
        """
        commands = [
            f"VALUE PARENT {parent_sg}",
            f"VALUE KPOINT {k_point}",
            f"VALUE IRREP {irrep_label}",
        ]
        if k_parameters:
            # 实测 iso 9.6.1（DISPLAY ISOTROPY 流程）：KVALUE 必须在选择 IR
            # **之后**设置；若在 IR 之前设置，DISPLAY ISOTROPY 会报
            # “parameters not selected for k vector”。
            # （list_irreps 的 DISPLAY IRREP 流程则相反，需在 IR 选择前设置，
            # 两者不要混用。）
            commands.append(self._kvalue_command(k_parameters))
        if opd_symbol:
            commands.append(f"VALUE DIRECTION {opd_symbol}")
        commands += [
            "SHOW SUBGROUP", "SHOW INDEX", "SHOW SIZE", "SHOW DIRECTION VEC",
            "SHOW BASIS", "SHOW ORIGIN", "SHOW MAXIMAL",
            "DISPLAY ISOTROPY",
        ]

        extra_input = ""
        timeout = None
        if generate_if_missing:
            # iso 的提示为 “Should the data base be added? Enter RETURN to
            # continue. Enter any character to stop”（实测 iso 9.6.1）：
            # **空行（回车）= 触发生成；任何非空字符 = 停止生成**。
            # 生成完成后 iso 会**自动显示一次子群表**（“Adding data base...”
            # 之后紧跟表格），若再发 DISPLAY ISOTROPY 会导致同一张表
            # 重复输出、解析出重复子群，因此此处只发空行应答。
            # 输入流：... DISPLAY ISOTROPY -> <空行> -> QUIT
            extra_input = "\n"
            timeout = float(get_config().generation_timeout)
        else:
            # 应答 “Enter any character to stop”：停止生成、回到命令提示符，
            # 避免程序等待更多输入导致 EOF 崩溃
            extra_input = "q\n"

        stdout = self._run_session("\n".join(commands) + "\n" + extra_input,
                                   timeout=timeout)

        # 先解析子群表：若在线生成成功（或数据库已存在），直接返回结果
        try:
            rows = parse_subgroup_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析子群表失败: {exc}") from exc

        if not rows:
            # 无子群行：区分“需要生成”与“参数错误”
            if detect_missing_subgroup_db(stdout):
                summary = "\n".join(stdout.splitlines()[-25:]) if stdout else "(无输出)"
                if generate_if_missing:
                    # 已应答 y 触发生成，但仍无子群表：给出 iso 输出末尾以定位原因
                    raise WrapperRunError(
                        "iso", 1,
                        f"k 点 {k_point} 的子群数据库缺失，按 generate_if_missing=True 应答 y "
                        f"生成后仍无子群表（iso 可能无法本地生成，或生成耗时超限）。iso 输出末尾：\n"
                        f"{summary}",
                    )
                raise WrapperRunError(
                    "iso", 1,
                    f"k 点 {k_point} 的子群数据库在本地不存在，需要在线生成"
                    f"（官网对应 “Generate isotropy subgroups”，可能耗时数分钟到数小时）。"
                    f"如确认需要，请以 generate_if_missing=True 重试。",
                )
            if detect_blocked_generation(stdout):
                summary = "\n".join(stdout.splitlines()[-30:]) if stdout else "(无输出)"
                raise WrapperRunError(
                    "iso", 1,
                    f"k 点 {k_point} 的参数未正确指定（iso 报 “parameters not selected for k vector”）。"
                    f"当前 KVALUE 命令为 {self._kvalue_command(k_parameters) if k_parameters else '(未设置)'}，"
                    f"k_parameters={k_parameters}。iso 输出末尾：\n{summary}",
                )

        subgroups: list[SubgroupInfo] = []
        for i, row in enumerate(rows):
            subgroups.append(self._subgroup_from_row(
                row,
                index=start_index + i,
                k_point=k_point,
                irrep_label=irrep_label,
                k_parameters=list(k_parameters) if k_parameters else [],
                k_coordinates=official_special_k_coords(
                    parent_sg, k_point, [], k_parameters
                ),
                parent_sg=parent_sg,
            ))
        return subgroups

    @staticmethod
    def _include_irrep(ir: IrrepInfo, distortion_types) -> bool:
        """Magnetic (m*) irreps stay out of the default Method 1/3 enumeration."""
        if not ir.label.startswith("m"):
            return True
        if distortion_types is None:
            return False
        names = (
            distortion_types
            if isinstance(distortion_types, (list, tuple, set))
            else [distortion_types]
        )
        return any(str(x).lower() == "magnetic" for x in names)

    @staticmethod
    def _subgroup_from_row(row: dict,
                           *,
                           index: int,
                           k_point: str,
                           irrep_label: str,
                           k_parameters: list[str] | None = None,
                           k_coordinates: list[str] | None = None,
                           parent_sg: int = 0) -> SubgroupInfo:
        coords = list(k_coordinates or [])
        dir_raw = row.get("opd_dir_raw") or ""
        basis_raw = row.get("basis_raw") or ""
        origin_raw = row.get("origin_raw") or ""
        k_active = format_k_active(
            dir_raw, coords or ["0", "0", "0"], parent_sg or None, None,
        )
        return SubgroupInfo(
            index=index,
            space_group_number=row["space_group_number"],
            space_group_symbol=row["space_group_symbol"],
            subgroup_index=row["subgroup_index"],
            size=row["size"],
            is_maximal=row["is_maximal"],
            opd_symbol=row["opd_symbol"],
            opd_vector=row["opd_vector"],
            basis_vectors=row["basis_vectors"],
            origin=row["origin"],
            k_point_label=k_point,
            irrep_label=irrep_label,
            k_parameters=list(k_parameters) if k_parameters else [],
            k_coordinates=coords,
            parent_sg=parent_sg,
            opd_dir_raw=dir_raw,
            basis_raw=basis_raw,
            origin_raw=origin_raw,
            k_active_raw=k_active,
        )

    def enumerate_all_special_subgroups(self, parent_sg: int,
                                        distortion_types=None) -> list[SubgroupInfo]:
        """
        枚举母相的全部特殊 k 点 × 全部 IR 的各向同性子群。

        对应官网 Method 1（“Search over all special k points”）的完整候选集。
        只包含无自由参数的特殊 k 点（有参数 k 点属于 Method 2 范畴）。

        实现：每个特殊 k 点一个 iso 会话，用 ``DISPLAY SETTING`` 的输出行
        （“Current setting is International ...”）作为各 IR 子群表之间的
        分隔标记，大幅减少 WSL 进程启动次数。

        Args:
            parent_sg: 母相空间群号
            distortion_types: 畸变类型。默认不含 magnetic，带 ``m`` 前缀的
                IR 不进入枚举（与官网默认 Types 一致）。

        Returns:
            List[SubgroupInfo]
        """
        subgroups: list[SubgroupInfo] = []
        kpoints = self.list_k_points(parent_sg)
        for kp in kpoints:
            if not kp.is_special:
                continue
            try:
                irreps = self.list_irreps(parent_sg, kp.label)
            except WrapperRunError:
                continue
            irreps = [ir for ir in irreps if self._include_irrep(ir, distortion_types)]
            if not irreps:
                continue

            commands = [
                f"VALUE PARENT {parent_sg}",
                f"VALUE KPOINT {kp.label}",
                "SHOW SUBGROUP", "SHOW INDEX", "SHOW SIZE",
                "SHOW DIRECTION VEC", "SHOW BASIS", "SHOW ORIGIN", "SHOW MAXIMAL",
            ]
            for ir in irreps:
                commands.append(f"VALUE IRREP {ir.label}")
                commands.append("DISPLAY ISOTROPY")
                commands.append("DISPLAY SETTING")
            stdout = self._run_session("\n".join(commands) + "\n")

            # 分隔符：每个 IR 的子群表后紧跟 DISPLAY SETTING 输出
            chunks = stdout.split("Current setting is International")
            # chunks[0] = 程序横幅（含起始的 setting 行之前部分）；
            # chunks[1:] 依次对应 irreps[0], irreps[1], ...
            for ir, chunk in zip(irreps, chunks[1:], strict=False):
                try:
                    rows = parse_subgroup_table(chunk)
                except (ValueError, IndexError) as exc:
                    raise OutputParseError(
                        "iso", f"解析 {kp.label}/{ir.label} 子群表失败: {exc}"
                    ) from exc
                for row in rows:
                    subgroups.append(self._subgroup_from_row(
                        row,
                        index=len(subgroups),
                        k_point=kp.label,
                        irrep_label=ir.label,
                        k_coordinates=official_special_k_coords(
                            parent_sg, kp.label, kp.coordinates,
                        ),
                        parent_sg=parent_sg,
                    ))
        return subgroups

    # ================================================================
    # 模式基矢（Method 2 的 Distortion Page 数据源）
    # ================================================================

    def calc_distortion_modes(self, parent_sg: int,
                              subgroup: SubgroupInfo,
                              wyckoff_letters: Sequence[str]) -> list[DistortionMode]:
        """
        计算指定子群路径下的畸变模式基矢（DISPLAY BUSH + SHOW MODES）。

        Args:
            parent_sg: 母相空间群号
            subgroup: 目标子群（须含 k_point_label / irrep_label / opd_symbol）
            wyckoff_letters: 母相结构中各 Wyckoff 位置字母（来自 findsym）

        Returns:
            List[DistortionMode]：每个模式含 BushMode 原子位移基矢
        """
        if not wyckoff_letters:
            raise WrapperRunError("iso", 1, "未提供任何 Wyckoff 位置，无法计算模式。")

        if subgroup.k_parameters:
            # iso 的 DISPLAY BUSH 仅支持对称 k 点（“Selected irrep must belong
            # to a k point of symmetry”）；参数 k 点（如 LD/DT）的模式计算
            # 依赖官网的 (3+d) 维超空间机制，本地二进制无法完成。
            raise WrapperRunError(
                "iso", 1,
                f"k 点 {subgroup.k_point_label}（参数 {subgroup.k_parameters}）为参数"
                f"（非对称）k 点：iso 二进制只能枚举其子群，无法计算原子位移模式；"
                f"官网对该场景使用 (3+d) 维超空间机制，本地暂不支持。",
            )

        commands = [
            f"VALUE PARENT {parent_sg}",
            f"VALUE KPOINT {subgroup.k_point_label}",
        ]
        commands.append(f"VALUE IRREP {subgroup.irrep_label}")
        if subgroup.k_parameters:
            # DISPLAY BUSH 流程：带参数 k 点须在选 IR 之后设置 KVALUE
            # （与 list_subgroups 的 DISPLAY ISOTROPY 流程一致；当前参数 k 点
            # 已在上方提前报错，此分支仅作防御性保留）
            commands.append(self._kvalue_command(subgroup.k_parameters))
        commands.append(f"VALUE DIRECTION {subgroup.opd_symbol}")
        commands += [
            "SHOW MODES", "SHOW MICROSCOPIC",
            "SHOW SUBGROUP", "SHOW DIRECTION VEC", "SHOW INDEX",
            "SHOW BASIS", "SHOW ORIGIN",
        ]
        # VALUE WYCKOFF replaces the previously selected orbit; it is not an
        # accumulating selector.  Emit one BUSH table per requested orbit in
        # the same ISO session and merge the parsed rows below.  Sending
        # ``VALUE WYCKOFF d`` followed by ``VALUE WYCKOFF e`` and only one
        # DISPLAY BUSH silently computed e alone, which removed complete-mode
        # copies on d (for example EuAl4 GM5+ C1: 2 modes instead of 6).
        for letter in dict.fromkeys(str(value).strip() for value in wyckoff_letters):
            if not letter:
                continue
            commands.extend((f"VALUE WYCKOFF {letter}", "DISPLAY BUSH"))
        stdout = self._run_session("\n".join(commands) + "\n")

        if detect_blocked_generation(stdout):
            raise WrapperRunError(
                "iso", 1,
                f"k 点 {subgroup.k_point_label} 的参数未正确指定，无法计算模式；"
                f"请确认子群包含有效的 k 参数（{subgroup.k_parameters}）。",
            )

        try:
            rows = parse_bush_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析模式基矢表失败: {exc}") from exc

        # DISPLAY BUSH 的每个 displacement 列是固定子空间
        # V^H = {u | D(h)u=u, h in H} 的一条独立基矢；续行则是同一
        # 基矢在其他代表原子上的分量。因此必须按
        # (IR, OPD, Wyckoff orbit, column) 拆成模式，不能把列相加，
        # 也不能把续行当成新模式。
        grouped: dict[tuple[str, str, str], list[BushMode]] = {}
        for row in rows:
            bush = BushMode(**row)
            key = (bush.irrep_label, bush.opd_symbol, bush.wyckoff_letter)
            grouped.setdefault(key, []).append(bush)

        modes: list[DistortionMode] = []
        for (ir, opd, letter), bushes in grouped.items():
            dimension = max((len(b.displacements) for b in bushes), default=0)
            if dimension <= 0:
                continue
            k_stem = re.match(r"^([A-Z]+)", ir)
            k_label = k_stem.group(1) if k_stem else subgroup.k_point_label
            k_coords = ",".join(str(value) for value in subgroup.k_coordinates)
            for column in range(dimension):
                component_rows: list[BushMode] = []
                for bush in bushes:
                    if column >= len(bush.displacements):
                        continue
                    component_rows.append(BushMode(
                        irrep_label=bush.irrep_label,
                        opd_symbol=bush.opd_symbol,
                        wyckoff_letter=bush.wyckoff_letter,
                        point=list(bush.point),
                        point_raw=list(bush.point_raw),
                        displacements=[list(bush.displacements[column])],
                    ))
                if not component_rows:
                    continue
                component = _component_label(column)
                # Prefixing with ``IR__`` preserves the public API's ability
                # to address/sum all copies by the bare irrep label.
                amplitude_key = f"{ir}__{letter}__{opd}__{component}"
                identity = ModeIdentity.unresolved(
                    parent_sg=int(parent_sg),
                    global_irrep=ir,
                    k_coordinates=tuple(str(value) for value in subgroup.k_coordinates),
                    wyckoff_letter=letter,
                    reason="display_bush_has_no_site_irrep_identity",
                )
                modes.append(DistortionMode(
                    irrep_label=ir,
                    dimension=dimension,
                    mode_type="displacive",
                    basis_vectors=[
                        list(bush.displacements[0]) for bush in component_rows
                    ],
                    wyckoff_site=letter,
                    k_point_label=k_label,
                    opd_symbol=opd,
                    bush_modes=component_rows,
                    amplitude_key=amplitude_key,
                    k_coords_label=k_coords,
                    opd_component=component,
                    mode_identity=identity,
                ))
        return modes

    # ================================================================
    # 畴变体
    # ================================================================

    def get_domains(self, parent_sg: int, subgroup: SubgroupInfo,
                    k_parameters: Sequence[str] | None = None) -> list[DomainInfo]:
        """
        获取指定子群的畴变体列表（SHOW DOMAIN）。

        畴总数等于子群在母相中的指数（subgroup_index），
        与官网 “Domains” 输出一致：每个畴含生成元、空间群、基矢与原点。

        Args:
            parent_sg: 母相空间群号
            subgroup: 目标子群（须含 k/IR/OPD）
            k_parameters: k 点参数

        Returns:
            List[DomainInfo]
        """
        commands = [
            f"VALUE PARENT {parent_sg}",
            f"VALUE KPOINT {subgroup.k_point_label}",
        ]
        if k_parameters is None:
            k_parameters = subgroup.k_parameters
        commands.append(f"VALUE IRREP {subgroup.irrep_label}")
        if k_parameters:
            # 顺序：先选 IR 再设 KVALUE（DISPLAY ISOTROPY 流程实测要求，
            # 见 list_subgroups 注释）
            commands.append(self._kvalue_command(k_parameters))
        commands += [
            f"VALUE DIRECTION {subgroup.opd_symbol}",
            "SHOW DOMAIN", "SHOW DOMAIN GENERATORS",
            "SHOW SUBGROUP", "SHOW INDEX", "SHOW BASIS", "SHOW ORIGIN",
            "SHOW DIRECTION VEC",
            "DISPLAY ISOTROPY",
        ]
        stdout = self._run_session("\n".join(commands) + "\n")
        try:
            rows = parse_domain_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析畴表失败: {exc}") from exc
        return [DomainInfo(**row) for row in rows]


def validate_microscopic_distortion_mode_source(
    mode: DistortionMode,
    subgroup: SubgroupInfo,
) -> None:
    """Bind one mutable working mode to its exact query and subgroup context."""

    identity = mode.mode_identity
    provenance = mode.microscopic_provenance
    if identity is None or provenance is None:
        raise ValueError("microscopic mode lacks verified source-column identity")
    validate_microscopic_mode_source(
        identity,
        provenance,
        mode_id=identity.stable_token,
    )
    if (
        mode.irrep_label != identity.global_irrep
        or mode.wyckoff_site != identity.wyckoff_letter
        or mode.wyckoff_orbit_id != identity.orbit_id
        or mode.site_irrep != identity.display_site_irrep
        or mode.opd_component != identity.component_label
        or mode.dimension != provenance.source_block_column_count
        or tuple(
            token.strip()
            for token in str(mode.k_coords_label).split(",")
            if token.strip()
        )
        != tuple(identity.k_coordinates)
    ):
        raise ValueError("mutable mode routing fields differ from ModeIdentity")
    try:
        direction_selector = IsoWrapper._direction_vector_selector(mode.opd_dir_raw)
    except ValueError as exc:
        raise ValueError("microscopic mode has an invalid exact OPD direction") from exc
    if direction_selector != provenance.source_direction_symbol:
        raise ValueError("microscopic mode direction differs from source VECTOR selector")
    context_kind = provenance.source_subgroup_context_kind
    if context_kind == "single_irrep":
        try:
            primary_direction_selector = IsoWrapper._direction_vector_selector(
                subgroup.opd_dir_raw
            )
        except ValueError as exc:
            raise ValueError(
                "target subgroup has an invalid exact OPD direction"
            ) from exc
        if (
            primary_direction_selector
            != provenance.source_subgroup_primary_direction_selector
        ):
            raise ValueError(
                "microscopic mode primary direction differs from target subgroup"
            )
        if mode.opd_symbol != provenance.source_subgroup_opd_symbol:
            raise ValueError(
                "microscopic mode OPD differs from source subgroup context"
            )
    elif (
        str(getattr(subgroup, "_method3_route_resolution", ""))
        != "exact_fixed_space"
        or any((
            str(subgroup.k_point_label or "").strip(),
            tuple(subgroup.k_coordinates or ()),
            tuple(subgroup.k_parameters or ()),
            str(subgroup.k_active_raw or "").strip(),
            str(subgroup.irrep_label or "").strip(),
            str(subgroup.opd_symbol or "").strip(),
            str(subgroup.opd_dir_raw or "").strip(),
            tuple(subgroup.opd_vector or ()),
            str(mode.opd_symbol or "").strip(),
        ))
    ):
        raise ValueError(
            "exact-fixed-space microscopic mode claims a primary subgroup route"
        )
    if any(
        bush.irrep_label != identity.global_irrep
        or bush.wyckoff_letter != identity.wyckoff_letter
        or (
            context_kind == "single_irrep"
            and bush.opd_symbol != provenance.source_subgroup_opd_symbol
        )
        or (
            context_kind == "exact_fixed_space"
            and str(bush.opd_symbol or "").strip()
        )
        for bush in mode.bush_modes
    ):
        raise ValueError("microscopic row routing fields differ from source context")
    basis, origin = microscopic_subgroup_embedding_tokens(subgroup)
    actual_context = (
        int(subgroup.parent_sg),
        int(subgroup.space_group_number),
        basis,
        origin,
        (
            str(subgroup.irrep_label).strip()
            if context_kind == "single_irrep" else None
        ),
        (
            str(subgroup.opd_symbol).strip()
            if context_kind == "single_irrep" else None
        ),
        context_kind,
    )
    source_context = (
        provenance.source_parent_sg,
        provenance.source_subgroup_space_group_number,
        provenance.source_subgroup_basis,
        provenance.source_subgroup_origin,
        provenance.source_subgroup_irrep_label,
        provenance.source_subgroup_opd_symbol,
        provenance.source_subgroup_context_kind,
    )
    if actual_context != source_context:
        raise ValueError("microscopic mode query context differs from target subgroup")
