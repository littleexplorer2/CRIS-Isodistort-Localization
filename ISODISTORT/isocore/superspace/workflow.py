"""(3+nmod) 端到端工作流、JSON/YAML 序列化、CLI 任务。"""
from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from ..io.result_serializer import ResultSerializer
from ..utils import DimensionMismatchError, InputError
from ..utils.text_parser import parse_fraction
from .constants import SERIAL_KIND, SERIAL_KIND_3P1, eps
from .group import SuperspaceGroup
from .kvector import KsVector, SuperspaceLattice, assert_ks_compatible, default_lattice_matrix
from .representation import (
    Irrep,
    OpdOperator,
    SuperspaceMode,
    build_irreps,
    filter_allowed_opds,
    generate_modes,
    generate_opds,
    little_group,
)
from .validate import dim_3pd, validate_nmod, validate_space_group_number, validate_vector_dim


def _parse_number_list(text: str) -> list[float]:
    parts = [p.strip() for p in text.replace(";", ",").split(",") if p.strip()]
    if not parts:
        raise InputError("empty numeric list")
    return [parse_fraction(p) for p in parts]


def parse_ks_text(text: str, nmod: int) -> list[float]:
    values = _parse_number_list(text)
    expected = dim_3pd(nmod)
    if len(values) == 3 and nmod > 0:
        extra = [0.0] * nmod
        extra[0] = 1.0
        values = values + extra
    validate_vector_dim(values, expected, name="k_s")
    return values


def parse_q_vectors_text(text: str | None, nmod: int) -> list[list[float]] | None:
    if not text or not str(text).strip():
        return None
    nmod = validate_nmod(nmod)
    chunks = [c.strip() for c in str(text).split(";") if c.strip()]
    vecs = [_parse_number_list(c) for c in chunks]
    for i, v in enumerate(vecs):
        validate_vector_dim(v, 3, name=f"q_vectors[{i}]")
    if len(vecs) != nmod:
        raise DimensionMismatchError(
            f"Need {nmod} q-vectors (semicolon-separated); got {len(vecs)}"
        )
    return vecs


@dataclass
class SuperspaceResult:
    """一次 (3+nmod) 计算的完整结果（可 JSON/YAML 往返）。"""

    nmod: int
    space_group_number: int
    group: SuperspaceGroup
    lattice: SuperspaceLattice
    ks: KsVector
    irreps: list[Irrep] = field(default_factory=list)
    opds: list[OpdOperator] = field(default_factory=list)
    modes: list[SuperspaceMode] = field(default_factory=list)
    little_group_order: int = 0

    def to_dict(self) -> dict:
        kind = SERIAL_KIND_3P1 if self.nmod == 1 else SERIAL_KIND
        return {
            "kind": kind,
            "nmod": self.nmod,
            "d": self.nmod,
            "space_group_number": self.space_group_number,
            "space_group_symbol": self.group.space_group_symbol,
            "schoenflies": self.group.schoenflies,
            "setting": self.group.setting,
            "group": self.group.to_dict(),
            "lattice": self.lattice.to_dict(),
            "ks": self.ks.to_dict(),
            "irreps": [ir.to_dict() for ir in self.irreps],
            "opds": [opd.to_dict() for opd in self.opds],
            "modes": [m.to_dict() for m in self.modes],
            "little_group_order": self.little_group_order,
        }

    @classmethod
    def from_dict(cls, data: dict) -> SuperspaceResult:
        nmod = validate_nmod(data.get("nmod", data.get("d", 0)))
        group = SuperspaceGroup.from_dict(data.get("group") or data)
        lattice = SuperspaceLattice.from_dict(
            data.get("lattice")
            or {
                "nmod": nmod,
                "lattice_3d": default_lattice_matrix(int(data.get("space_group_number", 1))),
                "space_group_number": int(data.get("space_group_number", 1)),
            }
        )
        ks = KsVector.from_dict(data["ks"]) if "ks" in data else KsVector(
            coords=np.zeros(dim_3pd(nmod)), nmod=nmod
        )
        return cls(
            nmod=nmod,
            space_group_number=int(data.get("space_group_number", group.space_group_number)),
            group=group,
            lattice=lattice,
            ks=ks,
            irreps=[Irrep.from_dict(x) for x in data.get("irreps", [])],
            opds=[OpdOperator.from_dict(x) for x in data.get("opds", [])],
            modes=[SuperspaceMode.from_dict(x) for x in data.get("modes", [])],
            little_group_order=int(data.get("little_group_order", 0)),
        )

    def equivalent(self, other: SuperspaceResult, tol: float | None = None) -> bool:
        tol = eps() if tol is None else tol
        if self.nmod != other.nmod:
            return False
        if self.space_group_number != other.space_group_number:
            return False
        if self.group.order != other.group.order:
            return False
        if not np.allclose(self.ks.coords, other.ks.coords, atol=tol):
            return False
        if len(self.irreps) != len(other.irreps):
            return False
        if len(self.modes) != len(other.modes):
            return False
        return True


def run_superspace_workflow(
    space_group_number: int,
    nmod: int,
    *,
    q_vectors: Sequence[Sequence[float]] | None = None,
    ks_coords: Sequence[float] | None = None,
    k_point_label: str = "",
    lattice_3d: Sequence[Sequence[float]] | None = None,
    check_ks: bool = True,
) -> SuperspaceResult:
    """空间群 → (3+nmod) 超空间群 → k_s → IR → OPD → 模式 → 三维投影。"""
    nmod = validate_nmod(nmod)
    sg = validate_space_group_number(space_group_number)
    group = SuperspaceGroup.from_space_group(sg, nmod=nmod, q_vectors=q_vectors)
    lat_mat = (
        np.asarray(lattice_3d, dtype=float)
        if lattice_3d is not None
        else default_lattice_matrix(sg)
    )
    lattice = SuperspaceLattice(lattice_3d=lat_mat, nmod=nmod, space_group_number=sg)
    if ks_coords is None:
        if nmod and group.q_vectors:
            ks = KsVector.from_3d(
                group.q_vectors[0],
                nmod,
                k_point_label=k_point_label,
                space_group_number=sg,
            )
        else:
            ks = KsVector.from_3d(
                [0.0, 0.0, 0.0],
                nmod,
                satellite=[0.0] * nmod if nmod else None,
                k_point_label=k_point_label or "GM",
                space_group_number=sg,
            )
    else:
        ks = KsVector(
            coords=ks_coords,
            nmod=nmod,
            k_point_label=k_point_label,
            space_group_number=sg,
        )
    ks = ks.reduce(group.centering)
    if check_ks and nmod and group.q_vectors:
        assert_ks_compatible(ks, group.q_vectors)
    irreps = build_irreps(group, ks, k_point_label=k_point_label or ks.k_point_label)
    opds = generate_opds(irreps, nmod=nmod)
    allowed = filter_allowed_opds(opds, irreps)
    modes = generate_modes(group, ks, irreps, allowed)
    lg = little_group(group, ks)
    return SuperspaceResult(
        nmod=nmod,
        space_group_number=sg,
        group=group,
        lattice=lattice,
        ks=ks,
        irreps=irreps,
        opds=allowed,
        modes=modes,
        little_group_order=len(lg),
    )


def save_superspace_result(result: SuperspaceResult, path: str | Path) -> Path:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    data = result.to_dict()
    suffix = dest.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        dest.write_text(
            yaml.safe_dump(data, allow_unicode=True, sort_keys=False),
            encoding="utf-8",
        )
        return dest
    if suffix != ".json":
        dest = dest.with_suffix(".json")
    ser = ResultSerializer(output_dir=dest.parent)
    written = ser.save(data, dest.stem)
    if written.resolve() != dest.resolve():
        dest.write_text(written.read_text(encoding="utf-8"), encoding="utf-8")
        return dest
    return written


def load_superspace_result(path: str | Path) -> SuperspaceResult:
    src = Path(path)
    if not src.is_file():
        # ResultSerializer 约定不带后缀
        json_path = src if src.suffix else src.with_suffix(".json")
        if json_path.is_file():
            src = json_path
        else:
            raise InputError(f"superspace result file not found: {path}")
    text = src.read_text(encoding="utf-8")
    if src.suffix.lower() in {".yaml", ".yml"}:
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    if not isinstance(data, dict):
        raise InputError("superspace file must contain a JSON/YAML object")
    return SuperspaceResult.from_dict(data)


def format_superspace_report(result: SuperspaceResult) -> str:
    lines = [
        f"Superspace (3+{result.nmod})  nmod={result.nmod}  setting={result.group.setting}",
        f"Parent space group: {result.space_group_number} "
        f"{result.group.space_group_symbol} {result.group.schoenflies}",
        f"Operations: {result.group.order}  dim={result.group.dim}",
        f"k_s: {result.ks.coords.tolist()}  label={result.ks.k_point_label or '-'}",
        f"Little-group order: {result.little_group_order}",
        f"IRs ({len(result.irreps)}):",
    ]
    for ir in result.irreps:
        lines.append(f"  {ir.label:<10s} dim={ir.dimension}")
    lines.append(f"OPDs ({len(result.opds)}):")
    for opd in result.opds:
        lines.append(f"  {opd.irrep_label:<10s} {opd.opd_symbol:<4s} {opd.opd_vector}")
    lines.append(f"Modes ({len(result.modes)}):")
    for mode in result.modes:
        lines.append(
            f"  {mode.irrep_label:<10s} OPD={mode.opd_symbol:<4s} "
            f"3D components={mode.project_to_3d()}"
        )
    return "\n".join(lines)


def run_superspace_cli(
    *,
    nmod: object,
    space_group: int = 139,
    ks_text: str | None = None,
    q_text: str | None = None,
    k_point_label: str = "",
    export: str | None = None,
    load_path: str | None = None,
    print_report: bool = True,
) -> tuple[int, str, SuperspaceResult | None]:
    """终端 ``--superspace-d`` 入口。返回 (exit_code, stdout, result)。"""
    if load_path:
        result = load_superspace_result(load_path)
        text = format_superspace_report(result)
        if export:
            save_superspace_result(result, export)
            text += f"\nWrote {export}"
        if print_report:
            print(text)
        return 0, text, result

    try:
        d = validate_nmod(nmod)
    except Exception as exc:  # noqa: BLE001 - CLI 边界：非法 d 转非 0 退出码
        msg = str(exc)
        print(msg)
        return 2, msg, None

    if d > 1:
        note = (
            f"Note: nmod={d} is a (3+{d}) task; (3+1)-only dump format is not used. "
            f"Writing generic superspace JSON/YAML (nmod field)."
        )
    else:
        note = ""

    q_vectors = parse_q_vectors_text(q_text, d) if d else None
    ks_coords = parse_ks_text(ks_text, d) if ks_text else None
    try:
        result = run_superspace_workflow(
            space_group,
            d,
            q_vectors=q_vectors,
            ks_coords=ks_coords,
            k_point_label=k_point_label,
            check_ks=bool(q_vectors),
        )
    except Exception as exc:  # noqa: BLE001 - CLI 边界
        msg = str(exc)
        print(msg)
        return 1, msg, None
    text = format_superspace_report(result)
    if note:
        text = note + "\n" + text
    if export:
        path = save_superspace_result(result, export)
        text += f"\nWrote {path}"
    if print_report:
        print(text)
    return 0, text, result
