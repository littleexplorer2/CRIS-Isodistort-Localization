"""Application workflows shared by transport adapters.

Owns the live session, atomic revisions, upload lifetime and export snapshots.
The HTTP adapter supplies decoded inputs and serializes the returned values.
"""
from __future__ import annotations

import contextlib
import copy
import re
import threading
import uuid
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from application.session import IsoDistort
from backend.export import parse_export_formats, parse_export_method
from backend.messages import MESSAGES
from backend.method1 import DISTORTION_TYPES
from backend.utils import get_config
from backend.utils.schoenflies import hm_symbol, schoenflies_symbol

_COUPLED_METHOD2_UNSUPPORTED = (
    "Method 2 with more than one superposed IR component is not implemented "
    "locally. Choose one component or use the official ISODISTORT website; "
    "independent single-IR subgroup tables cannot be combined."
)


def _require_single_method2_component(data: dict) -> None:
    """Reject coupled Method 2 requests before any candidate state is changed.

    ``nsup`` counts submitted IR components, not distinct k-point labels.  In
    particular, two components at the same k point are still a coupled request
    and must not be collapsed or implemented as a union of two single-IR
    inventories.
    """

    counts: dict[str, int] = {}
    if "nsup" in data:
        raw_nsup = data["nsup"]
        if type(raw_nsup) is not int:
            raise ValueError("nsup must be a JSON integer")
        nsup = raw_nsup
        if nsup < 1:
            raise ValueError("nsup must be a positive integer")
        counts["nsup"] = nsup

    for key in ("components", "kpoints"):
        if key not in data or data[key] is None:
            continue
        values = data[key]
        if not isinstance(values, list):
            raise ValueError(f"{key} must be a list of Method 2 components")
        if not values:
            raise ValueError(f"{key} must contain at least one Method 2 component")
        counts[key] = len(values)

    if any(count > 1 for count in counts.values()):
        raise ValueError(_COUPLED_METHOD2_UNSUPPORTED)
    if len(set(counts.values())) > 1:
        rendered = ", ".join(f"{key}={value}" for key, value in counts.items())
        raise ValueError(f"inconsistent Method 2 component counts: {rendered}")


class WebSession:
    """网页会话：持有唯一 IsoDistort 实例（单用户本地工具）。

    IsoDistort 构造较慢（底层 BaseWrapper 在 Windows 下会初始化 WSL 短路径
    暂存目录与 ISODATA 符号链接，约 3-4 秒），采用懒初始化：首次访问
    ``iso`` 属性时才创建，避免拖慢模块导入（否则 main_web.py 启动会阻塞数秒）。
    """

    def __init__(self) -> None:
        self._iso = None
        self.parent_upload_path: Path | None = None
        # Monotonic context version used to reject writes from an older tab.
        # It protects the shared single-user session without pretending that
        # each browser tab owns an independent IsoDistort instance.
        self.revision = 0
        # 官网默认（见 webpage_info 第 2 页 HTML）：includestrain 勾选，
        # Displacive 行的各物种复选框逐个勾选（Eu/Al，等价于全部物种），
        # Occupational/Magnetic/Rotational 整行不勾选
        self.distortion_types: list[str] = ["strain", "displacive"]
        # 各畸变类型的作用域物种（官网 all/none/Eu/Al 复选框），"*"=全部
        self.distortion_scope: dict[str, list[str]] = {
            "displacive": ["*"],
            "occupational": [],
            "strain": [],
            "magnetic": [],
            "rotational": [],
        }
        self.method1: list = []
        self.method2 = None
        self.method3: list = []
        # Method 2 k 点枚举得到的子群（Download all 只导出这份列表，不扫 output_dir）
        self.method2_subgroups: list = []

    @property
    def iso(self):
        if self._iso is None:
            self._iso = IsoDistort()
        return self._iso

    def replace_parent_upload(self, path: str | Path) -> None:
        previous = self.parent_upload_path
        self.parent_upload_path = Path(path)
        if previous is not None and previous != self.parent_upload_path:
            _discard_upload(previous)

    def cleanup_uploads(self) -> None:
        if self.parent_upload_path is not None:
            _discard_upload(self.parent_upload_path)
            self.parent_upload_path = None


_SESSION = WebSession()
# IsoDistort owns mutable structure, candidate and mode state.  Serialize live
# session reads/mutations and export snapshot creation so concurrent requests
# cannot observe a half-updated pool.  Long export work runs on that independent
# snapshot after this lock is released.  A monotonic revision rejects stale
# writes from another tab; tabs still share one live session.
_SESSION_LOCK = threading.RLock()
_EXPORT_JOB_LOCK = threading.Lock()
_EXPORT_PROGRESS_LOCK = threading.Lock()
_EXPORT_CANCEL = threading.Event()


def _initial_export_progress() -> dict:
    return {
        "active": False,
        "phase": "idle",
        "total": 0,
        "completed": 0,
        "successful": 0,
        "ineligible": 0,
        "failed": 0,
        "estimated_remaining_seconds": None,
        "cancellation_requested": False,
    }


_EXPORT_PROGRESS: dict = _initial_export_progress()


def _export_progress_snapshot() -> dict:
    with _EXPORT_PROGRESS_LOCK:
        return dict(_EXPORT_PROGRESS)


def _update_export_progress(**values) -> None:
    with _EXPORT_PROGRESS_LOCK:
        _EXPORT_PROGRESS.update(values)


def _reset_export_progress(**values) -> None:
    """Start a clean export status record without stale job-specific fields."""
    with _EXPORT_PROGRESS_LOCK:
        _EXPORT_PROGRESS.clear()
        _EXPORT_PROGRESS.update(_initial_export_progress())
        _EXPORT_PROGRESS.update(values)


class StaleContextError(RuntimeError):
    """A browser mutation was based on an older shared-session snapshot."""


def _assert_current_revision(value) -> None:
    """Validate an optional client revision while ``_SESSION_LOCK`` is held."""
    if value is None or value == "":
        return  # Backward compatibility for API clients predating revisions.
    try:
        expected = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("revision must be an integer") from exc
    if expected != _SESSION.revision:
        raise StaleContextError(
            "Stale context: this tab used revision "
            f"{expected}, but the current session is revision {_SESSION.revision}. "
            "Its cached results were discarded; retry from the current state."
        )

# 230 个空间群（序号 + HM 符号 + Schoenflies 符号），供 Method 1/3 下拉使用
_SPACE_GROUPS = [
    {
        "number": i,
        "symbol": hm_symbol(i),
        "schoenflies": schoenflies_symbol(i),
    }
    for i in range(1, 231)
]


def _space_groups() -> list[dict]:
    return _SPACE_GROUPS


def _state_summary() -> dict:
    """返回当前会话状态摘要（供前端展示）。"""
    iso = _SESSION.iso
    modes = list(iso.mode_displacements.keys()) + list(iso.mode_occupancies.keys())
    summary = {
        "revision": _SESSION.revision,
        "structure": None,
        "subgroups": len(iso.subgroups),
        "modes": modes,
        "distorted_atoms": len(iso.distorted_structure) if iso.distorted_structure else None,
        "distortion_types": _SESSION.distortion_types,
        "distortion_scope": _SESSION.distortion_scope,
        "species": iso.species(),
    }
    if iso.structure is not None:
        lattice = iso.structure.lattice
        summary["structure"] = {
            "space_group_number": iso.symmetry_info["space_group_number"],
            "space_group_symbol": iso.symmetry_info["space_group_symbol"],
            "space_group_schoenflies": schoenflies_symbol(
                iso.symmetry_info["space_group_number"]
            ),
            "atoms": len(iso.structure),
            "preferences": iso.space_group_preferences(),
            "lattice": {
                "a": round(float(lattice.a), 5),
                "b": round(float(lattice.b), 5),
                "c": round(float(lattice.c), 5),
                "alpha": round(float(lattice.alpha), 5),
                "beta": round(float(lattice.beta), 5),
                "gamma": round(float(lattice.gamma), 5),
            },
            "wyckoff": [
                {
                    "letter": s["wyckoff_letter"],
                    "multiplicity": s["multiplicity"],
                    "species": s["species"],
                    "coordinates": [
                        round(float(x), 6)
                        for x in iso.structure[s["representative_index"]].frac_coords
                    ],
                }
                for s in iso.symmetry_info["wyckoff_sites"]
            ],
            "wyckoff_display": iso.parent_wyckoff_display(),
        }
    return summary


def _subgroup_rows(subgroups) -> list[dict]:
    """把 SubgroupInfo 列表序列化为前端友好的 dict 列表。"""
    rows = []
    for sg in subgroups:
        identity_status = str(
            getattr(sg, "_mode_identity_status", "not_checked") or "not_checked"
        )
        identity_failures = tuple(
            getattr(sg, "_mode_identity_failures", ()) or ()
        )
        rows.append({
            "index": sg.index,
            "space_group_number": sg.space_group_number,
            "space_group_symbol": sg.space_group_symbol,
            "subgroup_index": sg.subgroup_index,
            "size": sg.size,
            "is_maximal": sg.is_maximal,
            "opd_symbol": sg.opd_symbol,
            "k_point_label": sg.k_point_label,
            "irrep_label": sg.irrep_label,
            "basis_vectors": sg.basis_vectors,
            "origin": sg.origin,
            "k_parameters": list(sg.k_parameters or []),
            "mode_identity_status": identity_status,
            "mode_identity_warning": "; ".join(
                str(getattr(failure, "reason", "") or "")
                for failure in identity_failures
                if str(getattr(failure, "reason", "") or "")
            ),
            "selectable": identity_status != "unresolved",
        })
    return rows


def _method1_rows(items) -> list[dict]:
    rows = []
    for item in items:
        sg = item.subgroup
        fields = sg.official_fields()
        rows.append({
            "index": sg.index,
            "space_group_number": sg.space_group_number,
            "space_group_symbol": sg.space_group_symbol,
            "crystal_system": item.crystal_system,
            "is_maximal": item.is_maximal,
            "k_point_label": sg.k_point_label,
            "irrep_label": fields["irrep"],
            "opd_symbol": fields["opd"],
            "opd_dir_raw": fields["dir"],
            "size": fields["s"],
            "subgroup_index": fields["i"],
            "basis": fields["basis"],
            "origin": fields["origin"],
            "k_active": fields["k_active"],
            "opd_line": sg.opd_line(),
            "k_parameters": list(sg.k_parameters or []),
        })
    return rows


def _method3_rows(items) -> list[dict]:
    """把 Method 3 的 Method3ResultItem 序列化为前端友好的 dict 列表。

    官网首屏的科学身份是 SG/basis/origin/s/i embedding，而不是某一个
    来源 IR/OPD。兼容代表 route 仍随行返回供当前模式计算使用，但界面把
    所有已知 route 单独标成诊断信息，不能拿它当候选身份。
    """
    rows = []
    for item in items:
        sg = item.subgroup
        fields = sg.official_fields()
        routes = item.routes if item.routes is not None else [sg]
        known_routes = [
            {
                "k_point_label": route.k_point_label,
                "irrep_label": route.irrep_label,
                "opd_symbol": route.opd_symbol,
                "k_parameters": list(route.k_parameters or []),
            }
            for route in routes
        ]
        rows.append({
            "index": sg.index,
            "space_group_number": sg.space_group_number,
            "space_group_symbol": sg.space_group_symbol,
            "subgroup_index": sg.subgroup_index,
            "size": sg.size,
            "is_maximal": sg.is_maximal,
            "opd_symbol": sg.opd_symbol,
            "k_point_label": sg.k_point_label,
            "irrep_label": sg.irrep_label,
            "point_group": item.point_group,
            "basis_vectors": sg.basis_vectors,
            "origin": sg.origin,
            "basis_display": fields["basis"],
            "origin_display": fields["origin"],
            "k_parameters": list(sg.k_parameters or []),
            "known_routes": known_routes,
            "known_route_count": len(known_routes),
            "route_resolution": getattr(item, "route_resolution", "known_single_ir"),
            "selectable": bool(known_routes) or getattr(
                item, "route_resolution", ""
            ) == "exact_fixed_space",
        })
    return rows


def _write_upload(filename: str, content: str) -> str:
    """把上传的 CIF 内容写入临时目录，返回文件路径。"""
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", Path(filename).name) or "upload.cif"
    cfg = get_config()
    upload_dir = cfg.temp_dir / "web_uploads"
    upload_dir.mkdir(parents=True, exist_ok=True)
    path = upload_dir / f"{uuid.uuid4().hex}-{safe}"
    path.write_text(content, encoding="utf-8")
    return str(path)


def _discard_upload(path: str | Path) -> None:
    """Remove one generated upload and its directory when it becomes empty."""
    upload_path = Path(path)
    with contextlib.suppress(FileNotFoundError):
        upload_path.unlink()
    with contextlib.suppress(FileNotFoundError, OSError):
        upload_path.parent.rmdir()


def _api_load_cif(data: dict) -> dict:
    content = data.get("content", "")
    filename = data.get("filename", "upload.cif")
    if not content.strip():
        raise ValueError("CIF 内容为空 / CIF content is empty")
    path = _write_upload(filename, content)
    try:
        _SESSION.iso.load_structure(path)
    except Exception:
        _discard_upload(path)
        raise
    _SESSION.replace_parent_upload(path)
    _SESSION.iso.set_distortion_scope(_SESSION.distortion_scope)
    _SESSION.iso.set_distortion_types(_SESSION.distortion_types)
    _SESSION.method1, _SESSION.method2, _SESSION.method3 = [], None, []
    _SESSION.method2_subgroups = []
    return {}


def _api_set_types(data: dict) -> dict:
    types = data.get("types", ["strain"])
    valid = set(DISTORTION_TYPES)
    _SESSION.distortion_types = [t for t in types if t in valid] or ["strain"]
    scope = data.get("scope")
    if scope is not None:
        _SESSION.distortion_scope = {
            tp: (["*"] if (v == "*" or v == "all") else
                 (v if isinstance(v, list) else []))
            for tp, v in scope.items() if tp in valid
        }
        # 同步到底层 IsoDistort（模式计算按作用域过滤）
        _SESSION.iso.set_distortion_scope(_SESSION.distortion_scope)
    _SESSION.iso.set_distortion_types(_SESSION.distortion_types)
    # Every result below the Types panel depends on these selections.
    # Clear both the web-owned tables and all derived engine state so a
    # later click cannot reuse a path computed under the previous scope.
    _SESSION.method1, _SESSION.method2, _SESSION.method3 = [], None, []
    _SESSION.method2_subgroups = []
    _SESSION.iso.set_subgroup_candidates([])
    _SESSION.iso.clear_selected_modes()
    return {}


def _api_method1(data: dict) -> dict:
    lattice = data.get("lattice")
    if lattice is not None:
        lattice = _SESSION.iso.lattice_in_conventional_frame(
            lattice, data.get("frame", "conventional")
        )
    result = _SESSION.iso.search_method_1(
        distortion_types=data.get("distortion_types", _SESSION.distortion_types),
        crystal_system=data.get("crystal_system") or None,
        subgroup_space_group=data.get("subgroup_space_group") or None,
        lattice=lattice,
        lattice_kind=data.get("lattice_kind", "conventional"),
        maximal_subgroup_only=bool(data.get("maximal_subgroup_only", False)),
    )
    _SESSION.method1 = result
    return {"candidates": _method1_rows(result)}


def _api_subgroups(data: dict) -> dict:
    # 对齐官网 Method 2 的单组分量：枚举指定 k 点
    # （+ 参数）下全部 IR 的子群。多组分量需要 coupled-IR
    # 稳定子的交集枚举，不能拼接各自的单-IR 子群表。
    # 候选池由应用层会话统一持有，后续模式计算必须携带同一查询上下文。
    _require_single_method2_component(data)
    gen = bool(data.get("generate", False))
    groups = data.get("kpoints")
    if groups:
        group = groups[0]
        subs = _SESSION.iso.list_subgroups_at_kpoint(
            group["k"],
            k_parameters=group.get("params"),
            generate_if_missing=gen,
        )
        _SESSION.iso.set_subgroup_candidates(subs)
        _SESSION.method2_subgroups = list(subs)
        return {"subgroups": _subgroup_rows(subs)}
    # 兼容旧版单 k 点（可带 ir 参数）路径
    if data.get("ir"):
        subs = _SESSION.iso.list_subgroups_at(
            data["k"], data["ir"],
            k_parameters=data.get("params"),
            opd_symbol=data.get("opd"),
            generate_if_missing=gen,
        )
    else:
        subs = _SESSION.iso.list_subgroups_at_kpoint(
            data["k"],
            k_parameters=data.get("params"),
            generate_if_missing=gen,
        )
    _SESSION.iso.set_subgroup_candidates(subs)
    _SESSION.method2_subgroups = list(subs)
    return {"subgroups": _subgroup_rows(subs)}


def _api_method2(data: dict) -> dict:
    _require_single_method2_component(data)
    idx = data.get("subgroup_idx")
    if idx is None:
        raise ValueError("subgroup_idx 缺失 / subgroup_idx missing")
    idx = int(idx)
    iso = _SESSION.iso
    # 按结果表来源显式传入候选池。不同 Method 的 index 会重复，
    # 因此不能再通过覆盖 iso.subgroups 来隐式切换上下文。
    source = data.get("source")
    if source == "method1":
        if not _SESSION.method1:
            raise ValueError("Method 1 candidate table is not available")
        candidates = [item.subgroup for item in _SESSION.method1]
    elif source == "method3":
        if not _SESSION.method3:
            raise ValueError("Method 3 candidate table is not available")
        selected_item = next(
            (
                item
                for item in _SESSION.method3
                if int(item.subgroup.index) == idx
            ),
            None,
        )
        if (
            selected_item is not None
            and selected_item.routes == []
            and getattr(selected_item, "route_resolution", "")
            != "exact_fixed_space"
        ):
            raise ValueError(
                "This affine embedding has no resolved single-IR or coupled-IR "
                "second-stage route; mode calculation is not implemented for it"
            )
        candidates = [item.subgroup for item in _SESSION.method3]
    elif source == "subgroups":
        if not _SESSION.method2_subgroups:
            raise ValueError("Method 2 subgroup table is not available")
        candidates = list(_SESSION.method2_subgroups)
    elif source in (None, ""):
        # Backward compatibility for older clients that did not send source.
        candidates = list(_SESSION.method2_subgroups) or None
    else:
        raise ValueError(f"Unknown Method 2 candidate source: {source}")
    iso.set_distortion_scope(_SESSION.distortion_scope)
    iso.set_distortion_types(_SESSION.distortion_types)
    nmod = data.get("nmod", data.get("number_of_independent_modulations", 0))
    try:
        nmod = int(nmod or 0)
    except (TypeError, ValueError):
        nmod = 0
    result = iso.search_method_2(
        subgroup_idx=idx,
        distortion_type=data.get("distortion_type", _SESSION.distortion_types),
        number_of_independent_modulations=nmod,
        number_of_superposed_irreps=int(data.get("nsup", 1) or 1),
        candidates=candidates,
    )
    _SESSION.method2 = result
    modes = []
    for m in result.modes:
        key = str(getattr(m, "amplitude_key", "") or m.irrep_label)
        pretty = iso._mode_label_overrides.get(key, "")
        identity = getattr(m, "mode_identity", None)
        modes.append({
            "irrep_label": key,
            "pretty_label": pretty or m.irrep_label,
            "opd_symbol": m.opd_symbol,
            "mode_type": m.mode_type,
            "wyckoff_sites": sorted({b.wyckoff_letter for b in m.bush_modes}),
            "n_representatives": len(m.bush_modes),
            "site_irrep": getattr(m, "site_irrep", "") or "",
            "k_coords": getattr(m, "k_coords_label", "") or "",
            "identity_status": getattr(identity, "status", "unresolved"),
            "identity_reason": getattr(identity, "reason", "") or "",
        })
    for label, entry in iso.mode_occupancies.items():
        om = entry["mode"]
        modes.append({
            "irrep_label": label,
            "opd_symbol": om.irrep_label or "",
            "mode_type": "occupational",
            "wyckoff_sites": [om.wyckoff_letter],
            "n_representatives": int(np.count_nonzero(om.pattern)),
            "validated": entry["validated"],
            "note": entry["note"],
        })
    metadata = getattr(result, "metadata", None) or {}
    identity_failures = list(metadata.get("mode_identity_failures") or [])
    return {
        "modes": modes,
        "export_ready": bool(metadata.get("export_ready", True)),
        "mode_identity_status": metadata.get(
            "mode_identity_status", "not_applicable"
        ),
        "mode_identity_failures": identity_failures,
        "export_warning": (
            "This candidate has a numerical mode basis, but its ISO microscopic "
            "identity is unresolved and authoritative subgroup files cannot be "
            "exported. "
            + "; ".join(
                str(failure.get("reason") or "identity unresolved")
                for failure in identity_failures
            )
            if identity_failures
            else ""
        ),
    }


def _api_method3(data: dict) -> dict:
    result = _SESSION.iso.search_method_3(
        distortion_types=data.get("distortion_types", _SESSION.distortion_types),
        point_group=data.get("point_group") or None,
        space_group_type=data.get("space_group_type") or None,
        supercell_basis=data.get("supercell_basis"),
        direct_sublattice_centering=data.get("direct_sublattice_centering") or None,
        lattice_type=data.get("lattice_type", "direct"),
        generate_if_missing=bool(data.get("generate", False)),
    )
    _SESSION.method3 = result
    return {"candidates": _method3_rows(result)}


def _api_method4(data: dict) -> dict:
    content = data.get("content", "")
    if not content.strip():
        raise ValueError("Daughter CIF 内容为空 / daughter CIF content is empty")
    path = _write_upload(data.get("filename", "daughter.cif"), content)
    try:
        result = _SESSION.iso.search_method_4(
            distorted_cif_path=path,
            atom_matching_method=data.get("atom_matching_method", "nearest-site"),
            robust_distance_threshold=float(data.get("robust_distance_threshold", 0.25)),
            provided_origin_shift=data.get("provided_origin_shift"),
        )
    finally:
        _discard_upload(path)
    ranked = sorted(result.amplitudes.items(), key=lambda kv: abs(kv[1]), reverse=True)
    return {
        "amplitudes": {k: float(v) for k, v in ranked},
        "parent_cell_amplitudes": {
            k: float(result.parent_cell_amplitudes[k]) for k, _ in ranked
        },
        "raw_coefficients": {
            k: float(result.raw_coefficients[k]) for k, _ in ranked
        },
        "mode_normfactors": {
            k: float(result.mode_normfactors[k]) for k, _ in ranked
        },
        "strain_mode_amplitudes": dict(result.strain_mode_amplitudes),
        "strain_modes": list(result.strain_modes),
        "strain_raw_coordinate_sum_parent_basis": (
            dict(result.strain_raw_coordinate_sum_parent_basis)
            if result.strain_raw_coordinate_sum_parent_basis is not None
            else None
        ),
        "strain_applied_engineering_q_parent_basis": dict(
            result.strain_applied_engineering_q_parent_basis
        ),
        "strain_tensor_parent_basis": result.strain_tensor_parent_basis,
        "strain_multiplier_parent_basis": result.strain_multiplier_parent_basis,
        # Deprecated JSON aliases; both retain the same parent-basis
        # semantics and never carry the discarded Cartesian contract.
        "strain_voigt_engineering": dict(
            result.strain_applied_engineering_q_parent_basis
        ),
        "strain_tensor": result.strain_tensor_parent_basis,
        "rms_residual": result.rms_residual,
        "max_abs_residual": result.max_abs_residual,
        "metadata": result.metadata,
    }


def _api_isotropy_cache_list() -> dict:
    entries = _SESSION.iso.list_isotropy_cache()
    return {
        "entries": [e.to_dict() for e in entries],
        "count": len(entries),
    }


def _api_isotropy_cache_delete(data: dict) -> dict:
    names = data.get("names") or []
    if not isinstance(names, list):
        raise ValueError("names must be a list of cache file names")
    result = _SESSION.iso.delete_isotropy_cache([str(n) for n in names])
    remaining = _api_isotropy_cache_list()
    return {**result, **remaining}


def _subgroups_for_method(method: int) -> list:
    """取出所选 Method 当次计算得到的子群列表（不含其它 Method）。"""
    if method == 1:
        return [item.subgroup for item in _SESSION.method1]
    if method == 3:
        return [item.subgroup for item in _SESSION.method3]
    subs = list(_SESSION.method2_subgroups)
    if not subs and _SESSION.method2 is not None:
        subs = [_SESSION.method2.subgroup]
    return subs


@dataclass(frozen=True)
class RuntimeSettings:
    web_port: int
    web_idle_timeout: float


@dataclass(frozen=True)
class Download:
    filename: str
    body: bytes


class ExportBusyError(RuntimeError):
    """A previous export still owns the single export worker."""


class ExportFailedError(RuntimeError):
    """The selected export worker failed after its input was accepted."""


def runtime_settings() -> RuntimeSettings:
    cfg = get_config()
    return RuntimeSettings(cfg.web_port, cfg.web_idle_timeout)


def state_summary() -> dict:
    with _SESSION_LOCK:
        return _state_summary()


def messages() -> dict:
    return MESSAGES


def space_groups() -> list[dict]:
    return _space_groups()


def export_progress_snapshot() -> dict:
    return _export_progress_snapshot()


def cancel_export() -> bool:
    progress = _export_progress_snapshot()
    if progress.get("active"):
        _EXPORT_CANCEL.set()
        _update_export_progress(
            phase="cancellation_requested", cancellation_requested=True,
        )
    return bool(progress.get("active"))


def reset_export_progress() -> None:
    _EXPORT_CANCEL.clear()
    _reset_export_progress()


def cleanup_uploads() -> None:
    with _SESSION_LOCK:
        _SESSION.cleanup_uploads()


def _api_kpoints(_data: dict) -> dict:
    return {"kpoints": [
        {"label": kp.label, "coordinates": kp.coordinates,
         "parameters": kp.parameters, "is_special": kp.is_special,
         "kovalev": kp.kovalev}
        for kp in _SESSION.iso.list_k_points()
    ]}


def _api_irreps(data: dict) -> dict:
    return {"irreps": [
        {"label": ir.label, "dimension": ir.dimension, "active": ir.active}
        for ir in _SESSION.iso.list_irreps(data.get("k", ""), data.get("params"))
    ]}


def _api_method1_options(_data: dict) -> dict:
    return {"options": _SESSION.iso.method1_options()}


_ACTIONS = {
    "load_cif": _api_load_cif,
    "set_types": _api_set_types,
    "method1": _api_method1,
    "subgroups": _api_subgroups,
    "method2": _api_method2,
    "method3": _api_method3,
    "method4": _api_method4,
    "isotropy_cache": lambda _data: _api_isotropy_cache_list(),
    "isotropy_cache/delete": _api_isotropy_cache_delete,
    "kpoints": _api_kpoints,
    "irreps": _api_irreps,
    "method1_options": _api_method1_options,
}
_MUTATIONS = frozenset({
    "load_cif", "set_types", "method1", "subgroups", "method2", "method3", "method4",
})


def execute(action: str, data: dict | None = None, *, revision=None) -> dict:
    """Apply an application action atomically and return its committed DTO.

    Loading a new parent intentionally replaces the context regardless of the
    sender's previous revision; all other mutations reject stale senders.
    """
    data = data or {}
    with _SESSION_LOCK:
        if action != "load_cif":
            _assert_current_revision(revision)
        result = _ACTIONS[action](data) or {}
        if action in _MUTATIONS:
            _SESSION.revision += 1
            result["state"] = _state_summary()
    result.setdefault("ok", True)
    return result


def read_download(filename: str) -> Download:
    """Read only a previously generated file inside the configured output root."""
    output_root = get_config().output_dir.resolve()
    path = (output_root / filename).resolve()
    if not path.is_relative_to(output_root):
        raise PermissionError("invalid file")
    with _SESSION_LOCK:
        if not path.is_file():
            raise FileNotFoundError("file not found")
        return Download(path.name, path.read_bytes())


def export_archive(values: dict[str, list[str]]) -> Download:
    """Export a frozen candidate pool, releasing the live-session lock first.

    Values are decoded form fields. This workflow has no HTTP or rendering
    dependency; a caller decides how to present validation and worker errors.
    """
    requested_revision = (values.get("revision") or [None])[0]
    with _SESSION_LOCK:
        _assert_current_revision(requested_revision)
    method_values = values.get("method") or ["2"]
    if len(method_values) > 1:
        raise ValueError("只能选择一个 Method 导出，不能多选 / select exactly one Method")
    method = parse_export_method(method_values[0])
    formats = parse_export_formats(
        (values.get("formats") or ["cif,isoviz,modes,topas"])[0]
    )
    if method == 2:
        nsup_values = values.get("nsup") or ["1"]
        if len(nsup_values) != 1:
            raise ValueError("nsup must be provided at most once")
        nsup_text = nsup_values[0]
        if re.fullmatch(r"[1-9]\d*", nsup_text) is None:
            raise ValueError("nsup query parameter must be a positive integer")
        _require_single_method2_component({"nsup": int(nsup_text)})
    indices_raw = (values.get("indices") or [""])[0].strip()
    wanted = None
    if indices_raw:
        try:
            wanted = {int(item.strip()) for item in indices_raw.split(",") if item.strip()}
        except ValueError as exc:
            raise ValueError("indices must be comma-separated integers") from exc
    with _SESSION_LOCK:
        # Recheck at the snapshot point, after parsing the potentially long input.
        _assert_current_revision(requested_revision)
        live_iso = _SESSION.iso
        if live_iso.structure is None:
            raise FileNotFoundError("请先加载母相 CIF / load a parent CIF first")
        subgroups = _subgroups_for_method(method)
        if wanted is not None:
            subgroups = [sg for sg in subgroups if sg.index in wanted]
        if not subgroups:
            raise FileNotFoundError(
                f"没有可导出的 Method {method} 子群；请先完成该 Method 的计算"
                f" / no Method {method} subgroups; run that Method first"
            )
        snapshot = getattr(live_iso, "snapshot_for_export", None)
        iso = snapshot() if callable(snapshot) else live_iso
        subgroups = copy.deepcopy(subgroups)
    need_modes = any(fmt != "cif" for fmt in formats)
    compute = (values.get("compute_modes") or ["1"])[0].strip().lower()
    compute_missing_modes = need_modes and compute not in ("0", "false", "no")
    nmod_raw = (values.get("nmod") or ["0"])[0].strip()
    try:
        export_nmod = int(nmod_raw or 0)
    except ValueError:
        export_nmod = 0
    if not _EXPORT_JOB_LOCK.acquire(blocking=False):
        raise ExportBusyError("Another subgroup export is already running")
    _EXPORT_CANCEL.clear()
    _reset_export_progress(
        active=True, phase="starting", method=method, total=len(subgroups),
        completed=0, successful=0, ineligible=0, failed=0,
        estimated_remaining_seconds=None, cancellation_requested=False,
    )

    def progress_callback(update: dict) -> None:
        _update_export_progress(active=True, **update)

    try:
        body = iso.export_subgroups_zip(
            formats=formats, subgroups=subgroups,
            compute_missing_modes=compute_missing_modes, wrapping=None,
            export_method=method,
            number_of_independent_modulations=export_nmod,
            progress_callback=progress_callback,
            cancel_check=_EXPORT_CANCEL.is_set,
        )
    except Exception as exc:
        _update_export_progress(active=False, phase="failed", error=str(exc))
        raise ExportFailedError(str(exc)) from exc
    finally:
        _EXPORT_JOB_LOCK.release()
    # Only the collector's candidate-boundary event proves cancellation.
    final_phase = (
        "cancelled" if _export_progress_snapshot().get("phase") == "cancelled"
        else "complete"
    )
    _update_export_progress(
        active=False, phase=final_phase, estimated_remaining_seconds=0.0,
    )
    return Download(f"isodistort_method{method}.zip", body)
