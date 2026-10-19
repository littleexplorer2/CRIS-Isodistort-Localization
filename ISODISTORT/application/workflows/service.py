"""Application workflows shared by transport adapters.

Owns the live session, atomic revisions, upload lifetime and export snapshots.
The HTTP adapter supplies decoded inputs and serializes the returned values.
"""
from __future__ import annotations
from backend.models.phase_path import validate_modulation_count

import contextlib
import copy
import re
import logging
import threading
import uuid
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from application.session import IsoDistort
from backend.export import parse_export_formats, parse_export_method
from backend.messages import MESSAGES
from backend.utils import get_config
from backend.utils.schoenflies import hm_symbol, schoenflies_symbol

from .inputs import action_input, export_job_id, nonnegative_integer, query_value
from .jobs import ExportJobs
from .wire import response_data

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

    def stage(self, *, replace_parent: bool = False) -> "WebSession":
        """Copy all web state together and isolate the mutable calculation engine."""
        staged = WebSession()
        if replace_parent:
            # A new parent starts a fresh runtime, including after a read-only
            # query retired the old wrapper on timeout. Retain only the user's
            # selections and revision; generated tables belong to the replaced
            # calculation. The old session stays live until the upload commits.
            staged.revision = self.revision
            staged.distortion_types, staged.distortion_scope = copy.deepcopy((
                self.distortion_types, self.distortion_scope,
            ))
            return staged
        fields = {key: value for key, value in vars(self).items() if key != "_iso"}
        staged.__dict__.update(copy.deepcopy(fields))
        if self._iso is not None:
            staged._iso = self._iso.snapshot_for_action()
        return staged

    def close(self) -> None:
        if self._iso is not None:
            _close_worker(self._iso)

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
_EXPORT_JOBS = ExportJobs()


def _close_worker(worker) -> None:
    """Housekeeping must neither hide an action failure nor undo its commit."""
    try:
        worker.close()
    except Exception:
        logging.getLogger(__name__).warning("Could not close calculation worker", exc_info=True)


class StaleContextError(RuntimeError):
    """A browser mutation was based on an older shared-session snapshot."""


def _assert_current_revision(value) -> None:
    """Validate an optional client revision while ``_SESSION_LOCK`` is held."""
    if value is None or value == "":
        return  # Backward compatibility for API clients predating revisions.
    expected = nonnegative_integer(value, "revision", text=isinstance(value, str))
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


def _state_summary(session: WebSession) -> dict:
    """返回当前会话状态摘要（供前端展示）。"""
    iso = session.iso
    modes = list(iso.mode_displacements.keys()) + list(iso.mode_occupancies.keys())
    summary = {
        "revision": session.revision,
        "structure": None,
        "subgroups": len(iso.subgroups),
        "modes": modes,
        "distorted_atoms": len(iso.distorted_structure) if iso.distorted_structure else None,
        "distortion_types": session.distortion_types,
        "distortion_scope": session.distortion_scope,
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
        identity_status = sg.mode_identity_status
        identity_failures = sg.mode_identity_failures
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
    try:
        path.write_text(content, encoding="utf-8")
    except BaseException:
        try:
            _discard_upload(path)
        except OSError:
            logging.getLogger(__name__).warning(
                "Could not remove interrupted upload", exc_info=True,
            )
        raise
    return str(path)


def _discard_upload(path: str | Path) -> None:
    """Remove one generated upload and its directory when it becomes empty."""
    upload_path = Path(path)
    with contextlib.suppress(FileNotFoundError):
        upload_path.unlink()
    with contextlib.suppress(FileNotFoundError, OSError):
        upload_path.parent.rmdir()


def _api_load_cif(session: WebSession, data: dict) -> dict:
    content = data.get("content", "")
    filename = data.get("filename", "upload.cif")
    if not content.strip():
        raise ValueError("CIF 内容为空 / CIF content is empty")
    path = _write_upload(filename, content)
    # Register ownership before constructing/loading the new engine so the
    # transaction cleans this upload for every failure, including cancellation.
    session.parent_upload_path = Path(path)
    session.iso.load_structure(path)
    session.iso.set_distortion_scope(session.distortion_scope)
    session.iso.set_distortion_types(session.distortion_types)
    session.method1, session.method2, session.method3 = [], None, []
    session.method2_subgroups = []
    return {}


def _api_set_types(session: WebSession, data: dict) -> dict:
    session.distortion_types = list(data["types"])
    scope = data.get("scope")
    if scope is not None:
        session.distortion_scope = copy.deepcopy(scope)
        # 同步到底层 IsoDistort（模式计算按作用域过滤）
        session.iso.set_distortion_scope(session.distortion_scope)
    session.iso.set_distortion_types(session.distortion_types)
    # Every result below the Types panel depends on these selections.
    # Clear both the web-owned tables and all derived engine state so a
    # later click cannot reuse a path computed under the previous scope.
    session.method1, session.method2, session.method3 = [], None, []
    session.method2_subgroups = []
    session.iso.set_subgroup_candidates([])
    session.iso.clear_selected_modes()
    return {}


def _api_method1(session: WebSession, data: dict) -> dict:
    lattice = data.get("lattice")
    if lattice is not None:
        lattice = session.iso.lattice_in_conventional_frame(
            lattice, data.get("frame", "conventional")
        )
    result = session.iso.search_method_1(
        distortion_types=data.get("distortion_types", session.distortion_types),
        crystal_system=data.get("crystal_system") or None,
        subgroup_space_group=data.get("subgroup_space_group") or None,
        lattice=lattice,
        lattice_kind=data.get("lattice_kind", "conventional"),
        maximal_subgroup_only=bool(data.get("maximal_subgroup_only", False)),
    )
    session.method1 = result
    return {"candidates": _method1_rows(result)}


def _api_subgroups(session: WebSession, data: dict) -> dict:
    # 对齐官网 Method 2 的单组分量：枚举指定 k 点
    # （+ 参数）下全部 IR 的子群。多组分量需要 coupled-IR
    # 稳定子的交集枚举，不能拼接各自的单-IR 子群表。
    # 候选池由应用层会话统一持有，后续模式计算必须携带同一查询上下文。
    gen = bool(data.get("generate", False))
    groups = data.get("kpoints")
    if groups:
        group = groups[0]
        subs = session.iso.list_subgroups_at_kpoint(
            group["k"],
            k_parameters=group.get("params"),
            generate_if_missing=gen,
        )
        session.iso.set_subgroup_candidates(subs)
        session.method2_subgroups = list(subs)
        return {"subgroups": _subgroup_rows(subs)}
    # 兼容旧版单 k 点（可带 ir 参数）路径
    if data.get("ir"):
        subs = session.iso.list_subgroups_at(
            data["k"], data["ir"],
            k_parameters=data.get("params"),
            opd_symbol=data.get("opd"),
            generate_if_missing=gen,
        )
    else:
        subs = session.iso.list_subgroups_at_kpoint(
            data["k"],
            k_parameters=data.get("params"),
            generate_if_missing=gen,
        )
    session.iso.set_subgroup_candidates(subs)
    session.method2_subgroups = list(subs)
    return {"subgroups": _subgroup_rows(subs)}


def _api_method2(session: WebSession, data: dict) -> dict:
    idx = data["subgroup_idx"]
    iso = session.iso
    # 按结果表来源显式传入候选池。不同 Method 的 index 会重复，
    # 因此不能再通过覆盖 iso.subgroups 来隐式切换上下文。
    source = data.get("source")
    if source == "method1":
        if not session.method1:
            raise ValueError("Method 1 candidate table is not available")
        candidates = [item.subgroup for item in session.method1]
    elif source == "method3":
        if not session.method3:
            raise ValueError("Method 3 candidate table is not available")
        selected_item = next(
            (
                item
                for item in session.method3
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
        candidates = [item.subgroup for item in session.method3]
    elif source == "subgroups":
        if not session.method2_subgroups:
            raise ValueError("Method 2 subgroup table is not available")
        candidates = list(session.method2_subgroups)
    elif source in (None, ""):
        # Backward compatibility for older clients that did not send source.
        candidates = list(session.method2_subgroups) or None
    else:
        raise ValueError(f"Unknown Method 2 candidate source: {source}")
    iso.set_distortion_scope(session.distortion_scope)
    iso.set_distortion_types(session.distortion_types)
    nmod = data["nmod"]
    result = iso.search_method_2(
        subgroup_idx=idx,
        distortion_type=data.get("distortion_type", session.distortion_types),
        number_of_independent_modulations=nmod,
        number_of_superposed_irreps=data.get("nsup", 1),
        candidates=candidates,
    )
    session.method2 = result
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


def _api_method3(session: WebSession, data: dict) -> dict:
    result = session.iso.search_method_3(
        distortion_types=data.get("distortion_types", session.distortion_types),
        point_group=data.get("point_group") or None,
        space_group_type=data.get("space_group_type") or None,
        supercell_basis=data.get("supercell_basis"),
        direct_sublattice_centering=data.get("direct_sublattice_centering") or None,
        lattice_type=data.get("lattice_type", "direct"),
        generate_if_missing=bool(data.get("generate", False)),
    )
    session.method3 = result
    return {"candidates": _method3_rows(result)}


def _api_method4(session: WebSession, data: dict) -> dict:
    content = data.get("content", "")
    if not content.strip():
        raise ValueError("Daughter CIF 内容为空 / daughter CIF content is empty")
    path = _write_upload(data.get("filename", "daughter.cif"), content)
    try:
        result = session.iso.search_method_4(
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


def _api_isotropy_cache_list(session: WebSession) -> dict:
    entries = session.iso.list_isotropy_cache()
    return {
        "entries": [e.to_dict() for e in entries],
        "count": len(entries),
    }


def _api_isotropy_cache_delete(session: WebSession, data: dict) -> dict:
    names = data.get("names") or []
    if not isinstance(names, list):
        raise ValueError("names must be a list of cache file names")
    result = session.iso.delete_isotropy_cache([str(n) for n in names])
    remaining = _api_isotropy_cache_list(session)
    return {**result, **remaining}


def _subgroups_for_method(session: WebSession, method: int) -> list:
    """取出所选 Method 当次计算得到的子群列表（不含其它 Method）。"""
    if method == 1:
        return [item.subgroup for item in session.method1]
    if method == 3:
        return [item.subgroup for item in session.method3]
    subs = list(session.method2_subgroups)
    if not subs and session.method2 is not None:
        subs = [session.method2.subgroup]
    return subs


@dataclass(frozen=True)
class RuntimeSettings:
    web_port: int
    web_idle_timeout: float


@dataclass(frozen=True)
class Download:
    filename: str
    body: bytes
    export_status: dict | None = None


class ExportFailedError(RuntimeError):
    """The selected export worker failed after its input was accepted."""


def runtime_settings() -> RuntimeSettings:
    cfg = get_config()
    return RuntimeSettings(cfg.web_port, cfg.web_idle_timeout)


def state_summary() -> dict:
    with _SESSION_LOCK:
        return copy.deepcopy(_state_summary(_SESSION))


def messages() -> dict:
    return copy.deepcopy(MESSAGES)


def space_groups() -> list[dict]:
    return copy.deepcopy(_space_groups())


def export_progress_snapshot() -> dict:
    return _EXPORT_JOBS.snapshot()


def cancel_export(job_id: str | None = None) -> bool:
    """Request cancellation only for the identified currently running export."""
    return _EXPORT_JOBS.cancel(job_id)


def reset_export_progress() -> None:
    _EXPORT_JOBS.reset()


def cleanup_uploads() -> None:
    with _SESSION_LOCK:
        _SESSION.cleanup_uploads()


def close_session() -> None:
    """Release the live session after transport request workers have stopped."""
    global _SESSION
    with _SESSION_LOCK:
        previous = _SESSION
        _SESSION = WebSession()
        try:
            previous.cleanup_uploads()
        finally:
            previous.close()


def _api_kpoints(session: WebSession, _data: dict) -> dict:
    return {"kpoints": [
        {"label": kp.label, "coordinates": kp.coordinates,
         "parameters": kp.parameters, "is_special": kp.is_special,
         "kovalev": kp.kovalev}
        for kp in session.iso.list_k_points()
    ]}


def _api_irreps(session: WebSession, data: dict) -> dict:
    return {"irreps": [
        {"label": ir.label, "dimension": ir.dimension, "active": ir.active}
        for ir in session.iso.list_irreps(data.get("k", ""), data.get("params"))
    ]}


def _api_method1_options(session: WebSession, _data: dict) -> dict:
    return {"options": session.iso.method1_options()}


_ACTIONS = {
    "load_cif": _api_load_cif,
    "set_types": _api_set_types,
    "method1": _api_method1,
    "subgroups": _api_subgroups,
    "method2": _api_method2,
    "method3": _api_method3,
    "method4": _api_method4,
    "isotropy_cache": lambda session, _data: _api_isotropy_cache_list(session),
    "isotropy_cache/delete": _api_isotropy_cache_delete,
    "kpoints": _api_kpoints,
    "irreps": _api_irreps,
    "method1_options": _api_method1_options,
}
_MUTATIONS = frozenset({
    "load_cif", "set_types", "method1", "subgroups", "method2", "method3", "method4",
})


def execute(action: str, data: dict | None = None, *, revision=None) -> dict:
    """Validate, compute against isolated state, then commit one new revision."""
    global _SESSION
    if action not in _ACTIONS:
        raise ValueError(f"Unknown action: {action}")
    values = action_input(action, data)
    if action in ("subgroups", "method2"):
        _require_single_method2_component(values)
    with _SESSION_LOCK:
        if action != "load_cif":
            _assert_current_revision(revision)
        if action not in _MUTATIONS:
            result = _ACTIONS[action](_SESSION, values) or {}
            result.setdefault("ok", True)
            result = response_data(result)
        else:
            previous = _SESSION
            staged = previous.stage(replace_parent=action == "load_cif")
            try:
                result = _ACTIONS[action](staged, values) or {}
                staged.revision = previous.revision + 1
                result["state"] = _state_summary(staged)
                result.setdefault("ok", True)
                result = response_data(result)
            except BaseException:
                try:
                    if staged.parent_upload_path != previous.parent_upload_path:
                        staged.cleanup_uploads()
                except OSError:
                    logging.getLogger(__name__).warning(
                        "Could not remove failed parent upload", exc_info=True,
                    )
                finally:
                    staged.close()
                raise
            _SESSION = staged
            previous.close()
            if previous.parent_upload_path != staged.parent_upload_path:
                try:
                    previous.cleanup_uploads()
                except OSError:
                    logging.getLogger(__name__).warning(
                        "Could not remove superseded parent upload", exc_info=True,
                    )
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
    """Validate a request and export a frozen pool under one identified job."""
    requested_revision = query_value(values, "revision")
    method = parse_export_method(query_value(values, "method", "2"))
    formats = parse_export_formats(query_value(values, "formats", "cif,isoviz,modes,topas"))
    if method == 2:
        nsup = nonnegative_integer(query_value(values, "nsup", "1"), "nsup", text=True)
        _require_single_method2_component({"nsup": nsup})
    export_nmod = validate_modulation_count(nonnegative_integer(query_value(values, "nmod", "0"), "nmod", text=True))
    job_id = export_job_id(query_value(values, "job_id"))
    compute = query_value(values, "compute_modes", "1").lower()
    if compute not in ("0", "false", "no", "1", "true", "yes"):
        raise ValueError("compute_modes must be 0/1, false/true or no/yes")
    compute_missing_modes = any(fmt != "cif" for fmt in formats) and compute in ("1", "true", "yes")
    indices_raw = query_value(values, "indices", "").strip()
    wanted = None
    if indices_raw:
        wanted = {nonnegative_integer(item.strip(), "indices", text=True)
                  for item in indices_raw.split(",")}
    worker = None
    job = None
    try:
        with _SESSION_LOCK:
            _assert_current_revision(requested_revision)
            live_iso = _SESSION.iso
            if live_iso.structure is None:
                raise FileNotFoundError("load a parent CIF first")
            subgroups = _subgroups_for_method(_SESSION, method)
            if wanted is not None:
                subgroups = [sg for sg in subgroups if sg.index in wanted]
            if not subgroups:
                raise FileNotFoundError(f"no Method {method} subgroups; run that Method first")
            job = _EXPORT_JOBS.start(job_id, revision=_SESSION.revision,
                                      method=method, total=len(subgroups))
            worker = live_iso.snapshot_for_export()
            subgroups = copy.deepcopy(subgroups)
        body = worker.export_subgroups_zip(
            formats=formats, subgroups=subgroups,
            compute_missing_modes=compute_missing_modes, wrapping=None,
            export_method=method, number_of_independent_modulations=export_nmod,
            progress_callback=lambda update: _EXPORT_JOBS.update(job, update),
            cancel_check=job.cancellation.is_set,
        )
        _close_worker(worker)
        worker = None
        status = _EXPORT_JOBS.finish(job)
        job = None
        return Download(f"isodistort_method{method}.zip", body, export_status=status)
    except BaseException as exc:
        if job is not None:
            _EXPORT_JOBS.finish(job, error=str(exc))
            if isinstance(exc, Exception):
                raise ExportFailedError(str(exc)) from exc
        raise
    finally:
        if worker is not None:
            _close_worker(worker)
