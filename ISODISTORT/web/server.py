"""
ISODISTORT local web UI (web/server.py)

Start a local HTTP server and open the browser:
    python main_web.py
    python web/server.py

Port defaults to 8000 (config runtime.web_port). If that port is taken the
server tries the next ports automatically.

The page is English only. It shares isocore.api.IsoDistort with the terminal
and the Python API.
"""

from __future__ import annotations

import contextlib
import json
import re
import sys
import threading
import time
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer, ThreadingHTTPServer
from pathlib import Path

import numpy as np

# 确保能导入 isocore（server.py 位于 ISODISTORT/web/）
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from isocore.api import IsoDistort  # noqa: E402
from isocore.distortion import DISTORTION_TYPES  # noqa: E402
from isocore.i18n import MESSAGES  # noqa: E402
from isocore.io import parse_export_formats, parse_export_method  # noqa: E402
from isocore.utils import get_config  # noqa: E402
from isocore.utils.schoenflies import hm_symbol, schoenflies_symbol  # noqa: E402

WEB_DIR = Path(__file__).resolve().parent


class WebSession:
    """网页会话：持有唯一 IsoDistort 实例（单用户本地工具）。

    IsoDistort 构造较慢（底层 BaseWrapper 在 Windows 下会初始化 WSL 短路径
    暂存目录与 ISODATA 符号链接，约 3-4 秒），采用懒初始化：首次访问
    ``iso`` 属性时才创建，避免拖慢模块导入（否则 main_web.py 启动会阻塞数秒）。
    """

    def __init__(self) -> None:
        self._iso = None
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


_SESSION = WebSession()

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


# ------------------------------------------------------------
# 生命周期管理：网页关闭 -> 自动停止服务并释放端口
# - 页面打开后周期性发送心跳（/api/ping）；关闭页面前发送
#   shutdown 信标（/api/shutdown）
# - 守护线程检测：收到 shutdown 请求，或“页面曾打开但心跳超时”，
#   则关闭 HTTPServer（释放端口）并退出进程
# ------------------------------------------------------------
class _Lifecycle:
    """页面生命周期状态（模块级单例，守护线程与请求处理器共享）。"""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.page_seen = False          # 页面是否至少打开过一次（未打开则服务常驻）
        self.page_heartbeat = 0.0       # 最近一次心跳/页面请求时间
        self.shutdown_requested = False
        self.in_flight = 0              # 进行中的长请求（ZIP 导出等）；>0 时不因心跳超时停服


_LIFE = _Lifecycle()


def _touch_heartbeat() -> None:
    """页面活动刷新（任何页面请求都算活动）。"""
    with _LIFE.lock:
        _LIFE.page_heartbeat = time.time()


def _mark_page_seen() -> None:
    with _LIFE.lock:
        _LIFE.page_seen = True


def _request_shutdown() -> None:
    """请求关闭服务（由 /api/shutdown 触发，稍后由守护线程执行）。"""
    with _LIFE.lock:
        _LIFE.shutdown_requested = True


def _begin_long_request() -> None:
    """ZIP 等长请求：刷新心跳并阻止看门狗在计算期间误杀服务。"""
    with _LIFE.lock:
        _LIFE.in_flight += 1
        _LIFE.page_heartbeat = time.time()


def _end_long_request() -> None:
    with _LIFE.lock:
        _LIFE.in_flight = max(0, _LIFE.in_flight - 1)
        _LIFE.page_heartbeat = time.time()


def _watchdog(server: HTTPServer, idle_timeout: float) -> None:
    """守护线程：页面关闭（心跳停止）或收到 shutdown 后关闭服务。"""
    while True:
        time.sleep(2)
        with _LIFE.lock:
            busy = _LIFE.in_flight > 0
            if busy:
                _LIFE.page_heartbeat = time.time()
            stale = (
                _LIFE.page_seen
                and not busy
                and (time.time() - _LIFE.page_heartbeat > idle_timeout)
            )
            stop = _LIFE.shutdown_requested or stale
        if stop:
            with contextlib.suppress(Exception):  # 关闭失败不影响退出
                server.shutdown()
            return


def _state_summary() -> dict:
    """返回当前会话状态摘要（供前端展示）。"""
    iso = _SESSION.iso
    modes = list(iso.mode_displacements.keys()) + list(iso.mode_occupancies.keys())
    summary = {
        "language": None,
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

    Method 3 结果项只含 subgroup / point_group / basis，没有 Method 1 的
    crystal_system / is_maximal，因此不能复用 _method1_rows（否则读不存在的
    属性抛 AttributeError）。Method 3 结果表格展示子群 + 点群。
    """
    rows = []
    for item in items:
        sg = item.subgroup
        rows.append({
            "index": sg.index,
            "space_group_number": sg.space_group_number,
            "space_group_symbol": sg.space_group_symbol,
            "k_point_label": sg.k_point_label,
            "irrep_label": sg.irrep_label,
            "opd_symbol": sg.opd_symbol,
            "point_group": item.point_group,
            "k_parameters": list(sg.k_parameters or []),
        })
    return rows


def _write_upload(filename: str, content: str) -> str:
    """把上传的 CIF 内容写入临时目录，返回文件路径。"""
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", Path(filename).name) or "upload.cif"
    cfg = get_config()
    upload_dir = cfg.temp_dir / "web_uploads"
    upload_dir.mkdir(parents=True, exist_ok=True)
    path = upload_dir / safe
    path.write_text(content, encoding="utf-8")
    return str(path)


class IsoHandler(BaseHTTPRequestHandler):
    """HTTP 请求处理器：JSON API + 静态页面。"""

    # ------------------------------------------------------------
    # 基础
    # ------------------------------------------------------------

    def log_message(self, fmt, *args):  # 精简日志
        sys.stdout.write("[web] " + fmt % args + "\n")

    def _send_json(self, data: dict, status: int = 200) -> None:
        body = json.dumps(data, ensure_ascii=False, default=str).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _read_json(self) -> dict:
        length = int(self.headers.get("Content-Length", 0))
        if length <= 0:
            return {}
        raw = self.rfile.read(length)
        return json.loads(raw.decode("utf-8"))

    def _run(self, fn) -> None:
        """执行 API 动作，统一错误处理（错误信息按当前语言输出）。"""
        try:
            result = fn()
            if result is None:
                result = {}
            result.setdefault("ok", True)
            self._send_json(result)
        except Exception as exc:  # noqa: BLE001 - web 边界：统一转为 JSON 错误
            self._send_json({"ok": False, "error": str(exc)}, status=200)

    # ------------------------------------------------------------
    # 路由
    # ------------------------------------------------------------

    def do_GET(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path in ("/", "/index.html"):
            _mark_page_seen()
            _touch_heartbeat()
            self._serve_index()
        elif path.startswith("/static/"):
            _touch_heartbeat()
            self._serve_static(path)
        elif path == "/api/ping":
            # 心跳：页面存活标记（关闭页面后心跳停止 -> 守护线程自动停服）
            _touch_heartbeat()
            self._send_json({"ok": True})
        elif path == "/api/state":
            _touch_heartbeat()
            self._send_json({"ok": True, "state": _state_summary()})
        elif path == "/api/i18n":
            _touch_heartbeat()
            self._send_json({
                "ok": True,
                "messages": MESSAGES,
            })
        elif path == "/api/kpoints":
            _touch_heartbeat()
            self._run(lambda: {
                "kpoints": [
                    {"label": kp.label, "coordinates": kp.coordinates,
                     "parameters": kp.parameters, "is_special": kp.is_special,
                     "kovalev": kp.kovalev}
                    for kp in _SESSION.iso.list_k_points()
                ],
            })
        elif path == "/api/space_groups":
            # 230 个空间群的 序号+HM 符号（Method 3 下拉，对齐官网表单）
            _touch_heartbeat()
            self._send_json({"ok": True, "space_groups": _space_groups()})
        elif path == "/api/irreps":
            qs = urllib.parse.parse_qs(parsed.query)
            k = (qs.get("k") or [""])[0]
            params = (qs.get("params") or [""])[0].split(",") if qs.get("params") else None
            self._run(lambda: {
                "irreps": [
                    {"label": ir.label, "dimension": ir.dimension, "active": ir.active}
                    for ir in _SESSION.iso.list_irreps(k, params)
                ],
            })
        elif path == "/api/method1_options":
            # Method 1 下拉数据（对齐官网：可达子群空间群 + Conventional/Primitive lattice）
            _touch_heartbeat()
            self._run(lambda: {"options": _SESSION.iso.method1_options()})
        elif path == "/api/isotropy_cache":
            _touch_heartbeat()
            self._run(lambda: self._api_isotropy_cache_list())
        elif path == "/api/download":
            self._serve_download(parsed.query)
        elif path == "/api/download_all":
            # 一键下载全部输出文件（打包为 ZIP）
            _touch_heartbeat()
            _begin_long_request()
            try:
                self._serve_download_all()
            finally:
                _end_long_request()
        else:
            self._send_json({"ok": False, "error": f"Unknown path: {path}"}, 404)

    def do_POST(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        data = self._read_json()

        if path == "/api/shutdown":
            # 网页关闭/用户点击“停止服务”：先应答，再由守护线程关闭服务释放端口
            _request_shutdown()
            self._send_json({"ok": True, "shutdown": True})
            return
        _touch_heartbeat()

        if path == "/api/load_cif":
            self._run(lambda: self._api_load_cif(data))
        elif path == "/api/set_types":
            self._run(lambda: self._api_set_types(data))
        elif path == "/api/method1":
            self._run(lambda: self._api_method1(data))
        elif path == "/api/subgroups":
            self._run(lambda: self._api_subgroups(data))
        elif path == "/api/method2":
            self._run(lambda: self._api_method2(data))
        elif path == "/api/method3":
            self._run(lambda: self._api_method3(data))
        elif path == "/api/method4":
            self._run(lambda: self._api_method4(data))
        elif path == "/api/isotropy_cache/delete":
            self._run(lambda: self._api_isotropy_cache_delete(data))
        else:
            self._send_json({"ok": False, "error": f"Unknown path: {path}"}, 404)

    # ------------------------------------------------------------
    # API 实现
    # ------------------------------------------------------------

    def _api_load_cif(self, data: dict) -> dict:
        content = data.get("content", "")
        filename = data.get("filename", "upload.cif")
        if not content.strip():
            raise ValueError("CIF 内容为空 / CIF content is empty")
        path = _write_upload(filename, content)
        _SESSION.iso.load_structure(path)
        _SESSION.iso.set_distortion_scope(_SESSION.distortion_scope)
        _SESSION.iso.set_distortion_types(_SESSION.distortion_types)
        _SESSION.method1, _SESSION.method2, _SESSION.method3 = [], None, []
        _SESSION.method2_subgroups = []
        return {"state": _state_summary()}

    def _api_set_types(self, data: dict) -> dict:
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
        return {"state": _state_summary()}

    def _api_method1(self, data: dict) -> dict:
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
            maximal_subgroup_only=bool(data.get("maximal_subgroup_only", False)),
        )
        _SESSION.method1 = result
        return {"candidates": _method1_rows(result), "state": _state_summary()}

    def _api_subgroups(self, data: dict) -> dict:
        # 对齐官网 Method 2：枚举指定 k 点（+ 参数）下全部 IR 的子群。
        # 多组 k 点（superposed IRs，nsup>1）由前端一次性以 kpoints 列表提交，
        # 后端在此合并全部 k 点组的子群并连续编号，供行点击正确回查子群。
        # 注意：_SESSION 是模块级单例（WebSession），BaseHTTPRequestHandler
        # 实例上并不存在该属性，误写 self._SESSION 会抛
        # "'IsoHandler' object has no attribute '_SESSION'"。
        gen = bool(data.get("generate", False))
        groups = data.get("kpoints")
        if groups:
            all_subs: list = []
            for grp in groups:
                all_subs.extend(_SESSION.iso.list_subgroups_at_kpoint(
                    grp["k"],
                    k_parameters=grp.get("params"),
                    generate_if_missing=gen,
                ))
            # 连续重编号：避免多 k 点组之间 index 冲突（行点击用 index 定位子群）
            for j, sg in enumerate(all_subs):
                sg.index = j
            _SESSION.iso.subgroups = all_subs
            _SESSION.method2_subgroups = list(all_subs)
            return {"subgroups": _subgroup_rows(all_subs), "state": _state_summary()}
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
        _SESSION.iso.subgroups = list(subs)
        _SESSION.method2_subgroups = list(subs)
        return {"subgroups": _subgroup_rows(subs), "state": _state_summary()}

    def _api_method2(self, data: dict) -> dict:
        idx = data.get("subgroup_idx")
        if idx is None:
            raise ValueError("subgroup_idx 缺失 / subgroup_idx missing")
        idx = int(idx)
        iso = _SESSION.iso
        # 按结果表来源选择候选列表（前端随行点击携带 source）。
        # 会话内 iso.subgroups 会被 method1/method3/subgroups 相互覆盖，
        # 不指定来源时点击旧表行可能选错子群。
        source = data.get("source")
        if source == "method1" and _SESSION.method1:
            iso.subgroups = [item.subgroup for item in _SESSION.method1]
        elif source == "method3" and _SESSION.method3:
            iso.subgroups = [item.subgroup for item in _SESSION.method3]
        # source 缺省/"subgroups"：沿用 _api_subgroups（k 点枚举）设置的列表
        iso.set_distortion_scope(_SESSION.distortion_scope)
        iso.set_distortion_types(_SESSION.distortion_types)
        result = iso.search_method_2(
            subgroup_idx=idx,
            distortion_type=data.get("distortion_type", _SESSION.distortion_types),
        )
        _SESSION.method2 = result
        modes = []
        for m in result.modes:
            modes.append({
                "irrep_label": m.irrep_label,
                "opd_symbol": m.opd_symbol,
                "mode_type": m.mode_type,
                "wyckoff_sites": sorted({b.wyckoff_letter for b in m.bush_modes}),
                "n_representatives": len(m.bush_modes),
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
        return {"modes": modes, "state": _state_summary()}

    def _api_method3(self, data: dict) -> dict:
        result = _SESSION.iso.search_method_3(
            distortion_types=data.get("distortion_types", _SESSION.distortion_types),
            point_group=data.get("point_group") or None,
            space_group_type=data.get("space_group_type") or None,
            supercell_basis=data.get("supercell_basis"),
            direct_sublattice_centering=data.get("direct_sublattice_centering") or None,
            lattice_type=data.get("lattice_type", "direct"),
        )
        _SESSION.method3 = result
        return {"candidates": _method3_rows(result), "state": _state_summary()}

    def _api_method4(self, data: dict) -> dict:
        content = data.get("content", "")
        if not content.strip():
            raise ValueError("Daughter CIF 内容为空 / daughter CIF content is empty")
        path = _write_upload(data.get("filename", "daughter.cif"), content)
        result = _SESSION.iso.search_method_4(
            distorted_cif_path=path,
            atom_matching_method=data.get("atom_matching_method", "nearest-site"),
            robust_distance_threshold=float(data.get("robust_distance_threshold", 0.25)),
            provided_origin_shift=data.get("provided_origin_shift"),
        )
        ranked = sorted(result.amplitudes.items(), key=lambda kv: abs(kv[1]), reverse=True)
        return {
            "amplitudes": {k: float(v) for k, v in ranked},
            "rms_residual": result.rms_residual,
            "max_abs_residual": result.max_abs_residual,
        }

    def _api_isotropy_cache_list(self) -> dict:
        from isocore.backend.isotropy_cache import list_isotropy_cache

        wrapper = _SESSION.iso._iso
        entries = list_isotropy_cache(wrapper)
        return {
            "entries": [e.to_dict() for e in entries],
            "count": len(entries),
        }

    def _api_isotropy_cache_delete(self, data: dict) -> dict:
        from isocore.backend.isotropy_cache import delete_isotropy_cache

        names = data.get("names") or []
        if not isinstance(names, list):
            raise ValueError("names must be a list of cache file names")
        wrapper = _SESSION.iso._iso
        result = delete_isotropy_cache(wrapper, [str(n) for n in names])
        remaining = self._api_isotropy_cache_list()
        return {**result, **remaining}

    # ------------------------------------------------------------
    # 静态文件 / 下载
    # ------------------------------------------------------------

    def _serve_index(self) -> None:
        index = WEB_DIR / "index.html"
        body = index.read_bytes() if index.exists() else b"<h1>index.html missing</h1>"
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _serve_static(self, path: str) -> None:
        """提供 web/static/ 下的静态资源（bootstrap.css / docs.css / help.jpg 等）。"""
        rel = path[len("/static/"):]
        file_path = (WEB_DIR / "static" / rel).resolve()
        if not str(file_path).startswith(str((WEB_DIR / "static").resolve())):
            self._send_json({"ok": False, "error": "invalid path"}, 403)
            return
        if not file_path.is_file():
            self._send_json({"ok": False, "error": "not found"}, 404)
            return
        mime = "image/jpeg" if file_path.suffix.lower() in (".jpg", ".jpeg") else \
            "text/css" if file_path.suffix.lower() == ".css" else \
            "application/octet-stream"
        body = file_path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", f"{mime}; charset=utf-8" if "text/" in mime else mime)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _serve_download(self, query: str) -> None:
        qs = urllib.parse.parse_qs(query)
        fname = (qs.get("file") or [""])[0]
        cfg = get_config()
        path = (cfg.output_dir / fname).resolve()
        # 只允许输出目录内的文件
        if not str(path).startswith(str(cfg.output_dir.resolve())):
            self._send_json({"ok": False, "error": "invalid file"}, 403)
            return
        if not path.exists():
            self._send_json({"ok": False, "error": "file not found"}, 404)
            return
        body = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Disposition", f'attachment; filename="{path.name}"')
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _subgroups_for_method(self, method: int) -> list:
        """取出所选 Method 当次计算得到的子群列表（不含其它 Method）。"""
        if method == 1:
            return [item.subgroup for item in _SESSION.method1]
        if method == 3:
            return [item.subgroup for item in _SESSION.method3]
        subs = list(_SESSION.method2_subgroups)
        if not subs and _SESSION.method2 is not None:
            subs = [_SESSION.method2.subgroup]
        return subs

    def _serve_download_all(self) -> None:
        """按用户所选的**一个** Method 的子群打包导出（不扫描 output_dir）。

        查询参数：
            ``method``：1 / 2 / 3（不可多选；缺省 2）
            ``formats``：cif,isoviz,modes,topas（官网第 6 页对应选项）
        ZIP 结构：每子群一个文件夹（Method 1 为完整 OPD 行；Method 2/3 为 ``IR OPD``），
        内含所选格式文件；外层 ZIP 文件名为 ``isodistort_methodN.zip``（``wrapping=None``）。
        查询参数 ``indices``：逗号分隔的子群 index；若提供，只打包这些子群
        （网页在当前 Method 有筛选时传入命中行）。
        """
        parsed = urllib.parse.urlparse(self.path)
        qs = urllib.parse.parse_qs(parsed.query)
        method_vals = qs.get("method") or ["2"]
        try:
            if len(method_vals) > 1:
                raise ValueError(
                    "只能选择一个 Method 导出，不能多选 / select exactly one Method"
                )
            method = parse_export_method(method_vals[0])
            fmts = parse_export_formats(
                (qs.get("formats") or ["cif,isoviz,modes,topas"])[0]
            )
        except ValueError as exc:
            self._send_json({"ok": False, "error": str(exc)}, 400)
            return
        iso = _SESSION.iso
        if iso.structure is None:
            self._send_json(
                {"ok": False, "error": "请先加载母相 CIF / load a parent CIF first"},
                404,
            )
            return
        subs = self._subgroups_for_method(method)
        indices_raw = (qs.get("indices") or [""])[0].strip()
        if indices_raw:
            try:
                want = {int(x.strip()) for x in indices_raw.split(",") if x.strip()}
            except ValueError:
                self._send_json(
                    {"ok": False, "error": "indices must be comma-separated integers"},
                    400,
                )
                return
            subs = [sg for sg in subs if sg.index in want]
        if not subs:
            self._send_json({
                "ok": False,
                "error": (
                    f"没有可导出的 Method {method} 子群；请先完成该 Method 的计算"
                    f" / no Method {method} subgroups; run that Method first"
                ),
            }, 404)
            return
        # 勾选了 isoviz / modes / topas 时，对子群补跑 Method 2 以填充模式
        # （长计算由网页 busy 进度条提示）。参数 k 点本地无法算位移模式。
        need_modes = any(fmt != "cif" for fmt in fmts)
        compute_q = (qs.get("compute_modes") or ["1"])[0].strip().lower()
        # 默认开启；显式 compute_modes=0 可跳过（仅结构骨架，速度快）
        want_compute = compute_q not in ("0", "false", "no")
        compute_missing_modes = need_modes and want_compute
        saved_subs = list(iso.subgroups)
        try:
            iso.subgroups = list(subs)
            body = iso.export_subgroups_zip(
                formats=fmts,
                subgroups=subs,
                compute_missing_modes=compute_missing_modes,
                wrapping=None,
                use_opd_line_folders=(method == 1),
            )
        except Exception as exc:  # noqa: BLE001 - web 边界：统一转为 JSON 错误
            self._send_json({"ok": False, "error": str(exc)}, 500)
            return
        finally:
            iso.subgroups = saved_subs
        self.send_response(200)
        self.send_header("Content-Type", "application/zip")
        self.send_header(
            "Content-Disposition",
            f'attachment; filename="isodistort_method{method}.zip"',
        )
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def _bind_server(host: str, preferred_port: int) -> ThreadingHTTPServer | None:
    """按“配置端口 → 顺延端口 → 系统空闲端口(0)”的顺序尝试绑定，保证成功。

    使用 ThreadingHTTPServer：在线生成子群数据库等长耗时请求会阻塞较久，
    单线程 HTTPServer 会同时阻塞心跳，导致看门狗误判“页面关闭”并停服；
    多线程版本可让心跳/状态请求在长请求期间继续正常响应。
    """
    candidates = list(range(preferred_port, preferred_port + 21))
    candidates.append(0)  # 交给系统分配空闲端口，兜底保证可启动
    for port in candidates:
        try:
            return ThreadingHTTPServer((host, port), IsoHandler)
        except OSError:
            continue
    return None


def _open_browser(url: str) -> None:
    """自动打开默认浏览器；失败时重试一次，并在控制台醒目输出网址。"""
    opened = False
    for opener in (webbrowser.open, webbrowser.open_new):
        try:
            if opener(url):
                opened = True
                break
        except Exception:  # noqa: BLE001, S112 - 浏览器异常不应影响服务运行，继续尝试
            continue
    if not opened:
        print(f"\n无法自动打开浏览器，请手动访问: {url}")


def main() -> int:
    cfg = get_config()
    host = "127.0.0.1"

    server = _bind_server(host, cfg.web_port)
    if server is None:
        print(f"无法绑定端口（{cfg.web_port} 及顺延端口均被占用），请检查网络环境。")
        return 1
    port = server.server_address[1]

    url = f"http://{host}:{port}/"
    # flush=True：即使输出被重定向/记录，网址也立即可见
    print("=" * 60, flush=True)
    print("ISODISTORT Local Web Console", flush=True)
    print(f"  URL: {url}", flush=True)
    idle = cfg.web_idle_timeout
    print(f"  Auto-stop: shuts down ~{idle}s after the page is closed", flush=True)
    print("  Press Ctrl+C to stop", flush=True)
    print("=" * 60, flush=True)

    # 延迟自动打开浏览器（等待服务就绪）；无论成败都会在控制台给出网址
    threading.Timer(1.0, _open_browser, args=[url]).start()

    # 守护线程：网页关闭/心跳停止后自动停服并释放端口
    watchdog = threading.Thread(
        target=_watchdog, args=(server, float(idle)), daemon=True,
        name="isodistort-web-watchdog",
    )
    watchdog.start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n服务已停止 / Server stopped.")
    finally:
        server.server_close()
    print("已退出并释放端口 / Exited, port released.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
