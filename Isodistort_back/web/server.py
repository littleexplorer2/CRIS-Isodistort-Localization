"""
ISODISTORT 本地网页交互程序（web/server.py）

在本地端口启动简易网页界面，并自动打开浏览器。
用法（两者等价）：
    python main_web.py
    python web/server.py
    （端口默认 8000，可用配置 runtime.web_port 修改；被占用时自动顺延）

依赖：仅 Python 标准库（http.server），无需额外安装。

与终端/API 的关系：
- 底层复用 isocore.api.IsoDistort（同一套真实 iso/findsym 计算）
- 界面语言：页面右上角 中/EN/中+EN 按钮切换（前端渲染）；服务器控制台输出
  跟随请求的 ?lang= 参数（未指定时用配置 runtime.language）
"""

from __future__ import annotations

import json
import re
import sys
import threading
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

# 确保能导入 isocore（server.py 位于 Isodistort_back/web/）
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from isocore.api import IsoDistort  # noqa: E402
from isocore.i18n import MESSAGES, TERMS_EN2ZH, set_language  # noqa: E402
from isocore.utils import get_config  # noqa: E402

WEB_DIR = Path(__file__).resolve().parent


class WebSession:
    """网页会话：持有唯一 IsoDistort 实例（单用户本地工具）。"""

    def __init__(self) -> None:
        self.iso = IsoDistort()
        self.distortion_types: list[str] = ["displacement", "strain"]
        self.method1: list = []
        self.method2 = None
        self.method3: list = []


_SESSION = WebSession()


def _state_summary() -> dict:
    """返回当前会话状态摘要（供前端展示）。"""
    iso = _SESSION.iso
    summary = {
        "language": None,
        "structure": None,
        "subgroups": len(iso.subgroups),
        "modes": list(iso.mode_displacements.keys()),
        "distorted_atoms": len(iso.distorted_structure) if iso.distorted_structure else None,
        "distortion_types": _SESSION.distortion_types,
    }
    if iso.structure is not None:
        summary["structure"] = {
            "space_group_number": iso.symmetry_info["space_group_number"],
            "space_group_symbol": iso.symmetry_info["space_group_symbol"],
            "atoms": len(iso.structure),
            "wyckoff": [
                {"letter": s["wyckoff_letter"], "species": s["species"]}
                for s in iso.symmetry_info["wyckoff_sites"]
            ],
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
        })
    return rows


def _method1_rows(items) -> list[dict]:
    rows = []
    for item in items:
        sg = item.subgroup
        rows.append({
            "index": sg.index,
            "space_group_number": sg.space_group_number,
            "space_group_symbol": sg.space_group_symbol,
            "crystal_system": item.crystal_system,
            "is_maximal": item.is_maximal,
            "k_point_label": sg.k_point_label,
            "irrep_label": sg.irrep_label,
            "opd_symbol": sg.opd_symbol,
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

    def _query_lang(self) -> str:
        """从查询参数取语言（?lang=zh|en|mixed），并应用到服务器端输出。

        未指定时使用配置 runtime.language（默认 en）。
        """
        parsed = urllib.parse.urlparse(self.path)
        qs = urllib.parse.parse_qs(parsed.query)
        default = get_config().language
        lang = (qs.get("lang") or [default])[0]
        if lang not in ("zh", "en", "mixed"):
            lang = default if default in ("zh", "en", "mixed") else "en"
        try:
            set_language(lang)
        except ValueError:
            pass
        return lang

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
        lang = self._query_lang()
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path in ("/", "/index.html"):
            self._serve_index()
        elif path == "/api/state":
            self._send_json({"ok": True, **{"state": _state_summary()},
                             "language": lang})
        elif path == "/api/i18n":
            # mixed 模式：返回中文文案 + 术语表，前端自行做专有名词替换
            messages = MESSAGES["zh"] if lang == "mixed" else MESSAGES.get(lang, MESSAGES["zh"])
            self._send_json({
                "ok": True,
                "language": lang,
                "messages": messages,
                "terms": TERMS_EN2ZH,
            })
        elif path == "/api/kpoints":
            self._run(lambda: {
                "kpoints": [
                    {"label": kp.label, "coordinates": kp.coordinates,
                     "parameters": kp.parameters, "is_special": kp.is_special}
                    for kp in _SESSION.iso.list_k_points()
                ],
            })
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
        elif path == "/api/download":
            self._serve_download(parsed.query)
        else:
            self._send_json({"ok": False, "error": f"Unknown path: {path}"}, 404)

    def do_POST(self) -> None:
        _ = self._query_lang()  # 按请求设置服务器端语言（含控制台输出语言）
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        data = self._read_json()

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
        elif path == "/api/generate":
            self._run(lambda: self._api_generate(data))
        elif path == "/api/export":
            self._run(lambda: self._api_export(data))
        elif path == "/api/mixed":
            self._run(lambda: self._api_mixed(data))
        elif path == "/api/domains":
            self._run(self._api_domains)
        elif path == "/api/set_language":
            lang = data.get("language", get_config().language)
            set_language(lang)
            self._send_json({"ok": True, "language": lang})
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
        _SESSION.method1, _SESSION.method2, _SESSION.method3 = [], None, []
        return {"state": _state_summary()}

    def _api_set_types(self, data: dict) -> dict:
        types = data.get("types", ["displacement", "strain"])
        valid = {"displacement", "order", "strain", "magnetic"}
        _SESSION.distortion_types = [t for t in types if t in valid] or ["displacement", "strain"]
        return {"state": _state_summary()}

    def _api_method1(self, data: dict) -> dict:
        result = _SESSION.iso.search_method_1(
            distortion_types=data.get("distortion_types", _SESSION.distortion_types),
            crystal_system=data.get("crystal_system") or None,
            subgroup_space_group=data.get("subgroup_space_group") or None,
            direct_sublattice=data.get("direct_sublattice"),
            maximal_subgroup_only=bool(data.get("maximal_subgroup_only", False)),
        )
        _SESSION.method1 = result
        return {"candidates": _method1_rows(result), "state": _state_summary()}

    def _api_subgroups(self, data: dict) -> dict:
        subs = _SESSION.iso.list_subgroups_at(
            data["k"], data["ir"],
            k_parameters=data.get("params"),
            opd_symbol=data.get("opd"),
            generate_if_missing=bool(data.get("generate", False)),
        )
        return {"subgroups": _subgroup_rows(subs), "state": _state_summary()}

    def _api_method2(self, data: dict) -> dict:
        idx = int(data["subgroup_idx"])
        result = _SESSION.iso.search_method_2(
            subgroup_idx=idx,
            distortion_type=data.get("distortion_type", "displacement"),
        )
        _SESSION.method2 = result
        modes = []
        for m in result.modes:
            modes.append({
                "irrep_label": m.irrep_label,
                "opd_symbol": m.opd_symbol,
                "wyckoff_sites": sorted({b.wyckoff_letter for b in m.bush_modes}),
                "n_representatives": len(m.bush_modes),
            })
        return {"modes": modes, "state": _state_summary()}

    def _api_method3(self, data: dict) -> dict:
        result = _SESSION.iso.search_method_3(
            distortion_types=data.get("distortion_types", _SESSION.distortion_types),
            point_group=data.get("point_group") or None,
            space_group_type=data.get("space_group_type") or None,
            supercell_basis=data.get("supercell_basis"),
            direct_sublattice_centering=data.get("direct_sublattice_centering") or None,
        )
        _SESSION.method3 = result
        return {"candidates": _method1_rows(result), "state": _state_summary()}

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

    def _api_generate(self, data: dict) -> dict:
        distorted = _SESSION.iso.generate_distortion(
            irrep_label=data.get("irrep_label"),
            amplitude=float(data.get("amplitude", 1.0)) if data.get("amplitude") is not None else None,
            supercell=data.get("supercell"),
        )
        return {
            "atoms": len(distorted),
            "volume": distorted.volume,
            "state": _state_summary(),
        }

    def _api_mixed(self, data: dict) -> dict:
        distorted = _SESSION.iso.generate_mixed_distortion(
            contributions={k: float(v) for k, v in (data.get("contributions") or {}).items()},
            supercell=data.get("supercell"),
        )
        return {"atoms": len(distorted), "volume": distorted.volume,
                "state": _state_summary()}

    def _api_export(self, data: dict) -> dict:
        filename = data.get("filename", "web_distorted")
        formats = data.get("formats", ["cif"])
        paths = _SESSION.iso.export(filename, formats=formats)
        return {"files": [Path(p).name for p in paths],
                "state": _state_summary()}

    def _api_domains(self) -> dict:
        domains = _SESSION.iso.generate_domains()
        return {
            "domains": [
                {"domain_number": d.domain_number, "generator": d.generator,
                 "space_group_number": d.space_group_number,
                 "space_group_symbol": d.space_group_symbol}
                for d in domains
            ],
            "state": _state_summary(),
        }

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


def main() -> int:
    cfg = get_config()
    port = cfg.web_port
    host = "127.0.0.1"

    server = None
    for _attempt in range(10):
        try:
            server = HTTPServer((host, port), IsoHandler)
            break
        except OSError:
            port += 1
    if server is None:
        print(f"无法绑定端口 {cfg.web_port}-{port}，请检查占用。")
        return 1

    url = f"http://{host}:{port}/"
    print("=" * 60)
    print("ISODISTORT Local Web Console")
    print(f"  URL: {url}")
    print(f"  Language: {cfg.language}  (click the top-right button to switch)")
    print("  Press Ctrl+C to stop")
    print("=" * 60)

    # 自动打开浏览器
    threading.Timer(1.0, webbrowser.open, args=[url]).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n服务已停止 / Server stopped.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
