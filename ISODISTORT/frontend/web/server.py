"""
ISODISTORT local web UI (frontend/web/server.py)

Start a local HTTP server and open the browser:
    python main_web.py
    python frontend/web/server.py

Port defaults to 8000 (config runtime.web_port). If that port is taken the
server tries the next ports automatically.

The page is English only. It calls application workflows shared with the Python API.
"""

from __future__ import annotations

import contextlib
import json
import shutil
import subprocess
import sys
import threading
import time
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer, ThreadingHTTPServer
from pathlib import Path

# Keep the package root available when this transport is launched directly.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from application import workflows  # noqa: E402

WEB_DIR = Path(__file__).resolve().parent

_PAGE_CLOSE_GRACE = 1.5
_WATCHDOG_POLL_SECONDS = 0.25


class _Lifecycle:
    """页面生命周期状态（模块级单例，守护线程与请求处理器共享）。"""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.page_seen = False          # 页面是否至少打开过一次（未打开则服务常驻）
        self.page_heartbeat = 0.0       # 最近一次心跳/页面请求时间
        self.shutdown_requested = False
        self.in_flight = 0              # 进行中的长请求（ZIP 导出等）；>0 时不因心跳超时停服
        self.clients: dict[str, float] = {}  # 活跃标签页 client_id -> 最近心跳
        self.last_page_closed = 0.0     # 最后一个标签显式关闭的时间；0 表示仍有页面


_LIFE = _Lifecycle()


def _reset_lifecycle() -> None:
    """重置单进程生命周期状态（正式启动及测试隔离共用）。"""
    with _LIFE.lock:
        _LIFE.page_seen = False
        _LIFE.page_heartbeat = 0.0
        _LIFE.shutdown_requested = False
        _LIFE.in_flight = 0
        _LIFE.clients.clear()
        _LIFE.last_page_closed = 0.0


def _touch_heartbeat(client_id: str | None = None) -> None:
    """刷新服务活动；带 client_id 的 ping 同时登记一个活跃标签页。"""
    now = time.time()
    with _LIFE.lock:
        _LIFE.page_heartbeat = now
        if client_id:
            _LIFE.page_seen = True
            _LIFE.clients[client_id] = now
            _LIFE.last_page_closed = 0.0


def _mark_page_seen() -> None:
    with _LIFE.lock:
        _LIFE.page_seen = True


def _request_shutdown() -> None:
    """请求关闭服务（由 /api/shutdown 触发，稍后由守护线程执行）。"""
    with _LIFE.lock:
        _LIFE.shutdown_requested = True


def _mark_page_closed(client_id: str | None) -> int:
    """注销一个标签页；返回仍活跃的页面数。重复信标是幂等的。"""
    with _LIFE.lock:
        if client_id:
            _LIFE.clients.pop(client_id, None)
        else:
            # 兼容没有 client_id 的旧页面：只有无法区分标签时才视为全部关闭。
            _LIFE.clients.clear()
        if not _LIFE.clients:
            _LIFE.last_page_closed = time.time()
        return len(_LIFE.clients)


def _begin_long_request() -> None:
    """ZIP 等长请求：刷新心跳并阻止看门狗在计算期间误杀服务。"""
    with _LIFE.lock:
        _LIFE.in_flight += 1
        _LIFE.page_heartbeat = time.time()


def _end_long_request() -> None:
    with _LIFE.lock:
        _LIFE.in_flight = max(0, _LIFE.in_flight - 1)
        _LIFE.page_heartbeat = time.time()


def _watchdog(
    server: HTTPServer,
    idle_timeout: float,
    close_grace: float = _PAGE_CLOSE_GRACE,
    poll_interval: float = _WATCHDOG_POLL_SECONDS,
) -> None:
    """最后页面关闭后停服；心跳超时负责兜底浏览器异常退出。"""
    while True:
        time.sleep(poll_interval)
        with _LIFE.lock:
            busy = _LIFE.in_flight > 0
            if busy:
                _LIFE.page_heartbeat = time.time()
            explicitly_closed = (
                _LIFE.page_seen
                and not busy
                and not _LIFE.clients
                and _LIFE.last_page_closed > 0
                and (time.time() - _LIFE.last_page_closed >= close_grace)
            )
            stale = (
                _LIFE.page_seen
                and not busy
                and (time.time() - _LIFE.page_heartbeat > idle_timeout)
            )
            stop = _LIFE.shutdown_requested or explicitly_closed or stale
        if stop:
            with contextlib.suppress(Exception):  # 关闭失败不影响退出
                server.shutdown()
            return


class IsoHandler(BaseHTTPRequestHandler):
    """HTTP transport: decode requests, invoke workflows, encode responses."""

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

    def _run(self, action: str, data: dict | None = None, *, revision=None) -> None:
        try:
            self._send_json(workflows.execute(action, data, revision=revision))
        except workflows.StaleContextError as exc:
            self._send_json({
                "ok": False, "stale_context": True, "error": str(exc),
                "state": workflows.state_summary(),
            }, status=409)
        except Exception as exc:  # noqa: BLE001 - HTTP error boundary
            self._send_json({"ok": False, "error": str(exc)}, status=200)

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
            # 每个标签页独立续租；最后一个标签关闭才会触发快速停服。
            qs = urllib.parse.parse_qs(parsed.query)
            client_id = (qs.get("client_id") or [""])[0].strip() or None
            _touch_heartbeat(client_id)
            with _LIFE.lock:
                active_clients = len(_LIFE.clients)
            self._send_json({"ok": True, "active_clients": active_clients})
        elif path == "/api/state":
            _touch_heartbeat()
            # Build one coherent snapshot while other ThreadingHTTPServer
            # handlers may be loading a CIF or replacing candidate pools.
            state = workflows.state_summary()
            self._send_json({"ok": True, "state": state})
        elif path == "/api/export_status":
            _touch_heartbeat()
            self._send_json({"ok": True, "export": workflows.export_progress_snapshot()})
        elif path == "/api/messages":
            _touch_heartbeat()
            self._send_json({
                "ok": True,
                "messages": workflows.messages(),
            })
        elif path == "/api/kpoints":
            _touch_heartbeat()
            self._run("kpoints")
        elif path == "/api/space_groups":
            # 230 个空间群的 序号+HM 符号（Method 3 下拉，对齐官网表单）
            _touch_heartbeat()
            self._send_json({"ok": True, "space_groups": workflows.space_groups()})
        elif path == "/api/irreps":
            qs = urllib.parse.parse_qs(parsed.query)
            k = (qs.get("k") or [""])[0]
            params = (qs.get("params") or [""])[0].split(",") if qs.get("params") else None
            self._run("irreps", {"k": k, "params": params})
        elif path == "/api/method1_options":
            # Method 1 下拉数据（对齐官网：可达子群空间群 + Conventional/Primitive lattice）
            _touch_heartbeat()
            self._run("method1_options")
        elif path == "/api/isotropy_cache":
            _touch_heartbeat()
            self._run("isotropy_cache")
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
        if path == "/api/export_cancel":
            requested = workflows.cancel_export()
            self._send_json({"ok": True, "cancellation_requested": requested})
            return
        if path == "/api/client/close":
            qs = urllib.parse.parse_qs(parsed.query)
            client_id = (qs.get("client_id") or [""])[0].strip() or None
            active_clients = _mark_page_closed(client_id)
            self._send_json({
                "ok": True,
                "closing": active_clients == 0,
                "active_clients": active_clients,
            })
            return
        _touch_heartbeat()

        if path == "/api/load_cif":
            # Loading a parent explicitly replaces the shared context, even
            # when another tab has advanced it since this tab's last snapshot.
            self._run("load_cif", data)
        elif path == "/api/set_types":
            self._run(
                "set_types", data,
                revision=data.get("revision"),
            )
        elif path == "/api/method1":
            self._run(
                "method1", data,
                revision=data.get("revision"),
            )
        elif path == "/api/subgroups":
            self._run(
                "subgroups", data,
                revision=data.get("revision"),
            )
        elif path == "/api/method2":
            self._run(
                "method2", data,
                revision=data.get("revision"),
            )
        elif path == "/api/method3":
            self._run(
                "method3", data,
                revision=data.get("revision"),
            )
        elif path == "/api/method4":
            self._run(
                "method4", data,
                revision=data.get("revision"),
            )
        elif path == "/api/isotropy_cache/delete":
            self._run("isotropy_cache/delete", data)
        else:
            self._send_json({"ok": False, "error": f"Unknown path: {path}"}, 404)

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
        """提供 frontend/web/static/ 下的静态资源。"""
        rel = path[len("/static/"):]
        static_root = (WEB_DIR / "static").resolve()
        file_path = (static_root / rel).resolve()
        if not file_path.is_relative_to(static_root):
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

    def _send_download(self, download: workflows.Download, content_type: str) -> None:
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Disposition", f'attachment; filename="{download.filename}"')
        self.send_header("Content-Length", str(len(download.body)))
        self.end_headers()
        self.wfile.write(download.body)

    def _serve_download(self, query: str) -> None:
        values = urllib.parse.parse_qs(query)
        try:
            download = workflows.read_download((values.get("file") or [""])[0])
        except PermissionError as exc:
            self._send_json({"ok": False, "error": str(exc)}, 403)
            return
        except FileNotFoundError as exc:
            self._send_json({"ok": False, "error": str(exc)}, 404)
            return
        self._send_download(download, "application/octet-stream")

    def _serve_download_all(self) -> None:
        values = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
        try:
            download = workflows.export_archive(values)
        except workflows.StaleContextError as exc:
            self._send_json({
                "ok": False, "stale_context": True, "error": str(exc),
                "state": workflows.state_summary(),
            }, 409)
            return
        except ValueError as exc:
            self._send_json({"ok": False, "error": str(exc)}, 400)
            return
        except FileNotFoundError as exc:
            self._send_json({"ok": False, "error": str(exc)}, 404)
            return
        except workflows.ExportBusyError as exc:
            self._send_json({
                "ok": False, "error": str(exc),
                "export": workflows.export_progress_snapshot(),
            }, 409)
            return
        except workflows.ExportFailedError as exc:
            self._send_json({"ok": False, "error": str(exc)}, 500)
            return
        self._send_download(download, "application/zip")

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

def _open_browser_from_windows_shell(url: str) -> bool:
    """Ask Windows Shell to open *url* through a hidden PowerShell helper."""
    if not sys.platform.startswith("win"):
        return False
    powershell = shutil.which("powershell.exe")
    if not powershell:
        return False
    quoted_url = "'" + url.replace("'", "''") + "'"
    try:
        completed = subprocess.run(  # noqa: S603 - fixed Windows PowerShell executable
            [
                powershell,
                "-NoProfile",
                "-WindowStyle",
                "Hidden",
                "-Command",
                f"Start-Process -FilePath {quoted_url}",
            ],
            check=False,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return completed.returncode == 0

def _open_browser(url: str) -> bool:
    """Request the default browser and report whether dispatch succeeded."""
    if _open_browser_from_windows_shell(url):
        print(f"Browser launch requested: {url}", flush=True)
        return True
    opened = False
    for opener in (webbrowser.open, webbrowser.open_new):
        try:
            if opener(url):
                opened = True
                break
        except Exception:  # noqa: BLE001, S112 - 浏览器异常不应影响服务运行，继续尝试
            continue
    if not opened:
        print(f"\n无法自动打开浏览器，请手动访问: {url}", flush=True)
        return False
    print(f"Browser launch requested: {url}", flush=True)
    return True

def main() -> int:
    cfg = workflows.runtime_settings()
    host = "127.0.0.1"
    _reset_lifecycle()
    workflows.reset_export_progress()

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
    print(
        "  Auto-stop: exits shortly after the last tab closes "
        f"({idle}s heartbeat fallback)",
        flush=True,
    )
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
        workflows.cleanup_uploads()
    print("已退出并释放端口 / Exited, port released.", flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
