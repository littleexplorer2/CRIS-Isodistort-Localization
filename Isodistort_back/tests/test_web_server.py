"""
测试：web 服务（web/server.py）基础端点冒烟。

不依赖 WSL 的端点（/api/i18n、/api/state、/）直接验证；
依赖真实计算的端点（load_cif 等）在 WSL 可用时验证。
"""
import json
import shutil
import subprocess
import sys
import threading
import time
import urllib.request
from pathlib import Path

import pytest

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))
from web import server as web_server  # noqa: E402  (导入即校验语法/路由)


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


class _ServerHandle:
    """在测试进程内启动/停止 HTTPServer。"""

    def __init__(self, port: int = 0) -> None:
        self.httpd = web_server.HTTPServer(("127.0.0.1", port), web_server.IsoHandler)
        self.port = self.httpd.server_address[1]

    def __enter__(self):
        self._thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self.httpd.shutdown()
        self.httpd.server_close()


@pytest.fixture(scope="module")
def server():
    with _ServerHandle() as handle:
        time.sleep(0.2)
        yield handle


def _get(port: int, path: str) -> dict:
    with urllib.request.urlopen(f"http://127.0.0.1:{port}{path}") as resp:
        return json.loads(resp.read().decode("utf-8"))


def test_index_served(server):
    with urllib.request.urlopen(f"http://127.0.0.1:{server.port}/") as resp:
        body = resp.read().decode("utf-8")
    assert resp.status == 200
    assert "ISODISTORT" in body
    assert "langBtn" in body  # 语言切换按钮存在


def test_i18n_endpoint(server):
    data = _get(server.port, "/api/i18n?lang=en")
    assert data["ok"]
    assert data["language"] == "en"
    assert data["messages"]["load.done"].startswith("[Loaded]")
    assert data["terms"]["space group"] == "空间群"

    data_zh = _get(server.port, "/api/i18n?lang=zh")
    assert data_zh["messages"]["load.done"].startswith("[加载完成]")

    # mixed：返回中文文案 + 术语表（前端做专有名词替换）
    data_mixed = _get(server.port, "/api/i18n?lang=mixed")
    assert data_mixed["language"] == "mixed"
    assert data_mixed["messages"]["load.done"].startswith("[加载完成]")
    assert data_mixed["terms"]["space group"] == "空间群"


def test_state_endpoint(server):
    data = _get(server.port, "/api/state")
    assert data["ok"]
    assert data["state"]["structure"] is None


@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_load_cif_endpoint(server):
    """上传 EuAl4 母相 CIF 并校验识别结果（真实计算）。"""
    cif = (Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码")
           / "EuAl4 Springer (parent).cif")
    if not cif.exists():
        pytest.skip("EuAl4 CIF 不存在")
    body = json.dumps({"filename": "eual4.cif", "content": cif.read_text(encoding="utf-8")}).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/load_cif?lang=zh", data=body,
        headers={"Content-Type": "application/json"},
    )
    # 仅访问本地测试服务（http scheme），无 file:// 等自定义 scheme
    with urllib.request.urlopen(req) as resp:  # noqa: S310
        data = json.loads(resp.read().decode("utf-8"))
    assert data["ok"]
    assert data["state"]["structure"]["space_group_number"] == 139
