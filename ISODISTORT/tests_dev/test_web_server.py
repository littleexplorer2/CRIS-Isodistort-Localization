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
import urllib.error
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


def _post(port: int, path: str, payload: dict | None = None) -> dict:
    body = json.dumps(payload or {}).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}", data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:  # noqa: S310 - 仅访问本地测试服务（http scheme）
        return json.loads(resp.read().decode("utf-8"))


def test_index_served(server):
    with urllib.request.urlopen(f"http://127.0.0.1:{server.port}/") as resp:
        body = resp.read().decode("utf-8")
    assert resp.status == 200
    assert "ISODISTORT" in body
    assert "langSel" not in body
    assert "generateDistortion" not in body
    assert "showDomains" not in body
    assert 'id="btnDlAll"' in body
    assert "sortMethod1" not in body
    assert "sortResult" in body
    assert "downloadFilteredTable" in body
    assert "downloadFilteredSubgroups" not in body
    assert 'id="dlMethodOpt4"' in body
    assert 'id="btnDlTxt"' in body and 'id="btnDlCsv"' in body
    assert 'id="genDbHelp"' in body
    assert "onResultFilter" in body
    assert "m1Items" in body and "m3Items" in body
    # 对齐官网搜索页布局：well 面板 + OK/Change 按钮 + help 图标
    assert 'class="well"' in body
    assert 'id="hM1"' in body and 'id="hM2"' in body
    assert 'id="btnM1"' in body and 'id="btnTypes"' in body
    assert "help.jpg" in body
    assert "onpagehide" in body or "pagehide" in body or "sendBeacon" in body or "/api/shutdown" in body
    # Distortion 批量导出：官网第 6 页四种格式勾选 + Download all
    assert 'id="fmtCif"' in body and 'id="fmtIsoviz"' in body
    assert 'id="fmtModes"' in body and 'id="fmtTopas"' in body
    assert "downloadAll()" in body
    assert "indices=" in body
    assert 'name="orderparam"' in body
    assert "Save interactive distortion" in body
    assert "Complete modes details" in body
    assert "TOPAS.STR" in body
    assert 'id="dlMethod"' in body
    assert 'id="dlMethodOpt1"' in body and 'id="dlMethodOpt2"' in body
    assert 'id="dlMethodOpt3"' in body
    assert 'id="dlMethodOpt4"' in body
    assert 'multiple' not in body.split('id="dlMethod"')[1].split(">")[0]


def test_i18n_endpoint(server):
    data = _get(server.port, "/api/i18n")
    assert data["ok"]
    assert "language" not in data
    assert "terms" not in data
    assert data["messages"]["load.done"].startswith("[Loaded]")
    assert "ui.menu.language" not in data["messages"]
    assert "Generate isotropy subgroups" in data["messages"]["m2.genDbHelp"]
    assert data["messages"]["m1.orderParam"] == "Order parameter:"
    assert data["messages"]["dist.method4"] == "Method 4"


def test_state_endpoint(server):
    data = _get(server.port, "/api/state")
    assert data["ok"]
    assert data["state"]["structure"] is None
    # 默认对齐官网：strain + displacive（全物种）
    assert data["state"]["distortion_types"] == ["strain", "displacive"]
    assert "distortion_scope" in data["state"]


def test_set_types_with_scope(server):
    """官网 per-species 作用域：类型 + 各类型物种范围一起提交。"""
    body = {
        "types": ["displacive", "occupational"],
        "scope": {
            "displacive": ["*"],
            "occupational": ["Eu"],
            "strain": [],
            "magnetic": [],
            "rotational": [],
        },
    }
    data = _post(server.port, "/api/set_types", body)
    assert data["ok"]
    assert data["state"]["distortion_types"] == ["displacive", "occupational"]
    assert data["state"]["distortion_scope"]["occupational"] == ["Eu"]
    assert data["state"]["distortion_scope"]["displacive"] == ["*"]


def test_ping_endpoint(server):
    data = _get(server.port, "/api/ping")
    assert data["ok"]
    assert web_server._LIFE.page_heartbeat > 0  # 心跳被记录（守护线程停服依据）


def test_shutdown_endpoint(server):
    """/api/shutdown 置位关闭请求（守护线程稍后执行 server.shutdown）。"""
    data = _post(server.port, "/api/shutdown")
    assert data["ok"] and data["shutdown"]
    assert web_server._LIFE.shutdown_requested is True
    web_server._LIFE.shutdown_requested = False  # 复位，避免影响模块内其他测试


def test_space_groups_endpoint(server):
    data = _get(server.port, "/api/space_groups")
    assert data["ok"]
    assert len(data["space_groups"]) == 230
    assert data["space_groups"][0]["number"] == 1
    assert data["space_groups"][0]["symbol"] == "P1"
    assert data["space_groups"][0]["schoenflies"] == "C1-1"
    assert data["space_groups"][-1]["number"] == 230


def test_static_css_served(server):
    with urllib.request.urlopen(f"http://127.0.0.1:{server.port}/static/bootstrap.css") as resp:
        body = resp.read()
    assert resp.status == 200
    assert len(body) > 1000  # 官网 Bootstrap 2.x 样式表已就位
    assert "text/css" in resp.headers["Content-Type"]


def test_download_all_without_method2_is_json_error(server):
    """未完成 Method 2 时 Download all 不得打包 output_dir，应返回 JSON 错误。"""
    url = f"http://127.0.0.1:{server.port}/api/download_all?method=2&formats=cif"
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(url)
    err = excinfo.value
    assert err.code == 404
    data = json.loads(err.read().decode("utf-8"))
    assert data["ok"] is False
    assert "Method 2" in data["error"] or "CIF" in data["error"]


def test_download_all_rejects_multiple_methods(server):
    """method=1,2 必须拒绝，不允许一次导出多个 Method。"""
    url = f"http://127.0.0.1:{server.port}/api/download_all?method=1,2&formats=cif"
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(url)
    err = excinfo.value
    assert err.code == 400
    data = json.loads(err.read().decode("utf-8"))
    assert data["ok"] is False
    assert "不能多选" in data["error"] or "exactly one" in data["error"]


@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_load_cif_endpoint(server):
    """上传 EuAl4 母相 CIF 并校验识别结果（真实计算）。"""
    cif = (Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码")
           / "EuAl4 Parent.cif")
    if not cif.exists():
        cif = (Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码")
               / "EuAl4 Springer (parent).cif")  # 旧命名兼容
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
    # 页头信息（对齐官网 "Reading CIF file..." 段落）：点阵参数 + Wyckoff 坐标
    lat = data["state"]["structure"]["lattice"]
    assert abs(lat["a"] - 4.402) < 0.01
    assert abs(lat["c"] - 11.163) < 0.01
    wy = data["state"]["structure"]["wyckoff"]
    assert len(wy) == 3
    assert any(w["coordinates"] == [0, 0, 0] for w in wy)  # Eu 2a (0,0,0)
    # 官网 "Default space-group preferences: ..." 行
    assert "monoclinic axes a(b)c" in data["state"]["structure"]["preferences"]
    assert data["state"]["species"] == ["Al", "Eu"]


@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_method1_options_endpoint(server):
    """Method 1 下拉数据（真实枚举）：可达子群空间群（<230）+ lattice 选项。"""
    cif = Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码\EuAl4 Parent.cif")
    if not cif.exists():
        pytest.skip("EuAl4 CIF 不存在")
    body = json.dumps({"filename": "eual4.cif", "content": cif.read_text(encoding="utf-8")}).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/load_cif", data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:  # noqa: S310 - 仅访问本地测试服务（http scheme）
        data = json.loads(resp.read().decode("utf-8"))
    assert data["ok"]

    data = _get(server.port, "/api/method1_options")
    assert data["ok"]
    opts = data["options"]
    # 官网只列出与母相结构相容的可达子群，而非全部 230 个
    assert 0 < len(opts["space_groups"]) < 230
    # 序号升序排列（官网下拉按序号排列，便于查找）；母相自身（139）包含在内
    numbers = [sg["number"] for sg in opts["space_groups"]]
    assert numbers == sorted(numbers)
    assert 139 in numbers
    assert opts["conventional_lattices"]
    assert opts["primitive_lattices"]
    # 官网风格标签，如 (1,0,0),(0,1,0),(0,0,1)
    assert "(" in opts["conventional_lattices"][0]["label"]
    # 按格点等价去重后的选项数应远少于旧实现（I4/mmm 旧 19 个，官网 12 个）
    assert len(opts["conventional_lattices"]) < 19
    assert len(opts["primitive_lattices"]) < 19
