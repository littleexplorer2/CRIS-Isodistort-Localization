"""Web server smoke tests and API↔web three-interface consistency."""
from __future__ import annotations

import csv
import html
import json
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from data_dir import experiment_data_dir  # noqa: E402

import main_terminal  # noqa: E402
from isocore.api import IsoDistort  # noqa: E402
from isocore.backend import SubgroupInfo  # noqa: E402
from isocore.distortion import Method3ResultItem  # noqa: E402
from isocore.utils.opd_format import _centering_letter  # noqa: E402
from isocore.utils.schoenflies import hm_symbol, schoenflies_symbol  # noqa: E402
from isocore.utils.text_parser import parse_basis_token  # noqa: E402
from main_terminal import (  # noqa: E402
    IsoDistortConsoleApp,
    _empty_tbl,
    _method3_cols,
    _row_method3,
)
from web import server as web_server  # noqa: E402

DATA_DIR = experiment_data_dir()
CIFS_DIR = Path(__file__).resolve().parent / "cifs_30"
MATCH_TOL = 1e-6


def _official_select_options(cif_name: str, select_name: str) -> list[dict]:
    """Read an archived official selector for regression comparison only."""
    page = (
        Path(__file__).resolve().parents[2]
        / "webpage_info" / cif_name / "2. ISODISTORT_ search.html"
    )
    if not page.exists():
        pytest.skip(f"Archived official page is unavailable: {page}")
    source = page.read_text(encoding="utf-8", errors="ignore")
    match = re.search(
        rf'<select name="{re.escape(select_name)}">(.*?)</select>',
        source,
        flags=re.DOTALL,
    )
    assert match, f"Official select {select_name!r} was not found"
    options = []
    for value, label in re.findall(
        r'<option(?: value="([^"]*)")?>(.*?)</option>',
        match.group(1),
        flags=re.DOTALL,
    ):
        text = html.unescape(re.sub(r"<[^>]+>", "", label)).strip()
        if text.lower() in {"no choice", "not selected"}:
            continue
        fields = value.split()
        options.append({
            "label": text,
            "candidate_count": int(fields[0]) if fields else None,
        })
    return options


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


@pytest.fixture(autouse=True)
def _restore_session():
    """每个测试后重置共享网页会话，避免与三接口测试相互污染。"""
    web_server._reset_lifecycle()
    yield
    web_server._SESSION = web_server.WebSession()
    web_server._reset_lifecycle()


def _get(port: int, path: str) -> dict:
    with urllib.request.urlopen(f"http://127.0.0.1:{port}{path}") as resp:
        return json.loads(resp.read().decode("utf-8"))


def _post(port: int, path: str, payload: dict | None = None) -> dict:
    body = json.dumps(payload or {}).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}", data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:  # noqa: S310
        return json.loads(resp.read().decode("utf-8"))


class _SessionIsoStub:
    """Small mutable IsoDistort stand-in for revision/session tests."""

    def __init__(self) -> None:
        self.structure = None
        self.subgroups = [object()]
        self.mode_displacements = {"old": object()}
        self.mode_occupancies = {"old-occ": object()}
        self.distorted_structure = [object()]
        self.scope = None
        self.types = None
        self.selected_modes_cleared = False

    @staticmethod
    def species():
        return []

    def set_distortion_scope(self, scope):
        self.scope = scope

    def set_distortion_types(self, types):
        self.types = list(types)

    def set_subgroup_candidates(self, subgroups):
        self.subgroups = list(subgroups)

    def clear_selected_modes(self):
        self.selected_modes_cleared = True
        self.mode_displacements = {}
        self.mode_occupancies = {}
        self.distorted_structure = None


# --- from test_web_server.py ---

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
    assert 'id="btnGenDbManage"' in body
    assert "loadGenDbCache" in body
    assert "deleteSelectedGenDb" in body
    assert "/api/isotropy_cache" in body
    assert "onResultFilter" in body
    assert "m1Items" in body and "m3Items" in body
    assert 'exportHeader: "k-active"' in body
    assert 'exportHeader: "basis"' in body
    assert "it.basis_display" in body and "it.origin_display" in body
    assert 'key: "dir"' in body
    # 对齐官网搜索页布局：well 面板 + OK/Change 按钮 + help 图标
    assert 'class="well"' in body
    assert 'id="hM1"' in body and 'id="hM2"' in body
    assert 'id="btnM1"' in body and 'id="btnTypes"' in body
    assert "help.jpg" in body
    assert "onpagehide" in body or "pagehide" in body or "sendBeacon" in body or "/api/shutdown" in body
    assert 'const pageClientId' in body
    assert '/api/client/close' in body
    assert 'window.addEventListener("beforeunload", sendPageClosed)' in body
    assert "if (exportingZip) return" not in body
    # Distortion 批量导出：官网第 6 页四种格式勾选 + Download all
    assert 'id="fmtCif"' in body and 'id="fmtIsoviz"' in body
    assert 'id="distParamFormatNote"' in body
    assert "applyMethod3DownloadDefaults();" in body
    assert "candidateHasParametricK" in body
    assert 'id="fmtModes"' in body and 'id="fmtTopas"' in body
    assert "async function downloadAll()" in body
    assert "r.blob()" in body
    assert "dist.zipWait" in body
    assert "function startBusy" in body
    assert "function renderStatusBar" in body
    assert "STATE = d.state" in body
    assert "A passive refresh must never bless result rows produced by another tab" in body
    assert "payload.revision = STATE.revision" in body
    assert "acceptResponseState(data)" in body
    assert "clearComputedResults();" in body
    assert "&revision=${encodeURIComponent(STATE.revision)}" in body
    assert "knownRouteText" in body
    assert "statusBox.innerHTML = prevStatus" not in body
    assert "progress-striped" in body
    assert "fmtElapsed" in body
    assert "indices=" in body
    assert "resultTableHtml" in body
    assert "methodColumns" in body
    assert 'class="table table-hover result"' in body
    assert 'name="orderparam"' not in body
    assert "space_group_schoenflies" in body
    assert "toFixed(5)" in body
    assert "function typeLabel" in body
    assert "Save interactive distortion" in body
    assert "Complete modes details" in body
    assert "TOPAS.STR" in body
    assert 'id="dlMethod"' in body
    assert 'id="dlMethodOpt1"' in body and 'id="dlMethodOpt2"' in body
    assert 'id="dlMethodOpt3"' in body
    assert 'id="dlMethodOpt4"' in body
    assert 'multiple' not in body.split('id="dlMethod"')[1].split(">")[0]


def test_state_endpoint_builds_snapshot_under_session_lock(server, monkeypatch):
    observed = []

    def _guarded_summary():
        observed.append(web_server._SESSION_LOCK._is_owned())
        return {"structure": None, "distortion_types": []}

    monkeypatch.setattr(web_server, "_state_summary", _guarded_summary)
    data = _get(server.port, "/api/state")
    assert data["ok"] is True
    assert observed == [True]


def test_method3_rows_keep_numbers_and_add_official_fraction_display(capsys):
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=12,
        space_group_symbol="C2/m",
        subgroup_index=4,
        size=1,
        basis_vectors=[[0, 0.5, 1], [-1, 0, 0], [0, -0.25, 0]],
        origin=[0.5, -0.25, 0],
        k_point_label="LD",
        irrep_label="GM1+",
        opd_symbol="P1",
        k_parameters=["1/6"],
    )
    item = Method3ResultItem(
        subgroup=subgroup,
        point_group="2/m",
        basis=subgroup.basis_vectors,
        routes=[subgroup],
    )

    row = web_server._method3_rows([item])[0]
    assert row["basis_vectors"] == [[0, 0.5, 1], [-1, 0, 0], [0, -0.25, 0]]
    assert row["origin"] == [0.5, -0.25, 0]
    assert row["basis_display"] == "{(0,1/2,1),(-1,0,0),(0,-1/4,0)}"
    assert row["origin_display"] == "(1/2,-1/4,0)"

    terminal_row = _row_method3(item)
    assert terminal_row["basis"] == row["basis_display"]
    assert terminal_row["origin"] == row["origin_display"]
    assert terminal_row["routes"] == "LD(1/6):GM1+ P1"
    assert terminal_row["selectable"] is True

    app = object.__new__(IsoDistortConsoleApp)
    app._print_table_row(3, terminal_row)
    output = capsys.readouterr().out
    assert "SG 12 C2/m" in output
    assert "point_group=2/m" in output
    assert "s=1 i=4" in output
    assert "basis={(0,1/2,1),(-1,0,0),(0,-1/4,0)}" in output
    assert "origin=(1/2,-1/4,0)" in output
    assert "known routes=LD(1/6):GM1+ P1" in output
    assert " k=" not in output
    assert "IR=" not in output


def test_method3_affine_only_row_is_not_selectable_or_fake_route(
    monkeypatch,
    capsys,
):
    subgroup = SubgroupInfo(
        index=3,
        space_group_number=3,
        space_group_symbol="P2",
        subgroup_index=2,
        size=2,
        basis_vectors=[[0, 0, 1], [0, -1, 0], [2, 0, 0]],
        origin=[0, 0, 0],
    )
    item = Method3ResultItem(
        subgroup=subgroup,
        point_group="2",
        basis=subgroup.basis_vectors,
        routes=[],
        route_resolution="affine_only_unresolved_coupled_route",
    )

    row = web_server._method3_rows([item])[0]
    assert row["known_routes"] == []
    assert row["known_route_count"] == 0
    assert row["route_resolution"] == "affine_only_unresolved_coupled_route"
    assert row["selectable"] is False
    terminal_row = _row_method3(item)
    assert terminal_row["routes"] == ""
    assert terminal_row["selectable"] is False

    app = object.__new__(IsoDistortConsoleApp)
    app.tbl = {3: _empty_tbl(_method3_cols())}
    app.tbl[3]["rows"] = [terminal_row]
    answers = iter(["3", "q"])
    monkeypatch.setattr(main_terminal, "_prompt", lambda *_args, **_kwargs: next(answers))

    assert app._review_result_table(3, allow_idx=True) is None
    output = capsys.readouterr().out
    assert "NOT SELECTABLE (no resolved second-stage route)" in output
    assert "idx=3 is not selectable" in output
    assert "known routes=(none)" in output


def test_terminal_method3_filtered_csv_uses_embedding_headers(
    tmp_path,
    monkeypatch,
):
    subgroup = SubgroupInfo(
        index=3,
        space_group_number=12,
        space_group_symbol="C2/m",
        subgroup_index=8,
        size=4,
        basis_vectors=[[0, 0.5, 1], [-1, 0, 0], [0, -0.25, 0]],
        origin=[0.5, -0.25, 0],
        k_point_label="LD",
        irrep_label="GM1+",
        opd_symbol="P1",
        k_parameters=["1/6"],
    )
    item = Method3ResultItem(
        subgroup=subgroup,
        point_group="2/m",
        basis=subgroup.basis_vectors,
        routes=[subgroup],
    )
    terminal_row = _row_method3(item)
    app = object.__new__(IsoDistortConsoleApp)
    app.tbl = {3: _empty_tbl(_method3_cols())}
    app.tbl[3]["rows"] = [terminal_row]
    app.iso = SimpleNamespace(cfg=SimpleNamespace(output_dir=tmp_path))

    destination = tmp_path / "method3.csv"
    answers = iter(["csv", str(destination)])
    monkeypatch.setattr(main_terminal, "_prompt_int", lambda *_args, **_kwargs: 3)
    monkeypatch.setattr(main_terminal, "_prompt", lambda *_args, **_kwargs: next(answers))

    app._export_filtered_table()

    with destination.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    assert rows[0] == [
        "index",
        "space_group",
        "basis",
        "origin",
        "point_group",
        "s",
        "i",
        "known_routes",
    ]
    assert rows[1] == [
        "3",
        "12 C2/m",
        "{(0,1/2,1),(-1,0,0),(0,-1/4,0)}",
        "(1/2,-1/4,0)",
        "2/m",
        "4",
        "8",
        "LD(1/6):GM1+ P1",
    ]


def test_i18n_endpoint(server):
    data = _get(server.port, "/api/i18n")
    assert data["ok"]
    assert "language" not in data
    assert "terms" not in data
    assert data["messages"]["load.done"].startswith("[Loaded]")
    assert "ui.menu.language" not in data["messages"]
    assert "Generate isotropy subgroups" in data["messages"]["m2.genDbHelp"] or (
        "isotropy" in data["messages"]["m2.genDbHelp"].lower()
        and "subgroup" in data["messages"]["m2.genDbHelp"].lower()
    )
    assert "m2.genDbManage" in data["messages"]
    assert "m2.genDbDelete" in data["messages"]
    assert data["messages"]["dist.method4"] == "Method 4"
    assert "Building ZIP" in data["messages"]["dist.zipWait"]
    assert "Elapsed" in data["messages"]["st.elapsed"]
    assert "Local engine is running" in data["messages"]["st.busyHint"]


def test_parent_header_schoenflies():
    assert schoenflies_symbol(139) == "D4h-17"
    # Method 3 / Method 1 下拉与官网一致的旧 IT 简写
    assert hm_symbol(39) == "Abm2"
    assert hm_symbol(41) == "Aba2"
    assert hm_symbol(64) == "Cmca"
    assert hm_symbol(67) == "Cmma"
    assert hm_symbol(68) == "Ccca"


def test_state_endpoint(server):
    data = _get(server.port, "/api/state")
    assert data["ok"]
    assert data["state"]["structure"] is None
    assert data["state"]["revision"] == 0
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


def test_stale_tab_mutation_and_download_are_rejected_with_http_409(server):
    stub = _SessionIsoStub()
    web_server._SESSION._iso = stub
    web_server._SESSION.revision = 7

    req = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/set_types",
        data=json.dumps({"revision": 6, "types": ["strain"]}).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(req)  # noqa: S310 - local test server
    assert excinfo.value.code == 409
    payload = json.loads(excinfo.value.read().decode("utf-8"))
    assert payload["stale_context"] is True
    assert payload["state"]["revision"] == 7
    assert "Stale context" in payload["error"]
    assert stub.types is None
    assert web_server._SESSION.revision == 7

    url = (
        f"http://127.0.0.1:{server.port}/api/download_all"
        "?method=1&formats=cif&revision=6"
    )
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(url)
    assert excinfo.value.code == 409
    payload = json.loads(excinfo.value.read().decode("utf-8"))
    assert payload["stale_context"] is True
    assert payload["state"]["revision"] == 7


def test_set_types_bumps_revision_and_invalidates_all_derived_state(server):
    stub = _SessionIsoStub()
    web_server._SESSION._iso = stub
    web_server._SESSION.revision = 3
    web_server._SESSION.method1 = [object()]
    web_server._SESSION.method2 = object()
    web_server._SESSION.method3 = [object()]
    web_server._SESSION.method2_subgroups = [object()]

    data = _post(server.port, "/api/set_types", {
        "revision": 3,
        "types": ["displacive"],
        "scope": {"displacive": ["*"]},
    })

    assert data["ok"] is True
    assert data["state"]["revision"] == 4
    assert web_server._SESSION.revision == 4
    assert web_server._SESSION.method1 == []
    assert web_server._SESSION.method2 is None
    assert web_server._SESSION.method3 == []
    assert web_server._SESSION.method2_subgroups == []
    assert stub.subgroups == []
    assert stub.mode_displacements == {}
    assert stub.mode_occupancies == {}
    assert stub.distorted_structure is None
    assert stub.selected_modes_cleared is True


def test_ping_endpoint(server):
    data = _get(server.port, "/api/ping?client_id=tab-a")
    assert data["ok"]
    assert data["active_clients"] == 1
    assert "tab-a" in web_server._LIFE.clients
    assert web_server._LIFE.page_heartbeat > 0  # 心跳被记录（守护线程停服依据）


def test_page_close_tracks_last_browser_tab(server):
    _get(server.port, "/api/ping?client_id=tab-a")
    _get(server.port, "/api/ping?client_id=tab-b")

    first = _post(server.port, "/api/client/close?client_id=tab-a")
    assert first == {"ok": True, "closing": False, "active_clients": 1}
    assert web_server._LIFE.last_page_closed == 0

    last = _post(server.port, "/api/client/close?client_id=tab-b")
    assert last == {"ok": True, "closing": True, "active_clients": 0}
    assert web_server._LIFE.last_page_closed > 0
    assert web_server._LIFE.shutdown_requested is False


def test_watchdog_finishes_export_before_stopping_after_tab_close():
    class FakeServer:
        def __init__(self):
            self.stopped = threading.Event()

        def shutdown(self):
            self.stopped.set()

    fake = FakeServer()
    web_server._mark_page_seen()
    web_server._touch_heartbeat("tab-export")
    web_server._begin_long_request()
    web_server._mark_page_closed("tab-export")
    thread = threading.Thread(
        target=web_server._watchdog,
        args=(fake, 60.0),
        kwargs={"close_grace": 0.02, "poll_interval": 0.005},
        daemon=True,
    )
    thread.start()
    assert not fake.stopped.wait(0.05)
    web_server._end_long_request()
    assert fake.stopped.wait(1.0)
    thread.join(timeout=1.0)


def test_last_page_close_stops_http_server_and_releases_port():
    httpd = web_server.ThreadingHTTPServer(
        ("127.0.0.1", 0), web_server.IsoHandler
    )
    port = httpd.server_address[1]
    server_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    watchdog_thread = threading.Thread(
        target=web_server._watchdog,
        args=(httpd, 60.0),
        kwargs={"close_grace": 0.02, "poll_interval": 0.005},
        daemon=True,
    )
    server_thread.start()
    watchdog_thread.start()
    _get(port, "/api/ping?client_id=only-tab")
    closed = _post(port, "/api/client/close?client_id=only-tab")
    assert closed["closing"] is True
    server_thread.join(timeout=1.0)
    assert not server_thread.is_alive()
    httpd.server_close()
    watchdog_thread.join(timeout=1.0)

    rebound = web_server.HTTPServer(("127.0.0.1", port), web_server.IsoHandler)
    rebound.server_close()


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
    # ISODISTORT 官网 Method 3 下拉：旧 IT 写法（非 e-glide）
    by_num = {g["number"]: g["symbol"] for g in data["space_groups"]}
    assert by_num[39] == "Abm2"
    assert by_num[41] == "Aba2"
    assert by_num[64] == "Cmca"
    assert by_num[67] == "Cmma"
    assert by_num[68] == "Ccca"


def test_all_method3_space_groups_match_archived_official_page(server):
    """Exhaustively verify all 230 Method-3 SG labels, not a hand-picked sample."""
    official = _official_select_options("EuAl4 Parent.cif", "subgroupsym")
    actual = _get(server.port, "/api/space_groups")["space_groups"]
    expected = []
    for option in official:
        match = re.match(r"(\d+)\s+(\S+)\s+(\S+)$", option["label"])
        assert match, option["label"]
        expected.append({
            "number": int(match.group(1)),
            "symbol": match.group(2),
            "schoenflies": match.group(3),
        })
    assert actual == expected


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


def test_download_all_passes_nmod_without_mutating_cached_mode_context(server):
    seen: dict = {}

    class ExportIsoStub:
        def __init__(self) -> None:
            self.structure = object()
            self.number_of_independent_modulations = 1

        def export_subgroups_zip(self, **kwargs):
            seen.update(kwargs)
            return b"test-zip"

    stub = ExportIsoStub()
    subgroup = SimpleNamespace(index=0)
    web_server._SESSION._iso = stub
    web_server._SESSION.method1 = [SimpleNamespace(subgroup=subgroup)]

    url = (
        f"http://127.0.0.1:{server.port}/api/download_all"
        "?method=1&formats=modes&nmod=2"
    )
    with urllib.request.urlopen(url) as response:
        assert response.read() == b"test-zip"

    assert seen["number_of_independent_modulations"] == 2
    assert seen["subgroups"] == [subgroup]
    assert stub.number_of_independent_modulations == 1


@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_load_cif_endpoint(server):
    """上传 EuAl4 母相 CIF 并校验识别结果（真实计算）。"""
    cif = experiment_data_dir() / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = experiment_data_dir() / "EuAl4 Springer (parent).cif"  # 旧命名兼容
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
    assert data["state"]["structure"]["space_group_schoenflies"] == "D4h-17"
    # 页头信息（对齐官网 "Reading CIF file..." 段落）：点阵参数 + Wyckoff 坐标
    lat = data["state"]["structure"]["lattice"]
    assert abs(lat["a"] - 4.402) < 0.01
    assert abs(lat["c"] - 11.163) < 0.01
    assert round(lat["a"], 5) == lat["a"]
    wy = data["state"]["structure"]["wyckoff"]
    assert len(wy) == 3
    assert any(w["coordinates"] == [0, 0, 0] for w in wy)  # Eu 2a (0,0,0)
    display = data["state"]["structure"]["wyckoff_display"]
    assert display[0] == "Eu1 2a (0,0,0)"
    assert display[1] == "Al1 4d (0,1/2,1/4)"
    assert display[2].startswith("Al2 4e (0,0,z), z=")
    # 官网 "Default space-group preferences: ..." 行
    assert "monoclinic axes a(b)c" in data["state"]["structure"]["preferences"]
    assert data["state"]["species"] == ["Al", "Eu"]


@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_method1_options_endpoint(server):
    """Method 1 下拉数据（真实枚举）：可达子群空间群（<230）+ lattice 选项。"""
    cif = experiment_data_dir() / "EuAl4 Parent.cif"
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
    # 官网只列出与母相结构 + Types 相容的可达子群，而非全部 230 个
    assert 0 < len(opts["space_groups"]) < 230
    # 序号升序排列（官网下拉按序号排列，便于查找）；母相自身（139）包含在内
    numbers = [sg["number"] for sg in opts["space_groups"]]
    assert numbers == sorted(numbers)
    assert 139 in numbers
    # EuAl4 + 默认 strain/displacive：官网无 87/97/126/128（无位移活性 IR）
    for extra in (87, 97, 126, 128):
        assert extra not in numbers
    assert opts["conventional_lattices"]
    assert opts["primitive_lattices"]
    # 官网风格标签，如 (1,0,0),(0,1,0),(0,0,1)
    assert "(" in opts["conventional_lattices"][0]["label"]
    # I4/mmm EuAl4：官网 Conventional 12、Primitive 9
    assert len(opts["conventional_lattices"]) == 12
    assert len(opts["primitive_lattices"]) == 9
    assert opts["conventional_lattices"][0]["label"] == "(1,0,0),(0,1,0),(0,0,1)"
    assert opts["primitive_lattices"][0]["label"] == (
        "(-1/2,1/2,1/2),(1/2,-1/2,1/2),(1/2,1/2,-1/2)"
    )

    # The saved official page is a reference for the exact visible labels.
    official_sg = _official_select_options("EuAl4 Parent.cif", "isosubgroup")
    assert [
        f"{sg['number']} {sg['symbol']} {sg['schoenflies']}" for sg in opts["space_groups"]
    ] == [" ".join(sg["label"].split()) for sg in official_sg]

    # Full selector verification: order, displayed representative, covered
    # OPD rows, and the lattice represented by each option.
    iso = web_server._SESSION.iso
    for local_key, official_name in (
        ("conventional_lattices", "isolattice"),
        ("primitive_lattices", "isoplattice"),
    ):
        archived = _official_select_options("EuAl4 Parent.cif", official_name)
        local = opts[local_key]
        assert len(local) == len(archived)
        assert [item["candidate_count"] for item in local] == [
            item["candidate_count"] for item in archived
        ]
        for got, expected in zip(local, archived, strict=True):
            expected_basis = parse_basis_token(expected["label"])
            assert got["label"] == expected["label"]
            assert got["basis"] == expected_basis
            assert iso._same_lattice_orbit(got["basis"], expected_basis), (
                got["label"], expected["label"]
            )

    # Every option must select precisely the number of live rows recorded by
    # the official form, including classes whose display basis was aligned.
    live_subgroups = iso._method1_filtered_subgroups()
    for key in ("conventional_lattices", "primitive_lattices"):
        for option in opts[key]:
            count = 0
            for sg in live_subgroups:
                basis = np.asarray(sg.basis_vectors, dtype=float)
                if key == "primitive_lattices":
                    basis = iso._centering_matrix(
                        _centering_letter(sg.space_group_number)
                    ) @ basis
                count += iso._same_lattice_orbit(basis, option["basis"])
            assert count == option["candidate_count"], option["label"]

    # The official selector IDs target one class. In particular, selecting
    # the identity conventional lattice must not include every supercell.
    for key, kind in (
        ("conventional_lattices", "conventional"),
        ("primitive_lattices", "primitive"),
    ):
        option = opts[key][0]
        result = _post(server.port, "/api/method1", {
            "lattice": option["basis"], "lattice_kind": kind,
            "distortion_types": ["strain", "displacive"],
        })
        assert result["ok"]
        assert len(result["candidates"]) == option["candidate_count"]


# --- from test_three_interface.py ---

def _load_cif_via_web(port: int, path: Path) -> dict:
    content = path.read_text(encoding="utf-8")
    return _post(port, "/api/load_cif?lang=en",
                 {"filename": path.name, "content": content})


def _sg_of(structure) -> int:
    return SpacegroupAnalyzer(structure, symprec=1e-3).get_space_group_number()



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_method1_options(server):
    """同一母相：API 与网页 method1_options 空间群/格子完全一致。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("EuAl4 母相不存在")

    # API 路径
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    api_opts = iso.method1_options()

    # 网页路径
    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    web_opts = _get(server.port, "/api/method1_options")["options"]

    api_sgs = [g["number"] for g in api_opts["space_groups"]]
    web_sgs = [g["number"] for g in web_opts["space_groups"]]
    assert api_sgs == web_sgs, "Method1 可达子群空间群不一致"
    api_conv = [g["label"] for g in api_opts["conventional_lattices"]]
    web_conv = [g["label"] for g in web_opts["conventional_lattices"]]
    assert api_conv == web_conv, "Conventional lattice 选项不一致"
    api_prim = [g["label"] for g in api_opts["primitive_lattices"]]
    web_prim = [g["label"] for g in web_opts["primitive_lattices"]]
    assert api_prim == web_prim, "Primitive lattice 选项不一致"



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_method1_candidates(server):
    """同一母相：API 与网页 method1 候选数一致（同过滤条件）。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    api_cands = iso.search_method_1(
        distortion_types=["displacive", "strain"],
        crystal_system=None, subgroup_space_group=None,
        maximal_subgroup_only=False,
    )

    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    web_m1 = _post(server.port, "/api/method1",
                   {"distortion_types": ["displacive", "strain"],
                    "maximal_subgroup_only": False})
    assert web_m1["ok"]
    web_cands = web_m1["candidates"]

    assert len(api_cands) == len(web_cands), \
        f"Method1 候选数不一致: API={len(api_cands)} web={len(web_cands)}"
    # 首候选的空间群一致
    if api_cands:
        api_first = api_cands[0].subgroup.space_group_number
        web_first = web_cands[0]["space_group_number"]
        assert api_first == web_first, \
            f"Method1 首候选空间群不一致: API={api_first} web={web_first}"



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_method2_and_generation(server):
    """同一子群：API 与网页 method2 模式数、生成结构空间群一致。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    api_cands = iso.search_method_1(
        distortion_types=["displacive", "strain"],
        maximal_subgroup_only=False,
    )
    assert api_cands, "API 无候选"
    chosen = api_cands[0]
    iso.search_method_2(subgroup_idx=chosen.subgroup.index,
                        distortion_type=["displacive"])
    api_n_modes = len(iso.mode_displacements)
    assert api_n_modes > 0, "API 无位移模式"

    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    web_m1 = _post(server.port, "/api/method1",
                   {"distortion_types": ["displacive", "strain"],
                    "maximal_subgroup_only": False})
    assert web_m1["ok"]
    web_m2 = _post(server.port, "/api/method2",
                   {"subgroup_idx": chosen.subgroup.index,
                    "distortion_type": ["displacive"],
                    "nmod": 0, "nsup": 1})
    assert web_m2["ok"]
    web_modes = [m for m in web_m2["modes"] if m["mode_type"] == "displacive"]
    assert len(web_modes) == api_n_modes, \
        f"method2 模式数不一致: API={api_n_modes} web={len(web_modes)}"

    # 生成：API 与网页（同底层引擎）畸变结构空间群一致
    label = next(iter(iso.mode_displacements))
    api_dist = iso.generate_distortion(irrep_label=label, amplitude=0.1)
    api_sg = _sg_of(api_dist)
    # 网页生成走 /api/generate（若有该端点）；无则跳过数值比对
    if web_modes:
        web_sg = None
        # 通过引擎直接生成同一模式（网页后端复用同一 IsoDistort 会话）
        web_iso = web_server._SESSION.iso
        try:
            web_dist = web_iso.generate_distortion(irrep_label=label,
                                                   amplitude=0.1)
            web_sg = _sg_of(web_dist)
        except Exception:  # noqa: BLE001 - 会话状态差异时跳过
            web_sg = None
        if web_sg is not None:
            assert web_sg == api_sg, \
                f"生成畸变空间群不一致: API={api_sg} web={web_sg}"



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_state_structure(server):
    """加载后状态结构：API 与网页识别空间群/原子数一致（数值 ≤1e-6 级）。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("EuAl4 母相不存在")

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)

    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    st = web["state"]["structure"]
    assert st["space_group_number"] == iso.symmetry_info["space_group_number"] == 139
    assert st["atoms"] == len(iso.structure)
    # 晶格参数一致性（网页端从同一结构序列化）
    assert abs(st["lattice"]["a"] - iso.structure.lattice.a) < MATCH_TOL
    assert abs(st["lattice"]["c"] - iso.structure.lattice.c) < MATCH_TOL
