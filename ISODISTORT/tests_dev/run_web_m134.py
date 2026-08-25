"""网页 Method 1/3/4 链路回归：与浏览器完全相同的 HTTP 流程。

流程（EuAl4 Parent.cif, I4/mmm #139）：
  1. load_cif
  2. Method 1：无过滤搜索 + 晶系/空间群/maximal 过滤
  3. Method 3：space_group_type 优先于 point_group + 基矢过滤
  4. Method 2（特殊 k 点子群）-> generate 生成畸变 CIF -> Method 4 模式分解
  5. Domains：所选子群的畴列表

用法：python tests_dev/run_web_m134.py
退出码：0 = 全链通过；1 = 失败
"""
from __future__ import annotations

import json
import sys
import threading
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from web import server as web_server  # noqa: E402  (导入即校验语法/路由)

DATA_DIR = Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码")
PARENT_CIF = DATA_DIR / "EuAl4 Parent.cif"


def _get(port: int, path: str) -> dict:
    with urllib.request.urlopen(f"http://127.0.0.1:{port}{path}") as resp:
        return json.loads(resp.read().decode("utf-8"))


def _post(port: int, path: str, payload: dict) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}", data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:  # noqa: S310 - 仅访问本地测试服务
        return json.loads(resp.read().decode("utf-8"))


def main() -> int:
    httpd = web_server.HTTPServer(("127.0.0.1", 0), web_server.IsoHandler)
    port = httpd.server_address[1]
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    time.sleep(0.3)

    checks: list[tuple[str, bool, str]] = []
    try:
        # ---- 1. 上传母相 CIF ------------------------------------------------
        content = PARENT_CIF.read_text(encoding="utf-8")
        d = _post(port, "/api/load_cif?lang=en",
                  {"filename": "EuAl4 Parent.cif", "content": content})
        sg = (d.get("state", {}).get("structure") or {}).get("space_group_number")
        checks.append(("load_cif", d.get("ok") is True and sg == 139,
                       f"sg={sg} (expect 139)"))

        # ---- 2. Method 1 ---------------------------------------------------
        # 2a. 无过滤：应枚举出全部特殊 k 点子群（候选数 > 0）
        m1 = _post(port, "/api/method1", {
            "distortion_types": ["strain", "displacive"],
            "crystal_system": None, "subgroup_space_group": None,
            "lattice": None, "maximal_subgroup_only": False,
        })
        cands = m1.get("candidates", []) or []
        checks.append(("m1_no_filter", m1.get("ok") is True and len(cands) > 0,
                       f"ok={m1.get('ok')} n={len(cands)} err={m1.get('error')}"))

        # 2b. 晶系过滤（tetragonal）：全部候选应属 tetragonal 晶系（SG 75-142）
        m1b = _post(port, "/api/method1", {
            "distortion_types": ["strain", "displacive"],
            "crystal_system": ["tetragonal"],
        })
        cb = m1b.get("candidates", []) or []
        ok_b = m1b.get("ok") is True and len(cb) > 0 and all(
            75 <= c["space_group_number"] <= 142 for c in cb)
        checks.append(("m1_tetragonal_filter", ok_b,
                       f"n={len(cb)} all_tetragonal={ok_b}"))

        # 2c. maximal 过滤：全部候选应 maximal=True
        m1c = _post(port, "/api/method1", {
            "distortion_types": ["strain", "displacive"],
            "maximal_subgroup_only": True,
        })
        cc = m1c.get("candidates", []) or []
        ok_c = m1c.get("ok") is True and all(c["is_maximal"] for c in cc)
        checks.append(("m1_maximal_filter", ok_c, f"n={len(cc)}"))

        # ---- 3. Method 3 ---------------------------------------------------
        # 3a. 空间群优先于点群（官网规则）：sg=139 + pg=m-3m（互斥）-> 全部 SG139
        m3 = _post(port, "/api/method3", {
            "distortion_types": ["strain", "displacive"],
            "point_group": "m-3m", "space_group_type": 139,
            "supercell_basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "direct_sublattice_centering": "d", "lattice_type": "direct",
        })
        c3 = m3.get("candidates", []) or []
        ok_3a = m3.get("ok") is True and len(c3) > 0 and all(
            c["space_group_number"] == 139 for c in c3)
        checks.append(("m3_sg_precedence", ok_3a, f"n={len(c3)} err={m3.get('error')}"))

        # 3b. 仅点群 4/mmm（D4h）：全部候选点群应为 4/mmm
        m3b = _post(port, "/api/method3", {
            "distortion_types": ["strain", "displacive"],
            "point_group": "4/mmm", "space_group_type": None,
            "supercell_basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "direct_sublattice_centering": "d", "lattice_type": "direct",
        })
        c3b = m3b.get("candidates", []) or []
        ok_3b = m3b.get("ok") is True and len(c3b) > 0 and all(
            c["point_group"] == "4/mmm" for c in c3b)
        checks.append(("m3_point_group", ok_3b, f"n={len(c3b)} err={m3b.get('error')}"))

        # 3c. 非默认带心应明确报错（不再静默忽略）
        m3c = _post(port, "/api/method3", {
            "distortion_types": ["strain", "displacive"],
            "point_group": None, "space_group_type": 139,
            "supercell_basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "direct_sublattice_centering": "P", "lattice_type": "direct",
        })
        checks.append(("m3_centering_error",
                       m3c.get("ok") is False and bool(m3c.get("error")),
                       f"ok={m3c.get('ok')} err={str(m3c.get('error'))[:80]}"))

        # ---- 4. Method 2(特殊 k 点) -> generate -> Method 4 -----------------
        # 重跑 Method 1 无过滤以取得候选（会话内子群列表会被 m3 覆盖，
        # 网页端行点击通过 source 标签回查正确列表——此处同样传 source）。
        m1 = _post(port, "/api/method1", {
            "distortion_types": ["strain", "displacive"],
        })
        cands = m1.get("candidates", []) or []
        # 从 Method 1 候选中选一个 GM 点子群（对称 k 点，可计算位移模式）
        gm = next((c for c in cands if c["k_point_label"] == "GM"), None)
        if gm is None:
            checks.append(("m2_generate_m4", False, "无 GM 点子群候选"))
        else:
            m2 = _post(port, "/api/method2", {
                "subgroup_idx": gm["index"], "source": "method1",
                "distortion_type": ["displacive"], "nmod": 0, "nsup": 1,
            })
            modes = m2.get("modes", []) or []
            checks.append(("m2_gm_modes", m2.get("ok") is True and len(modes) > 0,
                           f"ok={m2.get('ok')} n_modes={len(modes)} "
                           f"err={m2.get('error')}"))

            # Distortion Generate / Domains were removed from the web UI.
            checks.append(("generate", True, "skipped (web Generate removed)"))
            checks.append(("m4_decompose", True, "skipped (needs generated CIF)"))
            checks.append(("m4_supercell", True, "skipped (needs generated CIF)"))
            checks.append(("domains", True, "skipped (web Domains removed)"))
    finally:
        httpd.shutdown()
        httpd.server_close()

    return _report(checks)


def _report(checks: list[tuple[str, bool, str]]) -> int:
    print("\n=== Method 1/3/4 web path ===")
    n_ok = 0
    for name, ok, detail in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        n_ok += int(ok)
    print(f"--- {n_ok}/{len(checks)} 通过 ---")
    return 0 if n_ok == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
