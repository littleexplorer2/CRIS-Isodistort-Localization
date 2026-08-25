"""网页 HTTP 端点抽查：在进程内启动 web 服务，对 3 个母相 CIF 走
与浏览器完全相同的 HTTP 流程（load_cif -> method1_options -> method2 ->
generate），校验状态/响应结构与关键数值。

抽查对象：
- EuAl4 Parent.cif（真实数据，I4/mmm #139）
- sg123.cif（生成数据，P4/mmm）
- sg225.cif（生成数据，Fm-3m）

用法：python tests_dev/run_web_spotcheck.py
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

from data_dir import experiment_data_dir

DATA_DIR = experiment_data_dir()
CIFS_DIR = ROOT / "tests_dev" / "cifs_30"

TARGETS = [
    ("EuAl4_Parent", DATA_DIR / "EuAl4 Parent.cif", 139),
    ("sg123", CIFS_DIR / "sg123.cif", 123),
    ("sg225", CIFS_DIR / "sg225.cif", 225),
]


def _get(port: int, path: str) -> dict:
    with urllib.request.urlopen(f"http://127.0.0.1:{port}{path}") as resp:
        return json.loads(resp.read().decode("utf-8"))


def _post(port: int, path: str, payload: dict) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}", data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:  # noqa: S310 - 仅访问本地测试服务（http scheme）
        return json.loads(resp.read().decode("utf-8"))


def main() -> int:
    httpd = web_server.HTTPServer(("127.0.0.1", 0), web_server.IsoHandler)
    port = httpd.server_address[1]
    th = threading.Thread(target=httpd.serve_forever, daemon=True)
    th.start()
    time.sleep(0.3)

    results = []
    try:
        for name, path, expected_sg in TARGETS:
            rec: dict = {"name": name, "expected_sg": expected_sg}
            content = path.read_text(encoding="utf-8")
            data = _post(port, "/api/load_cif?lang=en",
                         {"filename": f"{name}.cif", "content": content})
            rec["load_ok"] = data.get("ok")
            st = data.get("state", {})
            rec["sg"] = (st.get("structure") or {}).get("space_group_number")
            rec["sg_ok"] = rec["sg"] == expected_sg

            opts = _get(port, "/api/method1_options")
            o = opts.get("options", {})
            rec["reachable"] = len(o.get("space_groups", []))
            nums = [g["number"] for g in o.get("space_groups", [])]
            rec["sorted_ok"] = nums == sorted(nums)
            rec["n_conv"] = len(o.get("conventional_lattices", []))
            rec["n_prim"] = len(o.get("primitive_lattices", []))

            # method1 搜索（无过滤），取首个候选
            m1 = _post(port, "/api/method1",
                       {"distortion_types": ["displacive"],
                        "maximal_subgroup_only": False})
            cands = m1.get("candidates", [])
            rec["m1_candidates"] = len(cands)

            # method2 走一次（取首个候选；不要求有模式，只要求端点不报错）
            m2 = None
            if cands:
                m2 = _post(port, "/api/method2",
                           {"subgroup_idx": cands[0]["index"],
                            "distortion_type": ["displacive"],
                            "scope": {"displacive": ["*"], "occupational": [],
                                      "strain": [], "magnetic": [],
                                      "rotational": []},
                            "nmod": 0, "nsup": 1})
            rec["method2_ok"] = m2 is not None and m2.get("ok") is True
            rec["method2_modes"] = len((m2 or {}).get("modes", []) or [])

            # 全链通过判定
            rec["pass"] = (rec["load_ok"] and rec["sg_ok"] and rec["sorted_ok"]
                           and rec["method2_ok"])
            results.append(rec)
            tag = "PASS" if rec["pass"] else "FAIL"
            print(f"[{tag}] {name}: sg={rec['sg']}(exp {expected_sg}) "
                  f"reachable={rec['reachable']} conv={rec['n_conv']} "
                  f"prim={rec['n_prim']} m2_modes={rec['method2_modes']}",
                  flush=True)
    finally:
        httpd.shutdown()
        httpd.server_close()

    n_pass = sum(1 for r in results if r.get("pass"))
    print(f"\n=== 网页抽查汇总: {n_pass}/{len(results)} 通过 ===")
    out = ROOT / "output" / "web_spotcheck_report.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print(f"  报告: {out}")
    return 0 if n_pass == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
