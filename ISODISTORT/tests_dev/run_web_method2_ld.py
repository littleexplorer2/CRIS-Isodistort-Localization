"""网页 Method 2 场景回归：精确复现官网交互流程。

场景（与浏览器操作一致）：
  1. 首页上传 EuAl4 Parent.cif（实验数据与GD代码/）
  2. 勾选畸变类型 strains + displacive（作用域 Eu, Al）
  3. Method 2 选择 LD（K10, (0,0,g)），g=1/6，点 OK 枚举子群
     —— 此处曾报错 "'IsoHandler' object has no attribute '_SESSION'"
  4. 选择首个子群计算畸变模式（/api/method2）
  5. 生成畸变结构（/api/generate，对齐官网 Distortion Page）

用法：python tests_dev/run_web_method2_ld.py
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

        # ---- 2. 勾选 strains + displacive(Eu, Al) ---------------------------
        d = _post(port, "/api/set_types", {
            "types": ["strain", "displacive"],
            "scope": {"strain": ["*"], "displacive": ["Eu", "Al"],
                      "occupational": [], "magnetic": [], "rotational": []},
        })
        st = d.get("state", {})
        checks.append(("set_types",
                       d.get("ok") is True
                       and st.get("distortion_types") == ["strain", "displacive"]
                       and st.get("distortion_scope", {}).get("displacive") == ["Eu", "Al"],
                       f"types={st.get('distortion_types')} "
                       f"scope={st.get('distortion_scope', {}).get('displacive')}"))

        # ---- 3. Method 2：LD (K10, (0,0,g)), g=1/6 枚举子群 -----------------
        # 先确认 k 点下拉包含 LD 且带参数 g
        kp = _get(port, "/api/kpoints").get("kpoints", [])
        ld = next((k for k in kp if k["label"] == "LD"), None)
        checks.append(("kpoints_LD",
                       ld is not None and ld.get("parameters") == ["g"],
                       f"LD={ld}"))

        # 3a. 未勾选 generate：不报错即可。若数据库已生成并缓存在 WSL 暂存
        #     目录，则直接返回子群（n>0）；首次使用时数据库缺失返回空列表，
        #     前端据此显示“本地计算 / 前往官网”两个选项（对齐官网交互）。
        #     （旧 bug1：_api_subgroups 误用 self._SESSION 抛 AttributeError；
        #       旧 bug2：KVALUE 在 IR 之前导致 parameters not selected；
        #       旧 bug3：生成应答误用 'y'（应为空行）导致 EOF 崩溃。）
        d = _post(port, "/api/subgroups", {
            "kpoints": [{"k": "LD", "params": ["1/6"]}],
            "generate": False,
        })
        subs = d.get("subgroups", []) or []
        checks.append(("subgroups_no_generate",
                       d.get("ok") is True,
                       f"ok={d.get('ok')} n={len(subs)}(缓存) err={d.get('error')}"))

        # 3b. 勾选 generate：本地 iso 生成 LD(g=1/6) 子群数据库 -> 全部 IR 子群
        d = _post(port, "/api/subgroups", {
            "kpoints": [{"k": "LD", "params": ["1/6"]}],
            "generate": True,
        })
        subs = d.get("subgroups", []) or []
        ir_labels = sorted({s["irrep_label"] for s in subs})
        ld5_opds = sorted({s["opd_symbol"] for s in subs if s["irrep_label"] == "LD5"})
        checks.append(("subgroups_generate", d.get("ok") is True and len(subs) > 0,
                       f"ok={d.get('ok')} n={len(subs)} irreps={ir_labels} "
                       f"LD5_opds={ld5_opds} err={d.get('error')}"))
        checks.append(("subgroups_official_ir_set",
                       ir_labels == ["LD1", "LD2", "LD5"],
                       f"irreps={ir_labels} (expect LD1,LD2,LD5)"))
        checks.append(("subgroups_ld5_opd_count",
                       len(ld5_opds) == 16,
                       f"LD5 n_opds={len(ld5_opds)} (expect 16)"))
        if not subs:
            raise SystemExit(_report(checks))

        # ---- 4. 选首个子群计算模式（/api/method2） --------------------------
        # 已知限制：LD 为参数 k 点，iso 二进制的 DISPLAY BUSH 仅支持对称 k 点，
        # 位移模式计算依赖官网 (3+d) 维超空间机制，本地应返回**明确错误信息**
        # 而非崩溃（畸变后“子群”本身的枚举在步骤 3 已完成）。
        idx = subs[0]["index"]
        d = _post(port, "/api/method2", {
            "subgroup_idx": idx,
            "distortion_type": ["strain", "displacive"],
            "nmod": 0, "nsup": 1,
        })
        modes = d.get("modes", []) or []
        err = d.get("error") or ""
        if d.get("ok"):
            checks.append(("method2_modes", True, f"n_modes={len(modes)}"))
        else:
            # 已知限制：应给出“参数 k 点”说明，而非 AttributeError 等崩溃
            known = ("参数" in err and "k 点" in err) or "k point of symmetry" in err
            checks.append(("method2_modes(已知限制提示)", known,
                           f"ok=False err={err[:160]}"))
        sub0 = subs[0]
        print(f"  子群#{idx}: SG {sub0['space_group_number']} "
              f"{sub0['space_group_symbol']} k={sub0['k_point_label']} "
              f"irrep={sub0['irrep_label']} OPD={sub0['opd_symbol']}",
              flush=True)

        # ---- 5. 生成畸变结构（Distortion Page；有位移模式时） ----------------
        disp = next((m for m in modes if m.get("mode_type") != "occupational"
                     and m.get("n_representatives")), None)
        if disp:
            d = _post(port, "/api/generate",
                      {"contributions": {disp["irrep_label"]: 0.05}})
            checks.append(("generate", d.get("ok") is True,
                           f"ok={d.get('ok')} file={d.get('filename')} "
                           f"err={d.get('error')}"))
    finally:
        httpd.shutdown()
        httpd.server_close()

    return _report(checks)


def _report(checks: list[tuple[str, bool, str]]) -> int:
    print("\n=== Method 2 (LD k10, g=1/6) 场景回归 ===")
    n_ok = 0
    for name, ok, detail in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        n_ok += int(ok)
    print(f"--- {n_ok}/{len(checks)} 通过 ---")
    return 0 if n_ok == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
