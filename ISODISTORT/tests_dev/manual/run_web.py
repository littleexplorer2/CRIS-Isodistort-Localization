"""网页 HTTP 手工回归。

用法：
  python tests_dev/manual/run_web.py spotcheck
  python tests_dev/manual/run_web.py m134
  python tests_dev/manual/run_web.py method2_ld
"""
from __future__ import annotations

import sys
from pathlib import Path

_TESTS_DEV = Path(__file__).resolve().parents[1]
ROOT = _TESTS_DEV.parent  # ISODISTORT
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(_TESTS_DEV))


import json
import threading
import time
import urllib.request


from web import server as web_server  # noqa: E402  (导入即校验语法/路由)

from data_dir import experiment_data_dir

DATA_DIR = experiment_data_dir()
PARENT_CIF = DATA_DIR / "EuAl4 Parent.cif"
CIFS_DIR = _TESTS_DEV / "cifs_30"

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


def cmd_spotcheck() -> int:
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




def cmd_m134() -> int:
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




def cmd_method2_ld() -> int:
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

        # Distortion Generate was removed from the web UI (parent CIF → subgroup files).
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





def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description="Web HTTP manual runners")
    p.add_argument(
        "command",
        choices=["spotcheck", "m134", "method2_ld"],
    )
    args = p.parse_args(argv)
    return {
        "spotcheck": cmd_spotcheck,
        "m134": cmd_m134,
        "method2_ld": cmd_method2_ld,
    }[args.command]()


if __name__ == "__main__":
    raise SystemExit(main())
