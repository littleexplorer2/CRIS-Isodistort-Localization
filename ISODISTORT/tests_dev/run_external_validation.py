"""外部结构源验证：COD 真实 CIF + 官网示例路径对照。

背景：用外部数据库（COD 真实晶体结构）驱动本地「网页版交互」同底层 API，
与 Isodistort 官网已知输出做差异对比，检查网页交互是否有 bug。
（数据库内容仅作为测试输入，不嵌入程序。）

三层验证：
1. API 全流程（与网页同一底层）：加载 → SG 识别 vs COD 记录 →
   method1_options 可达子群 → 前若干子群 method2 模式 → 生成畸变 →
   spglib 校验 SG == 目标子群 → 物理自洽自检。
2. 官网输出对照：SrTiO₃ 三条文献公认路径（R₄⁺→I4/mcm #140、
   M₃⁺→P4/mbm #127、Γ₄⁻→P4mm #99，来源 isodistortexample.php 官网示例）；
   其余结构做自洽校验（官网无公开逐条输出时以母相识别 + 子群自洽为准）。
3. 网页 HTTP 抽查：每个结构走与浏览器相同的 HTTP 流程
   （load_cif → method1_options → method1 → method2 → generate），
   检查网页端点无报错、响应与 API 结果一致。

用法：python tests_dev/run_external_validation.py [--only NAME]
"""
from __future__ import annotations

import argparse
import json
import sys
import threading
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: E402
from run_30cif_validation import _validate_parent  # noqa: E402

from isocore.api import IsoDistort  # noqa: E402
from web import server as web_server  # noqa: E402

CIFS = Path(__file__).resolve().parent / "cifs_external"

# name -> 预期空间群集合（与 fetch_cod_cifs.TARGETS 一致）
EXPECTED: dict[str, set[int]] = {
    "SrTiO3_Pm3m": {221}, "BaTiO3_Pm3m": {221}, "WO3_Pm3m": {221},
    "LaMnO3_Pbnm": {62}, "Cu_Fm3m": {225}, "Si_Fd3m": {227},
    "NaCl_Fm3m": {225}, "BaTiO3_P4mmm": {123, 99}, "BiFeO3_R3c": {161},
    "ZnO_P63mc": {186}, "GaN_P63mc": {186}, "CsCl_Pm3m": {221},
    "SiO2_quartz": {152, 154}, "MgO_Fm3m": {225},
    # 第二批（不同 COD 批次）
    "Fe2O3_hematite": {167}, "CaCO3_calcite": {167}, "TiO2_rutile": {136},
    "TiO2_anatase": {141}, "CaF2_fluorite": {225}, "BaSO4_barite": {62},
    "MgAl2O4_spinel": {227}, "KBr_Fm3m": {225}, "Al2O3_corundum": {167},
    "CaTiO3_Pbnm": {62},
}

# 官网示例已知路径（isodistortexample.php）：(k, IR, 目标子群)
OFFICIAL_SRTIO3_PATHS = [("R", "R4+", 140), ("M", "M3+", 127),
                         ("GM", "GM4-", 99)]


def _sg_of(structure) -> int:
    return SpacegroupAnalyzer(structure, symprec=1e-3).get_space_group_number()


def _official_paths_check(path: Path) -> list[dict]:
    """官网示例路径对照：SrTiO₃ 的三条已知子群路径。"""
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(path)
    recs: list[dict] = []
    for k, ir, expected in OFFICIAL_SRTIO3_PATHS:
        rec: dict = {"k": k, "ir": ir, "official_target": expected}
        try:
            subs = iso.list_subgroups_at(k, ir)
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"list_subgroups: {type(exc).__name__}: {exc}"
            recs.append(rec)
            continue
        found = None
        for sg in subs:
            try:
                iso.search_method_2(subgroup_idx=sg.index,
                                    distortion_type=["displacive"])
            except Exception:  # noqa: BLE001,S112 - 无模式/计算失败的子群跳过
                continue
            if not iso.mode_displacements:
                continue
            label = next(iter(iso.mode_displacements))
            try:
                d = iso.generate_distortion(irrep_label=label, amplitude=0.1)
                rec["local_sg"] = sg.space_group_number
                rec["local_symbol"] = sg.space_group_symbol
                rec["distorted_sg"] = _sg_of(d)
                found = sg
                break
            except Exception:  # noqa: BLE001,S112 - 生成失败的子群跳过
                continue
        rec["match_official"] = (found is not None
                                 and rec.get("distorted_sg") == expected)
        recs.append(rec)
    return recs


def _web_http_check(name: str, path: Path, expected: set[int],
                    port: int) -> dict:
    """网页 HTTP 全流程抽查（与浏览器相同的端点调用）。"""
    rec: dict = {"name": name}
    content = path.read_text(encoding="utf-8")

    def _get(p: str) -> dict:
        with urllib.request.urlopen(f"http://127.0.0.1:{port}{p}") as resp:
            return json.loads(resp.read().decode("utf-8"))

    def _post(p: str, payload: dict) -> dict:
        body = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            f"http://127.0.0.1:{port}{p}", data=body,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req) as resp:  # noqa: S310 - 本地测试服务
            return json.loads(resp.read().decode("utf-8"))

    data = _post("/api/load_cif?lang=en",
                 {"filename": f"{name}.cif", "content": content})
    rec["load_ok"] = bool(data.get("ok"))
    st = data.get("state", {})
    rec["web_sg"] = (st.get("structure") or {}).get("space_group_number")
    if not rec["load_ok"]:
        rec["error"] = data.get("error", "load_cif failed")
        rec["pass"] = False
        return rec

    opts = _get("/api/method1_options")
    o = opts.get("options", {})
    rec["reachable"] = len(o.get("space_groups", []))
    nums = [g["number"] for g in o.get("space_groups", [])]
    rec["sorted_ok"] = nums == sorted(nums)

    m1 = _post("/api/method1",
               {"distortion_types": ["displacive", "strain"],
                "maximal_subgroup_only": False})
    cands = m1.get("candidates", [])
    rec["m1_candidates"] = len(cands)
    rec["m1_ok"] = bool(m1.get("ok"))

    m2 = None
    if cands:
        m2 = _post("/api/method2",
                   {"subgroup_idx": cands[0]["index"],
                    "distortion_type": ["displacive"],
                    "scope": {"displacive": ["*"], "occupational": [],
                              "strain": [], "magnetic": [], "rotational": []},
                    "nmod": 0, "nsup": 1})
    rec["m2_ok"] = m2 is not None and bool(m2.get("ok"))
    rec["m2_modes"] = len((m2 or {}).get("modes", []) or [])
    if m2 and not m2.get("ok"):
        rec["error"] = f"method2: {m2.get('error')}"

    gen = None
    if rec["m2_ok"] and rec["m2_modes"] > 0:
        modes = [m for m in m2["modes"] if m.get("mode_type") == "displacive"]
        if modes:
            gen = _post("/api/generate",
                        {"irrep_label": modes[0]["irrep_label"],
                         "amplitude": 0.1})
            rec["gen_ok"] = bool(gen.get("ok"))
            rec["gen_atoms"] = gen.get("atoms")
            if not gen.get("ok"):
                rec["error"] = f"generate: {gen.get('error')}"
    rec["pass"] = (rec["load_ok"] and rec["web_sg"] in expected
                   and rec["sorted_ok"] and rec["m1_ok"]
                   and rec["m2_ok"] and (rec["m2_modes"] == 0
                                         or rec.get("gen_ok")))
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="只验证指定名称")
    args = ap.parse_args()

    httpd = None
    port = None
    try:
        httpd = web_server.HTTPServer(("127.0.0.1", 0), web_server.IsoHandler)
        port = httpd.server_address[1]
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        time.sleep(0.3)
    except Exception as exc:  # noqa: BLE001
        print(f"警告: 网页服务启动失败，跳过 HTTP 抽查: {exc}", flush=True)

    results = []
    for name, sgs in EXPECTED.items():
        if args.only and name != args.only:
            continue
        path = CIFS / f"{name}.cif"
        if not path.exists():
            print(f"[--] {name}: 无 CIF（COD 未获取到合格条目）", flush=True)
            continue
        rec: dict = {"name": name, "expected_sgs": sorted(sgs)}
        try:
            api = _validate_parent(name, path, next(iter(sgs)))
        except Exception as exc:  # noqa: BLE001
            api = {"error": f"api: {type(exc).__name__}: {exc}"}
        rec["api"] = {k: v for k, v in api.items()
                      if k in ("identified_sg", "reachable_sgs", "candidates",
                               "symmetry_check", "error", "attempts")}
        rec["sg_ok"] = api.get("identified_sg") in sgs

        if name == "SrTiO3_Pm3m":
            rec["official_paths"] = _official_paths_check(path)
            rec["official_paths_ok"] = all(
                p.get("match_official") for p in rec["official_paths"])

        if port is not None:
            try:
                web_server._SESSION = web_server.WebSession()
                rec["web"] = _web_http_check(name, path, sgs, port)
            except Exception as exc:  # noqa: BLE001
                rec["web"] = {"error": f"{type(exc).__name__}: {exc}",
                              "pass": False}

        # 汇总判定：symmetry_check=True 为通过；False 时进一步归因——
        # 全部尝试的畸变结构都被 spglib 判为 P1（对称性完全破碎）且非刚性
        # 平移、非共享字母时，归类为“取位/原点约定差异”：
        # iso 的 BUSH 符号化代表点（如 32e 的 'x,-x+1/4,...'）与 COD 记录的
        # 轨道设置不一致（MgAl2O4 尖晶石实测），任何参数取值都无法对齐两套
        # 轨道，映射结果错误但被明确捕获，不属于算法静默 bug。
        check = api.get("symmetry_check")
        rec["classification"] = check
        if check is False:
            attempts = api.get("attempts") or []
            all_p1 = bool(attempts) and all(
                e.get("distorted_sg") == 1 for e in attempts
                if e.get("distorted_sg") is not None)
            not_rigid = bool(attempts) and not all(
                e.get("rigid_translation") for e in attempts
                if not e.get("symmetry_ok"))
            if all_p1 and not_rigid and not api.get("shared_letter"):
                rec["classification"] = "limitation-origin-setting"
        rec["pass"] = bool(rec.get("sg_ok") and check is True)
        if name == "SrTiO3_Pm3m":
            rec["pass"] = rec["pass"] and rec.get("official_paths_ok")
        if "web" in rec:
            rec["pass"] = rec["pass"] and rec["web"].get("pass") is True
        results.append(rec)
        tag = "PASS" if rec["pass"] else "FAIL"
        print(f"[{tag}] {name}: sg={api.get('identified_sg')} "
              f"exp={sorted(sgs)} check={check} "
              f"class={rec.get('classification')} "
              f"web_pass={rec.get('web', {}).get('pass')} "
              f"official_ok={rec.get('official_paths_ok')}", flush=True)

    if httpd is not None:
        httpd.shutdown()
        httpd.server_close()

    out = ROOT / "output" / "external_validation_report.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    n_pass = sum(1 for r in results if r.get("pass"))
    n_tot = len(results)
    print(f"\n=== 外部结构源验证汇总: {n_pass}/{n_tot} 通过 ===")
    print(f"  报告: {out}")
    return 0 if n_pass == n_tot else 1


if __name__ == "__main__":
    raise SystemExit(main())
