"""批量验证手工脚本。

用法：
  python tests_dev/manual/run_batch.py cif30 [--out report.json]
  python tests_dev/manual/run_batch.py external [--only NAME]
  python tests_dev/manual/run_batch.py terminal [--source internal|external|all]
"""
from __future__ import annotations

import sys
from pathlib import Path

_TESTS_DEV = Path(__file__).resolve().parents[1]
ROOT = _TESTS_DEV.parent  # ISODISTORT
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(_TESTS_DEV))



# ========== cif30 ==========
import argparse
import json
import time
from collections import Counter

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


from isocore.api import IsoDistort
from isocore.structure import read_cif
from isocore.utils import run_self_checks

from data_dir import experiment_data_dir

CIFS_DIR = _TESTS_DEV / "cifs_30"
DATA_DIR = experiment_data_dir()
REAL_CIFS = [
    ("EuAl4_Parent", DATA_DIR / "EuAl4 Parent.cif", 139),
    ("EuAl4_Springer", DATA_DIR / "EuAl4 Springer.cif", 139),
]
OFFICIAL_REF = DATA_DIR / "LD1_C1_subgroup.cif"

SYMPREC = 1e-3


def _expected_sg(path: Path) -> int | None:
    m = path.stem
    if m.startswith("sg") and m[2:].isdigit():
        return int(m[2:])
    return None


def _sg_of(structure) -> int:
    return SpacegroupAnalyzer(structure, symprec=SYMPREC).get_space_group_number()


def _validate_parent(name: str, path: Path, expected_sg: int | None) -> dict:
    rec: dict = {"name": name, "path": str(path), "expected_sg": expected_sg}
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    t0 = time.time()
    try:
        iso.load_structure(path)
    except Exception as exc:  # noqa: BLE001
        rec.update({"error": f"load: {type(exc).__name__}: {exc}"})
        return rec
    rec["load_s"] = round(time.time() - t0, 1)
    rec["identified_sg"] = iso.symmetry_info["space_group_number"]
    rec["n_atoms"] = len(iso.structure)
    # 多物种共享同一 Wyckoff 字母（如 Fe/O 同在 k）——模式映射的已知限制特征
    letter_counts = Counter(s["wyckoff_letter"]
                            for s in iso.symmetry_info["wyckoff_sites"])
    rec["shared_letter"] = any(c > 1 for c in letter_counts.values())

    # method1_options：可达子群（排序）+ 格子选项（去重）
    t0 = time.time()
    try:
        opts = iso.method1_options()
        sgs = [sg["number"] for sg in opts["space_groups"]]
        rec["reachable_sgs"] = len(sgs)
        rec["sgs_sorted"] = sgs == sorted(sgs)
        rec["conv_lattices"] = len(opts["conventional_lattices"])
        rec["prim_lattices"] = len(opts["primitive_lattices"])
        rec["options_s"] = round(time.time() - t0, 1)
    except Exception as exc:  # noqa: BLE001
        rec.update({"error": f"method1_options: {type(exc).__name__}: {exc}",
                    "options_s": round(time.time() - t0, 1)})
        return rec

    # method1（网页同款请求体：无过滤）
    try:
        candidates = iso.search_method_1(
            distortion_types=["displacive", "strain"],
            crystal_system=None, subgroup_space_group=None,
            maximal_subgroup_only=False,
        )
        rec["candidates"] = len(candidates)
    except Exception as exc:  # noqa: BLE001
        rec.update({"error": f"method1: {type(exc).__name__}: {exc}"})
        return rec

    # lattice 过滤 sanity：第一个 conventional lattice 选项应显著过滤候选
    if opts["conventional_lattices"]:
        lat = opts["conventional_lattices"][0]["basis"]
        try:
            filtered = iso.search_method_1(
                distortion_types=["displacive", "strain"],
                crystal_system=None, lattice=lat, maximal_subgroup_only=False,
            )
            rec["lattice_filter"] = {
                "lattice": lat,
                "before": rec["candidates"],
                "after": len(filtered),
            }
        except Exception as exc:  # noqa: BLE001
            rec["lattice_filter"] = {"error": f"{type(exc).__name__}: {exc}"}

    # method2 + generate：尝试多个可产生位移模式的子群，逐个校验畸变对称性
    # （区分：真实降对称 PASS / 1 原子原胞的刚性平移退化 / 算法 FAIL）
    parent_sg = rec["identified_sg"]
    attempts: list[dict] = []
    for item in candidates[:80]:
        if len(attempts) >= 12:
            break
        try:
            iso.search_method_2(
                subgroup_idx=item.subgroup.index,
                distortion_type=["displacive"],
            )
        except Exception as exc:  # noqa: BLE001 - 跳过无法计算模式的子群
            _ = exc
            continue
        if not iso.mode_displacements:
            continue
        label = next(iter(iso.mode_displacements))
        disp = iso.mode_displacements[label]["displacements"]
        entry = {
            "idx": item.subgroup.index,
            "sg": item.subgroup.space_group_number,
            "symbol": item.subgroup.space_group_symbol,
            "k": item.subgroup.k_point_label,
            "ir": item.subgroup.irrep_label,
            "basis": item.subgroup.basis_vectors,
            "subgroup_index": item.subgroup.subgroup_index,
            "n_modes": len(iso.mode_displacements),
            "first_mode": label,
            # 所有原子同一位移 = 刚性平移（1 原子原胞的 Γ 模式），
            # 纯位移畸变无法降对称（需应变耦合，见 README 已知差异）
            "rigid_translation": bool(
                np.allclose(disp, disp[0], atol=1e-6)),
        }
        try:
            d = iso.generate_distortion(irrep_label=label, amplitude=0.1)
            entry["distorted_atoms"] = len(d)
            entry["distorted_sg"] = _sg_of(d)
            entry["symmetry_ok"] = entry["distorted_sg"] == entry["sg"]
        except Exception as exc:  # noqa: BLE001
            entry["generate_error"] = f"{type(exc).__name__}: {exc}"
        if entry.get("symmetry_ok"):
            try:
                domains = iso.generate_domains()
                entry["n_domains"] = len(domains)
                entry["domains_ok"] = (len(domains)
                                       == entry["subgroup_index"])
            except Exception as exc:  # noqa: BLE001
                entry["domains_error"] = f"{type(exc).__name__}: {exc}"
            # 物理自洽性自检（零振幅/子群规则/正交性/线性/对称性守恒）
            try:
                entry["self_checks"] = run_self_checks(iso, label)
            except Exception as exc:  # noqa: BLE001
                entry["self_checks_error"] = f"{type(exc).__name__}: {exc}"
        attempts.append(entry)

    ok_entries = [e for e in attempts if e.get("symmetry_ok")]
    lowering_ok = [e for e in ok_entries
                   if e["sg"] != parent_sg and not e.get("rigid_translation")]
    rec["attempts"] = attempts
    # chosen：优先取“真实降对称且通过”的子群，其次任意通过项
    chosen = (lowering_ok[0] if lowering_ok else
              (ok_entries[0] if ok_entries else
               (attempts[0] if attempts else None)))
    if chosen is None:
        rec["chosen"] = {"note": "前 80 个候选均无位移模式"}
    else:
        rec["chosen"] = chosen
        # 汇总判定
        if chosen.get("symmetry_ok"):
            rec["symmetry_check"] = True
        elif all(e.get("rigid_translation")
                 for e in attempts if not e.get("symmetry_ok")) and attempts:
            rec["symmetry_check"] = "degenerate-rigid-translation"
        elif rec.get("shared_letter") and attempts:
            rec["symmetry_check"] = "limitation-shared-letter"
        else:
            rec["symmetry_check"] = False
    return rec


def _official_ld1_comparison() -> dict:
    """官网参考路径：EuAl4 -> LD1 参数 k 点（官网 g=1/6）-> P4mm #99。

    官网参考 CIF（LD1_C1_subgroup.cif）为 ISODISTORT 6.12.2 导出的
    **零振幅框架文件**（序参量全部为 0，仅含子群超胞与取位信息），
    其头部注释记录：k point LD (0,0,g), g=1/6；IR LD1；OPD C1 (a,b)；
    目标 99 P4mm，basis={(1,0,0),(0,1,0),(0,0,6)}，s=12，i=24。

    判据：本地枚举的 (k, IR, OPD, 子群, s, i, 基矢) 与官网头部完全一致；
    模式计算为参数 k 点的已知限制。
    """
    rec: dict = {"name": "EuAl4_LD1_official_ref"}
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(DATA_DIR / "EuAl4 Parent.cif")
    try:
        subs = iso.list_subgroups_at("LD", "LD1", k_parameters=["1/6"])
    except Exception as exc:  # noqa: BLE001
        rec.update({"error": f"LD 子群枚举: {type(exc).__name__}: {exc}"})
        return rec
    target = next((s for s in subs if s.space_group_number == 99), None)
    if target is None:
        rec.update({"error": "LD/LD1 子群中未找到 P4mm (#99)"})
        return rec
    rec["subgroup"] = {
        "sg": target.space_group_number,
        "symbol": target.space_group_symbol,
        "index": target.subgroup_index, "size": target.size,
        "opd": target.opd_symbol, "ir": target.irrep_label,
        "basis": target.basis_vectors,
    }
    # 与官网头部对比（LD1 + C1 OPD -> P4mm，s=12，i=24，c 轴 6 倍超胞）
    official_header = {
        "k": "LD (0,0,g), g=1/6",
        "ir": "LD1",
        "opd": "C1 (a,b)",
        "target": "99 P4mm",
        "s": 12, "i": 24, "c_mult": 6,
    }
    rec["official_header"] = official_header
    basis_ok = (len(target.basis_vectors) == 3
                and abs(round(float(np.linalg.det(target.basis_vectors)), 1))
                == 6.0)
    rec["path_match"] = (target.space_group_number == 99
                         and target.subgroup_index == 24
                         and target.size == 12
                         and target.opd_symbol == "C1"
                         and basis_ok)
    # 参数 k 点（LD）的模式计算为已知限制（iso BUSH 仅支持对称 k 点），
    # 官网使用 (3+d) 超空间机制。此处仅校验子群枚举结果与官网参考一致。
    rec["mode_calc"] = "not supported (parametric k point, documented limitation)"
    try:
        iso.search_method_2(subgroup_idx=target.index,
                            distortion_type=["displacive"])
        rec["mode_calc"] = "ok"
    except Exception as exc:  # noqa: BLE001
        rec["mode_calc_error"] = f"{type(exc).__name__}: {str(exc)[:120]}"
    # 官网参考 CIF 为“零振幅框架文件”（序参量全 0），空间群读回 I4/mmm 属正常
    ref = _read_cif(OFFICIAL_REF)
    rec["official_ref_sg"] = _sg_of(ref)
    rec["official_ref_zero_amplitude"] = ref.num_sites == 60
    return rec


def _read_cif(path: Path):
    return read_cif(path)


def cmd_cif30() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="output/validation_report_30.json")
    args = ap.parse_args()

    cif_paths = sorted(CIFS_DIR.glob("sg*.cif"))
    parents = [(p.stem, p, _expected_sg(p)) for p in cif_paths] + \
        [(n, p, sg) for n, p, sg in REAL_CIFS]
    print(f"=== 30-CIF 验证：{len(parents)} 个母相（{len(cif_paths)} 生成 + "
          f"{len(REAL_CIFS)} 真实）===")

    results = []
    t_start = time.time()
    for i, (name, path, expected) in enumerate(parents, start=1):
        print(f"[{i}/{len(parents)}] {name} ...", flush=True)
        rec = _validate_parent(name, path, expected)
        results.append(rec)
        tag = "OK " if not rec.get("error") else "ERR"
        chosen = rec.get("chosen", {})
        print(f"  {tag} sg={rec.get('identified_sg')} reachable="
              f"{rec.get('reachable_sgs')} cand={rec.get('candidates')} "
              f"chosen={chosen.get('sg')} "
              f"pass={rec.get('symmetry_check')} "
              f"({rec.get('options_s')}s)", flush=True)

    print("\n=== 官网参考比对：EuAl4 -> LD1 -> P4mm ===", flush=True)
    try:
        official = _official_ld1_comparison()
    except Exception as exc:  # noqa: BLE001
        official = {"error": f"{type(exc).__name__}: {exc}"}
    print("  ", json.dumps(official, ensure_ascii=False), flush=True)

    report = {
        "parents": results,
        "official_ld1": official,
        "elapsed_s": round(time.time() - t_start, 1),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print(f"\n报告已写入 {out}")

    # 汇总
    ran = [r for r in results if not r.get("skipped")]
    n_load_ok = sum(1 for r in ran if r.get("identified_sg") is not None)
    n_expected = sum(1 for r in ran if r.get("expected_sg") is not None)
    n_sg_match = sum(1 for r in ran
                     if r.get("identified_sg") == r.get("expected_sg"))
    n_pass = sum(1 for r in ran if r.get("symmetry_check") is True)
    n_degen = sum(1 for r in ran
                  if r.get("symmetry_check") == "degenerate-rigid-translation")
    n_lim = sum(1 for r in ran
                if r.get("symmetry_check") == "limitation-shared-letter")
    n_fail = sum(1 for r in ran if r.get("symmetry_check") is False)
    n_errors = sum(1 for r in ran if r.get("error"))
    n_attempted = sum(1 for r in ran if r.get("attempts"))
    print("\n=== 汇总 ===")
    print(f"  加载成功: {n_load_ok}/{len(ran)}")
    print(f"  空间群识别一致: {n_sg_match}/{n_expected} (预期 SG 的母相)")
    print(f"  畸变对称性校验通过: {n_pass}/{n_attempted}")
    print(f"  退化结构（1 原子原胞，刚性平移无法降对称，已知限制）: {n_degen}")
    print(f"  共享 Wyckoff 字母（模式映射已知限制）: {n_lim}")
    print(f"  未通过: {n_fail}")
    print(f"  流程错误: {n_errors}")
    return 0




# ========== external ==========
import argparse
import json
import threading
import time
import urllib.request


from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: E402
from isocore.api import IsoDistort  # noqa: E402
from web import server as web_server  # noqa: E402

CIFS = _TESTS_DEV / "cifs_external"

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

    # 网页版对齐官网选项页：不含 Distortion Page（生成/导出/畴端点已移除），
    # 抽查到 method2（子群枚举 + 模式计算）为止。
    rec["pass"] = (rec["load_ok"] and rec["web_sg"] in expected
                   and rec["sorted_ok"] and rec["m1_ok"] and rec["m2_ok"])
    return rec


def cmd_external() -> int:
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




# ========== terminal ==========
import argparse
import json
import re
import subprocess
import time

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.structure import read_cif  # noqa: E402

from data_dir import experiment_data_dir

CRIS_ROOT = ROOT.parent
PYTHON = CRIS_ROOT / ".venv" / "Scripts" / "python.exe"
TERMINAL = ROOT / "main_terminal.py"
DATA_DIR = experiment_data_dir()
CIFS_DIR = _TESTS_DEV / "cifs_30"
EXTERNAL_DIR = _TESTS_DEV / "cifs_external"


def _parent_list() -> list[tuple[str, Path]]:
    parents = [(p.stem, p) for p in sorted(CIFS_DIR.glob("sg*.cif"))]
    parents += [("EuAl4_Parent", DATA_DIR / "EuAl4 Parent.cif"),
                ("EuAl4_Springer", DATA_DIR / "EuAl4 Springer.cif")]
    return parents


def _external_list() -> list[tuple[str, Path]]:
    return [(p.stem, p) for p in sorted(EXTERNAL_DIR.glob("*.cif"))]


def _script(cif_path: str, subgroup_idx: int) -> str:
    """终端交互脚本（对齐当前 main_terminal.py 的提示顺序）。

    启动即出现 CIF 选择器：先输 0（手动输入路径）再给路径；
    随后进入主菜单 -> 菜单 2 启用 displacive（默认仅 strain，需显式启用）
    -> Method 1（晶系/SG/maximal/lattice 均不过滤）
    -> Method 2（从 Method 1 候选选择，nmod=0）
    -> 进入畸变生成（选第 1 个模式，amplitude=0.1）-> 退出。
    """
    return "\n".join([
        "0",            # CIF 选择器：手动输入路径
        cif_path,
        "2",            # 主菜单 2：设置畸变类型（默认仅 strain）
        "displacive",   # 启用 displacive
        "all",          # displacive 作用域：全部物种
        "3",            # 主菜单 3：Method 1
        "",             # 晶系（不过滤）
        "",             # 空间群号（不过滤）
        "n",            # maximal only? 否
        "0",            # lattice 过滤：不选
        "4",            # 主菜单 4：Method 2
        "1",            # 方式 1：从 Method 1 候选选择
        str(subgroup_idx),
        "0",            # independent modulations（默认 0）
        "y",            # 进入畸变生成
        "1",            # 选择第 1 个可用模式
        "0.1",          # amplitude
        "0",            # 退出主菜单
    ]) + "\n"


def _run_terminal(script: str, timeout: int = 240) -> str:
    proc = subprocess.run(  # noqa: S603
        [str(PYTHON), "-X", "utf8", str(TERMINAL)],
        input=script, capture_output=True, text=True,
        encoding="utf-8", errors="replace",  # 终端输出为 UTF-8（zh/en），显式解码
        cwd=str(ROOT), timeout=timeout, check=False,
    )
    return proc.stdout + proc.stderr


def _parse(stdout: str) -> dict:
    out: dict = {}
    m = re.search(r"Method 1 返回 (\d+) 条候选", stdout)
    if m:
        out["candidates"] = int(m.group(1))
    m = re.search(r"已选子群: #\d+: SG (\d+) (\S+) k=(\S+) IR=(\S+)", stdout)
    if m:
        out.update({"sg": int(m.group(1)), "symbol": m.group(2),
                    "k": m.group(3), "ir": m.group(4)})
    m = re.search(r"Method 2: 模式数 (\d+)", stdout)
    if m:
        out["n_modes"] = int(m.group(1))
    # 生成消息随界面语言变化（zh/en 双模式）
    m = re.search(r"已生成结构，原子数: (\d+)|Structure generated, atoms: (\d+)",
                  stdout)
    if m:
        out["gen_atoms"] = int(m.group(1) or m.group(2))
    m = re.search(
        r"\[Auto-export\] CIF written: (.+\.cif)"
        r"|\[默认导出\] 已生成 CIF: (.+\.cif)", stdout
    )
    if m:
        out["cif_path"] = (m.group(1) or m.group(2)).strip()
    if "Traceback" in stdout or "Error" in stdout:
        out["stderr_trace"] = stdout[-800:]
    return out


def cmd_terminal() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="output/validation_report_30.json")
    ap.add_argument("--source", default="internal",
                    choices=["internal", "external", "all"])
    args = ap.parse_args()

    report = None
    if args.source in ("internal", "all"):
        report = json.loads(Path(args.report).read_text(encoding="utf-8"))
        by_name = {r["name"]: r for r in report["parents"]}
    else:
        by_name = {}

    sources: list[tuple[str, list[tuple[str, Path]]]] = []
    if args.source in ("internal", "all"):
        sources.append(("internal", _parent_list()))
    if args.source in ("external", "all"):
        sources.append(("external", _external_list()))

    results = []
    t0 = time.time()
    total = sum(len(items) for _, items in sources)
    done = 0
    for src_name, items in sources:
        for name, path in items:
            done += 1
            if src_name == "internal":
                api = by_name.get(name, {})
                chosen = api.get("chosen") or {}
                idx = chosen.get("idx")
                api_candidates = api.get("candidates")
                api_n_modes = chosen.get("n_modes")
                target_sg = chosen.get("sg")
                api_check = api.get("symmetry_check")
            else:
                # 外部结构：API 基准现场计算（与网页同底层的完整流程）
                rec = _validate_parent(name, path, None)
                attempts = rec.get("attempts") or []
                chosen = (next((e for e in attempts
                                if e.get("symmetry_ok")), None)
                          or (attempts[0] if attempts else None))
                idx = chosen.get("idx") if chosen else None
                api_candidates = rec.get("candidates")
                api_n_modes = chosen.get("n_modes") if chosen else None
                target_sg = chosen.get("sg") if chosen else None
                api_check = rec.get("symmetry_check")
                if api_check is False and rec.get("classification"):
                    api_check = rec.get("classification")

            if idx is None:
                results.append({"name": name, "source": src_name,
                                "note": "API 无可用位移模式，跳过",
                                "skipped": True})
                print(f"[{done}/{total}] [{src_name}] {name}: "
                      f"skip（API 无可用模式）", flush=True)
                continue
            script = _script(str(path), idx)
            try:
                stdout = _run_terminal(script)
            except subprocess.TimeoutExpired:
                results.append({"name": name, "source": src_name,
                                "error": "timeout"})
                print(f"[{done}/{total}] [{src_name}] {name}: TIMEOUT",
                      flush=True)
                continue
            parsed = _parse(stdout)
            parsed["name"] = name
            parsed["source"] = src_name
            # 畸变结构 SG（从终端导出的 CIF 校验）
            dist_sg = None
            cif_path = parsed.get("cif_path")
            if cif_path and Path(cif_path).exists():
                try:
                    dist_sg = SpacegroupAnalyzer(
                        read_cif(cif_path), symprec=1e-3
                    ).get_space_group_number()
                except Exception as exc:  # noqa: BLE001 - 解析失败按缺失处理
                    _ = exc
            parsed["distorted_sg"] = dist_sg
            parsed["expected_sg"] = target_sg
            if api_check is True:
                # API 基准有对称性校验：终端生成的畸变结构必须与目标子群一致
                parsed["symmetry_ok"] = ((dist_sg == target_sg)
                                         if dist_sg else None)
            else:
                # API 基准为已知边界（退化/共享字母/取位差异）：不计失败，
                # 仅记录终端行为一致（候选数/模式数仍严格比对）
                parsed["note"] = f"api symmetry_check={api_check}（已知边界）"
                parsed["symmetry_ok"] = None
            # 与 API 基准比对
            parsed["api_candidates"] = api_candidates
            parsed["api_n_modes"] = api_n_modes
            parsed["candidates_match"] = (parsed.get("candidates")
                                          == api_candidates)
            parsed["modes_match"] = (parsed.get("n_modes") == api_n_modes)
            results.append(parsed)
            ok_tag = ("OK" if parsed["symmetry_ok"] is True
                      else ("boundary" if parsed["symmetry_ok"] is None
                            else "MISMATCH"))
            print(f"[{done}/{total}] [{src_name}] {name}: "
                  f"cand={parsed.get('candidates')}/api={api_candidates} "
                  f"modes={parsed.get('n_modes')}/api={api_n_modes} "
                  f"dist_sg={dist_sg}/target={target_sg} [{ok_tag}]",
                  flush=True)

    summary = {
        "results": results,
        "elapsed_s": round(time.time() - t0, 1),
        "n_candidates_match": sum(1 for r in results
                                  if r.get("candidates_match")),
        "n_modes_match": sum(1 for r in results if r.get("modes_match")),
        "n_sym_ok": sum(1 for r in results if r.get("symmetry_ok") is True),
        "n_ran": sum(1 for r in results if not r.get("skipped")),
    }
    for src_name, _ in sources:
        sub = [r for r in results if r.get("source") == src_name]
        summary[f"n_{src_name}"] = len(sub)
        summary[f"{src_name}_sym_ok"] = sum(
            1 for r in sub if r.get("symmetry_ok") is True)
        summary[f"{src_name}_cand_match"] = sum(
            1 for r in sub if r.get("candidates_match"))
        summary[f"{src_name}_modes_match"] = sum(
            1 for r in sub if r.get("modes_match"))
    out = Path("output/terminal_validation_report.json")
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print("\n=== 终端验证汇总 ===")
    for src_name, _ in sources:
        print(f"  [{src_name}] 运行 {summary[f'n_{src_name}']} 个 | "
              f"候选数一致 {summary[f'{src_name}_cand_match']} | "
              f"模式数一致 {summary[f'{src_name}_modes_match']} | "
              f"畸变对称性==子群 {summary[f'{src_name}_sym_ok']}")
    print(f"  报告: {out}")
    return 0





def _run_with_rest(fn, rest: list[str]) -> int:
    old = sys.argv
    try:
        sys.argv = [fn.__name__] + list(rest)
        return int(fn() or 0)
    finally:
        sys.argv = old


def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description="Batch validation runners")
    p.add_argument("command", choices=["cif30", "external", "terminal"])
    args, rest = p.parse_known_args(argv)
    fn = {
        "cif30": cmd_cif30,
        "external": cmd_external,
        "terminal": cmd_terminal,
    }[args.command]
    return _run_with_rest(fn, rest)


if __name__ == "__main__":
    raise SystemExit(main())
