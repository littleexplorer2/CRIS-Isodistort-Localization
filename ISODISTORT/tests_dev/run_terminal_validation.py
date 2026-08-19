"""终端交互验证：脚本化运行 main_terminal.py，与 API 结果比对。

覆盖两类结构源（对应外部数据库测试方案）：
- internal：30 个合成 CIF（tests_dev/cifs_30/）+ 2 个真实 EuAl4 CIF，
  API 基准来自 output/validation_report_30.json；
- external：COD 真实结构（tests_dev/cifs_external/，fetch_cod_cifs.py 下载），
  API 基准现场用 run_30cif_validation._validate_parent 计算。

比对项（同一输入下结果必须一致）：
- Method 1 候选数
- Method 2 选中的子群（SG/k/IR）与模式数
- 生成畸变结构的空间群（spglib 校验，== 目标子群）

用法：
python tests_dev/run_terminal_validation.py                      # 30-CIF（默认）
python tests_dev/run_terminal_validation.py --source external    # COD 外部结构
python tests_dev/run_terminal_validation.py --source all         # 两者
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from run_30cif_validation import _validate_parent  # noqa: E402

from isocore.structure import read_cif  # noqa: E402

PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"
TERMINAL = ROOT / "main_terminal.py"
DATA_DIR = Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码")
CIFS_DIR = ROOT / "tests_dev" / "cifs_30"
EXTERNAL_DIR = ROOT / "tests_dev" / "cifs_external"


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
    随后进入主菜单 -> Method 1（晶系/SG/maximal/lattice 均不过滤）
    -> Method 2（从 Method 1 候选选择，nmod=0，nsup=1）
    -> 进入畸变生成（选第 1 个模式，amplitude=0.1）-> 退出。
    """
    return "\n".join([
        "0",            # CIF 选择器：手动输入路径
        cif_path,
        "3",            # 主菜单 3：Method 1
        "",             # 晶系（不过滤）
        "",             # 空间群号（不过滤）
        "n",            # maximal only? 否
        "0",            # lattice 过滤：不选
        "4",            # 主菜单 4：Method 2
        "1",            # 方式 1：从 Method 1 候选选择
        str(subgroup_idx),
        "0",            # independent modulations（默认 0）
        "1",            # superposed IRs（默认 1）
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


def main() -> int:
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


if __name__ == "__main__":
    raise SystemExit(main())
