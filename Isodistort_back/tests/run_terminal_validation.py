"""终端交互验证：对 30-CIF 母相逐个脚本化运行 main_terminal.py，并与
API/网页验证报告（output/validation_report_30.json）比对。

比对项（同一输入下结果必须一致）：
- Method 1 候选数
- Method 2 选中的子群（SG/k/IR）与模式数
- 生成畸变结构的空间群（spglib 校验，== 目标子群）

用法：python tests/run_terminal_validation.py [--report output/validation_report_30.json]
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
from isocore.structure import read_cif  # noqa: E402

PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"
TERMINAL = ROOT / "main_terminal.py"
DATA_DIR = Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码")
CIFS_DIR = ROOT / "tests" / "cifs_30"


def _parent_list() -> list[tuple[str, Path]]:
    parents = [(p.stem, p) for p in sorted(CIFS_DIR.glob("sg*.cif"))]
    parents += [("EuAl4_Parent", DATA_DIR / "EuAl4 Parent.cif"),
                ("EuAl4_Springer", DATA_DIR / "EuAl4 Springer.cif")]
    return parents


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
    args = ap.parse_args()
    report = json.loads(Path(args.report).read_text(encoding="utf-8"))
    by_name = {r["name"]: r for r in report["parents"]}

    results = []
    t0 = time.time()
    for i, (name, path) in enumerate(_parent_list(), start=1):
        api = by_name.get(name, {})
        chosen = api.get("chosen") or {}
        idx = chosen.get("idx")
        if idx is None:
            results.append({"name": name, "note": "API 无可用位移模式，跳过",
                            "skipped": True})
            print(f"[{i}/32] {name}: skip（API 无可用模式）", flush=True)
            continue
        script = _script(str(path), idx)
        try:
            stdout = _run_terminal(script)
        except subprocess.TimeoutExpired:
            results.append({"name": name, "error": "timeout"})
            print(f"[{i}/32] {name}: TIMEOUT", flush=True)
            continue
        parsed = _parse(stdout)
        parsed["name"] = name
        # 畸变结构 SG（从终端导出的 CIF 校验）
        dist_sg = None
        cif_path = parsed.get("cif_path")
        if cif_path and Path(cif_path).exists():
            try:
                dist_sg = SpacegroupAnalyzer(
                    read_cif(cif_path), symprec=1e-3).get_space_group_number()
            except Exception as exc:  # noqa: BLE001 - 解析失败按缺失处理
                _ = exc
        parsed["distorted_sg"] = dist_sg
        parsed["expected_sg"] = chosen.get("sg")
        parsed["symmetry_ok"] = (dist_sg == chosen.get("sg")) if dist_sg else None
        # 与 API 报告比对
        parsed["api_candidates"] = api.get("candidates")
        parsed["api_n_modes"] = chosen.get("n_modes")
        parsed["candidates_match"] = (parsed.get("candidates")
                                      == api.get("candidates"))
        parsed["modes_match"] = (parsed.get("n_modes")
                                 == chosen.get("n_modes"))
        results.append(parsed)
        print(f"[{i}/32] {name}: cand={parsed.get('candidates')}"
              f"/api={api.get('candidates')} modes={parsed.get('n_modes')}"
              f"/api={chosen.get('n_modes')} dist_sg={dist_sg}"
              f"/target={chosen.get('sg')}", flush=True)

    summary = {
        "results": results,
        "elapsed_s": round(time.time() - t0, 1),
        "n_candidates_match": sum(1 for r in results if r.get("candidates_match")),
        "n_modes_match": sum(1 for r in results if r.get("modes_match")),
        "n_sym_ok": sum(1 for r in results if r.get("symmetry_ok") is True),
        "n_ran": sum(1 for r in results if not r.get("skipped")),
    }
    out = Path("output/terminal_validation_report.json")
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print("\n=== 终端验证汇总 ===")
    print(f"  运行: {summary['n_ran']}  候选数一致: "
          f"{summary['n_candidates_match']}  模式数一致: "
          f"{summary['n_modes_match']}  畸变对称性==子群: "
          f"{summary['n_sym_ok']}")
    print(f"  报告: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
