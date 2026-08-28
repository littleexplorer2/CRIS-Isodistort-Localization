"""Compare official Method1 X/M (and remaining) vs local export folders."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3] / "output_compare"
OFF = ROOT / "官网" / "Method1"
LOC = ROOT / "现有网页版交互" / "Method1"


def _cif_keys(text: str) -> dict:
    out: dict = {"sites": []}
    for ln in text.splitlines():
        if ln.startswith("_symmetry_space_group_name_H-M"):
            out["hm"] = ln
        elif ln.startswith("_symmetry_Int_Tables_number"):
            out["sg"] = ln
        elif ln.startswith("_iso_parent-to-child.transform_Pp_abc"):
            out["p2c"] = ln
        elif ln.startswith("_iso_displacivemode_number"):
            out["ndisp"] = ln
        elif ln.startswith("_iso_strainmode_number"):
            out["nstrain"] = ln
        elif ln.startswith("# k point:"):
            out["k"] = ln
        elif ln.startswith("# IR:"):
            out["ir"] = ln
        elif re.match(r"^# (P|C)\d", ln) or (
            ln.startswith("# ") and "basis={" in ln and "k-active=" in ln and not ln.startswith("# GM")
            and not ln.startswith("# X") and not ln.startswith("# M") and not ln.startswith("# N")
            and not ln.startswith("# P") and not ln.startswith("# LD")
        ):
            # OPD body line (no irrep) OR legacy with irrep
            if "basis={" in ln:
                out["opd"] = ln
        elif re.match(r"^# (GM|X|M|N|P|LD)\S*\s+(P|C)\d", ln) and "basis={" in ln:
            out["opd"] = ln
        elif re.match(r"^(Eu|Al)\S*\s+(Eu|Al)\s+\d+\s+[a-z]\s+", ln):
            if len(out["sites"]) < 6:
                out["sites"].append(ln.strip())
    return out


def _prefix(name: str) -> str:
    # "X1+ P3 (a;0) ..." or "M5+ C1 (a,b) ..."
    m = re.match(r"^(\S+)\s+(\S+)", name)
    return f"{m.group(1)} {m.group(2)}" if m else name[:20]


def _match_loc(od: Path) -> Path | None:
    pref = _prefix(od.name)
    cands = [p for p in LOC.iterdir() if p.is_dir() and p.name.startswith(pref)]
    if not cands:
        # try IR only + OPD token
        parts = pref.split()
        if len(parts) >= 2:
            cands = [
                p for p in LOC.iterdir()
                if p.is_dir() and p.name.startswith(parts[0] + " ") and parts[1] in p.name.split()[:3]
            ]
    if not cands:
        return None
    # prefer same SG number
    sm = re.search(r"\b(\d{1,3})\b", od.name)
    if sm:
        for c in cands:
            if f" {sm.group(1)} " in f" {c.name} " or c.name.split(",")[0].endswith(sm.group(1)):
                # also check OPD dir fragment
                return c
        for c in cands:
            if sm.group(1) in c.name:
                return c
    return cands[0]


def compare(od: Path, ld: Path | None) -> list[str]:
    lines = [f"\n===== {od.name[:70]} ====="]
    if ld is None:
        lines.append("LOC: MISSING")
        return lines
    lines.append(f"LOC: {ld.name[:70]}")
    # folder naming
    if "  " in ld.name.replace("  ", "X"):  # noop
        pass
    if re.search(r"\d \d", ld.name) and "/" not in ld.name:
        # possible slash→space leftover
        if re.search(r"(-?\d) (\d)", ld.name):
            lines.append("FOLDER: possible slash→space fraction")
    if re.search(r"[IPCF]4 mmm|[IPCF]2 m\b", ld.name):
        lines.append("FOLDER: HM has spaces from slash strip")
    oc = od / "subgroup.cif"
    lc = ld / "subgroup.cif"
    if not oc.is_file() or not lc.is_file():
        lines.append("CIF missing")
        return lines
    ok = _cif_keys(oc.read_text(encoding="utf-8", errors="replace"))
    lk = _cif_keys(lc.read_text(encoding="utf-8", errors="replace"))
    for k in ("hm", "sg", "p2c", "ndisp", "nstrain", "k", "ir", "opd"):
        ov, lv = ok.get(k), lk.get(k)
        if ov != lv:
            lines.append(f"  DIFF {k}:")
            lines.append(f"    OFF: {ov}")
            lines.append(f"    LOC: {lv}")
        else:
            lines.append(f"  OK {k}")
    osites, lsites = ok.get("sites") or [], lk.get("sites") or []
    if osites != lsites:
        lines.append("  DIFF sites:")
        for a, b in zip(osites, lsites):
            if a != b:
                lines.append(f"    OFF: {a}")
                lines.append(f"    LOC: {b}")
    else:
        lines.append(f"  OK sites ({len(osites)})")
    return lines


def main() -> None:
    pairs = []
    for od in sorted(OFF.iterdir()):
        if not od.is_dir():
            continue
        ir = od.name.split()[0]
        if ir.startswith(("X", "M")) or ir in {"N1+", "N2+", "P3", "P5"}:
            pairs.append((od, _match_loc(od)))
    print(f"Comparing {len(pairs)} official folders (X/M + sample N/P)")
    all_lines: list[str] = []
    for od, ld in pairs:
        all_lines.extend(compare(od, ld))
    out = Path(__file__).resolve().parent / "_cmp_xm_report.txt"
    out.write_text("\n".join(all_lines), encoding="utf-8")
    print(f"Wrote {out}")
    # summary counts
    diff_hm = sum(1 for ln in all_lines if ln.startswith("  DIFF hm"))
    diff_p2c = sum(1 for ln in all_lines if ln.startswith("  DIFF p2c"))
    diff_sites = sum(1 for ln in all_lines if ln.startswith("  DIFF sites"))
    diff_ndisp = sum(1 for ln in all_lines if ln.startswith("  DIFF ndisp"))
    diff_opd = sum(1 for ln in all_lines if ln.startswith("  DIFF opd"))
    ok_hm = sum(1 for ln in all_lines if ln.startswith("  OK hm"))
    print(f"summary: DIFF hm={diff_hm} p2c={diff_p2c} sites={diff_sites} ndisp={diff_ndisp} opd={diff_opd}; OK hm={ok_hm}")
    # print first 120 lines
    print("\n".join(all_lines[:120]))


if __name__ == "__main__":
    main()
