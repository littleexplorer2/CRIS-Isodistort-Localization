"""Compare official vs local output_compare Method1/2 exports."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3] / "output_compare"
OFF = ROOT / "官网"
LOC = ROOT / "现有网页版交互"


def _cif_keys(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
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
        elif "k-active=" in ln and ln.startswith("#"):
            out["opd"] = ln
        elif re.match(r"^(Eu|Al)\S*\s+(Eu|Al)\s+\d+\s+[a-z]\s+", ln):
            out.setdefault("sites", [])
            if isinstance(out["sites"], list) and len(out["sites"]) < 8:
                out["sites"].append(ln.strip())
    return out


def _pair_method1() -> list[tuple[Path, Path | None]]:
    pairs = []
    for od in sorted((OFF / "Method1").iterdir()):
        if not od.is_dir():
            continue
        # match by leading IR OPD token
        key = od.name.split(",")[0].strip()  # up to first comma-ish
        # better: match irrep+opd prefix
        m = re.match(r"^(\S+\s+\S+)", od.name)
        prefix = m.group(1) if m else od.name[:20]
        cands = [
            p
            for p in (LOC / "Method1").iterdir()
            if p.is_dir() and p.name.startswith(prefix.split()[0])
            and prefix.split()[1] in p.name.split()[1:3]
        ]
        # refine: same irrep and opd at start
        ir, opd = prefix.split()[:2]
        cands = [
            p
            for p in (LOC / "Method1").iterdir()
            if p.is_dir() and p.name.startswith(f"{ir} {opd}")
        ]
        # pick closest by space group number in name
        best = None
        for c in cands:
            if od.name[:30] in c.name or c.name[:30] in od.name:
                best = c
                break
        if best is None and len(cands) == 1:
            best = cands[0]
        elif best is None and cands:
            # same SG number token
            sm = re.search(r"\b(\d{1,3})\b", od.name)
            if sm:
                for c in cands:
                    if f" {sm.group(1)} " in f" {c.name} ":
                        best = c
                        break
            if best is None:
                best = cands[0]
        pairs.append((od, best))
    return pairs


def _pair_method2() -> list[tuple[Path, Path | None]]:
    pairs = []
    for od in sorted((OFF / "Method2").iterdir()):
        if not od.is_dir():
            continue
        loc = LOC / "Method2" / od.name
        pairs.append((od, loc if loc.is_dir() else None))
    return pairs


def compare_cif(off_dir: Path, loc_dir: Path | None, label: str) -> None:
    print(f"\n===== {label} =====")
    print(f"OFF folder: {off_dir.name[:90]}")
    if loc_dir is None:
        print("LOC: MISSING")
        return
    print(f"LOC folder: {loc_dir.name[:90]}")
    # fraction slash in folder?
    if "1 2" in loc_dir.name or "1 4" in loc_dir.name:
        print("FOLDER: local uses space instead of slash in fractions")
    if "I4 mmm" in loc_dir.name or "P4 mmm" in loc_dir.name:
        print("FOLDER: local HM has spaces (official often I4mmm / P4nmm)")
    oc = off_dir / "subgroup.cif"
    lc = loc_dir / "subgroup.cif"
    if not oc.is_file() or not lc.is_file():
        print("CIF missing")
        return
    ok = _cif_keys(oc.read_text(encoding="utf-8", errors="replace"))
    lk = _cif_keys(lc.read_text(encoding="utf-8", errors="replace"))
    for k in ("hm", "sg", "p2c", "ndisp", "nstrain", "k", "opd"):
        ov, lv = ok.get(k), lk.get(k)
        mark = "OK" if ov == lv else "DIFF"
        if ov != lv:
            print(f"  {mark} {k}:")
            print(f"    OFF: {ov}")
            print(f"    LOC: {lv}")
        else:
            print(f"  {mark} {k}: {ov}")
    osites = ok.get("sites") or []
    lsites = lk.get("sites") or []
    if osites != lsites:
        print("  DIFF sites (first rows):")
        for a, b in zip(osites, lsites):
            if a != b:
                print(f"    OFF: {a}")
                print(f"    LOC: {b}")
        if len(osites) != len(lsites):
            print(f"    count OFF={len(osites)} LOC={len(lsites)}")
    else:
        print(f"  OK sites preview: {osites[:3]}")
    # topas / isoviz markers
    ot = (off_dir / "topas.str").read_text(encoding="utf-8", errors="replace")[:400]
    lt = (loc_dir / "topas.str").read_text(encoding="utf-8", errors="replace")[:400]
    print(f"  TOPAS OFF has transformPp: {'transformPp' in ot}; LOC: {'transformPp' in lt}")
    print(f"  TOPAS OFF prm count ~ {ot.count('prm')}; LOC head has mode: {'!a1' in lt or 'a1' in lt}")
    oi = (off_dir / "data.isoviz").read_text(encoding="utf-8", errors="replace")
    li = (loc_dir / "data.isoviz").read_text(encoding="utf-8", errors="replace")
    print(f"  ISOVIZ OFF !isoversion: {'!isoversion' in oi}; LOC: {'!isoversion' in li}")
    print(f"  ISOVIZ sizes OFF={len(oi)} LOC={len(li)}")
    om = list(off_dir.glob("*modes*"))
    lm = list(loc_dir.glob("*modes*"))
    print(f"  MODES files OFF={[p.name for p in om]} LOC={[p.name for p in lm]}")


def main() -> None:
    print("Method1 pairs:", len(_pair_method1()))
    for od, ld in _pair_method1():
        if od.name.startswith("GM"):
            compare_cif(od, ld, f"M1 {od.name[:40]}")
    print("\n\n#### Method2 ####")
    for od, ld in _pair_method2():
        compare_cif(od, ld, f"M2 {od.name}")


if __name__ == "__main__":
    main()
