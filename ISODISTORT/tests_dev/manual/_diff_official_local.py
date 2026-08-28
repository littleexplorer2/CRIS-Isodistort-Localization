"""Byte/structure diff: output_compare official vs local web dumps."""
from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

OC = Path(r"c:\Users\devou\OneDrive\Desktop\CRIS\output_compare")
OFF, LOC = sorted(OC.iterdir(), key=lambda p: len(list((p / "Method1").iterdir()) if (p / "Method1").is_dir() else []))
# Prefer: shorter Chinese name 官网 vs longer 现有网页版交互 — identify by Method1 folder count equality and name length
cands = [p for p in OC.iterdir() if p.is_dir() and (p / "Method1").is_dir()]
cands.sort(key=lambda p: p.name)
# 官网 is typically the shorter dirname
OFF = min(cands, key=lambda p: len(p.name.encode("utf-8", "replace")))
LOC = max(cands, key=lambda p: len(p.name.encode("utf-8", "replace")))
print("OFFICIAL:", OFF.name)
print("LOCAL:   ", LOC.name)

COMPARE_FILES = ("subgroup.cif", "data.isoviz", "topas.str")


def folder_map(method: str) -> dict[str, Path]:
    root = OFF / method
    return {p.name: p for p in root.iterdir() if p.is_dir()}


def loc_map(method: str) -> dict[str, Path]:
    root = LOC / method
    return {p.name: p for p in root.iterdir() if p.is_dir()}


def first_diff(a: bytes, b: bytes, context: int = 40) -> str:
    n = min(len(a), len(b))
    for i in range(n):
        if a[i] != b[i]:
            return (
                f"@{i} off={a[max(0,i-context):i+context]!r} "
                f"loc={b[max(0,i-context):i+context]!r}"
            )
    if len(a) != len(b):
        return f"length off={len(a)} loc={len(b)}"
    return "identical"


def cif_fields(text: str) -> dict:
    out = {}
    m = re.search(r"(?i)_symmetry_Int_Tables_number\s+(\d+)", text)
    if m:
        out["sg"] = m.group(1)
    m = re.search(r"(?i)_symmetry_space_group_name_H-M\s+\"([^\"]+)\"", text)
    if m:
        out["hm"] = m.group(1)
    m = re.search(r"(?i)_space_group_name_Hall\s+\"([^\"]+)\"", text)
    if m:
        out["hall"] = m.group(1)
    m = re.search(r"(?i)_cell_length_a\s+([^\s]+)", text)
    if m:
        out["a"] = m.group(1)
    sites = re.findall(
        r"(?m)^([A-Za-z]+\d+)\s+[A-Za-z]+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)",
        text,
    )
    # better: _atom_site loop
    loop = re.search(
        r"_atom_site_label\s*\n(?:_atom_site[^\n]+\n)*((?:(?!loop_|_[a-z]).*\n)+)",
        text,
        re.I,
    )
    labels = re.findall(r"(?m)^([A-Za-z]+\d+)\s+", text.split("_atom_site")[-1] if "_atom_site" in text else "")
    out["n_atom_lines"] = len(re.findall(r"(?m)^[A-Za-z]+\d+\s+", text))
    # OPD comment
    for line in text.splitlines():
        if "basis={" in line and "k-active" in line:
            out["opd"] = line.lstrip("# ").strip()
            break
    out["n_lines"] = text.count("\n")
    return out


def audit_method(method: str) -> None:
    om, lm = folder_map(method), loc_map(method)
    print(f"\n===== {method}: official={len(om)} local={len(lm)} =====")
    missing = sorted(set(om) - set(lm))
    extra = sorted(set(lm) - set(om))
    print(f"local missing folders: {len(missing)}")
    for n in missing[:15]:
        print(f"  - {n[:100]}")
    print(f"local extra folders: {len(extra)}")
    for n in extra[:15]:
        print(f"  - {n[:100]}")

    # Match by exact name first; else by short key IR+OPD+SG
    def short(name: str) -> str:
        m = re.match(r"(\S+)\s+(\S+)\s+\(([^)]*)\)\s+(\d+)", name)
        if m:
            return f"{m.group(1)}|{m.group(2)}|{m.group(3)}|{m.group(4)}"
        parts = name.split()
        if len(parts) >= 2:
            return f"{parts[0]}|{parts[1]}"
        return name

    o_by_short = defaultdict(list)
    for n in om:
        o_by_short[short(n)].append(n)
    l_by_short = defaultdict(list)
    for n in lm:
        l_by_short[short(n)].append(n)

    exact = set(om) & set(lm)
    print(f"exact-name overlap: {len(exact)}")

    # For non-exact, try short match
    paired = []
    used_l = set()
    for on in om:
        if on in lm:
            paired.append((on, on, "exact"))
            used_l.add(on)
            continue
        sk = short(on)
        cands = [x for x in l_by_short.get(sk, []) if x not in used_l]
        if len(cands) == 1:
            paired.append((on, cands[0], "short"))
            used_l.add(cands[0])
        elif len(cands) > 1:
            paired.append((on, cands[0], f"short-ambig:{len(cands)}"))
            used_l.add(cands[0])
        else:
            paired.append((on, None, "unpaired"))

    unpaired = sum(1 for _, l, _ in paired if l is None)
    print(f"paired: {len(paired)-unpaired} unpaired-off: {unpaired}")

    stats = Counter()
    field_diffs = Counter()
    samples = []

    for on, ln, how in paired:
        if ln is None:
            stats["unpaired"] += 1
            continue
        if how != "exact":
            stats[f"name_diff:{how}"] += 1
            if len(samples) < 8 and how.startswith("short"):
                samples.append(("NAME", on[:80], ln[:80]))
        for fn in COMPARE_FILES:
            op, lp = om[on] / fn, lm[ln] / fn
            if not op.is_file():
                stats[f"off_missing_{fn}"] += 1
                continue
            if not lp.is_file():
                stats[f"loc_missing_{fn}"] += 1
                continue
            ob, lb = op.read_bytes(), lp.read_bytes()
            if ob == lb:
                stats[f"{fn}_byte_ok"] += 1
            else:
                stats[f"{fn}_byte_diff"] += 1
                if fn == "subgroup.cif":
                    of, lf = cif_fields(ob.decode("utf-8", "replace")), cif_fields(
                        lb.decode("utf-8", "replace")
                    )
                    for k in sorted(set(of) | set(lf)):
                        if of.get(k) != lf.get(k):
                            field_diffs[k] += 1
                    if len(samples) < 20:
                        samples.append(
                            (
                                "CIF",
                                on[:70],
                                first_diff(ob, lb),
                                {k: (of.get(k), lf.get(k)) for k in of if of.get(k) != lf.get(k)},
                            )
                        )

    print("stats:", dict(stats))
    print("CIF field mismatch counts:", dict(field_diffs))
    print("\nsamples:")
    for s in samples[:15]:
        print(" ", s)


def main():
    audit_method("Method1")
    audit_method("Method2")
    # LD5 P4 official modes check
    m2 = OFF / "Method2" / "LD5 P4"
    if m2.is_dir():
        print("\nOFF LD5 P4 files:", sorted(x.name for x in m2.iterdir()))
    m2l = LOC / "Method2" / "LD5 P4"
    if m2l.is_dir():
        print("LOC LD5 P4 files:", sorted(x.name for x in m2l.iterdir()))


if __name__ == "__main__":
    main()
