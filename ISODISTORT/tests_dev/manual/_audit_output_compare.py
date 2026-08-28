"""Audit output_compare/官网 Method1|2 vs webpage_info/a. (and Method2 OPD page).

Compares after normalizing Windows download folder naming:
  '/' deleted from HM (I4/mmm -> I4mmm), fractions compacted (1/2 -> 12).
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from html import unescape
from pathlib import Path

ROOT = Path(r"c:\Users\devou\OneDrive\Desktop\CRIS")
WPI = ROOT / "webpage_info"
OC = ROOT / "output_compare" / "官网"
A_HTML = WPI / "a. ISODISTORT_ order parameter direction.html"

MODES_PAT = re.compile(r"(?i)complete\s*modes|ISODISTORT_ complete modes")

# Capture one OPD radio/label line; k-active stops at quote/tag/newline.
OPD_RE = re.compile(
    r"(?P<head>"
    r"(?:GM|X|M|N|P|H|LD|SM|DT|Z|A|R|T|W|V|U|S|Y|B|F|G|K|L|Q)\S*"
    r"\s+\S+\s+\([^)]*\)\s+\d+\s+[A-Za-z0-9_/\-]+)"
    r",\s*basis=\{(?P<basis>[^}]*)\},\s*"
    r"origin=(?P<origin>\([^)]*\)),\s*"
    r"s=(?P<s>\d+),\s*i=(?P<i>\d+),\s*"
    r"k-active=\s*(?P<kactive>(?:\([^)]*\)(?:,\s*)?)+)",
    re.I,
)

FOLDER_RE = re.compile(
    r"^(?P<head>.+?),\s*basis=\{(?P<basis>[^}]*)\},\s*"
    r"origin=(?P<origin>\([^)]*\)),\s*"
    r"s=(?P<s>\d+),\s*i=(?P<i>\d+),\s*"
    r"k-active=\s*(?P<kactive>.+)$"
)

# CIF comment omits irrep: "P1 (a) 139 I4/mmm, basis=..."
CIF_OPD_RE = re.compile(
    r"(?P<head>\S+\s+\([^)]*\)\s+\d+\s+[A-Za-z0-9_/\-]+),\s*"
    r"basis=\{(?P<basis>[^}]*)\},\s*"
    r"origin=(?P<origin>\([^)]*\)),\s*"
    r"s=(?P<s>\d+),\s*i=(?P<i>\d+),\s*"
    r"k-active=\s*(?P<kactive>.+)",
)


def compact_fractions(s: str) -> str:
    """1/2->12, 1/4->14, 3/4->34, 7/4->74, -1/2->-12, etc. Also strip HM '/'."""
    s = s.replace("\u2212", "-")
    # Compact signed fractions first: -1/2, 3/4, 1/2
    s = re.sub(
        r"(?<![A-Za-z0-9_])(-?\d+)/(\d+)(?![A-Za-z0-9_])",
        lambda m: f"{m.group(1)}{m.group(2)}",
        s,
    )
    # Remove remaining '/' (HM symbols I4/mmm -> I4mmm, P4_2/nmc -> P4_2nmc)
    s = s.replace("/", "")
    s = re.sub(r"\s+", "", s)
    return s


def fields_key(head: str, basis: str, origin: str, s: str, i: str, kactive: str) -> str:
    return compact_fractions(
        f"{head}|{basis}|{origin}|{s}|{i}|{kactive}"
    )


def display_line(head, basis, origin, s, i, kactive) -> str:
    return (
        f"{head}, basis={{{basis}}}, origin={origin}, "
        f"s={s}, i={i}, k-active= {kactive.strip()}"
    )


def parse_opd_text(text: str) -> list[dict]:
    text = unescape(text)
    out = []
    seen = set()
    for m in OPD_RE.finditer(text):
        key = fields_key(
            m.group("head"),
            m.group("basis"),
            m.group("origin"),
            m.group("s"),
            m.group("i"),
            m.group("kactive"),
        )
        if key in seen:
            continue
        seen.add(key)
        head = re.sub(r"\s+", " ", m.group("head").strip())
        out.append(
            {
                "key": key,
                "head": head,
                "basis": m.group("basis"),
                "origin": m.group("origin"),
                "s": m.group("s"),
                "i": m.group("i"),
                "kactive": m.group("kactive").strip(),
                "display": display_line(
                    head,
                    m.group("basis"),
                    m.group("origin"),
                    m.group("s"),
                    m.group("i"),
                    m.group("kactive"),
                ),
                "irrep": head.split()[0],
                "sg": re.search(r"(\d+)\s+\S+$", head).group(1)
                if re.search(r"(\d+)\s+\S+$", head)
                else None,
            }
        )
    return out


def parse_folder(name: str) -> dict | None:
    m = FOLDER_RE.match(name.strip())
    if not m:
        return None
    head = re.sub(r"\s+", " ", m.group("head").strip())
    return {
        "key": fields_key(
            head,
            m.group("basis"),
            m.group("origin"),
            m.group("s"),
            m.group("i"),
            m.group("kactive"),
        ),
        "head": head,
        "basis": m.group("basis"),
        "origin": m.group("origin"),
        "s": m.group("s"),
        "i": m.group("i"),
        "kactive": m.group("kactive").strip(),
        "irrep": head.split()[0],
        "sg": re.search(r"(\d+)\s+\S+$", head).group(1)
        if re.search(r"(\d+)\s+\S+$", head)
        else None,
        "name": name,
    }


def head_without_irrep(head: str) -> str:
    parts = head.split(None, 1)
    return parts[1] if len(parts) == 2 else head


def audit_method(
    method_dir: Path,
    expected: list[dict],
    *,
    label: str,
    allow_empty_expected: bool = False,
):
    expected_keys = {e["key"] for e in expected}
    expected_disp = {e["key"]: e["display"] for e in expected}

    folders = [p for p in method_dir.iterdir() if p.is_dir()]
    print(f"\n===== {label}: {len(folders)} folders | expected {len(expected)} =====")

    metas = []
    bad = []
    for p in sorted(folders, key=lambda x: x.name.lower()):
        meta = parse_folder(p.name)
        if not meta:
            bad.append(p.name)
            continue
        metas.append((p, meta))

    got = {m["key"] for _, m in metas}
    missing = sorted(expected_keys - got)
    extra = sorted(got - expected_keys) if expected_keys or not allow_empty_expected else []

    print(f"parsable: {len(metas)}; unparsable: {len(bad)}")
    if bad:
        for n in bad[:10]:
            print(f"  BAD NAME: {n}")

    print(f"MISSING vs reference: {len(missing)}")
    for k in missing[:25]:
        print(f"  - {expected_disp.get(k, k)[:130]}")
    if len(missing) > 25:
        print(f"  ... +{len(missing) - 25}")

    print(f"EXTRA (not in reference): {len(extra)}")
    key_to_name = {m["key"]: p.name for p, m in metas}
    for k in extra[:15]:
        print(f"  - {key_to_name.get(k, k)[:130]}")
    if len(extra) > 15:
        print(f"  ... +{len(extra) - 15}")

    # Duplicate keys
    key_counts = Counter(m["key"] for _, m in metas)
    dups = [k for k, c in key_counts.items() if c > 1]
    if dups:
        print(f"DUPLICATE keys: {len(dups)}")
        for k in dups[:10]:
            print(f"  - {key_to_name.get(k, k)[:100]}")

    incomplete = []
    mismatches = []
    suspect_wrong_download = []
    file_freq = Counter()

    for p, meta in metas:
        files = [f for f in p.iterdir() if f.is_file()]
        for f in files:
            file_freq[f.name] += 1
        names = {f.name for f in files}
        miss_pay = []
        if "subgroup.cif" not in names:
            miss_pay.append("subgroup.cif")
        if "data.isoviz" not in names:
            miss_pay.append("data.isoviz")
        if "topas.str" not in names:
            miss_pay.append("topas.str")
        modes = [f for f in files if MODES_PAT.search(f.name)]
        if not modes:
            miss_pay.append("Complete modes details")
        if miss_pay:
            incomplete.append((p.name, miss_pay, sorted(names)))

        issues = []
        cif_path = p / "subgroup.cif"
        if cif_path.is_file():
            cif = cif_path.read_text(encoding="utf-8", errors="replace")
            sg_m = re.search(r"(?i)_space_group_IT_number\s+(\d+)", cif)
            if not sg_m:
                sg_m = re.search(r"(?i)_symmetry_Int_Tables_number\s+(\d+)", cif)
            if sg_m and meta["sg"] and sg_m.group(1) != meta["sg"]:
                issues.append(f"CIF SG#{sg_m.group(1)} != folder SG#{meta['sg']}")
                suspect_wrong_download.append(p.name)

            # Find OPD comment line
            opd_line = None
            for line in cif.splitlines():
                if "basis={" in line and "k-active" in line:
                    opd_line = line.lstrip("# ").strip()
                    break
            if opd_line:
                cm = CIF_OPD_RE.search(opd_line)
                if cm:
                    cif_key_body = fields_key(
                        cm.group("head"),
                        cm.group("basis"),
                        cm.group("origin"),
                        cm.group("s"),
                        cm.group("i"),
                        cm.group("kactive"),
                    )
                    folder_body = fields_key(
                        head_without_irrep(meta["head"]),
                        meta["basis"],
                        meta["origin"],
                        meta["s"],
                        meta["i"],
                        meta["kactive"],
                    )
                    if cif_key_body != folder_body:
                        issues.append(
                            "CIF OPD fields != folder name "
                            f"(cif_head={cm.group('head')})"
                        )
                        suspect_wrong_download.append(p.name)

        if modes:
            html = modes[0].read_text(encoding="utf-8", errors="replace")
            # Prefer full OPD with irrep from modes page
            found = parse_opd_text(html)
            if found:
                # pick best: matching irrep or first
                pick = None
                for it in found:
                    if it["irrep"] == meta["irrep"]:
                        pick = it
                        break
                if pick is None:
                    pick = found[0]
                if pick["key"] != meta["key"]:
                    # Compare without requiring exact if only HM slash differs — keys already normed
                    issues.append(
                        f"modes OPD != folder "
                        f"(modes={pick['display'][:80]})"
                    )
                    suspect_wrong_download.append(p.name)
            else:
                # at least check irrep string appears
                if meta["irrep"] not in html and meta["irrep"].replace("+", "") not in html:
                    # weak signal
                    if f" {meta['sg']} " not in html and f">{meta['sg']}" not in html:
                        issues.append("modes HTML lacks folder SG/irrep markers")

        if issues:
            mismatches.append((p.name, issues))

    print(f"\npayload files present counts:")
    for n in ("subgroup.cif", "data.isoviz", "topas.str", "ISODISTORT_ complete modes details.html"):
        print(f"  {file_freq.get(n, 0):4d}/{len(metas)}  {n}")
    other = [n for n in file_freq if n not in {
        "subgroup.cif", "data.isoviz", "topas.str",
        "ISODISTORT_ complete modes details.html",
    }]
    if other:
        print(f"  other file names: {sorted(other)[:20]}")

    print(f"\nincomplete payloads: {len(incomplete)}")
    for name, miss, have in incomplete[:20]:
        print(f"  {name[:110]}")
        print(f"    missing={miss} have={have}")

    print(f"content mismatches: {len(mismatches)}")
    for name, iss in mismatches[:25]:
        print(f"  {name[:110]}")
        for i in iss:
            print(f"    - {i}")
    if len(mismatches) > 25:
        print(f"  ... +{len(mismatches) - 25}")

    # Unique wrong-download suspects
    suspects = sorted(set(suspect_wrong_download))
    print(f"LIKELY wrong-download / misnamed: {len(suspects)}")
    for n in suspects[:30]:
        print(f"  ! {n[:120]}")

    return {
        "n": len(metas),
        "missing": missing,
        "extra": extra,
        "incomplete": incomplete,
        "mismatches": mismatches,
        "suspects": suspects,
        "bad": bad,
    }


def main() -> int:
    a_html = A_HTML.read_text(encoding="utf-8", errors="replace")
    a_items = parse_opd_text(a_html)
    print(f"webpage_info/a. unique OPD lines: {len(a_items)}")
    fam = Counter(it["irrep"] for it in a_items)
    print("irrep counts:", dict(sorted(fam.items())))

    m1 = OC / "Method1"
    r1 = audit_method(m1, a_items, label="Method1 vs a.")

    # Method2: page 4 is OPD for LD path in webpage_info numbered flow
    m2_pages = [
        WPI / "4. ISODISTORT_ order parameter direction.html",
    ]
    m2_items: list[dict] = []
    src = None
    for p in m2_pages:
        if p.exists():
            m2_items = parse_opd_text(p.read_text(encoding="utf-8", errors="replace"))
            src = p.name
            break
    print(f"\nMethod2 OPD reference: {src} -> {len(m2_items)} lines")
    if m2_items:
        print("Method2 irreps:", dict(Counter(it["irrep"] for it in m2_items)))

    m2 = OC / "Method2"
    if m2.is_dir():
        r2 = audit_method(
            m2,
            m2_items,
            label=f"Method2 vs {src}",
            allow_empty_expected=not m2_items,
        )
    else:
        print("NO Method2 directory")
        r2 = None

    # Also: Method1 folders that look like Method2-only (LD) by irrep
    ld_in_m1 = [p.name for p in m1.iterdir() if p.is_dir() and p.name.startswith("LD")]
    if ld_in_m1:
        print(f"\nLD folders wrongly under Method1: {len(ld_in_m1)}")
        for n in ld_in_m1:
            print(f"  {n}")

    print("\n===== SUMMARY =====")
    print(
        f"Method1: have={r1['n']} expect={len(a_items)} "
        f"missing={len(r1['missing'])} extra={len(r1['extra'])} "
        f"incomplete={len(r1['incomplete'])} mismatch={len(r1['mismatches'])} "
        f"suspect={len(r1['suspects'])}"
    )
    if r2:
        print(
            f"Method2: have={r2['n']} expect={len(m2_items)} "
            f"missing={len(r2['missing'])} extra={len(r2['extra'])} "
            f"incomplete={len(r2['incomplete'])} mismatch={len(r2['mismatches'])} "
            f"suspect={len(r2['suspects'])}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
