from isocore.data.irreps_cdml import lookup_irrep_kovalev

cases = [
    ("GM1+", "k14", "k14t1"),
    ("GM2+", "k14", None),
    ("GM5-", "k14", "k14t10"),
    ("M3-", "k15", "k15t4"),
    ("M4-", "k15", "k15t8"),
    ("M5+", "k15", "k15t9"),
    ("A1+", "k20", "k20t1"),
    ("A2+", "k20", "k20t5"),
    ("A3+", "k20", "k20t3"),
    ("A3-", "k20", "k20t4"),
    ("Y1", "k6", "k6t1"),
    ("Y2", "k6", "k6t2"),
]
for ir, kk, want in cases:
    got = lookup_irrep_kovalev(ir, k_kovalev=kk)
    print(ir, kk, "->", got, "want", want, "OK" if want is None or got == want else "BAD")

# From official EuAl4 Method1 samples in compare
from pathlib import Path
import re
off = Path(r"output_compare/EuAl4 Parent.cif/官网/Method1")
pairs = {}
for p in off.rglob("subgroup.cif"):
    t = p.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"# IR:\s*(\S+),\s*(\S+)", t)
    kp = re.search(r"# k point:\s*([A-Za-z]+),\s*(k\d+)", t)
    if m and kp:
        key = (m.group(1), kp.group(2))
        pairs[key] = m.group(2)

bad = 0
for (ir, kk), want in sorted(pairs.items()):
    got = lookup_irrep_kovalev(ir, k_kovalev=kk)
    if got != want:
        bad += 1
        if bad <= 15:
            print("MISMATCH", ir, kk, "got", got, "want", want)
print("EuAl4 M1 unique IR tags checked", len(pairs), "mismatches without companions", bad)

# NdNiO2
off2 = Path(r"output_compare/NdNiO2 own.cif/官网/Method1")
pairs2 = {}
for p in off2.rglob("subgroup.cif"):
    t = p.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"# IR:\s*(\S+),\s*(\S+)", t)
    kp = re.search(r"# k point:\s*([A-Za-z]+),\s*(k\d+)", t)
    if m and kp:
        pairs2[(m.group(1), kp.group(2))] = m.group(2)
bad2 = 0
for (ir, kk), want in sorted(pairs2.items()):
    got = lookup_irrep_kovalev(ir, k_kovalev=kk)
    if got != want:
        bad2 += 1
        if bad2 <= 15:
            print("Nd MISMATCH", ir, kk, "got", got, "want", want)
print("NdNiO2 M1 unique", len(pairs2), "mismatches", bad2)
