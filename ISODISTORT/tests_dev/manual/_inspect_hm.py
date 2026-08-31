from pathlib import Path
import re

for label, root in [
    ("official", Path(r"output_compare/NdNiO2 own.cif/官网/Method1")),
    ("local", Path(r"output_compare/NdNiO2 own.cif/现有网页版交互/Method1")),
]:
    p = next(root.rglob("subgroup.cif"))
    t = p.read_text(encoding="utf-8")
    print("===", label, p.parent.name[:60])
    for ln in t.splitlines():
        if "space_group" in ln or "H-M" in ln or ln.startswith("# k") or ln.startswith("# IR"):
            print(repr(ln))

print("\n=== bad HM scan ===")
bad = []
for p in Path("output_compare").rglob("subgroup.cif"):
    if "现有" not in str(p):
        continue
    t = p.read_text(encoding="utf-8", errors="replace")
    m = re.search(r'_space_group_name_H-M_alt\s+"([^"]+)"', t)
    if not m:
        continue
    hm = m.group(1)
    # glued letter+digit that isn't a known screw axis token
    tokens = hm.split()
    for tok in tokens:
        if re.fullmatch(r"[a-zA-Z]+\d+", tok) and tok not in {
            # none expected; screws are like 21, 42 alone
        }:
            bad.append((hm, str(p)[-90:]))
            break
print("bad count", len(bad))
for x in bad[:25]:
    print(x)
