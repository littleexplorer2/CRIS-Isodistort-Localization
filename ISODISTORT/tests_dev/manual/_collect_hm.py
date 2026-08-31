from pathlib import Path
import re
from isocore.io.isodistort_cif import _HM_FULL_1992, _space_hm_symbol

seen = {}
for p in Path("output_compare").rglob("subgroup.cif"):
    if "官网" not in str(p):
        continue
    t = p.read_text(encoding="utf-8", errors="replace")
    m = re.search(r'_symmetry_space_group_name_H-M\s+"([^"]+)"', t)
    n = re.search(r"_symmetry_Int_Tables_number\s+(\d+)", t)
    if m and n:
        num = int(n.group(1))
        if num not in seen:
            seen[num] = m.group(1)

print("official unique HM", len(seen))
missing = []
wrong = []
for num, hm in sorted(seen.items()):
    tab = _HM_FULL_1992.get(num)
    if tab is None:
        missing.append((num, hm))
    elif tab != hm:
        wrong.append((num, tab, hm))
print("missing from table", missing)
print("table differs", wrong[:20], "n=", len(wrong))
