from isocore.data.irreps_cdml import _section_tokens, _ISO_DATA
from pathlib import Path
import re

text = _ISO_DATA.read_text(encoding="latin-1", errors="replace")

def section_ints(start, end):
    lines = text.splitlines()
    s = e = None
    for i, line in enumerate(lines):
        if line.strip() == start:
            s = i + 1
        elif line.strip() == end and s is not None:
            e = i
            break
    blob = "\n".join(lines[s:e])
    return [int(x) for x in re.findall(r"-?\d+", blob)]

sgs = section_ints("irrep_space_group", "irrep_space_group_pointer")
ptrs = section_ints("irrep_space_group_pointer", "irrep_image")
print("len sgs", len(sgs), "len ptrs", len(ptrs), "max ptr", max(ptrs), "n sg unique", len(set(sgs)))

# Check SG 139 range
idxs = [i for i, g in enumerate(sgs) if g == 139]
print("SG139 count", len(idxs), "range", min(idxs), max(idxs))

# pointer for 139: pointers are 1-based starts; SG number i has pointer ptrs[i-1]?
# Actually ptrs has 230 entries for SG 1..230?
print("n pointers", len(ptrs))
# Find start of 139
# If ptrs[k] is 1-based start index for space group (k+1)
start139 = ptrs[138] - 1  # 0-based
end139 = (ptrs[139] - 1) if len(ptrs) > 139 else len(sgs)
print("ptr start139", ptrs[138], "end", ptrs[139] if len(ptrs)>139 else None)
print("slice check sgs in range", set(sgs[start139:end139]))

from isocore.data.irreps_cdml import _ml_kov_tables
ml, kov = _ml_kov_tables()
print("ml len", len(ml), "sgs len", len(sgs))
# In SG139 slice look up GM2+ k14 and M3- k15
for i in range(start139, end139):
    if ml[i] in {"GM2+", "M3-", "LD1", "Y1"} and kov[i].startswith(("k14", "k15", "k10", "k6")):
        print(i, ml[i], kov[i])
