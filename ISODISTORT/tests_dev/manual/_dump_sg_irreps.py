from isocore.data.irreps_cdml import _ISO_DATA
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
from isocore.data.irreps_cdml import _ml_kov_tables
ml, kov = _ml_kov_tables()

for sg in (123, 139):
    start = ptrs[sg - 1] - 1
    end = ptrs[sg] - 1
    print(f"\n=== SG {sg} [{start}:{end}] n={end-start} ===")
    for i in range(start, end):
        print(f"  {ml[i]:10s} {kov[i]}")
