from pathlib import Path
import re

root = Path(r"output_compare/EuAl4 Parent.cif/官网/Method2")
print("exists", root.exists(), "n", len(list(root.iterdir())) if root.exists() else 0)
if root.exists():
    for d in sorted(root.iterdir())[:30]:
        print(d.name)
    for p in sorted(root.rglob("subgroup.cif")):
        t = p.read_text(encoding="utf-8", errors="replace")
        ir = next((ln for ln in t.splitlines() if ln.startswith("# IR:")), None)
        kp = next((ln for ln in t.splitlines() if ln.startswith("# k point:")), None)
        print(p.parent.name[:40], "|", kp, "|", ir)
