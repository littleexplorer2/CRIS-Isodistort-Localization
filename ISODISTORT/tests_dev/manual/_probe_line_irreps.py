from isocore.data.irreps_cdml import _ml_kov_tables

ml, kov = _ml_kov_tables()
for prefix in ["SM", "DT", "LD", "Q", "W", "Y", "GP", "B", "C", "A"]:
    hits = sorted({n for n in ml if n.startswith(prefix) and n[len(prefix):len(prefix)+1].isdigit()})
    print(prefix, hits[:20], "n=", len(hits))

# Official NdNiO2 Y IR tags
from pathlib import Path
root = Path(r"output_compare/NdNiO2 own.cif/官网/Method2")
seen = {}
for p in root.rglob("subgroup.cif"):
    t = p.read_text(encoding="utf-8", errors="replace")
    ir = next((ln for ln in t.splitlines() if ln.startswith("# IR:")), "")
    # IR: Y1, k6t1
    m = __import__("re").search(r"# IR:\s*(\S+),\s*(\S+)", ir)
    if m:
        seen[m.group(1)] = m.group(2)
print("NdNiO2 Method2 IR tags", seen)
