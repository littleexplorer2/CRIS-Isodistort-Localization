"""Independent P3 special-position serialization counterexample, no reference answers."""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace as NS

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from cris_rea.backend.structure import (  # noqa: E402
    Lattice, AtomicSite, CrystalStructure, SymmetryOperation, parse_cif, expand_structure,
)
from cris_rea.backend.modes.core import (  # noqa: E402
    construct_modes, distort, close_group, rational_matrix, IDENTITY,
)
from cris_rea.backend.modes.export import cif_text, checked_cif  # noqa: E402

identity = rational_matrix(IDENTITY)
zero = (F(0),) * 3
rotation = SymmetryOperation(((0, -1, 0), (1, -1, 0), (0, 0, 1)), zero)
group = close_group([rotation])
structure = CrystalStructure(
    "threefold", Lattice(4, 4, 6, 90, 90, 120),
    (AtomicSite("Si1", "Si", (F(1, 3), F(2, 3), F(1, 7))),), group, 143,
)
parent = NS(operations=group, primitive_basis=identity)
subgroup = NS(
    operations=group, basis=identity, origin=zero, index=1, size=1,
    subgroup_number=143, subgroup_symbol="P3",
)
model = construct_modes(structure, parent, subgroup, parent)
checked = "--checked" in sys.argv[1:]
if checked:
    output, serialization = checked_cif(model, distort(model), subgroup, "roundtrip")
else:
    output, serialization = cif_text(model, distort(model), subgroup, "roundtrip"), None
result = {
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "role": "Algebraic review counterexample; not a website case or full regression",
    "source_sha256": {
        path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        for path in (
            "cris_rea/backend/modes/core.py", "cris_rea/backend/modes/export.py",
            "cris_rea/backend/structure/model.py",
        )
    },
    "input_definition": {
        "group_generator": "-y,x-y,z", "position": ["1/3", "2/3", "1/7"],
        "cell": [4, 4, 6, 90, 90, 120], "all_mode_amplitudes": 0,
    },
    "pre_export_atom_count": len(model["atoms"]),
    "exported_cif": output,
    "checked_serialization": serialization,
}
try:
    reimported = parse_cif(output)
    result["reimport_atom_count"] = len(expand_structure(reimported))
    result["status"] = "roundtrip_succeeded"
except Exception as error:
    result["status"] = "roundtrip_rejected"
    result["error"] = {
        "type": type(error).__name__, "code": getattr(error, "code", None),
        "message": str(error), "details": getattr(error, "details", None),
    }
destinations = [argument for argument in sys.argv[1:] if argument != "--checked"]
path = Path(destinations[0]) if destinations else Path(__file__).with_suffix(".json")
if path.exists():
    raise SystemExit(f"Refusing to overwrite retained evidence: {path}")
path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(result["status"], path)
