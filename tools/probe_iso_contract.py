"""Frozen, independent-input probes of the official ISO command contract.

Research only. This script calls the production isolation boundary and never
reads any website calculation result or uses expected answers as input.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cris_rea.application import RunService
from cris_rea.backend.iso_protocol import (
    build_direction_query, build_irrep_query, build_isotropy_query, build_parent_query,
)

CASES = {
    "parent-221-wyckoff-parameters": """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 221
SHOW PARENT
SHOW WYCKOFF VECTOR
DISPLAY PARENT
QUIT
""",
    "frequency-221-x1plus-p1": """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 221
VALUE IRREP X1+
VALUE DIRECTION P1
SHOW PARENT
SHOW IRREP
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW BASIS
SHOW ORIGIN
SHOW SIZE
SHOW INDEX
SHOW FREQUENCY DIR
SHOW ELEMENTS
DISPLAY ISOTROPY
QUIT
""",
    "frequency-221-gm4minus-p1": """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
VALUE PARENT 221
VALUE IRREP GM4-
VALUE DIRECTION P1
SHOW PARENT
SHOW IRREP
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW BASIS
SHOW ORIGIN
SHOW FREQUENCY DIR
DISPLAY ISOTROPY
QUIT
""",
    "direction-221-to99": """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
VALUE PARENT 221
VALUE SUBGROUP 99
VALUE BASIS 0,1,0 0,0,1 1,0,0
VALUE ORIGIN 0,0,0
SHOW PARENT
SHOW IRREP
SHOW KPOINT
SHOW DIRECTION VECTOR
DISPLAY DIRECTION
QUIT
""",
    "elements-221-gm4minus-p1": """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 221
VALUE IRREP GM4-
VALUE DIRECTION P1
SHOW PARENT
SHOW IRREP
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW BASIS
SHOW ORIGIN
SHOW SIZE
SHOW INDEX
SHOW ELEMENTS
DISPLAY ISOTROPY
QUIT
""",
    "dt-221-quarter-subgroups": """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 221
VALUE IRREP DT1
VALUE KVALUE 1,1/4
SHOW PARENT
SHOW IRREP
SHOW KPOINT
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW BASIS
SHOW ORIGIN
SHOW SIZE
SHOW INDEX
SHOW ELEMENTS
DISPLAY ISOTROPY
Y
QUIT
""",
    "gp-1-rational-subgroups": """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 1
VALUE IRREP GP1GQ1
VALUE KVALUE 3,1/3,1/4,0
SHOW PARENT
SHOW IRREP
SHOW KPOINT
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW BASIS
SHOW ORIGIN
SHOW SIZE
SHOW INDEX
SHOW ELEMENTS
DISPLAY ISOTROPY
Y
QUIT
""",
    "discover-221-dt": """PAGE NOBREAK
VALUE PARENT 221
VALUE KPOINT DT
SHOW PARENT
SHOW IRREP
SHOW KPOINT
SHOW DIMENSION
DISPLAY IRREP
QUIT
""",
    "discover-1-gp": """PAGE NOBREAK
VALUE PARENT 1
VALUE KPOINT GP
SHOW PARENT
SHOW IRREP
SHOW KPOINT
SHOW DIMENSION
DISPLAY IRREP
QUIT
""",
    "gamma-1-free-sites-cartesian": """PAGE NOBREAK
LABEL VECTOR CONVENTIONAL
VALUE PARENT 1
VALUE IRREP GM1
VALUE DIRECTION P1
VALUE LATTICE PARAMETER 4 5 6 80 85 75
VALUE CELL 1,0,0 0,1,0 0,0,1
SHOW PARENT
SHOW IRREP
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW BASIS
SHOW ORIGIN
SHOW SIZE
SHOW INDEX
DISPLAY ISOTROPY
VALUE WYCKOFF a
VALUE WYCKOFF XYZ 0.123 0.234 0.345
SHOW WYCKOFF
SHOW MICROSCOPIC VECTOR
SHOW CARTESIAN
DISPLAY DISTORTION
VALUE WYCKOFF XYZ 0.267 0.389 0.491
DISPLAY DISTORTION
QUIT
""",
    "gamma-221-multisite": """PAGE NOBREAK
LABEL VECTOR CONVENTIONAL
VALUE PARENT 221
VALUE IRREP GM4-
VALUE DIRECTION P1
VALUE LATTICE PARAMETER 4 4 4 90 90 90
VALUE CELL 1,0,0 0,1,0 0,0,1
SHOW PARENT
SHOW IRREP
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW BASIS
SHOW ORIGIN
SHOW SIZE
SHOW INDEX
DISPLAY ISOTROPY
VALUE WYCKOFF a b c
SHOW WYCKOFF
SHOW MICROSCOPIC VECTOR
SHOW CARTESIAN
DISPLAY DISTORTION
QUIT
""",
}
CASES["dt-221-quarter-generate"] = CASES["dt-221-quarter-subgroups"].replace(
    "DISPLAY ISOTROPY\nY\n", "DISPLAY ISOTROPY\n\n"
)
CASES["gp-1-rational-generate"] = CASES["gp-1-rational-subgroups"].replace(
    "DISPLAY ISOTROPY\nY\n", "DISPLAY ISOTROPY\n\n"
)
CASES["frequency-dt-221-quarter-p1"] = CASES["dt-221-quarter-generate"].replace(
    "SHOW PARENT\n", "VALUE DIRECTION P1\nSHOW FREQUENCY DIR\nSHOW PARENT\n"
).replace("DISPLAY ISOTROPY\n\nQUIT", "DISPLAY ISOTROPY\n" + "\n" * 12 + "QUIT")
CASES["frequency-gp-1-rational-c1"] = CASES["gp-1-rational-generate"].replace(
    "SHOW PARENT\n", "VALUE DIRECTION C1\nSHOW FREQUENCY DIR\nSHOW PARENT\n"
).replace("DISPLAY ISOTROPY\n\nQUIT", "DISPLAY ISOTROPY\n" + "\n" * 12 + "QUIT")
CASES["secondary-221-gm3plus-domain3-general"] = """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
VALUE PARENT 221
VALUE IRREP GM3+
VALUE DIRECTION P1
VALUE DOMAIN 3
VALUE WYCKOFF n
VALUE WYCKOFF XYZ 0.123 0.234 0.345
VALUE LATTICE PARAMETER 4 4 4 90 90 90
VALUE CELL 1,0,0 0,1,0 0,0,1
SHOW PARENT
SHOW IRREP
SHOW DIRECTION
SHOW WYCKOFF
SHOW MICROSCOPIC VECTOR
DISPLAY DISTORTION
QUIT
"""
CASES["gp-1-rational-two-sites-modes"] = """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
VALUE PARENT 1
VALUE IRREP GP1GQ1
VALUE KVALUE 3,1/3,1/4,0
VALUE DIRECTION C1
VALUE CELL 0,0,1 0,4,0 -3,0,0
VALUE LATTICE PARAMETER 4 5 6 80 85 75
VALUE WYCKOFF a
VALUE WYCKOFF XYZ 0.123 0.234 0.345
SHOW PARENT
SHOW IRREP
SHOW DIRECTION
SHOW WYCKOFF
SHOW MICROSCOPIC VECTOR
DISPLAY DISTORTION

VALUE WYCKOFF XYZ 0.267 0.389 0.491
DISPLAY DISTORTION
QUIT
"""
CASES["parent-225-primitive-basis-elements"] = """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 225
SHOW PARENT
SHOW BASIS
SHOW ELEMENTS
DISPLAY PARENT
QUIT
"""
CASES["subgroup-221-gm4minus-p3-centered"] = CASES["elements-221-gm4minus-p1"].replace(
    "VALUE DIRECTION P1", "VALUE DIRECTION P3"
)
CASES["parent-160-primitive-basis-elements"] = CASES[
    "parent-225-primitive-basis-elements"
].replace("VALUE PARENT 225", "VALUE PARENT 160")
CASES["direction-221-to123-quarter"] = CASES["direction-221-to99"].replace(
    "VALUE SUBGROUP 99", "VALUE SUBGROUP 123"
).replace("VALUE BASIS 0,1,0 0,0,1 1,0,0", "VALUE BASIS 1,0,0 0,0,1 0,-4,0"
).replace("DISPLAY DIRECTION\nQUIT", "DISPLAY DIRECTION\n" + "\n" * 12 + "QUIT")
CASES["direction-1-to1-twelve-cell"] = CASES["direction-221-to99"].replace(
    "VALUE PARENT 221", "VALUE PARENT 1"
).replace("VALUE SUBGROUP 99", "VALUE SUBGROUP 1"
).replace("VALUE BASIS 0,1,0 0,0,1 1,0,0", "VALUE BASIS 0,0,1 0,4,0 -3,0,0"
).replace("DISPLAY DIRECTION\nQUIT", "DISPLAY DIRECTION\n" + "\n" * 12 + "QUIT")
CASES["gamma-221-forbidden-scalar-irrep"] = """PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
VALUE PARENT 221
VALUE IRREP GM1+
VALUE DIRECTION P1
VALUE WYCKOFF a
SHOW PARENT
SHOW IRREP
SHOW DIRECTION
SHOW WYCKOFF
SHOW MICROSCOPIC VECTOR
DISPLAY DISTORTION
QUIT
"""
CASES["gamma-1-cartesian-vs-fractional"] = CASES["gamma-1-free-sites-cartesian"].replace(
    "SHOW CARTESIAN\nDISPLAY DISTORTION", "DISPLAY DISTORTION\nSHOW CARTESIAN\nDISPLAY DISTORTION"
)
CASES["negative-invalid-cell-gp"] = CASES["gp-1-rational-two-sites-modes"].replace(
    "VALUE CELL 0,0,1 0,4,0 -3,0,0", "VALUE CELL 1,0,0 0,1,0 0,0,1"
).replace("VALUE WYCKOFF XYZ 0.267 0.389 0.491\nDISPLAY DISTORTION\n", "")
CASES["builder-parent-225"] = build_parent_query(225).decode()
CASES["builder-primary-gamma221"] = build_isotropy_query(221, "GM4-", "P1").decode()
CASES["builder-primary-gp-k"] = build_isotropy_query(
    1, "GP1GQ1", "C1", ["1/3", "1/4", 0]
).decode()
CASES["builder-decline-gp-k"] = build_isotropy_query(
    1, "GP1GQ1", "C1", ["1/3", "1/4", 0], allow_generate=False
).decode()
CASES["builder-direction-twelve-cell"] = build_direction_query(
    1, 1, ((0, 0, -3), (0, 4, 0), (1, 0, 0)), (0, 0, 0)
).decode()
CASES["screen-1000-primary"] = build_isotropy_query(221, "GM4-", "P1").decode().replace(
    "SCREEN 240", "SCREEN 1000"
)
for builder_case in ("builder-parent-225", "builder-primary-gamma221",
                     "builder-primary-gp-k", "builder-decline-gp-k",
                     "builder-direction-twelve-cell", "screen-1000-primary"):
    CASES[builder_case + "-default-setting"] = CASES[builder_case]
CASES["builder-parent-225-wide"] = build_parent_query(225).decode()
CASES["builder-dt-quarter-p2-wide"] = build_isotropy_query(
    221, "DT1", "P2", ["1/4"]
).decode()
for metadata_parent, metadata_irrep in (
    (221, "GM4-"), (221, "DT1"), (221, "X1+"), (1, "GP1GQ1"), (225, "LD1"),
    (221, "C1")
):
    CASES[f"irrep-metadata-{metadata_parent}-{metadata_irrep.lower()}"] = build_irrep_query(
        metadata_parent, metadata_irrep
    ).decode()
CASES["discover-221-all-kpoints"] = """PAGE NOBREAK
SCREEN 1000
LABEL VECTOR CONVENTIONAL
VALUE PARENT 221
SHOW PARENT
SHOW KPOINT
SHOW KDEGREE
DISPLAY KPOINT
QUIT
"""


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    names = sys.argv[1:] or list(CASES)
    if any(name not in CASES for name in names):
        raise SystemExit("Unknown case; choose: " + ", ".join(CASES))
    base = ROOT / "evidence" / "phase2" / "iso-contract"
    base.mkdir(parents=True, exist_ok=True)
    # Freeze every input before the first run. Existing cases are immutable.
    for name in names:
        target = base / name
        target.mkdir(exist_ok=False)
        input_path = target / "input.txt"
        input_path.write_bytes(CASES[name].encode("ascii"))
        (target / "input-provenance.json").write_text(json.dumps({
            "case": name,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "source": "Independent synthetic crystallographic inputs; official ISO manual syntax",
            "manual": "evidence/database/official/isotropy_doc.php.html",
            "manual_sha256": sha(ROOT / "evidence/database/official/isotropy_doc.php.html"),
            "probe_script_sha256": sha(Path(__file__)),
            "input_sha256": sha(input_path),
            "website_case_answers_read": False,
            "acceptance": "Research probe; no general scientific validation implied",
        }, indent=2) + "\n", encoding="utf-8")
    service = RunService(ROOT)
    for name in names:
        target = base / name
        result = service.run("iso", (target / "input.txt").read_bytes(), timeout=60)
        run_dir = Path(result["run_directory"])
        for leaf in ("evidence.json", "stdout.txt", "stderr.txt"):
            shutil.copyfile(run_dir / leaf, target / leaf)
        for path in (run_dir / "work").glob("*"):
            if path.is_file() and not path.is_symlink() and path.stat().st_size < 1024 * 1024:
                (target / "generated-work").mkdir(exist_ok=True)
                shutil.copyfile(path, target / "generated-work" / path.name)
        (target / "result.json").write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        print(json.dumps({
            "case": name, "run_directory": str(run_dir),
            "status": result["evidence"]["status"],
            "stdout": result["stdout"][:20000], "stderr": result["stderr"],
        }), flush=True)


if __name__ == "__main__":
    main()
