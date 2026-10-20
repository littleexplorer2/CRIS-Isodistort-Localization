"""Targeted real-binary integration of the new finite-cell workflow.

Authored algebraic structures are explicit independent inputs, not experimental
claims or downloaded case answers. Official parent tables are generic inputs.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import traceback
import os

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from cris_rea.backend.iso_protocol import build_parent_query, parse_parent_output
from cris_rea.backend.runtime import execute_program
from cris_rea.backend.modes.core import complete_parent_group
from cris_rea.backend.modes import seitz_expression
from cris_rea.backend.structure import parse_cif, expand_structure
from cris_rea.application.modes import ModeWorkflow

CASES = [
    {"name": "cubic_polar_axis", "parent": 221, "irrep": "GM4-", "opd": "P1", "k": [],
     "cell": [4, 4, 4, 90, 90, 90],
     "sites": [["A", "Ba", "0", "0", "0"], ["B", "Ti", "0.5", "0.5", "0.5"],
               ["O", "O", "0", "0.5", "0.5"]]},
    {"name": "cubic_rhombohedral_centered", "parent": 221, "irrep": "GM4-", "opd": "P3", "k": [],
     "cell": [4, 4, 4, 90, 90, 90],
     "sites": [["A", "Ba", "0", "0", "0"], ["B", "Ti", "0.5", "0.5", "0.5"],
               ["O", "O", "0", "0.5", "0.5"]]},
    {"name": "face_centered_parent", "parent": 225, "irrep": "GM4-", "opd": "P1", "k": [],
     "cell": [5.6, 5.6, 5.6, 90, 90, 90],
     "sites": [["A", "Na", "0", "0", "0"], ["B", "Cl", "0.5", "0.5", "0.5"]]},
    {"name": "quarter_k_supercell", "parent": 221, "irrep": "DT1", "opd": "P1", "k": ["1/4"],
     "cell": [4, 4, 4, 90, 90, 90],
     "sites": [["A", "Ba", "0", "0", "0"], ["B", "Ti", "0.5", "0.5", "0.5"],
               ["O", "O", "0", "0.5", "0.5"]]},
    {"name": "triclinic_general_rational_k", "parent": 1, "irrep": "GP1GQ1", "opd": "C1",
     "k": ["1/3", "1/4", "0"], "cell": [4, 5, 6, 80, 85, 75],
     "sites": [["A", "Si", "0.123", "0.234", "0.345"],
               ["B", "O", "0.287", "0.419", "0.531"]]},
]

_audit_active = False
_audit_reads = set()


def _audit(event, args):
    if _audit_active and event == "open" and isinstance(args[0], (str, bytes, os.PathLike)):
        mode = args[1]
        if mode is None or (isinstance(mode, str) and ("r" in mode or "+" in mode)):
            _audit_reads.add(os.fsdecode(args[0]))


sys.addaudithook(_audit)


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    global _audit_active
    selected = [c for c in CASES if not sys.argv[1:] or c["name"] in sys.argv[1:]]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = ROOT / "evidence" / "phase2" / "workflow" / stamp
    out.mkdir(parents=True)
    # Freeze authored inputs before invoking the native program.
    write(out / "input_recipes.json", selected)
    record = {"input_origin": "Authored mathematical test structures; not experimental data.",
              "input_recipe_sha256": sha(out / "input_recipes.json"),
              "reference_case_reads": [], "cases": [], "started_utc": stamp}
    parent_tables = {}
    for case in selected:
        item = {"name": case["name"]}
        record["cases"].append(item)
        dest = out / case["name"]; dest.mkdir()
        try:
            number = case["parent"]
            if number not in parent_tables:
                ev = execute_program(ROOT, "iso", build_parent_query(number), 60)
                if ev["status"] != "completed_uninterpreted":
                    raise RuntimeError(f"parent database query: {ev['status']}, {ev.get('error')}")
                source = ROOT / "runs" / ev["run_id"]
                pg = parse_parent_output((source / "stdout.txt").read_text(), expected_parent=number)
                parent_tables[number] = (pg, ev["run_id"])
            pg, table_run = parent_tables[number]
            ops = complete_parent_group(pg.operations, pg.primitive_basis)
            lines = [f"data_{case['name']}", f"_space_group_IT_number {number}"]
            for tag, value in zip(("length_a", "length_b", "length_c", "angle_alpha",
                                   "angle_beta", "angle_gamma"), case["cell"]):
                lines.append(f"_cell_{tag} {value}")
            lines += ["loop_", "_space_group_symop_operation_xyz"]
            lines += ["'" + seitz_expression(op) + "'" for op in ops]
            lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol",
                      "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
            lines += [" ".join(site) for site in case["sites"]]
            cif = ("\n".join(lines)+"\n").encode()
            (dest / "parent.cif").write_bytes(cif)
            item["parent_generic_table_run"] = table_run
            item["input_sha256"] = sha(dest / "parent.cif")
            # Use zero state first; its frozen result fixes amplitude length independently.
            _audit_reads.clear()
            _audit_active = True
            try:
                result = ModeWorkflow(ROOT).run(cif, case["irrep"], case["opd"], case["k"], timeout=90)
            finally:
                _audit_active = False
                item["python_reads_during_production"] = sorted(_audit_reads)
            forbidden = [p for p in _audit_reads if str(ROOT / "evidence").lower() in p.lower()
                         or p.lower().startswith("e:\\cris")]
            if forbidden:
                raise AssertionError("Production read reference/source project paths: " + repr(forbidden))
            item["workflow_id"] = result["workflow_id"]
            run = Path(result["directory"])
            shutil.copyfile(run / "workflow.json", dest / "workflow.json")
            for file in run.iterdir():
                if file.suffix in (".in", ".out"):
                    shutil.copyfile(file, dest / file.name)
            shutil.copytree(run / "exports", dest / "exports")
            # Freeze production output before independent read-only checks.
            item["frozen_exports"] = {p.name: sha(p) for p in (dest / "exports").iterdir()}
            data = json.loads((dest / "exports" / "modes.json").read_text(encoding="utf-8"))
            basis = np.array(data["displacement_basis"])
            dimension = len(data["mode_names"])
            # Handles zero-dimensional spaces without pretending an absent output is failure.
            basis = basis.reshape(3*len(data["atoms"]), dimension)
            coeff = np.linspace(-0.023, 0.017, dimension)
            displacement = basis@coeff
            np.testing.assert_allclose(basis.T@displacement, coeff, atol=1e-11, rtol=0)
            np.testing.assert_allclose(displacement@displacement, coeff@coeff, atol=1e-11, rtol=0)
            parsed = parse_cif((dest / "exports" / "reference.cif").read_bytes())
            count = len(expand_structure(parsed))
            if count != len(data["atoms"]):
                raise AssertionError(f"CIF roundtrip atom count {count} != {len(data['atoms'])}")
            item["roundtrip_atom_count"] = count
            item["summary"] = result["summary"]
            item["validation"] = result["validation"]
            item["status"] = "passed_local_contract"
            print(json.dumps({"case": case["name"], "status": item["status"],
                              "atoms": count, "modes": dimension,
                              "subgroup": result["summary"]["subgroup"]}), flush=True)
        except Exception as exc:
            item["status"] = "failed"
            item["error"] = f"{type(exc).__name__}: {exc}"
            (dest / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
            print(json.dumps(item), flush=True)
        write(out / "record.json", record)
    record["finished_utc"] = datetime.now(timezone.utc).isoformat()
    record["status"] = "passed" if all(c["status"] == "passed_local_contract" for c in record["cases"]) else "failed"
    write(out / "record.json", record)
    print(str(out / "record.json"), flush=True)
    return 0 if record["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
