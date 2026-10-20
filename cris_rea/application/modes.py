"""CIF + user-selected IR/OPD to independently checked finite-cell mode exports."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import uuid

import numpy as np

from ..backend.runtime import InputError, execute_program, MAX_INPUT_BYTES
from ..backend.structure import parse_cif
from ..backend.iso_protocol import (
    build_parent_query, build_isotropy_query, build_direction_query,
    parse_parent_output, parse_isotropy_output, build_irrep_query, parse_irrep_output,
)
from ..backend.modes import construct_modes, distort, ModeError
from ..backend.modes.export import export_modes, write_json
from .locking import project_lock

FILES = {"modes.json", "modes.csv", "reference.cif", "distorted.cif", "report.json"}
UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}")


def fingerprint(path):
    data = path.read_bytes()
    return {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}


class ModeWorkflow:
    def __init__(self, project_root):
        self.root = Path(project_root).resolve()
        self._lock = project_lock(self.root)

    def run(self, cif_bytes, irrep, opd, k_values=(), amplitudes=None,
            strain_amplitudes=None, timeout=60):
        if not isinstance(cif_bytes, bytes) or not 0 < len(cif_bytes) <= MAX_INPUT_BYTES:
            raise InputError("CIF input must contain 1 to 1048576 bytes.")
        if not self._lock.acquire(blocking=False):
            raise InputError("Another computation is already running in this project.")
        try:
            return self._run(cif_bytes, irrep, opd, k_values, amplitudes, strain_amplitudes, timeout)
        finally:
            self._lock.release()

    def _run(self, cif_bytes, irrep, opd, k_values, amplitudes, strain_amplitudes, timeout):
        runs = self.root / "runs"
        if runs.is_symlink() or runs.resolve() != runs:
            raise InputError("runs must be a real directory inside this project.")
        runs.mkdir(exist_ok=True)
        workflow_id = str(uuid.uuid4())
        directory = runs / workflow_id
        directory.mkdir()
        (directory / "parent.cif").write_bytes(cif_bytes)
        record = {
            "workflow_id": workflow_id, "started_utc": datetime.now(timezone.utc).isoformat(),
            "status": "preparing", "directory": str(directory), "input": fingerprint(directory / "parent.cif"),
            "request": {"irrep": irrep, "opd": opd, "k_values": list(k_values),
                        "amplitudes": amplitudes, "strain_amplitudes": strain_amplitudes},
            "binary_runs": [], "versions": {"python": platform.python_version(), "numpy": np.__version__},
            "external_resources": "Only vendor/isobyu generic programs/databases; no reference-result reads.",
            "source_code": {str(p.relative_to(self.root)): fingerprint(p)
                            for p in sorted((self.root / "cris_rea").rglob("*.py"))}}
        write_json(directory / "workflow.json", record)

        def native(label, commands):
            (directory / f"{label}.in").write_bytes(commands)
            evidence = execute_program(self.root, "iso", commands, timeout)
            rid = evidence["run_id"]
            record["binary_runs"].append({"role": label, "run_id": rid,
                                           "evidence": str(runs / rid / "evidence.json")})
            raw = (runs / rid / "stdout.txt").read_bytes() if (runs / rid / "stdout.txt").exists() else b""
            (directory / f"{label}.out").write_bytes(raw)
            if evidence["status"] != "completed_uninterpreted":
                raise ModeError("native_execution_failed", f"Official ISO {label}: {evidence['status']}",
                                run_id=rid, error=evidence.get("error"))
            return raw.decode("utf-8", errors="replace")

        try:
            structure = parse_cif(cif_bytes)
            if structure.space_group_number is None:
                raise ModeError("missing_parent_number", "CIF must declare the parent space-group number.")
            parent_number = structure.space_group_number
            metadata = parse_irrep_output(native("irrep_metadata", build_irrep_query(parent_number, irrep)),
                                          expected_parent=parent_number, expected_irrep=irrep)
            nominal_k = metadata.k_vector(k_values)
            record["irrep_metadata"] = {
                "wavevector_expression": metadata.wavevector_expression,
                "parameter_symbols": list(metadata.parameter_symbols),
                "dimension": metadata.dimension,
                "nominal_k": [str(v) for v in nominal_k],
                "note": "Nominal reciprocal conventional star representative; active arms depend on OPD. "
                        "No claim of independent full irrep-matrix/phase verification."}
            parent = parse_parent_output(native("parent", build_parent_query(parent_number)),
                                         expected_parent=parent_number)
            primary = parse_isotropy_output(native("subgroup", build_isotropy_query(
                parent_number, irrep, opd, k_parameters=k_values)),
                expected_parent=parent_number, expected_irrep=irrep, expected_opd=opd)
            h_parent = parse_parent_output(native("subgroup_lattice", build_parent_query(
                primary.subgroup_number)), expected_parent=primary.subgroup_number)
            model = construct_modes(structure, parent, primary, h_parent)
            state = distort(model, amplitudes, strain_amplitudes)
            # This independently generated listing is annotation, not a source of mode coefficients.
            secondary = native("compatible_irreps", build_direction_query(
                parent_number, primary.subgroup_number, primary.basis, primary.origin))
            record["compatible_irreps_note"] = (
                "Official DISPLAY DIRECTION raw listing saved; rounded vectors never enter the basis.")
            record["structure_notices"] = list(structure.notices)
            record["summary"] = {
                "parent": parent_number, "subgroup": primary.subgroup_number,
                "subgroup_symbol": primary.subgroup_symbol, "atoms": len(model["atoms"]),
                "displacement_modes": model["basis"].shape[1],
                "rigid_translations": model["rigid_translations"],
                "strain_modes": model["strain_basis"].shape[1],
                "basis_convention": model["basis_convention"],
                "compatible_irreps_raw": secondary}
            export_modes(directory / "exports", model, state, primary, record)
            record["status"] = "validated_local_subspace"
            record["files"] = {name: str(directory / "exports" / name) for name in sorted(FILES)}
            record["export_fingerprints"] = {
                name: fingerprint(directory / "exports" / name) for name in sorted(FILES)}
            record["validation"] = model["validation"]
            return record
        except Exception as exc:
            record["status"] = "rejected_or_inconclusive"
            record["error"] = {"type": type(exc).__name__, "code": getattr(exc, "code", None),
                               "message": str(exc), "details": getattr(exc, "details", {})}
            raise
        finally:
            record["finished_utc"] = datetime.now(timezone.utc).isoformat()
            write_json(directory / "workflow.json", record)

    def artifact(self, workflow_id, filename):
        if not isinstance(workflow_id, str) or not UUID.fullmatch(workflow_id) or filename not in FILES:
            raise InputError("Invalid mode artifact.")
        directory = self.root / "runs" / workflow_id
        path = directory / "exports" / filename
        if (directory.is_symlink() or (directory / "exports").is_symlink() or path.is_symlink()
                or not path.is_file() or not path.resolve().is_relative_to(self.root / "runs")):
            raise InputError("Artifact not found.")
        manifest = directory / "workflow.json"
        if manifest.is_symlink() or not manifest.is_file():
            raise InputError("Workflow record not found.")
        try:
            record = json.loads(manifest.read_text(encoding="utf-8"))
            expected = record["export_fingerprints"][filename]
        except (ValueError, KeyError, TypeError) as exc:
            raise InputError("Workflow record is incomplete or invalid.") from exc
        if record.get("status") != "validated_local_subspace":
            raise InputError("This workflow has no validated export.")
        data = path.read_bytes()
        if {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)} != expected:
            raise InputError("Artifact fingerprint changed after validation.")
        return path, data
