import json
from pathlib import Path
import re
from .locking import project_lock

from ..backend import PROGRAMS, InputError, build_iso_query, execute_program, inspect_runtime

RUN_ID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}")


class RunService:
    def __init__(self, project_root):
        self.root = Path(project_root).resolve()
        self._run_lock = project_lock(self.root)

    def programs(self):
        return {name: dict(details) for name, details in PROGRAMS.items()}

    def doctor(self):
        return inspect_runtime(self.root)

    def prepare_iso_query(self, parent, irrep, opd, wyckoff=None):
        return {
            "program": "iso",
            "input": build_iso_query(parent, irrep, opd, wyckoff).decode("ascii"),
            "scope": "Raw official ISO subgroup query and optional Wyckoff projection; "
                     "no complete CIF structure workflow or normalized mode export.",
        }

    def run(self, program, data, timeout=30):
        if not self._run_lock.acquire(blocking=False):
            raise InputError("Another computation is already running in this server.")
        try:
            evidence = execute_program(self.root, program, data, timeout)
            return self.result(evidence["run_id"])
        finally:
            self._run_lock.release()

    def result(self, run_id):
        if not isinstance(run_id, str) or not RUN_ID.fullmatch(run_id):
            raise InputError("Invalid run identifier.")
        run_dir = self.root / "runs" / run_id
        if run_dir.is_symlink() or not run_dir.resolve().is_relative_to(self.root / "runs"):
            raise InputError("Run path is outside this project's runs directory.")
        manifest = run_dir / "evidence.json"
        if not manifest.is_file() or manifest.is_symlink() or manifest.resolve().parent != run_dir.resolve():
            raise InputError("Run not found.")
        evidence = json.loads(manifest.read_text(encoding="utf-8"))
        result = {"evidence": evidence, "run_directory": str(run_dir)}
        for name in ("stdout", "stderr"):
            path = run_dir / (name + ".txt")
            # The UI reads a bounded preview; the full bytes stay in the run.
            if path.is_file() and not path.is_symlink():
                with path.open("rb") as stream:
                    raw = stream.read(1024 * 1024 + 1)
                result[name] = raw[:1024 * 1024].decode("utf-8", errors="replace")
                result[name + "_preview_truncated"] = len(raw) > 1024 * 1024
            else:
                result[name] = ""
                result[name + "_preview_truncated"] = False
        return result
