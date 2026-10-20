"""Engineering boundaries; these checks are not scientific acceptance."""
import ast
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from cris_rea.application.locking import project_lock
from cris_rea.application.modes import ModeWorkflow
from cris_rea.application.service import RunService


class ModeBoundaryTests(unittest.TestCase):
    def test_layers_do_not_depend_on_upper_layers(self):
        root = Path(__file__).resolve().parents[1] / "cris_rea"
        for layer, forbidden in (("backend", ("frontend", "application")),
                                 ("application", ("frontend",)),
                                 ("frontend", ("backend",))):
            for path in (root / layer).rglob("*.py"):
                tree = ast.parse(path.read_text(encoding="utf-8"))
                for node in ast.walk(tree):
                    if isinstance(node, ast.ImportFrom):
                        parts = (node.module or "").split(".")
                        self.assertFalse(set(parts) & set(forbidden), str(path))
                    elif isinstance(node, ast.Import):
                        for name in node.names:
                            self.assertFalse(set(name.name.split(".")) & set(forbidden), str(path))

    def test_raw_and_modes_share_same_lock(self):
        with tempfile.TemporaryDirectory() as folder:
            raw, modes = RunService(folder), ModeWorkflow(folder)
            self.assertIs(raw._run_lock, modes._lock)
            self.assertIs(modes._lock, project_lock(Path(folder) / "."))
            with modes._lock:
                with self.assertRaisesRegex(ValueError, "already running"):
                    raw.run("iso", b"QUIT\n")
                with self.assertRaisesRegex(ValueError, "already running"):
                    modes.run(b"cif", "GM1", "P1")

    def test_download_requires_frozen_validated_export(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            wid = "11111111-1111-4111-8111-111111111111"
            directory = root / "runs" / wid
            exports = directory / "exports"; exports.mkdir(parents=True)
            path = exports / "modes.json"; path.write_bytes(b"{}")
            record = {"status": "rejected_or_inconclusive", "export_fingerprints": {
                "modes.json": {"sha256": hashlib.sha256(b"{}").hexdigest(), "bytes": 2}}}
            (directory / "workflow.json").write_text(json.dumps(record), encoding="utf-8")
            service = ModeWorkflow(root)
            with self.assertRaises(ValueError):
                service.artifact(wid, "modes.json")
            record["status"] = "validated_local_subspace"
            (directory / "workflow.json").write_text(json.dumps(record), encoding="utf-8")
            self.assertEqual(service.artifact(wid, "modes.json")[1], b"{}")
            with patch.object(Path, "read_bytes", side_effect=[b"{}", b"unverified"]) as read:
                self.assertEqual(service.artifact(wid, "modes.json")[1], b"{}")
                self.assertEqual(read.call_count, 1)
            path.write_bytes(b'{"changed":true}')
            with self.assertRaisesRegex(ValueError, "fingerprint"):
                service.artifact(wid, "modes.json")
            for invalid in ("../parent.cif", "parent.cif", "workflow.json"):
                with self.assertRaises(ValueError):
                    service.artifact(wid, invalid)
            (directory / "workflow.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "incomplete or invalid"):
                service.artifact(wid, "modes.json")


if __name__ == "__main__":
    unittest.main()
