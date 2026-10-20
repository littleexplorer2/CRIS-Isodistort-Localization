"""Engineering tests only: fake resources never represent scientific evidence."""

import contextlib
import hashlib
import http.client
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

from cris_rea.application import RunService
from cris_rea.backend import runtime
from cris_rea.backend import build_iso_query
from cris_rea.backend.worker import supervise
from cris_rea.frontend.cli import main
from cris_rea.frontend.web import create_server


class InputBoundaryTests(unittest.TestCase):
    def test_rejects_unknown_and_path_programs(self):
        for program in ("../iso", "/bin/sh", "findsym; echo hi", "", None):
            with self.subTest(program=program), self.assertRaises(runtime.InputError):
                runtime.validate_request(program, b"input\n", 1)

    def test_rejects_invalid_bytes_and_timeouts(self):
        for data in (b"", b"\0", b"x" * (runtime.MAX_INPUT_BYTES + 1), "text"):
            with self.subTest(data_type=type(data)), self.assertRaises(runtime.InputError):
                runtime.validate_request("iso", data, 1)
        for timeout in (0, -1, 301, float("nan"), float("inf"), True, "30"):
            with self.subTest(timeout=timeout), self.assertRaises(runtime.InputError):
                runtime.validate_request("iso", b"QUIT\n", timeout)

    def test_rejects_run_traversal(self):
        with tempfile.TemporaryDirectory() as directory:
            service = RunService(directory)
            for value in ("../../README.md", "../", "/tmp", "0123", None):
                with self.subTest(value=value), self.assertRaises(ValueError):
                    service.result(value)

    def test_query_uses_user_labels_and_rejects_injection(self):
        commands = build_iso_query(123, "GM1+", "P1", "b").decode("ascii")
        self.assertIn("VALUE PARENT 123\nVALUE IRREP GM1+\nVALUE DIRECTION P1\n", commands)
        self.assertIn("VALUE WYCKOFF b\n", commands)
        self.assertTrue(commands.endswith("QUIT\n"))
        for args in (
            (0, "GM1+", "P1", None), (231, "GM1+", "P1", None),
            (221, "GM1+\nQUIT", "P1", None), (221, "GM1+", "P1\nQUIT", None),
            (221, "GM1+", "P1", "../"),
        ):
            with self.subTest(args=args), self.assertRaises(ValueError):
                build_iso_query(*args)

    def test_doctor_requires_vendor_even_when_linux_probe_passes(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(runtime.sys, "platform", "linux"), \
                patch.object(runtime.platform, "platform", return_value="engineering Linux stub"), \
                patch.object(runtime.subprocess, "run", return_value=
                             subprocess.CompletedProcess([], 0, b"Linux", b"")):
            result = runtime.inspect_runtime(directory)
            self.assertFalse(result["execution_supported"])
            self.assertFalse(result["vendor_ready"])
            self.assertIn("iso", result["missing_vendor_files"])


class ExecutionEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        source = self.root / "vendor" / "isobyu"
        source.mkdir(parents=True)
        for name in ("iso", "findsym"):
            (source / name).write_bytes(b"\x7fELFengineering-test-only")
        (source / "const.dat").write_bytes(b"fake constants")
        (source / "data_engineering.txt").write_bytes(b"fake database")
        (source / "findsym_sample.log").write_bytes(b"REFERENCE MUST NOT BE COPIED")
        (self.root / "evidence").mkdir()
        (self.root / "evidence" / "answer.txt").write_text("REFERENCE")
        self.status = "completed_uninterpreted"
        self.stdout = b"engineering stub"
        self.code = 0
        self.seen = []

    def tearDown(self):
        self.temp.cleanup()

    def fake_supervisor(self, command, **kwargs):
        spec = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        self.seen.append(spec)
        Path(spec["stdout_path"]).write_bytes(self.stdout)
        Path(spec["stderr_path"]).write_bytes(b"")
        Path(spec["result_path"]).write_text(json.dumps({
            "status": self.status, "exit_code": self.code,
        }), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, b"", b"")

    def run_stub(self, program="iso", data=b"raw bytes\r\n"):
        # Explicit control-flow stub: no fake binary is actually executed.
        with patch.object(runtime.sys, "platform", "linux"), \
                patch.object(runtime.platform, "platform", return_value="engineering Linux stub"), \
                patch.object(runtime.subprocess, "run", side_effect=self.fake_supervisor):
            return runtime.execute_program(self.root, program, data, 1)

    def test_separate_runs_preserve_bytes_and_vendor(self):
        original = (self.root / "vendor/isobyu/data_engineering.txt").read_bytes()
        first = self.run_stub()
        second = self.run_stub()
        self.assertNotEqual(first["run_id"], second["run_id"])
        for evidence in (first, second):
            folder = self.root / "runs" / evidence["run_id"]
            self.assertEqual((folder / "input.dat").read_bytes(), b"raw bytes\r\n")
            self.assertEqual(
                evidence["input"]["sha256"], hashlib.sha256(b"raw bytes\r\n").hexdigest()
            )
            self.assertEqual(
                set(evidence["vendor_resources_before"]),
                {"iso", "const.dat", "data_engineering.txt"},
            )
            self.assertFalse((folder / "isodata/findsym_sample.log").exists())
            self.assertEqual(evidence["scientific_validation"], "not_performed")
            self.assertIn("stdout.txt", evidence["artifacts"])
            (folder / "isodata/data_engineering.txt").write_bytes(b"modified copy")
        self.assertEqual((self.root / "vendor/isobyu/data_engineering.txt").read_bytes(), original)

    def test_file_argument_and_isodata_remain_short(self):
        self.run_stub("findsym")
        self.assertEqual(self.seen[-1]["command"], ["../isodata/findsym", "../input.dat"])
        self.assertEqual(self.seen[-1]["isodata"], "../isodata/")
        self.assertEqual(self.seen[-1]["input_mode"], "file")
        self.assertEqual(Path(self.seen[-1]["cwd"]).name, "work")

    def test_timeout_and_failure_never_become_completion(self):
        for status, code in (("timed_out", -9), ("process_failed", 2), ("output_limit", -9)):
            self.status, self.code = status, code
            evidence = self.run_stub()
            self.assertEqual(evidence["status"], status)
            self.assertEqual(evidence["process"]["exit_code"], code)

    def test_zero_exit_syntax_error_is_diagnostic(self):
        self.stdout = b"Isotropy\n*Syntax error: BAD_COMMAND\n*"
        result = self.run_stub()
        self.assertEqual(result["status"], "program_diagnostic")
        self.assertEqual(result["process"]["exit_code"], 0)
        self.assertEqual(result["diagnostics"][0]["code"], "syntax_error")

    def test_consecutive_iso_prompts_do_not_hide_syntax_error(self):
        self.stdout = b"Isotropy\n***Syntax error: BAD_COMMAND\n*"
        self.assertEqual(self.run_stub()["status"], "program_diagnostic")

    def test_zero_exit_explicit_converter_error_is_diagnostic(self):
        self.stdout = b"!useKeyWords\nError reading symmetry operators from CIF file.\n" \
            b"This program has bombed: read_cifv 93\n"
        result = self.run_stub()
        self.assertEqual(result["status"], "program_diagnostic")
        self.assertEqual(result["process"]["exit_code"], 0)
        self.assertEqual(
            [entry["code"] for entry in result["diagnostics"]],
            ["explicit_error", "program_bombed"],
        )

    def test_missing_vendor_records_runtime_error(self):
        (self.root / "vendor/isobyu/iso").unlink()
        result = self.run_stub()
        self.assertEqual(result["status"], "runtime_error")
        self.assertEqual(len(self.seen), 0)
        self.assertTrue((self.root / "runs" / result["run_id"] / "evidence.json").is_file())

    def test_unavailable_supervisor_does_not_make_output(self):
        with patch.object(runtime.sys, "platform", "linux"), \
                patch.object(runtime.platform, "platform", return_value="engineering Linux stub"), \
                patch.object(runtime.subprocess, "run", side_effect=FileNotFoundError("no runtime")):
            result = runtime.execute_program(self.root, "iso", b"QUIT\n", 1)
        self.assertEqual(result["status"], "runtime_error")
        self.assertNotIn("stdout.txt", result["artifacts"])


class CliTests(unittest.TestCase):
    def test_doctor_unavailable_exit_code(self):
        with patch("cris_rea.frontend.cli.RunService.doctor", return_value={
            "execution_supported": False, "scientific_validation": "not_performed"
        }), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(["doctor"]), 2)

    def test_cli_preserves_bytes_and_diagnostic_exit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.txt"
            path.write_bytes(b"raw\r\n")
            with patch("cris_rea.frontend.cli.RunService.run", return_value={
                "evidence": {"status": "program_diagnostic"}
            }) as run, contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["run", "iso", "--input", str(path)]), 1)
                self.assertEqual(run.call_args.args[:2], ("iso", b"raw\r\n"))

    def test_query_preview_does_not_execute(self):
        stream = io.StringIO()
        with patch("cris_rea.frontend.cli.RunService.run") as run, \
                contextlib.redirect_stdout(stream):
            self.assertEqual(main([
                "query", "--parent", "221", "--irrep", "GM4-", "--opd", "P1",
            ]), 0)
            self.assertIn("VALUE PARENT 221", stream.getvalue())
            run.assert_not_called()


class HttpTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.service = RunService(self.temp.name)
        self.server = create_server(self.service, 0)
        self.thread = threading.Thread(target=self.server.serve_forever)
        self.thread.start()
        self.port = self.server.server_port
        self.client = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        self.client.request("GET", "/")
        response = self.client.getresponse()
        self.page = response.read().decode()
        self.assertEqual(response.status, 200)
        self.token = re.search(r'const token="([0-9a-f]+)"', self.page).group(1)

    def tearDown(self):
        self.client.close()
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temp.cleanup()

    def request(self, method, path, body=None, headers=None):
        self.client.request(method, path, body=body, headers=headers or {})
        response = self.client.getresponse()
        raw = response.read()
        return response.status, json.loads(raw)

    def test_local_host_and_token_are_required(self):
        status, _ = self.request("GET", "/api/programs")
        self.assertEqual(status, 403)
        status, _ = self.request("GET", "/", headers={"Host": "attacker.example"})
        self.assertEqual(status, 403)
        status, result = self.request("GET", "/api/programs", headers={"X-CRIS-Token": self.token})
        self.assertEqual(status, 200)
        self.assertIn("iso", result)
        self.assertIn("Parent structure + k / IR / OPD", self.page)

    def test_cross_origin_run_rejected_before_service(self):
        with patch.object(self.service, "run") as run:
            status, _ = self.request("POST", "/api/run", "{}", {
                "X-CRIS-Token": self.token, "Origin": "https://attacker.example",
                "Content-Type": "application/json",
            })
            self.assertEqual(status, 403)
            run.assert_not_called()

    def test_file_bytes_are_passed_and_unknown_program_rejected(self):
        headers = {"X-CRIS-Token": self.token, "Content-Type": "application/json"}
        with patch.object(self.service, "run", return_value={"evidence": {"status": "timed_out"}}) as run:
            status, result = self.request("POST", "/api/run", json.dumps({
                "program": "iso", "input_b64": "cmF3DQo=", "timeout": 0.2,
            }), headers)
            self.assertEqual(status, 200)
            self.assertEqual(result["evidence"]["status"], "timed_out")
            self.assertEqual(run.call_args.args, ("iso", b"raw\r\n", 0.2))
        status, _ = self.request("POST", "/api/run", json.dumps({
            "program": "../../bin/sh", "input": "bad",
        }), headers)
        self.assertEqual(status, 400)
        self.assertFalse((Path(self.temp.name) / "runs").exists())

    def test_rejects_oversized_body_and_invalid_run_id(self):
        self.client.putrequest("POST", "/api/run")
        self.client.putheader("X-CRIS-Token", self.token)
        self.client.putheader("Content-Type", "application/json")
        self.client.putheader("Content-Length", "2000000")
        self.client.endheaders()
        response = self.client.getresponse()
        response.read()
        self.assertEqual(response.status, 413)
        status, _ = self.request("GET", "/api/runs/../../README.md", headers={
            "X-CRIS-Token": self.token
        })
        self.assertEqual(status, 404)

    def test_query_preparation_is_separate_from_execution(self):
        with patch.object(self.service, "run") as run:
            status, result = self.request("POST", "/api/iso-query", json.dumps({
                "parent": 221, "irrep": "GM4-", "opd": "P1", "wyckoff": "a",
            }), {"X-CRIS-Token": self.token, "Content-Type": "application/json"})
            self.assertEqual(status, 200)
            self.assertIn("VALUE WYCKOFF a\n", result["input"])
            run.assert_not_called()
        for invalid in (0, False, []):
            status, _ = self.request("POST", "/api/iso-query", json.dumps({
                "parent": 221, "irrep": "GM4-", "opd": "P1", "wyckoff": invalid,
            }), {"X-CRIS-Token": self.token, "Content-Type": "application/json"})
            self.assertEqual(status, 400)

    def test_partial_request_body_has_a_finite_read_timeout(self):
        self.server.RequestHandlerClass.read_timeout = 0.1
        partial = http.client.HTTPConnection("127.0.0.1", self.port, timeout=2)
        try:
            partial.putrequest("POST", "/api/run")
            partial.putheader("X-CRIS-Token", self.token)
            partial.putheader("Content-Type", "application/json")
            partial.putheader("Content-Length", "100")
            partial.endheaders()
            partial.send(b"{")
            response = partial.getresponse()
            response.read()
            self.assertEqual(response.status, 408)
        finally:
            partial.close()


@unittest.skipUnless(sys.platform == "linux", "Symlink creation needs Linux/WSL.")
class LinkBoundaryTests(unittest.TestCase):
    def test_runs_link_is_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "project"
            outside = Path(directory) / "outside"
            root.mkdir()
            outside.mkdir()
            (root / "runs").symlink_to(outside, target_is_directory=True)
            with self.assertRaises(ValueError):
                runtime.execute_program(root, "iso", b"QUIT\n")
            self.assertEqual(list(outside.iterdir()), [])

    def test_manifest_link_is_not_read(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run_id = "11111111-1111-4111-8111-111111111111"
            folder = root / "runs" / run_id
            folder.mkdir(parents=True)
            target = root / "outside.json"
            target.write_text('{"private": true}')
            (folder / "evidence.json").symlink_to(target)
            with self.assertRaises(ValueError):
                RunService(root).result(run_id)


class BrowserInputContractTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("node"), "Optional Node runtime is unavailable.")
    def test_file_selection_races(self):
        script = Path(__file__).with_name("test_webpage_races.js")
        result = subprocess.run(
            [shutil.which("node"), str(script)], capture_output=True, timeout=10,
        )
        self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))


@unittest.skipUnless(sys.platform == "linux", "Process groups require Linux/WSL.")
class LinuxSupervisorTests(unittest.TestCase):
    def run_worker(self, code, timeout=2, limit=10000):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder / "input.dat").write_bytes(b"input\r\n")
            spec = {
                "command": [sys.executable, "-c", code], "cwd": directory,
                "isodata": directory, "input_mode": "stdin",
                "input_path": str(folder / "input.dat"),
                "stdout_path": str(folder / "stdout.txt"),
                "stderr_path": str(folder / "stderr.txt"),
                "result_path": str(folder / "process.json"),
                "timeout_seconds": timeout, "max_output_bytes": limit,
            }
            path = folder / "execution.json"
            path.write_text(json.dumps(spec), encoding="utf-8")
            self.assertEqual(supervise(path), 0)
            return json.loads((folder / "process.json").read_text()), \
                (folder / "stdout.txt").read_bytes()

    def test_child_input_and_nonzero_exit(self):
        result, stdout = self.run_worker(
            "import sys; sys.stdout.buffer.write(sys.stdin.buffer.read()); sys.exit(7)"
        )
        self.assertEqual(stdout, b"input\r\n")
        self.assertEqual(result["exit_code"], 7)
        self.assertEqual(result["status"], "process_failed")

    def test_real_timeout(self):
        result, _ = self.run_worker("import time; time.sleep(5)", timeout=0.1)
        self.assertEqual(result["status"], "timed_out")
        self.assertLess(result["duration_seconds"], 2)
        self.assertEqual(result["exit_code"], -9)

    def test_real_output_limit(self):
        result, _ = self.run_worker("print('x'*20000)", limit=1000)
        self.assertEqual(result["status"], "output_limit")


if __name__ == "__main__":
    unittest.main()
