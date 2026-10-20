"""Frontend transport tests use explicit stubs; they are not scientific output evidence."""
import base64
import contextlib
import http.client
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import threading
import unittest
from unittest.mock import Mock, patch

from cris_rea.application import RunService
from cris_rea.frontend.cli import main
from cris_rea.frontend.web import create_server
from cris_rea.frontend.webpage import PAGE

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_ID = "0acdabef-1412-4111-8111-0123456789ab"


def result(status="validated_local_subspace"):
    return {"workflow_id": WORKFLOW_ID, "status": status, "directory": "engineering-stub",
            "summary": {"parent": 1, "subgroup": 1, "atoms": 2, "displacement_modes": 3,
                        "rigid_translations": 3, "strain_modes": 6,
                        "basis_convention": "engineering transport stub"},
            "files": {"modes.json": "engineering/modes.json"}}


class ModesCliTests(unittest.TestCase):
    def test_comma_k_parameters_preserve_order_and_reject_empty_entries(self):
        (ROOT / "runs/frontend-test-temp").mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT / "runs/frontend-test-temp") as directory:
            path = Path(directory) / "parent.cif"
            path.write_bytes(b"data_independent\n")
            for tokens, expected in (
                    (["--k-values=1/3,-1/4,0"], ["1/3", "-1/4", "0"]),
                    (["--k-values=-1/3,-1/4,0"], ["-1/3", "-1/4", "0"]),
                    (["--k-values", "1/3", "1/4", "0"], ["1/3", "1/4", "0"])):
                with self.subTest(tokens=tokens):
                    mode = Mock()
                    mode.run.return_value = result()
                    with patch("cris_rea.frontend.cli._mode_workflow", return_value=mode), \
                            contextlib.redirect_stdout(io.StringIO()):
                        code = main(["modes", "--cif", str(path), "--irrep", "GP1",
                                     "--opd", "P1", *tokens])
                    self.assertEqual(code, 0)
                    self.assertEqual(mode.run.call_args.kwargs["k_values"], expected)
            for token in ("--k-values=1/3,,0", "--k-values=,1/3", "--k-values=1/3,", "--k-values="):
                with self.subTest(token=token), \
                        patch("cris_rea.frontend.cli._mode_workflow") as factory, \
                        contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(main(["modes", "--cif", str(path), "--irrep", "GP1",
                                           "--opd", "P1", token]), 2)
                    factory.assert_not_called()

    def test_mode_bytes_parameters_and_status(self):
        (ROOT / "runs/frontend-test-temp").mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT / "runs/frontend-test-temp") as directory:
            path = Path(directory) / "parent.cif"
            path.write_bytes(b"data_original\r\n# unchanged\r\n")
            mode = Mock()
            mode.run.return_value = result()
            stream = io.StringIO()
            with patch("cris_rea.frontend.cli._mode_workflow", return_value=mode), contextlib.redirect_stdout(stream):
                code = main(["modes", "--cif", str(path), "--irrep", "LD1", "--opd", "P1",
                             "--k-values", "1/3", "--amplitudes", ".1", "-.2",
                             "--strain-amplitudes", ".01", "--timeout", "15"])
            self.assertEqual(code, 0)
            self.assertEqual(json.loads(stream.getvalue())["status"], "validated_local_subspace")
            mode.run.assert_called_once_with(
                b"data_original\r\n# unchanged\r\n", "LD1", "P1",
                k_values=["1/3"], amplitudes=[.1, -.2], strain_amplitudes=[.01], timeout=15)
            mode.run.return_value = result("rejected_or_inconclusive")
            with patch("cris_rea.frontend.cli._mode_workflow", return_value=mode), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["modes", "--cif", str(path), "--irrep", "GM1+", "--opd", "P1"]), 1)

    def test_cli_missing_or_oversize_input_does_not_start(self):
        (ROOT / "runs/frontend-test-temp").mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT / "runs/frontend-test-temp") as directory:
            path = Path(directory) / "large.cif"
            path.write_bytes(b"x" * (1024 * 1024 + 1))
            with patch("cris_rea.frontend.cli._mode_workflow") as factory, contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(main(["modes", "--cif", str(path), "--irrep", "GM1+", "--opd", "P1"]), 2)
                factory.assert_not_called()

    def test_nonfinite_cli_amplitudes_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
            main(["modes", "--cif", "not-read.cif", "--irrep", "GM1+", "--opd", "P1", "--amplitudes", "nan"])
        self.assertEqual(raised.exception.code, 2)


class ModesHttpTests(unittest.TestCase):
    def setUp(self):
        (ROOT / "runs/frontend-test-temp").mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / "runs/frontend-test-temp")
        self.service = RunService(self.temp.name)
        self.mode = Mock()
        self.mode.run.return_value = result()
        self.mode.artifact.return_value = (Path("modes.json"), b'{"transport":"stub"}\r\n')
        self.server = create_server(self.service, 0, mode_workflow=self.mode)
        self.thread = threading.Thread(target=self.server.serve_forever)
        self.thread.start()
        self.client = http.client.HTTPConnection("127.0.0.1", self.server.server_port, timeout=5)
        self.client.request("GET", "/")
        response = self.client.getresponse()
        page = response.read().decode()
        self.token = re.search(r'const token="([0-9a-f]+)"', page).group(1)
        self.headers = {"X-CRIS-Token": self.token, "Content-Type": "application/json"}

    def tearDown(self):
        self.client.close()
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temp.cleanup()

    def send(self, path, payload, headers=None):
        self.client.request("POST", path, json.dumps(payload),
                            headers=self.headers if headers is None else headers)
        response = self.client.getresponse()
        return response.status, json.loads(response.read())

    def test_cif_file_bytes_and_optional_values_reach_workflow(self):
        raw = b"data_frozen\r\n# engineering stub\r\n"
        status, body = self.send("/api/modes", {
            "cif_b64": base64.b64encode(raw).decode(), "irrep": "LD1", "opd": "P1",
            "k_values": ["1/3"], "amplitudes": [.2], "strain_amplitudes": [0], "timeout": 12})
        self.assertEqual(status, 200)
        self.assertEqual(body["status"], "validated_local_subspace")
        self.mode.run.assert_called_once_with(raw, "LD1", "P1", k_values=["1/3"],
                                             amplitudes=[.2], strain_amplitudes=[0], timeout=12)

    def test_text_and_default_values_reach_workflow(self):
        status, _ = self.send("/api/modes", {"cif": "data_text\n", "irrep": "GM1+", "opd": "P1"})
        self.assertEqual(status, 200)
        self.mode.run.assert_called_once_with(b"data_text\n", "GM1+", "P1", k_values=[],
                                             amplitudes=None, strain_amplitudes=None, timeout=60)

    def test_bad_shape_and_nonfinite_requests_never_dispatch(self):
        base = {"cif": "data_text\n", "irrep": "GM1+", "opd": "P1"}
        bad = [
            {"cif_b64": "ZGF0YQ=="}, {"amplitudes": [float("nan")]},
            {"amplitudes": [True]}, {"amplitudes": [10 ** 400]},
            {"strain_amplitudes": "1 2"}, {"k_values": "1/3"},
            {"k_values": [1]}, {"timeout": float("inf")}, {"timeout": 0},
            {"irrep": 1}, {"opd": ""}, {"cif": ""}, {"cif": 1},
        ]
        for changed in bad:
            with self.subTest(changed=changed):
                status, body = self.send("/api/modes", {**base, **changed})
                self.assertEqual(status, 400)
                self.assertIn("error", body)
        self.mode.run.assert_not_called()

    def test_token_origin_and_host_still_apply(self):
        payload = {"cif": "data_text\n", "irrep": "GM1+", "opd": "P1"}
        for headers in ({}, {**self.headers, "Origin": "https://attacker.example"},
                        {**self.headers, "Host": "attacker.example"}):
            self.assertEqual(self.send("/api/modes", payload, headers)[0], 403)
        self.mode.run.assert_not_called()

    def test_application_rejection_or_exception_is_not_success(self):
        payload = {"cif": "data_text\n", "irrep": "GM1+", "opd": "P1"}
        self.mode.run.return_value = result("rejected_or_inconclusive")
        status, value = self.send("/api/modes", payload)
        self.assertEqual(status, 200)
        self.assertEqual(value["status"], "rejected_or_inconclusive")
        for error, expected in ((ValueError("unsupported setting"), 400),
                                (RuntimeError("local failure"), 500)):
            self.mode.run.side_effect = error
            status, body = self.send("/api/modes", payload)
            self.assertEqual(status, expected)
            self.assertIn(str(error), body["error"])

    def test_downloads_are_authenticated_and_exact(self):
        path = f"/api/modes/{WORKFLOW_ID}/modes.json"
        self.client.request("GET", path)
        response = self.client.getresponse()
        self.assertEqual(response.status, 403)
        response.read()
        self.mode.artifact.assert_not_called()
        self.client.request("GET", path, headers={"X-CRIS-Token": self.token})
        response = self.client.getresponse()
        self.assertEqual(response.status, 200)
        self.assertEqual(response.getheader("Content-Disposition"), 'attachment; filename="modes.json"')
        self.assertEqual(response.getheader("Cache-Control"), "no-store")
        self.assertEqual(response.read(), b'{"transport":"stub"}\r\n')
        self.mode.artifact.assert_called_once_with(WORKFLOW_ID, "modes.json")

    def test_traversal_and_unknown_artifacts_rejected(self):
        for suffix in ("../parent.cif", "modes.json/extra", "%2e%2e%2fparent.cif", "workflow.json", "x%0d%0aHeader"):
            self.client.request("GET", f"/api/modes/{WORKFLOW_ID}/{suffix}", headers={"X-CRIS-Token": self.token})
            response = self.client.getresponse()
            self.assertEqual(response.status, 404)
            response.read()
        self.mode.artifact.assert_not_called()
        self.mode.artifact.side_effect = ValueError("bad UUID or changed fingerprint")
        self.client.request("GET", "/api/modes/not-a-uuid/modes.json", headers={"X-CRIS-Token": self.token})
        response = self.client.getresponse()
        self.assertEqual(response.status, 404)
        response.read()


class ModesBrowserContractTests(unittest.TestCase):
    def test_scientific_scope_is_present_and_old_console_preserved(self):
        for text in ("ordinary 3D", "commensurate", "standard setting", "subgroup-fixed",
                     "mix parent irreps", "per-irrep column convention", 'id="program"', 'id="run"',
                     'id="mode-run"', "ISO independent k parameters"):
            self.assertIn(text, PAGE)

    @unittest.skipUnless(shutil.which("node"), "Node is unavailable for DOM adapter checks.")
    def test_mode_file_races_dispatch_and_download_names(self):
        # Executed JS uses the actual page script and a deliberately minimal DOM.
        harness = r"""
const assert=require("node:assert/strict"),vm=require("node:vm");
const source=JSON.parse(process.argv[1]),elements=new Map(),sent=[];
function node(){return {value:"",files:[],disabled:false,textContent:"",className:"",children:[],handlers:{},
 addEventListener(k,f){this.handlers[k]=f;},appendChild(x){this.children.push(x);},
 replaceChildren(){this.children=[];}};}
function el(id){if(!elements.has(id))elements.set(id,node());return elements.get(id);}
let finishRun;const held=new Promise(resolve=>finishRun=resolve);
const context={document:{getElementById:el,createElement:node},TextDecoder,Uint8Array,
 btoa:s=>Buffer.from(s,"binary").toString("base64"),
 fetch:async(path,options)=>{
 if(path==="/api/programs")return{ok:true,json:async()=>({iso:{description:"stub"}})};
 assert.equal(path,"/api/modes");sent.push(JSON.parse(options.body));await held;
 return{ok:true,json:async()=>({status:"validated_local_subspace",workflow_id:"id",directory:"stub",
 summary:{basis_convention:"local"},files:{"modes.json":"/unused/path","evil.html":"/unused"}})};}};
vm.createContext(context);vm.runInContext(source,context);
function deferred(text){let done;const p=new Promise(resolve=>done=resolve);
 return{file:{name:text,size:text.length,arrayBuffer:()=>p},finish:()=>done(Uint8Array.from(Buffer.from(text)).buffer)};}
async function choose(x){el("mode-file").files=[x.file];return el("mode-file").handlers.change();}
(async()=>{
 await new Promise(resolve=>setImmediate(resolve));
 el("mode-irrep").value="GM1+";el("mode-opd").value="P1";el("mode-timeout").value="60";
 const a=deferred("FIRST"),b=deferred("SECOND"),ap=choose(a);
 assert.equal(el("mode-run").disabled,true);await el("mode-run").handlers.click();assert.equal(sent.length,0);
 const bp=choose(b);b.finish();await bp;a.finish();await ap;assert.equal(el("mode-source").value,"SECOND");
 const c=deferred("STALE"),cp=choose(c);el("mode-source").value="EDITED";el("mode-source").handlers.input();c.finish();await cp;
 assert.equal(el("mode-source").value,"EDITED");
 el("mode-k-values").value="1/3";el("mode-amplitudes").value=".2 -.1";
 const run=el("mode-run").handlers.click();assert.equal(el("run").disabled,true);
 assert.equal(el("mode-run").disabled,true);assert.equal(sent[0].cif,"EDITED");assert.equal("cif_b64" in sent[0],false);
 assert.deepEqual(sent[0].amplitudes,[.2,-.1]);assert.deepEqual(sent[0].k_values,["1/3"]);
 finishRun();await run;assert.equal(el("mode-run").disabled,false);
 assert.equal(el("mode-downloads").children.length,1);
 assert.equal(el("mode-downloads").children[0].textContent,"下载 modes.json");
 const d=deferred("ORIGINAL\r\n"),dp=choose(d);d.finish();await dp;
 await el("mode-run").handlers.click();assert.equal(Buffer.from(sent[1].cif_b64,"base64").toString(),"ORIGINAL\r\n");
 el("mode-amplitudes").value="nan";await el("mode-run").handlers.click();assert.equal(sent.length,2);
 assert.match(el("mode-status").textContent,/finite/);assert.equal(el("mode-downloads").children.length,0);
 process.stdout.write("Mode input race, pending guard, shared UI busy state, dispatch and export allowlist: passed.\n");
})().catch(e=>{console.error(e);process.exitCode=1;});
"""
        script = re.search(r'<script nonce="__TOKEN__">([\s\S]*?)</script>', PAGE).group(1)
        completed = subprocess.run([shutil.which("node"), "-e", harness, json.dumps(script)],
                                   capture_output=True, text=True, timeout=15)
        self.assertEqual(completed.returncode, 0, completed.stderr)


if __name__ == "__main__":
    unittest.main()
