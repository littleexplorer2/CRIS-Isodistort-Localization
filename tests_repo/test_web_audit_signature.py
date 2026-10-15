"""Manual browser audit signatures must cover every production layer."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "ISODISTORT/tests/manual/download_local_web_outputs.mjs"
NODE = shutil.which("node")
pytestmark = pytest.mark.skipif(NODE is None, reason="Node.js is required for the browser audit helper")


def run_node(code: str, *args: Path) -> str:
    result = subprocess.run(
        [NODE, "--input-type=module", "--eval", code, "signature-test", *(str(arg) for arg in args)],
        cwd=ROOT, capture_output=True, text=True, encoding="utf-8", timeout=30,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_source_signature_detects_changes_in_each_layer(tmp_path: Path) -> None:
    fixture_root = tmp_path / "workspace with spaces"
    helper = fixture_root / SCRIPT.relative_to(ROOT)
    helper.parent.mkdir(parents=True)
    shutil.copyfile(SCRIPT, helper)
    code = """
        import { pathToFileURL } from 'node:url';
        const { sourceSignature } = await import(pathToFileURL(process.argv[2]));
        console.log(JSON.stringify(sourceSignature()));
    """
    before = json.loads(run_node(code, helper))
    for layer in ("backend", "application", "frontend"):
        module = fixture_root / "ISODISTORT" / layer / "example" / "__init__.py"
        module.parent.mkdir(parents=True)
        module.write_text("VALUE = 1\n", encoding="utf-8")
        added = json.loads(run_node(code, helper))
        assert added["sha256"] != before["sha256"]
        assert module.relative_to(fixture_root).as_posix() in {row["path"] for row in added["entries"]}
        module.write_text("VALUE = 2\n", encoding="utf-8")
        changed = json.loads(run_node(code, helper))
        assert changed["sha256"] != added["sha256"]
        cache = module.parent / "__pycache__" / "example.pyc"
        cache.parent.mkdir()
        cache.write_bytes(b"ignored bytecode")
        assert json.loads(run_node(code, helper))["sha256"] == changed["sha256"]
        before = changed


def test_frozen_contracts_cover_canonical_modules_and_reject_drift(tmp_path: Path) -> None:
    helper = tmp_path / SCRIPT.relative_to(ROOT)
    helper.parent.mkdir(parents=True)
    shutil.copyfile(SCRIPT, helper)
    manifest = tmp_path / "ISODISTORT/docs/manifests/local_web_output_refresh_manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("{}\n", encoding="utf-8")
    run_node("""
        import assert from 'node:assert/strict';
        import fs from 'node:fs';
        import path from 'node:path';
        import crypto from 'node:crypto';
        import { pathToFileURL } from 'node:url';
        const { assertFrozenContracts, WEB_CONTRACT_PATHS } = await import(pathToFileURL(process.argv[2]));
        const root = process.argv[3];
        const read = rel => fs.readFileSync(path.join(root, rel));
        const manifestPath = 'ISODISTORT/docs/manifests/local_web_output_refresh_manifest.json';
        const manifest = JSON.parse(read(manifestPath));
        const revision = JSON.parse(fs.readFileSync(process.argv[4]));
        assert.deepEqual(new Set(Object.keys(revision.web_contract_files)), new Set(WEB_CONTRACT_PATHS));
        assert(WEB_CONTRACT_PATHS.includes('ISODISTORT/application/workflows/service.py'));
        assert(WEB_CONTRACT_PATHS.includes('ISODISTORT/application/session/core_api.py'));
        const contracts = { ...revision.web_contract_files };
        revision.parent_manifest_sha256 = crypto.createHash('sha256').update(read(manifestPath)).digest('hex');
        assertFrozenContracts(manifest, { web_contracts: contracts }, revision);
        for (const rel of WEB_CONTRACT_PATHS) {
            assert.throws(() => assertFrozenContracts(manifest, {
                web_contracts: { ...contracts, [rel]: '0'.repeat(64) },
            }, revision), /Frozen web contracts drifted/);
        }
        assert.throws(() => assertFrozenContracts(manifest, { web_contracts: contracts }, {
            ...revision, parent_manifest_sha256: '0'.repeat(64),
        }), /different frozen manifest/);
    """, helper, tmp_path, ROOT / "ISODISTORT/docs/manifests/local_web_output_refresh_runtime_revision_20261015.json")
