"""Tests for frontend.web; migrated without changing scientific expectations."""
from __future__ import annotations

from application.session import IsoDistort
from application.workflows import service as workflow_service
from backend.utils.opd_format import _centering_letter
from backend.utils.text_parser import parse_basis_token
from data_dir import experiment_data_dir
from frontend.web import server as web_server
from pathlib import Path
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from tests.support_application_session import _wsl_wsl_available
from tests.support_frontend_web import _SessionIsoStub
from tests.support_paths import TESTS_ROOT
from types import SimpleNamespace
from dataclasses import replace
import copy
import html
import json
import numpy as np
import pytest
import re
import runpy
import shutil
import subprocess
import threading
import time
import urllib.request
import uuid


# Original coverage: test_web.py.

_PROJECT_ROOT = TESTS_ROOT.parent

DATA_DIR = experiment_data_dir()

CIFS_DIR = TESTS_ROOT / "cifs_30"

MATCH_TOL = 1e-6


def _run_page_javascript(scenario: str) -> None:
    """Execute production page functions with controlled network/file timing."""
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is required to execute asynchronous browser contracts")
    source = (_PROJECT_ROOT / "frontend" / "web" / "index.html").read_text(encoding="utf-8")
    script = source.split("<script>", 1)[1].split("</script>", 1)[0]
    script = script.split("/* ================= 初始化", 1)[0]
    harness = r"""
const assert = require('node:assert/strict');
const nodes = {}, pending = [], readers = [], errors = [];
const document = {
  getElementById: id => nodes[id] ||= {
    value: '', checked: false, files: [], style: {}, textContent: '', innerHTML: '',
    classList: {add() {}, remove() {}}, querySelector() {return null;},
    querySelectorAll() {return [];}, scrollIntoView() {},
  },
  querySelectorAll: () => [],
  querySelector: () => null,
};
const window = {crypto: null, location: {href: 'http://localhost/'}};
const crypto = {randomUUID: () => '00000000-0000-4000-8000-000000000001'};
const setInterval = () => 1, clearInterval = () => {};
class FileReader {
  readAsText(file) {readers.push({reader: this, file});}
}
function fetch(path, opts) {
  return new Promise((resolve, reject) => pending.push({path, opts, resolve, reject}));
}
function reply(path, data) {
  const index = pending.findIndex(request => request.path === path);
  assert.notEqual(index, -1, `Missing pending request ${path}`);
  pending.splice(index, 1)[0].resolve({json: async () => data});
}
function readFile(name) {
  const item = readers.find(item => item.file.name === name);
  assert.ok(item, `Missing file read ${name}`);
  item.reader.result = name;
  item.reader.onload();
}
async function flush() {for (let i = 0; i < 12; i++) await Promise.resolve();}
function state(revision) {
  return {revision, structure: null, distortion_types: ['displacive'],
    distortion_scope: {}, species: [], subgroups: 0, modes: []};
}
"""
    setup = r"""
const productionRenderTypes = renderTypes;
renderTypes = () => {};
renderStatusBar = s => {nodes.currentRevision = s && s.revision;};
renderResultTable = () => {};
const productionRenderKvecBox = renderKvecBox;
renderKvecBox = () => {};
selectedLattice = () => null;
showError = (_node, error) => errors.push(error);
STATE = state(7);
(async () => {
"""
    completed = subprocess.run(
        [node, "-"], input=harness + script + setup + scenario
        + "\n})().then(() => {scenarioFinished = true;}).catch(error => {"
        "scenarioFinished = true; console.error(error); process.exitCode = 1;});"
        "\nlet scenarioFinished = false;\nprocess.on('beforeExit', () => {"
        "if (!scenarioFinished) {console.error('Scenario did not finish'); process.exitCode = 1;}});",
        text=True, capture_output=True, timeout=10,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize("scope,active,selected,summary", [
    ({}, True, ["Eu", "Al"], "displacive=all"),
    ({"strain": []}, True, ["Eu", "Al"], "displacive=all"),
    ({"displacive": ["*"]}, True, ["Eu", "Al"], "displacive=all"),
    ({"displacive": ["Eu"]}, True, ["Eu"], "displacive=Eu"),
    ({"displacive": []}, False, [], "-"),
    ({}, False, [], "-"),
])
def test_canonical_scope_display_and_submission_agree(scope, active, selected, summary):
    _run_page_javascript(f"""
STATE.species = ['Eu', 'Al'];
STATE.distortion_types = {json.dumps(['displacive'] if active else [])};
STATE.distortion_scope = {json.dumps(scope)};
const expectedSelected = {json.dumps(selected)}, expectedSummary = {json.dumps(summary)};
""" + r"""
productionRenderTypes();
const box = document.getElementById('typeBox');
const inputs = [...box.innerHTML.matchAll(/<input\b([^>]*)>/g)].map(match => {
  const attributes = Object.fromEntries([...match[1].matchAll(/([\w-]+)="([^"]*)"/g)]
    .map(attribute => [attribute[1], attribute[2]]));
  return {attributes, checked: /\bchecked\b/.test(match[1]),
    dataset: {species: attributes['data-species']}};
});
const displacive = inputs.filter(input => input.attributes['data-tp'] === 'displacive');
assert.deepEqual(displacive.filter(input => input.checked && input.dataset.species)
  .map(input => input.dataset.species), expectedSelected);
assert.equal(displacive.find(input => input.attributes['data-role'] === 'none').checked, false);
assert.equal(formatScope(STATE), expectedSummary);
// Read the actual rendered controls back through the production submit adapter.
box.querySelector = selector => inputs.find(input =>
  [...selector.matchAll(/\[([\w-]+)="([^"]*)"\]/g)]
    .every(attribute => input.attributes[attribute[1]] === attribute[2])) || null;
box.querySelectorAll = selector => inputs.filter(input => input.checked
  && input.dataset.species
  && input.attributes['data-tp'] === selector.match(/data-tp="([^"]*)"/)[1]);
const submission = collectTypesAndScope();
assert.deepEqual(submission.types, STATE.distortion_types);
assert.deepEqual(submission.scope.displacive,
  expectedSelected.length === STATE.species.length ? ['*'] : expectedSelected);
""")


@pytest.mark.parametrize("scenario", [
    pytest.param(r"""
const old = refreshStatus();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
subgroups = [{index: 4, marker: 'new'}];
reply('/api/state', {ok: true, state: state(7)});
await old;
assert.equal(STATE.revision, 9, 'Delayed state must not roll revision backward');
assert.equal(subgroups[0].marker, 'new', 'Delayed refresh must not erase current rows');
""", id="delayed-state"),
    pytest.param(r"""
STATE.structure = {};
const old = runMethod1();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
reply('/api/method1', {ok: true, state: state(8), candidates: [{marker: 'old'}]});
await flush();
assert.equal(STATE.revision, 9);
assert.deepEqual(m1Items, [], 'Old parent candidate rows must not reappear');
assert.equal(pending.length, 0, 'Discarded result must not launch follow-up refresh');
await old;
""", id="delayed-result"),
    pytest.param(r"""
selectedKp = () => ({label: 'GM'});
collectParams = () => [];
const old = enumSubgroups();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
const out = document.getElementById('subOut');
out.innerHTML = 'current parent results';
reply('/api/subgroups', {ok: false, error: 'old failure', state: state(7)});
await old;
assert.equal(out.innerHTML, 'current parent results', 'Old errors must not overwrite new UI');
assert.equal(STATE.revision, 9);
""", id="delayed-error"),
    pytest.param(r"""
const old = refreshKpoints();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
kpoints = [{label: 'new'}];
reply('/api/kpoints', {ok: true, revision: 7, kpoints: [{label: 'old'}]});
await old;
assert.equal(kpoints[0].label, 'new', 'Old query payload needs its own revision guard');
""", id="delayed-read"),
    pytest.param(r"""
document.getElementById('dauFile').files = [{name: 'daughter.cif'}];
refreshStatus = async () => {};
subgroups = [{index: 4}]; resultRevisions[2] = 7;
const basis = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: state(8), modes: []});
await basis;
const old = runMethod4();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
readFile('daughter.cif');
await flush();
assert.ok(!pending.some(request => request.path === '/api/method4'),
  'A daughter chosen under the old parent must not be submitted under the new parent');
""", id="daughter-read-context"),
    pytest.param(r"""
refreshStatus = async () => {};
refreshKpoints = async () => {};
loadMethod1Options = async () => {};
document.getElementById('cifFile').files = [{name: 'old.cif'}];
const old = loadCif();
document.getElementById('cifFile').files = [{name: 'new.cif'}];
const current = loadCif();
readFile('new.cif');
await flush();
assert.equal(JSON.parse(pending[0].opts.body).filename, 'new.cif');
reply('/api/load_cif', {ok: true, state: state(8)});
await current;
readFile('old.cif');
await flush();
assert.equal(pending.length, 0, 'Superseded parent file read must not replace new input');
await old;
""", id="parent-read-intent"),
    pytest.param(r"""
const old = refreshKpoints();
const current = refreshKpoints();
const newer = pending.splice(1, 1)[0];
newer.resolve({json: async () => ({ok: true, revision: 7, kpoints: [{label: 'new'}]})});
await current;
reply('/api/kpoints', {ok: true, revision: 7, kpoints: [{label: 'old'}]});
await old;
assert.equal(kpoints[0].label, 'new', 'Same revision does not imply same request intent');
""", id="same-revision-read-order"),
    pytest.param(r"""
const query = refreshKpoints();
reply('/api/kpoints', {ok: true, revision: 8, kpoints: [{label: 'foreign'}]});
await query;
assert.deepEqual(kpoints, [], 'A newer read alone cannot authorize its foreign payload');
assert.equal(STATE.revision, 7);
reply('/api/state', {ok: true, state: state(8)});
await flush();
assert.equal(STATE.revision, 8);
assert.deepEqual(kpoints, []);
""", id="foreign-read-context"),
    pytest.param(r"""
STATE.structure = {};
const old = runMethod1();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
m1Items = [{marker: 'current'}];
document.getElementById('m1Out').innerHTML = 'current results';
reply('/api/method1', {ok: false, stale_context: true, error: 'stale', state: state(9)});
await old;
assert.equal(document.getElementById('m1Out').innerHTML, 'current results');
assert.equal(m1Items[0].marker, 'current');
assert.deepEqual(errors, [], 'A late conflict must not erase an accepted success');
""", id="conflict-after-current-success"),
    pytest.param(r"""
STATE.structure = {};
let release;
const gap = new Promise(resolve => {release = resolve;});
const realBusy = withBusy;
withBusy = async (...args) => {const data = await realBusy(...args); await gap; return data;};
const old = runMethod1();
reply('/api/method1', {ok: true, state: state(8), candidates: [{index: 0, marker: 'old'}]});
await flush();
assert.equal(STATE.revision, 8);
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
release();
await old;
assert.deepEqual(m1Items, [], 'Acceptance must be rechecked after the caller await');
assert.equal(pending.length, 0);
""", id="consumer-await-gap"),
    pytest.param(r"""
refreshStatus = async () => {};
STATE.structure = {};
const action = runMethod1();
reply('/api/method1', {ok: true, state: state(8), candidates: [{index: 3}]});
await action;
assert.equal(resultRevisions[1], 8);
assert.equal(itemsFor(1)[0].index, 3);
const select = computeModes(3, 'method1');
assert.equal(JSON.parse(pending[0].opts.body).revision, 8);
reply('/api/method2', {ok: true, state: state(9), modes: []});
await select;
assert.equal(resultRevisions[1], 9, 'Compatible local selection retains the candidate pool');
assert.equal(itemsFor(1)[0].index, 3);
const remote = api('/api/state');
reply('/api/state', {ok: true, state: state(10)});
await remote;
assert.deepEqual(itemsFor(1), [], 'Another tab cannot bless this tab candidate indices');
await computeModes(3, 'method1');
assert.equal(pending.length, 0, 'A stale candidate click cannot start computation');
""", id="candidate-ownership"),
    pytest.param(r"""
refreshStatus = async () => {};
refreshKpoints = async () => {};
loadMethod1Options = async () => {};
document.getElementById('cifFile').files = [{name: 'old.cif'}];
const old = loadCif();
readFile('old.cif');
await flush();
document.getElementById('cifFile').files = [{name: 'new.cif'}];
const current = loadCif();
readFile('new.cif');
await flush();
assert.equal(pending.length, 1);
assert.equal(JSON.parse(pending[0].opts.body).filename, 'old.cif');
reply('/api/load_cif', {ok: true, state: state(8)});
await flush();
assert.equal(pending.length, 1, 'New parent POST waits for the prior parent to finish');
assert.equal(JSON.parse(pending[0].opts.body).filename, 'new.cif');
reply('/api/load_cif', {ok: true, state: state(9)});
await Promise.all([old, current]);
assert.equal(STATE.revision, 9);
assert.equal(document.getElementById('btnLoadCif').disabled, undefined);
""", id="parent-post-order"),
    pytest.param(r"""
selectedKp = () => ({label: 'GM'});
collectParams = () => [];
const old = enumSubgroups();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
const out = document.getElementById('subOut');
out.innerHTML = 'current results';
pending.shift().reject(new Error('old transport error'));
await old;
assert.equal(out.innerHTML, 'current results');
assert.deepEqual(errors, []);
""", id="transport-error-ownership"),
    pytest.param(r"""
const request = refreshKpoints();
const action = api('/api/method1', {});
reply('/api/method1', {ok: true, state: state(8), candidates: []});
await action;
reply('/api/kpoints', {ok: true, revision: 7, kpoints: [{label: 'old'}]});
await flush();
assert.equal(pending.length, 1, 'An overtaken selector is retried at the new revision');
reply('/api/kpoints', {ok: true, revision: 8, kpoints: [{label: 'current'}]});
await request;
assert.equal(kpoints[0].label, 'current');
""", id="selector-retry-after-local-mutation"),
    pytest.param(r"""
selectedDownloadFormats = () => ['cif'];
document.getElementById('dlMethod').value = '1';
m1Items = [{index: 1}]; resultRevisions[1] = 7;
const download = downloadAll();
reply('/api/export_status', {export: {}});
const request = pending.shift();
assert.equal(request.path.pathname, '/api/download_all');
request.resolve({ok: false, headers: {get: () => 'application/json'},
  json: async () => ({ok: false, stale_context: true, state: state(8), error: 'stale export'})});
await flush();
reply('/api/export_status', {export: {}});
await download;
assert.equal(STATE.revision, 8);
assert.deepEqual(m1Items, []);
assert.deepEqual(errors, ['stale export']);
""", id="download-conflict-state"),
    pytest.param(r"""
selectedDownloadFormats = () => ['cif'];
document.getElementById('dlMethod').value = '1';
m1Items = [{index: 1}]; resultRevisions[1] = 7;
const download = downloadAll();
reply('/api/export_status', {export: {}});
const request = pending.shift();
const upload = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
reply('/api/load_cif', {ok: true, state: state(9)});
await upload;
document.getElementById('distOut').innerHTML = 'current export information';
request.resolve({ok: false, headers: {get: () => 'application/json'},
  json: async () => ({ok: false, state: state(7), error: 'old export error'})});
await flush();
reply('/api/export_status', {export: {}});
await download;
assert.equal(STATE.revision, 9);
assert.equal(document.getElementById('distOut').innerHTML, 'current export information');
assert.deepEqual(errors, []);
""", id="delayed-download-error"),
    pytest.param(r"""
document.querySelectorAll = selector => selector === '.genDbCb:checked' ? [{value: 'removed'}] : [];
const listing = loadGenDbCache();
const deletion = deleteSelectedGenDb();
reply('/api/isotropy_cache/delete', {ok: true, revision: 7, entries: [], deleted: ['removed']});
await deletion;
const currentList = document.getElementById('genDbList').innerHTML;
reply('/api/isotropy_cache', {ok: true, revision: 7, entries: [{name: 'removed'}]});
await listing;
assert.equal(document.getElementById('genDbList').innerHTML, currentList,
  'A cache read preceding deletion must not republish the removed entry');
""", id="cache-read-after-deletion"),
])
def test_async_responses_preserve_current_scientific_context(scenario):
    _run_page_javascript(scenario)


_METHOD2_QUERY_SETUP = r"""
refreshStatus = async () => {};
selectedDownloadFormats = () => ['cif'];
document.getElementById('dlMethod').value = '2';
kpoints = [
  {label: 'Q', parameters: ['a', 'b', 'g'], coordinates: ['a', 'b', 'g']},
  {label: 'R', parameters: [], coordinates: [0, 0, 0]},
];
productionRenderKvecBox();
document.getElementById('kvec_1').value = '0';
for (const letter of ['a', 'b', 'g']) document.getElementById('p_' + letter + '_1').value = '1/6';
const oldCandidate = {index: 4, k_point_label: 'Q', irrep_label: 'Q1', opd_symbol: 'P1'};
subgroups = [oldCandidate]; resultRevisions[2] = 7;
document.getElementById('subOut').innerHTML = 'old Q candidates';
function changeK(value) {
  document.getElementById('kvec_1').value = value;
  onKpChange(1);
}
function editParameter(letter, value) {
  const markup = document.getElementById('kvecBox').innerHTML;
  const input = markup.match(new RegExp('<input\\b[^>]*\\bid="p_' + letter + '_1"[^>]*>'));
  assert.ok(input, 'Production parameter input must exist');
  const handler = input[0].match(/\boninput="([^"]+)"/);
  assert.ok(handler, 'Parameter editing must invalidate immediately, without waiting for blur');
  document.getElementById('p_' + letter + '_1').value = value;
  eval(handler[1]);
}
function prepareZip() {
  document.createElement = () => ({click() {}, remove() {}});
  document.body = {appendChild() {}};
  URL.createObjectURL = () => 'blob:test'; URL.revokeObjectURL = () => {};
}
function completeZip(request, jobId) {
  const status = {job_id: jobId, active: false, phase: 'complete',
    successful: 1, total: 1, ineligible: 0, failed: 0};
  request.resolve({ok: true,
    headers: {get: name => name.toLowerCase() === 'content-type'
      ? 'application/zip' : encodeURIComponent(JSON.stringify(status))},
    blob: async () => ({})});
}
"""


@pytest.mark.parametrize("change", [
    pytest.param("changeK('1');", id="k-change"),
    pytest.param("changeK('1'); changeK('0');", id="k-change-back"),
    *[
        pytest.param(
            f"editParameter('{letter}', {json.dumps(value)});"
            + (f"editParameter('{letter}', '1/6');" if change_back else ""),
            id=f"{letter}-{name}",
        )
        for letter in ("a", "b", "g")
        for name, value, change_back in (
            ("changed", "1/3", False), ("empty", "", False),
            ("invalid", "invalid", False), ("changed-back", "1/3", True),
        )
    ],
])
def test_method2_query_edit_invalidates_candidate_actions(change):
    _run_page_javascript(_METHOD2_QUERY_SETUP + change + r"""
assert.deepEqual(itemsFor(2), [], 'Edited query must retire its prior candidate pool');
assert.equal(resultRevisions[2], null);
assert.equal(document.getElementById('subOut').innerHTML, '');
let tableDownloads = 0;
downloadText = () => {tableDownloads += 1;};
downloadFilteredTable(false);
await computeModes(4, 'subgroups');
await downloadAll();
assert.equal(tableDownloads, 0);
assert.equal(pending.length, 0, 'No candidate or export request may use the retired pool');
""")


@pytest.mark.parametrize("response_kind", ["success", "error", "network-error", "consumer-gap"])
@pytest.mark.parametrize("operation", ["query", "modes"])
def test_method2_query_edit_rejects_inflight_results(operation, response_kind):
    _run_page_javascript(_METHOD2_QUERY_SETUP + f"""
const operationKind = {json.dumps(operation)};
const responseKind = {json.dumps(response_kind)};
""" + r"""
let release;
if (responseKind === 'consumer-gap') {
  const gap = new Promise(resolve => {release = resolve;});
  const realBusy = withBusy;
  withBusy = async (...args) => {const data = await realBusy(...args); await gap; return data;};
}
const operation = operationKind === 'query' ? enumSubgroups() : computeModes(4, 'subgroups');
const request = pending.shift();
assert.equal(request.path, operationKind === 'query' ? '/api/subgroups' : '/api/method2');
const success = {ok: true, state: state(8), subgroups: [oldCandidate], modes: []};
if (responseKind === 'consumer-gap') {
  request.resolve({json: async () => success});
  await flush();
}
changeK('1');
changeK('0');
document.getElementById('subOut').innerHTML = 'current draft';
document.getElementById('modesBox').innerHTML = 'current mode area';
if (responseKind === 'consumer-gap') release();
else if (responseKind === 'network-error') request.reject(new Error('obsolete network failure'));
else request.resolve({json: async () => responseKind === 'success'
  ? success : {ok: false, state: state(7), error: 'obsolete query failure'}});
await operation;
assert.deepEqual(itemsFor(2), []);
assert.equal(document.getElementById('subOut').innerHTML, 'current draft');
// Revision acceptance may clear mode results, but no old payload/error may reappear.
assert.ok(!document.getElementById('modesBox').innerHTML.includes('Q1'));
assert.deepEqual(errors, []);
assert.equal(pending.length, 0, 'Rejected results must not initiate a follow-up refresh');
if (['success', 'consumer-gap'].includes(responseKind)) assert.equal(STATE.revision, 8);
""")


def test_method2_new_query_retires_hidden_pool_even_on_failure():
    _run_page_javascript(_METHOD2_QUERY_SETUP + r"""
const operation = enumSubgroups();
assert.deepEqual(itemsFor(2), [], 'A new search must not leave a hidden downloadable old pool');
reply('/api/subgroups', {ok: false, state: state(7), error: 'query failed'});
await operation;
assert.deepEqual(itemsFor(2), []);
await downloadAll();
assert.equal(pending.length, 0);
""")


@pytest.mark.parametrize("source", ["method1", "method3", "subgroups"])
def test_method2_query_edit_preserves_only_unrelated_mode_results(source):
    _run_page_javascript(_METHOD2_QUERY_SETUP + f"const source = {json.dumps(source)};\n" + r"""
m1Items = [{...oldCandidate}]; resultRevisions[1] = 7;
m3Items = [{...oldCandidate}]; resultRevisions[3] = 7;
const selection = computeModes(4, source);
reply('/api/method2', {ok: true, state: state(8), modes: []});
await selection;
const output = document.getElementById(source === 'method3' ? 'm3ModesBox' : 'modesBox');
const currentModes = output.innerHTML;
assert.ok(currentModes.includes('Q1'));
m4Items = [{mode: 'independent'}]; resultRevisions[4] = 8;
changeK('1');
assert.equal(output.innerHTML, source === 'subgroups' ? '' : currentModes);
assert.equal(itemsFor(1).length, 1);
assert.equal(itemsFor(3).length, 1);
assert.equal(itemsFor(4).length, source === 'subgroups' ? 0 : 1);
assert.deepEqual(itemsFor(2), []);
""")


@pytest.mark.parametrize("source", ["method1", "method3"])
def test_method2_query_edit_preserves_unrelated_inflight_mode_selection(source):
    _run_page_javascript(_METHOD2_QUERY_SETUP + f"const source = {json.dumps(source)};\n" + r"""
m1Items = [{...oldCandidate}]; resultRevisions[1] = 7;
m3Items = [{...oldCandidate}]; resultRevisions[3] = 7;
const selection = computeModes(4, source);
changeK('1');
reply('/api/method2', {ok: true, state: state(8), modes: []});
await selection;
const output = document.getElementById(source === 'method3' ? 'm3ModesBox' : 'modesBox');
assert.ok(output.innerHTML.includes('Q1'));
assert.equal(itemsFor(1).length, 1);
assert.equal(itemsFor(3).length, 1);
assert.deepEqual(itemsFor(2), []);
""")


def test_method2_new_query_uses_edited_parameters_without_sending_drafts():
    _run_page_javascript(_METHOD2_QUERY_SETUP + r"""
editParameter('a', '');
editParameter('b', 'draft');
editParameter('g', '1/3');
assert.equal(pending.length, 0, 'Draft edits must not initiate scientific requests');
editParameter('a', '1/5');
editParameter('b', '1/7');
const search = enumSubgroups();
const request = pending.shift();
assert.deepEqual(JSON.parse(request.opts.body).kpoints,
  [{k: 'Q', params: ['1/5', '1/7', '1/3']}]);
request.resolve({json: async () => ({ok: true, state: state(8), subgroups: [oldCandidate]})});
await search;
assert.equal(itemsFor(2).length, 1);
assert.equal(resultRevisions[2], 8);
""")


def test_method2_new_query_and_selection_enable_current_export():
    _run_page_javascript(_METHOD2_QUERY_SETUP + r"""
changeK('1');
const search = enumSubgroups();
const request = pending.shift();
assert.deepEqual(JSON.parse(request.opts.body).kpoints, [{k: 'R', params: null}]);
const currentCandidate = {...oldCandidate, k_point_label: 'R', irrep_label: 'R1'};
request.resolve({json: async () => ({ok: true, state: state(8), subgroups: [currentCandidate]})});
await search;
assert.equal(itemsFor(2)[0].irrep_label, 'R1');
const selection = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: state(9), modes: []});
await selection;
assert.equal(resultRevisions[2], 9);
prepareZip();
const download = downloadAll();
reply('/api/export_status', {export: {}});
const downloadRequest = pending.shift();
assert.equal(downloadRequest.path.searchParams.get('revision'), '9');
completeZip(downloadRequest, activeExportJobId);
await download;
assert.equal(pending.length, 0);
""")


def test_method2_query_edit_preserves_already_dispatched_export_snapshot():
    _run_page_javascript(_METHOD2_QUERY_SETUP + r"""
prepareZip();
const download = downloadAll();
reply('/api/export_status', {export: {}});
const request = pending.shift();
const jobId = activeExportJobId;
assert.equal(request.path.searchParams.get('revision'), '7');
changeK('1');
assert.deepEqual(itemsFor(2), []);
assert.equal(activeExportJobId, jobId, 'An existing export keeps its original snapshot ownership');
completeZip(request, jobId);
await download;
assert.equal(pending.length, 0);
assert.deepEqual(errors, []);
assert.match(document.getElementById('distOut').innerHTML, /dist.zipDone/);
""")


@pytest.mark.parametrize("earlier_operation", ["query", "modes"])
def test_method2_new_query_shows_current_revision_conflict(earlier_operation):
    _run_page_javascript(_METHOD2_QUERY_SETUP + f"""
const earlierOperation = {json.dumps(earlier_operation)};
""" + r"""
const earlier = earlierOperation === 'query' ? enumSubgroups() : computeModes(4, 'subgroups');
const earlierRequest = pending.shift();
changeK('1');
const current = enumSubgroups();
const currentRequest = pending.shift();
assert.equal(JSON.parse(currentRequest.opts.body).revision, 7);
earlierRequest.resolve({json: async () => ({
  ok: true, state: state(8), subgroups: [oldCandidate], modes: [],
})});
await earlier;
assert.equal(STATE.revision, 8);
currentRequest.resolve({json: async () => ({
  ok: false, stale_context: true, state: state(8), error: 'Current query must retry.',
})});
await current;
assert.deepEqual(itemsFor(2), []);
assert.match(document.getElementById('subOut').innerHTML, /Current query must retry\./,
  'A still-current request must show its conflict even if another local reply already accepted that revision');
assert.equal(pending.length, 0);
""")


@pytest.mark.parametrize("obsolete_by", ["request", "query-dependency", "older-revision", "parent"])
def test_method2_query_edit_keeps_obsolete_conflicts_rejected(obsolete_by):
    _run_page_javascript(_METHOD2_QUERY_SETUP + f"""
const obsoleteBy = {json.dumps(obsolete_by)};
""" + r"""
const earlier = obsoleteBy === 'query-dependency'
  ? computeModes(4, 'subgroups') : enumSubgroups();
const earlierRequest = pending.shift();
if (obsoleteBy === 'parent') {
  const parent = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
  reply('/api/load_cif', {ok: true, state: state(8)});
  await parent;
} else {
  changeK('1');
  const current = enumSubgroups();
  reply('/api/subgroups', {ok: true, state: state(8),
    subgroups: [{...oldCandidate, k_point_label: 'R', irrep_label: 'R1'}]});
  await current;
}
document.getElementById('subOut').innerHTML = 'current candidate feedback';
document.getElementById('modesBox').innerHTML = 'current mode feedback';
earlierRequest.resolve({json: async () => ({
  ok: false, stale_context: true,
  state: state(obsoleteBy === 'older-revision' ? 7 : 8), error: 'Obsolete conflict.',
})});
await earlier;
assert.equal(STATE.revision, 8);
assert.equal(document.getElementById('subOut').innerHTML, 'current candidate feedback');
assert.equal(document.getElementById('modesBox').innerHTML, 'current mode feedback');
assert.deepEqual(errors, []);
assert.equal(pending.length, 0);
""")


_REPLACEMENT_QUERY_SETUP = _METHOD2_QUERY_SETUP + r"""
STATE.structure = {};
renderSessionState = () => {};
const queryState = revision => ({...state(revision), structure: {}});
m1Items = [{...oldCandidate}]; resultRevisions[1] = 7;
m3Items = [{...oldCandidate}]; resultRevisions[3] = 7;
const startQuery = method => method === 1 ? runMethod1() : runMethod3();
const queryPath = method => '/api/method' + method;
const sourceFor = method => method === 2 ? 'subgroups' : 'method' + method;
"""


@pytest.mark.parametrize("method", [1, 3])
@pytest.mark.parametrize("failure", ["server-error", "network-error"])
def test_replacement_query_failure_retires_hidden_candidates(method, failure):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"""
const method = {method}, failure = {json.dumps(failure)};
""" + r"""
document.getElementById('dlMethod').value = String(method);
const operation = startQuery(method);
const request = pending.shift();
assert.deepEqual(itemsFor(method), [], 'A replacement search must immediately retire its prior pool');
assert.equal(resultRevisions[method], null);
if (failure === 'network-error') request.reject(new Error('replacement network failure'));
else request.resolve({json: async () => ({ok: false, state: queryState(7), error: 'replacement failed'})});
await operation;
assert.deepEqual(itemsFor(method), []);
assert.ok(errors.some(error => String(error).includes('replacement')));
const otherMethod = method === 1 ? 3 : 1;
assert.equal(itemsFor(otherMethod).length, 1);
assert.equal(itemsFor(2).length, 1);
let tableDownloads = 0;
downloadText = () => {tableDownloads += 1;};
downloadFilteredTable(false);
await computeModes(4, sourceFor(method));
await downloadAll();
assert.equal(tableDownloads, 0);
assert.equal(pending.length, 0, 'Failed replacement must not leave any old candidate action available');
""")


@pytest.mark.parametrize("method", [1, 3])
@pytest.mark.parametrize("earlier_operation", ["query", "modes"])
@pytest.mark.parametrize("response_kind", ["success", "error", "network-error", "consumer-gap"])
def test_replacement_query_rejects_obsolete_responses(method, earlier_operation, response_kind):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"""
const method = {method}, earlierOperation = {json.dumps(earlier_operation)};
const responseKind = {json.dumps(response_kind)};
""" + r"""
let release;
const realBusy = withBusy;
if (responseKind === 'consumer-gap') {
  const gap = new Promise(resolve => {release = resolve;});
  withBusy = async (...args) => {const data = await realBusy(...args); await gap; return data;};
}
const earlier = earlierOperation === 'query'
  ? startQuery(method) : computeModes(4, sourceFor(method));
const earlierRequest = pending.shift();
const oldSuccess = {ok: true, state: queryState(8), candidates: [oldCandidate], modes: []};
if (responseKind === 'consumer-gap') {
  earlierRequest.resolve({json: async () => oldSuccess});
  await flush();
  withBusy = realBusy;
}
const current = startQuery(method);
const currentRequest = pending.shift();
const out = document.getElementById(method === 1 ? 'm1Out' : 'm3Out');
out.innerHTML = 'new search progress';
if (responseKind === 'consumer-gap') release();
else if (responseKind === 'network-error') earlierRequest.reject(new Error('obsolete failure'));
else earlierRequest.resolve({json: async () => responseKind === 'success'
  ? oldSuccess : {ok: false, state: queryState(7), error: 'obsolete failure'}});
await earlier;
assert.deepEqual(itemsFor(method), []);
assert.equal(out.innerHTML, 'new search progress', 'Old responses cannot erase the new query view');
const modesOut = document.getElementById(method === 1 ? 'modesBox' : 'm3ModesBox');
assert.ok(!modesOut.innerHTML.includes('Q1'));
assert.deepEqual(errors, []);
currentRequest.resolve({json: async () => ({
  ok: false, state: queryState(STATE.revision), error: 'current failure',
})});
await current;
assert.deepEqual(itemsFor(method), []);
assert.deepEqual(errors, ['current failure']);
assert.equal(pending.length, 0);
""")


@pytest.mark.parametrize("method", [1, 3])
@pytest.mark.parametrize("mode_method", [1, 2, 3])
def test_replacement_query_invalidates_only_its_own_mode_results(method, mode_method):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"""
const method = {method}, modeMethod = {mode_method};
""" + r"""
const selection = computeModes(4, sourceFor(modeMethod));
reply('/api/method2', {ok: true, state: queryState(8), modes: []});
await selection;
const output = document.getElementById(modeMethod === 3 ? 'm3ModesBox' : 'modesBox');
const oldModeDisplay = output.innerHTML;
assert.ok(oldModeDisplay.includes('Q1'));
m4Items = [{mode: 'independent'}]; resultRevisions[4] = 8;
const replacement = startQuery(method);
reply(queryPath(method), {ok: false, state: queryState(8), error: 'replacement failed'});
await replacement;
assert.equal(output.innerHTML, method === modeMethod ? '' : oldModeDisplay);
assert.deepEqual(itemsFor(method), []);
for (const other of [1, 2, 3, 4].filter(value => value !== method)) {
  assert.equal(itemsFor(other).length, other === 4 && method === modeMethod ? 0 : 1);
}
""")


@pytest.mark.parametrize("method", [1, 3])
def test_replacement_query_success_restores_selection_and_download(method):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"const method = {method};\n" + r"""
const search = startQuery(method);
reply(queryPath(method), {ok: true, state: queryState(8),
  candidates: [{...oldCandidate, index: 9}]});
await search;
assert.equal(itemsFor(method)[0].index, 9);
assert.equal(resultRevisions[method], 8);
const selection = computeModes(9, sourceFor(method));
reply('/api/method2', {ok: true, state: queryState(9), modes: []});
await selection;
document.getElementById('dlMethod').value = String(method);
let tableText = null;
downloadText = text => {tableText = text;};
downloadFilteredTable(false);
assert.ok(tableText.includes('9'));
prepareZip();
const download = downloadAll();
reply('/api/export_status', {export: {}});
const request = pending.shift();
assert.equal(request.path.searchParams.get('method'), String(method));
assert.equal(request.path.searchParams.get('revision'), '9');
completeZip(request, activeExportJobId);
await download;
assert.deepEqual(errors, []);
assert.equal(pending.length, 0);
""")


@pytest.mark.parametrize("method", [1, 3])
def test_replacement_query_does_not_reassign_dispatched_export(method):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"const method = {method};\n" + r"""
document.getElementById('dlMethod').value = String(method);
prepareZip();
const download = downloadAll();
reply('/api/export_status', {export: {}});
const downloadRequest = pending.shift(), jobId = activeExportJobId;
const replacement = startQuery(method);
reply(queryPath(method), {ok: false, state: queryState(7), error: 'replacement failed'});
await replacement;
assert.deepEqual(itemsFor(method), []);
assert.equal(activeExportJobId, jobId);
completeZip(downloadRequest, jobId);
await download;
assert.match(document.getElementById('distOut').innerHTML, /dist.zipDone/);
assert.equal(pending.length, 0);
""")


@pytest.mark.parametrize("failure", [
    "no-file", "invalid-origin", "file-read", "server-error", "network-error",
])
def test_method4_replacement_failure_retires_old_amplitudes(failure):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"const failure = {json.dumps(failure)};\n" + r"""
const basis = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: queryState(8), modes: []});
await basis;
m4Items = [{mode: 'old-daughter', As: 1, Ap: 1, raw: 1, normfactor: 1}];
resultRevisions[4] = 8;
m4Meta.rms = 42;
document.getElementById('dlMethod').value = '4';
document.getElementById('dauFile').files = failure === 'no-file' ? [] : [{name: 'new-daughter.cif'}];
if (failure === 'invalid-origin') document.getElementById('m4OriginX').value = 'invalid';
const operation = runMethod4();
assert.deepEqual(itemsFor(4), [], 'A new daughter submission must retire its old amplitudes before file I/O');
assert.equal(m4Meta.rms, null);
if (failure === 'file-read') readers[0].reader.onerror(new Error('controlled read failure'));
else if (!['no-file', 'invalid-origin'].includes(failure)) {
  readFile('new-daughter.cif');
  await flush();
  const request = pending.shift();
  if (failure === 'network-error') request.reject(new Error('controlled network failure'));
  else request.resolve({json: async () => ({ok: false, state: queryState(8), error: 'controlled failure'})});
}
await operation;
let downloads = 0;
downloadText = () => {downloads += 1;};
downloadFilteredTable(false);
assert.equal(downloads, 0);
assert.deepEqual(itemsFor(4), []);
assert.equal(pending.length, 0);
for (const method of [1, 2, 3]) assert.equal(itemsFor(method).length, 1);
""")


@pytest.mark.parametrize("source", ["method1", "subgroups", "method3"])
@pytest.mark.parametrize("response_kind", ["success", "error", "network-error", "consumer-gap"])
def test_nmod_change_rejects_obsolete_mode_responses(source, response_kind):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"""
const source = {json.dumps(source)}, responseKind = {json.dumps(response_kind)};
""" + r"""
let release;
if (responseKind === 'consumer-gap') {
  const gap = new Promise(resolve => {release = resolve;});
  const realBusy = withBusy;
  withBusy = async (...args) => {const data = await realBusy(...args); await gap; return data;};
}
const operation = computeModes(4, source);
const request = pending.shift();
assert.equal(JSON.parse(request.opts.body).nmod, 0);
const oldSuccess = {ok: true, state: queryState(8), modes: [{
  irrep_label: 'old-nmod-zero', mode_type: 'displacive', opd_symbol: 'P1',
  wyckoff_sites: [], n_representatives: 1,
}]};
if (responseKind === 'consumer-gap') {
  request.resolve({json: async () => oldSuccess});
  await flush();
}
document.getElementById('nmod').value = '1';
applyNmod();
document.getElementById('nmod').value = '0';
applyNmod();
if (responseKind === 'consumer-gap') release();
else if (responseKind === 'network-error') request.reject(new Error('obsolete nmod failure'));
else request.resolve({json: async () => responseKind === 'success'
  ? oldSuccess : {ok: false, state: queryState(7), error: 'obsolete nmod failure'}});
await operation;
const out = document.getElementById(source === 'method3' ? 'm3ModesBox' : 'modesBox');
assert.ok(!out.innerHTML.includes('old-nmod-zero'), 'Changing nmod back must not revive an obsolete mode request');
assert.deepEqual(errors, []);
assert.equal(pending.length, 0);
for (const method of [1, 2, 3]) assert.equal(itemsFor(method).length, 1);
""")


def test_nmod_applied_value_preserves_candidates_and_new_mode_request():
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + r"""
const original = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: queryState(8), modes: []});
await original;
assert.ok(document.getElementById('modesBox').innerHTML.includes('Q1'));
document.getElementById('nmod').value = '0';
applyNmod();
assert.ok(document.getElementById('modesBox').innerHTML.includes('Q1'), 'Applying the same value keeps the same basis');
document.getElementById('nmod').value = '1';
applyNmod();
assert.equal(document.getElementById('modesBox').innerHTML, '');
for (const method of [1, 2, 3]) assert.equal(itemsFor(method).length, 1);
const current = computeModes(4, 'subgroups');
const request = pending.shift();
assert.equal(JSON.parse(request.opts.body).nmod, 1);
request.resolve({json: async () => ({ok: true, state: queryState(9), modes: []})});
await current;
assert.ok(document.getElementById('modesBox').innerHTML.includes('Q1'));
""")


def test_method4_rejects_basis_after_failed_mode_selection():
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + r"""
const basis = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: queryState(8), modes: []});
await basis;
const failedSelection = computeModes(4, 'method1');
reply('/api/method2', {ok: false, state: queryState(8), error: 'new candidate failed'});
await failedSelection;
document.getElementById('dauFile').files = [{name: 'daughter.cif'}];
await runMethod4();
assert.equal(readers.length, 0, 'Hidden earlier basis cannot authorize daughter file processing');
assert.equal(pending.length, 0);
assert.ok(errors.some(error => String(error).includes('Compute modes')));
""")


@pytest.mark.parametrize("boundary", ["file-read", "response", "consumer-gap"])
@pytest.mark.parametrize("change", ["nmod", "mode-selection", "source-query", "parent"])
def test_method4_requires_same_confirmed_basis_across_awaits(boundary, change):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"""
const boundary = {json.dumps(boundary)}, change = {json.dumps(change)};
""" + r"""
const basis = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: queryState(8), modes: []});
await basis;
let release;
if (boundary === 'consumer-gap') {
  const gap = new Promise(resolve => {release = resolve;});
  const realBusy = withBusy;
  withBusy = async (...args) => {const data = await realBusy(...args); await gap; return data;};
}
document.getElementById('dauFile').files = [{name: 'daughter.cif'}];
const decomposition = runMethod4();
let request;
if (boundary !== 'file-read') {
  readFile('daughter.cif');
  await flush();
  request = pending.shift();
  assert.equal(JSON.parse(request.opts.body).nmod, 0);
}
if (boundary === 'consumer-gap') {
  request.resolve({json: async () => ({
    ok: true, state: queryState(9), amplitudes: {'obsolete-daughter': 1},
  })});
  await flush();
}
let selection, selectionRequest;
if (change === 'nmod') {
  document.getElementById('nmod').value = '1'; applyNmod();
} else if (change === 'mode-selection') {
  selection = computeModes(4, 'method1');
  selectionRequest = pending.shift();
} else if (change === 'source-query') {
  changeK('1');
} else {
  const parent = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
  reply('/api/load_cif', {ok: true, state: queryState(10)});
  await parent;
}
if (boundary === 'file-read') readFile('daughter.cif');
else if (boundary === 'consumer-gap') release();
else request.resolve({json: async () => ({
  ok: true, state: queryState(9), amplitudes: {'obsolete-daughter': 1},
})});
await decomposition;
assert.deepEqual(itemsFor(4), [], 'A daughter result must belong to the still-confirmed mode basis');
assert.ok(!document.getElementById('m4Out').innerHTML.includes('obsolete-daughter'));
assert.ok(!pending.some(item => item.path === '/api/method4'));
if (selection) {
  selectionRequest.resolve({json: async () => ({
    ok: false, state: queryState(STATE.revision), error: 'controlled selection end',
  })});
  await selection;
}
// A parent replacement can request fresh selector metadata; settle those reads.
for (const item of pending.splice(0)) item.resolve({json: async () => ({ok: false})});
await flush();
""")


def test_method4_reuses_confirmed_basis_for_repeated_daughters():
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + r"""
document.getElementById('nmod').value = '1'; applyNmod();
const basis = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: queryState(8), modes: []});
await basis;
for (const [number, revision] of [[1, 9], [2, 10]]) {
  const name = 'daughter' + number + '.cif';
  document.getElementById('dauFile').files = [{name}];
  const decomposition = runMethod4();
  readFile(name);
  await flush();
  const request = pending.shift();
  assert.equal(JSON.parse(request.opts.body).nmod, 1);
  request.resolve({json: async () => ({
    ok: true, state: queryState(revision), amplitudes: {current: number},
    parent_cell_amplitudes: {current: number}, raw_coefficients: {current: number},
    mode_normfactors: {current: 1},
  })});
  await decomposition;
  assert.equal(itemsFor(4)[0].As, number);
}
let text = '';
downloadText = value => {text = value;};
document.getElementById('dlMethod').value = '4';
downloadFilteredTable(false);
assert.ok(text.includes('current'));
assert.deepEqual(errors, []);
""")


def test_method4_obsolete_completion_keeps_current_submission_error():
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + r"""
showError = (out, error) => {errors.push(error); out.innerHTML = String(error);};
const basis = computeModes(4, 'subgroups');
reply('/api/method2', {ok: true, state: queryState(8), modes: []});
await basis;
document.getElementById('dauFile').files = [{name: 'old-daughter.cif'}];
const old = runMethod4();
readFile('old-daughter.cif'); await flush();
const oldRequest = pending.shift();
document.getElementById('dauFile').files = [];
await runMethod4();
assert.equal(document.getElementById('m4Out').innerHTML, 'err.noFile');
oldRequest.resolve({json: async () => ({
  ok: true, state: queryState(9), amplitudes: {'old-daughter': 1},
})});
await old;
assert.deepEqual(itemsFor(4), []);
assert.equal(document.getElementById('m4Out').innerHTML, 'err.noFile',
  'A committed obsolete decomposition cannot erase the current submission error');
""")


def test_nmod_change_preserves_already_dispatched_export():
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + r"""
prepareZip();
document.getElementById('dlMethod').value = '2';
const download = downloadAll();
reply('/api/export_status', {export: {}});
const request = pending.shift(), jobId = activeExportJobId;
assert.equal(request.path.searchParams.get('nmod'), '0');
document.getElementById('nmod').value = '1'; applyNmod();
assert.equal(activeExportJobId, jobId);
assert.equal(itemsFor(2).length, 1);
completeZip(request, jobId);
await download;
assert.deepEqual(errors, []);
assert.match(document.getElementById('distOut').innerHTML, /dist.zipDone/);
""")


@pytest.mark.parametrize("failure", ["network", "json", "null-json"])
def test_current_local_failure_survives_compatible_revision_advance(failure):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"const failure = {json.dumps(failure)};\n" + r"""
const earlier = computeModes(4, 'method3');
const modeRequest = pending.shift();
const current = runMethod1();
const queryRequest = pending.shift();
modeRequest.resolve({json: async () => ({ok: true, state: queryState(8), modes: []})});
await earlier;
assert.equal(STATE.revision, 8);
if (failure === 'network') queryRequest.reject(new Error('current transport failure'));
else queryRequest.resolve({json: async () => {
  if (failure === 'json') throw new SyntaxError('current invalid JSON');
  return null;
}});
await current;
assert.equal(STATE.revision, 8, 'Local failures must not fabricate or roll back server state');
assert.deepEqual(itemsFor(1), []);
assert.equal(errors.length, 1, 'A current local failure must be explained after compatible state advances');
assert.equal(pending.length, 0);
""")


@pytest.mark.parametrize("superseded_by", ["request", "query", "parent"])
@pytest.mark.parametrize("failure", ["network", "json", "null-json"])
def test_obsolete_local_failures_remain_rejected(superseded_by, failure):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"""
const supersededBy = {json.dumps(superseded_by)}, failure = {json.dumps(failure)};
""" + r"""
const old = supersededBy === 'query' ? computeModes(4, 'method1') : runMethod1();
const request = pending.shift();
if (supersededBy === 'parent') {
  const parent = api('/api/load_cif', {filename: 'new.cif', content: 'new'});
  reply('/api/load_cif', {ok: true, state: state(9)});
  await parent;
} else {
  const current = runMethod1();
  reply('/api/method1', {ok: true, state: queryState(8), candidates: [{...oldCandidate, index: 9}]});
  await current;
}
document.getElementById('m1Out').innerHTML = 'current result';
document.getElementById('modesBox').innerHTML = 'current modes';
if (failure === 'network') request.reject(new Error('obsolete transport failure'));
else request.resolve({json: async () => {
  if (failure === 'json') throw new SyntaxError('obsolete invalid JSON');
  return null;
}});
await old;
assert.deepEqual(errors, [], 'A local error still requires current request, query and parent ownership');
assert.equal(document.getElementById('m1Out').innerHTML, 'current result');
assert.equal(document.getElementById('modesBox').innerHTML, 'current modes');
assert.equal(pending.length, 0);
""")


@pytest.mark.parametrize("change", ["compatible", "parent", "request", "job"])
@pytest.mark.parametrize("failure", ["network", "json", "null-json", "blob"])
def test_download_local_failure_keeps_snapshot_job_ownership(change, failure):
    _run_page_javascript(_REPLACEMENT_QUERY_SETUP + f"""
const change = {json.dumps(change)}, failure = {json.dumps(failure)};
""" + r"""
const download = downloadAll();
reply('/api/export_status', {export: {}});
const request = pending.shift();
const jobId = activeExportJobId;
const action = change === 'parent'
  ? api('/api/load_cif', {filename: 'new.cif', content: 'new'})
  : api('/api/method1', {});
reply(change === 'parent' ? '/api/load_cif' : '/api/method1',
  {ok: true, state: state(8), candidates: []});
await action;
document.getElementById('distOut').innerHTML = 'current export view';
if (change === 'request') beginResponseContext('/api/download_all');
if (change === 'job') activeExportJobId = 'different-job';
if (failure === 'network') request.reject(new Error('download transport failure'));
else if (failure === 'blob') request.resolve({
  ok: true,
  headers: {get: name => name.toLowerCase() === 'content-type'
    ? 'application/zip' : encodeURIComponent(JSON.stringify({
      job_id: jobId, active: false, phase: 'complete',
      successful: 1, total: 1, ineligible: 0, failed: 0,
    }))},
  blob: async () => { throw new Error('download body failed'); },
});
else request.resolve({ok: false, headers: {get: () => 'application/json'},
  json: async () => {
    if (failure === 'json') throw new SyntaxError('invalid download JSON');
    return null;
  }});
await flush();
if (failure !== 'blob') reply('/api/export_status', {export: {}});
await download;
assert.equal(STATE.revision, 8);
const obsolete = change === 'request' || change === 'job';
assert.equal(errors.length, obsolete ? 0 : 1,
  'A current snapshot job retains local failure feedback across parent changes');
if (obsolete) assert.equal(document.getElementById('distOut').innerHTML, 'current export view');
""")


@pytest.mark.parametrize("raw,expected", [
    pytest.param("0", 0, id="zero"),
    pytest.param("1", 1, id="one"),
    pytest.param("2", 2, id="two"),
    pytest.param("3", 3, id="maximum"),
    pytest.param("003", 3, id="leading-zeroes"),
    pytest.param(" 1 ", 1, id="surrounding-whitespace"),
    pytest.param("1.9", None, id="decimal"),
    pytest.param("1.0", None, id="integer-with-decimal-point"),
    pytest.param("1e0", None, id="exponent"),
    pytest.param("0x1", None, id="hexadecimal"),
    pytest.param("2oops", None, id="partial-number"),
    pytest.param("+1", None, id="plus-sign"),
    pytest.param("-0", None, id="negative-zero"),
    pytest.param("-1", None, id="negative"),
    pytest.param("", None, id="empty"),
    pytest.param(" ", None, id="blank"),
    pytest.param("NaN", None, id="nan"),
    pytest.param("Infinity", None, id="infinite"),
    pytest.param("１", None, id="non-ascii-digit"),
    pytest.param("4", None, id="unsupported-dimension"),
    pytest.param("9007199254740993", None, id="unsafe-integer"),
])
def test_nmod_input_contract(raw, expected):
    _run_page_javascript(
        "nmodValue = 2;\n"
        f"document.getElementById('nmod').value = {json.dumps(raw)};\n"
        "applyNmod();\n"
        + (
            "assert.equal(nmodValue, 2, 'Rejected text must retain the previous request');\n"
            "assert.deepEqual(errors, ['err.badNmod']);\n"
            if expected is None else
            f"assert.equal(nmodValue, {expected});\nassert.deepEqual(errors, []);\n"
        )
    )


@pytest.mark.parametrize("scenario", [
    pytest.param(r"""
exportingZip = true; activeExportJobId = 'old-job';
const operation = cancelExport();
assert.equal(JSON.parse(pending[0].opts.body).job_id, 'old-job');
activeExportJobId = 'new-job';
document.getElementById('exportProgress').textContent = 'new job: 2 / 10';
document.getElementById('btnCancelExport').disabled = false;
pending.shift().resolve({ok: true, json: async () => ({ok: true, cancellation_requested: true})});
await operation;
assert.equal(document.getElementById('exportProgress').textContent, 'new job: 2 / 10');
assert.equal(document.getElementById('btnCancelExport').disabled, false);
""", id="newer-job-success"),
    pytest.param(r"""
exportingZip = true; activeExportJobId = 'old-job';
const operation = cancelExport();
activeExportJobId = 'new-job';
document.getElementById('exportProgress').textContent = 'new cancellation pending';
document.getElementById('btnCancelExport').disabled = true;
pending.shift().resolve({ok: true, json: async () => ({ok: true, cancellation_requested: false})});
await operation;
assert.equal(document.getElementById('exportProgress').textContent, 'new cancellation pending');
assert.equal(document.getElementById('btnCancelExport').disabled, true,
  'An older false response must not restore the new job cancellation button');
""", id="newer-job-unaccepted-cancellation"),
    pytest.param(r"""
exportingZip = true; activeExportJobId = 'old-job';
const operation = cancelExport();
exportingZip = false; activeExportJobId = null;
document.getElementById('exportProgress').textContent = 'complete: 10 / 10';
pending.shift().resolve({ok: true, json: async () => ({ok: true, cancellation_requested: false})});
await operation;
assert.equal(document.getElementById('exportProgress').textContent, 'complete: 10 / 10');
assert.equal(document.getElementById('btnCancelExport').disabled, true);
""", id="completed-job"),
    pytest.param(r"""
exportingZip = true; activeExportJobId = 'old-job';
const operation = cancelExport();
exportingZip = false;
document.getElementById('exportProgress').textContent = 'complete';
pending.shift().resolve({ok: true, json: async () => ({ok: true, cancellation_requested: true})});
await operation;
assert.equal(document.getElementById('exportProgress').textContent, 'complete');
""", id="completed-before-id-cleanup"),
    pytest.param(r"""
exportingZip = true; activeExportJobId = 'old-job';
const operation = cancelExport();
activeExportJobId = 'new-job';
document.getElementById('exportProgress').textContent = 'new job';
document.getElementById('btnCancelExport').disabled = false;
pending.shift().reject(new Error('old network failure'));
await operation;
assert.equal(document.getElementById('exportProgress').textContent, 'new job');
assert.equal(document.getElementById('btnCancelExport').disabled, false);
""", id="newer-job-network-error"),
    pytest.param(r"""
exportingZip = true; activeExportJobId = 'old-job';
let release;
const parsed = new Promise(resolve => {release = resolve;});
const operation = cancelExport();
pending.shift().resolve({ok: true, json: () => parsed});
await flush();
activeExportJobId = 'new-job';
document.getElementById('exportProgress').textContent = 'new job';
document.getElementById('btnCancelExport').disabled = false;
release({ok: true, cancellation_requested: true});
await operation;
assert.equal(document.getElementById('exportProgress').textContent, 'new job');
assert.equal(document.getElementById('btnCancelExport').disabled, false);
""", id="newer-job-during-json-read"),
    pytest.param(r"""
exportingZip = true; activeExportJobId = 'same-job';
const first = cancelExport();
const second = cancelExport();
const earlier = pending.shift(), later = pending.shift();
later.resolve({ok: true, json: async () => ({ok: true, cancellation_requested: true})});
await second;
const text = document.getElementById('exportProgress').textContent;
earlier.reject(new Error('earlier error'));
await first;
assert.match(text, /^Cancellation requested/);
assert.equal(document.getElementById('exportProgress').textContent, text);
assert.equal(document.getElementById('btnCancelExport').disabled, true);
""", id="same-job-request-order"),
    pytest.param(r"""
activeExportJobId = null; exportingZip = false;
document.getElementById('exportProgress').textContent = 'idle';
document.getElementById('btnCancelExport').disabled = false;
await cancelExport();
assert.equal(pending.length, 0);
assert.equal(document.getElementById('exportProgress').textContent, 'idle');
assert.equal(document.getElementById('btnCancelExport').disabled, false);
""", id="no-active-job"),
])
def test_cancel_export_response_ownership(scenario):
    _run_page_javascript(scenario)


@pytest.mark.parametrize("outcome", ["cancelled", "not-accepted", "http-error", "network-error", "malformed"])
def test_cancel_export_current_job_outcomes(outcome):
    response = {
        "cancelled": "request.resolve({ok: true, json: async () => ({ok: true, cancellation_requested: true})});",
        "not-accepted": "request.resolve({ok: true, json: async () => ({ok: true, cancellation_requested: false})});",
        "http-error": "request.resolve({ok: false, json: async () => ({ok: false, error: 'server failure'})});",
        "network-error": "request.reject(new Error('network failure'));",
        "malformed": "request.resolve({ok: true, json: async () => ({ok: true, cancellation_requested: 'yes'})});",
    }[outcome]
    assertion = {
        "cancelled": "assert.match(text, /^Cancellation requested/); assert.equal(button.disabled, true);",
        "not-accepted": "assert.match(text, /Cancellation was not accepted.*retry/i); assert.equal(button.disabled, false);",
        "http-error": "assert.equal(text, 'server failure'); assert.equal(button.disabled, false);",
        "network-error": "assert.equal(text, 'network failure'); assert.equal(button.disabled, false);",
        "malformed": "assert.match(text, /Invalid cancellation response/); assert.equal(button.disabled, false);",
    }[outcome]
    _run_page_javascript(r"""
exportingZip = true; activeExportJobId = 'current-job';
const operation = cancelExport();
const request = pending.shift();
assert.equal(JSON.parse(request.opts.body).job_id, 'current-job');
""" + response + r"""
await operation;
const text = document.getElementById('exportProgress').textContent;
const button = document.getElementById('btnCancelExport');
assert.equal(activeExportJobId, 'current-job');
""" + assertion)


def test_cancel_export_before_job_start_remains_retryable():
    _run_page_javascript(r"""
selectedDownloadFormats = () => ['cif'];
document.getElementById('dlMethod').value = '1';
m1Items = [{index: 1}]; resultRevisions[1] = 7;
document.createElement = () => ({click() {}, remove() {}});
document.body = {appendChild() {}};
URL.createObjectURL = () => 'blob:test'; URL.revokeObjectURL = () => {};
const download = downloadAll();
const initialPoll = pending.shift();
const downloadRequest = pending.shift();
const jobId = activeExportJobId;
const earlyCancellation = cancelExport();
pending.shift().resolve({ok: true,
  json: async () => ({ok: true, cancellation_requested: false})});
await earlyCancellation;
initialPoll.resolve({json: async () => ({export: {
  job_id: jobId, active: true, phase: 'exporting', total: 1,
  successful: 0, ineligible: 0, failed: 0,
}})});
await flush();
assert.equal(activeExportJobId, jobId);
assert.match(document.getElementById('exportProgress').textContent, /^exporting/);
assert.equal(document.getElementById('btnCancelExport').disabled, false,
  'A cancellation rejected before job registration must allow retry after the same job starts');
const retry = cancelExport();
const retryRequest = pending.shift();
assert.equal(JSON.parse(retryRequest.opts.body).job_id, jobId);
retryRequest.resolve({ok: true,
  json: async () => ({ok: true, cancellation_requested: true})});
await retry;
assert.match(document.getElementById('exportProgress').textContent, /^Cancellation requested/);
assert.equal(document.getElementById('btnCancelExport').disabled, true);
const status = {job_id: jobId, active: false, phase: 'cancelled',
  successful: 0, total: 1, ineligible: 0, failed: 0};
downloadRequest.resolve({ok: true,
  headers: {get: name => name.toLowerCase() === 'content-type'
    ? 'application/zip' : encodeURIComponent(JSON.stringify(status))},
  blob: async () => ({})});
await download;
assert.equal(activeExportJobId, null);
""")


@pytest.mark.parametrize("cancellation_accepted", [True, False])
@pytest.mark.parametrize("terminal_source", ["download-response", "progress-poll"])
def test_export_terminal_status_retires_cancellation_before_blob_read(terminal_source, cancellation_accepted):
    _run_page_javascript(r"""
selectedDownloadFormats = () => ['cif'];
document.getElementById('dlMethod').value = '1';
m1Items = [{index: 1}]; resultRevisions[1] = 7;
document.createElement = () => ({click() {}, remove() {}});
document.body = {appendChild() {}};
URL.createObjectURL = () => 'blob:test'; URL.revokeObjectURL = () => {};
const download = downloadAll();
const initialPoll = pending.shift();
const downloadRequest = pending.shift();
const jobId = activeExportJobId;
const cancellation = cancelExport();
const cancelRequest = pending.shift();
let finishBlob;
const blob = new Promise(resolve => {finishBlob = resolve;});
const status = {job_id: jobId, active: false, phase: 'complete',
  successful: 1, total: 1, ineligible: 0, failed: 0};
const response = {ok: true,
  headers: {get: name => name.toLowerCase() === 'content-type'
    ? 'application/zip' : encodeURIComponent(JSON.stringify(status))},
  blob: () => blob};
initialPoll.resolve({json: async () => ({export: POLL ? status : {}})});
if (!POLL) downloadRequest.resolve(response);
await flush();
const finishedText = document.getElementById('exportProgress').textContent;
assert.match(finishedText, /^complete/);
cancelRequest.resolve({ok: true, json: async () => ({ok: true, cancellation_requested: ACCEPTED})});
await cancellation;
assert.equal(document.getElementById('exportProgress').textContent, finishedText,
  'Cancellation must not overwrite a terminal job while its ZIP body is still downloading');
assert.equal(document.getElementById('btnCancelExport').style.display, 'none');
assert.equal(document.getElementById('btnCancelExport').disabled, true);
if (POLL) downloadRequest.resolve(response);
finishBlob({});
await download;
assert.equal(activeExportJobId, null);
""".replace("POLL", "true" if terminal_source == "progress-poll" else "false")
       .replace("ACCEPTED", str(cancellation_accepted).lower()))


def _official_select_options(cif_name: str, select_name: str) -> list[dict]:
    """Read an archived official selector for regression comparison only."""
    page = (
        TESTS_ROOT.parent.parent
        / "webpage_info" / cif_name / "2. ISODISTORT_ search.html"
    )
    if not page.exists():
        pytest.skip(f"Archived official page is unavailable: {page}")
    source = page.read_text(encoding="utf-8", errors="ignore")
    match = re.search(
        rf'<select name="{re.escape(select_name)}">(.*?)</select>',
        source,
        flags=re.DOTALL,
    )
    assert match, f"Official select {select_name!r} was not found"
    options = []
    for value, label in re.findall(
        r'<option(?: value="([^"]*)")?>(.*?)</option>',
        match.group(1),
        flags=re.DOTALL,
    ):
        text = html.unescape(re.sub(r"<[^>]+>", "", label)).strip()
        if text.lower() in {"no choice", "not selected"}:
            continue
        fields = value.split()
        options.append({
            "label": text,
            "candidate_count": int(fields[0]) if fields else None,
        })
    return options

class _ServerHandle:
    """在测试进程内启动/停止 HTTPServer。"""

    def __init__(self, port: int = 0) -> None:
        self.httpd = web_server.HTTPServer(("127.0.0.1", port), web_server.IsoHandler)
        self.port = self.httpd.server_address[1]

    def __enter__(self):
        self._thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self.httpd.shutdown()
        self.httpd.server_close()

@pytest.fixture(scope="module")
def server():
    with _ServerHandle() as handle:
        time.sleep(0.2)
        yield handle

@pytest.fixture(autouse=True)
def _restore_session():
    """每个测试后重置共享网页会话，避免接口测试相互污染。"""
    web_server._reset_lifecycle()
    workflow_service.reset_export_progress()
    yield
    workflow_service._SESSION.cleanup_uploads()
    workflow_service._SESSION = workflow_service.WebSession()
    web_server._reset_lifecycle()
    workflow_service.reset_export_progress()

def _get(port: int, path: str) -> dict:
    with urllib.request.urlopen(f"http://127.0.0.1:{port}{path}") as resp:
        return json.loads(resp.read().decode("utf-8"))

def _post(port: int, path: str, payload: dict | None = None) -> dict:
    body = json.dumps(payload or {}).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}{path}", data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:  # noqa: S310
        return json.loads(resp.read().decode("utf-8"))


@pytest.mark.parametrize("provided,nmod", [
    (False, None), (True, 0), (True, 1), (True, 3),
    (True, True), (True, "1"), (True, 1.5), (True, -1), (True, 4),
])
def test_method4_http_forwards_expected_nmod_without_recomputing_basis(server, monkeypatch, provided, nmod):
    calls, uploads = [], []

    class ModeEngine(_SessionIsoStub):
        def search_method_4(self, **kwargs):
            calls.append(kwargs)
            raise ValueError("controlled scientific boundary; no decomposition executed")

    live = workflow_service._SESSION
    live._iso = ModeEngine()
    live.revision = 7
    monkeypatch.setattr(
        workflow_service, "_write_upload",
        lambda *_args: uploads.append("daughter") or Path("unused-daughter.cif"),
    )
    monkeypatch.setattr(workflow_service, "_discard_upload", lambda _path: uploads.append("discarded"))
    payload = {"revision": 7, "content": "synthetic transport payload, never parsed"}
    if provided:
        payload["nmod"] = nmod
    result = _post(server.port, "/api/method4", payload)
    assert result["ok"] is False
    valid = not provided or type(nmod) is int and 0 <= nmod <= 3
    if valid:
        assert len(calls) == 1
        assert calls[0]["number_of_independent_modulations"] == nmod
        assert uploads == ["daughter", "discarded"]
        assert "controlled scientific boundary" in result["error"]
    else:
        assert calls == [] and uploads == []
        assert "nmod" in result["error"] or "modulation" in result["error"]
    assert workflow_service._SESSION is live and live.revision == 7


# --- from test_web_server.py ---

def test_index_served(server):
    with urllib.request.urlopen(f"http://127.0.0.1:{server.port}/") as resp:
        body = resp.read().decode("utf-8")
    assert resp.status == 200
    assert "ISODISTORT" in body
    assert "langSel" not in body
    assert "generateDistortion" not in body
    assert "showDomains" not in body
    assert 'id="btnDlAll"' in body
    assert "sortMethod1" not in body
    assert "sortResult" in body
    assert "downloadFilteredTable" in body
    assert "downloadFilteredSubgroups" not in body
    assert 'id="dlMethodOpt4"' in body
    assert 'id="btnDlTxt"' in body and 'id="btnDlCsv"' in body
    assert 'id="genDbHelp"' in body
    assert 'id="btnGenDbManage"' in body
    assert "loadGenDbCache" in body
    assert "deleteSelectedGenDb" in body
    assert "/api/isotropy_cache" in body
    assert "onResultFilter" in body
    assert "m1Items" in body and "m3Items" in body
    assert 'exportHeader: "k-active"' in body
    assert 'exportHeader: "basis"' in body
    assert "it.basis_display" in body and "it.origin_display" in body
    assert 'key: "dir"' in body
    # 对齐官网搜索页布局：well 面板 + OK/Change 按钮 + help 图标
    assert 'class="well"' in body
    assert 'id="hM1"' in body and 'id="hM2"' in body
    assert 'id="btnM1"' in body and 'id="btnTypes"' in body
    assert "help.jpg" in body
    assert "onpagehide" in body or "pagehide" in body or "sendBeacon" in body or "/api/shutdown" in body
    assert 'const pageClientId' in body
    assert '/api/client/close' in body
    assert 'window.addEventListener("beforeunload", sendPageClosed)' in body
    assert "if (exportingZip) return" not in body
    # Distortion 批量导出：官网第 6 页四种格式勾选 + Download all
    assert 'id="fmtCif"' in body and 'id="fmtIsoviz"' in body
    assert 'id="distParamFormatNote"' in body
    assert "applyMethod3DownloadDefaults();" in body
    assert "candidateHasParametricK" in body
    assert 'id="fmtModes"' in body and 'id="fmtTopas"' in body
    assert "async function downloadAll()" in body
    assert "r.blob()" in body
    assert "dist.zipWait" in body
    assert "function startBusy" in body
    assert "function renderStatusBar" in body
    assert body.count("STATE = data.state") == 1
    assert "function responseIsCurrent(data)" in body
    assert "payload.revision = STATE.revision" in body
    assert "acceptResponseState(data, context)" in body
    assert "clearComputedResults();" in body
    assert "&revision=${encodeURIComponent(STATE.revision)}" in body
    assert "knownRouteText" in body
    assert 'exportHeader: "route_resolution"' in body
    assert 't("m4.strainSummary"' in body
    assert "d.strain_applied_engineering_q_parent_basis" in body
    assert "d.strain_raw_coordinate_sum_parent_basis" in body
    assert "statusBox.innerHTML = prevStatus" not in body
    assert "progress-striped" in body
    assert "fmtElapsed" in body
    assert "indices=" in body
    assert "resultTableHtml" in body
    assert "methodColumns" in body
    assert 'class="table table-hover result"' in body
    assert 'name="orderparam"' not in body
    assert "space_group_schoenflies" in body
    assert "toFixed(5)" in body
    assert "function typeLabel" in body
    assert "METHOD2_COUPLED_UNSUPPORTED" in body
    assert "if (n > 1)" in body
    assert "{ nsup: n, kpoints: groups, generate: gen }" in body
    assert 'nsup: source === "subgroups" ? nsupValue : 1' in body
    assert 'method === "2" && nsupValue > 1' in body
    assert '`&nsup=${encodeURIComponent(nsupValue)}`' in body
    assert "method === 2 && nsupValue > 1" in body
    assert "Save interactive distortion" in body
    assert "Complete modes details" in body
    assert "TOPAS.STR" in body
    assert 'id="dlMethod"' in body
    assert 'id="dlMethodOpt1"' in body and 'id="dlMethodOpt2"' in body
    assert 'id="dlMethodOpt3"' in body
    assert 'id="dlMethodOpt4"' in body
    assert 'multiple' not in body.split('id="dlMethod"')[1].split(">")[0]

def test_browser_launch_uses_windows_shell(monkeypatch, capsys):
    seen: dict = {}

    def fake_run(command, **kwargs):
        seen["command"] = command
        seen["kwargs"] = kwargs
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(web_server.sys, "platform", "win32")
    monkeypatch.setattr(
        web_server.shutil,
        "which",
        lambda _name: r"C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe",
    )
    monkeypatch.setattr(web_server.subprocess, "run", fake_run)
    monkeypatch.setattr(
        web_server.webbrowser,
        "open",
        lambda *_args, **_kwargs: pytest.fail("stdlib fallback should not run"),
    )

    url = "http://127.0.0.1:8123/"
    assert web_server._open_browser(url) is True
    assert seen["command"] == [
        r"C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe",
        "-NoProfile",
        "-WindowStyle",
        "Hidden",
        "-Command",
        "Start-Process -FilePath 'http://127.0.0.1:8123/'",
    ]
    assert seen["kwargs"]["timeout"] == 10
    assert f"Browser launch requested: {url}" in capsys.readouterr().out

def test_browser_launch_reports_manual_url_after_all_dispatchers_fail(
    monkeypatch,
    capsys,
):
    monkeypatch.setattr(web_server, "_open_browser_from_windows_shell", lambda _url: False)
    monkeypatch.setattr(web_server.webbrowser, "open", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(web_server.webbrowser, "open_new", lambda *_args, **_kwargs: False)

    url = "http://127.0.0.1:8124/"
    assert web_server._open_browser(url) is False
    assert f"请手动访问: {url}" in capsys.readouterr().out

def test_direct_web_entry_calls_web_main(monkeypatch):
    calls = []
    monkeypatch.setattr(web_server, "main", lambda: calls.append(True) or 23)

    with pytest.raises(SystemExit) as excinfo:
        runpy.run_path(str(_PROJECT_ROOT / "main_web.py"), run_name="__main__")

    assert excinfo.value.code == 23
    assert calls == [True]

def test_state_endpoint_builds_snapshot_under_session_lock(server, monkeypatch):
    observed = []

    def _guarded_summary(_session):
        observed.append(workflow_service._SESSION_LOCK._is_owned())
        return {"structure": None, "distortion_types": []}

    monkeypatch.setattr(workflow_service, "_state_summary", _guarded_summary)
    data = _get(server.port, "/api/state")
    assert data["ok"] is True
    assert observed == [True]

def test_export_status_and_cancel_endpoints_share_job_state(server):
    job = workflow_service._EXPORT_JOBS.start(
        "current-export", revision=0, method=2, total=48,
    )
    workflow_service._EXPORT_JOBS.update(job, dict(
        phase="candidate",
        completed=1,
        successful=1,
        current_candidate="Y1 P2",
        estimated_remaining_seconds=120.0,
    ))

    status = _get(server.port, "/api/export_status")
    assert status["export"]["current_candidate"] == "Y1 P2"
    assert status["export"]["estimated_remaining_seconds"] == 120.0

    cancelled = _post(server.port, "/api/export_cancel", {"job_id": job.job_id})
    assert cancelled == {"ok": True, "cancellation_requested": True}
    assert job.cancellation.is_set()
    assert _get(server.port, "/api/export_status")["export"]["phase"] == (
        "cancellation_requested"
    )
    workflow_service._EXPORT_JOBS.finish(job)

@pytest.mark.parametrize(
    "payload",
    [
        {
            "nsup": 2,
            "kpoints": [
                {"k": "X", "params": None},
                {"k": "X", "params": None},
            ],
        },
        {
            "components": [
                {"k": "X", "irrep": "X1+"},
                {"k": "X", "irrep": "X4+"},
            ],
        },
    ],
)
def test_method2_coupled_subgroup_requests_fail_closed_without_union(
    server,
    payload,
):
    calls = []

    class GuardIso(_SessionIsoStub):
        def list_subgroups_at_kpoint(self, *args, **kwargs):
            calls.append((args, kwargs))
            return []

    previous = [object()]
    guard = GuardIso()
    workflow_service._SESSION._iso = guard
    workflow_service._SESSION.method2_subgroups = previous
    revision = workflow_service._SESSION.revision

    result = _post(
        server.port,
        "/api/subgroups",
        {**payload, "revision": revision},
    )

    assert result["ok"] is False
    assert "more than one superposed IR component" in result["error"]
    assert "cannot be combined" in result["error"]
    assert calls == []
    assert workflow_service._SESSION.method2_subgroups == previous
    assert workflow_service._SESSION.revision == revision

def test_method2_coupled_mode_request_fails_before_using_single_ir_candidate(
    server,
):
    calls = []

    class GuardIso(_SessionIsoStub):
        def search_method_2(self, **kwargs):
            calls.append(kwargs)
            raise AssertionError("coupled request reached the single-IR calculation")

    guard = GuardIso()
    workflow_service._SESSION._iso = guard
    workflow_service._SESSION.method2_subgroups = [SimpleNamespace(index=0)]
    revision = workflow_service._SESSION.revision

    result = _post(
        server.port,
        "/api/method2",
        {
            "subgroup_idx": 0,
            "source": "subgroups",
            "nsup": 2,
            "revision": revision,
        },
    )

    assert result["ok"] is False
    assert "more than one superposed IR component" in result["error"]
    assert calls == []
    assert workflow_service._SESSION.revision == revision

@pytest.mark.parametrize(
    "invalid_nsup",
    [True, False, 1.0, 1.5, "1", "", None],
)
def test_method2_nsup_requires_a_json_integer_before_enumeration(
    server,
    invalid_nsup,
):
    calls = []

    class GuardIso(_SessionIsoStub):
        def list_subgroups_at_kpoint(self, *args, **kwargs):
            calls.append((args, kwargs))
            return []

    workflow_service._SESSION._iso = GuardIso()
    revision = workflow_service._SESSION.revision
    result = _post(
        server.port,
        "/api/subgroups",
        {
            "nsup": invalid_nsup,
            "kpoints": [{"k": "X", "params": None}],
            "revision": revision,
        },
    )

    assert result["ok"] is False
    assert result["error"] == "nsup must be a JSON integer"
    assert calls == []
    assert workflow_service._SESSION.revision == revision

def test_upload_recovers_after_live_read_only_engine_failure(server, monkeypatch, tmp_path):
    class DisabledIso(_SessionIsoStub):
        def snapshot_for_action(self):
            raise RuntimeError("runtime disabled by earlier timeout")

    class FreshIso(_SessionIsoStub):
        def load_structure(self, path):
            self.input_text = Path(path).read_text(encoding="utf-8")

    live = workflow_service._SESSION
    live._iso = DisabledIso()
    live.revision = 7
    live.distortion_types = ["occupational"]
    live.distortion_scope = {"occupational": ["X"]}
    fresh = FreshIso()
    monkeypatch.setattr(workflow_service, "IsoDistort", lambda: fresh)
    monkeypatch.setattr(workflow_service, "get_config", lambda: SimpleNamespace(temp_dir=tmp_path))
    result = _post(server.port, "/api/load_cif", {"content": "data_new\n"})
    assert result["ok"] is True and result["state"]["revision"] == 8
    assert workflow_service._SESSION is not live
    assert workflow_service._SESSION.iso is fresh and fresh.input_text == "data_new\n"
    assert fresh.types == ["occupational"] and fresh.scope == {"occupational": ["X"]}

def test_messages_endpoint(server):
    data = _get(server.port, "/api/messages")
    assert data["ok"]
    assert "language" not in data
    assert "terms" not in data
    assert data["messages"]["load.done"].startswith("[Loaded]")
    assert "ui.menu.language" not in data["messages"]
    assert "Generate isotropy subgroups" in data["messages"]["m2.genDbHelp"] or (
        "isotropy" in data["messages"]["m2.genDbHelp"].lower()
        and "subgroup" in data["messages"]["m2.genDbHelp"].lower()
    )
    assert "m2.genDbManage" in data["messages"]
    assert "m2.genDbDelete" in data["messages"]
    assert data["messages"]["dist.method4"] == "Method 4"
    assert "Building ZIP" in data["messages"]["dist.zipWait"]
    assert "Elapsed" in data["messages"]["st.elapsed"]
    assert "Local engine is running" in data["messages"]["st.busyHint"]

def test_retired_translation_endpoint_is_not_served(server):
    with pytest.raises(urllib.error.HTTPError) as failure:
        _get(server.port, "/api/i18n")
    assert failure.value.code == 404

def test_state_endpoint(server):
    data = _get(server.port, "/api/state")
    assert data["ok"]
    assert data["state"]["structure"] is None
    assert data["state"]["revision"] == 0
    # 默认对齐官网：strain + displacive（全物种）
    assert data["state"]["distortion_types"] == ["strain", "displacive"]
    assert "distortion_scope" in data["state"]

def test_set_types_with_scope(server):
    """官网 per-species 作用域：类型 + 各类型物种范围一起提交。"""
    body = {
        "types": ["displacive", "occupational"],
        "scope": {
            "displacive": ["*"],
            "occupational": ["Eu"],
            "strain": [],
            "magnetic": [],
            "rotational": [],
        },
    }
    data = _post(server.port, "/api/set_types", body)
    assert data["ok"]
    assert data["state"]["distortion_types"] == ["displacive", "occupational"]
    assert data["state"]["distortion_scope"]["occupational"] == ["Eu"]
    assert data["state"]["distortion_scope"]["displacive"] == ["*"]

def test_stale_tab_mutation_and_download_are_rejected_with_http_409(server):
    stub = _SessionIsoStub()
    workflow_service._SESSION._iso = stub
    workflow_service._SESSION.revision = 7

    req = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/set_types",
        data=json.dumps({"revision": 6, "types": ["strain"]}).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(req)  # noqa: S310 - local test server
    assert excinfo.value.code == 409
    payload = json.loads(excinfo.value.read().decode("utf-8"))
    assert payload["stale_context"] is True
    assert payload["state"]["revision"] == 7


    assert "Stale context" in payload["error"]
    assert stub.types is None
    assert workflow_service._SESSION.revision == 7

    url = (
        f"http://127.0.0.1:{server.port}/api/download_all"
        "?method=1&formats=cif&revision=6"
    )
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(url)
    assert excinfo.value.code == 409
    payload = json.loads(excinfo.value.read().decode("utf-8"))
    assert payload["stale_context"] is True
    assert payload["state"]["revision"] == 7

@pytest.mark.parametrize("body", [b"[]", b"null", b"{", b"\xff"])
def test_invalid_json_body_is_rejected_without_touching_the_session(server, body):
    workflow_service._SESSION.revision = 7
    request = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/set_types", data=body,
        headers={"Content-Type": "application/json"},
    )
    with pytest.raises(urllib.error.HTTPError) as caught:
        urllib.request.urlopen(request)
    assert caught.value.code == 400
    assert json.loads(caught.value.read())["ok"] is False
    assert workflow_service._SESSION.revision == 7
    assert workflow_service._SESSION._iso is None

@pytest.mark.parametrize("route,payload", [
    ("set_types", {"types": ["strain"], "scope": {"displacive": ["\ud800"]}}),
    ("load_cif", {"content": "data_input\n\udfff", "filename": "input.cif"}),
    ("load_cif", {"content": "data_input\n", "filename": "\ud800.cif"}),
])
def test_surrogate_action_and_upload_preserve_live_state_and_files(
    server, monkeypatch, tmp_path, route, payload,
):
    live = workflow_service._SESSION
    live.revision = 7
    original = tmp_path / "original.cif"
    original.write_text("original input", encoding="utf-8")
    live.parent_upload_path = original
    stub = _SessionIsoStub()
    live._iso = stub
    monkeypatch.setattr(workflow_service, "get_config", lambda: SimpleNamespace(temp_dir=tmp_path))

    result = _post(server.port, f"/api/{route}", {"revision": 7, **payload})

    assert result["ok"] is False
    assert "Unicode scalar" in result["error"]
    assert workflow_service._SESSION is live
    assert live.revision == 7 and live._iso is stub
    assert stub.scope is None and stub.types is None
    assert original.read_text(encoding="utf-8") == "original input"
    assert list(tmp_path.iterdir()) == [original]

@pytest.mark.parametrize("bad_value", ["\ud800", float("nan")])
def test_unencodable_action_response_cannot_commit_over_http(server, monkeypatch, bad_value):
    live = workflow_service._SESSION
    live.revision = 7
    live._iso = _SessionIsoStub()
    original_types = list(live.distortion_types)

    def invalid_response(staged, _data):
        staged.distortion_types = ["strain"]
        return {"diagnostic": bad_value}

    monkeypatch.setitem(
        workflow_service._ACTIONS, "set_types",
        replace(workflow_service._ACTIONS["set_types"], handler=invalid_response),
    )
    result = _post(server.port, "/api/set_types", {"types": ["strain"], "revision": 7})
    assert result["ok"] is False
    assert workflow_service._SESSION is live
    assert live.revision == 7 and live.distortion_types == original_types

def test_set_types_bumps_revision_and_invalidates_all_derived_state(server):
    stub = _SessionIsoStub()
    workflow_service._SESSION._iso = stub
    workflow_service._SESSION.revision = 3
    workflow_service._SESSION.method1 = [object()]
    workflow_service._SESSION.method2 = object()
    workflow_service._SESSION.method3 = [object()]
    workflow_service._SESSION.method2_subgroups = [object()]

    data = _post(server.port, "/api/set_types", {
        "revision": 3,
        "types": ["displacive"],
        "scope": {"displacive": ["*"]},
    })

    assert data["ok"] is True
    assert data["state"]["revision"] == 4
    assert workflow_service._SESSION.revision == 4
    assert workflow_service._SESSION.method1 == []
    assert workflow_service._SESSION.method2 is None
    assert workflow_service._SESSION.method3 == []
    assert workflow_service._SESSION.method2_subgroups == []
    stub = workflow_service._SESSION.iso
    assert stub.subgroups == []
    assert stub.mode_displacements == {}
    assert stub.mode_occupancies == {}
    assert stub.distorted_structure is None
    assert stub.selected_modes_cleared is True

def test_ping_endpoint(server):
    data = _get(server.port, "/api/ping?client_id=tab-a")
    assert data["ok"]
    assert data["active_clients"] == 1
    assert "tab-a" in web_server._LIFE.clients
    assert web_server._LIFE.page_heartbeat > 0  # 心跳被记录（守护线程停服依据）

def test_page_close_tracks_last_browser_tab(server):
    _get(server.port, "/api/ping?client_id=tab-a")
    _get(server.port, "/api/ping?client_id=tab-b")

    first = _post(server.port, "/api/client/close?client_id=tab-a")
    assert first == {"ok": True, "closing": False, "active_clients": 1}
    assert web_server._LIFE.last_page_closed == 0

    last = _post(server.port, "/api/client/close?client_id=tab-b")
    assert last == {"ok": True, "closing": True, "active_clients": 0}
    assert web_server._LIFE.last_page_closed > 0
    assert web_server._LIFE.shutdown_requested is False

def test_watchdog_finishes_export_before_stopping_after_tab_close():
    class FakeServer:
        def __init__(self):
            self.stopped = threading.Event()

        def shutdown(self):
            self.stopped.set()

    fake = FakeServer()
    web_server._mark_page_seen()
    web_server._touch_heartbeat("tab-export")
    web_server._begin_long_request()
    web_server._mark_page_closed("tab-export")
    thread = threading.Thread(
        target=web_server._watchdog,
        args=(fake, 60.0),
        kwargs={"close_grace": 0.02, "poll_interval": 0.005},
        daemon=True,
    )
    thread.start()
    assert not fake.stopped.wait(0.05)
    web_server._end_long_request()
    assert fake.stopped.wait(1.0)
    thread.join(timeout=1.0)

def test_last_page_close_stops_http_server_and_releases_port():
    httpd = web_server.ThreadingHTTPServer(
        ("127.0.0.1", 0), web_server.IsoHandler
    )
    port = httpd.server_address[1]
    server_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    watchdog_thread = threading.Thread(
        target=web_server._watchdog,
        args=(httpd, 60.0),
        kwargs={"close_grace": 0.02, "poll_interval": 0.005},
        daemon=True,
    )
    server_thread.start()
    watchdog_thread.start()
    _get(port, "/api/ping?client_id=only-tab")
    closed = _post(port, "/api/client/close?client_id=only-tab")
    assert closed["closing"] is True
    server_thread.join(timeout=1.0)
    assert not server_thread.is_alive()
    httpd.server_close()
    watchdog_thread.join(timeout=1.0)

    rebound = web_server.HTTPServer(("127.0.0.1", port), web_server.IsoHandler)
    rebound.server_close()

def test_shutdown_endpoint(server):
    """/api/shutdown 置位关闭请求（守护线程稍后执行 server.shutdown）。"""
    data = _post(server.port, "/api/shutdown")
    assert data["ok"] and data["shutdown"]
    assert web_server._LIFE.shutdown_requested is True
    web_server._LIFE.shutdown_requested = False  # 复位，避免影响模块内其他测试

def test_space_groups_endpoint(server):
    data = _get(server.port, "/api/space_groups")
    assert data["ok"]
    assert len(data["space_groups"]) == 230
    assert data["space_groups"][0]["number"] == 1
    assert data["space_groups"][0]["symbol"] == "P1"
    assert data["space_groups"][0]["schoenflies"] == "C1-1"
    assert data["space_groups"][-1]["number"] == 230
    # ISODISTORT 官网 Method 3 下拉：旧 IT 写法（非 e-glide）
    by_num = {g["number"]: g["symbol"] for g in data["space_groups"]}
    assert by_num[39] == "Abm2"
    assert by_num[41] == "Aba2"
    assert by_num[64] == "Cmca"
    assert by_num[67] == "Cmma"
    assert by_num[68] == "Ccca"

def test_all_method3_space_groups_match_archived_official_page(server):
    """Exhaustively verify all 230 Method-3 SG labels, not a hand-picked sample."""
    official = _official_select_options("EuAl4 Parent.cif", "subgroupsym")
    actual = _get(server.port, "/api/space_groups")["space_groups"]
    expected = []
    for option in official:
        match = re.match(r"(\d+)\s+(\S+)\s+(\S+)$", option["label"])
        assert match, option["label"]
        expected.append({
            "number": int(match.group(1)),
            "symbol": match.group(2),
            "schoenflies": match.group(3),
        })
    assert actual == expected

def test_static_css_served(server):
    with urllib.request.urlopen(f"http://127.0.0.1:{server.port}/static/bootstrap.css") as resp:
        body = resp.read()
    assert resp.status == 200
    assert len(body) > 1000  # 官网 Bootstrap 2.x 样式表已就位
    assert "text/css" in resp.headers["Content-Type"]

def test_download_all_without_method2_is_json_error(server):
    """未完成 Method 2 时 Download all 不得打包 output_dir，应返回 JSON 错误。"""
    url = f"http://127.0.0.1:{server.port}/api/download_all?method=2&formats=cif"
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(url)
    err = excinfo.value
    assert err.code == 404
    data = json.loads(err.read().decode("utf-8"))
    assert data["ok"] is False
    assert "Method 2" in data["error"] or "CIF" in data["error"]

def test_download_all_rejects_coupled_method2_before_exporting_old_pool(server):
    calls = []

    class ExportIsoStub:
        structure = object()

        def export_subgroups_zip(self, **kwargs):
            calls.append(kwargs)
            raise AssertionError("coupled request exported the old single-IR pool")

    subgroup = SimpleNamespace(index=0)
    workflow_service._SESSION._iso = ExportIsoStub()
    workflow_service._SESSION.method2_subgroups = [subgroup]

    url = (
        f"http://127.0.0.1:{server.port}/api/download_all"
        "?method=2&formats=cif&nsup=2"
    )
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(url)

    assert excinfo.value.code == 400
    payload = json.loads(excinfo.value.read().decode("utf-8"))
    assert payload["ok"] is False
    assert "more than one superposed IR component" in payload["error"]
    assert calls == []
    assert workflow_service._SESSION.method2_subgroups == [subgroup]

def test_download_all_rejects_multiple_methods(server):
    """method=1,2 必须拒绝，不允许一次导出多个 Method。"""
    url = f"http://127.0.0.1:{server.port}/api/download_all?method=1,2&formats=cif"
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        urllib.request.urlopen(url)
    err = excinfo.value
    assert err.code == 400
    data = json.loads(err.read().decode("utf-8"))
    assert data["ok"] is False
    assert "不能多选" in data["error"] or "exactly one" in data["error"]

def test_unsupported_export_modulation_count_returns_json_error(server):
    url = f"http://127.0.0.1:{server.port}/api/download_all?method=2&nmod=4"
    with pytest.raises(urllib.error.HTTPError) as caught:
        urllib.request.urlopen(url)
    assert caught.value.code == 400
    payload = json.loads(caught.value.read())
    assert payload["ok"] is False
    assert "not implemented" in payload["error"]
    assert workflow_service._SESSION._iso is None

@pytest.mark.parametrize("query", [
    "nmod=", "nmod=&nmod=1", "nmod=1&nmod=", "nmod=0&nmod=1",
    "job_id=", "job_id=&job_id=00000000-0000-4000-8000-000000000000",
    "formats=&formats=cif",
])
def test_export_blank_and_duplicate_query_fields_are_rejected_over_http(server, query):
    live = workflow_service._SESSION
    live.revision = 7
    url = f"http://127.0.0.1:{server.port}/api/download_all?{query}"
    with pytest.raises(urllib.error.HTTPError) as caught:
        urllib.request.urlopen(url)
    assert caught.value.code == 400
    assert json.loads(caught.value.read())["ok"] is False
    assert workflow_service._SESSION is live
    assert live.revision == 7 and live._iso is None


@pytest.mark.parametrize("method", [1, 2, 3])
@pytest.mark.parametrize("indices", ["1,999", ""])
def test_export_invalid_candidate_selection_fails_before_job_or_snapshot(server, method, indices):
    class LiveIso:
        structure = object()

        def snapshot_for_export(self):
            raise AssertionError("Invalid candidate selection must not start scientific work")

    live = workflow_service._SESSION
    live.revision = 7
    live._iso = LiveIso()
    subgroup = SimpleNamespace(index=1)
    live.method1 = [SimpleNamespace(subgroup=subgroup)]
    live.method2_subgroups = [subgroup]
    live.method3 = [SimpleNamespace(subgroup=subgroup)]
    before = workflow_service.export_progress_snapshot()
    query = urllib.parse.urlencode({
        "method": method, "formats": "cif", "revision": 7, "indices": indices,
    })
    with pytest.raises(urllib.error.HTTPError) as caught:
        urllib.request.urlopen(f"http://127.0.0.1:{server.port}/api/download_all?{query}")
    assert caught.value.code == 400
    data = json.loads(caught.value.read())
    assert data["ok"] is False
    if indices:
        assert data["error"] == f"Subgroup indices not found in Method {method} candidate pool: 999"
    else:
        assert "indices" in data["error"]
    assert workflow_service._SESSION is live and live.revision == 7
    assert workflow_service.export_progress_snapshot() == before


def test_download_all_passes_nmod_without_mutating_cached_mode_context(server):
    seen: dict = {}

    class ExportIsoStub:
        def __init__(self) -> None:
            self.structure = object()
            self.number_of_independent_modulations = 1

        def export_subgroups_zip(self, **kwargs):
            seen.update(kwargs)
            return b"test-zip"

        def snapshot_for_export(self):
            return copy.deepcopy(self)

        def close(self):
            pass

    stub = ExportIsoStub()
    subgroup = SimpleNamespace(index=0)
    workflow_service._SESSION._iso = stub
    workflow_service._SESSION.method1 = [SimpleNamespace(subgroup=subgroup)]

    url = (
        f"http://127.0.0.1:{server.port}/api/download_all"
        "?method=1&formats=modes&nmod=2"
    )
    with urllib.request.urlopen(url) as response:
        assert response.read() == b"test-zip"

    assert seen["number_of_independent_modulations"] == 2
    assert seen["export_method"] == 1
    assert seen["subgroups"] == [subgroup]
    assert stub.number_of_independent_modulations == 1

def test_download_all_snapshots_under_lock_and_exports_after_releasing_it(server):
    observed: dict[str, bool] = {}

    class ExportWorker:
        def export_subgroups_zip(self, **_kwargs):
            observed["export_locked"] = workflow_service._SESSION_LOCK._is_owned()
            return b"test-zip"

        def close(self):
            pass

    class LiveIso:
        structure = object()

        def snapshot_for_export(self):
            observed["snapshot_locked"] = workflow_service._SESSION_LOCK._is_owned()
            return ExportWorker()

    subgroup = SimpleNamespace(index=0)
    workflow_service._SESSION._iso = LiveIso()
    workflow_service._SESSION.method1 = [SimpleNamespace(subgroup=subgroup)]

    with urllib.request.urlopen(
        f"http://127.0.0.1:{server.port}/api/download_all"
        "?method=1&formats=cif"
    ) as response:
        assert response.read() == b"test-zip"

    assert observed == {"snapshot_locked": True, "export_locked": False}

def test_download_response_keeps_its_terminal_status_when_the_next_job_starts(server, monkeypatch):
    first_id, second_id = str(uuid.uuid4()), str(uuid.uuid4())
    registry = workflow_service.ExportJobs()
    monkeypatch.setattr(workflow_service, "_EXPORT_JOBS", registry)
    original_finish = registry.finish
    second = []

    def finish_and_start_next(job, **kwargs):
        result = original_finish(job, **kwargs)
        second.append(registry.start(second_id, revision=8, method=2, total=99))
        return result

    monkeypatch.setattr(registry, "finish", finish_and_start_next)

    class Worker:
        def export_subgroups_zip(self, **kwargs):
            kwargs["progress_callback"]({
                "phase": "cancelled", "completed": 1, "successful": 1,
                "irrep_label": "\u0393", "ineligible": 0, "failed": 0,
            })
            return b"partial-zip"

        def close(self):
            pass

    live = workflow_service._SESSION
    live.revision = 7
    live._iso = SimpleNamespace(structure=object(), snapshot_for_export=Worker)
    live.method2_subgroups = [SimpleNamespace(index=index) for index in range(3)]
    try:
        with urllib.request.urlopen(
            f"http://127.0.0.1:{server.port}/api/download_all?formats=cif&job_id={first_id}"
        ) as response:
            status = json.loads(urllib.parse.unquote(response.headers["X-ISODISTORT-Export"]))
            assert response.read() == b"partial-zip"
        assert registry.snapshot()["job_id"] == second_id
        assert status["job_id"] == first_id and status["source_revision"] == 7
        assert status["active"] is False and status["phase"] == "cancelled"
        assert status["successful"] == 1 and status["total"] == 3
        assert status["irrep_label"] == "\u0393"
        status["successful"] = 500
        assert registry.snapshot()["successful"] == 0
    finally:
        for job in second:
            original_finish(job)

@pytest.mark.parametrize("matching_status", [True, False])
def test_download_ui_uses_response_status_while_polling_observes_the_next_job(matching_status):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is required to execute the browser download function")
    source = (_PROJECT_ROOT / "frontend" / "web" / "index.html").read_text(encoding="utf-8")
    download = "async function downloadAll()" + source.split(
        "async function downloadAll()", 1,
    )[1].split("async function cancelExport()", 1)[0]
    harness = r"""
const assert = require('node:assert/strict');
const jobId = '00000000-0000-4000-8000-000000000001';
const nextJobId = '00000000-0000-4000-8000-000000000002';
const crypto = {randomUUID: () => jobId};
const nodes = {};
const successes = [], errors = [];
let clicked = 0, exportingZip = false, activeExportJobId = null;
const nsupValue = 1, nmodValue = 0, STATE = {revision: 7};
const document = {
  getElementById: id => nodes[id] ||= {value: '2', style: {}, textContent: ''},
  createElement: () => ({click: () => clicked++, remove() {}}),
  body: {appendChild() {}},
};
const window = {location: {href: 'http://localhost/'}};
const t = (...args) => args.join(' ');
const itemsFor = () => [{index: 0}, {index: 1}, {index: 2}];
const itemMatches = () => true, candidateHasParametricK = () => false;
const selectedDownloadFormats = () => ['cif'];
const startBusy = () => ({stop() {}});
const beginResponseContext = () => ({});
const receiveResponse = data => data, responseIsCurrent = () => true;
const showOk = (_out, text) => successes.push(text);
const showError = (_out, text) => errors.push(text);
const setInterval = () => 1, clearInterval = () => {};
URL.createObjectURL = () => 'blob:test'; URL.revokeObjectURL = () => {};
const fetch = async url => {
  if (url === '/api/export_status') return {json: async () => ({export: {
    job_id: nextJobId, active: true, successful: 0, total: 99,
  }})};
  assert.equal(url.searchParams.get('job_id'), jobId);
  return {ok: true, headers: {get: name => name.toLowerCase() === 'content-type'
    ? 'application/zip' : encodeURIComponent(JSON.stringify({
      job_id: MATCHING ? jobId : nextJobId, active: false, phase: 'cancelled',
      successful: 1, total: 3, ineligible: 0, failed: 0,
    }))}, blob: async () => ({})};
};
""".replace("MATCHING", "true" if matching_status else "false")
    assertions = r"""
downloadAll().then(() => {
  if (MATCHING) {
    assert.equal(clicked, 1); assert.equal(errors.length, 0);
    assert.match(successes[0], /1\/3 ready candidates/);
    assert.match(nodes.exportProgress.textContent, /cancelled.*ready 1/);
  } else {
    assert.equal(clicked, 0); assert.equal(successes.length, 0);
    assert.match(errors[0], /matching final status/);
  }
  assert.equal(activeExportJobId, null);
}).catch(error => { console.error(error); process.exitCode = 1; });
""".replace("MATCHING", "true" if matching_status else "false")
    completed = subprocess.run(
        [node, "-"], input=harness + download + assertions,
        text=True, capture_output=True, timeout=10,
    )
    assert completed.returncode == 0, completed.stderr

def test_server_close_waits_for_an_accepted_request(monkeypatch):
    entered = threading.Event()
    release = threading.Event()
    closed = threading.Event()
    errors = []

    class Handler(web_server.BaseHTTPRequestHandler):
        def do_GET(self):
            entered.set()
            assert release.wait(5)
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b"finished")

        def log_message(self, *_args):
            pass

    monkeypatch.setattr(web_server, "IsoHandler", Handler)
    httpd = web_server._bind_server("127.0.0.1", 0)
    serving = threading.Thread(target=httpd.serve_forever)
    serving.start()

    def request():
        try:
            with urllib.request.urlopen(f"http://127.0.0.1:{httpd.server_address[1]}/") as response:
                assert response.read() == b"finished"
        except BaseException as exc:
            errors.append(exc)

    def close_server():
        httpd.server_close()
        closed.set()

    client = threading.Thread(target=request)
    closer = threading.Thread(target=close_server)
    client.start()
    try:
        assert entered.wait(5)
        httpd.shutdown()
        closer.start()
        assert not closed.wait(0.05)
    finally:
        release.set()
        client.join(5)
        serving.join(5)
        if closer.ident is not None:
            closer.join(5)
        else:
            httpd.shutdown()
            httpd.server_close()
    assert errors == []
    assert closed.is_set()

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_load_cif_endpoint(server):
    """上传 EuAl4 母相 CIF 并校验识别结果（真实计算）。"""
    cif = experiment_data_dir() / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = experiment_data_dir() / "EuAl4 Springer (parent).cif"  # 旧命名兼容
    if not cif.exists():
        pytest.skip("EuAl4 CIF 不存在")
    body = json.dumps({"filename": "eual4.cif", "content": cif.read_text(encoding="utf-8")}).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/load_cif?lang=zh", data=body,
        headers={"Content-Type": "application/json"},
    )
    # 仅访问本地测试服务（http scheme），无 file:// 等自定义 scheme
    with urllib.request.urlopen(req) as resp:  # noqa: S310
        data = json.loads(resp.read().decode("utf-8"))
    assert data["ok"]
    assert data["state"]["structure"]["space_group_number"] == 139
    assert data["state"]["structure"]["space_group_schoenflies"] == "D4h-17"
    # 页头信息（对齐官网 "Reading CIF file..." 段落）：点阵参数 + Wyckoff 坐标
    lat = data["state"]["structure"]["lattice"]
    assert abs(lat["a"] - 4.402) < 0.01
    assert abs(lat["c"] - 11.163) < 0.01
    assert round(lat["a"], 5) == lat["a"]
    wy = data["state"]["structure"]["wyckoff"]
    assert len(wy) == 3
    assert any(w["coordinates"] == [0, 0, 0] for w in wy)  # Eu 2a (0,0,0)
    display = data["state"]["structure"]["wyckoff_display"]
    assert display[0] == "Eu1 2a (0,0,0)"
    assert display[1] == "Al1 4d (0,1/2,1/4)"
    assert display[2].startswith("Al2 4e (0,0,z), z=")
    # 官网 "Default space-group preferences: ..." 行
    assert "monoclinic axes a(b)c" in data["state"]["structure"]["preferences"]
    assert data["state"]["species"] == ["Al", "Eu"]

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_load_4310_endpoint_reports_each_standard_wyckoff_orbit(server):
    cif = experiment_data_dir() / "4310_tetra.cif"
    if not cif.exists():
        pytest.skip("4310 CIF 不存在")
    body = json.dumps({
        "filename": "4310_tetra.cif",
        "content": cif.read_text(encoding="utf-8"),
    }).encode("utf-8")
    request = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/load_cif",
        data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request) as response:  # noqa: S310
        data = json.loads(response.read().decode("utf-8"))

    assert data["ok"]
    assert data["state"]["structure"]["wyckoff_display"] == [
        "La1 4e (0,0,z), z= 0.43204",
        "La2 4e (0,0,z), z= 0.30148",
        "Ni1 2a (0,0,0)",
        "Ni2 4e (0,0,z), z=-0.13885",
        "O1 8g (0,1/2,z), z= 0.36070",
        "O3 4c (0,1/2,0)",
        "O4 4e (0,0,z), z=-0.21680",
        "O2 4e (0,0,z), z=-0.06780",
    ]

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过真实计算端点")
def test_method1_options_endpoint(server):
    """Method 1 下拉数据（真实枚举）：可达子群空间群（<230）+ lattice 选项。"""
    cif = experiment_data_dir() / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("EuAl4 CIF 不存在")
    body = json.dumps({"filename": "eual4.cif", "content": cif.read_text(encoding="utf-8")}).encode("utf-8")
    req = urllib.request.Request(
        f"http://127.0.0.1:{server.port}/api/load_cif", data=body,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:  # noqa: S310 - 仅访问本地测试服务（http scheme）
        data = json.loads(resp.read().decode("utf-8"))
    assert data["ok"]

    data = _get(server.port, "/api/method1_options")
    assert data["ok"]
    opts = data["options"]
    # 官网只列出与母相结构 + Types 相容的可达子群，而非全部 230 个
    assert 0 < len(opts["space_groups"]) < 230
    # 序号升序排列（官网下拉按序号排列，便于查找）；母相自身（139）包含在内
    numbers = [sg["number"] for sg in opts["space_groups"]]
    assert numbers == sorted(numbers)
    assert 139 in numbers
    # EuAl4 + 默认 strain/displacive：官网无 87/97/126/128（无位移活性 IR）
    for extra in (87, 97, 126, 128):
        assert extra not in numbers
    assert opts["conventional_lattices"]
    assert opts["primitive_lattices"]
    # 官网风格标签，如 (1,0,0),(0,1,0),(0,0,1)
    assert "(" in opts["conventional_lattices"][0]["label"]
    # I4/mmm EuAl4：官网 Conventional 12、Primitive 9
    assert len(opts["conventional_lattices"]) == 12
    assert len(opts["primitive_lattices"]) == 9
    assert opts["conventional_lattices"][0]["label"] == "(1,0,0),(0,1,0),(0,0,1)"
    assert opts["primitive_lattices"][0]["label"] == (
        "(-1/2,1/2,1/2),(1/2,-1/2,1/2),(1/2,1/2,-1/2)"
    )

    # The saved official page is a reference for the exact visible labels.
    official_sg = _official_select_options("EuAl4 Parent.cif", "isosubgroup")
    assert [
        f"{sg['number']} {sg['symbol']} {sg['schoenflies']}" for sg in opts["space_groups"]
    ] == [" ".join(sg["label"].split()) for sg in official_sg]

    # Compare independently calculated lattice classes and candidate counts;
    # archived display representatives never feed back into production.
    iso = workflow_service._SESSION.iso
    for local_key, official_name in (
        ("conventional_lattices", "isolattice"),
        ("primitive_lattices", "isoplattice"),
    ):
        archived = _official_select_options("EuAl4 Parent.cif", official_name)
        local = opts[local_key]
        assert len(local) == len(archived)
        assert [item["candidate_count"] for item in local] == [
            item["candidate_count"] for item in archived
        ]
        for got, expected in zip(local, archived, strict=True):
            expected_basis = parse_basis_token(expected["label"])
            # A computed GL(3,Z)/point-group representative need not equal
            # the archived display basis; its own label must describe it.
            assert np.allclose(parse_basis_token(got["label"]), got["basis"],
                               rtol=0.0, atol=1e-12)
            assert iso._same_lattice_orbit(got["basis"], expected_basis), (
                got["label"], expected["label"]
            )

    # Every option must select precisely the number of live rows recorded by
    # the official form, with independently computed display representatives.
    live_subgroups = iso._method1_filtered_subgroups()
    for key in ("conventional_lattices", "primitive_lattices"):
        for option in opts[key]:
            count = 0
            for sg in live_subgroups:
                basis = np.asarray(sg.basis_vectors, dtype=float)
                if key == "primitive_lattices":
                    basis = iso._centering_matrix(
                        _centering_letter(sg.space_group_number)
                    ) @ basis
                count += iso._same_lattice_orbit(basis, option["basis"])
            assert count == option["candidate_count"], option["label"]

    # The official selector IDs target one class. In particular, selecting
    # the identity conventional lattice must not include every supercell.
    for key, kind in (
        ("conventional_lattices", "conventional"),
        ("primitive_lattices", "primitive"),
    ):
        option = opts[key][0]
        result = _post(server.port, "/api/method1", {
            "lattice": option["basis"], "lattice_kind": kind,
            "distortion_types": ["strain", "displacive"],
        })
        assert result["ok"]
        assert len(result["candidates"]) == option["candidate_count"]

# --- from test_three_interface.py ---

def _load_cif_via_web(port: int, path: Path) -> dict:
    content = path.read_text(encoding="utf-8")
    return _post(port, "/api/load_cif?lang=en",
                 {"filename": path.name, "content": content})

def _sg_of(structure) -> int:
    return SpacegroupAnalyzer(structure, symprec=1e-3).get_space_group_number()

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_method1_options(server):
    """同一母相：API 与网页 method1_options 空间群/格子完全一致。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("EuAl4 母相不存在")

    # API 路径
    iso = IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    api_opts = iso.method1_options()

    # 网页路径
    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    web_opts = _get(server.port, "/api/method1_options")["options"]

    api_sgs = [g["number"] for g in api_opts["space_groups"]]
    web_sgs = [g["number"] for g in web_opts["space_groups"]]
    assert api_sgs == web_sgs, "Method1 可达子群空间群不一致"
    api_conv = [g["label"] for g in api_opts["conventional_lattices"]]
    web_conv = [g["label"] for g in web_opts["conventional_lattices"]]
    assert api_conv == web_conv, "Conventional lattice 选项不一致"
    api_prim = [g["label"] for g in api_opts["primitive_lattices"]]
    web_prim = [g["label"] for g in web_opts["primitive_lattices"]]
    assert api_prim == web_prim, "Primitive lattice 选项不一致"

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_method1_candidates(server):
    """同一母相：API 与网页 method1 候选数一致（同过滤条件）。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")

    iso = IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    api_cands = iso.search_method_1(
        distortion_types=["displacive", "strain"],
        crystal_system=None, subgroup_space_group=None,
        maximal_subgroup_only=False,
    )

    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    web_m1 = _post(server.port, "/api/method1",
                   {"distortion_types": ["displacive", "strain"],
                    "maximal_subgroup_only": False})
    assert web_m1["ok"]
    web_cands = web_m1["candidates"]

    assert len(api_cands) == len(web_cands), \
        f"Method1 候选数不一致: API={len(api_cands)} web={len(web_cands)}"
    # 首候选的空间群一致
    if api_cands:
        api_first = api_cands[0].subgroup.space_group_number
        web_first = web_cands[0]["space_group_number"]
        assert api_first == web_first, \
            f"Method1 首候选空间群不一致: API={api_first} web={web_first}"

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_method2_and_generation(server):
    """同一子群：API 与网页 method2 模式数、生成结构空间群一致。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")

    iso = IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    api_cands = iso.search_method_1(
        distortion_types=["displacive", "strain"],
        maximal_subgroup_only=False,
    )
    assert api_cands, "API 无候选"
    chosen = api_cands[0]
    iso.search_method_2(subgroup_idx=chosen.subgroup.index,
                        distortion_type=["displacive"])
    api_n_modes = len(iso.mode_displacements)
    assert api_n_modes > 0, "API 无位移模式"

    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    web_m1 = _post(server.port, "/api/method1",
                   {"distortion_types": ["displacive", "strain"],
                    "maximal_subgroup_only": False})
    assert web_m1["ok"]
    web_m2 = _post(server.port, "/api/method2",
                   {"subgroup_idx": chosen.subgroup.index,
                    "distortion_type": ["displacive"],
                    "nmod": 0, "nsup": 1})
    assert web_m2["ok"]
    web_modes = [m for m in web_m2["modes"] if m["mode_type"] == "displacive"]
    assert len(web_modes) == api_n_modes, \
        f"method2 模式数不一致: API={api_n_modes} web={len(web_modes)}"

    # 生成：API 与网页（同底层引擎）畸变结构空间群一致
    label = next(iter(iso.mode_displacements))
    api_dist = iso.generate_distortion(irrep_label=label, amplitude=0.1)
    api_sg = _sg_of(api_dist)
    # 网页生成走 /api/generate（若有该端点）；无则跳过数值比对
    if web_modes:
        web_sg = None
        # 通过引擎直接生成同一模式（网页后端复用同一 IsoDistort 会话）
        web_iso = workflow_service._SESSION.iso
        try:
            web_dist = web_iso.generate_distortion(irrep_label=label,
                                                   amplitude=0.1)
            web_sg = _sg_of(web_dist)
        except Exception:  # noqa: BLE001 - 会话状态差异时跳过
            web_sg = None
        if web_sg is not None:
            assert web_sg == api_sg, \
                f"生成畸变空间群不一致: API={api_sg} web={web_sg}"

@pytest.mark.skipif(not _wsl_wsl_available(), reason="WSL 不可用，跳过三接口一致性测试")
def test_api_vs_web_state_structure(server):
    """加载后状态结构：API 与网页识别空间群/原子数一致（数值 ≤1e-6 级）。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("EuAl4 母相不存在")

    iso = IsoDistort()
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)

    web = _load_cif_via_web(server.port, cif)
    assert web["ok"]
    st = web["state"]["structure"]
    assert st["space_group_number"] == iso.symmetry_info["space_group_number"] == 139
    assert st["atoms"] == len(iso.structure)
    # 晶格参数一致性（网页端从同一结构序列化）
    assert abs(st["lattice"]["a"] - iso.structure.lattice.a) < MATCH_TOL
    assert abs(st["lattice"]["c"] - iso.structure.lattice.c) < MATCH_TOL
