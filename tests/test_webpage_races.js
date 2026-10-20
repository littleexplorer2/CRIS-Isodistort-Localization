// DOM-adapter engineering test, not a full browser or scientific acceptance test.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const {TextDecoder} = require("node:util");
const file = path.join(__dirname, "../cris_rea/frontend/webpage.py");
const script = fs.readFileSync(file, "utf8").match(/<script nonce="__TOKEN__">([\s\S]*?)<\/script>/)[1];
const elements = new Map();
function element(id) {
  if (!elements.has(id)) elements.set(id, {
    value: "", files: [], disabled: true, textContent: "", className: "",
    handlers: {}, addEventListener(event, callback) {this.handlers[event] = callback;},
    appendChild() {},
  });
  return elements.get(id);
}
const submissions = [];
const context = {
  document: {getElementById: element, createElement() {return {};}},
  TextDecoder, Uint8Array,
  btoa: value => Buffer.from(value, "binary").toString("base64"),
  fetch: async (url, options) => {
    if (url === "/api/programs") return {ok: true, json: async () => ({iso: {description: "raw"}})};
    if (url === "/api/iso-query") return {
      ok: true, json: async () => ({program: "iso", input: "PREPARED\nQUIT\n", scope: "query"}),
    };
    assert.equal(url, "/api/run");
    submissions.push(JSON.parse(options.body));
    return {ok: true, json: async () => ({
      evidence: {status: "completed_uninterpreted", run_id: "engineering"},
      run_directory: "engineering", stdout: "", stderr: "",
    })};
  },
};
vm.createContext(context);
vm.runInContext(script, context);
function deferredFile(name, text) {
  let resolve;
  const pending = new Promise(done => {resolve = done;});
  return {
    file: {name, size: text.length, arrayBuffer: () => pending},
    finish: () => resolve(Uint8Array.from(Buffer.from(text)).buffer),
  };
}
function choose(pending) {
  element("file").files = [pending.file];
  return element("file").handlers.change();
}
(async () => {
  await new Promise(resolve => setImmediate(resolve));
  element("program").value = "iso";
  element("timeout").value = "30";
  element("source").value = "OLD TEXT";
  const first = deferredFile("first.in", "FIRST");
  const firstChange = choose(first);
  assert.equal(element("run").disabled, true);
  await element("run").handlers.click();
  assert.equal(submissions.length, 0, "Pending file must not submit the old textarea.");
  const second = deferredFile("second.in", "SECOND");
  const secondChange = choose(second);
  second.finish();
  await secondChange;
  assert.equal(element("source").value, "SECOND");
  assert.equal(element("run").disabled, false);
  first.finish();
  await firstChange;
  assert.equal(element("source").value, "SECOND", "An older file cannot replace a newer selection.");

  const third = deferredFile("third.in", "THIRD");
  const thirdChange = choose(third);
  element("source").value = "EDITED";
  element("source").handlers.input();
  third.finish();
  await thirdChange;
  assert.equal(element("source").value, "EDITED");
  await element("run").handlers.click();
  assert.equal(submissions.at(-1).input, "EDITED");
  assert.equal("input_b64" in submissions.at(-1), false);

  const fourth = deferredFile("fourth.in", "FOURTH");
  const fourthChange = choose(fourth);
  element("parent").value = "221";
  element("irrep").value = "GM4-";
  element("opd").value = "P1";
  await element("prepare").handlers.click();
  fourth.finish();
  await fourthChange;
  assert.equal(element("source").value, "PREPARED\nQUIT\n");
  assert.equal(element("run").disabled, false);
  process.stdout.write("File-selection ordering, pending-run guard, edit and query cancellation: passed.\n");
})().catch(error => {console.error(error); process.exitCode = 1;});
