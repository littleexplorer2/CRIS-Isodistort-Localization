#!/usr/bin/env node

/**
 * Drive the real local ISODISTORT web page through Microsoft Edge CDP and
 * persist every browser download required by the frozen refresh manifest.
 *
 * No third-party Node package is required. The script launches the existing
 * Python web entry through run_cris.ps1, opens that page in headless Edge,
 * calls the page's own API/download functions, and waits for the resulting
 * browser files before marking a case complete.
 */

import { spawn, spawnSync } from "node:child_process";
import crypto from "node:crypto";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import process from "node:process";
import { fileURLToPath } from "node:url";

const SCRIPT_PATH = fileURLToPath(import.meta.url);
const PACKAGE_ROOT = path.resolve(path.dirname(SCRIPT_PATH), "..", "..");
const REPO_ROOT = path.resolve(PACKAGE_ROOT, "..");
const MANIFEST_PATH = path.join(PACKAGE_ROOT, "docs", "manifests", "local_web_output_refresh_manifest.json");
const METHOD3_MANIFEST_PATH = path.join(PACKAGE_ROOT, "docs", "manifests", "method3_download_manifest.json");
const METHOD4_MANIFEST_PATH = path.join(PACKAGE_ROOT, "docs", "manifests", "method4_download_manifest.json");
const RUNTIME_REVISION_PATH = path.join(
  PACKAGE_ROOT,
  "docs",
  "manifests",
  "local_web_output_refresh_runtime_revision_20261015.json",
);
const EDGE_PATHS = [
  "C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe",
  "C:\\Program Files\\Microsoft\\Edge\\Application\\msedge.exe",
];
const WEB_CONTRACT_PATHS = [
  "ISODISTORT/main_web.py",
  "ISODISTORT/frontend/web/server.py",
  "ISODISTORT/frontend/web/index.html",
  "ISODISTORT/backend/messages/__init__.py",
  "ISODISTORT/backend/messages/messages.py",
  "ISODISTORT/application/session/__init__.py",
  "ISODISTORT/application/session/core_api.py",
  "ISODISTORT/application/workflows/__init__.py",
  "ISODISTORT/application/workflows/service.py",
  "ISODISTORT/resources/config/settings.yaml",
];
const SOURCE_TREE_ROOTS = [
  "ISODISTORT/backend",
  "ISODISTORT/application",
  "ISODISTORT/features",
  "ISODISTORT/frontend",
];
const SOURCE_TREE_FILES = [
  "ISODISTORT/main_web.py",
  "ISODISTORT/pyproject.toml",
  "ISODISTORT/requirements.txt",
  "ISODISTORT/requirements-dev.txt",
  "ISODISTORT/resources/config/settings.yaml",
  "ISODISTORT/docs/manifests/local_web_output_refresh_manifest.json",
  "ISODISTORT/docs/manifests/method3_download_manifest.json",
  "ISODISTORT/docs/manifests/method4_download_manifest.json",
  "ISODISTORT/docs/manifests/local_web_output_refresh_runtime_revision_20261015.json",
  "ISODISTORT/tests/manual/download_local_web_outputs.mjs",
  "run_cris.ps1",
];
let activeServerUrl = null;

function nowIso() {
  return new Date().toISOString();
}

function log(message, extra = undefined) {
  const suffix = extra === undefined ? "" : ` ${typeof extra === "string" ? extra : JSON.stringify(extra)}`;
  process.stdout.write(`[${nowIso()}] ${message}${suffix}\n`);
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

function readJson(file) {
  return JSON.parse(fs.readFileSync(file, "utf8"));
}

function writeJson(file, value) {
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, `${JSON.stringify(value, null, 2)}\n`, "utf8");
}

function sha256Buffer(buffer) {
  return crypto.createHash("sha256").update(buffer).digest("hex");
}

function sha256File(file) {
  return sha256Buffer(fs.readFileSync(file));
}

function relativePosix(file) {
  return path.relative(REPO_ROOT, file).split(path.sep).join("/");
}

function listFiles(root) {
  if (!fs.existsSync(root)) return [];
  const stat = fs.statSync(root);
  if (stat.isFile()) return [root];
  const files = [];
  const stack = [root];
  while (stack.length) {
    const current = stack.pop();
    for (const entry of fs.readdirSync(current, { withFileTypes: true })) {
      const child = path.join(current, entry.name);
      if (entry.isDirectory()) stack.push(child);
      else if (entry.isFile()) files.push(child);
    }
  }
  return files.sort((a, b) => relativePosix(a).localeCompare(relativePosix(b)));
}

function hashFileSet(files) {
  const entries = files.map((file) => ({
    path: relativePosix(file),
    bytes: fs.statSync(file).size,
    sha256: sha256File(file),
  }));
  const digest = crypto.createHash("sha256");
  for (const entry of entries) {
    digest.update(`${entry.path}\0${entry.bytes}\0${entry.sha256}\n`, "utf8");
  }
  return { sha256: digest.digest("hex"), file_count: entries.length, entries };
}

function sourceSignature() {
  const files = [];
  for (const rel of SOURCE_TREE_ROOTS) {
    for (const file of listFiles(path.join(REPO_ROOT, rel))) {
      if (file.endsWith(".pyc") || file.includes(`${path.sep}__pycache__${path.sep}`)) continue;
      files.push(file);
    }
  }
  for (const rel of SOURCE_TREE_FILES) {
    const file = path.join(REPO_ROOT, rel);
    if (fs.existsSync(file)) files.push(file);
  }
  const unique = [...new Set(files.map((file) => path.resolve(file)))].sort();
  return hashFileSet(unique);
}

function runText(command, args, options = {}) {
  const result = spawnSync(command, args, {
    cwd: options.cwd || REPO_ROOT,
    encoding: "utf8",
    windowsHide: true,
    maxBuffer: 64 * 1024 * 1024,
  });
  if (result.error) throw result.error;
  if (result.status !== 0 && !options.allowFailure) {
    throw new Error(`${command} exited ${result.status}: ${(result.stderr || result.stdout || "").trim()}`);
  }
  return {
    status: result.status,
    stdout: (result.stdout || "").trim(),
    stderr: (result.stderr || "").trim(),
  };
}

function pythonEnvironmentSignature() {
  const code = [
    "import hashlib, importlib.metadata, json, platform, sys",
    "d=[]",
    "for x in importlib.metadata.distributions():",
    " n=(x.metadata.get('Name') or '').strip()",
    " v=(x.version or '').strip()",
    " if n: d.append([n.lower(),v])",
    "d=sorted(set(map(tuple,d)))",
    "raw='\\n'.join(f'{n}=={v}' for n,v in d)",
    "print(json.dumps({'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'distributions':d,'distributions_sha256':hashlib.sha256(raw.encode()).hexdigest()}))",
  ].join("\n");
  const result = runText("powershell.exe", [
    "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
    path.join(REPO_ROOT, "run_cris.ps1"), "-c", code,
  ]);
  return JSON.parse(result.stdout.split(/\r?\n/).filter(Boolean).at(-1));
}

function gitSignature() {
  const head = runText("git", ["rev-parse", "HEAD"], { allowFailure: true });
  const status = runText("git", ["status", "--porcelain=v1", "--untracked-files=all"], { allowFailure: true });
  return {
    head: head.status === 0 ? head.stdout : null,
    status_sha256: sha256Buffer(Buffer.from(status.stdout, "utf8")),
    status_lines: status.stdout ? status.stdout.split(/\r?\n/).length : 0,
  };
}

function binarySignature() {
  const result = {};
  for (const name of ["iso", "findsym", "smodes"]) {
    const file = path.join(PACKAGE_ROOT, "resources", "isobyu", name);
    result[name] = fs.existsSync(file)
      ? { path: relativePosix(file), bytes: fs.statSync(file).size, sha256: sha256File(file) }
      : { path: relativePosix(file), missing: true };
  }
  return result;
}

function captureRuntimeSignature(edgePath, includeDependencies = false) {
  const webContracts = {};
  for (const rel of WEB_CONTRACT_PATHS) {
    const file = path.join(REPO_ROOT, rel);
    webContracts[rel] = sha256File(file);
  }
  const edgeVersion = runText(edgePath, ["--version"], { allowFailure: true });
  return {
    captured_at: nowIso(),
    source_tree: sourceSignature(),
    web_contracts: webContracts,
    git: gitSignature(),
    binaries: binarySignature(),
    edge: { path: edgePath, version: edgeVersion.stdout || edgeVersion.stderr || null },
    node: process.version,
    python_environment: includeDependencies ? pythonEnvironmentSignature() : undefined,
  };
}

function assertFrozenContracts(manifest, signature, runtimeRevision) {
  const expected = runtimeRevision?.web_contract_files
    || manifest.source_snapshot.web_contract_files;
  const mismatches = [];
  for (const [rel, hash] of Object.entries(expected)) {
    const actual = signature.web_contracts[rel];
    if (String(actual).toLowerCase() !== String(hash).toLowerCase()) {
      mismatches.push({ path: rel, expected: hash, actual });
    }
  }
  if (mismatches.length) {
    throw new Error(`Frozen web contracts drifted: ${JSON.stringify(mismatches)}`);
  }
  if (runtimeRevision) {
    const parentHash = sha256File(MANIFEST_PATH);
    if (runtimeRevision.parent_manifest_sha256 !== parentHash) {
      throw new Error(
        `Runtime revision points to a different frozen manifest: ${runtimeRevision.parent_manifest_sha256} != ${parentHash}`,
      );
    }
  }
}

function findEdge() {
  const edge = EDGE_PATHS.find((candidate) => fs.existsSync(candidate));
  if (!edge) throw new Error("Microsoft Edge was not found in the supported installation paths");
  return edge;
}

function parseArgs(argv) {
  const args = { cases: [], stopAfter: null };
  for (let i = 0; i < argv.length; i += 1) {
    const token = argv[i];
    if (token === "--case") args.cases.push(argv[++i]);
    else if (token === "--stop-after") args.stopAfter = Number(argv[++i]);
    else if (token === "--help") {
      process.stdout.write(
        "Usage: node download_local_web_outputs.mjs [--case CASE_ID] [--stop-after N]\n",
      );
      process.exit(0);
    } else throw new Error(`Unknown argument: ${token}`);
  }
  if (args.stopAfter !== null && (!Number.isInteger(args.stopAfter) || args.stopAfter < 1)) {
    throw new Error("--stop-after must be a positive integer");
  }
  return args;
}

function flattenCases(manifest, method3Manifest, method4Manifest) {
  const cases = [];
  for (const item of manifest.method1_cases) {
    cases.push({ ...item, method: 1, parent: manifest.parents[item.parent_ref] });
  }
  for (const item of manifest.method2_cases) {
    cases.push({ ...item, method: 2, parent: manifest.parents[item.parent_ref] });
  }
  const method3ById = new Map();
  for (const parent of method3Manifest.parents) {
    for (const item of parent.cases) method3ById.set(item.id, { ...item, parentName: parent.parent_cif });
  }
  for (const id of manifest.method3_cases.case_ids) {
    const source = method3ById.get(id);
    if (!source) throw new Error(`Method 3 source case missing: ${id}`);
    const parentRef = source.parentName === manifest.parents.EU.name ? "EU" : "ND";
    cases.push({
      id,
      method: 3,
      parent_ref: parentRef,
      parent: manifest.parents[parentRef],
      query: source,
      output_dir: manifest.method3_cases.output_dir_template
        .replace("{parent_cif}", source.parentName)
        .replace("{case_id}", id),
      artifact_contract: manifest.method3_cases.artifact_contract,
    });
  }
  for (const item of manifest.method4_cases.cases) {
    const parentSource = method4Manifest.parents[item.source_parent];
    if (!parentSource) throw new Error(`Method 4 parent missing: ${item.source_parent}`);
    const source = parentSource.cases[item.source_case_id];
    if (!source) throw new Error(`Method 4 source case missing: ${item.source_parent}/${item.source_case_id}`);
    const parentRef = item.source_parent === manifest.parents.EU.name ? "EU" : "ND";
    cases.push({
      ...item,
      method: 4,
      parent_ref: parentRef,
      parent: manifest.parents[parentRef],
      query: source,
      output_dir: manifest.method4_cases.output_dir_template
        .replace("{parent_cif}", item.source_parent)
        .replace("{id}", item.id),
    });
  }
  return cases;
}

class CdpClient {
  constructor(url) {
    this.url = url;
    this.ws = null;
    this.nextId = 1;
    this.pending = new Map();
  }

  async connect() {
    this.ws = new WebSocket(this.url);
    await new Promise((resolve, reject) => {
      const timer = setTimeout(() => reject(new Error("Timed out opening CDP websocket")), 15000);
      this.ws.addEventListener("open", () => { clearTimeout(timer); resolve(); }, { once: true });
      this.ws.addEventListener("error", (event) => { clearTimeout(timer); reject(new Error(`CDP websocket error: ${event.message || "unknown"}`)); }, { once: true });
    });
    this.ws.addEventListener("message", (event) => {
      const message = JSON.parse(String(event.data));
      if (!message.id) return;
      const pending = this.pending.get(message.id);
      if (!pending) return;
      this.pending.delete(message.id);
      if (message.error) pending.reject(new Error(`${pending.method}: ${message.error.message}`));
      else pending.resolve(message.result || {});
    });
    this.ws.addEventListener("close", () => {
      for (const pending of this.pending.values()) pending.reject(new Error("CDP websocket closed"));
      this.pending.clear();
    });
  }

  send(method, params = {}) {
    const id = this.nextId++;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject, method });
      this.ws.send(JSON.stringify({ id, method, params }));
    });
  }

  async evaluate(expression) {
    const result = await this.send("Runtime.evaluate", {
      expression,
      awaitPromise: true,
      returnByValue: true,
      userGesture: true,
    });
    if (result.exceptionDetails) {
      const detail = result.exceptionDetails.exception?.description
        || result.exceptionDetails.text
        || "JavaScript evaluation failed";
      throw new Error(detail);
    }
    return result.result?.value;
  }

  close() {
    if (this.ws && this.ws.readyState <= WebSocket.OPEN) this.ws.close();
  }
}

async function waitForFile(file, timeoutMs = 300000) {
  const started = Date.now();
  let previousSize = -1;
  let stableCount = 0;
  while (Date.now() - started < timeoutMs) {
    if (fs.existsSync(file) && !fs.existsSync(`${file}.crdownload`)) {
      const size = fs.statSync(file).size;
      stableCount = size > 0 && size === previousSize ? stableCount + 1 : 0;
      previousSize = size;
      if (stableCount >= 2) return file;
    }
    await sleep(250);
  }
  throw new Error(`Timed out waiting for browser download: ${file}`);
}

async function waitFor(predicate, description, timeoutMs = 30000) {
  const started = Date.now();
  while (Date.now() - started < timeoutMs) {
    const value = await predicate();
    if (value) return value;
    await sleep(250);
  }
  throw new Error(`Timed out waiting for ${description}`);
}

async function startServer(batchRuntimeDir) {
  fs.mkdirSync(batchRuntimeDir, { recursive: true });
  const stdoutPath = path.join(batchRuntimeDir, "server.stdout.log");
  const stderrPath = path.join(batchRuntimeDir, "server.stderr.log");
  const stdoutLog = fs.createWriteStream(stdoutPath, { flags: "a" });
  const stderrLog = fs.createWriteStream(stderrPath, { flags: "a" });
  const code = "from frontend.web import server as s; s._open_browser=lambda url: False; raise SystemExit(s.main())";
  const child = spawn("powershell.exe", [
    "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
    path.join(REPO_ROOT, "run_cris.ps1"), "-c", code,
  ], { cwd: REPO_ROOT, windowsHide: true, stdio: ["ignore", "pipe", "pipe"] });
  let combined = "";
  let url = null;
  const onChunk = (chunk, stream) => {
    const text = chunk.toString("utf8");
    stream.write(text);
    combined += text;
    const match = combined.match(/URL:\s*(http:\/\/127\.0\.0\.1:\d+\/)/);
    if (match) url = match[1];
  };
  child.stdout.on("data", (chunk) => onChunk(chunk, stdoutLog));
  child.stderr.on("data", (chunk) => onChunk(chunk, stderrLog));
  child.on("exit", (codeValue) => log(`Local web server exited (${codeValue})`));
  await waitFor(() => url, "local web server URL", 60000);
  log("Local web server started", url);
  return { child, url, stdoutLog, stderrLog };
}

async function startEdge(edgePath, serverUrl, batchRuntimeDir) {
  const profileDir = fs.mkdtempSync(path.join(os.tmpdir(), "cris-edge-refresh-"));
  const stdoutPath = path.join(batchRuntimeDir, "edge.stdout.log");
  const stderrPath = path.join(batchRuntimeDir, "edge.stderr.log");
  const stdoutLog = fs.createWriteStream(stdoutPath, { flags: "a" });
  const stderrLog = fs.createWriteStream(stderrPath, { flags: "a" });
  const child = spawn(edgePath, [
    "--headless=new",
    "--remote-debugging-port=0",
    "--remote-allow-origins=*",
    `--user-data-dir=${profileDir}`,
    "--no-first-run",
    "--no-default-browser-check",
    "--disable-extensions",
    "--disable-background-networking",
    "--disable-background-timer-throttling",
    "--disable-renderer-backgrounding",
    serverUrl,
  ], { cwd: REPO_ROOT, windowsHide: true, stdio: ["ignore", "pipe", "pipe"] });
  child.stdout.pipe(stdoutLog);
  child.stderr.pipe(stderrLog);
  const portFile = path.join(profileDir, "DevToolsActivePort");
  await waitFor(() => fs.existsSync(portFile), "Edge DevToolsActivePort", 30000);
  const [portText] = fs.readFileSync(portFile, "utf8").split(/\r?\n/);
  const debugBase = `http://127.0.0.1:${Number(portText)}`;
  const target = await waitFor(async () => {
    try {
      const response = await fetch(`${debugBase}/json/list`);
      const targets = await response.json();
      return targets.find((item) => item.type === "page" && item.url.startsWith(serverUrl))
        || targets.find((item) => item.type === "page");
    } catch {
      return null;
    }
  }, "Edge page target", 30000);
  const cdp = new CdpClient(target.webSocketDebuggerUrl);
  await cdp.connect();
  await cdp.send("Page.enable");
  await cdp.send("Runtime.enable");
  await waitFor(async () => {
    try {
      return await cdp.evaluate("document.readyState === 'complete' && typeof api === 'function'");
    } catch {
      return false;
    }
  }, "local web page initialization", 30000);
  log("Headless Edge connected to the real local page", { debug_port: Number(portText), target: target.url });
  return { child, cdp, profileDir, stdoutLog, stderrLog };
}

async function setDownloadDirectory(cdp, directory) {
  fs.mkdirSync(directory, { recursive: true });
  try {
    await cdp.send("Browser.setDownloadBehavior", {
      behavior: "allow",
      downloadPath: directory,
      eventsEnabled: true,
    });
  } catch {
    await cdp.send("Page.setDownloadBehavior", { behavior: "allow", downloadPath: directory });
  }
}

async function downloadFromPage(cdp, caseDir, method, kind, finalName, expression, timeoutMs) {
  const attemptDir = path.join(caseDir, "_browser_downloads", `${Date.now()}-${kind}`);
  await setDownloadDirectory(cdp, attemptDir);
  const originalName = kind === "zip"
    ? `isodistort_method${method}.zip`
    : `isodistort_method${method}_filtered.${kind}`;
  log(`Triggering browser ${kind.toUpperCase()} download`, { method, final: path.join(caseDir, finalName) });
  let lastProgress = null;
  const progressTimer = kind === "zip" && activeServerUrl
    ? setInterval(async () => {
        try {
          const response = await fetch(`${activeServerUrl}api/export_status`);
          const progress = (await response.json()).export || {};
          const marker = `${progress.completed}/${progress.total}/${progress.phase}`;
          if (marker !== lastProgress) {
            lastProgress = marker;
            log("Web ZIP progress", {
              method,
              phase: progress.phase,
              completed: progress.completed,
              total: progress.total,
              successful: progress.successful,
              ineligible: progress.ineligible,
              failed: progress.failed,
              candidate_index: Number.isInteger(progress.candidate_index)
                ? progress.candidate_index : null,
              candidate_status: progress.candidate_status || null,
              candidate: progress.irrep_label
                ? `${progress.irrep_label}/${progress.opd_symbol || "-"}` : null,
              candidate_error_type: progress.candidate_error_type || null,
              candidate_error_message: progress.candidate_error_message || null,
              estimated_remaining_seconds: Number.isFinite(progress.estimated_remaining_seconds)
                ? Math.round(progress.estimated_remaining_seconds) : null,
            });
          }
        } catch { /* the browser download remains authoritative */ }
      }, 30000)
    : null;
  try {
    await cdp.evaluate(expression);
  } finally {
    if (progressTimer) clearInterval(progressTimer);
  }
  const downloaded = await waitForFile(path.join(attemptDir, originalName), timeoutMs);
  const destination = path.join(caseDir, finalName);
  if (fs.existsSync(destination)) {
    const oldHash = sha256File(destination);
    const newHash = sha256File(downloaded);
    if (oldHash !== newHash) {
      throw new Error(`Existing artifact conflicts with new browser download: ${destination}`);
    }
    return destination;
  }
  fs.renameSync(downloaded, destination);
  log("Browser download persisted", { file: relativePosix(destination), bytes: fs.statSync(destination).size, sha256: sha256File(destination) });
  return destination;
}

function artifactInfo(file) {
  return {
    relative_path: path.relative(path.dirname(file), file).split(path.sep).join("/"),
    repository_path: relativePosix(file),
    byte_size: fs.statSync(file).size,
    sha256: sha256File(file),
  };
}

function archiveExistingAttemptFiles(caseDir) {
  const names = [
    "complete_unfiltered_result_table.txt",
    "complete_unfiltered_result_table.csv",
    "raw_all_candidates_four_format.zip",
    "complete_result_table.txt",
    "complete_result_table.csv",
    "case_execution.json",
    "raw_web_or_api_error.json",
    "raw_web_response.json",
    "zip_inventory.json",
  ];
  const existing = names.map((name) => path.join(caseDir, name)).filter((file) => fs.existsSync(file));
  for (const name of ["raw_all_candidates_four_format", "_browser_downloads"]) {
    const candidate = path.join(caseDir, name);
    if (fs.existsSync(candidate)) existing.push(candidate);
  }
  if (!existing.length) return;
  const archiveDir = path.join(caseDir, "attempts", nowIso().replace(/[:.]/g, "-"));
  fs.mkdirSync(archiveDir, { recursive: true });
  for (const file of existing) fs.renameSync(file, path.join(archiveDir, path.basename(file)));
  log("Preserved incomplete prior attempt", relativePosix(archiveDir));
}

function isCaseComplete(caseDir, method, currentSignature) {
  const recordPath = path.join(caseDir, "case_execution.json");
  if (!fs.existsSync(recordPath)) return false;
  try {
    const record = readJson(recordPath);
    if (record.status !== "complete") return false;
    if (
      record.end_runtime_signature?.source_tree?.sha256
      !== currentSignature.source_tree.sha256
    ) return false;
    for (const [rel, hash] of Object.entries(currentSignature.web_contracts)) {
      if (record.end_runtime_signature?.web_contracts?.[rel] !== hash) return false;
    }
    const required = method <= 3
      ? ["complete_unfiltered_result_table.txt", "complete_unfiltered_result_table.csv", "raw_all_candidates_four_format.zip", "zip_inventory.json"]
      : record.result_kind === "success"
        ? ["complete_result_table.txt", "complete_result_table.csv", "raw_web_response.json"]
        : ["raw_web_or_api_error.json"];
    if (!required.every((name) => {
      const file = path.join(caseDir, name);
      return fs.existsSync(file) && fs.statSync(file).size > 0;
    })) return false;
    for (const artifact of record.artifacts || []) {
      const file = path.join(REPO_ROOT, artifact.repository_path);
      if (!fs.existsSync(file)) return false;
      if (fs.statSync(file).size !== artifact.byte_size) return false;
      if (sha256File(file) !== artifact.sha256) return false;
    }
    if (method <= 3) {
      const inventory = readJson(path.join(caseDir, "zip_inventory.json"));
      if (inventory.status_report?.status !== "success") return false;
      if (Number(inventory.status_report?.ineligible_candidate_count) !== 0) return false;
      if (Number(inventory.candidate_count) !== Number(inventory.status_report?.requested_candidate_count)) return false;
      if ((inventory.incomplete_candidates || []).length) return false;
      const extractedRoot = path.join(REPO_ROOT, inventory.extracted_root);
      for (const entry of inventory.entries || []) {
        const file = path.join(extractedRoot, ...entry.relative_path.split("/"));
        if (!fs.existsSync(file) || fs.statSync(file).size !== entry.byte_size) return false;
        if (sha256File(file) !== entry.sha256) return false;
      }
    }
    return true;
  } catch {
    return false;
  }
}

function loadParentExpression(caseSpec) {
  const parentPath = path.join(REPO_ROOT, caseSpec.parent.path);
  const content = fs.readFileSync(parentPath, "utf8");
  const payload = { filename: path.basename(parentPath), content };
  const types = caseSpec.parent.distortion_types;
  const scope = caseSpec.parent.distortion_scope;
  return `(async()=>{
    const loaded=await api('/api/load_cif',${JSON.stringify(payload)});
    if(!loaded.ok) throw new Error(loaded.error||'load_cif failed');
    const typed=await api('/api/set_types',${JSON.stringify({ types, scope })});
    if(!typed.ok) throw new Error(typed.error||'set_types failed');
    return {revision:STATE&&STATE.revision,structure:STATE&&STATE.structure};
  })()`;
}

function setupMethod1Expression(caseSpec) {
  const query = {
    distortion_types: caseSpec.parent.distortion_types,
    crystal_system: caseSpec.query.crystal_system,
    subgroup_space_group: caseSpec.query.subgroup_space_group,
    maximal_subgroup_only: Boolean(caseSpec.query.maximal_subgroup_only),
  };
  if (caseSpec.query.lattice) query.lattice = caseSpec.query.lattice;
  return `(async()=>{
    const d=await api('/api/method1',${JSON.stringify(query)});
    if(!d.ok) throw new Error(d.error||'method1 failed');
    m1Items=d.candidates||[]; tblState[1]=emptyTblState();
    document.getElementById('dlMethod').value='1';
    for(const id of ['fmtCif','fmtIsoviz','fmtModes','fmtTopas']) document.getElementById(id).checked=true;
    return {candidate_count:m1Items.length,revision:STATE&&STATE.revision};
  })()`;
}

function setupMethod2Expression(caseSpec) {
  const parameters = caseSpec.query.k_parameters || {};
  const params = Object.keys(parameters).length ? Object.values(parameters) : null;
  const payload = {
    nsup: caseSpec.query.number_of_superposed_irs || 1,
    kpoints: [{ k: caseSpec.query.k_label, params }],
    generate: Boolean(caseSpec.query.generate_if_missing),
  };
  return `(async()=>{
    const d=await api('/api/subgroups',${JSON.stringify(payload)});
    if(!d.ok) throw new Error(d.error||'method2 subgroup enumeration failed');
    subgroups=d.subgroups||[]; tblState[2]=emptyTblState();
    nsupValue=${Number(caseSpec.query.number_of_superposed_irs || 1)};
    nmodValue=${Number(caseSpec.query.nmod || 0)};
    document.getElementById('dlMethod').value='2';
    for(const id of ['fmtCif','fmtIsoviz','fmtModes','fmtTopas']) document.getElementById(id).checked=true;
    return {candidate_count:subgroups.length,revision:STATE&&STATE.revision};
  })()`;
}

function setupMethod3Expression(caseSpec) {
  const query = caseSpec.query;
  const payload = {
    distortion_types: query.distortion_types,
    point_group: null,
    space_group_type: query.space_group_type,
    supercell_basis: query.supercell_basis,
    direct_sublattice_centering: query.direct_sublattice_centering,
    lattice_type: query.lattice_type,
    generate: false,
  };
  return `(async()=>{
    const d=await api('/api/method3',${JSON.stringify(payload)});
    if(!d.ok) throw new Error(d.error||'method3 failed');
    m3Items=d.candidates||[]; tblState[3]=emptyTblState();
    document.getElementById('dlMethod').value='3';
    for(const id of ['fmtCif','fmtIsoviz','fmtModes','fmtTopas']) document.getElementById(id).checked=true;
    return {candidate_count:m3Items.length,revision:STATE&&STATE.revision};
  })()`;
}

function setupMethod4Expression(caseSpec) {
  const daughterPath = path.join(REPO_ROOT, caseSpec.query.daughter_cif);
  const payload = {
    filename: path.basename(daughterPath),
    content: fs.readFileSync(daughterPath, "utf8"),
    atom_matching_method: caseSpec.query.atom_matching_method || "nearest-site",
    robust_distance_threshold: Number(caseSpec.query.robust_distance_threshold ?? 0.25),
    provided_origin_shift: caseSpec.query.provided_origin_shift || null,
  };
  return `(async()=>{
    const d=await api('/api/method4',${JSON.stringify(payload)});
    if(d.ok){
      m4Items=Object.entries(d.amplitudes||{}).map(([mode,As])=>({
        mode,As:+As,Ap:+(d.parent_cell_amplitudes||{})[mode],
        raw:+(d.raw_coefficients||{})[mode],normfactor:+(d.mode_normfactors||{})[mode]
      }));
      m4Meta={rms:d.rms_residual,maxAbs:d.max_abs_residual,
        strainApplied:d.strain_applied_engineering_q_parent_basis||null,
        strainRaw:d.strain_raw_coordinate_sum_parent_basis||null,
        strainTensor:d.strain_tensor_parent_basis||null,
        strainMultiplier:d.strain_multiplier_parent_basis||null,
        strainModes:d.strain_modes||[],metadata:d.metadata||null};
      tblState[4]=emptyTblState(); document.getElementById('dlMethod').value='4';
    }
    return d;
  })()`;
}

async function inventoryZip(zipPath, caseDir) {
  const extractDir = path.join(caseDir, "raw_all_candidates_four_format");
  if (fs.existsSync(extractDir)) {
    throw new Error(`ZIP extraction directory already exists: ${extractDir}`);
  }
  fs.mkdirSync(extractDir, { recursive: true });
  const result = runText("tar.exe", ["-xf", zipPath, "-C", extractDir], { allowFailure: true });
  if (result.status !== 0) throw new Error(`tar.exe could not extract ${zipPath}: ${result.stderr || result.stdout}`);
  const entries = listFiles(extractDir).map((file) => ({
    relative_path: path.relative(extractDir, file).split(path.sep).join("/"),
    byte_size: fs.statSync(file).size,
    sha256: sha256File(file),
    candidate_identity: path.relative(extractDir, file).split(path.sep)[0] || null,
    status: "present",
  }));
  const basenames = new Set(entries.map((entry) => path.posix.basename(entry.relative_path)));
  for (const required of ["export_candidate_status.json", "export_candidate_status.txt"]) {
    if (!basenames.has(required)) throw new Error(`ZIP is missing required root report: ${required}`);
  }
  const reportEntry = entries.find(
    (entry) => path.posix.basename(entry.relative_path) === "export_candidate_status.json",
  );
  const report = readJson(path.join(extractDir, ...reportEntry.relative_path.split("/")));
  const reportPrefix = path.posix.dirname(reportEntry.relative_path) === "."
    ? "" : `${path.posix.dirname(reportEntry.relative_path)}/`;
  const candidateFiles = new Map();
  for (const entry of entries) {
    if (!entry.relative_path.startsWith(reportPrefix)) continue;
    const remainder = entry.relative_path.slice(reportPrefix.length);
    const parts = remainder.split("/");
    if (parts.length !== 2) continue;
    const [folder, filename] = parts;
    if (!candidateFiles.has(folder)) candidateFiles.set(folder, new Set());
    candidateFiles.get(folder).add(filename);
  }
  const requiredCandidateFiles = new Set([
    "subgroup.cif", "data.isoviz", "Complete modes details.txt", "topas.str",
  ]);
  const incompleteCandidates = [];
  for (const [folder, files] of candidateFiles) {
    const missing = [...requiredCandidateFiles].filter((name) => !files.has(name));
    if (missing.length) incompleteCandidates.push({ folder, missing });
  }
  if (incompleteCandidates.length) {
    throw new Error(`ZIP has incomplete candidate folders: ${JSON.stringify(incompleteCandidates)}`);
  }
  if (candidateFiles.size !== Number(report.successful_candidate_count)) {
    throw new Error(
      `ZIP candidate count ${candidateFiles.size} != status report successful count ${report.successful_candidate_count}`,
    );
  }
  if (
    report.status !== "success"
    || Number(report.ineligible_candidate_count) !== 0
    || Number(report.successful_candidate_count) !== Number(report.requested_candidate_count)
  ) {
    throw new Error(
      `ZIP is not a complete successful batch: status=${report.status} requested=${report.requested_candidate_count} successful=${report.successful_candidate_count} ineligible=${report.ineligible_candidate_count}`,
    );
  }
  const inventory = {
    schema: "local-web-output-zip-inventory/v1",
    generated_at: nowIso(),
    raw_zip: { byte_size: fs.statSync(zipPath).size, sha256: sha256File(zipPath) },
    extracted_root: relativePosix(extractDir),
    entry_count: entries.length,
    candidate_count: candidateFiles.size,
    status_report: report,
    incomplete_candidates: incompleteCandidates,
    entries,
  };
  writeJson(path.join(caseDir, "zip_inventory.json"), inventory);
  return inventory;
}

async function executeCase(cdp, caseSpec, outputRoot, startSignature) {
  const caseDir = path.join(outputRoot, ...caseSpec.output_dir.split("/"));
  fs.mkdirSync(caseDir, { recursive: true });
  if (isCaseComplete(caseDir, caseSpec.method, startSignature)) {
    log(`Skipping already complete case ${caseSpec.id}`, relativePosix(caseDir));
    return { id: caseSpec.id, status: "complete", skipped: true };
  }
  archiveExistingAttemptFiles(caseDir);
  const startedAt = nowIso();
  const record = {
    schema: "local-web-output-case-execution/v1",
    case_id: caseSpec.id,
    method: caseSpec.method,
    status: "running",
    started_at: startedAt,
    parent: {
      path: caseSpec.parent.path,
      expected_sha256: caseSpec.parent.sha256,
      actual_sha256: sha256File(path.join(REPO_ROOT, caseSpec.parent.path)),
    },
    query: caseSpec.query,
    start_runtime_signature: startSignature,
    artifacts: [],
  };
  if (caseSpec.method === 4) {
    const daughterPath = path.join(REPO_ROOT, caseSpec.query.daughter_cif);
    record.daughter = {
      path: caseSpec.query.daughter_cif,
      expected_sha256: caseSpec.query.sha256,
      actual_sha256: sha256File(daughterPath),
    };
  }
  writeJson(path.join(caseDir, "case_execution.json"), record);
  log(`Starting ${caseSpec.id}`, { method: caseSpec.method, output: relativePosix(caseDir) });

  try {
    if (record.parent.actual_sha256 !== record.parent.expected_sha256) {
      throw new Error(`Parent hash mismatch for ${caseSpec.id}`);
    }
    if (record.daughter && record.daughter.actual_sha256 !== record.daughter.expected_sha256) {
      throw new Error(`Daughter hash mismatch for ${caseSpec.id}`);
    }
    await cdp.evaluate(loadParentExpression(caseSpec));
    let responseSummary;
    if (caseSpec.method === 1) responseSummary = await cdp.evaluate(setupMethod1Expression(caseSpec));
    else if (caseSpec.method === 2) responseSummary = await cdp.evaluate(setupMethod2Expression(caseSpec));
    else if (caseSpec.method === 3) responseSummary = await cdp.evaluate(setupMethod3Expression(caseSpec));
    else responseSummary = await cdp.evaluate(setupMethod4Expression(caseSpec));

    if (caseSpec.method <= 3) {
      if (!responseSummary || responseSummary.candidate_count < 1) {
        throw new Error(`Method ${caseSpec.method} returned no candidates`);
      }
      record.web_response_summary = responseSummary;
      const txt = await downloadFromPage(
        cdp, caseDir, caseSpec.method, "txt", "complete_unfiltered_result_table.txt",
        "(()=>{tblState[Number(document.getElementById('dlMethod').value)]=emptyTblState();downloadFilteredTable(false);return true})()",
        60000,
      );
      const csv = await downloadFromPage(
        cdp, caseDir, caseSpec.method, "csv", "complete_unfiltered_result_table.csv",
        "(()=>{tblState[Number(document.getElementById('dlMethod').value)]=emptyTblState();downloadFilteredTable(true);return true})()",
        60000,
      );
      const zip = await downloadFromPage(
        cdp, caseDir, caseSpec.method, "zip", "raw_all_candidates_four_format.zip",
        "downloadAll()",
        6 * 60 * 60 * 1000,
      );
      const inventory = await inventoryZip(zip, caseDir);
      record.artifacts = [artifactInfo(txt), artifactInfo(csv), artifactInfo(zip), artifactInfo(path.join(caseDir, "zip_inventory.json"))];
      record.zip_entry_count = inventory.entry_count;
      record.result_kind = "success";
    } else {
      const expectsFailure = String(caseSpec.query.expect || "").startsWith("reject_");
      if (!responseSummary || !responseSummary.ok) {
        const errorFile = path.join(caseDir, "raw_web_or_api_error.json");
        writeJson(errorFile, {
          captured_at: nowIso(),
          response: responseSummary,
          submitted_parameters: {
            atom_matching_method: caseSpec.query.atom_matching_method || "nearest-site",
            robust_distance_threshold: Number(caseSpec.query.robust_distance_threshold ?? 0.25),
            provided_origin_shift: caseSpec.query.provided_origin_shift || null,
          },
          parent_sha256: record.parent.actual_sha256,
          daughter_sha256: record.daughter.actual_sha256,
        });
        record.artifacts = [artifactInfo(errorFile)];
        record.web_response_summary = responseSummary;
        record.result_kind = "failure_response";
        if (!expectsFailure) {
          throw new Error(
            `Method 4 returned an unexpected failure: ${responseSummary?.error || "missing response"}`,
          );
        }
      } else {
        const responseFile = path.join(caseDir, "raw_web_response.json");
        writeJson(responseFile, responseSummary);
        const txt = await downloadFromPage(
          cdp, caseDir, 4, "txt", "complete_result_table.txt",
          "(()=>{tblState[4]=emptyTblState();document.getElementById('dlMethod').value='4';downloadFilteredTable(false);return true})()",
          60000,
        );
        const csv = await downloadFromPage(
          cdp, caseDir, 4, "csv", "complete_result_table.csv",
          "(()=>{tblState[4]=emptyTblState();document.getElementById('dlMethod').value='4';downloadFilteredTable(true);return true})()",
          60000,
        );
        record.artifacts = [artifactInfo(txt), artifactInfo(csv), artifactInfo(responseFile)];
        record.web_response_summary = {
          ok: true,
          amplitude_count: Object.keys(responseSummary.amplitudes || {}).length,
          rms_residual: responseSummary.rms_residual,
          max_abs_residual: responseSummary.max_abs_residual,
        };
        record.result_kind = "success";
        if (expectsFailure) {
          throw new Error("Method 4 unexpectedly accepted a case whose frozen expectation is rejection");
        }
      }
    }
    const endSignature = captureRuntimeSignature(findEdge(), false);
    record.end_runtime_signature = endSignature;
    if (endSignature.source_tree.sha256 !== startSignature.source_tree.sha256) {
      throw new Error("Source tree drifted while the case was running; artifacts are not valid for the frozen runtime signature");
    }
    record.status = "complete";
    record.finished_at = nowIso();
    writeJson(path.join(caseDir, "case_execution.json"), record);
    log(`Completed ${caseSpec.id}`, { result_kind: record.result_kind, artifacts: record.artifacts.length });
    return { id: caseSpec.id, status: "complete", result_kind: record.result_kind, candidate_count: responseSummary?.candidate_count };
  } catch (error) {
    record.status = "blocked";
    record.finished_at = nowIso();
    record.error = { message: error.message, stack: error.stack };
    record.end_runtime_signature = captureRuntimeSignature(findEdge(), false);
    writeJson(path.join(caseDir, "case_execution.json"), record);
    log(`Blocked ${caseSpec.id}`, error.message);
    return { id: caseSpec.id, status: "blocked", error: error.message };
  }
}

async function shutdown(server, edge) {
  try {
    if (server?.url) await fetch(`${server.url}api/shutdown`, { method: "POST" });
  } catch { /* best effort */ }
  try { edge?.cdp?.close(); } catch { /* best effort */ }
  if (edge?.child && !edge.child.killed) edge.child.kill();
  if (server?.child && !server.child.killed) {
    await Promise.race([
      new Promise((resolve) => server.child.once("exit", resolve)),
      sleep(5000),
    ]);
    if (server.child.exitCode === null) server.child.kill();
  }
  for (const stream of [server?.stdoutLog, server?.stderrLog, edge?.stdoutLog, edge?.stderrLog]) {
    try { stream?.end(); } catch { /* best effort */ }
  }
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const manifest = readJson(MANIFEST_PATH);
  const method3Manifest = readJson(METHOD3_MANIFEST_PATH);
  const method4Manifest = readJson(METHOD4_MANIFEST_PATH);
  const runtimeRevision = fs.existsSync(RUNTIME_REVISION_PATH)
    ? readJson(RUNTIME_REVISION_PATH) : null;
  if (manifest.status !== "frozen") throw new Error(`Manifest is not frozen: ${manifest.status}`);
  const allCases = flattenCases(manifest, method3Manifest, method4Manifest);
  if (allCases.length !== manifest.case_count) {
    throw new Error(`Manifest count mismatch: expected ${manifest.case_count}, resolved ${allCases.length}`);
  }
  let selected = allCases;
  if (args.cases.length) {
    const wanted = new Set(args.cases);
    selected = allCases.filter((item) => wanted.has(item.id));
    const found = new Set(selected.map((item) => item.id));
    const missing = args.cases.filter((id) => !found.has(id));
    if (missing.length) throw new Error(`Unknown case ID(s): ${missing.join(", ")}`);
  }
  if (args.stopAfter !== null) selected = selected.slice(0, args.stopAfter);

  const edgePath = findEdge();
  const outputRoot = path.resolve(REPO_ROOT, manifest.output_root);
  const runtimeDir = path.join(outputRoot, "_runtime");
  fs.mkdirSync(runtimeDir, { recursive: true });
  const startSignature = captureRuntimeSignature(edgePath, true);
  assertFrozenContracts(manifest, startSignature, runtimeRevision);
  writeJson(path.join(outputRoot, "batch_start_runtime_signature.json"), startSignature);
  const protectedSnapshots = {};
  for (const rel of Object.values(manifest.parents).map((parent) => parent.path)) {
    protectedSnapshots[rel] = sha256File(path.join(REPO_ROOT, rel));
  }
  writeJson(path.join(outputRoot, "protected_input_hashes.json"), protectedSnapshots);

  let server;
  let edge;
  const results = [];
  try {
    server = await startServer(runtimeDir);
    activeServerUrl = server.url;
    edge = await startEdge(edgePath, server.url, runtimeDir);
    for (let index = 0; index < selected.length; index += 1) {
      log(`Batch progress ${index + 1}/${selected.length}`, selected[index].id);
      const result = await executeCase(edge.cdp, selected[index], outputRoot, startSignature);
      results.push(result);
      const summary = {
        schema: "local-web-output-refresh-progress/v1",
        manifest: relativePosix(MANIFEST_PATH),
        batch_id: manifest.batch_id,
        updated_at: nowIso(),
        requested_case_count: selected.length,
        completed: results.filter((item) => item.status === "complete").length,
        blocked: results.filter((item) => item.status === "blocked").length,
        results,
      };
      writeJson(path.join(outputRoot, "progress.json"), summary);
    }
  } finally {
    await shutdown(server, edge);
  }

  const endSignature = captureRuntimeSignature(edgePath, false);
  writeJson(path.join(outputRoot, "batch_end_runtime_signature.json"), endSignature);
  const inputDrift = [];
  for (const [rel, original] of Object.entries(protectedSnapshots)) {
    const actual = sha256File(path.join(REPO_ROOT, rel));
    if (actual !== original) inputDrift.push({ path: rel, start: original, end: actual });
  }
  const batchSummary = {
    schema: "local-web-output-refresh-summary/v1",
    batch_id: manifest.batch_id,
    finished_at: nowIso(),
    selected_case_count: selected.length,
    complete_count: results.filter((item) => item.status === "complete").length,
    blocked_count: results.filter((item) => item.status === "blocked").length,
    skipped_complete_count: results.filter((item) => item.skipped).length,
    source_signature_stable: startSignature.source_tree.sha256 === endSignature.source_tree.sha256,
    protected_inputs_stable: inputDrift.length === 0,
    protected_input_drift: inputDrift,
    results,
  };
  writeJson(path.join(outputRoot, "batch_summary.json"), batchSummary);
  log("Batch invocation finished", batchSummary);
  if (batchSummary.blocked_count || !batchSummary.source_signature_stable || !batchSummary.protected_inputs_stable) {
    process.exitCode = 2;
  }
}

// Importing the signature helpers must never launch a browser or change outputs.
export { sourceSignature, assertFrozenContracts, WEB_CONTRACT_PATHS };

if (process.argv[1] && path.resolve(process.argv[1]) === SCRIPT_PATH) {
  main().catch((error) => {
    process.stderr.write(`${error.stack || error.message}\n`);
    process.exitCode = 1;
  });
}
