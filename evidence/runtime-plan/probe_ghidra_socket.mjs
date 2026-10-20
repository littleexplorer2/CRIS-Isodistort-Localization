// Transport-only probe. Does not start Ghidra or analyze an executable.
import {createServer} from "node:net";
import {lstat, mkdir, writeFile} from "node:fs/promises";
import {fileURLToPath} from "node:url";
import {dirname, join, resolve} from "node:path";
import {createGhidraEndpoint} from "./rea/node_modules/rea-agents/dist/ghidra/GhidraTransport.js";

const runtime = dirname(fileURLToPath(import.meta.url));
const root = resolve(runtime, "../..");
const candidate = join(runtime, "temp/ghidra-projects-and-analysis-artifacts/rea-ghidra-PROBE0");
const record = {
  purpose: "Test REA's existing Unix socket fallback only; no Ghidra analysis.",
  started_at_utc: new Date().toISOString(),
  rea_version: "6.7.0",
  runtime_root_candidate: candidate,
  unshortened_socket_bytes: Buffer.byteLength(join(candidate, "bridge.sock")),
  platform: process.platform,
  node: process.version,
};
let allocated;
let server;
try {
  allocated = await createGhidraEndpoint(candidate, "unix-socket", process.platform);
  record.endpoint = allocated.endpoint;
  record.socket_root = allocated.socketRoot?.path ?? null;
  if (!allocated.socketRoot || !allocated.endpoint.path.startsWith("/tmp/rea-ghidra-socket-")) {
    throw new Error("Expected the supported private /tmp socket fallback.");
  }
  record.socket_root_mode = ((await lstat(allocated.socketRoot.path)).mode & 0o777).toString(8);
  server = createServer();
  await new Promise((ok, fail) => {
    server.once("error", fail);
    server.listen(allocated.endpoint.path, ok);
  });
  record.bind_succeeded = true;
  await new Promise((ok, fail) => server.close(error => error ? fail(error) : ok()));
  server = null;
  await allocated.socketRoot.close();
  record.cleaned = await lstat(allocated.socketRoot.path).then(() => false, error => error.code === "ENOENT");
  if (!record.cleaned) throw new Error("Private probe socket directory remains.");
  record.status = "passed";
} catch (error) {
  record.status = "failed";
  record.error = String(error);
  process.exitCode = 1;
} finally {
  if (server?.listening) await new Promise(ok => server.close(ok));
  if (allocated?.socketRoot) await allocated.socketRoot.close();
  record.finished_at_utc = new Date().toISOString();
  const directory = join(root, "evidence/environment");
  await mkdir(directory, {recursive: true});
  await writeFile(join(directory, "ghidra-socket-transport-probe.json"), JSON.stringify(record, null, 2) + "\n");
  console.log(JSON.stringify(record));
}
