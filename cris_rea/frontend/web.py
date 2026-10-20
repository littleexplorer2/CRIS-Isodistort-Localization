import base64
import binascii
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import secrets
from urllib.parse import urlsplit

from .webpage import PAGE

MAX_BODY_BYTES = 1500000


def make_handler(service, token):
    class Handler(BaseHTTPRequestHandler):
        server_version = "CRIS-REA/0.1"
        read_timeout = 10

        def setup(self):
            super().setup()
            self.connection.settimeout(self.read_timeout)

        def log_message(self, fmt, *args):
            # Request inputs and scientific results are not logged to console.
            return

        def _reply(self, status, payload, content_type="application/json; charset=utf-8"):
            raw = payload.encode("utf-8") if isinstance(payload, str) else \
                json.dumps(payload, ensure_ascii=False).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(raw)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "no-referrer")
            self.send_header("X-Frame-Options", "DENY")
            self.send_header(
                "Content-Security-Policy",
                "default-src 'none'; script-src 'nonce-" + token +
                "'; style-src 'nonce-" + token +
                "'; connect-src 'self'; base-uri 'none'; frame-ancestors 'none'; form-action 'none'",
            )
            self.end_headers()
            self.wfile.write(raw)

        def _host_allowed(self):
            allowed = {
                f"127.0.0.1:{self.server.server_port}",
                f"localhost:{self.server.server_port}",
            }
            return self.headers.get("Host") in allowed

        def _authorized(self):
            value = self.headers.get("X-CRIS-Token", "")
            origin = self.headers.get("Origin")
            expected = {
                f"http://127.0.0.1:{self.server.server_port}",
                f"http://localhost:{self.server.server_port}",
            }
            return secrets.compare_digest(value, token) and \
                (origin is None or origin in expected)

        def do_GET(self):
            if not self._host_allowed():
                self._reply(403, {"error": "Invalid local host."})
                return
            path = urlsplit(self.path).path
            if path == "/":
                self._reply(200, PAGE.replace("__TOKEN__", token), "text/html; charset=utf-8")
                return
            if not self._authorized():
                self._reply(403, {"error": "A local session token is required."})
                return
            if path == "/api/programs":
                self._reply(200, service.programs())
            elif path.startswith("/api/runs/"):
                try:
                    self._reply(200, service.result(path.removeprefix("/api/runs/")))
                except (ValueError, OSError):
                    self._reply(404, {"error": "Run not found."})
            else:
                self._reply(404, {"error": "Not found."})

        def do_POST(self):
            if not self._host_allowed() or not self._authorized():
                self._reply(403, {"error": "Invalid local session or origin."})
                return
            path = urlsplit(self.path).path
            if path not in ("/api/run", "/api/iso-query"):
                self._reply(404, {"error": "Not found."})
                return
            if self.headers.get("Transfer-Encoding"):
                self._reply(400, {"error": "Transfer-Encoding is not supported."})
                return
            if self.headers.get("Content-Type", "").split(";")[0] != "application/json":
                self._reply(415, {"error": "Expected application/json."})
                return
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if length <= 0 or length > MAX_BODY_BYTES:
                    self._reply(413, {"error": "Request is empty or exceeds the size limit."})
                    return
                request = json.loads(self.rfile.read(length))
                if not isinstance(request, dict):
                    raise ValueError("Expected a JSON object.")
                if path == "/api/iso-query":
                    self._reply(200, service.prepare_iso_query(
                        request.get("parent"), request.get("irrep"),
                        request.get("opd"), request.get("wyckoff"),
                    ))
                    return
                if ("input" in request) == ("input_b64" in request):
                    raise ValueError("Provide either text input or one uploaded file.")
                if "input_b64" in request:
                    data = base64.b64decode(request["input_b64"], validate=True)
                else:
                    data = request["input"].encode("utf-8")
                result = service.run(
                    request.get("program"), data, request.get("timeout", 30)
                )
                self._reply(200, result)
            except (ValueError, TypeError, AttributeError, binascii.Error) as exc:
                self._reply(400, {"error": str(exc)})
            except TimeoutError:
                self.close_connection = True
                self._reply(408, {"error": "Request body was not received in time."})
            except OSError as exc:
                self._reply(500, {"error": f"Local file operation failed: {exc}"})

    return Handler


def create_server(service, port=8765):
    if not isinstance(port, int) or not 0 <= port <= 65535:
        raise ValueError("Port must be between 0 and 65535.")
    token = secrets.token_hex(24)
    server = ThreadingHTTPServer(("127.0.0.1", port), make_handler(service, token))
    server.daemon_threads = False
    return server


def serve(service, port=8765):
    server = create_server(service, port)
    print(f"CRIS REA: http://127.0.0.1:{server.server_port}", flush=True)
    print("Local official-tool console. Ctrl+C stops the server.", flush=True)
    try:
        server.serve_forever()
    finally:
        server.server_close()
