"""Independent reference capture after freezing a completed local calculation.

Never imported by the product. Each stage is inspectable, append-only, and uses
only forms returned by the official site. No reference values enter local modes.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import Request, urlopen
from urllib.parse import urljoin

import capture_web_method2 as protocol

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "evidence/phase2/workflow/20261020T145121147627Z/cubic_polar_axis"
OUT = ROOT / "evidence/phase2/web-comparison/20261020_cubic_polar_axis"
protocol.OUT = OUT


def fingerprint(path):
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def write_json(path, value):
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def freeze():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = sorted(p for p in LOCAL.rglob("*") if p.is_file())
    if not paths or not (LOCAL / "parent.cif").is_file():
        raise RuntimeError("Completed local evidence must exist before reference access.")
    write_json(OUT / "local_freeze_manifest.json", {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "role": "Independent completed local calculation frozen before this website comparison",
        "input_origin": "Independently constructed cubic BaTiO3 geometry, a=4 angstrom; not experimental",
        "requested_inputs": {"irrep": "GM4-", "opd": "P1", "k": [0, 0, 0]},
        "local_directory": LOCAL.relative_to(ROOT).as_posix(),
        "files": {p.relative_to(ROOT).as_posix(): fingerprint(p) for p in paths},
        "capture_script": fingerprint(Path(__file__)),
        "reference_role": "comparison_only_never_read_by_local_calculation",
    })
    print("Frozen", len(paths), "local files before reference requests.")


def verify_freeze():
    record = json.loads((OUT / "local_freeze_manifest.json").read_text())
    for path, expected in record["files"].items():
        if fingerprint(ROOT / path) != expected:
            raise RuntimeError(f"Frozen local file changed: {path}")


def capture_http(name, url, data=None, content_type=None):
    if (OUT / f"{name}.html").exists():
        raise FileExistsError(name)
    if data is not None:
        (OUT / f"{name}.request.bin").write_bytes(data)
    started = datetime.now(timezone.utc).isoformat()
    headers = {"User-Agent": "CRIS-REA-independent-comparison/0.1"}
    if content_type:
        headers["Content-Type"] = content_type
    with urlopen(Request(url, data=data, headers=headers), timeout=60) as response:
        raw = response.read()
        metadata = {
            "url": url, "final_url": response.url, "method": "POST" if data else "GET",
            "http_status": response.status, "content_type": response.headers.get("Content-Type"),
            "started_at_utc": started, "finished_at_utc": datetime.now(timezone.utc).isoformat(),
            "response_sha256": hashlib.sha256(raw).hexdigest(), "response_bytes": len(raw),
            "request_sha256": hashlib.sha256(data).hexdigest() if data else None,
            "request_bytes": len(data) if data else 0,
            "role": "official-reference-only; never computation input",
            "visibility": "Decoded HTTPS response; no server-internal process observation",
        }
    (OUT / f"{name}.html").write_bytes(raw)
    html = raw.decode("utf-8")
    write_json(OUT / f"{name}.forms.json", protocol.Forms(html).forms)
    write_json(OUT / f"{name}.metadata.json", metadata)
    return html


def upload():
    html = capture_http("00_home", protocol.BASE + "isodistort.php")
    form = protocol.select_form(html, "uploadparentcif")
    action = urljoin(protocol.BASE, form["attributes"]["action"])
    if action != protocol.BASE + "isodistortuploadfile.php":
        raise ValueError("Unexpected upload action.")
    file_controls = [x for x in form["inputs"] if x.get("type") == "file"]
    if len(file_controls) != 1 or file_controls[0].get("name") != "toProcess":
        raise ValueError("Unexpected file field.")
    values = protocol.fields(form)
    boundary = "CRIS_REA_phase2_synthetic_cif_20261020"
    body = bytearray()
    for key, value in values.items():
        body += f'--{boundary}\r\nContent-Disposition: form-data; name="{key}"\r\n\r\n{value}\r\n'.encode()
    cif = (LOCAL / "parent.cif").read_bytes()
    body += (
        f'--{boundary}\r\nContent-Disposition: form-data; name="toProcess"; filename="parent.cif"\r\n'
        "Content-Type: chemical/x-cif\r\n\r\n"
    ).encode() + cif + f"\r\n--{boundary}--\r\n".encode()
    write_json(OUT / "01_upload.fields.json", {
        "form_fields": values, "file_field": "toProcess", "file_sha256": hashlib.sha256(cif).hexdigest(),
        "file_source": (LOCAL / "parent.cif").relative_to(ROOT).as_posix(),
    })
    capture_http("01_upload", action, bytes(body), f"multipart/form-data; boundary={boundary}")


def read(name):
    return (OUT / f"{name}.html").read_text(encoding="utf-8")


def main():
    command = sys.argv[1]
    if command == "freeze":
        freeze()
        return
    verify_freeze()
    if command == "upload":
        upload()
    elif command == "process":
        protocol.post(protocol.select_form(read("01_upload"), "uploadparentcif"), "01_processed")
    elif command == "k":
        form = protocol.select_form(read("01_processed"), "kvector")
        choices = [o for s in form["selects"] if s["name"] == "kvec1"
                   for o in s["options"] if o["text"].startswith("GM,")]
        if len(choices) != 1:
            raise ValueError("Unique GM option not found.")
        protocol.post(form, "02_irrep", {"kvec1": choices[0].get("value", choices[0]["text"])})
    elif command == "irrep":
        form = protocol.select_form(read("02_irrep"), "irrep")
        choices = [o for s in form["selects"] if s["name"] == "irrep1"
                   for o in s["options"] if o["text"].startswith("GM4-")]
        if len(choices) != 1:
            raise ValueError("Unique GM4- option not found.")
        protocol.post(form, "03_opd", {"irrep1": choices[0].get("value", choices[0]["text"])})
    elif command == "opd":
        form = protocol.select_form(read("03_opd"), "distort")
        if "P1      (a,0,0)" not in protocol.fields(form)["orderparam"]:
            raise ValueError("Inspect current OPD page; P1 must be explicitly selected.")
        protocol.post(form, "04_distortion")
    elif command == "complete":
        protocol.post(protocol.select_form(read("04_distortion"), "displaydistort"),
                      "05_complete_modes", {"origintype": "completemodesdetails"})
    else:
        raise ValueError("Unknown stage.")
    verify_freeze()
    write_json(OUT / f"stage_{command}.json", {
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "capture_script_sha256": fingerprint(Path(__file__))["sha256"],
        "frozen_local_files_reverified": True,
    })
    print("Completed", command, "with local frozen files unchanged.")


if __name__ == "__main__":
    main()
