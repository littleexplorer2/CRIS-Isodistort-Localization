"""Read-only authenticated artifact checks for a real UI-created workflow."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
from urllib.request import Request, urlopen
from urllib.error import HTTPError

ROOT = Path(__file__).resolve().parents[1]


def main():
    wid = sys.argv[1]
    if not re.fullmatch(r"[0-9a-f-]{36}", wid):
        raise ValueError("Use a workflow UUID.")
    base = "http://127.0.0.1:8765"
    html = urlopen(base, timeout=10).read().decode()
    token = re.search(r'const token="([a-f0-9]+)"', html).group(1)
    headers = {"X-CRIS-Token": token, "Origin": base}
    dest = ROOT / "evidence" / "phase2" / "browser" / wid
    dest.mkdir()
    source = ROOT / "runs" / wid
    workflow = json.loads((source / "workflow.json").read_text(encoding="utf-8"))
    record = {"workflow_id": wid, "status": workflow["status"], "summary": workflow["summary"],
              "checks": [], "session_token_saved": False}
    for name, expected in workflow["export_fingerprints"].items():
        request = Request(f"{base}/api/modes/{wid}/{name}", headers=headers)
        with urlopen(request, timeout=10) as response:
            data = response.read()
            actual = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
            if actual != expected:
                raise AssertionError("Downloaded artifact differs from frozen output.")
            (dest / name).write_bytes(data)
            record["checks"].append({"name": name, "status": response.status, **actual,
                                     "disposition": response.headers.get("Content-Disposition")})
    for name in ("parent.cif", "workflow.json"):
        shutil.copyfile(source / name, dest / name)
    try:
        urlopen(f"{base}/api/modes/{wid}/modes.json", timeout=10)
    except HTTPError as exc:
        record["missing_token_status"] = exc.code
    assert record["missing_token_status"] == 403
    payload = json.dumps({"cif": (source / "parent.cif").read_text(), "irrep": "GP1GQ1",
                          "opd": "C1", "k_values": []}).encode()
    try:
        urlopen(Request(base+"/api/modes", data=payload,
                        headers={**headers, "Content-Type": "application/json"}), timeout=120)
        raise AssertionError("Missing k parameters were silently accepted.")
    except HTTPError as exc:
        record["missing_k_status"] = exc.code
        record["missing_k_diagnostic"] = json.loads(exc.read())
    assert record["missing_k_status"] == 400
    (dest / "download_record.json").write_text(json.dumps(record, ensure_ascii=False, indent=2)+"\n",
                                              encoding="utf-8")
    print(json.dumps({"record": str(dest / "download_record.json"), "downloads": len(record["checks"]),
                      "missing_token_status": record["missing_token_status"],
                      "missing_k_status": record["missing_k_status"]}))


if __name__ == "__main__":
    main()
