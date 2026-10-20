"""Capture public theory/interface definitions only; no website calculation requests."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
SOURCES = (
    ("isodisplace_2006.pdf", "https://iso.byu.edu/jac39_607reprint.pdf",
     "Campbell et al. J. Appl. Cryst. 39, 607-614 (2006)",
     "Sections 2.1-2.4, 4.3-4.5; historical method, not current normalization"),
    ("isotheory.pdf", "https://iso.byu.edu/isotheory.pdf",
     "Stokes and Hatch, Introduction to Isotropy Subgroups and Displacive Phase Transitions",
     "Definitions of representations, displacement space and stabilizer; no case-result reuse"),
    ("versions.html", "https://landau3.byu.edu/isodistortversions.php",
     "ISODISTORT version history", "6.1 orthogonality; 5.5 strain normalization; 5.6 strain mapping"),
    ("gsas2_G2phase_CIF.py", "https://raw.githubusercontent.com/AdvancedPhotonSource/GSAS-II/main/GSASII/imports/G2phase_CIF.py",
     "GSAS-II author-maintained CIF consumer source",
     "Only ISODISTORT parsing routine; third-party interoperability, not official BYU contract"),
)


def main():
    if (HERE / "sources.json").exists():
        raise FileExistsError("This capture already has a ledger; preserve it and use a new evidence directory.")
    records = []
    for filename, url, title, scope in SOURCES:
        request = Request(url, headers={"User-Agent": "CRIS_REA theory-source audit"})
        started = datetime.now(timezone.utc).isoformat()
        try:
            with urlopen(request, timeout=60) as response:
                raw = response.read()
                record = {
                    "id": filename, "url": url, "final_url": response.url,
                    "retrieved_at_utc": started, "status": response.status,
                    "media_type": response.headers.get("Content-Type"),
                    "etag": response.headers.get("ETag"),
                    "last_modified": response.headers.get("Last-Modified"),
                    "path": (HERE / filename).relative_to(ROOT).as_posix(),
                    "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                    "title": title, "intended_read_scope": scope,
                    "role": "public_definition_or_method_source_not_case_answer",
                }
            (HERE / filename).write_bytes(raw)
        except Exception as error:
            record = {"id": filename, "url": url, "retrieved_at_utc": started,
                      "error_type": type(error).__name__, "error": str(error)}
        records.append(record)
    for filename in ("isodistorthelp.php", "isodistorthelp.php.txt"):
        path = ROOT / "evidence/web" / filename
        raw = path.read_bytes()
        records.append({
            "id": filename, "url": "https://landau3.byu.edu/isodistorthelp.php",
            "path": path.relative_to(ROOT).as_posix(), "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "role": "existing_public_help_archive_not_case_answer",
            "original_capture_ledger": "evidence/web/public_pages_manifest.json",
        })
    (HERE / "sources.json").write_text(
        json.dumps({"schema": 1, "sources": records}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps([{"id": r["id"], "status": r.get("status"),
                       "bytes": r.get("bytes"), "error": r.get("error")}
                      for r in records], indent=2))


if __name__ == "__main__":
    main()
