"""Capture only explicitly listed public ISODISTORT pages; no user inputs."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence" / "web"
BASE = "https://landau3.byu.edu/"
FILES = (
    "isodistort.php", "isodistorthelp.php", "isodistorttutorials.php",
    "isodistortversions.php", "isodistortparent.html", "isodistortdistortion.html",
    "isodistort.js",
)


class StructureParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.forms = []
        self.links = []
        self.text = []
        self.form = None

    def handle_starttag(self, tag, attrs):
        attr = dict(attrs)
        if tag == "form":
            self.form = {"attributes": attr, "controls": []}
            self.forms.append(self.form)
        if self.form is not None and tag in ("input", "select", "option", "textarea"):
            self.form["controls"].append({"tag": tag, **attr})
        if tag == "a" and "href" in attr:
            self.links.append(attr)
        if tag in ("p", "br", "h1", "h2", "h3", "li", "tr", "div"):
            self.text.append("\n")

    def handle_endtag(self, tag):
        if tag == "form":
            self.form = None

    def handle_data(self, data):
        self.text.append(data)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    records = []
    for name in FILES:
        url = BASE + name
        with urlopen(Request(url, headers={"User-Agent": "CRIS-REA-research/0.1"}), timeout=60) as response:
            body = response.read()
            metadata = {
                "url": url, "final_url": response.url, "http_status": response.status,
                "captured_at_utc": datetime.now(timezone.utc).isoformat(),
                "sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body),
                "headers": {k: v for k, v in response.headers.items()
                            if k.lower() in {"content-type", "last-modified", "etag", "date"}},
                "path": f"evidence/web/{name}",
                "role": "official-reference-only",
            }
        (OUT / name).write_bytes(body)
        if not name.endswith(".js"):
            parser = StructureParser()
            parser.feed(body.decode("utf-8", errors="replace"))
            (OUT / f"{name}.structure.json").write_text(
                json.dumps({"forms": parser.forms, "links": parser.links}, indent=2), encoding="utf-8")
            (OUT / f"{name}.txt").write_text("".join(parser.text), encoding="utf-8")
        records.append(metadata)
    (OUT / "public_pages_manifest.json").write_text(
        json.dumps({"capture_method": "Python urllib HTTPS GET; decoded response body",
                    "records": records}, indent=2), encoding="utf-8")
    print(json.dumps([{"file": r["path"], "bytes": r["bytes"]} for r in records], indent=2))


if __name__ == "__main__":
    main()
