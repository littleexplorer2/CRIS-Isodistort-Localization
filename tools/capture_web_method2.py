"""Auditable official-web protocol replay. Never imported by local computation."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence" / "web" / "method2"
BASE = "https://landau3.byu.edu/"


class Forms(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.forms = []
        self.current = None
        self.select = None
        self.option = None
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "form":
            self.current = {"attributes": a, "inputs": [], "selects": []}
            self.forms.append(self.current)
        elif self.current is not None:
            if tag == "input":
                self.current["inputs"].append(a)
            elif tag == "select":
                self.select = {**a, "options": []}
                self.current["selects"].append(self.select)
            elif tag == "option" and self.select is not None:
                self.option = {**a, "text": ""}
                self.select["options"].append(self.option)

    def handle_endtag(self, tag):
        if tag == "form":
            self.current = None
        elif tag == "select":
            self.select = self.option = None
        elif tag == "option":
            self.option = None

    def handle_data(self, data):
        if self.option is not None:
            self.option["text"] += data


def fields(form):
    result = {}
    for field in form["inputs"]:
        kind = field.get("type", "").lower()
        if not field.get("name") or kind in {"submit", "button", "file"}:
            continue
        if kind in {"checkbox", "radio"} and "checked" not in field:
            continue
        result[field["name"]] = field.get("value", "")
    for select in form["selects"]:
        chosen = next((o for o in select["options"] if "selected" in o), select["options"][0])
        result[select["name"]] = chosen.get("value", chosen["text"])
    return result


def select_form(html, input_value):
    forms = Forms(html).forms
    matches = [f for f in forms
               if any(x.get("name") == "input" and x.get("value") == input_value for x in f["inputs"])]
    if len(matches) != 1:
        raise ValueError(f"Expected one form for {input_value}; found {len(matches)}")
    return matches[0]


def post(form, name, overrides=None):
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / f"{name}.html").exists():
        raise FileExistsError(f"Do not overwrite captured response: {name}")
    url = urljoin(BASE, form["attributes"]["action"])
    if url != BASE + "isodistortform.php":
        raise ValueError("Replay restricted to observed official form endpoint")
    values = fields(form)
    values.update(overrides or {})
    request_bytes = urlencode(values).encode()
    (OUT / f"{name}.request.txt").write_bytes(request_bytes)
    (OUT / f"{name}.fields.json").write_text(json.dumps(values, indent=2), encoding="utf-8")
    started = datetime.now(timezone.utc).isoformat()
    with urlopen(Request(url, data=request_bytes,
                         headers={"Content-Type": "application/x-www-form-urlencoded",
                                  "User-Agent": "CRIS-REA-research/0.1"}), timeout=60) as response:
        raw = response.read()
        metadata = {"url": url, "final_url": response.url, "method": "POST", "status": response.status,
                    "content_type": response.headers.get("Content-Type"),
                    "started_at_utc": started, "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                    "request_sha256": hashlib.sha256(request_bytes).hexdigest(),
                    "response_sha256": hashlib.sha256(raw).hexdigest(),
                    "request_bytes": len(request_bytes), "response_bytes": len(raw),
                    "role": "official-reference-only",
                    "limitations": ["Decoded HTTP response; no server-internal process visibility."]}
    (OUT / f"{name}.html").write_bytes(raw)
    html = raw.decode("utf-8")
    (OUT / f"{name}.forms.json").write_text(json.dumps(Forms(html).forms, indent=2), encoding="utf-8")
    (OUT / f"{name}.metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return html


def main():
    # Public pedagogical SrTiO3 parent; never uploads private or experimental data.
    html = (ROOT / "evidence/web/isodistortparent.html").read_text(encoding="utf-8")
    html = post(select_form(html, "autoloadparent"), "01_search")
    html = post(select_form(html, "kvector"), "02_irrep")
    form = select_form(html, "irrep")
    irrep = next(o for s in form["selects"] if s["name"] == "irrep1"
                 for o in s["options"] if o["text"].startswith("GM4-"))
    html = post(form, "03_opd", {"irrep1": irrep["value"]})
    form = select_form(html, "distort")
    if "P1      (a,0,0)" not in fields(form)["orderparam"]:
        raise ValueError("Official default OPD changed; inspect before continuing")
    html = post(form, "04_distortion")
    post(select_form(html, "displaydistort"), "05_complete_modes",
         {"origintype": "completemodesdetails"})
    print("Captured public parent > GM > GM4- > P1 > complete modes details.")


if __name__ == "__main__":
    main()
