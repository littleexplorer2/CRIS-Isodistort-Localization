"""Produce token-free derived browser evidence; retain raw bundles only locally.

This extracts screenshot bytes without editing images. Raw HTML and scripts are
deliberately excluded. No production computation imports this utility.
"""
from __future__ import annotations

import base64
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "evidence/validation/live-http"
SCENARIOS = (
    ("success", "rea-local-browser-bundle.json",
     "ev_5308c0d89e7d80cb537e045f28c06e98e49df50bc000e63b8e3064cffbda39d4",
     "8a69195a-4f22-4574-b980-451275c41fb5"),
    ("negative", "rea-local-negative-bundle.json",
     "ev_c2a71bb109f6829f44b92d418433920119e504aca8c4519deebfbcd4404a784c",
     "6a34fad6-bf5d-4157-87ca-086bc24adc13"),
)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def run():
    results = []
    for label, filename, evidence_id, run_id in SCENARIOS:
        source = BASE / filename
        raw_bundle = source.read_bytes()
        bundle = json.loads(raw_bundle)
        record = next(r for r in bundle["records"] if r["evidence_id"] == evidence_id)
        normalized = record["normalized_result"]
        tokens = set()
        token_pattern = r"""const\s+token\s*=\s*["']([^"']+)"""
        for candidate in bundle["records"]:
            captured = candidate.get("normalized_result", {})
            for step in captured.get("steps", []):
                dom = step["artifacts"]["dom"]
                if dom["state"] == "captured":
                    tokens.update(re.findall(token_pattern, dom["value"]["text"]))
            for event in captured.get("events", {}).get("items", []):
                body = event.get("body", {})
                if body.get("state") == "captured" and body.get("media_type", "").startswith("text/html"):
                    html = base64.b64decode(body["content"], validate=True).decode("utf-8")
                    tokens.update(re.findall(token_pattern, html))
        if not tokens:
            raise ValueError("Expected a known local session token to remove")

        def scrub(text):
            for token in tokens:
                text = text.replace(token, "[LOCAL_SESSION_TOKEN_REMOVED]")
            return text

        output = BASE / label
        output.mkdir(exist_ok=True)
        exported = []

        def write(name, raw, original_sha=None, transformation="none"):
            if any(token.encode() in raw for token in tokens):
                raise ValueError("Refusing to export a local session token")
            destination = output / name
            destination.write_bytes(raw)
            exported.append({
                "path": destination.relative_to(ROOT).as_posix(),
                "bytes": len(raw), "sha256": sha(raw),
                "source_sha256": original_sha or sha(raw),
                "transformation": transformation,
            })

        final = normalized["steps"][-1]
        for kind in ("accessibility", "screenshot"):
            artifact = final["artifacts"][kind]
            if artifact["state"] != "captured":
                continue
            value = artifact["value"]
            if kind == "screenshot":
                raw = base64.b64decode(value["data_base64"], validate=True)
                assert sha(raw) == value["sha256"]
                assert value["media_type"] == "image/png"
                assert raw.startswith(b"\x89PNG\r\n\x1a\n")
                write("final.png", raw, value["sha256"], "base64 decode only; image unedited")
            else:
                raw = value["text"].encode("utf-8")
                assert sha(raw) == value["sha256"]
                write("accessibility.txt", scrub(value["text"]).encode("utf-8"),
                      value["sha256"], "remove known session-token strings if present")

        events = normalized["events"]["items"]
        transactions = {}
        for event in events:
            tx = event.get("transaction_id")
            if tx and event["kind"] in ("request", "response", "request-finished", "request-failed"):
                item = transactions.setdefault(tx, {"transaction_id": tx})
                item.update({key: event[key] for key in
                             ("method", "url", "resource_type") if key in event})
                if event.get("status") is not None:
                    item["status"] = event["status"]
                if event.get("failure") is not None:
                    item["failure"] = event["failure"]
        unavailable = []
        for event in events:
            if event["kind"] != "network-content":
                continue
            tx, phase, body = event["transaction_id"], event["phase"], event["body"]
            item = transactions.setdefault(tx, {"transaction_id": tx})
            body_summary = {k: v for k, v in body.items() if k != "content"}
            if body.get("state") == "captured":
                data = base64.b64decode(body["content"], validate=True)
                assert sha(data) == body["sha256"]
                media = body.get("media_type", "")
                # Keep only observed API JSON, never HTML, scripts, CSS, or headers.
                path = urlsplit(item.get("url", {}).get("url", "")).path
                if media.startswith("application/json") and path in {
                    "/api/programs", "/api/iso-query", "/api/run"
                }:
                    body_summary["json"] = json.loads(scrub(data.decode("utf-8")))
                else:
                    body_summary["content_omitted"] = "HTML/scripts/non-API data excluded"
            else:
                unavailable.append({
                    "transaction_id": tx, "phase": phase,
                    "url": item.get("url"), "body_state": body,
                })
            item[f"{phase}_body"] = body_summary

        run_dir = ROOT / "runs" / run_id
        evidence = json.loads((run_dir / "evidence.json").read_text(encoding="utf-8"))
        native_input = (run_dir / "input.dat").read_bytes()
        run_tx = next(t for t in transactions.values()
                      if urlsplit(t.get("url", {}).get("url", "")).path == "/api/run")
        posted = run_tx["request_body"]["json"]
        posted_bytes = (base64.b64decode(posted["input_b64"], validate=True)
                        if "input_b64" in posted else posted["input"].encode("utf-8"))
        response = run_tx["response_body"]["json"]
        assert posted_bytes == native_input
        assert evidence["input"]["sha256"] == sha(native_input)
        assert response["evidence"] == evidence
        assert response["stdout"] == (run_dir / "stdout.txt").read_bytes().decode("utf-8")
        assert response["stderr"] == (run_dir / "stderr.txt").read_bytes().decode("utf-8")
        source_matches = all(sha((ROOT / name).read_bytes()) == value["sha256"]
                             for name, value in evidence["source_code"].items())
        for name in ("input.dat", "stdout.txt", "stderr.txt", "evidence.json",
                     "execution.json", "process.json",
                     "supervisor.stdout.txt", "supervisor.stderr.txt"):
            raw = (run_dir / name).read_bytes()
            write(name, raw)

        summary = {
            "scope": "Real local browser/HTTP/Windows/WSL/official ISO execution; not scientific validation",
            "label": label, "run_id": run_id, "evidence_id": evidence_id,
            "derived_at_utc": datetime.now(timezone.utc).isoformat(),
            "raw_bundle": {"path": source.relative_to(ROOT).as_posix(),
                           "sha256": sha(raw_bundle), "bytes": len(raw_bundle),
                           "retained_locally_only": True, "contains_session_token": True},
            "redaction": {"known_token_count": len(tokens), "tokens_exported": False,
                          "raw_dom_exported": False, "full_scripts_exported": False,
                          "header_values_exported": False,
                          "screenshot": "exact captured PNG bytes, not edited",
                          "derivative_is_original_REA_record": False},
            "provider": record["provider"], "browser": normalized["browser"],
            "steps": [{k: s[k] for k in ("step_id", "action", "status", "elapsed_ms", "error")}
                      for s in normalized["steps"]],
            "completeness": normalized["completeness"],
            "events_retained": normalized["events"]["retained"],
            "events_dropped": normalized["events"]["dropped"],
            "event_omissions": normalized["events"]["omissions"],
            "noncaptured_network_bodies": unavailable,
            "network": list(transactions.values()),
            "limitations": record["limitations"],
            "cross_checks": {"posted_bytes_equal_native_input": True,
                             "input_hash_matches_manifest": True,
                             "http_evidence_equals_saved_manifest": True,
                             "http_outputs_equal_saved_raw_utf8_outputs": True,
                             "source_fingerprints_match_current_disk": source_matches,
                             "source_fingerprints_prove_loaded_code": False},
            "status": evidence["status"],
            "official_exit_code": evidence["process"]["exit_code"],
            "scientific_validation": evidence["scientific_validation"],
            "exports": exported,
            "omitted_run_files": "Official executables, databases, work directory and worker copy are not duplicated",
        }
        encoded = json.dumps(summary, ensure_ascii=False, indent=2).encode("utf-8")
        assert not any(token.encode() in encoded for token in tokens)
        (output / "summary.json").write_bytes(encoded)
        results.append({"label": label, "status": evidence["status"],
                        "official_exit_code": evidence["process"]["exit_code"],
                        "source_fingerprints_match": source_matches,
                        "noncaptured_bodies": len(unavailable),
                        "events_dropped": normalized["events"]["dropped"],
                        "export_count": len(exported)})
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    run()
