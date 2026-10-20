"""Extract browser captures with Evidence-ID and transaction provenance."""
import base64
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence" / "web"


def main():
    bundle = OUT / "rea_method2_bundle.json"
    data = json.loads(bundle.read_text(encoding="utf-8"))
    manifest = []
    for record in data["records"]:
        if record["operation"] != "capture_browser_scenario":
            continue
        evidence_id = record["evidence_id"]
        target = OUT / evidence_id
        target.mkdir(exist_ok=True)
        (target / "record.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
        normalized = record["normalized_result"]
        for step in normalized["steps"]:
            for kind in ("dom", "accessibility"):
                artifact = step["artifacts"][kind]
                if artifact["state"] == "captured":
                    name = f"step-{step['step_index']}-{kind}.{'html' if kind == 'dom' else 'txt'}"
                    (target / name).write_text(artifact["value"]["text"], encoding="utf-8")
        for event in normalized["events"]["items"]:
            body = event.get("body", {})
            if event.get("kind") != "network-content" or body.get("state") != "captured":
                continue
            if body["encoding"] != "base64":
                raise ValueError(body["encoding"])
            raw = base64.b64decode(body["content"], validate=True)
            assert hashlib.sha256(raw).hexdigest() == body["sha256"]
            suffix = ".html" if body.get("media_type", "").startswith("text/html") else ".bin"
            name = f"{event['transaction_id']}-{event['phase']}{suffix}"
            (target / name).write_bytes(raw)
            manifest.append({"evidence_id": evidence_id, "transaction_id": event["transaction_id"],
                             "phase": event["phase"], "sha256": body["sha256"],
                             "path": str((target / name).relative_to(ROOT)).replace("\\", "/"),
                             "media_type": body.get("media_type"), "bytes": len(raw),
                             "role": "official-reference-only"})
    (OUT / "rea_extraction_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Extracted {len(manifest)} bodies with verified SHA-256.")


if __name__ == "__main__":
    main()
