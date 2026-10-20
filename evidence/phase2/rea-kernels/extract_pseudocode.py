"""Extract exact decompiler text from saved REA records; never execute it."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]


def main():
    entries = []
    for path in sorted(BASE.glob("*.json")):
        if not re.match(r"(?:0[4-9]|1[0-2])-", path.name):
            continue
        raw = path.read_bytes()
        outer = json.loads(raw)
        tool = outer["response"]["result"]
        record = json.loads(tool["content"][0]["text"])
        text = record["normalized_result"]
        if not isinstance(text, str):
            raise TypeError((path.name, type(text).__name__))
        destination = BASE / (path.stem + ".pseudocode.txt")
        encoded = text.encode("utf-8")
        destination.write_bytes(encoded)
        entries.append({
            "raw_path": path.relative_to(ROOT).as_posix(),
            "raw_sha256": hashlib.sha256(raw).hexdigest(),
            "extracted_path": destination.relative_to(ROOT).as_posix(),
            "extracted_sha256": hashlib.sha256(encoded).hexdigest(),
            "extracted_bytes": len(encoded),
            "evidence_id": record["evidence_id"],
            "operation": record["operation"],
            "parameters": record["parameters"],
            "subject": record["subject"],
            "provider": record["provider"],
            "tool_is_error": bool(tool.get("isError")),
            "confidence": record["confidence"],
            "authority": record["authority"],
            "limitations": record["limitations"],
            "raw_normalized_text_equal": text == record["raw_result"]["value"],
        })
    result = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "role": "exact text extraction of Ghidra pseudocode, not original source or proof",
        "entries": entries,
    }
    (BASE / "pseudocode_manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps([{"procedure": e["parameters"]["procedure"],
                       "bytes": e["extracted_bytes"], "is_error": e["tool_is_error"]}
                      for e in entries], indent=2))


if __name__ == "__main__":
    main()
