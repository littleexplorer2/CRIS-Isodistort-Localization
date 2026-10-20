"""Extract readable derived pseudocode with links to its original MCP Evidence."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "evidence/native/rea-linux-iso-recovered"
DEST = ROOT / "evidence/native/pseudocode"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    DEST.mkdir(exist_ok=True)
    manifest = []
    for path in sorted(SOURCE.glob("[0-9][0-9]-*.json")):
        raw = path.read_bytes()
        envelope = json.loads(raw)
        result = envelope.get("response", {}).get("result", {})
        if result.get("isError"):
            continue
        evidence = result.get("structuredContent")
        if evidence is None:
            texts = [x["text"] for x in result.get("content", []) if x.get("type") == "text"]
            if not texts:
                continue
            evidence = json.loads(texts[0])
        if not isinstance(evidence, dict) or not evidence.get("evidence_id"):
            continue
        normalized = evidence.get("normalized_result")
        code = normalized if isinstance(normalized, str) else (
            normalized.get("pseudocode") if isinstance(normalized, dict) else None
        )
        if not isinstance(code, str):
            continue
        output = DEST / (envelope["label"] + ".c")
        header = (
            "/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.\n"
            f" * Evidence: {evidence['evidence_id']}\n"
            f" * Raw response: {path.relative_to(ROOT).as_posix()}\n"
            " * Types, array bounds and control flow may be imperfectly recovered.\n"
            " */\n"
        )
        output.write_bytes((header + code).encode("utf-8"))
        manifest.append({
            "source": path.relative_to(ROOT).as_posix(),
            "source_sha256": digest(raw),
            "derived_file": output.relative_to(ROOT).as_posix(),
            "derived_sha256": digest(output.read_bytes()),
            "evidence_id": evidence["evidence_id"],
            "limitations": evidence.get("limitations", []),
        })
    (DEST / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"extracted": len(manifest), "directory": str(DEST)}))


if __name__ == "__main__":
    main()
