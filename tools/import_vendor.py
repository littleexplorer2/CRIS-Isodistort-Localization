"""Copy explicit local ISOTROPY distribution into this independent project."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target = root / "vendor/isobyu"
    source = args.source.resolve(strict=True)
    if source == target.resolve():
        parser.error("source and destination must differ")
    expected = json.loads((root / "evidence/environment/source_manifest.json").read_text(encoding="utf-8-sig"))
    for entry in expected:
        path = source / entry["name"]
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            parser.error("missing or different artifact: " + entry["name"])
        dest = target / entry["name"]
        if dest.exists() and hashlib.sha256(dest.read_bytes()).hexdigest() != entry["sha256"]:
            parser.error("refusing to overwrite different artifact: " + str(dest))
    target.mkdir(parents=True, exist_ok=True)
    for entry in expected:
        dest = target / entry["name"]
        if not dest.exists():
            shutil.copy2(source / entry["name"], dest)
    print(f"Verified {len(expected)} artifacts in {target}")

if __name__ == "__main__":
    main()
