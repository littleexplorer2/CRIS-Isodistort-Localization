"""Read-only verification of the original distribution and Git status.

This tool writes only the destination evidence directory. It never updates the
source repository, its index, regression counter, or official resource files.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence" / "environment"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((OUT / "source_manifest.json").read_text(encoding="utf-8-sig"))
    results = []
    for item in manifest:
        source = Path(item["source"])
        copy = ROOT / "vendor" / "isobyu" / item["name"]
        source_sha, copy_sha = sha(source), sha(copy)
        results.append({
            "name": item["name"],
            "source_sha256": source_sha,
            "copy_sha256": copy_sha,
            "expected_sha256": item["sha256"],
            "matches_frozen_source_and_copy": source_sha == copy_sha == item["sha256"],
        })
    source_root = Path(manifest[0]["source"]).parents[3]
    status = subprocess.run(
        ["git", "-C", str(source_root), "status", "--porcelain=v1"],
        capture_output=True, check=True, encoding="utf-8",
    ).stdout
    (OUT / "source_git_status_after.txt").write_text(status, encoding="utf-8")
    before = (OUT / "source_git_status_before.txt").read_text(encoding="utf-8-sig")
    counter = source_root / ".regression_validation_counter.json"
    frozen_counter = OUT / "source_counter_readonly.json"
    # Snapshot may have newline/encoding normalization. Compare parsed JSON.
    counter_equal = (
        json.loads(counter.read_text(encoding="utf-8-sig"))
        == json.loads(frozen_counter.read_text(encoding="utf-8-sig"))
    )
    record = {
        "at_utc": datetime.now(timezone.utc).isoformat(),
        "source_root": str(source_root),
        "source_file_count": len(results),
        "resources": results,
        "source_status_unchanged": before.splitlines() == status.splitlines(),
        "source_counter_semantically_unchanged": counter_equal,
        "scope": "Resource byte identity and status/counter comparison, not a hash of every source-worktree file.",
    }
    target = OUT / "source_preservation.json"
    target.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    passed = all(x["matches_frozen_source_and_copy"] for x in results)
    print(json.dumps({
        "resources_unchanged": passed,
        "files": len(results),
        "status_unchanged": record["source_status_unchanged"],
        "counter_unchanged": counter_equal,
        "record": str(target),
    }))
    return 0 if passed and record["source_status_unchanged"] and counter_equal else 1


if __name__ == "__main__":
    raise SystemExit(main())
