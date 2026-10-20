"""Read-only byte inventory of an explicitly selected ISOTROPY data directory.

This is a lexical evidence scanner, NOT a crystallographic database decoder.
Candidate headers match one lowercase identifier on an otherwise empty line.
No pointer origin, array shape, unit, or physical meaning is inferred.
Python standard library only; source files are opened in binary read mode.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys


HEADER = re.compile(rb"[a-z][a-z0-9_]*")
TOKEN = re.compile(rb'"(?:[^"]|"")*"|\'(?:[^\']|\'\')*\'|[^\s,]+')
INTEGER = re.compile(rb"[+-]?\d+")
REAL = re.compile(rb"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def byte_sample(data: bytes, offset: int, length: int = 96) -> dict:
    raw = data[offset : offset + length]
    return {
        "offset": offset,
        "length": len(raw),
        "hex": raw.hex(" "),
        "ascii_escaped": raw.decode("ascii", errors="backslashreplace"),
    }


def lexical_payload(raw: bytes) -> dict:
    counts = Counter()
    samples = []
    int_min = None
    int_max = None
    for match in TOKEN.finditer(raw):
        token = match.group()
        if token.startswith((b'"', b"'")):
            kind = "quoted"
        elif INTEGER.fullmatch(token):
            kind = "integer"
            value = int(token)
            int_min = value if int_min is None else min(int_min, value)
            int_max = value if int_max is None else max(int_max, value)
        elif REAL.fullmatch(token):
            kind = "real_lexeme"
        elif token in (b"T", b"F", b".TRUE.", b".FALSE.", b".true.", b".false."):
            kind = "logical_lexeme"
        else:
            kind = "other"
        counts[kind] += 1
        if len(samples) < 12:
            samples.append(
                {
                    "relative_byte_offset": match.start(),
                    "kind": kind,
                    "raw": token[:96].decode("ascii", errors="backslashreplace"),
                    "sample_truncated": len(token) > 96,
                }
            )
    return {
        "tokenizer": "quoted strings or non-whitespace/comma runs; lexical only",
        "token_counts": dict(sorted(counts.items())),
        "integer_min": int_min,
        "integer_max": int_max,
        "first_tokens": samples,
    }


def scan_bytes(data: bytes) -> dict:
    histogram = Counter(data)
    unusual = {
        f"0x{byte:02x}": count
        for byte, count in sorted(histogram.items())
        if byte not in (9, 10, 13) and not 32 <= byte <= 126
    }
    headers = []
    offset = 0
    for line_number, line in enumerate(data.splitlines(keepends=True), 1):
        if HEADER.fullmatch(line.strip()):
            headers.append(
                {
                    "name": line.strip().decode("ascii"),
                    "line": line_number,
                    "header_offset": offset,
                    "payload_offset": offset + len(line),
                }
            )
        offset += len(line)
    sections = []
    for index, header in enumerate(headers):
        stop = headers[index + 1]["header_offset"] if index + 1 < len(headers) else len(data)
        raw = data[header["payload_offset"] : stop]
        sections.append(
            {
                **header,
                "payload_end_exclusive": stop,
                "payload_bytes": len(raw),
                "payload_sha256": digest(raw),
                "sample": byte_sample(data, header["payload_offset"], min(96, len(raw))),
                **lexical_payload(raw),
            }
        )
    unusual_runs = []
    for match in re.finditer(rb"[^\x09\x0a\x0d\x20-\x7e]+", data):
        if len(unusual_runs) < 32:
            unusual_runs.append(
                {
                    "offset": match.start(),
                    "length": len(match.group()),
                    "hex_prefix": match.group()[:64].hex(" "),
                    "context": byte_sample(data, max(0, match.start() - 48), 144),
                }
            )
    return {
        "size_bytes": len(data),
        "sha256": digest(data),
        "ascii_text_bytes_only": not unusual,
        "unusual_byte_histogram": unusual,
        "nul_count": histogram[0],
        "line_break_counts": {
            "lf": data.count(b"\n"),
            "crlf": data.count(b"\r\n"),
            "bare_cr": data.count(b"\r") - data.count(b"\r\n"),
        },
        "first_96_bytes": byte_sample(data, 0),
        "unusual_runs_first_32": unusual_runs,
        "candidate_header_count": len(headers),
        "candidate_headers": sections,
        "limits": [
            "Headers are lexical candidates, not a decoded schema.",
            "Token counts do not establish array dimensions or pointer semantics.",
            "No scientific calculation or executable invocation is performed.",
        ],
    }


def self_test() -> None:
    fixture = b'alpha_count\r\n2\r\nalpha_label\r\n"A B" "C"\r\nflags\r\nT F\r\nnul_label\r\n"x\x00y"\r\nend\r\n'
    result = scan_bytes(fixture)
    sections = {item["name"]: item for item in result["candidate_headers"]}
    assert result["nul_count"] == 1
    assert not result["ascii_text_bytes_only"]
    assert result["line_break_counts"] == {"lf": 9, "crlf": 9, "bare_cr": 0}
    assert sections["alpha_label"]["sample"]["offset"] == fixture.index(b'"A B"')
    assert sections["alpha_label"]["token_counts"] == {"quoted": 2}
    assert sections["flags"]["token_counts"] == {"logical_lexeme": 2}
    assert sections["end"]["payload_bytes"] == 0
    assert scan_bytes(b"\x7f\xff")["unusual_byte_histogram"] == {"0x7f": 1, "0xff": 1}
    assert scan_bytes(b"")["size_bytes"] == 0
    print("PASS: lexical boundaries, CRLF offsets, NUL, logical rows, empty input")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.source is None or args.output is None:
        parser.error("--source and --output are required unless --self-test is used")
    source = args.source.resolve(strict=True)
    output = args.output.resolve()
    if output == source or source in output.parents:
        parser.error("output must not be inside the source directory")
    files = sorted(source.glob("data_*.txt"))
    if (source / "const.dat").is_file():
        files.insert(0, source / "const.dat")
    if not files:
        parser.error("no data_*.txt or const.dat files found")
    output.mkdir(parents=True, exist_ok=True)
    manifest = {
        "format": "cris-rea-database-byte-inventory-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(source),
        "scanner_sha256": digest(Path(__file__).read_bytes()),
        "scope": "all data_*.txt and const.dat directly inside source",
        "files": [],
    }
    for path in files:
        data = path.read_bytes()
        result = scan_bytes(data)
        result.update({"file": path.name, "source": str(path.resolve())})
        result["source_unchanged_during_scan"] = digest(path.read_bytes()) == result["sha256"]
        if not result["source_unchanged_during_scan"]:
            raise RuntimeError(f"source changed during scan: {path}")
        evidence_name = path.name + ".inventory.json"
        (output / evidence_name).write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        manifest["files"].append(
            {
                key: result[key]
                for key in (
                    "file", "size_bytes", "sha256", "ascii_text_bytes_only", "nul_count",
                    "candidate_header_count", "source_unchanged_during_scan",
                )
            }
            | {"evidence": evidence_name}
        )
    (output / "inventory.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
