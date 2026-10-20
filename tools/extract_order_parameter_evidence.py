"""Extract readable views of frozen REA evidence; no analysis or target execution."""

import hashlib
import json
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "evidence/native/rea-linux-iso-recovered"
OUT = ROOT / "evidence/native/order-parameter-audit"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    binary = (ROOT / "vendor/isobyu/iso").read_bytes()
    binary_sha = sha(binary)
    index = {
        "purpose": "Offline readable extraction of already retained REA results.",
        "target_sha256": binary_sha,
        "records": [],
    }
    for stem in ("05-orderparam-check", "06-orderparam-to-subgroup", "09-ncmp"):
        source_path = SOURCE / (stem + ".json")
        raw = source_path.read_bytes()
        packet = json.loads(raw)
        result = packet["response"]["result"]
        if result.get("isError", False):
            raise ValueError("Cannot extract an unsuccessful analysis: " + stem)
        evidence = result["structuredContent"]
        if evidence["subject"]["digest"]["sha256"] != binary_sha:
            raise ValueError("REA evidence belongs to different target bytes.")
        value = evidence["raw_result"]["value"]
        pseudo = value if isinstance(value, str) else value["pseudocode"]
        output = OUT / (stem + ".pseudo.txt")
        encoded = pseudo.encode("utf-8")
        output.write_bytes(encoded)
        index["records"].append({
            "source_packet": source_path.relative_to(ROOT).as_posix(),
            "source_sha256": sha(raw),
            "evidence_id": evidence["evidence_id"],
            "provider": evidence["provider"],
            "analysis_profile_digest": evidence["analysis_profile"]["digest"],
            "json_pointer": "/response/result/structuredContent/raw_result/value" +
                            ("" if isinstance(value, str) else "/pseudocode"),
            "extracted_path": output.relative_to(ROOT).as_posix(),
            "extracted_sha256": sha(encoded),
            "extraction_changes": "none; decoded JSON string encoded as UTF-8",
        })
    # Map the known read-only call-argument constants through the ELF PT_LOAD map.
    if binary[:6] != b"\x7fELF\x02\x01":
        raise ValueError("Expected ELF64 little-endian target.")
    phoff = struct.unpack_from("<Q", binary, 32)[0]
    phentsize, phnum = struct.unpack_from("<HH", binary, 54)
    virtual_address = 0x5AF834
    for offset in range(phoff, phoff + phentsize * phnum, phentsize):
        kind, flags, file_offset, vaddr, _, filesz, _, _ = struct.unpack_from(
            "<IIQQQQQQ", binary, offset
        )
        if kind == 1 and vaddr <= virtual_address and virtual_address + 16 <= vaddr + filesz:
            mapped = file_offset + virtual_address - vaddr
            data = binary[mapped:mapped + 16]
            constants = {
                "source": "vendor/isobyu/iso",
                "source_sha256": binary_sha,
                "virtual_address": hex(virtual_address),
                "file_offset": hex(mapped),
                "pt_load_flags": flags,
                "hex": data.hex(),
                "interpretation": "four little-endian signed int32 values used by call sites",
                "int32_values": list(struct.unpack("<4i", data)),
                "scope": "Call-site constants only; does not decode common_little database fields.",
            }
            (OUT / "call-argument-constants.json").write_text(
                json.dumps(constants, indent=2) + "\n", encoding="utf-8"
            )
            break
    else:
        raise ValueError("Call-argument constants are not in a mapped file-backed segment.")
    (OUT / "index.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"target_sha256": binary_sha, "records": len(index["records"]),
                      "output": str(OUT)}))


if __name__ == "__main__":
    main()
