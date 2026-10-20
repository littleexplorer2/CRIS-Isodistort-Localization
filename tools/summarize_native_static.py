"""Extract observed direct call sites from retained GNU objdump text."""
import hashlib
import json
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
static = root / "evidence/native/static"
records = []
for path in sorted(static.glob("*.asm.txt")):
    text = path.read_text()
    heads = re.findall(r"^([0-9a-f]+) <([^>]+)>:$", text, re.M)
    calls = [
        dict(site="0x"+site,target_address="0x"+address,target_symbol=symbol)
        for site,address,symbol in re.findall(
            r"^\s*([0-9a-f]+):.*\bcall\s+([0-9a-f]+)\s+<([^>]+)>",text,re.M)
    ]
    records.append(dict(file=path.relative_to(root).as_posix(),
                        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        functions=[dict(address="0x"+a,symbol=s) for a,s in heads],
                        direct_calls=calls,
                        limitation="Static call instructions only; no runtime execution or indirect-call completeness"))
out = root / "evidence/native/direct-calls.json"
out.write_text(json.dumps(records,indent=2),encoding="utf-8")
print(json.dumps(dict(files=len(records),direct_calls=sum(len(r["direct_calls"]) for r in records))))
