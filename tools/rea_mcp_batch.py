"""Run a real, isolated stdio REA MCP session and retain complete responses.

No client registration changes. Supply an explicit public REA launcher command.
This is analysis tooling and is never imported by the calculation application.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import queue
import subprocess
import threading
import time

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("an explicit REA MCP command is required")
    args.output.mkdir(parents=True, exist_ok=False)
    plan = json.loads(args.plan.read_text(encoding="utf-8-sig"))
    events = queue.Queue()
    stderr_file = (args.output / "server.stderr.txt").open("wb")
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=stderr_file, text=True, encoding="utf-8",
                            cwd=Path(__file__).resolve().parents[1])
    def reader():
        for line in proc.stdout:
            try:
                events.put(json.loads(line))
            except json.JSONDecodeError:
                events.put({"non_json_stdout": line})
        events.put({"eof": True})
    threading.Thread(target=reader, daemon=True).start()
    sequence = 0
    journal = []
    def send(obj):
        proc.stdin.write(json.dumps(obj) + "\n")
        proc.stdin.flush()
    def request(method, params, label, timeout=600):
        nonlocal sequence
        sequence += 1
        started = time.monotonic()
        req = dict(jsonrpc="2.0",id=sequence,method=method,params=params)
        send(req)
        notifications = []
        while True:
            event = events.get(timeout=max(0.1, timeout-(time.monotonic()-started)))
            if event.get("eof"):
                raise RuntimeError(f"REA exited unexpectedly: {proc.poll()}")
            if event.get("id") == sequence:
                break
            notifications.append(event)
        entry = dict(label=label,at_utc=datetime.now(timezone.utc).isoformat(),
                     elapsed_seconds=time.monotonic()-started,request=req,
                     response=event,notifications=notifications)
        file = args.output / f"{sequence:02d}-{label}.json"
        file.write_text(json.dumps(entry,ensure_ascii=False,indent=2),encoding="utf-8")
        result = event.get("result", {})
        journal.append(dict(label=label,path=file.name,elapsed_seconds=entry["elapsed_seconds"],
                            error=event.get("error"),isError=result.get("isError",False)))
        (args.output/"index.json").write_text(json.dumps(journal,indent=2),encoding="utf-8")
        print(json.dumps(journal[-1]),flush=True)
        return result
    try:
        request("initialize",dict(protocolVersion="2024-11-05",capabilities={},
                clientInfo=dict(name="cris-rea-evidence-client",version="0.1")),"initialize",60)
        send(dict(jsonrpc="2.0",method="notifications/initialized"))
        inventory = request("tools/list",{},"tools-list",60)
        advertised = {t["name"] for t in inventory.get("tools",[])}
        for step in plan:
            if step["name"] not in advertised:
                raise RuntimeError("Tool not advertised: " + step["name"])
            request("tools/call",dict(name=step["name"],arguments=step.get("arguments",{})),
                    step["label"],step.get("timeout",600))
    finally:
        proc.stdin.close()
        try:
            proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            proc.terminate()
            proc.wait(timeout=15)
        stderr_file.close()

if __name__ == "__main__":
    main()
