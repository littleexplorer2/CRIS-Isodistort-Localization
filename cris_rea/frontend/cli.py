import argparse
import json
from pathlib import Path
import sys

from ..application import RunService


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="python -m cris_rea",
        description="Local official ISOTROPY tools with per-run evidence.",
    )
    actions = parser.add_subparsers(dest="action", required=True)
    actions.add_parser("doctor", help="Check local vendor files and Linux/WSL runtime")
    run_parser = actions.add_parser("run", help="Run a trusted raw input file")
    run_parser.add_argument("program", help="iso, findsym, findsym_cifinput, smodes, comsubs")
    run_parser.add_argument("--input", type=Path, required=True)
    run_parser.add_argument("--timeout", type=float, default=30)
    query_parser = actions.add_parser("query", help="Prepare or run a restricted ISO query")
    query_parser.add_argument("--parent", type=int, required=True)
    query_parser.add_argument("--irrep", required=True)
    query_parser.add_argument("--opd", required=True)
    query_parser.add_argument("--wyckoff")
    query_parser.add_argument("--timeout", type=float, default=30)
    query_parser.add_argument("--run", action="store_true", help="Execute the prepared commands")
    serve_parser = actions.add_parser("serve", help="Start the local browser interface")
    serve_parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[2]
    service = RunService(root)
    try:
        if args.action == "doctor":
            result = service.doctor()
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if result["execution_supported"] else 2
        if args.action == "run":
            # Bound file reads before accepting arbitrarily large local files.
            with args.input.open("rb") as stream:
                data = stream.read(1024 * 1024 + 1)
            result = service.run(args.program, data, args.timeout)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if result["evidence"]["status"] == "completed_uninterpreted" else 1
        if args.action == "query":
            query = service.prepare_iso_query(args.parent, args.irrep, args.opd, args.wyckoff)
            if not args.run:
                print(query["input"], end="")
                return 0
            result = service.run("iso", query["input"].encode("ascii"), args.timeout)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if result["evidence"]["status"] == "completed_uninterpreted" else 1
        from .web import serve
        serve(service, args.port)
        return 0
    except (ValueError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130
