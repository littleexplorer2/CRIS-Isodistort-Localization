import argparse
import json
import math
from pathlib import Path
import sys

from ..application import RunService


def _finite_float(value):
    result = float(value)
    if not math.isfinite(result):
        raise argparse.ArgumentTypeError("A finite number is required.")
    return result


def _mode_workflow(root):
    from ..application.modes import ModeWorkflow
    return ModeWorkflow(root)


def _k_values(arguments):
    # Split transport delimiters only; the official protocol validates values and order.
    values = [value.strip() for argument in arguments for value in argument.split(",")]
    if any(not value for value in values):
        raise ValueError("k parameters cannot contain empty comma-separated entries.")
    return values


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
    modes_parser = actions.add_parser(
        "modes", help="Build an ordinary 3D commensurate subgroup-fixed mode basis from a parent CIF",
        description="Local orthonormal subgroup-fixed displacement/strain basis in the supported "
                    "standard setting. This basis is not the website's per-irrep column convention.",
    )
    modes_parser.add_argument("--cif", type=Path, required=True)
    modes_parser.add_argument("--irrep", required=True)
    modes_parser.add_argument("--opd", required=True)
    modes_parser.add_argument("--k-values", nargs="*", default=[],
                              help="ISO independent k parameters in order; spaces or --k-values=1/3,-1/4,0")
    modes_parser.add_argument("--amplitudes", nargs="*", type=_finite_float, default=None,
                              help="Displacement amplitudes in exported basis order (angstrom)")
    modes_parser.add_argument("--strain-amplitudes", nargs="*", type=_finite_float, default=None,
                              help="Dimensionless strain amplitudes in exported basis order")
    modes_parser.add_argument("--timeout", type=_finite_float, default=60,
                              help="Timeout in seconds per official ISO invocation (up to five calls)")
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
        if args.action == "modes":
            k_values = _k_values(args.k_values)
            with args.cif.open("rb") as stream:
                data = stream.read(1024 * 1024 + 1)
            if not data or len(data) > 1024 * 1024:
                raise ValueError("CIF input must be nonempty and at most 1 MiB.")
            result = _mode_workflow(root).run(
                data, args.irrep, args.opd, k_values=k_values,
                amplitudes=args.amplitudes, strain_amplitudes=args.strain_amplitudes,
                timeout=args.timeout,
            )
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if result.get("status") == "validated_local_subspace" else 1
        from .web import serve
        serve(service, args.port)
        return 0
    except (ValueError, OSError, ImportError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130
