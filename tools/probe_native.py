"""Bounded Linux-only experiments; all inputs fixed before any web comparison."""
from pathlib import Path
import hashlib
import json
import os
import platform
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "isobyu"
OUT = ROOT / "evidence" / "native" / "runtime"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(name, program, text, file_argument=False, timeout=30):
    dest = OUT / name
    dest.mkdir(parents=True, exist_ok=False)
    data = dest / "isodata"
    shutil.copytree(VENDOR, data)
    work = dest / "work"
    work.mkdir()
    inp = dest / "input.txt"
    inp.write_text(text, encoding="ascii")
    exe = data / program
    exe.chmod(exe.stat().st_mode | 0o100)
    command = ["../isodata/" + program] + (["../input.txt"] if file_argument else [])
    env = dict(os.environ, ISODATA="../isodata/")
    before = {p.name: digest(p) for p in data.iterdir() if p.is_file()}
    start = time.time()
    timeout_hit = False
    try:
        result = subprocess.run(command, cwd=work, env=env,
                                input=None if file_argument else text.encode("ascii"),
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                timeout=timeout, check=False)
        stdout, stderr, code = result.stdout, result.stderr, result.returncode
    except subprocess.TimeoutExpired as exc:
        stdout, stderr, code = exc.stdout or b"", exc.stderr or b"", None
        timeout_hit = True
    (dest / "stdout.txt").write_bytes(stdout)
    (dest / "stderr.txt").write_bytes(stderr)
    after = {p.name: digest(p) for p in data.iterdir() if p.is_file()}
    metadata = dict(
        name=name, program=program, command=command, cwd=str(work),
        platform=platform.platform(), input_sha256=digest(inp),
        binary_sha256=before[program], isodata_sha256=before,
        elapsed_seconds=time.time()-start, returncode=code,
        timeout=timeout_hit, timeout_seconds=timeout,
        stdout_sha256=digest(dest / "stdout.txt"),
        stderr_sha256=digest(dest / "stderr.txt"),
        changed_data=[k for k in sorted(set(before)|set(after)) if before.get(k)!=after.get(k)],
        generated_work_files=[str(p.relative_to(work)) for p in work.rglob("*") if p.is_file()],
        scientific_status="unreviewed",
        input_source="official command syntax; analyst-selected generic cubic case, no website answer read",
    )
    (dest / "record.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({k: metadata[k] for k in ("name","returncode","timeout","elapsed_seconds")}))
    print(stdout[:12000].decode("utf-8", errors="replace"))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-name", default="runtime")
    parser.add_argument("--repair-probes", action="store_true")
    args = parser.parse_args()
    if not args.output_name.replace("-", "").replace("_", "").isalnum():
        parser.error("output-name must be a simple directory name")
    OUT = ROOT / "evidence" / "native" / args.output_name
    if args.repair_probes:
        run("findsym-short-filename", "findsym", (VENDOR / "findsym_sample.in").read_text(), True)
        run("findsym-cif-short-filename", "findsym_cifinput",
            "data_generic_cubic\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n"
            "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
            "_symmetry_space_group_name_H-M 'P m -3 m'\n_symmetry_Int_Tables_number 221\n"
            "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n"
            "_atom_site_fract_y\n_atom_site_fract_z\nX1 Si 0 0 0\n", True)
        run("comsubs-shipped-input", "comsubs", (VENDOR / "comsubs_sample.in").read_text(), timeout=30)
        raise SystemExit(0)
    run("iso-gm4minus-p1", "iso",
        "PAGE NOBREAK\nVALUE PARENT 221\nVALUE IRREP GM4-\nVALUE DIRECTION P1\n"
        "SHOW PARENT\nSHOW IRREP\nSHOW DIRECTION\nSHOW SUBGROUP\nSHOW BASIS\n"
        "SHOW ORIGIN\nSHOW SIZE\nSHOW INDEX\nDISPLAY ISOTROPY\n"
        "VALUE WYCKOFF a\nSHOW MICROSCOPIC VECTOR\nDISPLAY DISTORTION\nQUIT\n")
    run("iso-invalid-keyword", "iso", "THIS_IS_NOT_A_COMMAND\nQUIT\n")
    run("smodes-shipped-input", "smodes", (VENDOR / "smodes_sample.in").read_text())
    run("findsym-shipped-input", "findsym", (VENDOR / "findsym_sample.in").read_text(), True)
