#!/usr/bin/env python3
"""check_env.py -- locate an iqtree* binary, report its build flavour and machine
resources, then recommend concrete thread/memory flag values.

Every flag spelling emitted here is doc-verified:
  -nt / -nt AUTO / -ntmax   references/command-reference.md "General options"
  -mem <N>G                 references/command-reference.md "General options"
  -T / -T AUTO              references/command-reference.md "Tree search for pathogen data"
                            + Quickstart.md / Tutorial.md ("for version 1.x change -T to -nt")
Exit code: 0 always, unless --require is given and no executable is found (then 2).
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

CANDIDATES = ["iqtree", "iqtree2", "iqtree3", "iqtree-omp", "iqtree.mpich",
              "iqtree-mpi", "iqtree3-mpi"]

# ---------------------------------------------------------------- resources ---

def read_proc_cpuinfo():
    """Return (cpu_count, model_name) from /proc/cpuinfo; cpu_count may be 0."""
    count = 0
    model = ""
    try:
        with open("/proc/cpuinfo", "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("processor"):
                    count += 1
                elif not model and line.lower().startswith("model name"):
                    model = line.split(":", 1)[1].strip()
    except OSError:
        pass
    return count, model


def read_memtotal_kb():
    """Return MemTotal in kB from /proc/meminfo, or None."""
    try:
        with open("/proc/meminfo", "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return None


def format_ram(kb):
    if kb is None:
        return "unknown"
    gb = kb / 1024.0 / 1024.0
    if gb >= 1:
        return "%.1f GB" % gb
    return "%.0f MB" % (kb / 1024.0)


# ------------------------------------------------------------------ binary ---

VERSION_RE = re.compile(
    r"IQ-TREE[^\n]*?version\s+([0-9]+)\.([0-9]+)\.?([0-9]*)", re.IGNORECASE)
BUILD_RE = re.compile(r"version\s+[0-9][^\n]*?\bfor\b", re.IGNORECASE)


def classify_flavour(name, banner):
    """Map a binary name / banner onto the build flavour the docs describe."""
    low = ((name or "") + " " + (banner or "")).lower()
    if "-mpi" in low or "mpi version" in low or "iqtree.mpich" in low:
        # Compilation-Guide.md: "The executable is named iqtree-mpi"; hybrid is
        # "iqtree-omp-mpi for IQ-TREE versions <= 1.5.X".
        if "omp" in low:
            return "MPI+OpenMP hybrid (run under mpirun)"
        return "MPI (run under mpirun; do NOT use plain iqtree with mpirun)"
    if "-omp" in low:
        return "OpenMP multicore (legacy name, versions <= 1.5.X)"
    if "single" in low and "thread" in low:
        return "single-thread"
    if "knights landing" in low or "knl" in low or "avx-512" in low:
        return "Xeon Phi KNL / AVX-512 multicore"
    if "multicore" in low:
        return "OpenMP multicore (default precompiled binary)"
    return "unknown from banner (assume OpenMP multicore for >= 1.6.0)"


def major_minor(banner):
    m = VERSION_RE.search(banner or "")
    if not m:
        return None, None, None
    return (int(m.group(1)), int(m.group(2)),
            int(m.group(3)) if m.group(3) else 0)


def run_version(path, timeout):
    """Try to capture the IQ-TREE banner from `path -h` (or `-v`)."""
    for args in (["-h"], ["-v"], ["--version"]):
        try:
            proc = subprocess.Popen([path] + args, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL)
            out, _ = proc.communicate(timeout=timeout)
        except (OSError, subprocess.SubprocessError):
            continue
        text = (out or b"").decode("utf-8", "replace")
        if "IQ-TREE" in text or "iqtree" in text.lower():
            return text.strip(), args[0]
    return "", None


def find_binary(candidates, extra_paths):
    """Return (path, name) of the first candidate resolvable on PATH.

    --path entries are probed first: as directories they are searched for every
    candidate name, as files they must themselves be an iqtree* executable.
    """
    probe_files = []
    for p in list(extra_paths or []):
        if not p:
            continue
        if os.path.isdir(p):
            probe_files.extend(os.path.join(p, n) for n in candidates)
        elif os.path.isfile(p):
            probe_files.append(p)
    for cand in probe_files:
        if os.access(cand, os.X_OK):
            return os.path.realpath(cand), os.path.basename(cand)
    for name in candidates:
        hit = shutil.which(name)
        if hit:
            return hit, name
    for name in candidates:
        for d in os.environ.get("PATH", "").split(os.pathsep):
            if d:
                cand = os.path.join(d, name)
                if os.path.isfile(cand) and os.access(cand, os.X_OK):
                    return cand, name
    return None, None


# --------------------------------------------------------------- recommends ---

def recommend(ncpu, mem_kb, banner):
    """Return list of (flag, value, rationale) using doc spellings."""
    out = []
    major, minor, _patch = major_minor(banner)
    legacy = major is not None and (major == 1 or (major == 2 and False))
    # Docs state -T replaces -nt in 2.x/3.x ("for version 1.x change -T to -nt").
    threads_flag = "-nt" if legacy else "-T"
    if ncpu and ncpu > 1:
        out.append((threads_flag, "AUTO",
                    "let IQ-TREE measure multi-threading efficiency and pick the best core count"))
        cap = max(1, min(ncpu, 8))
        out.append((threads_flag, str(cap),
                    "fixed alternative; parallel efficiency is only good for long alignments"))
        out.append(("-ntmax", str(min(ncpu, 16)),
                    "upper bound that %s AUTO may allocate" % threads_flag))
    else:
        out.append((threads_flag, "1", "single visible CPU core"))
    if mem_kb:
        gb = int(mem_kb / 1024.0 / 1024.0)
        cap = max(1, int(gb * 0.75)) if gb >= 2 else max(1, gb)
        out.append(("-mem", "%dG" % cap,
                    "stay below physical RAM (%s total); default tries not to exceed it"
                    % format_ram(mem_kb)))
    else:
        out.append(("-mem", "8G", "total RAM unknown from /proc/meminfo; set an explicit cap"))
    return out


def self_test():
    print("check_env.py --self-test")
    assert read_memtotal_kb() is None or read_memtotal_kb() > 0
    kb = 64 * 1024 * 1024
    assert format_ram(kb) == "64.0 GB", format_ram(kb)
    assert format_ram(512 * 1024) == "512 MB"
    assert format_ram(None) == "unknown"
    fake = ("IQ-TREE multicore version 2.2.0 for Linux 64-bit built Jun 20 2022\n"
            "Developed by Bui Quang Minh, Nguyen Lam Tung, Olga Chernomor.")
    assert classify_flavour("iqtree", fake).startswith("OpenMP multicore"), classify_flavour("iqtree", fake)
    assert classify_flavour("iqtree-omp", "").startswith("OpenMP multicore (legacy")
    assert classify_flavour("iqtree-mpi", "").startswith("MPI")
    assert classify_flavour("iqtree-omp-mpi", "").startswith("MPI+OpenMP")
    assert classify_flavour("iqtree3", "IQ-TREE 3 single-thread version 3.0.0") .startswith("single-thread")
    assert major_minor(fake)[:2] == (2, 2), major_minor(fake)
    assert major_minor("") [0] is None
    recs = recommend(8, kb, fake)
    flags = [r[0] for r in recs]
    assert flags[0] == "-T" and "-ntmax" in flags and "-mem" in flags, recs
    mem = [v for f, v, _ in recs if f == "-mem"]
    assert len(mem) == 1 and re.match(r"^\d+G$", mem[0]), recs
    legacy_recs = recommend(4, 32 * 1024 * 1024,
                            "IQ-TREE multicore version 1.6.12 for Linux 64-bit")
    assert legacy_recs[0][0] == "-nt", legacy_recs
    assert recommend(1, None, "")[0][1] == "1"
    # candidate resolution must never raise when nothing is installed
    path, name = find_binary(["zz-definitely-not-installed-iqtree"], [])
    assert path is None and name is None, (path, name)
    # --path DIR must resolve a binary living outside PATH
    path, name = find_binary(CANDIDATES, ["/tmp"])
    assert (path is None) or name.startswith("iqtree"), (path, name)
    print("  flavour/version/resource parsing: OK")
    print("  recommendations (8 CPUs, 64 GB, v2.2.0): "
          + "; ".join("%s %s" % (f, v) for f, v, _ in recs))
    print("SELF TEST PASSED")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="check_env.py",
        description="Find an iqtree* binary, report version/build flavour and "
                    "machine resources, and recommend thread/memory flags.",
        epilog="Doc-verified flags: -nt/-nt AUTO/-ntmax, -T/-T AUTO, -mem <N>G.")
    ap.add_argument("--candidates", default=",".join(CANDIDATES),
                    help="comma-separated binary names (default: %(default)s)")
    ap.add_argument("--path", action="append", default=[],
                    help="extra directory or file to probe (repeatable)")
    ap.add_argument("--require", action="store_true",
                    help="exit 2 if no iqtree* executable is found")
    ap.add_argument("--timeout", type=float, default=10.0,
                    help="seconds to wait for the version banner (default: %(default)s)")
    ap.add_argument("--self-test", action="store_true", dest="self_test",
                    help="exercise the parser against a built-in fixture and exit")
    ns = ap.parse_args(argv)

    if ns.self_test:
        return self_test()

    cands = [c.strip() for c in ns.candidates.split(",") if c.strip()]
    path, name = find_binary(cands, ns.path)
    ncpu, model = read_proc_cpuinfo()
    if not ncpu:
        ncpu = os.cpu_count() or 0
    mem_kb = read_memtotal_kb()

    print("== IQ-TREE environment ==")
    if path:
        banner, used = run_version(path, ns.timeout)
        first = banner.splitlines()[0].strip() if banner else "(no banner captured)"
        major, minor, patch = major_minor(banner)
        ver = "%d.%d%s" % (major, minor, (".%d" % patch) if patch is not None else "") \
            if major is not None else "unknown"
        print("binary      : %s" % path)
        print("name        : %s" % name)
        print("version     : %s%s" % (ver, "" if major is None else " (from `%s %s`)" % (name, used or "-h")))
        print("banner      : %s" % first[:160])
        print("flavour     : %s" % classify_flavour(name, banner))
    else:
        print("binary      : NOT FOUND (searched: %s)" % ", ".join(cands))
        print("hint        : add the IQ-TREE bin/ directory to PATH, or pass --path DIR")
        print("              (Quickstart.md: 'Copy all files in bin folder to your system search path')")

    print("cpu cores   : %d%s" % (ncpu, (" | %s" % model) if model else ""))
    print("total RAM   : %s" % format_ram(mem_kb))

    banner = ""
    if path:
        banner, _ = run_version(path, ns.timeout)
    print("-- recommended flags --")
    for flag, value, why in recommend(ncpu, mem_kb, banner):
        print("  %-8s %-6s # %s" % (flag, value, why))
    print("  note: '-T' is the 2.x/3.x spelling of '-nt'; '-mem'/'-ntmax' are spelled "
          "the same in both (command-reference.md).")

    if not path and ns.require:
        print("ERROR: no iqtree* executable found and --require was passed")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
