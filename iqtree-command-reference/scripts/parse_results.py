#!/usr/bin/env python3
"""parse_results.py -- summarise the output files of an IQ-TREE run.

Reads PREFIX.iqtree, PREFIX.log, PREFIX.treefile (plus a handful of optional
siblings) and prints a compact, LLM-friendly digest: best-fit model with its
BIC/AICc, log-likelihood, tree length, taxa/sites, topology-test table, run
time, thread/RAM usage and the trailing NOTE/WARNING/ERROR messages.

Missing inputs are reported as ABSENT and still exit 0 (a crashed or partial
run is exactly when this is most useful).  Exit 1 only if PREFIX matches no
output file at all.

Section names parsed here appear verbatim in the documentation mirror:
  MAXIMUM LIKELIHOOD TREE / USER TREES (Advanced-Tutorial.md, Tutorial.md)
  TOTAL CPU TIME / TOTAL WALL-CLOCK TIME / FINAL LOG-LIKELIHOOD VALUE
  BEST-FIT MODEL / BEST-FIT PARTITION SCHEME / LIKELIHOOD MAPPING STATISTICS
  (Command-Reference.md, Tutorial.md, Assessing-Phylogenetic-Assumptions.md)
"""

import argparse
import glob
import os
import re
import sys

# prefix -> (suffix, doc source that mentions it, one-line purpose)
OUTPUT_FILES = [
    (".iqtree",   "main report: model, tree, supports, tests (Tutorial.md)"),
    (".log",      "full run log; send first/last 10 lines with bug reports (FAQ)"),
    (".treefile", "ML tree in NEWICK (Tutorial.md)"),
    (".ckp.gz",   "checkpoint for resuming an interrupted run (Command-Reference.md)"),
    (".model.gz", "ModelFinder checkpoint, >=1.6.X (Command-Reference.md)"),
    (".model",    "ModelFinder checkpoint, <=1.5.X (Command-Reference.md)"),
    (".contree",  "consensus tree with branch supports (Tutorial.md)"),
    (".splits.nex", "split/bipartition support values for SplitsTree (Tutorial.md)"),
    (".splits",   "same as .splits.nex in star-dot format, needs -wsplits (Tutorial.md)"),
    (".ufboot",   "bootstrap trees, written with -wbt/-wbtl (Command-Reference.md)"),
    (".state",    "ancestral state reconstruction per node/site, -asr (Command-Reference.md)"),
    (".rate",     "per-site rates, --rate / -wsr (Command-Reference.md)"),
    (".mlrate",   "ML site rates, --mlrate, IQ-TREE 3 (Advanced-Tutorial.md)"),
    (".sitelh",   "site log-likelihoods, -wsl* (Command-Reference.md)"),
    (".siteprob", "site posterior probabilities, -wsp* (Command-Reference.md)"),
    (".sitefreq", "PMSF site frequency profiles from -ft (Command-Reference.md)"),
    (".partlh",   "partition log-likelihoods, -wpl (Command-Reference.md)"),
    (".alninfo",  "alignment site statistics, -alninfo (Command-Reference.md)"),
    (".lmap.svg", "likelihood-mapping plot, -lmap (Command-Reference.md)"),
    (".lmap.eps", "likelihood-mapping plot, -lmap (Command-Reference.md)"),
    (".varsites", "variable-sites alignment written when +ASC needs it (FAQ)"),
    (".uniqueseq.phy", "deduplicated alignment; may keep 2 identical seqs (FAQ)"),
    (".boottrees", "standard bootstrap trees, combine with cat (FAQ)"),
    (".suptree",  "target tree with mapped support values, -sup (Command-Reference.md)"),
    (".conaln",   "concatenated alignment for multi-file NEXUS partitions (Advanced-Tutorial.md)"),
    (".treels",   "all locally optimal trees, -wt (Command-Reference.md)"),
]

ABSENT_OK = {".model", ".model.gz", ".ckp.gz", ".splits", ".sitefreq", ".mlrate"}


def read_text(path):
    try:
        with open(path, "r", errors="replace") as fh:
            return fh.read()
    except OSError:
        return None


NUM = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"


def grab(text, pattern, cast=None):
    """Return group 1 of `pattern`, or the whole match when the pattern has no group."""
    m = re.search(pattern, text or "", re.MULTILINE)
    if not m:
        return None
    val = m.group(1) if m.groups() else m.group(0)
    if cast:
        try:
            return cast(val)
        except ValueError:
            return val
    return val


DASH_RULE = r"^[ \t]*[-=#]{4,}[ \t]*$"


def section(text, name):
    """Return the body of an ALL-CAPS heading block up to the next real heading.

    A "real heading" is an all-caps line immediately followed by a dashed rule,
    which is how IQ-TREE formats its report sections (Tutorial.md,
    Advanced-Tutorial.md); legend lines such as "deltaL : logL difference ..."
    must not terminate the block.
    """
    pat = r"^%s[^\n]*\n%s\n(.*?)(?=^(?:[A-Z][A-Z0-9 ,()/._-]{3,}[^\n]*\n%s)|\Z)" \
        % (re.escape(name), DASH_RULE, DASH_RULE)
    m = re.search(pat, text or "", re.M | re.S)
    if not m:
        return None
    body = m.group(1)
    # trailing blank line before the next heading is noise
    return body.rstrip("\n") + "\n"


def parse_iqtree(text):
    """Extract the interesting bits from a .iqtree report."""
    d = {}
    # real IQ-TREE writes "BEST-FIT MODEL BY BIC: TIM2+I+G4" and
    # "BEST-FIT MODEL: ..."; tolerate both plus the "(BY x)" variant.
    m = re.search(r"BEST-FIT MODEL(?:\s+BY\s+([A-Za-z]+)|\s*\(BY\s+([A-Za-z]+)\))?:"
                  r"\s*(\S+)", text or "")
    if m:
        d["best_model"] = m.group(3)
        d["best_criterion"] = m.group(1) or m.group(2)
    pm = re.search(r"BEST-FIT PARTITION SCHEME:\s*(\S+)", text or "")
    if pm:
        d["best_scheme"] = pm.group(1)
    d["loglik"] = grab(text, r"FINAL LOG-LIKELIHOOD VALUE:\s*(" + NUM + r")", float)
    d["tree_length"] = grab(text, r"Sum of branch lengths:\s*\(\s*" + NUM + r"\)\s*(" + NUM + ")", float) \
        or grab(text, r"Total tree length[^\n]*?\):\s*(" + NUM + ")", float)
    d["ntaxa"] = grab(text, r"Number of taxa:\s*(\d+)", int)
    d["nsites"] = grab(text, r"(?:Number of (?:constant )?sites|Alignment length):\s*(\d+)", int)
    d["npatterns"] = grab(text, r"Number of (?:distinct )?patterns:\s*(\d+)", int)
    d["cpu_time"] = grab(text, r"TOTAL CPU TIME:\s*(.+)")
    d["wall_time"] = grab(text, r"TOTAL WALL-CLOCK TIME:\s*(.+)")
    d["threads"] = grab(text, r"Best-scoring model using (\d+) CPU", int) \
        or grab(text, r"BEST NUMBER OF THREADS:\s*(\d+)", int) \
        or grab(text, r"Using\s+(\d+)\s+CPU", int)
    d["ram"] = grab(text, r"(?:Maximum RAM|RAM)\s*(?:usage)?:?\s*(\d[\d,.]*\s*[KMGT]?B)", None)
    d["tree_note"] = grab(text, r"^(NOTE: Tree is UNROOTED[^\n]*)", None) \
        or grab(text, r"^(NOTE: Tree is ROOTED[^\n]*)", None)
    d["support_order"] = grab(text, r"(Numbers in parentheses are[^\n]*)")

    crit = {}
    # ModelFinder ranking table, when the report contains one
    table = modelfinder_table(text)
    if table:
        d["mf_table"] = table
        top = table[0]
        for k in ("AIC", "AICc", "BIC", "LogL", "NumPr"):
            if k in top:
                crit[k] = top[k]
    d["criteria"] = crit
    # keep the raw numeric string: full precision matters more than a float cast
    for label, pat in (("AIC", r"\bAIC:\s*(" + NUM + ")"),
                       ("AICc", r"\bAICc:\s*(" + NUM + ")"),
                       ("BIC", r"\bBIC(?:_SCORE)?:\s*(" + NUM + ")")):
        if label not in d["criteria"]:
            v = grab(text, pat)
            if v is not None:
                d["criteria"][label] = v
    au = grab(text, r"USER TREES", None)
    d["user_trees"] = topology_table(text)
    lm = section(text, "LIKELIHOOD MAPPING STATISTICS")
    if lm:
        d["lmap"] = re.sub(r"\s+", " ", lm.strip())[:300]
    return d


def _split_header(line):
    """Split a multi-space-aligned report header into column names."""
    return [h for h in re.split(r"\s{2,}", line.strip()) if h]


def modelfinder_table(text):
    """Parse the ModelFinder ranking table (rows starting with a rank number)."""
    lines = (text or "").splitlines()
    hdr_idx = None
    for i, l in enumerate(lines[:600]):
        if re.search(r"Best-?model", l) and re.search(r"\bAIC\b", l) and re.search(r"\bBIC\b", l):
            hdr_idx = i
            break
    if hdr_idx is None:
        return []
    header = _split_header(lines[hdr_idx])
    rows = []
    for l in lines[hdr_idx + 1:]:
        if set(l.strip()) <= {"-", "=", " "} and l.strip():
            continue
        if not l.strip():
            if rows:
                break
            continue
        cells = l.split()
        if not cells or not re.match(r"^\d+$", cells[0]):
            break
        rows.append(dict(zip(header, cells)))
        if len(rows) >= 5:
            break
    return rows


def topology_table(text):
    """Parse the USER TREES block: logL, deltaL, bp-RELL, p-KH, p-SH, c-ELW, p-AU."""
    body = section(text, "USER TREES")
    if not body:
        return None
    lines = body.splitlines()
    hdr = None
    rows = []
    for l in lines:
        if re.match(r"^\s*(Tree|Rank)\s+logL\s+deltaL\s", l) or \
                (re.match(r"^\s*logL\s+deltaL\s+\S+", l) and ":" not in l):
            hdr = _split_header(l)
            continue
        if hdr and re.match(r"^\s*\d+\s+[-+]?" + NUM, l):
            cells = re.findall(NUM + r"\s*[-+]?", l)
            rows.append([c.strip() for c in cells])
    legend = "\n".join(l for l in lines if re.match(
        r"^(Plus signs|Minus signs|All tests|deltaL|bp-RELL|p-AU|p-KH|p-SH|c-ELW)",
        l.strip()))
    if not rows:
        return None
    return {"header": hdr, "rows": rows, "legend": legend}


def parse_log(text):
    d = {}
    d["version_line"] = (text or "").splitlines()[0].strip() if text else None
    for line in (text or "").splitlines():
        m = re.search(r"IQ-TREE .*version\s+(\S+)", line)
        if m:
            d["version"] = m.group(1)
            break
    d["threads_line"] = grab(text, r"(THREADS: .+|Best number of threads: .+)")
    d["auto_threads"] = grab(text, r"BEST NUMBER OF THREADS:\s*(\d+)", int)
    d["max_ram"] = grab(text, r"Maximum RAM used:\s*([\d.,]+\s*\w+)", None)
    d["finished"] = bool(re.search(r"Total Wall-Clock Time| iqtree completed|CONCLUSION",
                                   text or "", re.I))
    msgs = []
    for line in (text or "").splitlines():
        s = line.strip()
        if re.match(r"^(WARNING|ERROR|NOTE|STOP RULE|CONCLUSION)[: ]", s):
            msgs.append(s)
    d["messages"] = msgs
    d["stop_rules"] = [m for m in msgs if m.startswith("STOP RULE")]
    d["tail"] = [l.rstrip() for l in (text or "").splitlines() if l.strip()][-6:]
    return d


def count_newick_taxa(text):
    # Tokenise instead of regex-matching delimiters: a regex that consumes the
    # delimiter after each label breaks the next label's own leading match
    # (e.g. "Cow,Whale" lost "Whale").
    tokens = re.split(r"[(),;:\[\]]|\s", text or "")
    return len(set(t for t in tokens if t and not re.match(r"^\d", t)))


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="parse_results.py",
        description="Summarise IQ-TREE output files for a given run prefix.")
    ap.add_argument("--prefix", "-p", help="output prefix, e.g. myrun for myrun.iqtree")
    ap.add_argument("--list-files", dest="list_files", action="store_true",
                    help="only inventory which output files exist")
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="include extra tables (ModelFinder top-5 rows)")
    ap.add_argument("--self-test", action="store_true", dest="self_test",
                    help="run against a built-in fixture set in a temp dir")
    ns = ap.parse_args(argv)
    if ns.self_test:
        return self_test()
    if not ns.prefix:
        ap.error("--PREFIX is required (or use --self-test)")
    return run(ns.prefix, ns.list_files, ns.verbose)


def run(prefix, list_files=False, verbose=False):
    found = [(suf, desc) for suf, desc in OUTPUT_FILES
             if os.path.exists(prefix + suf)]
    print("== IQ-TREE results: %s ==" % prefix)
    if not found:
        cand = glob.glob(os.path.dirname(prefix + ".*") + "/" +
                         os.path.basename(prefix) + "*")
        print("NO OUTPUT FILES FOUND for prefix '%s'" % prefix)
        if cand:
            print("did you mean one of:")
            for c in sorted(cand)[:10]:
                print("  %s" % c)
        else:
            print("(no files starting with '%s*' in this directory either)"
                  % os.path.basename(prefix))
        print("hint: the prefix defaults to the alignment or partition file name "
              "(-pre/--prefix overrides it, command-reference.md)")
        return 1
    present = {s for s, _ in found}
    print("-- file inventory --")
    for suf, desc in OUTPUT_FILES:
        if suf in present:
            size = os.path.getsize(prefix + suf)
            print("  PRESENT %-16s %8d B  # %s" % (suf, size, desc))
    missing = [suf for suf, _ in OUTPUT_FILES if suf not in present]
    notable = [m for m in missing if m not in ABSENT_OK]
    if notable:
        print("  absent  : %s" % ", ".join(notable))
    if list_files:
        return 0

    iq = read_text(prefix + ".iqtree")
    lg = read_text(prefix + ".log")
    tf = read_text(prefix + ".treefile")

    print("-- summary --")
    if iq is None:
        print("  .iqtree : ABSENT (cannot report model/logL/times)")
    else:
        d = parse_iqtree(iq)
        def show(label, val, fmt="%s"):
            if val is not None and val != "" and val != {}:
                print("  %-11s %s" % (label + ":", fmt % val))
        show("model", d.get("best_model") + (" (by %s)" % d["best_criterion"]
                                            if d.get("best_criterion") else "")
             if d.get("best_model") else None)
        if d.get("best_scheme"):
            show("scheme", d["best_scheme"])
        show("log-lik", d.get("loglik"), "%.3f")
        crit = d.get("criteria") or {}
        if crit:
            print("  %-11s %s" % ("criteria:", "  ".join(
                "%s=%s" % (k, crit[k]) for k in ("AIC", "AICc", "BIC", "LogL", "NumPr")
                if k in crit)))
        elif any(k in iq for k in ("AIC", "BIC")):
            print("  criteria  : (no parsimony-criterion table found in .iqtree)")
        show("tree length", d.get("tree_length"), "%.6f")
        if d.get("ntaxa") or d.get("nsites"):
            show("taxa/sites", "%s / %s" % (d.get("ntaxa", "?"), d.get("nsites", "?")))
        if d.get("npatterns"):
            show("patterns", d["npatterns"])
        if d.get("threads"):
            show("threads", d["threads"])
        if d.get("cpu_time") or d.get("wall_time"):
            show("time", "CPU %s | wall-clock %s" % (d.get("cpu_time", "?"),
                                                     d.get("wall_time", "?")))
        if d.get("tree_note"):
            show("rooting", d["tree_note"].replace("NOTE: ", ""))
        if d.get("support_order"):
            show("supports", d["support_order"])
        if d.get("lmap"):
            show("lmap", d["lmap"][:200])
        ut = d.get("user_trees")
        if ut:
            print("  TOPOLOGY TESTS (USER TREES):")
            print("    %s" % " ".join(str(h) for h in ut["header"]))
            for r in ut["rows"]:
                print("    %s" % " ".join(str(c) for c in r))
            if ut["legend"]:
                keys = [l for l in ut["legend"].splitlines()
                        if l.startswith(("Plus signs", "Minus signs", "All tests"))]
                for l in keys:
                    print("    %s" % l.strip())
        if verbose and d.get("mf_table"):
            print("  MODELFINDER TOP ROWS:")
            for r in d["mf_table"]:
                print("    " + " ".join("%s=%s" % kv for kv in sorted(r.items())))
    if lg is None:
        print("  .log      : ABSENT")
    else:
        L = parse_log(lg)
        if L.get("version"):
            print("  version   : %s | %s" % (L["version"], L.get("version_line", "")[:70]))
        if L.get("auto_threads") is not None:
            print("  AUTO thr. : %s (from 'BEST NUMBER OF THREADS')" % L["auto_threads"])
        if L.get("max_ram"):
            print("  RAM used  : %s" % L["max_ram"])
        print("  finished  : %s" % ("yes" if L.get("finished") else
                                   "NO - run looks incomplete; rerun the SAME command line "
                                   "to resume from .ckp.gz"))
        msgs = L.get("messages") or []
        stop = L.get("stop_rules") or []
        if stop:
            print("  STOP RULES:")
            for m in stop[-6:]:
                print("    %s" % m)
        bad = [m for m in msgs if m.startswith(("WARNING", "ERROR"))]
        notes = [m for m in msgs if m.startswith("NOTE")]
        if bad:
            print("  WARNING/ERROR (%d, last 5):" % len(bad))
            for m in bad[-5:]:
                print("    %s" % m[:200])
        else:
            print("  WARNING/ERROR: none")
        if notes:
            print("  NOTES (last 3):")
            for m in notes[-3:]:
                print("    %s" % m[:200])
    if tf is None:
        print("  .treefile : ABSENT")
    else:
        nleaf = count_newick_taxa(tf)
        nl = len([l for l in tf.splitlines() if l.strip()])
        print("  treefile  : %d Newick line(s), ~%d distinct taxon label(s)" % (nl, nleaf))
        last = [l for l in tf.splitlines() if l.strip()][-1].strip()
        print("  tree head : %s%s" % (last[:150], " ..." if len(last) > 150 else ""))
    return 0


FIXTURE_IQTREE = """
-----------------------------------------------------------------------------------------
IQ-TREE 2.2.0 built 2022-06-01

DATA SET SUMMARY
----------------
Number of taxa:                  7
Number of sites:                 28
Number of distinct patterns:     6

MAXIMUM LIKELIHOOD TREE
-----------------------
NOTE: Tree is UNROOTED although outgroup taxon 'Frog' is drawn at root
Numbers in parentheses are SH-aLRT support (%) / ultrafast bootstrap support (%)

 +-------------Frog
 |
 +--| (100/100)
    +------Turtle

Sum of branch lengths: ( 0.01234567)      0.01234567

BEST-FIT MODEL BY BIC: TIM2+I+G4
AIC: -123.456  AICc: -110.111  BIC: -100.222

USER TREES
----------

See example.trees for trees with branch lengths.

Tree      logL    deltaL  bp-RELL    p-KH     p-SH    c-ELW     p-AU
-------------------------------------------------------------------------
  1   -21152.617   0.000  0.7110 + 0.7400 + 1.0000 + 0.6954 + 0.7939 +
  2   -21156.802   4.185  0.2220 + 0.2600 + 0.5910 + 0.2288 + 0.3079 +
  4   -21339.596 186.980  0.0000 - 0.0000 - 0.0000 - 0.0000 - 0.0000 -

deltaL  : logL difference from the maximal logl in the set.
Plus signs denote the 95% confidence sets.
Minus signs denote significant exclusion.
All tests performed 1000 resamplings using the RELL method.

TOTAL CPU TIME: 0.5 CPU sec, 0.3 wall sec
TOTAL WALL-CLOCK TIME: 0.3 wall sec
FINAL LOG-LIKELIHOOD VALUE: -123.4567
Maximum RAM used: 123.45 MB
"""

FIXTURE_LOG = """IQ-TREE multicore version 2.2.0 for Linux 64-bit built Jun  1 2022
Computing pairwise distances...
Measuring multi-threading efficiency up to 8 CPU cores
Threads: 1 / Time: 8.001 sec / Speedup: 1.000 / Efficiency: 100% / LogL: -22217
BEST NUMBER OF THREADS: 3
NOTE: Checkpointing is enabled
WARNING: Site 12 contains only gaps and is ignored
ERROR: Invalid option -zzz ignored
STOP RULE: Total computational time exceeded 10 seconds
Maximum RAM used: 123.45 MB
Total Wall-Clock Time: 0.3 wall sec
"""

FIXTURE_TREEFILE = "(Frog,(Turtle,(Bird,(Human,(Cow,Whale)))),Mouse);\n"


def self_test():
    import tempfile
    print("parse_results.py --self-test")
    tmp = tempfile.mkdtemp(prefix="iqtree-parse-selftest-")
    pre = os.path.join(tmp, "tiny")
    with open(pre + ".iqtree", "w") as fh:
        fh.write(FIXTURE_IQTREE)
    with open(pre + ".log", "w") as fh:
        fh.write(FIXTURE_LOG)
    with open(pre + ".treefile", "w") as fh:
        fh.write(FIXTURE_TREEFILE)

    d = parse_iqtree(FIXTURE_IQTREE)
    assert d["best_model"] == "TIM2+I+G4", d["best_model"]
    assert abs(d["loglik"] - (-123.4567)) < 1e-6, d["loglik"]
    assert abs(d["tree_length"] - 0.01234567) < 1e-9, d["tree_length"]
    assert (d["ntaxa"], d["nsites"], d["npatterns"]) == (7, 28, 6), d
    assert d["criteria"].get("BIC") == "-100.222", d["criteria"]
    assert d["criteria"].get("AICc") == "-110.111", d["criteria"]
    assert d["cpu_time"].startswith("0.5 CPU sec"), d["cpu_time"]
    ut = d["user_trees"]
    assert ut and len(ut["rows"]) == 3, ut
    assert "p-AU" in ut["header"], ut["header"]
    assert "Minus signs" in ut["legend"], ut["legend"]

    L = parse_log(FIXTURE_LOG)
    assert L["version"] == "2.2.0", L
    assert L["auto_threads"] == 3, L
    assert L["max_ram"] == "123.45 MB", L
    assert len(L["stop_rules"]) == 1, L["stop_rules"]
    assert any(x.startswith("ERROR") for x in L["messages"]), L["messages"]
    assert L["finished"], L
    assert count_newick_taxa(FIXTURE_TREEFILE) == 7, count_newick_taxa(FIXTURE_TREEFILE)

    rc = run(pre)
    assert rc == 0, rc
    print("  full-prefix digest: OK")
    # graceful degradation: only .treefile present
    pre2 = os.path.join(tmp, "partial")
    with open(pre2 + ".treefile", "w") as fh:
        fh.write(FIXTURE_TREEFILE)
    assert run(pre2) == 0, "partial run must exit 0"
    print("  partial run (.treefile only): exits 0 and reports absence")
    # nothing there at all
    assert run(os.path.join(tmp, "nothing")) == 1
    print("  unknown prefix: exits 1 with a hint")
    print("SELF TEST PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
