#!/usr/bin/env python3
"""build_command.py -- assemble a validated IQ-TREE command line.

Every flag in FLAG_MAP below was verified by grepping the documentation mirror;
the `src` key names the file and section it came from.  Two spellings are tracked
for the options that IQ-TREE 2.x/3.x renamed (see the rename table in SKILL.md):
`legacy` is the <=1.x spelling, `modern` is the >=2.x spelling.  The emitted
command uses one style only (--style legacy|modern, default modern).

Nothing here invents flags: anything not represented in FLAG_MAP must be passed
through --extra verbatim.
"""

import argparse
import os
import re
import shlex
import subprocess
import sys

# arg name -> (modern flag, legacy flag, doc source, description)
FLAG_MAP = {
    "alignment":  ("-s", "-s", "command-reference.md 'General options'",
                   "input alignment (PHYLIP/FASTA/NEXUS/CLUSTAL/MSF)"),
    "model":      ("-m", "-m", "command-reference.md 'Automatic model selection' + 'Specifying substitution models'", "substitution model or TEST/MF/MFP/TESTMERGE... keyword"),
    "seqtype":    ("-st", "-st", "command-reference.md 'General options' (-st DNA/AA/BIN/MORPH/CODON/NT2AA)", "DNA | AA | BIN | MORPH | CODON | NT2AA[<code>]"),
    "threads":    ("-T", "-nt", "command-reference.md 'General options' (-nt) / 'Tree search for pathogen data' (-T); Quickstart.md 'for version 1.x change -T to -nt'", "CPU cores, or AUTO"),
    "threads_max": ("-ntmax", "-ntmax", "command-reference.md 'General options' (-ntmax)", "max cores -T/-nt AUTO may allocate"),
    "memory":     ("-mem", "-mem", "command-reference.md 'General options' (-mem 64G)", "RAM cap, e.g. 8G"),
    "prefix":     ("--prefix", "-pre", "command-reference.md 'General options' (-pre); Tutorial.md 'for version 1.x change --prefix to -pre'", "output-file prefix"),
    "seed":       ("-seed", "-seed", "command-reference.md 'General options' (-seed)", "random number seed"),
    "outgroup":   ("-o", "-o", "command-reference.md 'General options' (-o)", "outgroup taxon name(s)"),
    "user_tree":  ("-t", "-t", "command-reference.md 'General options' (-t)", "starting/constraint/user tree file"),
    "fix_tree":   ("-te", "-te", "command-reference.md 'General options' (-te)", "fixed user tree (no tree search)"),
    "bootstrap":  ("-B", "-bb", "command-reference.md 'Ultrafast bootstrap parameters' (-bb); Tutorial.md 'for version 1.x change -B to -bb'", "ultrafast bootstrap replicates (>=1000)"),
    "bootstrap_nni": ("-bnni", "-bnni", "command-reference.md 'Ultrafast bootstrap parameters' (-bnni)", "NNI-optimise UFBoot trees (severe model violation)"),
    "alrt":       ("-alrt", "-alrt", "command-reference.md 'Single branch tests' (-alrt)", "SH-aLRT replicates (>=1000)"),
    "partition":  ("-p", "-spp", "command-reference.md 'Partition model options' (-spp or -p); Advanced-Tutorial.md", "edge-proportional partition file"),
    "partition_equal":  ("-q", "-q", "command-reference.md 'Partition model options' (-q)", "edge-equal partition file"),
    "partition_unlinked": ("-Q", "-sp", "command-reference.md 'Partition model options' (-sp or -Q)", "edge-unlinked partition file"),
    "redo":       ("-redo", "-redo", "command-reference.md 'Checkpointing to resume stopped run' (-redo)", "overwrite outputs and redo from scratch"),
    "constraint": ("-g", "-g", "command-reference.md 'Tree search parameters' (-g)", "topological constraint tree"),
    "site_rates": ("--rate", "-wsr", "command-reference.md 'Rate heterogeneity' (-wsr); Advanced-Tutorial.md 'for version 1.x change --rate to -wsr'", "write per-site rates to .rate"),
    "alninfo":    ("-alninfo", "-alninfo", "command-reference.md 'Miscellaneous options' (-alninfo)", "write alignment site stats"),
    "asr":        ("-asr", "-asr", "command-reference.md 'Ancestral sequence reconstruction' (-asr)", "ancestral state reconstruction -> .state"),
    "stop_iter":  ("-n", "-n", "command-reference.md 'Tree search parameters' (-n)", "number of iterations (0 = skip tree search)"),
    "nstop":      ("-nstop", "-nstop", "command-reference.md 'Tree search parameters' (-nstop, DEFAULT 100)", "unsuccessful iterations before stop (default 100)"),
    "pers":       ("-pers", "-pers", "command-reference.md 'Tree search parameters' (-pers, DEFAULT 0.5)", "NNI perturbation strength 0..1 (default 0.5)"),
}

# model keywords accepted by -m for automatic model selection
MODEL_KEYWORDS = {"TEST", "TESTONLY", "TESTNEW", "TESTNEWONLY", "MF", "MFP",
                  "TESTMERGE", "TESTMERGEONLY", "TESTNEWMERGE", "TESTNEWMERGEONLY",
                  "MF+MERGE", "MFP+MERGE"}
SEQTYPES = {"DNA", "AA", "BIN", "MORPH", "CODON", "NT2AA"}
MEM_RE = re.compile(r"^\d+(\.\d+)?[KMGT]B?$", re.IGNORECASE)
TREE_SUFFIXES = (".treefile", ".nwk", ".nex", ".nexus", ".treels", ".contree",
                 ".suptree", ".tree", ".newick")

# Options that appear in the documentation mirror but are NOT in FLAG_MAP.
# Kept so --extra can be sanity-checked without inventing anything new.
KNOWN_LONG_OPTIONS = {
    "--sampling", "--prefix", "--seed", "--model-joint", "--init-model",
    "--link-exchange", "--init-exchange", "--pathogen", "--pathogen-force",
    "--tree-search", "--sprta", "--sprta-zero-branch", "--sprta-other-places",
    "--disable-local-ref", "--estimate-MAT", "--robust-phy", "--mlrate",
    "--out-aln", "--out-format", "--scf", "--scfl", "--gcf", "--symtest",
    "--symtest-only", "--symtest-remove-bad", "--symtest-remove-good",
    "--symtest-type", "--symtest-pval", "--symtest-keep-zero", "--root-test",
    "--dating", "--date", "--date-tip", "--date-root", "--date-ci",
    "--date-outlier", "--date-options", "--clock-sd", "--alisim",
}
KNOWN_SHORT_OPTIONS = {
    "-h", "-?", "-s", "-st", "-t", "-te", "-o", "-pre", "-nt", "-ntmax", "-seed",
    "-v", "-quiet", "-keep-ident", "-safe", "-mem", "-redo", "-cptime", "-lmap",
    "-lmclust", "-n", "-wql", "-m", "-mwopt", "-a", "-gmedian", "-i", "-wsr",
    "-q", "-spp", "-p", "-sp", "-Q", "-ft", "-fs", "-fmax", "-allnni", "-djc",
    "-fast", "-g", "-ninit", "-ntop", "-nbest", "-nstop", "-pers", "-sprrad",
    "-bb", "-bcor", "-beps", "-bnni", "-bsam", "-nm", "-nstep", "-wbt", "-wbtl",
    "-b", "-bc", "-bo", "-alrt", "-abayes", "-lbp", "-asr", "-asr-min", "-z",
    "-zb", "-zw", "-au", "-con", "-net", "-minsup", "-bi", "-sup", "-suptag",
    "-scale", "-prec", "-rf", "-rf_all", "-rf_adj", "-r", "-ru", "-rcat",
    "-rbal", "-rcsg", "-rlen", "-alninfo", "-blfix", "-blmin", "-blmax", "-czb",
    "-me", "-wpl", "-wspr", "-wspm", "-wspmr", "-wsl", "-wslr", "-wslm",
    "-wslmr", "-wt", "-fconst", "-T", "--alrt", "-mset", "-msub", "-mfreq",
    "-mrate", "-cmin", "-cmax", "-merit", "-mtree", "-mredo", "-madd", "-mdef",
    "-rcluster", "-rclusterf", "-rcluster-max", "-wsplits", "-AIC", "-AICc",
    "-cmax", "--opt-gamma-inv", "-lrt", "-qmax",
}


def sniff(path):
    """Return (kind, detail) for a plain-text input: 'alignment', 'partition' or None.

    kind/detail feed the warnings about datatype and partition syntax; the rules
    come from Tutorial.md (input formats), Complex-Models.md (partition file
    format) and Advanced-Tutorial.md.
    """
    try:
        with open(path, "r", errors="replace") as fh:
            head = [next(fh, "") for _ in range(6)]
    except OSError:
        return None, None
    joined = "".join(head)
    stripped = [l.strip() for l in head if l.strip()]
    low = joined.lower().lstrip()
    if low.startswith("#nexus") or low.startswith("#nex"):
        if re.search(r"\bcharset\b|\bcharpartition\b", joined, re.IGNORECASE):
            return "partition", "nexus"
        return "alignment", "nexus"
    if re.match(r"^\s*\d+\s+\d+\s*$", head[0] or ""):
        return "alignment", "phylip"
    if stripped and all(re.match(r"^(DNA|AA|RNA|Protein|[A-Za-z]+)\s*,\s*[^=]+=\s*[^;]+$", l)
                        for l in stripped[:3]):
        return "partition", "raxml"
    if stripped and stripped[0].startswith(">"):
        return "alignment", "fasta"
    return "alignment", "unknown"


def sequence_length(path):
    """Best-effort site count of an alignment (PHYLIP header, FASTA or NEXUS)."""
    kind, detail = sniff(path)
    try:
        with open(path, "r", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return None
    if detail == "phylip":
        try:
            return int(text.splitlines()[0].split()[1])
        except (IndexError, ValueError):
            return None
    if detail == "fasta":
        seqs, cur = [], []
        for line in text.splitlines():
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur)); cur = []
            elif line.strip() and not line.strip().lower().startswith("#nexus"):
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
        return max((len(s) for s in seqs), default=None)
    if detail == "nexus":
        m = re.search(r"dimensions\s+ntax\s*=\s*(\d+)\s+nchar\s*=\s*(\d+)",
                      text, re.IGNORECASE)
        if m:
            return int(m.group(2))
    return None


RAXML_PART_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9]*)\s*,\s*([^=]+?)\s*=\s*([^;]+?)\s*;?\s*$")


def parse_partition_sites(spec):
    """Parse '1-100 250-384' / '1-100, 250-384' / '101-249\\3' / '*' into ranges.

    Syntax per Complex-Models.md 'Partition file format' and Advanced-Tutorial.md:
    NEXUS separates sub-ranges with whitespace, RAxML-style with commas, and a
    trailing '\\step' extracts every step-th site.  '*' means the whole alignment.
    Return list of (start, end, step, raw) with 1-based inclusive bounds.
    """
    spec = spec.strip()
    if spec in ("*", "ALL", "all"):
        return [(1, None, 1, spec)]
    out = []
    for token in re.split(r"[,\s]+", spec):
        if not token:
            continue
        m = re.match(r"^(\d+)(?:-(\d+))?(?:\\(\d+))?$", token)
        if not m:
            out.append((None, None, None, token))
            continue
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        step = int(m.group(3)) if m.group(3) else 1
        out.append((start, end, step, token))
    return out


def check_partition_ranges(path, n_sites):
    """Return error strings for ranges that fall outside an n_sites-wide alignment."""
    errs = []
    try:
        with open(path, "r", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return errs
    specs = []
    kind, detail = sniff(path)
    if detail == "raxml":
        for line in text.splitlines():
            m = RAXML_PART_RE.match(line)
            if m:
                specs.append((m.group(2).strip(), m.group(3)))
    elif detail == "nexus":
        for m in re.finditer(r"charset\s+(\S+)\s*=\s*([^;]+);", text, re.IGNORECASE):
            body = m.group(2)
            if ":" in body:                       # aln.phy: 1-100 -> drop filename
                body = body.split(":", 1)[1]
            body = re.sub(r"^(CODON|DNA|AA|MORPH|BIN)\b", "", body.strip(), flags=re.IGNORECASE)
            specs.append((m.group(1), body))
    for name, spec in specs:
        for start, end, step, raw in parse_partition_sites(spec):
            if start is None:
                errs.append("cannot parse site specification %r for partition %s "
                            "(expected START-END[,START-END] or *)" % (raw, name))
                continue
            if end is not None and end < start:
                errs.append("partition %s: range %s ends before it starts" % (name, raw))
                continue
            hi = end if end is not None else start
            if start < 1 or hi > n_sites:
                errs.append("partition %s: range %s lies outside the alignment "
                            "(%d sites)" % (name, raw, n_sites))
            if step < 1:
                errs.append("partition %s: \step must be >= 1 in %s" % (name, raw))
    return errs



def build(ns):
    """Return (cmd_list, warnings, errors)."""
    warns, errs = [], []
    style = ns.style
    parts = []

    def flag(name):
        modern, legacy, _src, _desc = FLAG_MAP[name]
        return legacy if style == "legacy" else modern

    exe = ns.binary
    if style == "legacy" and exe == "iqtree3":
        warns.append("--style legacy targets IQ-TREE <=1.x but --binary is iqtree3; "
                     "use 'iqtree' or switch to --style modern")

    # alignment / partition input ------------------------------------------------
    aln = ns.alignment
    if aln:
        if not os.path.exists(aln):
            if ns.allow_missing:
                warns.append("alignment not found: %s (--allow-missing given, building anyway)" % aln)
            else:
                errs.append("alignment not found: %s (use --allow-missing to build a command "
                            "for a file that is created later)" % aln)
        parts.extend(["-s", aln])
        _kind, detail = (None, None)
        if os.path.exists(aln):
            _kind, detail = sniff(aln)
        if _kind == "partition":
            errs.append("--alignment %s looks like a %s-style PARTITION file, not an "
                        "alignment; pass it via --partition/--partition-equal/--partition-unlinked "
                        "(command-reference.md 'Partition model options')" % (aln, detail))
        elif detail == "unknown":
            warns.append("could not recognise the format of %s; docs support PHYLIP, FASTA, "
                         "NEXUS, CLUSTAL and MSF (Tutorial.md). Force the datatype with -st if "
                         "IQ-TREE mis-detects it" % aln)
        elif ns.seqtype and detail == "unknown":
            warns.append("-st is only needed when IQ-TREE mis-detects the sequence type "
                         "(command-reference.md)")
    elif not ns.partition and not ns.fix_tree and not ns.user_tree:
        errs.append("no input given: IQ-TREE needs an alignment (-s) or a partition file "
                    "(-q/-p/-Q; version 1.x: -spp/-sp)")

    if ns.seqtype:
        token = ns.seqtype.upper()
        base = re.match(r"^(NT2AA|CODON|MORPH|DNA|AA|BIN)(.*)$", token)
        if not base:
            errs.append("-st %r is not one of DNA, AA, BIN, MORPH, CODON, NT2AA "
                        "(command-reference.md 'General options')" % ns.seqtype)
        else:
            parts.extend(["-st", token])
            if base.group(1) == "CODON":
                warns.append("-st CODON is mandatory for codon models, otherwise IQ-TREE "
                             "applies DNA models (command-reference.md)")

    if ns.model:
        m = ns.model
        up = m.upper()
        kw = up.split("+")[0]
        sel_kw = {"TEST", "TESTONLY", "TESTNEW", "TESTNEWONLY", "MF", "MFP",
                  "TESTMERGE", "TESTMERGEONLY", "TESTNEWMERGE", "TESTNEWMERGEONLY"}
        if kw in sel_kw or kw.startswith("TEST"):
            if "MERGE" in up and not (ns.partition or ns.partition_equal or ns.partition_unlinked):
                errs.append("-m %s needs a partition file (-p/-q/-Q; version 1.x -spp/-sp); "
                            "MERGE merges existing partitions (command-reference.md "
                            "'Automatic model selection')" % m)

            elif "MERGE" not in up and (ns.partition or ns.partition_equal or ns.partition_unlinked):
                warns.append("-m %s fits a model per partition but does NOT look for a better "
                             "partitioning scheme; add +MERGE (MFP+MERGE, TESTMERGEONLY) for that "
                             "(Advanced-Tutorial.md 'Choosing the right partitioning scheme')" % m)
        elif kw.startswith("MIX") and "{" not in m and "+MF" not in up:
            warns.append("-m MIX{A,B} takes the components in braces; MixtureFinder uses "
                         "-m MIX+MF / MIX+MFP (Complex-Models.md 'MixtureFinder')")
        parts.extend(["-m", m])
    elif aln:
        warns.append("no -m given: the default is ModelFinder (-m MFP since 1.5.4, Tutorial.md), "
                     "so this run tests models before searching")

    # threads / memory ----------------------------------------------------------
    if ns.threads is not None:
        val = str(ns.threads)
        if val.upper() == "AUTO":
            pass
        elif re.match(r"^[1-9]\d*$", val):
            if int(val) == 1:
                warns.append("%s 1 forces the serial search path; the single-thread binary "
                             "needs no thread flag at all (command-reference.md)" % flag("threads"))
        else:
            errs.append("%s %r must be a positive integer or AUTO" % (flag("threads"), ns.threads))
            val = None
        if val is not None:
            parts.extend([flag("threads"), val])
        if val is not None and val.upper() == "AUTO" and ns.threads_max is None:
            ncpu = os.cpu_count() or 0
            if ncpu > 16:
                warns.append("%s AUTO can probe all %d visible cores; consider "
                             "-ntmax <N> (command-reference.md)" % (flag("threads"), ncpu))
    if ns.threads_max is not None:
        if ns.threads is None or str(ns.threads).upper() != "AUTO":
            warns.append("-ntmax only bounds '%s AUTO'; it has no effect with "
                         "%s %s" % (flag("threads"), flag("threads"), ns.threads))
        parts.extend(["-ntmax", str(ns.threads_max)])
    if ns.memory:
        if not MEM_RE.match(ns.memory):
            errs.append("-mem %r should look like 8G / 512MB / 64gb "
                        "(command-reference.md example: -mem 64G)" % ns.memory)
        else:
            parts.extend(["-mem", ns.memory])

    # output --------------------------------------------------------------------
    if ns.prefix:
        parts.extend([flag("prefix"), ns.prefix])
    if ns.seed is not None:
        if not re.match(r"^-?\d+$", str(ns.seed)):
            errs.append("-seed must be an integer, got %r" % ns.seed)
        else:
            parts.extend(["-seed", str(ns.seed)])
    if ns.outgroup:
        og = ns.outgroup.replace(" ", ",")
        parts.extend(["-o", og])
    if ns.redo:
        parts.append("-redo")

    # trees ---------------------------------------------------------------------
    if ns.user_tree and ns.fix_tree:
        errs.append("-t (starting tree) and -te (fixed tree) are mutually exclusive")
    if ns.user_tree:
        parts.extend(["-t", ns.user_tree])
    if ns.fix_tree:
        parts.extend(["-te", ns.fix_tree])
    if ns.constraint:
        parts.extend(["-g", ns.constraint])
    for name, val in (("-t", ns.user_tree), ("-te", ns.fix_tree), ("-g", ns.constraint)):
        if val and not os.path.exists(val):
            warns.append("%s file not found: %s" % (name, val))

    # partitions ----------------------------------------------------------------
    pflags = [k for k in ("partition", "partition_equal", "partition_unlinked")
              if getattr(ns, k)]
    if len(pflags) > 1:
        errs.append("give at most one partition option, got: %s"
                    % ", ".join(flag(k) for k in pflags))
    for k in pflags:
        val = getattr(ns, k)
        parts.extend([flag(k), val])
        if not os.path.exists(val):
            warns.append("partition file not found: %s" % val)
        else:
            kind, detail = sniff(val)
            if kind != "partition":
                warns.append("%s does not look like a partition file (detected: %s/%s); "
                             "expected RAxML-style 'DNA, part1 = 1-100' lines or a NEXUS "
                             "sets block (Complex-Models.md)" % (flag(k), kind, detail))
            elif detail == "raxml":
                warns.append("%s with a RAxML-style file applies the single -m rate-heterogeneity "
                             "type to every partition; use a NEXUS partition file to mix +G and "
                             "+I+G (Advanced-Tutorial.md)" % flag(k))
                if k == "partition_unlinked":
                    warns.append("RAxML-style files carry no per-partition models; that is the "
                                 "point of the NEXUS charpartition command (Complex-Models.md)")
            elif detail == "nexus" and ns.model and re.match(
                    r"^[A-Za-z0-9]+(\+[A-Za-z0-9]+)*$", ns.model) \
                    and ns.model.upper() not in MODEL_KEYWORDS:
                warns.append("-m %s may override the models declared in the NEXUS charpartition; "
                             "for per-partition models leave -m out or use -m TEST/MFP "
                             "(Complex-Models.md)" % ns.model)
            n_sites = sequence_length(aln) if aln else None
            if n_sites:
                bad = check_partition_ranges(val, n_sites)
                if bad:
                    errs.extend(bad)
        if pflags and not aln:
            kind, detail = sniff(val) if os.path.exists(val) else (None, None)
            if detail != "nexus":
                warns.append("without -s, IQ-TREE can read the alignment names only from a NEXUS "
                             "partition file (Advanced-Tutorial.md: 'iqtree -p example.nex')")

    # bootstrap / tests ---------------------------------------------------------
    if ns.bootstrap is not None:
        try:
            b = int(ns.bootstrap)
        except (TypeError, ValueError):
            errs.append("bootstrap replicates must be an integer, got %r" % ns.bootstrap)
        else:
            if b < 1000:
                warns.append("%s %d < 1000; the docs specify >=1000 replicates for "
                             "ultrafast bootstrap (command-reference.md 'Ultrafast bootstrap "
                             "parameters')" % (flag("bootstrap"), b))
            parts.extend([flag("bootstrap"), str(b)])
            if ns.bootstrap_nni:
                warns.append("-bnni roughly doubles computing time; recommended only under "
                             "severe model violation (command-reference.md)")
    elif ns.bootstrap_nni:
        errs.append("-bnni requires the ultrafast bootstrap option (%s)" % flag("bootstrap"))
    if ns.alrt is not None:
        parts.extend(["-alrt", str(int(ns.alrt))])

    # misc ----------------------------------------------------------------------
    if ns.site_rates:
        parts.append(flag("site_rates"))
        if ns.stop_iter is None and not ns.fix_tree:
            warns.append("--rate runs a full tree search; add '-n 0' or '-te <tree>' "
                         "to save time (Advanced-Tutorial.md)")
    if ns.alninfo:
        parts.append("-alninfo")
    if ns.asr:
        parts.append("-asr")
    if ns.stop_iter is not None:
        parts.extend(["-n", str(ns.stop_iter)])
    if ns.nstop is not None:
        parts.extend(["-nstop", str(ns.nstop)])
    if ns.pers is not None:
        try:
            p = float(ns.pers)
        except (TypeError, ValueError):
            errs.append("-pers must be numeric, got %r" % ns.pers)
        else:
            if not 0.0 <= p <= 1.0:
                errs.append("-pers %s is outside the documented 0..1 range" % ns.pers)
            else:
                parts.extend(["-pers", ("%g" % p)])

    if ns.extra:
        toks = shlex.split(ns.extra)
        parts.extend(toks)
        known = set()
        for modern, legacy, _s, _d in FLAG_MAP.values():
            known.add(modern)
            known.add(legacy)
        for i, t in enumerate(toks):
            if t.startswith("--") and t not in KNOWN_LONG_OPTIONS:
                warns.append("%s is not documented in this mirror (kept verbatim); verify with "
                             "`iqtree -h` before relying on it" % t)
            elif re.match(r"^-[A-Za-z]+$", t) and t not in known \
                    and t not in KNOWN_SHORT_OPTIONS:
                warns.append("%s is not documented in this mirror (kept verbatim); verify with "
                             "`iqtree -h` before relying on it" % t)

    cmd = [exe] + parts
    return cmd, warns, errs


def quote(cmd):
    out = []
    for tok in cmd:
        if tok and re.match(r"^[A-Za-z0-9_./:@+=,\-{}\\]+$", tok):
            out.append(tok)
        else:
            out.append("'" + tok.replace("'", "'\\''") + "'")
    return " ".join(out)


SELF_TEST_FIXTURE = """:
7 28
Frog AAATTTGGTCCTGTGATTCAGCAGTGAT
Turtle CTTCCACACCCCAGGACTCAGCAGTGAT
Bird CTACCACACCCCAGGACTCAGCAGTAAT
Human CTACCACACCCCAGGAAACAGCAGTGAT
Cow CTACCACACCCCAGGAAACAGCAGTGAC
Whale CTACCACGCCCCAGGACACAGCAGTGAT
Mouse CTACCACACCCCAGGACTCAGCAGTGAT
"""


def self_test():
    import tempfile
    print("build_command.py --self-test")
    assert set(FLAG_MAP) == {
        "alignment", "model", "seqtype", "threads", "threads_max", "memory", "prefix",
        "seed", "outgroup", "user_tree", "fix_tree", "bootstrap", "bootstrap_nni", "alrt",
        "partition", "partition_equal", "partition_unlinked", "redo", "constraint",
        "site_rates", "alninfo", "asr", "stop_iter", "nstop", "pers"}, "FLAG_MAP drift"
    for k, v in FLAG_MAP.items():
        assert len(v) == 4 and v[0] and v[1], k
    tmp = tempfile.mkdtemp(prefix="iqtree-build-selftest-")
    aln = os.path.join(tmp, "tiny.phy")
    with open(aln, "w") as fh:
        fh.write(SELF_TEST_FIXTURE.lstrip(":"))
    nex = os.path.join(tmp, "tiny.nex")
    with open(nex, "w") as fh:
        fh.write("#nexus\nbegin sets;\n charset part1 = 1-14;\n"
                 " charset part2 = 15-28;\n"
                 " charpartition mine = HKY+G:part1, GTR+I+G:part2;\nend;\n")
    tr = os.path.join(tmp, "tiny.treefile")
    with open(tr, "w") as fh:
        fh.write("(Frog,Turtle,Bird,Human,Cow,Whale,Mouse);\n")

    def parse(argv):
        ap = make_parser()
        return ap.parse_args(argv)

    cmd, w, e = build(parse(["--alignment", aln]))
    assert not e, e
    assert any("ModelFinder" in x for x in w), w
    assert cmd[:2] == ["iqtree3", "-s"] and "-m" not in cmd, cmd
    print("  happy path (bare alignment): %s" % quote(cmd))

    cmd, w, e = build(parse(["--alignment", aln, "--model", "GTR+I+G", "--threads", "4",
                             "--memory", "8G", "--prefix", "run1", "--seed", "42",
                             "--outgroup", "Frog", "--bootstrap", "1000", "--alrt", "1000"]))
    assert not e, e
    line = quote(cmd)
    for frag in ["iqtree3", "-s %s" % aln, "-m GTR+I+G", "-T 4", "-mem 8G",
                 "--prefix run1", "-seed 42", "-o Frog", "-B 1000", "-alrt 1000"]:
        assert frag in line, (frag, line)
    print("  modern full command: %s" % line)

    cmd, w, e = build(parse(["--alignment", aln, "--style", "legacy", "--threads", "AUTO",
                             "--bootstrap", "1000", "--prefix", "run1", "--partition", nex]))
    assert not e, e
    line = quote(cmd)
    for frag in ["-nt AUTO", "-bb 1000", "-pre run1", "-spp %s" % nex]:
        assert frag in line, (frag, line)
    assert "-T" not in line and "-B" not in line and "--prefix" not in line, line
    print("  legacy full command: %s" % line)

    # error cases
    _, _, e = build(parse(["--alignment", aln, "--threads", "0"]))
    assert e, e
    _, _, e = build(parse(["--alignment", aln, "--memory", "eightG"]))
    assert e, e
    _, _, e = build(parse(["--alignment", aln, "--pers", "1.5"]))
    assert e, e
    _, _, e = build(parse(["--alignment", aln, "--seqtype", "PROTEIN"]))
    assert e, e
    _, _, e = build(parse(["--alignment", aln, "--user-tree", tr, "--fix-tree", tr]))
    assert e, e
    _, _, e = build(parse(["--alignment", aln, "--bootstrap-nni"]))
    assert e, e
    _, _, e = build(parse([]))
    assert e, e
    print("  rejects bad threads/mem/pers/seqtype/-t+-te/bare -bnni/no input: OK")

    # warnings
    _, w, _ = build(parse(["--alignment", aln, "--bootstrap", "100"]))
    assert any("1000" in x for x in w), w
    _, w, _ = build(parse(["--alignment", aln, "--threads", "AUTO"]))
    _, w2, _ = build(parse(["--alignment", aln, "--threads", "4", "--threads-max", "8"]))
    assert any("-ntmax" in x for x in w2), w2
    _, w3, _ = build(parse(["--alignment", aln, "--model", "MFP", "--partition", nex]))
    assert any("MERGE" in x for x in w3), w3
    _, w4, _ = build(parse(["--alignment", aln, "--site-rates"]))
    assert any("-n 0" in x for x in w4), w4
    print("  warns on low bootstrap, -ntmax without AUTO, missing MERGE, --rate cost: OK")

    for extra_val in ("-wsplits -quiet", "-m LG+F+G"):
        cmd, w, e = build(parse(["--alignment", aln, "--extra=" + extra_val]))
        assert not e, (e, extra_val)
        toks = shlex.split(extra_val)
        for t in toks:
            assert t in cmd, (t, cmd)
        print("  --extra=%-22s -> %s" % (extra_val, quote(cmd)))
    print("SELF TEST PASSED")
    return 0


def make_parser():
    ap = argparse.ArgumentParser(
        prog="build_command.py",
        description="Assemble a doc-verified IQ-TREE command line. Dry-run by default.",
        epilog="Only flags listed in FLAG_MAP are emitted; use --extra for anything else.")
    ap.add_argument("--alignment", "-a", help="input alignment (-s): PHYLIP/FASTA/NEXUS/CLUSTAL/MSF")
    ap.add_argument("--model", "-m", help="model string, e.g. GTR+I+G, LG+F+G, MFP, TESTMERGEONLY")
    ap.add_argument("--seqtype", help="DNA|AA|BIN|MORPH|CODON|NT2AA (only when auto-detection fails)")
    ap.add_argument("--threads", help="'%s' CPU cores or AUTO" % "-T/-nt")
    ap.add_argument("--threads-max", dest="threads_max", help="-ntmax cap for AUTO")
    ap.add_argument("--memory", help="-mem RAM cap, e.g. 8G")
    ap.add_argument("--prefix", help="output prefix (--prefix / -pre)")
    ap.add_argument("--seed", help="-seed random seed")
    ap.add_argument("--outgroup", help="-o outgroup taxon name(s), comma separated")
    ap.add_argument("--user-tree", dest="user_tree", help="-t starting/user tree")
    ap.add_argument("--fix-tree", dest="fix_tree", help="-te fixed tree (no search)")
    ap.add_argument("--bootstrap", help="ultrafast bootstrap replicates (-B / legacy -bb)")
    ap.add_argument("--bootstrap-nni", dest="bootstrap_nni", action="store_true", help="-bnni")
    ap.add_argument("--alrt", help="-alrt SH-aLRT replicates")
    ap.add_argument("--partition", help="edge-proportional partition file (-p / legacy -spp)")
    ap.add_argument("--partition-equal", dest="partition_equal",
                    help="edge-equal partition file (-q)")
    ap.add_argument("--partition-unlinked", dest="partition_unlinked",
                    help="edge-unlinked partition file (-Q / legacy -sp)")
    ap.add_argument("--constraint", help="-g topological constraint tree")
    ap.add_argument("--redo", action="store_true", help="-redo (overwrites all outputs!)")
    ap.add_argument("--site-rates", dest="site_rates", action="store_true",
                    help="per-site rates (--rate / legacy -wsr)")
    ap.add_argument("--alninfo", action="store_true", help="-alninfo")
    ap.add_argument("--asr", action="store_true", help="-asr ancestral states")
    ap.add_argument("--stop-iter", dest="stop_iter", help="-n iterations (0 skips tree search)")
    ap.add_argument("--nstop", help="-nstop unsuccessful iterations before stopping")
    ap.add_argument("--pers", help="-pers NNI perturbation strength 0..1")
    ap.add_argument("--extra", default="", help="verbatim extra arguments (quoted string)")
    ap.add_argument("--binary", default="iqtree3", help="executable name (default: %(default)s)")
    ap.add_argument("--style", choices=["modern", "legacy"], default="modern",
                    help="modern = 2.x/3.x spellings, legacy = 1.x spellings")
    ap.add_argument("--allow-missing", dest="allow_missing", action="store_true",
                    help="do not error when the alignment file does not exist yet")
    ap.add_argument("--run", action="store_true", help="execute instead of dry-run")
    ap.add_argument("--dry-run", action="store_true", dest="dry_run",
                    help="print only (this is already the default)")
    ap.add_argument("--self-test", action="store_true", dest="self_test",
                    help="run built-in fixture checks and exit")
    return ap


def main(argv=None):
    ap = make_parser()
    ns = ap.parse_args(argv)
    if ns.self_test:
        return self_test()
    cmd, warns, errs = build(ns)
    line = quote(cmd)
    for w in warns:
        print("WARN: %s" % w)
    if errs:
        for e in errs:
            print("ERROR: %s" % e)
        print("NOTHING RUN. Fix the errors above.")
        return 1
    print("COMMAND: %s" % line)
    if ns.run:
        print("--- executing ---")
        try:
            proc = subprocess.Popen(cmd[1:], executable=cmd[0])
        except OSError as exc:
            print("ERROR: cannot launch %s: %s" % (cmd[0], exc))
            return 127
        return proc.wait()
    print("(dry run: nothing executed; add --run to execute)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
