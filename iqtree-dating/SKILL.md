---
name: iqtree-dating
description: Build time-calibrated phylogenies in IQ-TREE with either Bayesian MCMCtree (IQ2MC) or Least-Squares Dating LSD2. Use for molecular dating, time trees, tip dates, sampling dates from pathogen outbreaks (HIV, influenza, SARS-CoV-2), fossil/ancestral-node calibrations, relaxed clocks, gradient/Hessian generation for MCMCtree, and confidence intervals on node dates. Triggers include --date, --dating mcmctree, .mcmctree.ctl, .hessian, .timetree.nex/.nwk, mcmctree.ctl, FigTree date display, --date-ci, --date-outlier, and a dated tree as output.
---

# IQ-TREE Molecular Dating

## Overview

Convert a substitution-scale ML tree into a time tree using one of two engines built into IQ-TREE: **LSD2** (least-squares dating, since IQ-TREE 2.0.3; fast, runs inside a normal `iqtree` call via `--date`) or **MCMCtree/PAML via the IQ2MC workflow** (Bayesian dating, since IQ-TREE 3.0.1; `--dating mcmctree` writes a Hessian + control file that you then feed to a modified `mcmctree` binary). Pick LSD2 for routine pathogen/phylogeographic tip-dating; pick IQ2MC when you need a Bayesian posterior, fossil calibration priors, or mixture-model dating.

## Decision table

| Goal | Engine | Core options |
| --- | --- | --- |
| Tip sampling dates -> time tree in one command | LSD2 | `iqtree -s ALN --date DATE_FILE` |
| Dates embedded after `\|` in sequence names | LSD2 | `--date TAXNAME` |
| Fossil / ancestral-node (MRCA) calibrations | LSD2 | `--date DATE_FILE --date-tip 0` |
| Date an existing fixed topology | LSD2 | `-te TREE_FILE --date DATE_FILE` |
| Confidence intervals on dates | LSD2 | `--date-ci NUM` (+ `--clock-sd`) |
| Remove long-branch outliers before dating | LSD2 | `--date-outlier 3` |
| Any native LSD2 flag not exposed by IQ-TREE | LSD2 | `--date-options "-e 2"` |
| Bayesian posterior node ages, birth-death prior | IQ2MC | `--dating mcmctree` then `mcmctree X.mcmctree.ctl` |
| Clock model choice (equal / independent / correlated) | IQ2MC | `--mcmc-clock EQUAL\|IND\|CORR` |
| Partitioned or mixture-model dating | both | `-Q`/`-p` partitions, `-m "MIX{...}"`, profile mixtures |

## Part A — LSD2 least-squares dating

### Required version note
The docs carry a disclaimer: download **version 2.0.6 or later** to get newer options such as `--date-ci`; the CI method itself is documented as unpublished.

### The IQ-TREE LSD2 front-end (full list from the docs)

```
TIME TREE RECONSTRUCTION:
  --date FILE          Dates of tips or ancestral nodes
  --date TAXNAME       Extract dates from taxon names after last '|'
  --date-tip STRING    Tip dates as a real number or YYYY-MM-DD
  --date-root STRING   Root date as a real number or YYYY-MM-DD
  --date-ci NUM        Number of replicates to compute confidence interval
  --clock-sd NUM       Std-dev for lognormal relaxed clock (default: 0.2)
  --date-outlier NUM   Z-score cutoff to exclude outlier nodes (e.g. 3)
  --date-options ".."  Extra options passing directly to LSD2
```

### Case A1 — tip dates (the pathogen scenario: HIV, influenza, SARS-CoV-2)

Prepare a date file: one line per taxon, name and date separated by spaces/tabs. You do **not** need dates for every tip.

```
hCoV-19/Wuhan-Hu-1         2019-12-31
hCoV-19/China/WF0028       2020-02
hCoV-19/USA/WA-S88         2020-03-01
hCoV-19/USA/CA-CDPH-UC1	   2020
hCoV-19/bat/Yunnan         2013-07-24
hCoV-19/pangolin/Guangdong 2019-02-01:2019-12-31
```

Date precision is flexible (`YYYY`, `YYYY-MM`, `YYYY-MM-DD`) and ranges use `lower:upper` with `NA` for an open bound:

```
TaxonA  2018-02-01:NA     # sampled after 1 Feb 2018
TaxonB  NA:2018-03-31     # sampled before 31 Mar 2018
```

Run one command that does ModelFinder + ML tree (substitutions/site) + rescaling to time:

```
iqtree -s ALN_FILE --date DATE_FILE
```

Outputs (prefix = alignment file name unless `--prefix` given):

| File | Content |
| --- | --- |
| `ALN_FILE.timetree.lsd` | LSD report |
| `ALN_FILE.timetree.nex` | Time tree in NEXUS — open in FigTree, Node Labels → Display = "date" |
| `ALN_FILE.timetree.nwk` | Time tree in NEWICK |

Root position is auto-detected by the LSD criterion, and a wrong root gives wrong dates. Supply outgroups explicitly:

```
iqtree -s ALN_FILE --date DATE_FILE -o "hCoV-19/bat/Yunnan,hCoV-19/pangolin/Guangdong"
```

Alternative without a separate date file — append dates to sequence names with `|` and use the keyword `TAXNAME`:

```
>hCoV-19/Wuhan-Hu-1|2019-12-31
iqtree -s ALN_FILE --date TAXNAME -o "hCoV-19/bat/Yunnan,hCoV-19/pangolin/Guangdong"
```

### Case A2 — ancestral-node / fossil calibrations

Left column is a **comma-separated MRCA taxon list with no spaces** (space is the taxon-list/date separator); right column is the calibration date. Negative numbers mean "million years ago" relative to present.

```
taxon1,taxon2	      -50
taxon3,taxon4,taxon5  -100
taxon6                -10
```

```
iqtree -s ALN_FILE --date DATE_FILE --date-tip 0
```

`--date-tip 0` sets all undated tips to the present. If the root age is known, pin it with `--date-root` (real number or `YYYY-MM-DD`).

### Case A3 — dating a tree you already have

```
iqtree -s ALN_FILE --date DATE_FILE -te TREE_FILE
```

Add `-m MODEL` (or a partition file with models) to skip ModelFinder and speed the run up.

### Case A4 — confidence intervals and outliers

```
iqtree -s ALN_FILE --date DATE_FILE --date-ci 100          # resample branch lengths 100x
iqtree -s ALN_FILE --date DATE_FILE --date-ci 100 --clock-sd 0.4
iqtree -s ALN_FILE --date DATE_FILE --date-outlier 3       # z-score cutoff
```

- `--date-ci` is **not** a bootstrap; it is much faster and based on a mixture of Poisson and lognormal distributions for a relaxed clock (docs note it is unpublished).
- `--clock-sd` = lognormal SD (default 0.2); larger → wider intervals.
- `--date-outlier` = z-score threshold; higher value → more outliers removed. Long branches bias dates, so check this on uneven datasets.

### Case A5 — passthrough to native LSD2 (`--date-options`)

Any LSD2 flag can be injected. Verified mapping between the IQ-TREE wrapper and native LSD2 flags:

| IQ-TREE option | Equivalent LSD2 flag (lsd2-help.md) | Native semantics |
| --- | --- | --- |
| `--date-tip` | `-z tipsDate` | single date applied to all tips |
| `--date-root` | `-a rootDate` | root date; quote strings like `"2020-01-10"` or `"b(2019,2020)"` |
| `--date-ci` | `-f samplingNumberCI` | simulated-branch-length CIs; written into the NEXUS as label `CI_height` (FigTree "Node bar"); uses `-s` sequence length, default 1000 |
| `--clock-sd` | `-q standardDeviationRelaxedClock` | lognormal SD for CI simulation, default 0.2 |
| `--date-outlier` | `-e ZscoreOutlier` | normalize branch residuals; ~3 typical; cannot be combined with multi-rate estimation or imprecise constraints |
| `--date FILE` | `-d inputDateFile` | first line = number of constraints; supports `l(v)`, `u(v)`, `v`, `b(v1,v2)`, `mrca(a,b,c)`, node labels |

Other useful native-only flags to pass through `--date-options "..."`:

| Flag | Purpose |
| --- | --- |
| `-s sequenceLength` | sequence length used to draw integer branch lengths for CIs (default 1000); smaller → wider CIs under non-strict clocks |
| `-r rootingMethod` | `l` local re-root, `a` all branches, `as` constrained on all branches, `k` keep current branch |
| `-g outgroupFile` / `-G` | root on outgroups; `-G` deletes them from the result |
| `-v variance` | 0 none, 1 use input branch lengths (default), 2 two-pass |
| `-b varianceParameter` | pseudo-count `b` in `Vi = (Bi+b)`; small → strict clock, large → relaxed clock |
| `-p partitionFile` | branch subsets with different rates: `group1 1 {n1} {n5 n4}` |
| `-R round_time` | round minimum branch length to a time unit (default 365 ⇒ days) |
| `-u minBlen` / `-U minExBlen` | force minimum internal / external branch length in time units |
| `-l nullBlen` | informative-branch threshold (default `0.5/seq_length`) |
| `-S minSupport` | collapse internal branches with support ≤ value |
| `-w givenRte` | file of pre-set substitution rates (one per input tree) |
| `-D outDateFormat` | 1 = real number, 2 = year-month-day |
| `-F` | drop the ancestor-descendant date ordering constraint (faster) |
| `-t rateLowerBound` | lower bound on estimated rate (default 1e-10) |
| `-m samplingNumberOutlier` | dated nodes sampled for outlier detection (default 10) |
| `-n datasetNumber` | how many trees to read from the input tree file |
| `-j` | verbose messages |

Example from the docs:

```
iqtree -s ALN_FILE --date DATE_FILE --date-options "-e 2"
```

### Model guidance for pathogen dating

The LSD2 examples in the docs are nucleotide alignments (`-s ALN_FILE` with ModelFinder choosing e.g. GTR+G), and the dating step only rescales the ML substitution-scale tree, so what matters is that the **substitution model fits the marker**:

- Use DNA models for the usual RNA-virus genome / gene fragments; let ModelFinder choose (`-m MFP` behaviour is the default when `-m` is omitted) or fix it with `-m GTR+G4` to save time.
- For amino-acid or codon inputs, specify the corresponding model class with `-m` (see `iqtree-model-selection` / `iqtree-complex-models`); the IQ2MC control file documents the matching `seqtype` values used by MCMCtree: `0` nucleotides, `1` codons, `2` amino acids.
- Always give outgroups (`-o`) or a root date (`--date-root`) for outbreak datasets — an unrooted/auto-rooted tree is the main source of badly wrong tMRCA estimates.

## Part B — Bayesian dating with MCMCtree (IQ2MC)

Available from **IQ-TREE 3.0.1**. Cite Demotte et al. (2025), *IQ2MC: A New Framework to Infer Phylogenetic Time Trees Using IQ-TREE 3 and MCMCTree with Mixture Models*.

### Three-step workflow

| Step | What happens | Mandatory? |
| --- | --- | --- |
| 1 | IQ-TREE searches the ML tree and (optionally) the best-fit model | Skippable if you supply rooted tree + MSA + model |
| 2 | With `--dating mcmctree`, IQ-TREE computes the **gradient vector and Hessian** at ML and writes all MCMCtree input files | Required |
| 3 | Run `mcmctree` on the generated control file to sample node ages | Required |

MCMCtree accepts **rooted trees only** — root the tree yourself (FigTree/iTOL) if your input tree is unrooted, and add fossil/tip calibrations before step 3 if the tree lacks them.

### Step 2 command

```
iqtree3 -s example.phy -m GTR+G4 -te example_tree.nwk --dating mcmctree --prefix example
```

Files produced:

| File | Role |
| --- | --- |
| `example.mcmctree.hessian` | gradient vector + Hessian for approximate likelihood |
| `example.rooted.nwk` | rooted tree **compatible with the Hessian ordering** — must be the tree you hand to MCMCtree |
| `example.mcmctree.ctl` | PAML control file, runnable as-is |
| `example.dummy.aln` | placeholder alignment (current MCMCtree format still demands one); saves compute vs. passing the real MSA |

Workflow-shaping options written into the control file:

| Option | Arguments | Meaning |
| --- | --- | --- |
| `--mcmc-iter` | `burnin,samplefreq,nsample` | e.g. `20000,200,50000` |
| `--mcmc-bds` | `birth-rate,death-rate,sampling-fraction` | birth-death prior, e.g. `1,1,0.5` |
| `--mcmc-clock` | `EQUAL` \| `IND` \| `CORR` | global equal rates / independent rates / autocorrelated rates |

```
iqtree3 -s example.phy -m GTR+G4 -te example_tree.nwk \
  --mcmc-iter 20000,200,50000 --mcmc-bds 1,1,0.5 --mcmc-clock IND
```

Total MCMC iterations = `burnin + sampfreq × nsample`.

### Partitions and mixture models with `--dating mcmctree`

| Model type | Option | Hessian content |
| --- | --- | --- |
| Edge-unlinked partitions (EUL) | `-Q example.nex` | separate gradient + Hessian **per partition** |
| Edge-linked partitions (EL) | `-p example.nex` | one gradient vector + one Hessian (branches shared) |

```
iqtree3 -s example.phy -Q example.nex -m GTR+G4 -te example_tree.nwk --dating mcmctree
```

DNA/AA mixture models and profile mixtures are supported directly in `-m`:

```
iqtree3 -s example.phy -m "MIX{GTR,HKY}+G4" -te example_tree.nwk --dating mcmctree
iqtree3 -s example_aa.phy -m LG+G4+C60 -te example_aa_tree.nwk --dating mcmctree
```

Let MixtureFinder pick the mixture with `-m MIX+MFP`.

### Step 3 — running MCMCtree

Download the **modified PAML/MCMCtree** from `https://github.com/iqtree/paml` (it adds the `ckpfile` and `hessianfile` keywords). Then:

```
mcmctree example.mcmctree.ctl
```

Control-file skeleton (edit burnin/sampfreq for convergence before running):

```
seed = -1                        * current time used when seed < 0
seqfile = example.dummy.phy      * dummy alignment, only to let MCMCtree run
treefile = example.rooted.nwk    * rooted NEWICK with annotated fossil/tip dates
mcmcfile = example.mcmctree.log  * parameter traces, inspect in Tracer
outfile = example.mcmctree.out   * summarized MCMCtree results
ckpfile = example..mcmctree.ckp  * checkpoint file
hessianfile = example.mcmctree.hessian

checkpoint = 1     * 0 nothing; 1 save; 2 resume
ndata = 1          * number of partitions
seqtype = 0        * 0 nucleotides; 1 codons; 2 AAs (unused under approximate likelihood)
usedata = 2        * 0 priors only; 1 exact slow likelihood; 2 approximate likelihood
clock = 2          * 1 equal rates; 2 independent rates; 3 correlated rates
RootAge = <1.0     * safe constraint if the rooted tree has no root fossil

BDparas = 1 1 0.5  * birth-rate, death-rate, sampling fraction
finetune = 1: 0.1  0.1  0.1  0.01 .5  * auto : times, musigma2, rates, mixing, paras, FossilErr
print = 1          * 1 normal; 2 verbose

rgene_gamma = 2 2    * Dirichlet-gamma prior for mean rate across loci (clock=2 or 3, ndata>1)
sigma2_gamma = 1 10  * Dirichlet-gamma prior for rate variance across loci

burnin = 20000
sampfreq = 200
nsample = 50000

model = 0        * only for exact likelihood (usedata = 1): 0 JC69,1 K80,2 F81,3 F84,4 HKY85
alpha = 0        * only for usedata = 1
ncatG = 0        * only for usedata = 1
cleandata = 0    * only for usedata = 1
kappa_gamma = 6 2
alpha_gamma = 1 1
```

Compatibility fallback: with an unmodified MCMCtree release, delete the `ckpfile` and `hessianfile` lines from the control file and rename the Hessian file to `in.BV`.

## End-to-end worked examples

### Example 1 — SARS-CoV-2-style tip dating with LSD2

Inputs: `ncov.fasta` (alignment) and `dates.tsv`:

```
hCoV-19/Wuhan-Hu-1	2019-12-31
hCoV-19/Italy/SPL1	2020-01-29
hCoV-19/bat/Yunnan	2013-07-24
```

Command:

```
iqtree -s ncov.fasta --date dates.tsv -o "hCoV-19/bat/Yunnan" --date-ci 100 --date-outlier 3 -pre sarscov2
```

Check in the output:
1. `sarscov2.timetree.nwk` exists and its branch lengths are now in **time units**, not substitutions/site.
2. `sarscov2.timetree.lsd` reports a plausible substitution rate; root-to-tip distances should track sampling dates.
3. Open `sarscov2.timetree.nex` in FigTree → Node Labels → Display = `date`; with `--date-ci` also enable the `CI_height` node bar.
4. Confirm the excluded taxa listed by `--date-outlier 3` are genuinely long-branch artefacts, not real fast-evolving lineages.

### Example 2 — deep phylogeny with fossil (MRCA) calibrations

`calib.txt` (no spaces inside taxon lists):

```
taxon1,taxon2	-50
taxon3,taxon4,taxon5	-100
```

```
iqtree -s mammals.phy --date calib.txt --date-tip 0 --date-root -120 -m GTR+G4 -te mammals_ml.treefile
```

Check: node ages in the `.timetree.nwk` are older than their descendant tips, the root matches the imposed `--date-root`, and the calibration nodes correspond to the intended MRCA (verify with the branch-ID tree if unsure).

### Example 3 — Bayesian dating of a partitioned dataset (IQ2MC)

```
# Step 2: gradients + Hessian, edge-unlinked partitions, independent rates
iqtree3 -s example.phy -Q example.nex -m GTR+G4 -te example_tree.nwk \
  --dating mcmctree --prefix example \
  --mcmc-clock IND --mcmc-bds 1,1,0.5 --mcmc-iter 20000,200,50000

# Inspect/edit the control file (increase burnin if needed), then Step 3:
mcmctree example.mcmctree.ctl
```

Check:
1. `example.rooted.nwk` (not your original tree) is what `treefile` points at — the Hessian is ordered against it.
2. `example.mcmctree.log` traces look stationary in Tracer; raise `burnin`/`sampfreq` otherwise.
3. Summarized node-age posteriors with HPDs are in `example.mcmctree.out`.
4. Total iterations = 20000 + 200×50000; scale down for a test run first.

## Common pitfalls

| Symptom | Cause / fix |
| --- | --- |
| Dates wildly off | Wrong root — add `-o OUTGROUP` or `--date-root`; or use LSD2 `-r` via `--date-options` |
| `--date-ci` unknown | IQ-TREE < 2.0.6 — upgrade |
| Spaces in calibration file ignored | Taxon list must be comma-separated with **no** spaces; space separates list from date |
| MCMCtree errors on `ckpfile`/`hessianfile` | Using stock PAML — remove those lines and rename the Hessian to `in.BV` |
| MCMCtree rejects tree | Tree is unrooted, or you passed your own tree instead of `*.rooted.nwk` |
| Outlier removal fails with imprecise dates | Native LSD2 `-e` cannot combine with multiple rates or interval constraints |
| Dummy alignment name mismatch | Control file ships `seqfile = example.dummy.phy` while IQ-TREE lists `example.dummy.aln` — point `seqfile` at the file that actually exists |

## Loading only what you need

```bash
# LSD2 section of the dating manual (tip dates, calibrations, CI, outliers)
grep -n -A 60 "^## Inferring time tree with tip dates" references/dating.md

# IQ2MC / MCMCtree sections
grep -n -A 30 "^## Estimating the gradients" references/dating.md
grep -n -A 30 "^## Using partitions and mixture" references/dating.md
grep -n -A 55 "^## How to run MCMCtree" references/dating.md

# One specific native LSD2 flag (options are tab-indented in this file)
grep -n -A 8 "samplingNumberCI" references/lsd2-help.md      # -f : confidence intervals
grep -n -A 10 "rootingMethod" references/lsd2-help.md        # -r : rooting methods
grep -n -A 8 "ZscoreOutlier" references/lsd2-help.md         # -e : outlier exclusion
grep -n "^	-" references/lsd2-help.md                        # index of every flag with its line number

# IQ-TREE's own LSD2 front-end and worked examples
grep -n -A 20 "^## Least Square Dating (LSD2)" references/dating.md
```

## Resources

- `references/dating.md` — official "Phylogenetic Dating" manual: IQ2MC three-step workflow, `--dating mcmctree` output files, `--mcmc-*` options, full MCMCtree control file, and all LSD2 scenarios (tip dates, ancestral dates, dating an existing tree, CIs, outliers). Load for any dating task.
- `references/lsd2-help.md` — verbatim `lsd2 -h` output (LSD v1.8): every native LSD2 flag with defaults. Load when constructing `--date-options "..."` strings or debugging rate/variance/rooting behaviour.

## Related skills

`iqtree` (quickstart/orchestration), `iqtree-model-selection` and `iqtree-complex-models` (the `-m`/partition/mixture specifications dating depends on), `iqtree-tree-search` (getting the ML tree to date), `iqtree-rooting` (correct root before LSD2/MCMCtree), `iqtree-alisim` (simulate dated-tree benchmarks), `iqtree-command-reference`, `iqtree-recipes`.
