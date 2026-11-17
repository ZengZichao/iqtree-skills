---
name: iqtree-tree-search
description: Run and control IQ-TREE's ML tree search, fix or constrain a topology, score a given tree, test alternative trees with KH/SH/AU/ELW/deltaL, build consensus or graft bootstrap values onto a fixed tree, compute Robinson-Foulds distances, reconstruct ancestral states and site-specific rates, trim misleading sites, and resume interrupted runs. Triggers include tree search, NNI and SPR, starting tree, stopping criteria, .treefile, .treels, .contree, .suptree, user tree, -g constraint, topology tests, AU test, deltaL, RF distance, -asr, site rates, and checkpoint restart with -redo.
---

# Tree search and topology testing in IQ-TREE

## Overview

IQ-TREE builds an ML tree from 100 parsimony starting trees plus a BIONJ tree, refines it with randomized NNI over a small candidate set, and stops after `nstop` unsuccessful iterations. Everything downstream of that search — fixing a user topology, constraining it, comparing candidate topologies, summarising a set of trees, measuring their distances, and reading out per-site/per-node quantities — is driven by the flags below.

All flags here were verified against `references/command-reference.md` (CR) and `references/advanced-tutorial.md` (AT). Short single-dash forms are the documented ones; where AT shows a long form (`--rate`, `--mlrate`, `--robust-phy`, `--prefix`) that form is used too.

## The search algorithm as documented

| Stage | What the docs say | Knobs |
| --- | --- | --- |
| Starting trees | `-t` default = **100 parsimony trees + BIONJ tree**; `-t BIONJ` starts from BIONJ only, `-t RANDOM` from a completely random tree (CR:68) | `-t`, `-djc` (skip ML distances + BIONJ), `-ninit` (default 100 initial parsimony trees, uses randomized stepwise addition + parsimony SPR), `-sprrad` (SPR radius for that parsimony stage, default 6) |
| Candidate set | Top `-ntop` (default 20) parsimony trees are optimized by ML NNI to seed the candidate set; `-nbest` (default 5) trees are kept during the search (CR:406-407) | `-ntop`, `-nbest` |
| Rearrangements | Randomized NNI with perturbation strength `-pers` (0–1, default 0.5); `-allnni` considers **all** possible NNIs instead of only those near previously applied NNIs — slower but more thorough (default OFF) (CR:400, 409) | `-pers`, `-allnni` |
| Fast mode | `-fast`: only two starting trees (maximum parsimony + BIONJ), optimized by NNI; also raises the final-estimation epsilon to 0.05 (CR:402, 701) | `-fast`, `-me` |
| Stopping | `-nstop` = number of unsuccessful iterations (default 100); `-n` = number of iterations to stop and **overrides `-nstop`**; `-n 0` skips the tree search entirely (CR:405, 408, 130/593) | `-n`, `-nstop` |
| Local optima | Docs' own advice: repeat the analysis with **at least 10 IQ-TREE runs**; `-pers` and `-nstop` are "the most relevant options to change". For many short sequences use smaller perturbation and more patience, e.g. `-pers 0.2 -nstop 500` (CR:412) | `-pers`, `-nstop` |

```
iqtree -s data.phy -m TEST -nstop 500 -pers 0.2      # CR example: harder search
iqtree -s data.phy -m TEST -allnni                   # exhaustive-ish NNI pass
```

## Supplying, fixing and scoring a user tree

| Option | Documented meaning |
| --- | --- |
| `-t FILE` | **Starting** tree — search continues and may change the topology. Also accepts a *set* of trees for `-con`, `-net`, `-rf`, `-sup` post-analysis. |
| `-te FILE` | "Like `-t` but fixing user tree. That means, no tree search is performed and IQ-TREE computes the log-likelihood of the fixed user tree." (CR:69) |
| `-blfix` | Fix branch lengths of the tree passed via `-t`/`-te`; "useful to evaluate the log-likelihood of an input tree with fixed topology and branch lengths" (default OFF) (CR:697) |
| `-n 0` | Estimate model parameters on an initial parsimony tree only, skip the full search (CR:593) |
| `-o NAME` | Outgroup used to root the `.treefile` output; default = first taxon in the alignment (CR:70) |
| `-blmin` / `-blmax` | Min/max branch length; defaults: smaller of 0.000001 and 0.1/alignment_length, and 10 (CR:698-699) |
| `-czb` | Collapse near-zero branches so the final tree may be multifurcating (CR:700) |

Likelihood of one fixed tree, topology **and** branch lengths frozen:

```
iqtree -s example.phy -te my.treefile -blfix -m GTR+G --prefix fixed
```

Ancestral-state reconstruction along a user tree instead of the ML tree: `-asr` combined with `-te` (CR:548, 554).

## Constrained (topology-restricted) search

`-g FILE` gives a topological constraint tree in NEWICK; it "can be a multifurcating tree and need not to include all taxa" (CR:403, AT:187).

```
iqtree -s data.phy -m TEST -g constraint.tree          # CR:425
iqtree -s example.phy -m TIM2+I+G -g example.constr1 --prefix example.constr1   # AT:228
```

Two documented caveats (AT:220-250):

- On an **unrooted** tree a constraint like `((Human,Seal),(Cow,Whale));` only forces the split — Human and Seal need not be monophyletic. To force them as a clade, put the extra taxa inside the constraint, e.g. `((Human,Seal),(Cow,Whale),Mouse);`.
- Always check a constrained result against the unconstrained one with topology tests (next section); AT calls this "advised" because the constrained tree may be significantly worse.

Placement variant: give an incomplete tree with `-t` and let the search fill in the rest. With CMAPLE this is explicit: `--tree-search FAST` = "placement only, no SPR moves" (CR:436, 481).

## Tree topology tests (`-z`, `-zb`, `-zw`, `-au`)

| Option | Meaning (CR:587-594) |
| --- | --- |
| `-z FILE` | File containing a **set of trees**; IQ-TREE computes the log-likelihood of every tree. |
| `-zb NUM` | Number of RELL replicates (**>=1000**) for BP, KH, SH and ELW over the `-z` trees. |
| `-zw` | With `-zb`: additionally weighted-KH and weighted-SH. |
| `-au` | With `-zb`: additionally the AU test. The replicate count still comes from `-zb`. |
| `-n 0` / `-te` | Skip the search; `-te` estimates parameters on your tree instead of a parsimony tree. |
| `-wsl` | Write site log-likelihoods to `.sitelh` in TREE-PUZZLE format, for external CONSEL runs (CR:706). |

Reading the `USER TREES` block of `.iqtree` (AT:354-375):

| Column | Meaning | How to judge |
| --- | --- | --- |
| `logL` | Log-likelihood of that tree | Higher is better |
| `deltaL` | Log-likelihood difference from the **best** tree in the set | 0.000 marks the best tree |
| `bp-RELL` | Bootstrap proportion by the RELL method | A **posterior weight, not a p-value**; weights sum to 1 across trees |
| `p-KH` | One-sided Kishino-Hasegawa p-value | Designed for exactly 2 trees → no multiple-testing correction |
| `p-SH` | Shimodaira-Hasegawa p-value | Fixes KH's problem but gets conservative when many trees are tested |
| `c-ELW` | Expected likelihood weight | Weight, not a p-value |
| `p-AU` | Approximately unbiased p-value | Recommended replacement for KH and SH |
| `+` / `-` flag | `+` = inside the 95% confidence set; `-` = significantly rejected | Reject a tree when its p-value < 0.05 (marked `-`) |

Other documented facts: duplicate topologies in the `-z` set are detected and dropped automatically (AT:289); IQ-TREE's AU uses the **least-squares** p-value estimate whereas CONSEL uses MLE, so values differ slightly (CR:585); IQ-TREE's AU supports SSE/AVX/multicore and resamples sites *within* partitions, unlike CONSEL (CR:583). Trees with optimized branch lengths land in `<prefix>.<zfile>.trees` — AT:259 names `example.phy.treels.trees` for a `-z example.treels` run, while the sample report at AT:357 prints `See example.phy.trees`; check both.

## Consensus, supertree-style summary and grafting support values

| Option | Meaning (CR:614-624) |
| --- | --- |
| `-t FILE` | The set of trees to summarize (NEWICK or NEXUS, e.g. from MrBayes). |
| `-con` | Consensus tree of the `-t` trees → written to `.contree`. |
| `-minsup FRAC` | Minimum fraction (0–1) of trees a branch must appear in to be kept. Default **0** = extended majority-rule; `-minsup 0.5` = majority-rule. |
| `-bi NUM` | Burn-in: discard the first NUM trees of the `-t` file. |
| `-net` | Consensus network of the trees passed via `-t` → written to `.nex`. AT:405 also shows `iqtree -net mytrees`. |
| `-sup TARGET SET` | Support values are taken from the tree set and mapped onto the "target" tree; result → `.suptree`. AT:411 shows the two files as positional arguments (`-sup input_tree set_of_trees`), the FAQ example uses `-sup input_alignment.treefile -t alltrees`. |
| `-suptag NAME\|ALL` | Tag nodes of `.suptree` with the IDs of the trees in which they occur. |
| `-scale NUM` / `-prec NUM` | Scaling factor (default 100, i.e. percentages) and decimal precision (default 0) of `-sup` values. |
| `-wbt` / `-wbtl` | During UFBoot, write bootstrap trees to `.ufboot` (with branch lengths for `-wbtl`) (CR:499-500). |

```
iqtree -con mytrees                          # extended majority-rule
iqtree -con mytrees -minsup 0.5 -bi 100      # majority-rule, drop 100 burn-in trees
iqtree -net mytrees                          # consensus network -> .nex
iqtree -sup input_tree set_of_trees          # support onto a fixed topology -> .suptree
```

Documented multi-run bootstrap merge (FAQ "How do I save time for standard bootstrap?"): run `-bo` per machine with distinct `-pre`, `cat boot*.boottrees > alltrees`, then `iqtree -con -t alltrees`, then re-estimate lengths with `iqtree -s input_alignment -te alltrees.contree -pre alltrees.contree`, and finally `iqtree -sup input_alignment.treefile -t alltrees`.

## Robinson-Foulds distance

| Option | Meaning (CR:630-635) |
| --- | --- |
| `-t FILE` | First set of trees. |
| `-rf FILE2` | Second set; computes **all pairwise** RF distances between the `-t` set and the `-rf` set. |
| `-rf_all FILE` | All-to-all RF distances within one set. |
| `-rf_adj FILE` | RF distances between **adjacent** trees in the file (MCMC-style trace inspection). |

```
iqtree -rf tree_set1 tree_set2
iqtree -rf_all tree_set
```

## Ancestral state reconstruction

| Option | Meaning (CR:546-554) |
| --- | --- |
| `-asr` | Write empirical-Bayesian ancestral sequences for **all** nodes to `.state`. |
| `-asr-min P` | Minimum posterior probability to call the best ancestral state; default = observed state frequency from the alignment. |
| `-te FILE` | Reconstruct along this user tree. Named nodes keep their names; unnamed nodes become `Node1`, `Node2`, … |

```
iqtree -s example.phy -m JC+G -asr
```

`example.phy.state` is tab-separated with columns `Node`, `Site`, `State`, then `p_A p_C p_G p_T` (posterior probability per state, empirical Bayes) — loadable straight into R with `read.table(..., header=TRUE)`.

## Site-specific rates

| Option | Meaning | Output |
| --- | --- | --- |
| `-wsr` (CR:364) / `--rate` (AT:456) | Empirical-Bayes per-site rate = posterior-weighted mean over rate categories; requires a Gamma/FreeRate-type heterogeneity model | `<aln>.rate` with columns `Site`, `Rate`, `Category`, `Categorized_rate` (category 1 = slow … 4 = fast) |
| `--mlrate` (IQ-TREE 3, AT:493) | Per-site rates estimated by **maximum likelihood** | `<aln>.mlrate` with columns `Site`, `Rate` |
| `-wspr` / `-wslr` | Site posterior probabilities / site log-likelihoods per rate category | `.siteprob` / `.sitelh` |

```
iqtree -s example.phy --rate                                   # full search
iqtree -s example.phy -t ml.treefile -n 0 --rate               # reuse a known ML tree
iqtree -s example.phy -m GTR+R10 -n 0 --rate                   # known model, no selection
iqtree -s example.phy -n 0 --mlrate                            # ML site rates
```

`-n 0` uses a parsimony tree only (fast, coarser); `-t ml.treefile` improves the estimate; `-m` skips ModelFinder (AT:477-491).

## Trimming alignment sites by likelihood ("trimmed log-likelihood")

Requires IQ-TREE ≥ 3.0 and is invoked as `--robust-phy` (AT:521-545): during the search it repeatedly recomputes per-site likelihood under the current tree/model and removes a user-defined proportion of the lowest-likelihood sites.

```
iqtree3 -s <MY_ALIGNMENT> --robust-phy <PROPORTION_TO_RETAIN> -m <MODEL>
iqtree3 -s data.phy --robust-phy 0.98 -m GTR+G -wsl     # keep 98%, also write site LHs
```

Add `-wsl` to get `.sitelh`, then identify the removed sites from their log-likelihoods. Reference: Q. Liu, B.Q. Minh, R. Lanfear, M.A. Charleston, S.A. Richards & B.R. Holland (2025) *Robust Phylogenetics*, bioRxiv, doi:10.1101/2025.04.01.646540.

## Checkpointing and resuming

| Item | Documented behaviour |
| --- | --- |
| `.ckp.gz` | Checkpoint file (gzip), written periodically during a run. Rerunning the **exact same command line** resumes from the last stopped point and appends to the previous log. |
| Finished run | Rerunning after success raises `ERROR: Checkpoint (...) indicates that a previous run successfully finished`. |
| `-redo` | Redo everything regardless of state; **overwrites all existing output files**. |
| `-cptime SEC` | Minimum checkpoint interval in seconds, default **20**. |
| Deleting the file | Destroys all checkpointing information (CR:118). |
| `.model` / `.model.gz` | ModelFinder checkpoint of every model tested so far; reloading reuses that computation (CR:206). |

```
iqtree -s example.phy -nt AUTO -mem 8G -pre myrun     # interrupted -> rerun identical line
iqtree -s example.phy -redo                           # force a fresh run
```

## Pathogen / very low-divergence search (CMAPLE)

`--pathogen` uses CMAPLE when sequence divergence is low, otherwise the normal IQ-TREE engine; `--pathogen-force` always uses it (CR:428-445). `--tree-search FAST|NORMAL|EXHAUSTIVE` chooses placement-only / SPR at newly-added nodes / SPR at all nodes (default NORMAL). Branch supports: `--alrt NUM`, `--sprta`, `--sprta-zero-branch`, `--sprta-other-places`; threads `-T`; mutations `--estimate-MAT`. **When CMAPLE is active only these options plus `-s`, `-t`, `-m`, `-blfix`, `--prefix`, `-redo`, `-st` apply — everything else is ignored.** Handoff pattern: run the fast search, then feed `pathogen.treefile` to a normal run with `-te`.

## Worked example 1 — constrained vs unconstrained, settled by the AU test

Follows AT:303-378 exactly (four searches, one concatenated tree file, one test run).

```
# 1. unconstrained ML tree
iqtree -s example.phy -m TIM2+I+G --prefix example.unconstr
# 2-4. three constrained searches, one constraint per file
iqtree -s example.phy -m TIM2+I+G -g example.constr1 --prefix example.constr1
iqtree -s example.phy -m TIM2+I+G -g example.constr2 --prefix example.constr2
iqtree -s example.phy -m TIM2+I+G -g example.constr3 --prefix example.constr3
# 5. stack the four .treefiles into one tree set (Windows: use `type` instead of `cat`)
cat example.unconstr.treefile example.constr1.treefile example.constr2.treefile \
    example.constr3.treefile > example.treels
# 6. score all four and run the tests
iqtree -s example.phy -m TIM2+I+G -z example.treels -n 0 -zb 1000 -au
```

Reported `USER TREES` table (AT:359-364):

```
Tree      logL    deltaL  bp-RELL    p-KH     p-SH    c-ELW     p-AU
  1   -21152.617   0.000  0.7110 + 0.7400 + 1.0000 + 0.6954 + 0.7939 +
  2   -21156.802   4.185  0.2220 + 0.2600 + 0.5910 + 0.2288 + 0.3079 +
  3   -21158.579   5.962  0.0670 + 0.1330 + 0.5130 + 0.0758 + 0.1452 +
  4   -21339.596 186.980  0.0000 - 0.0000 - 0.0000 - 0.0000 - 0.0000 -
```

Read-out: trees 1–3 carry `+` under `p-AU`, so all stay inside the 95% confidence set — the three alternatives are statistically indistinguishable and the position of `Seal` is undecided. Tree 4 carries `-` in every test column and sits `deltaL = 186.98` below the best tree, so the `(Human,Mouse)` / `(Cow,Rat)` grouping is rejected. Conclusion matches the SH-aLRT/UFBoot supports from the basic tutorial, and justifies keeping (or dropping) the constraint.

## Worked example 2 — bootstrap trees → consensus → support on the ML tree → RF distance

Follows the FAQ recipe for splitting standard bootstrap across machines.

```
# 1. five independent 20-replicate runs on 5 machines (distinct prefixes, no overwriting)
iqtree -nt 4 -s input_alignment -bo 20 ... -pre boot1
iqtree -nt 4 -s input_alignment -bo 20 ... -pre boot2
iqtree -nt 4 -s input_alignment -bo 20 ... -pre boot3
iqtree -nt 4 -s input_alignment -bo 20 ... -pre boot4
iqtree -nt 4 -s input_alignment -bo 20 ... -pre boot5
# 2. pool the bootstrap tree files
cat boot*.boottrees > alltrees
# 3. consensus of the pooled trees            -> alltrees.contree
iqtree -con -t alltrees
# 4. re-estimate branch lengths of the consensus on the real alignment
iqtree -s input_alignment -te alltrees.contree -pre alltrees.contree
# 5. ordinary ML analysis on the original alignment
iqtree -nt 4 -s input_alignment ...
# 6. graft the bootstrap frequencies onto that ML tree   -> .suptree
iqtree -sup input_alignment.treefile -t alltrees
# 7. how different are the two tree sets, in splits?
iqtree -rf input_alignment.treefile alltrees
```

Read-out: step 3 writes an extended majority-rule consensus (add `-minsup 0.5` for strict majority-rule, `-bi N` to drop MCMC burn-in trees); step 4 turns that topology into a dated-by-lengths tree without letting the search move the topology (`-te` = no tree search); step 6 reports, per branch of the ML tree, how many of the pooled trees contain it, scaled by `-scale` (default 100 → percent) and rounded to `-prec` decimals; step 7 gives the pairwise Robinson-Foulds counts between the ML tree and the bootstrap set — large RF values mean the bootstrap set is not sampling around your ML topology.

## Fetching details from the references

```bash
cd "$(dirname "$0")"   # run from this skill folder

# Search-algorithm knobs and the local-optima advice
grep -n -A 20 "^## Tree search parameters" references/command-reference.md
# Fixed/user trees, starting-tree defaults, branch-length fixing
grep -n -A 30 "^## General options" references/command-reference.md
grep -n -A 22 "^## Miscellaneous options" references/command-reference.md
# Constraint tree semantics and worked constraints
grep -n -A 62 "^## Constrained tree search" references/advanced-tutorial.md
# Topology tests: option table, then the worked USER TREES reading
grep -n -A 28 "^## Tree topology tests" references/command-reference.md
grep -n -A 50 "^## Tree topology tests" references/advanced-tutorial.md
grep -n -A 46 "^## Testing constrained tree" references/advanced-tutorial.md
# Consensus / network / support mapping
grep -n -A 16 "^## Constructing consensus tree" references/command-reference.md
grep -n -A 33 "^## Consensus construction" references/advanced-tutorial.md
# Robinson-Foulds
grep -n -A 24 "^## Computing Robinson-Foulds distance" references/command-reference.md
# Ancestral states
grep -n -A 34 "^## Ancestral sequence reconstruction" references/command-reference.md
# Site rates and site trimming
grep -n -A 70 "^## Inferring site-specific rates" references/advanced-tutorial.md
grep -n -A 26 "^## Trimming alignment sites" references/advanced-tutorial.md
# Checkpointing, pathogen/CMAPLE search
grep -n -A 10 "^## Checkpointing to resume" references/command-reference.md
grep -n -A 60 "^## Tree search for pathogen" references/command-reference.md
# Single-flag lookups
grep -nE "^\| (-ninit|-ntop|-nbest|-nstop|-pers|-sprrad|-allnni|-fast|-g|-te|-z|-zb|-zw|-au|-asr|-con|-rf|-sup|-blfix|-redo|-cptime) " references/command-reference.md
```

Every pattern above returns hits. Note that `^## Tree topology tests` exists in **both** reference files, so name the file you want.

## Flags searched for but NOT documented

Searched the whole mirror (`grep -rin`) and found zero occurrences — do not use them: `--nni`, `--i-quartet`, `--n-quart`, `--epsilon`, `--n-spr`, `--n-best-spr`, `--start-spr`, `--num-rand-add`, `--constrain-branch-length`, `--tree-into-trees`, `--sup-tree`, `--rf-trees`, `--rft`, `--consensus`, `--map-anc-state`, `--anc-freq`, `--site-rates`, `--trimalias`, `--pretrim`, `--posttrim`. Real equivalents: `-allnni`/`-pers`, `-n`/`-nstop`/`-ninit`/`-ntop`/`-nbest`/`-sprrad`, `-me`, `-g`, `-te`/`-blfix`, `-z`/`-zb`/`-zw`/`-au`, `-t`+`-con`/`-net`/`-sup`, `-rf`/`-rf_all`/`-rf_adj`, `-asr`/`-asr-min`, `-wsr`(`--rate`)/`--mlrate`, `--robust-phy`. (`--site-rate` exists but belongs to AliSim, not to tree search.)

## Related skills

`iqtree-model-selection` (`-m TEST`/`MFP` before the search), `iqtree-branch-support` (`-bb`, `-alrt`, `-abayes`, `-lbp` values that the topology tests complement), `iqtree-complex-models` (partition/mixture models used in constrained and `-z` runs), `iqtree-concordance-factor` (`--gcf`/`--scf` on tree sets), `iqtree-rooting` (root the `.treefile` after a search), `iqtree-dating` (time trees from a fixed topology), `iqtree-alisim` (simulated benchmarks for search tests), `iqtree-command-reference`, `iqtree-recipes`.
