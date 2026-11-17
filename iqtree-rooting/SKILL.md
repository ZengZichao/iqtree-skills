---
name: iqtree-rooting
description: "Root, re-root and measure root support on IQ-TREE phylogenies. Use when a tree must be rooted without an outgroup, when asking for \"rootstrap\" or root support values, when running non-reversible substitution models (UNREST/12.12 and other Lie Markov models with RY/WS/MK symmetries, NONREV / NONREV+FO amino-acid models, NQ.* non-reversible matrices), when interpreting the report lines \"Tree is UNROOTED although outgroup taxon ... is drawn at root\" versus \"Tree is ROOTED at virtual root __root__\" in .treefile/.iqtree, when using -o, --model-joint, --root-test, -B/-bb bootstrap, -st NT2AA, .best_scheme.nex or .rootstrap.nex files, or when testing which branch a root should sit on with an AU/KH/SH test (.roottest.csv)."
---

# IQ-TREE: Rooting Trees and Rootstrap Support

## Overview

IQ-TREE infers **unrooted** trees by default because time-reversible models cannot locate a root; this skill covers the two supported ways to get a rooted tree — the outgroup approach (`-o`, plus manual re-rooting in a viewer) and the non-reversible **Rootstrap** approach (`--model-joint UNREST`/`NONREV` with `-B 1000`), including how to read `.rootstrap.nex` values and how to statistically test root positions with `--root-test`.

Use IQ-TREE **2.1.3 or later**; the Rootstrap features below are documented as requiring it.

## The rooting problem in one table

| Situation | What IQ-TREE does | Command core |
| --- | --- | --- |
| Reversible model (JC, HKY, GTR, LG, ...) | Infers an **unrooted** ML tree; the report prints `NOTE: Tree is UNROOTED although outgroup taxon '<first taxon>' is drawn at root` | `iqtree3 -s aln.phy -B 1000 -T AUTO --prefix rev` |
| Reversible model + named outgroup | `-o` specifies an outgroup taxon name to root the tree; "The output tree in `.treefile` will be rooted accordingly". Default outgroup = first taxon in the alignment | `... -o Pig` |
| Non-reversible model | Infers a **rooted** tree; report prints `NOTE: Tree is ROOTED at virtual root '__root__'` and a virtual `__root__` taxon appears in the ASCII/Newick tree | `... --model-joint UNREST` |
| Non-reversible model + UFBoot | Additionally writes `<prefix>.rootstrap.nex` with rootstrap support on every branch | `... --model-joint NONREV -B 1000` |

Key conceptual point: with reversible models the *placement of the name at the left of the printed tree carries no information* — IQ-TREE simply draws the first sequence of the alignment there. Never interpret that as a biological root.

## Recommended workflow (no outgroup available)

Two steps, because re-doing partition/model search under a parameter-rich non-reversible model is slow and over-parameterised.

```bash
# STEP 1 — unrooted tree + best partitioning scheme with REVERSIBLE models
iqtree3 -s bovidae.phy -p bovidae.nex -B 1000 -T AUTO --prefix rev_dna
#   writes rev_dna.treefile, rev_dna.iqtree, rev_dna.best_scheme.nex, ...

# STEP 2 — re-use the best scheme, force ONE LINKED NON-REVERSIBLE model, get a rooted tree
iqtree3 -s bovidae.phy -p rev_dna.best_scheme.nex --model-joint UNREST \
        -B 1000 -T AUTO --prefix nonrev_dna
#   writes nonrev_dna.treefile (ROOTED at __root__) and nonrev_dna.rootstrap.nex
```

Amino-acid version of the same two steps (protein-coding DNA translated on the fly):

```bash
# step 1: infer unrooted tree with reversible models
iqtree3 -s bovidae.phy -p bovidae.nex -B 1000 -T AUTO -st NT2AA --prefix rev_aa

# step 2: infer rooted tree with linked non-reversible models
iqtree3 -s bovidae.phy -p rev_aa.best_scheme.nex --model-joint NONREV \
        -B 1000 -T AUTO -st NT2AA --prefix nonrev_aa
```

Option meanings used above (all verified in the references):

| Option | Meaning |
| --- | --- |
| `-s` | Input alignment (PHYLIP, FASTA, NEXUS, CLUSTAL, MSF). |
| `-p` | Partition file, edge-proportional partition model (partitions may have different evolutionary speeds). Legacy spelling: `-spp`. |
| `-m TEST` / `-m MFP` | Model selection (`MFP` = like `-m MF` — extended candidate set incl. FreeRate — but immediately followed by tree reconstruction with the best-fit model). |
| `--model-joint <MODEL>` | Use one linked instance of `<MODEL>` across **all** partitions. Needed to avoid over-parameterisation: `UNREST` has 12 parameters, `NONREV` has 379. |
| `-B 1000` / `-bb 1000` | Ultrafast bootstrap, 1000 replicates (`-bb` is the legacy spelling). Rootstrap needs bootstrap replicates to count root positions. |
| `-T AUTO` / `-nt AUTO` | Auto-detect the best number of CPU cores. |
| `-st NT2AA` | Translate protein-coding DNA into AA and analyse the translated alignment. |
| `--prefix` / `-pre` | Prefix for all output files; use a distinct prefix per step so nothing is overwritten. |
| `-o <taxon>` | Outgroup taxon name used to root the tree in `.treefile`. |

## Non-reversible model catalogue

### DNA: Lie Markov models

Lie Markov models (since v1.6) form a nested series and many are non-reversible; they lack the inconsistency problem of plain GTR. In the names below the number before the dot is the number of rate parameters. Full list (Rev? = reversible?, Freq = number of free base frequencies: 0 = equal, 1 = f(A)=f(G) & f(C)=f(T), 2 = f(A)+f(G)=0.5=f(C)+f(T), 3 = unconstrained):

| Model | Rev? | Freq | Note |
| --- | --- | --- | --- |
| 1.1 | Yes | 0 | equiv. to JC |
| 2.2b | Yes | 0 | equiv. to K2P |
| 3.3a | Yes | 0 | equiv. to K3P |
| 3.3b | No | 0 | |
| 3.3c | Yes | 0 | equiv. to TNe |
| 3.4 | Yes | 1 | |
| 4.4a | Yes | 3 | equiv. to F81 |
| 4.4b | Yes | 1 | |
| 4.5a / 4.5b | No | 1 | |
| 5.6a | No | 0 | |
| 5.6b | No | 3 | |
| 5.7a | No | 2 | |
| 5.7b / 5.7c | No | 0 | |
| 5.11a | No | 2 | |
| 5.11b / 5.11c | No | 0 | |
| 5.16 | No | 1 | |
| 6.6 | No | 1 | equiv. to STRSYM, strand-symmetric model |
| 6.7a | No | 3 | F81+K3P |
| 6.7b / 6.8a | No | 3 | |
| 6.8b | No | 1 | |
| 6.17a / 6.17b | No | 1 | |
| 8.8 / 8.10a | No | 3 | |
| 8.10b | No | 1 | |
| 8.16 / 8.17 / 8.18 | No | 3 | |
| 9.20a | No | 2 | |
| 9.20b | No | 0 | Doubly stochastic |
| 10.12 / 10.34 | No | 3 | |
| **12.12** | **No** | 3 | **equiv. to UNREST (unrestricted model)** — the most general non-reversible DNA model, 12 parameters |

Every Lie Markov model can take a symmetry prefix:

| Prefix | Meaning |
| --- | --- |
| `RY` | purine-pyrimidine pairing (default) |
| `WS` | weak-strong pairing |
| `MK` | aMino-Keto pairing |

Include Lie Markov models in automatic model selection (v1.6+) by appending a keyword to `-m` (table reproduced in `../iqtree-recipes/references/command-reference.md`, section "Automatic model selection"):

| Option | Meaning |
| --- | --- |
| `-m TEST+LM` | additionally consider all Lie Markov models |
| `-m TEST+LMRY` | only those with RY symmetry |
| `-m TEST+LMWS` | only those with WS symmetry |
| `-m TEST+LMMK` | only those with MK symmetry |
| `-m TEST+LMSS` | only strand-symmetric Lie Markov models (`6.6` = STRSYM) |

### Amino acids

| Model string | Meaning |
| --- | --- |
| `NONREV` | Most general non-reversible amino-acid model, 379 parameters. Link it across partitions with `--model-joint NONREV`. |
| `NONREV+FO` | Same, with ML-optimised amino-acid frequencies (`+FO` = optimised, letter O not zero). Used when estimating a new non-reversible matrix (nQMaker workflow). |
| `NQ.pfam`, `NQ.mammal`, `NQ.bird`, `NQ.insect`, `NQ.plant`, `NQ.yeast` | Pre-existing **non-reversible** empirical Q matrices (Dang et al., 2022). Usable directly with `-m`. |
| `GTR20+FO` | Reversible GTR20 amino-acid model with ML-optimised frequencies — the documented step-2 `--model-joint` model when estimating a new non-reversible matrix; step 1 builds a guide tree with a reversible model such as `LG` (`../iqtree-model-selection/references/estimating-aa-models.md`). |

Frequency modifiers: `+F` empirical, `+FO` ML-optimised, `+FQ` equal.

### Naming hazard: do not confuse these

| String | Refers to | Not |
| --- | --- | --- |
| `GTR2` | General time-reversible model for **binary** data (unequal state frequencies) | a PoMo model, and not "GTR version 2" |
| `GTRX` (also invoked as `GTR`) | GTR with unequal rates for **morphological / multistate** data | DNA GTR |
| `MK` | Jukes-Cantor-type model with equal rates for **morphological** data | the Lie Markov aMino-Keto symmetry prefix |
| `GTR20` | Reversible amino-acid model with 190 rate parameters | a non-reversible DNA model |

If you need a PoMo (polymorphism-aware) `GTR2.0`/`GTR6.0` style model, that is a population-model feature — switch to the `iqtree-poms` skill; those names are **not** documented in the rooting references.

## Reading the outputs

| File | Content |
| --- | --- |
| `<prefix>.iqtree` | Main report. Look for the `NOTE: Tree is UNROOTED ...` / `NOTE: Tree is ROOTED at virtual root '__root__'` line and the ASCII tree with UFBoot supports. |
| `<prefix>.treefile` | ML tree in Newick. Rooted according to `-o`; carries the `__root__` node when a non-reversible model was used. Open in FigTree/iTOL to re-root on another branch. |
| `<prefix>.rootstrap.nex` | Written by the non-reversible + bootstrap run: the ML tree with **rootstrap** support annotated on every branch. In FigTree choose Branch Labels → Display → `rootstrap`. |
| `<prefix>.best_scheme.nex` | Best partitioning scheme from step 1; feed back to step 2 with `-p`. |
| `<prefix>.roottest.csv` | Written by `--root-test`: log-likelihood and topology-test p-values for every candidate rooting. |

### What rootstrap values mean

Rootstrap: a bootstrap analysis yields a set of **rooted** bootstrap trees under the non-reversible model; the rootstrap support of a branch of the ML tree is the **proportion of those rooted bootstrap trees whose root falls on that branch**. Values are computed for all branches, including external branches, and **sum to ≤ 1** along the tree (it is a distribution over branches, not a per-branch confidence). The sum is < 1 when some bootstrap replicates are rooted on a branch absent from the ML tree.

Worked reading of the documented Bovidae example:

| Analysis | Result | Interpretation |
| --- | --- | --- |
| DNA `UNREST` (`nonrev_dna`) | best branch (leading to `Cow`) = 72%, second-best branch (separating Cow+Yak from rest) = 17.9% | DNA model cannot reliably locate the root; it only narrows candidates. |
| AA `NONREV` (`nonrev_aa`) | branch separating Yak+Cow from the rest = 99.9% | High confidence; agrees with the independent outgroup rooting. |

General guidance from the source study: amino-acid non-reversible models have much more power to detect the root than DNA non-reversible models. Run both and prefer the one with concentrated support, but check that the winning position is also compatible with any outgroup/molecular-clock evidence.

## Testing specific root positions (`--root-test`)

Rootstrap is a support measure, not a statistical test. To compare log-likelihoods of the tree rooted on each branch, combine the non-reversible run with the tree-topology tests:

```bash
iqtree3 -s bovidae.phy -p rev_aa.best_scheme.nex --model-joint NONREV -st NT2AA \
        --root-test -zb 1000 -au -te nonrev_aa.treefile --prefix nonrev_aa_test
```

| Option | Meaning |
| --- | --- |
| `--root-test` | Re-root the tree on **every** branch and run the topology tests for you (v2.1.3+). |
| `-te <tree>` | Fix this user tree; no tree search, just likelihood computation. |
| `-zb 1000` | Number of RELL replicates (>= 1000) for the topology tests. |
| `-au` | Additionally perform the approximately unbiased (AU) test; requires `-zb`. |

`<prefix>.roottest.csv` columns (branches sorted by descending log-likelihood; comment lines start with `#`, so read in R with `dat=read.csv('nonrev_aa_test.roottest.csv',comment.char='#')`):

| Column | Meaning |
| --- | --- |
| `ID` | Branch ID |
| `logL` | Log-likelihood of the tree rooted at this branch |
| `deltaL` | logL difference from the maximal logL |
| `bp-RELL` | Bootstrap proportion using the RELL method |
| `p-KH` | One-sided Kishino-Hasegawa p-value |
| `p-SH` | Shimodaira-Hasegawa p-value |
| `c-ELW` | Expected likelihood weight |
| `p-AU` | Approximately unbiased test p-value |

Interpretation rule used in the docs: the top branch had `p-AU = 0.9976` while all others had `p-AU < 0.01`, so every alternative rooting is rejected — the root sits on branch ID 1. Map IDs back onto the tree by opening `<prefix>.rootstrap.nex` in FigTree and selecting Branch Labels → Display → `id`.

Tip from the docs: `--root-test -zb 1000 -au` can be added to the rootstrap command itself, giving rootstrap supports **and** rooting-test p-values in one run.

## Pitfalls

- Do not trust the leftmost taxon of a reversible-model tree as a root; re-root explicitly with `-o` or in a viewer.
- `UNREST`/`NONREV` are extremely parameter-rich — always link them with `--model-joint` on partitioned data rather than letting each partition get its own copy.
- Non-reversible runs are slower than reversible ones; keep step 1 (scheme + tree search) separate and reuse `.best_scheme.nex` and `-te`.
- Low, spread-out rootstrap supports (e.g. 72% / 17.9%) mean "undetermined", not "the 72% branch is right". Report the whole distribution.
- Rootstrap sums to ≤ 1; never rescale it to look like a percentage-support value comparable to UFBoot.
- If the non-reversible root contradicts a well-chosen outgroup, treat it as evidence of model misspecification rather than a new root.

## Related skills

- `iqtree-model-selection` — `-m MFP`/`TEST`, frequency and rate modifiers, Lie Markov inclusion in the model set.
- `iqtree-complex-models` — partition files, `.best_scheme.nex`, linked vs unlinked partition parameters.
- `iqtree-tree-search` — `-te`, `-g` constraints, starting trees.
- `iqtree-branch-support` — UFBoot (`-B`/`-bb`), SH-aLRT, AU/topology tests behind `--root-test`.
- `iqtree-concordance-factor` — mapping supports onto a fixed reference tree.
- `iqtree-poms` — polymorphism-aware models (the place for `POM`/`GTR2.0`-style population models).
- `iqtree-command-reference` — full option tables.
- `iqtree` — install, datatypes, output-file overview, troubleshooting.

## Resources

- `references/rootstrap.md` — the rooting guide: outgroup workflow with `bovidae_outgroup.phy`/`bovidae.nex`, the two-step non-reversible Rootstrap commands (`--model-joint UNREST`/`NONREV`), the rootstrap definition box, example ASCII trees with `__root__`, and the `--root-test` `.roottest.csv` walkthrough. Load first for any rooting or root-support question.
- `references/substitution-models.md` — all DNA models with 6-digit rate codes, frequency/rate-heterogeneity modifiers, the complete Lie Markov table (Rev?/Freq, incl. `12.12` = UNREST) with `RY`/`WS`/`MK` prefixes, amino-acid matrices including the non-reversible `NQ.*` set, binary/morphological models (`JC2`, `GTR2`, `MK`, `GTRX`) and `+ASC`. Load when choosing or spelling a model string.
- Cross-check outside this folder (read-only mirrors of the same docs): `../iqtree-recipes/references/command-reference.md` for the `+LM`/`+LMRY`/`+LMWS`/`+LMMK`/`+LMSS` model-selection keywords and `-T`/`-B`/`--prefix` spellings, and `../iqtree-model-selection/references/estimating-aa-models.md` for the `--model-joint NONREV+FO` estimation workflow.
- Targeted greps (avoid loading whole files):

```bash
cd "$(dirname "$0")"   # or: cd /path/to/iqtree-rooting
grep -n -A 45 "^### Lie Markov models" references/substitution-models.md
grep -n -A 8  "^All Lie Markov models can have" references/substitution-models.md
grep -n -A 3  "^| Model | Rev? | Freq" references/substitution-models.md
grep -n -A 12 "^## Rate heterogeneity across sites" references/substitution-models.md
grep -n -B 2 -A 30 "^## Inferring rooted trees without outgroup" references/rootstrap.md
grep -n -A 30 "^## Testing root positions" references/rootstrap.md
```
