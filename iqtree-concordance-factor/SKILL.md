---
name: iqtree-concordance-factor
description: Quantify gene-tree/species-tree discordance with IQ-TREE's gene concordance factor (gCF) and site concordance factor (sCF), and build full concordance vectors with ASTRAL quartet frequencies. Use when bootstrap support looks high but loci disagree, when checking whether a clade is genuinely conflicting or merely uninformative, or when annotating a reference tree with per-branch gCF/sCF/gDF values. Triggers include --gcf, --scf, --scfl, loci.treefile, .trees/.treefile gene tree sets, .cf.stat, quartet_psi columns, and adding new sequences to an existing tree.
---

# Concordance Factors: gCF, sCF and Discordance

## Overview

For every branch of a reference (usually species) tree, **gCF** is the percentage of *decisive* gene trees containing that branch and **sCF** is the percentage of *decisive alignment sites* supporting it. They complement bootstrap: bootstrap measures support *given the data*, concordance factors measure how much the loci/sites actually *agree*, so a branch can carry 100% UFBoot and 5% gCF. Available since IQ-TREE 2; `--scfl` (likelihood-based sCF) since v2.2.2, and `-S`/directory inputs described here need IQ-TREE 3.

## Option quick reference (all verified in the references)

| Option | Meaning |
| --- | --- |
| `--gcf TREES_FILE` | Compute gCF from a set of gene/locus trees against the reference tree given with `-t`. |
| `--scf NUM` | Original parsimony-based sCF; `NUM` = number of random quartets sampled per internal branch (≥100 recommended). |
| `--scfl NUM` | Likelihood-based sCF (IQ-TREE ≥2.2.2), more accurate; same quartet count argument. Recommended over `--scf`. |
| `-t TREE_FILE` | Reference tree for `--gcf`; also the documented meaning of `-t` in the tree-file option groups ("file containing a set of trees"), so a multi-tree Newick set such as `loci.treefile` is a valid argument. |
| `-te TREE_FILE` | Fixed user tree — no tree search. Required with `--scfl` if you want to keep the reference topology. |
| `-s ALN_FILE` | Input alignment (PHYLIP, FASTA, NEXUS, CLUSTAL, MSF). |
| `-p PARTITION_FILE \| ALN_DIR` | Edge-linked partition file, or (IQ-TREE 3) a directory of locus alignments auto-concatenated into a supermatrix. |
| `-S PARTITION_FILE \| ALN_DIR` | IQ-TREE 3: model selection + tree inference **separately per partition/alignment**; output `.treefile` contains one tree per locus. |
| `-m MODEL` / `-m MFP` | Fix the model (skips ModelFinder) or let ModelFinder choose. |
| `-B 1000` | Ultrafast bootstrap (v2 syntax; use `-bb 1000` in v1.x) for comparison against CFs. |
| `-blfix` | Keep the branch lengths of the `-t`/`-te` tree — reuse step-1 estimates instead of re-optimising. |
| `-seed NUM` | Make quartet sampling reproducible (sCF changes slightly between runs otherwise). |
| `-T AUTO` / `-T 10` | Number of CPU cores. |
| `--prefix NAME`, `-pre NAME` | Prefix for all output files. Omitting it writes outputs as `concat.treefile.*` / after the partition file name. |
| `-g CONSTRAINT_TREE` | Topological constraint tree — the basis of the placement recipe below. |

## Workflow 1 — gCF + sCF on a phylogenomic dataset

### Step 1: infer the reference (species) tree

```
# single concatenated alignment + partition file
iqtree3 -s ALN_FILE -p PARTITION_FILE --prefix concat -B 1000 -T AUTO

# or IQ-TREE 3: point -p straight at a folder of locus alignments
iqtree3 -p ALN_DIR --prefix concat -B 1000 -T AUTO
```

Any reference tree works: concatenation/supermatrix, coalescent (ASTRAL), reconciliation, supertree, or a published tree.

### Step 2: infer the gene/locus trees

```
iqtree3 -s ALN_FILE -S PARTITION_FILE --prefix loci -T AUTO
# or
iqtree3 -S ALN_DIR --prefix loci -T AUTO
```

`-S` runs ModelFinder + tree search independently per locus. Outputs match a normal run except that `loci.treefile` holds **a set of trees**. Also produced and useful later: `loci.best_model.nex` (every estimated model with all parameter values, in NEXUS) and `loci.iqtree`.

### Gene-tree input requirements for `--gcf`

| Requirement | Documented behaviour |
| --- | --- |
| Multi-tree file | A plain multi-Newick text file, one tree per line — exactly what `-S` writes to `loci.treefile`. The `.trees` suffix used by ASTRAL-era pipelines is just this same format; the references never mandate an extension. |
| Taxon coverage | Gene trees may each contain a **subset** of the reference taxa; variable coverage is accounted for correctly. |
| Tree type | Any set of trees works — locus trees or bootstrap trees. |
| Branch labels / IDs | gCF is appended to existing branch labels after a `/`; use `*.cf.branch` to map table rows onto branches. |
| Polytomies and `#N` partition labels | **Not covered by these references** (unverified). AliSim documents multifurcating Newick generally, but nothing here states how `--gcf` treats unresolved gene trees or whether `#N` site-mapping labels are required — test on a small subset and check `gN` counts before trusting the numbers. |

### Step 3: gCF

```
iqtree3 -t concat.treefile --gcf loci.treefile --prefix concord
```

| File | Content |
| --- | --- |
| `concord.cf.tree` | NEWICK reference tree with gCF on each internal branch. If a branch already carries a label (e.g. bootstrap), gCF is **appended after a `/`**. |
| `concord.cf.branch` | NEWICK reference tree with internal **branch IDs**, matching the `ID` column of the stat table. |
| `concord.cf.stat` | Tab-separated table: gCF and gDF (gene discordance factors) per internal branch; read directly in R. |

`-t` accepts any reference tree and `--gcf` accepts any set of trees (locus trees **or** bootstrap trees); gene trees may contain a subset of the reference taxa, and variable taxon coverage is handled correctly.

### Step 4: sCF

```
# IQ-TREE >= 2.2.2 (recommended)
iqtree3 -te concat.treefile -s ALN_FILE --scfl 100 --prefix concord -seed 42
# older versions (parsimony sCF)
iqtree2 -t concat.treefile -s ALN_FILE --scf 100 --prefix concord -T 10
```

Instead of `-s` you can pass a directory or partition file (`-p ALN_DIR`), which computes sCF over the concatenated alignment. Combining both in one run:

```
iqtree3 -t concat.treefile --gcf loci.treefile -p ALN_DIR --scf 100 --prefix concord -T 10
```

Then `concord.cf.tree` branches are labelled `gCF/sCF` and `concord.cf.stat` carries both columns. Note the combination example uses the original `--scf` ("only for the original sCF").

### Speeding up large datasets

Step 1 already estimated the model, its parameters and the branch lengths — reuse them. Extract the model string from the `ALISIM COMMAND` block at the end of the report:

```
tail concat.iqtree
```

```
ALISIM COMMAND
--------------
--alisim simulated_MSA -t concat.treefile -m "Q.plant+I{0.177536}+R8{0.147295,0.0935335,...}" --length 432014
```

Take everything after `-m` (drop the `--length` part) and add it plus `-blfix`:

```
iqtree3 -te concat.treefile -p ALN_DIR --scfl 100 --prefix concord2 -blfix \
  -m "Q.plant+I{0.177536}+R8{0.147295,0.0935335,0.114418,0.190578,0.108376,0.538389,0.113777,0.804005,0.0898871,1.30004,0.137297,1.95653,0.0958285,3.48597,0.0155849,6.09904}"
```

Use **your own** model string. Without `-te`, `--scfl` triggers ModelFinder and a tree search; `-m` alone skips ModelFinder.

## Workflow 2 — full concordance vectors (IQ-TREE + ASTRAL)

A concordance vector describes all four outcomes for a branch, not just the winning one:

| Element | Definition |
| --- | --- |
| Ψ1 | concordance factor for the branch as it appears in the species tree |
| Ψ2 | larger of the two discordance factors from a **single NNI** of that branch |
| Ψ3 | smaller of those two NNI discordance factors |
| Ψ4 | remaining discordance (sum of the other topologies); site and quartet CFs assume this is 0 |

Requires IQ-TREE ≥ 2.3, ASTRAL-III, and R (`tidyverse`, `boot`, `ape`, `ggtext`) with the scripts from `https://github.com/roblanf/concordance_vectors` (`concordance_vector.R`, `concordance_table.R`, `change_labels.R`).

```
# 1. gene trees (per-locus models recorded in loci.best_model.nex)
iqtree3 -S bird_400 --prefix loci -T 128

# 2. species tree
astral -i loci.treefile -o astral_species.tree 2> astral_species.log

# 3. quartet concordance vector + posterior probabilities on the fixed topology
astral -q astral_species.tree -i loci.treefile -t 2 -o astral_species_annotated.tree 2> astral_species_annotated.log

# 4. site concordance vectors (reuses the per-locus models -> big time saving)
iqtree3 -te astral_species_annotated.tree -p loci.best_model.nex --scfl 100 --prefix scfl -T 128

# 5. gene concordance vectors (note: input is scfl.cf.tree so labels accumulate)
iqtree3 -te scfl.cf.tree --gcf loci.treefile --prefix gcf -T 128

# 6. dummy run only to export ASTRAL's coalescent-unit branch lengths in an IQ-TREE table
#    -blfix keeps original lengths and --scfl 1 makes these sCFs meaningless: never report them
iqtree3 -te astral_species_annotated.tree -blfix -p loci.best_model.nex --scfl 1 --prefix coalescent_bl -T 128

# 7. merge into one CSV, relabel the tree, print a table for a branch of interest
Rscript concordance_vector.R                 # -> concordance_vectors.csv
Rscript change_labels.R                      # -> id_gcf_scf_qcf.nex
Rscript concordance_table.R concordance_vectors.csv 598   # -> concordance_table_598.{csv,pdf}
```

Key outputs: `gcf.cf.stat` (gCF, gDF1, gDF2, gDFP, counts, plus every ASTRAL label), `scfl.cf.stat` (same for sites), `coalescent_bl.cf.stat`, `gcf.cf.tree` (annotated tree), `gcf.cf.branch` (branch IDs matching the tables).

### Reading the annotated branch label

```
'[q1=0.570241231975882;q2=0.17602481596980715;q3=0.25373395205431093;f1=207.567808439221;f2=64.0730330130098;f3=92.3591585477691 7;pp1=1.0;pp2=1.575119358351017E-20;pp3=2.952157354351003E-20;QC=8496;EN=364.0]'/25.4/47.0:0.0055956557
```

| Token | Source | Meaning |
| --- | --- | --- |
| `q1`,`q2`,`q3` | ASTRAL | quartet concordance vector (fractions; **not** support values) |
| `f1`,`f2`,`f3` | ASTRAL | quartet counts |
| `pp1`,`pp2`,`pp3` | ASTRAL | posterior probabilities (roughly, P that q1 is the largest) |
| `QC`, `EN` | ASTRAL | quartet count / effective number of gene trees |
| `/25.4/47.0` | IQ-TREE | sCF then gCF appended to the label |
| `:0.0055956557` | IQ-TREE | branch length in substitutions per site |

Open `gcf.cf.tree` in Dendroscope and declare the labels as **edge labels**. Labels apply to branches, not nodes.

### `concordance_vectors.csv` columns

| Column group | Columns | Derived from |
| --- | --- | --- |
| identifier | `ID` | branch ID in `gcf.cf.branch` |
| genes | `gene_psi1..4` (%) , `gene_psi1_N..4_N` (counts), `gene_N` | `gCF`, `gDF1`/`gDF2` (larger→Ψ2, smaller→Ψ3), `gDFP` (Ψ4), `*_N`, `gN` in `gcf.cf.stat` |
| sites | `site_psi1..3`, `site_psi1_N..3_N`, `site_N`; `site_psi4(_N)` always 0 | `sCF`, `sDF1`/`sDF2`, `sCF_N`, `sDF1_N`/`sDF2_N`, `sN` in `scfl.cf.stat` |
| quartets | `quartet_psi1..3`, `quartet_psi1_N..3_N`, `quartet_N`, `quartet_psi1_pp..3_pp` | `q1..q3`, `f1..f3`, `EN`, `pp1..pp3` parsed from the `Label` column of `scfl.cf.stat` |
| lengths | `length_subs_per_site`, `length_coalescent` | IQ-TREE concatenated analysis / ASTRAL |

Ψ2 is always the **larger** and Ψ3 the **smaller** of the two NNI discordance values; sort accordingly rather than assuming gDF1 > gDF2. Site/quartet counts can be non-integers because of how they are computed; the bootstrap CIs in `concordance_table.R` output come from 1000 resamples of those counts and get unreliable at very low counts.

## Workflow 3 — adding new sequences to an existing tree (placement)

Verified spelling: the reference tree enters as a **constraint tree via `-g`**, with the new sequences supplied in the alignment. There is **no `--epa` or `--place` option in these references** (unverified — do not use).

```
iqtree -s sequences.fa -g existing_tree.nex
```

- `sequences.fa` — alignment of the 17 established taxa (`T1`–`T17`) **plus** the newcomers (`NEW_1`–`NEW_3`).
- `existing_tree.nex` — tree containing only the 17 established taxa.

IQ-TREE fixes the T1–T17 topology, places each new sequence onto the constraint tree by Maximum Parsimony, then runs standard ML optimisation that moves only the new sequences. This corrects suboptimal parsimony placements and finds the best placement under the full likelihood model. Last tested with IQ-TREE 2.2.0.3.

Optional global re-optimisation once many sequences have been added — drop the constraint and use the placed tree as the starting tree:

```
iqtree -s sequences.fa -t sequences.treefile -pre re_optimised
```

Reproducible build of the example from scratch (needs `newick_utils`):

```
iqtree3 --alisim sequences -t RANDOM{yh/20} -seed 438579088 -af fasta -m JC
grep '>' sequences.fa
sed -i -e 's/T18/NEW_1/g' -e 's/T19/NEW_2/g' -e 's/T20/NEW_3/g' sequences.fa
nw_prune sequences.treefile T18 T19 T20 > existing_tree.nex
iqtree -g existing_tree.nex -s sequences.fa -pre seqs_added -seed 438579088
iqtree -s sequences.fa -t seqs_added.treefile -pre re_optimised -seed 438579088
```

For datasets beyond a few tens of thousands of sequences, consider UShER/matOptimize or MAPLE instead.

## Interpreting gCF and sCF

| Pattern | Interpretation |
| --- | --- |
| High UFBoot (e.g. 95) + high gCF/sCF | Loci and sites agree; the split is well supported. |
| High UFBoot + low gCF and low sCF | Real, strongly estimated **conflict** — the concatenated signal is dominated by many loci telling one story while most loci tell another. Investigate alternatives (ILS, hybridisation, hidden paralogy, contamination). |
| Low UFBoot + low gCF + high Ψ4-style spread across vectors | Little information, not necessarily conflict; look at counts (`gene_N`, `site_N`). |
| gCF ≫ sCF or sCF ≫ gCF | Few loci drive the pattern vs. many sites within loci; short/low-information loci make gene trees noisy, so sCF is informative when individual loci are uninformative. |
| gCF low but one alternative high (Ψ2 ≈ 100 − Ψ1) | Two competing topologies — a classic ILS/hybrid signature. Large Ψ4 means genuinely unresolved extra topologies. |

Caveats to state in any write-up: sCF is computed from randomly sampled quartets, so repeat runs differ slightly (fix `-seed`); sCF assumes Ψ4 = 0; and CFs describe proportions of genes/sites, not confidence in the branch.

## End-to-end worked example

Inputs: `loci/` — 200 FASTA alignments, one per locus.

```
# 1. species tree with bootstrap for comparison
iqtree3 -p loci --prefix concat -B 1000 -T AUTO
# 2. one ML tree per locus
iqtree3 -S loci --prefix loci_trees -T AUTO
# 3. gCF and likelihood sCF
iqtree3 -t concat.treefile --gcf loci_trees.treefile --prefix concord
iqtree3 -te concat.treefile -p loci --scfl 100 --prefix concord_scfl -seed 1 -T 10
```

Checks, in order:
1. `loci_trees.treefile` contains ~200 Newick lines (`grep -c ';' loci_trees.treefile`) — if only one, you used `-p` instead of `-S`.
2. `concord.cf.stat` has one row per internal branch of the species tree and an `ID` column; `concord.cf.branch` shows the same IDs on the tree.
3. In `concord.cf.tree`, internal branches carry labels such as `95/80` (bootstrap/gCF) — appended after `/`, not replacing the support value.
4. Per branch, gCF + gDF1 + gDF2 + gDFP = 100%; if not, some gene trees were excluded (taxon coverage) — inspect `gN`.
5. `concord_scfl.cf.stat` supplies sCF; compare with gCF using the interpretation table above.
6. For publication-quality vectors, rerun the ASTRAL workflow and report Ψ1–Ψ4 with bootstrap CIs rather than Ψ1 alone.

## Common pitfalls

| Symptom | Cause / fix |
| --- | --- |
| `--scfl` unknown option | IQ-TREE < 2.2.2 — upgrade or use `--scf` |
| `-S`/`-p DIR` rejected | Directory input needs IQ-TREE 3; otherwise supply a partition file |
| Tree search ran again during sCF | `--scfl` without `-te`; add `-te REF.treefile` (and `-m` to skip ModelFinder, `-blfix` to freeze lengths) |
| sCF differs between identical runs | Random quartet sampling — set `-seed` |
| Bootstrap values overwritten | They are not; CFs are appended after `/`. Read `.cf.stat` for clean numbers |
| Cannot link table rows to the tree | Use the branch-ID tree `*.cf.branch`; IDs match the `ID` column |
| Very unstable sCF | Below 100 quartets per branch — increase the `--scfl`/`--scf` number |
| Quartet labels missing from `.cf.stat` | The ASTRAL annotation step (`astral -q ... -t 2`) was skipped |
| Placement changed existing relationships | `-g` constrains the topology; check you did not pass the tree with `-t` (starting tree) instead |

## Loading only what you need

```bash
# Concordance-factor manual sections
grep -n -A 22 "^## Inferring species tree" references/concordance-factor.md
grep -n -A 18 "^## Gene concordance factor" references/concordance-factor.md
grep -n -A 35 "^## Site concordance factor" references/concordance-factor.md
grep -n -A 40 "^## Calculating concordance factors on very large" references/concordance-factor.md

# Concordance-vector recipe (long file: load only the needed stage)
grep -n -A 34 "^## Estimating concordance vectors and support" references/concordance-vectors.md
grep -n -A 62 "^## Generate the concordance vectors for each branch" references/concordance-vectors.md   # full column table
grep -n -A 26 "^## Put concordance factors" references/concordance-vectors.md
grep -n -A 25 "^### Confidence intervals on the concordance" references/concordance-vectors.md

# Phylogenetic placement recipe
grep -n -A 20 "^### Command line" references/recipe-phylogenetic-placement.md
grep -n -A 25 "^### How to build this example yourself" references/recipe-phylogenetic-placement.md
```

## Resources

- `references/concordance-factor.md` — official Concordance Factor manual: species-tree and locus-tree inference commands, `--gcf`/`--scf`/`--scfl` usage, the three `.cf.*` output files, combined gCF+sCF runs, and the `-blfix` + reused-model speed-up for large datasets. Load first for any CF task.
- `references/concordance-vectors.md` — worked ASTRAL + IQ-TREE recipe (~4k words): Ψ1–Ψ4 definitions, exact command sequence, annotated-label parsing, the complete `concordance_vectors.csv` column dictionary, tree relabelling and bootstrap CIs. Load when you need full vectors rather than single CF percentages.
- `references/recipe-phylogenetic-placement.md` — adding new sequences to an existing tree with `-g` constraint trees, optional unconstrained global re-optimisation, and a self-contained reproduction script. Load for placement/online-phylogenetics questions.

## Related skills

`iqtree-tree-search` (reference and locus tree inference), `iqtree-branch-support` (UFBoot/SH-aLRT values that CFs complement), `iqtree-model-selection` and `iqtree-complex-models` (`-m MFP`, edge-linked partitions, `loci.best_model.nex`), `iqtree-alisim` (simulate known-discordance benchmarks; `ALISIM COMMAND` block), `iqtree-rooting` (root before comparing topologies), `iqtree-dating`, `iqtree-command-reference`, `iqtree-recipes`.
