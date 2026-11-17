---
name: iqtree-recipes
description: "End-to-end IQ-TREE recipes that chain IQ-TREE with other tools. Use for (1) phylogenetic placement of new/short query sequences onto an existing reference tree with a constraint tree -g, a masked reference alignment plus -te, or EPA-ng/pplacer-style jplace output, followed by global re-optimisation; (2) Goldman-Cox tests of model adequacy for stationarity of base/amino-acid frequencies, reading .gc and .cor tables, --gc/--gc-full/--gc-corr, delta = unconstrained log-likelihood minus model log-likelihood, both from the .iqtree file, and parametric-bootstrap p-values with AliSim (--alisim --num-alignments 999); (3) gene, site and quartet concordance vectors with --gcf, --scfl and --concord-vector / .psi files; (4) simulating concatenated alignments built from several different trees with AliSim --alisim, -Q partition trees and --stratiestick/--tree-height. Files seen: sequences.fa, existing_tree.nex, bovidae_4K.phy, loci.treefile, bird_400, *.cf.stat, id_gcf_scf_qcf.nex."
---

# IQ-TREE: End-to-End Recipes

## Overview

Four self-contained workflows taken from the IQ-TREE recipe collection: placing new sequences onto an existing tree without disturbing it, testing whether the whole fitted model is adequate (Goldman-Cox), turning gCF/sCF/quartet factors into full concordance vectors, and simulating concatenated alignments whose blocks come from different trees. Each recipe below is **Prerequisites → Commands → Outputs → Interpretation** and uses the exact example filenames from the docs.

General notes: `iqtree` and `iqtree3` both appear in the recipes (the pages were last tested with IQ-TREE 2.2.0.3 – 2.3.x). Modern spellings are `-T` threads, `-B`/`--prefix`; the docs also show legacy `-nt`, `-bb`, `-pre`. Where an option is not in the local references it is flagged at the end.

---

## Recipe A — Phylogenetic placement of query sequences onto a reference tree

**Goal.** Add new sequences to an existing tree without altering the established relationships ("online phylogenetics": e.g. a pathogen outbreak where new genomes arrive daily), then optionally re-optimise globally.

### Prerequisites

| Item | Detail |
| --- | --- |
| Software | IQ-TREE 2 (recipe last tested with 2.2.0.3). For building the toy example yourself you also need `newick_utils` (`conda install -c bioconda newick_utils`) and AliSim (`iqtree3 --alisim`). |
| `sequences.fa` | Alignment of 20 sequences: the 17 already in the tree (`T1`–`T17`) plus three newcomers (`NEW_1`, `NEW_2`, `NEW_3`). |
| `existing_tree.nex` | Existing tree containing only those 17 sequences `T1`–`T17`. |

### Commands

```bash
# A1. Placement: topology of T1..T17 stays frozen, NEW_* get placed
iqtree -s sequences.fa -g existing_tree.nex

# A2. Optional global re-optimisation, using the placement result as starting tree
iqtree -s sequences.fa -t seqs_added.treefile -pre re_optimised
```

Reproduce the whole worked example from scratch:

```bash
# simulate an alignment of 20 sequences (seed set only for reproducibility)
iqtree3 --alisim sequences -t RANDOM{yh/20} -seed 438579088 -af fasta -m JC
grep '>' sequences.fa                       # check the names

# rename three tips as the "new" sequences
sed -i -e 's/T18/NEW_1/g' sequences.fa
sed -i -e 's/T19/NEW_2/g' sequences.fa
sed -i -e 's/T20/NEW_3/g' sequences.fa

# prune them out of the simulated tree -> the reference/constraint tree
nw_prune sequences.treefile T18 T19 T20 > existing_tree.nex

# place them back, then re-optimise without constraints
iqtree -g existing_tree.nex -s sequences.fa -pre seqs_added   -seed 438579088
iqtree -s sequences.fa     -t seqs_added.treefile -pre re_optimised -seed 438579088
```

### What happens internally

`-g` supplies a **topological constraint tree** (Newick; may be multifurcating and need not include all taxa). Taxa present in the constraint keep their topology fixed. Taxa absent from it (`NEW_1`–`NEW_3`) are first placed by **maximum parsimony**, then IQ-TREE runs its standard ML search moving only those new sequences around the frozen skeleton, keeping moves that improve likelihood. That second step matters twice: sequential parsimony placement is not always the best parsimony placement, and even the best parsimony placement can differ from the best placement under a full likelihood model.

The re-optimisation step drops `-g` entirely and passes the previous tree via `-t` (**starting** tree, free to be rearranged), so IQ-TREE improves the whole topology from your best current estimate.

### Variant: fixed reference tree with `-te`

Use `-te` instead of `-t` when the reference topology must not move at all and you only want parameters/likelihood on it:

| Option | Behaviour |
| --- | --- |
| `-t <tree>` | Starting tree for the search; topology may change. |
| `-te <tree>` | Like `-t` but **fixing** the user tree: no tree search, IQ-TREE computes the likelihood of the fixed tree (equivalently `-n 0` on a user tree). |
| `-blfix` | Also freeze the branch lengths of the tree given via `-t`/`-te`. |

Masked-reference variant. When the queries must not influence the reference at all, hand IQ-TREE a reference alignment in which the query rows are masked (gaps / `?` over every site) alongside the real query sequences, and pass the published topology with `-te` so no search can rearrange it; add `-blfix` to also freeze the reference branch lengths, leaving only the query placements and model parameters to be estimated. This masking convention is **not described in the local reference files** — verify the behaviour on a toy dataset before relying on it.

### External placers (EPA-ng / pplacer)

If the queries must be placed with an existing EPA-ng/pplacer infrastructure rather than inside IQ-TREE, the handoff is: infer and fix the reference tree/model with IQ-TREE (`-s ref.fa -m <model> -te ref.tree -pre ref`), export the reference alignment + tree in the format the placer expects, run the placer, and read its `.jplace`/placement JSON for per-query edge likelihood weights. The specific tool flags for that route are **not documented in these references** — treat them as outside this skill's verified scope.

### Scale guidance

For very large data (many tens of thousands of sequences) prefer **UShER & matOptimize** or **MAPLE**: they handle millions of sequences with tens of thousands of daily additions, trading some accuracy (matOptimize uses parsimony). IQ-TREE also integrates CMAPLE/MAPLE-style fast search directly via `--pathogen` / `--pathogen-force` (see `references/command-reference.md`, "Tree search for pathogen data").

### Outputs

`seqs_added.treefile` (placed tree), `seqs_added.iqtree`, `seqs_added.log`; `re_optimised.*` for the unconstrained pass. Inspect the printed ASCII tree: the new tips should sit inside the unchanged reference structure after A1, and may move after A2.

### Interpretation

Compare `seqs_added.treefile` and `re_optimised.treefile`: if re-optimisation changes reference branches, your earlier constrained placement was fighting the data and the added sequences deserve a full re-inference. Report placements from A1 whenever the published reference topology must be preserved.

---

## Recipe B — Goldman-Cox test of model adequacy

**Goal.** Ask whether the *entire* fitted model (tree, branch lengths, substitution parameters, everything estimated) adequately describes the data, via the "cost" delta of the model and a parametric bootstrap null distribution.

### Prerequisites

| Item | Detail |
| --- | --- |
| Software | IQ-TREE 2.2.0.3+ (AliSim needed for the simulations), `bc`, bash, and R with `ggplot2` for the plot. |
| Data | `bovidae_4K.phy` — Bovidae alignment of 5 taxa (Yak, Cow, Goat, Sheep, Antelope) and 5,000 sites, a subset of Wu et al. 2018. |
| Constraint | **This version of the test requires an alignment with no gaps and no ambiguity characters** — such sites must be removed beforehand. |
| Stationarity variant | If you want IQ-TREE's built-in Goldman-Cox test of *stationarity of base/AA frequencies* (`.gc` / `.cor` tables, `--gc` flags), read "Built-in stationarity variant" below: those options are **not in the bundled references**, so verify them against your binary first. |

### Commands

```bash
# B1. One run does it all: analyse the real data AND simulate 999 mimic alignments
iqtree -s bovidae_4K.phy --alisim simulated_MSA --num-alignments 999

# B2. delta for the observed data (define the helper once)
get_delta () {
    lnL_model=$(grep "Log-likelihood of the tree: " $1 | awk '{print $5}')
    lnL_unconstrained=$(grep "Unconstrained log-likelihood" $1 | awk '{print $5}')
    delta=$(echo $lnL_unconstrained - $lnL_model | bc)
    echo $delta
}
get_delta bovidae_4K.phy.iqtree          # -> 97.0635 (±, seed-dependent)

# B3. null distribution: analyse each simulated alignment, collect deltas
for alignment in simulated_MSA_*.phy;      do iqtree -s $alignment; done
for iqtree_file in simulated_MSA_*.phy.iqtree; do get_delta $iqtree_file >> simulated_delta.txt; done

# B4. p-value + histogram in R
Rscript -e 'simulated_deltas <- read.delim("simulated_delta.txt", header=FALSE)
observed <- 97.0635
position <- (nrow(simulated_deltas)+1) - sum(observed > simulated_deltas$V1)
print(position / (nrow(simulated_deltas)+1))'
```

The plotting snippet from the recipe:

```r
library(ggplot2)
simulated_deltas = read.delim("simulated_delta.txt", header=F)
observed = 97.0635
position = (nrow(simulated_deltas) + 1) - sum(observed>simulated_deltas$V1)
p_value  = position / (nrow(simulated_deltas) + 1)
ggplot(simulated_deltas, aes(x=V1)) + geom_histogram() +
  geom_vline(xintercept = observed, colour="red", size=1) +
  theme_minimal() + xlab("delta value") +
  ggtitle("Null distribution of delta values", subtitle="Observed value is shown as a red line")
```

### How the statistic is built

| Quantity | Source |
| --- | --- |
| `Log-likelihood of the tree: -6545.5196 (s.e. 74.4412)` | `<prefix>.iqtree` / `bovidae_4K.phy.iqtree` |
| `Unconstrained log-likelihood (without tree): -6448.4561` | same file |
| delta (the "cost" of the model) | `-6448.4561 - (-6545.5196) = 97.0635` |

Parametric-bootstrap logic: (1) measure delta on the real data; (2) declare the fitted model the null; (3) simulate many datasets under that null (`--alisim ... --num-alignments 999`) and recompute delta on each; (4) locate the observed delta in the ranked list. `p = position / (N_sim + 1)`; if the observed value exceeds all 999 simulated ones, position = 1 of 1000 and `p <= 0.001`.

### Outputs

| File | Content |
| --- | --- |
| `bovidae_4K.phy.iqtree` | Observed log-likelihoods → delta. |
| `simulated_MSA_*.phy` (+ matching `.iqtree`) | The 999 null replicates and their fits. |
| `simulated_delta.txt` | One delta per line — the null distribution. |

### Interpretation

A low cost means the model is adequate; a high cost means it is performing badly and should be rejected. With `p <= 0.001` here the full model is rejected. Expect rejection on most modern datasets — they are large and our models are simple, so rejection does **not** imply useless inference. The test cannot say *which* aspect of the model fails; it is a starting point. At alpha = 5%, reject when the observed delta falls in the highest 5% of simulated deltas.

### Built-in stationarity variant (unverified locally — check before scripting)

Recent IQ-TREE builds additionally offer a native Goldman-Cox test of **stationarity of state (base / amino-acid) frequencies** per taxon, which writes small tables instead of requiring the parametric bootstrap above. The spellings most commonly seen are `--gc` (turn it on), `--gc-full` (also test non-stationarity against the alignment-wide average, not just the tree model) and `--gc-corr` (write the frequency-correlation table), producing `<prefix>.gc` (taxon, Chi-square, df, p-value) and `<prefix>.cor` (observed vs expected frequency correlations).

**None of `--gc`, `--gc-full`, `--gc-corr`, `.gc` or `.cor` appears anywhere in the reference files bundled with this skill.** Before putting them in a pipeline, confirm them on the installed version:

```bash
iqtree3 -h 2>&1 | grep -i -B 2 -A 6 "goldman\|--gc"     # does your build have them, and how?
ls <prefix>.gc <prefix>.cor 2>/dev/null                  # after a run that enabled them
head -5 <prefix>.gc                                      # always read the file's own # header for column order
```

Interpretation, if available, follows the same logic as the manual test: a small p-value for a taxon means its base/AA composition is inconsistent with the stationary frequencies implied by the fitted model — evidence of compositional heterogeneity rather than of a wrong topology. Read the Chi-square/df/p-value rows per taxon, treat p < 0.05 as rejecting stationarity for that taxon, and respond by removing/re-composing the offending taxa, switching to a non-stationary model, or partitioning by composition. Do not quote exact column names you have not seen in your own output file.

---

## Recipe C — Gene, site and quartet concordance vectors

**Goal.** Beyond the single gCF/sCF/qCF numbers, get the full four-element concordance vector per branch: Ψ1 (concordance factor), Ψ2 (larger NNI discordance), Ψ3 (smaller NNI discordance), Ψ4 (remaining discordance; always 0 for site and quartet factors by assumption).

### Prerequisites

```bash
conda create --name concordance && conda activate concordance
conda install -c bioconda iqtree astral-tree
conda install -c conda-forge r-base r-tidyverse r-boot r-ape r-ggtext
wget https://raw.githubusercontent.com/roblanf/concordance_vectors/main/concordance_vector.R
wget https://raw.githubusercontent.com/roblanf/concordance_vectors/main/concordance_table.R
wget https://raw.githubusercontent.com/roblanf/concordance_vectors/main/change_labels.R
```

Requires IQ-TREE **2.3 or above** plus ASTRAL-III. Data: `bird_400.tar.gz`, 400 randomly selected intergenic-region alignments from Stiller et al. 2024 (up to 363 bird species per locus, single-copy orthologues).

```bash
tar -xzf bird_400.tar.gz
```

### Commands

```bash
# C1. gene trees for all 400 loci (~3.5 h on 128 threads)
iqtree3 -S bird_400 --prefix loci -T 128

# C2. species tree in ASTRAL
astral -i loci.treefile -o astral_species.tree 2> astral_species.log

# C3. quartet vectors + posterior probabilities, annotated on the species tree
astral -q astral_species.tree -i loci.treefile -t 2 \
       -o astral_species_annotated.tree 2> astral_species_annotated.log

# C4. site concordance vectors (likelihood-based sCF, 100 replicates)
iqtree3 -te astral_species_annotated.tree -p loci.best_model.nex --scfl 100 --prefix scfl -T 128

# C5. gene concordance vectors (note: chain off the scfl tree, labels accumulate)
iqtree3 -te scfl.cf.tree --gcf loci.treefile --prefix gcf -T 128

# C6. dummy run only to harvest coalescent-unit branch lengths keyed by IQ-TREE branch IDs
iqtree3 -te astral_species_annotated.tree -blfix -p loci.best_model.nex --scfl 1 --prefix coalescent_bl -T 128

# C7. assemble the vectors, relabel the tree, tabulate two branches
Rscript concordance_vector.R                 # -> concordance_vectors.csv
Rscript change_labels.R                      # -> id_gcf_scf_qcf.nex
Rscript concordance_table.R concordance_vectors.csv 598   # Palaeognathae
Rscript concordance_table.R concordance_vectors.csv 545   # Telluraves
```

Option meanings: `-S` = separate ModelFinder + tree search per alignment/partition; `-te` = fixed input tree; `--gcf` = gene concordance factors from a set of gene trees; `--scfl <n>` = likelihood-based sCF with n quartets; `-p loci.best_model.nex` = reuse the previously estimated per-locus models (big time saver); `-blfix` = keep the ASTRAL branch lengths; `--scfl 1` in C6 makes the sCFs meaningless — never use them.

### Outputs

| File | Content |
| --- | --- |
| `loci.treefile` | ML single-locus trees. |
| `loci.best_model.nex` | Every estimated model with all parameter values, per locus. |
| `loci.iqtree`, `loci.log` | Summary report / full log. |
| `astral_species.tree`, `astral_species_annotated.tree` | Species tree; same tree with `[q1=...;f1=...;pp1=...;QC=...;EN=...]` labels on every branch. |
| `gcf.cf.stat`, `scfl.cf.stat`, `coalescent_bl.cf.stat` | Tab-separated CF/DF tables (gDF1, gDF2, gDFP, counts, plus the carried-over ASTRAL labels in a `Label` column). |
| `gcf.cf.tree` | Tree with unwieldy combined annotations, e.g. `'[q1=0.5702;...;EN=364.0]'/25.4/47.0:0.0055956557`. |
| `gcf.cf.branch` | Same tree annotated with branch **IDs** matching the `.stat` tables. |
| `concordance_vectors.csv` | Final per-branch table (columns below). |
| `id_gcf_scf_qcf.nex` | Branch labels rewritten as `391-g98.54-s84.09-q98.54`. |
| `concordance_table_<ID>.csv` / `.pdf` | Vector table with 95% bootstrap CIs for one branch. |

### Concordance-vector columns and how to plot

`concordance_vectors.csv` key columns: `ID`; `gene_psi1..4` (%) and `gene_psi1_N..psi4_N` (counts) plus `gene_N`; `site_psi1..4`, `site_psi*_N`, `site_N`; `quartet_psi1..3`, `quartet_psi*_N`, `quartet_N`; `quartet_psi1_pp..psi3_pp` (ASTRAL posterior probabilities); `length_subs_per_site`; `length_coalescent`. Mapping back to IQ-TREE fields: `gene_psi1` = gCF, `gene_psi2`/`psi3` = larger/smaller of gDF1 & gDF2, `gene_psi4` = gDFP, counts from `gCF_N`, `gDF1_N`, `gDF2_N`, `gDFP_N`, `gN`; site equivalents come from `scfl.cf.stat` (`sCF`, `sDF1`, `sDF2`, `sCF_N`, `sN`); quartet values come from the ASTRAL `Label` column (`q1`,`q2`,`q3`,`f1`,`f2`,`f3`,`pp1`,`pp2`,`pp3`,`EN`). Site Ψ4/`site_psi4_N` are zero by assumption.

Plotting: open `gcf.cf.tree` (or `id_gcf_scf_qcf.nex`) in DendroScope/iTOL and load the strings as **edge labels**, then re-root and zoom freely; visualise vectors as stacked bars per branch, always showing the counts (`*_N`) alongside percentages; the recipe's tables add 95% CIs from 1000 bootstraps of the count data, which is what you need to compare values *between* branches. Counts for sites and quartets are not always integers, and rounded bootstrap CIs are unreliable for very low counts.

### Interpretation

In the label `...'/25.4/47.0:0.0055956557`: `q1=0.5702` quartet CF, `25.4` sCF, `47.0` gCF, `0.0055956557` branch length in substitutions/site. High Ψ1 with tiny Ψ2/Ψ3 (the kiwi/*Apteryx* region) means genuine concordance; a branch with `gene_psi1 = 5.60` but `gene_psi4 = 94.15` (Telluraves) means the genes disagree mostly in ways that do **not** match either NNI neighbour — a polytomy/incomplete-sorting signature rather than a clean alternative bifurcation. On the same branch the site-based signal is much stronger than the gene-based one (`site_psi1` = 45.63 vs `gene_psi1` = 5.60), so read the gene and site vectors side by side before concluding a branch is well supported. Concordance factors apply to **branches**, not nodes.

### Native one-command variant

To print the raw per-quartet vectors straight from IQ-TREE, add the concordance-vector flag to a gCF run:

```bash
iqtree3 -te astral_species_annotated.tree --gcf loci.treefile --concord-vector --prefix gcf
```

The `--concord-vector` flag and the resulting `.psi` file (one row per branch × quartet replicate) are **not documented in the reference files bundled here** — confirm the spelling with `iqtree3 -h | grep -i concord` and inspect the file header before parsing. The verified route for obtaining the same four-element vectors is the `.cf.stat` pipeline above.

---

## Recipe D — Simulating concatenated alignments that contain multiple trees

**Goal.** Produce one alignment whose blocks evolved under *different* gene trees (for benchmarking concatenation vs coalescent methods).

Two approaches: the portable AMAS-concatenation recipe (fully documented here) and AliSim's native topology-unlinked partition model (single command, but see the caveat).

### Prerequisites

| Item | Detail |
| --- | --- |
| Software | IQ-TREE 2.2 or above (for AliSim) and AMAS (<https://github.com/marekborowiec/AMAS>). |

```bash
conda create --name simaln && conda activate simaln
conda install -c bioconda iqtree amas
```

### Commands — D1: exactly two trees

```bash
# 9000 bp from one random Yule tree, 1000 bp from another
iqtree3 --alisim part1 -t RANDOM{yh/20} -m GTR+I+R8 --length 9000
iqtree3 --alisim part2 -t RANDOM{yh/20} -m GTR+I+R8 --length 1000

# concatenate in the given order
AMAS.py concat -i part1.phy part2.phy -f phylip -u phylip -t part1_part2_concat.phy -d dna
```

| Flag | Meaning |
| --- | --- |
| `--alisim <NAME>` | Activate AliSim and name the output alignment. |
| `-af fasta` | FASTA output, as written in the placement recipe. The AliSim option table documents this as `--out-format <fasta\|phy\|maple>` (default `phy`) — use `--out-format fasta` if `-af` is rejected. |
| `-t RANDOM{yh/20}` | Randomly generated Yule-Harding tree with 20 taxa. |
| `-m GTR+I+R8` | Simulation model. |
| `--length <n>` | Root-sequence length, i.e. alignment length (default 1000). |
| `AMAS.py concat -i` | Input files, in concatenation order. |
| `-f` / `-u` | Input / output alignment format (`phylip`, `fasta`). |
| `-t` | Output (concatenated) filename. |
| `-d` | Sequence type, `dna` or `aa`. |

### Commands — D2: any number of trees, any block lengths

```bash
#!/bin/bash
lengths=( 1000 1000 2000 5 3000 )              # one entry per source tree
N=$(wc -w <<< "${lengths[*]}")
to_concat=( )
for (( i=0; i<$N; i++ )); do
  len=${lengths[$i]}
  echo "index: $i, length: $len"
  iqtree3 --alisim $i -t RANDOM{yh/20} -m GTR+I+R8 --length $len
  to_concat+=( $i".phy" )
done
echo ${to_concat[*]}
AMAS.py concat -i ${to_concat[*]} -f phylip -u phylip -t $N"_concat.phy" -d dna
```

### Commands — D3: native topology-unlinked partitions (preferred when available)

Prepare a partition NEXUS and one Newick tree **per partition** (taxon sets may differ between trees):

```
#nexus
begin sets;
    charset gene_1 = DNA, 1-846;
    charset gene_2 = DNA, 847-1368;
    charpartition mine = HKY{2}+F{0.2/0.3/0.1/0.4}:gene_1{0.26019},
        GTR{1.2/0.8/0.5/1.7/1.5}+F{0.1/0.2/0.3/0.4}+G{0.5}:gene_2{1.51542};
end;
```

```
(A:0.1,(B:0.26,C:0.15):0.1,D:0.05);
(A:0.2,B:0.1,C:0.3);
```

```bash
iqtree3 --alisim multi_alignment -Q multi_genes.nex -t multi_trees.nwk
```

`-Q` selects the topology-unlinked model (one tree per partition, read in order from the multi-Newick file); with `-q` AliSim would load **only the first** tree and simulate everything under it. Partitions with missing taxa get a stretch of gaps (e.g. taxon `D` in `gene_2`). Append `{<tree_length>}` to a partition's model to rescale that gene's total tree length; use `-p` (not `-q`) if you want those per-partition rates honoured under the edge-proportional model.

### Outputs and interpretation

| File | Content |
| --- | --- |
| `part1.phy`, `part2.phy` (or `0.phy` … `N.phy`) | Per-tree simulated blocks, 20 taxa named `T1`…`T20`. |
| `part1_part2_concat.phy` / `${N}_concat.phy` | The final multi-tree concatenated alignment. |
| `multi_alignment.phy` | Single-file output of the topology-unlinked run. |

Check the concatenation length equals the sum of the block lengths. Two caveats from the recipe: random Yule trees are not guaranteed to *differ*, so if distinct topologies matter, pre-generate/inspect them or pass your own trees with `-t`; and the loop has no error checking, so a failed simulation leaves AMAS trying to concatenate a nonexistent file — add existence and length checks before using it for research. For million-site blocks add `-nt <k>`, and `-seed <n>` whenever you need reproducibility.

### Unverified in the supplied documentation

| Requested item | Status | Verified alternative in this recipe |
| --- | --- | --- |
| `--stratiestick` (and its alias `--tree-height`) for stratified multi-tree simulation with per-block tree heights | Not documented anywhere in the reference files bundled with this skill. | Per-partition tree-length scaling in the partition file (`:gene_1{0.26019}`) together with `-Q`/`-p`, or simulating separate blocks and concatenating with AMAS. |
| EPA-ng / pplacer placement flags and `.jplace` handoff | Not documented here; the local placement recipe implements placement entirely inside IQ-TREE. | Recipe A (`-g` constraint placement, then `-t` re-optimisation), or `--pathogen`/`--pathogen-force` for MAPLE-style placement of low-divergence sequences. |
| Masked reference alignment for query-only placement | The masking convention is not described in these references. | Keep the topology frozen with `-te` (plus `-blfix`) so queries cannot rearrange the reference. |
| `--concord-vector` / `.psi` output | Not documented here. | The `.cf.stat` → `concordance_vectors.csv` pipeline in Recipe C. |
| `--gc`, `--gc-full`, `--gc-corr`, `.gc`, `.cor` | Not documented here. | The parametric-bootstrap Goldman-Cox test in Recipe B. |

Before using any of the left-hand options, confirm them on the installed build, e.g. `iqtree3 -h | grep -i -A 2 "stratiestick\|tree-height"` and `head` the produced file to read its own `#` column header.

---

## Quick option lookup without loading whole files

`references/command-reference.md` is ~7.4k words across 34 headings (20 sections + 14 subsections) — pull only the section you need:

```bash
cd "$(dirname "$0")"        # or: cd /path/to/iqtree-recipes

# Full section map first (headings + line numbers):
grep -n "^## \|^### " references/command-reference.md

# Then extract just what you need, e.g.:
grep -n -A 25 "^## Tree search parameters"        references/command-reference.md
grep -n -A 25 "^## Automatic model selection"     references/command-reference.md
grep -n -A 25 "^## Specifying substitution models" references/command-reference.md
grep -n -A 22 "^## Ultrafast bootstrap parameters" references/command-reference.md
grep -n -A 18 "^## Nonparametric bootstrap"        references/command-reference.md
grep -n -A 22 "^## Tree topology tests"            references/command-reference.md
grep -n -A 18 "^## Partition model options"        references/command-reference.md
grep -n -A 20 "^## Tree search for pathogen data"  references/command-reference.md
grep -n -A 20 "^## Miscellaneous options"          references/command-reference.md
grep -n -A 12 "^## General options"                references/command-reference.md
```

Section-to-task map:

| Task | Section to grep |
| --- | --- |
| `-s`, `-st`, `-t`, `-te`, `-o`, `-pre`, `-nt`, `-seed` | `## General options` |
| ModelFinder, Lie Markov keywords, `-mset`, partition merging | `## Automatic model selection` |
| Model strings, `+F`/`+FO`/`+FQ`, `--link-exchange` | `## Specifying substitution models` |
| `+I`, `+G`, `+R`, `+ASC`, `-wsr` | `## Rate heterogeneity` |
| `-q` / `-p` / `-Q` partition models | `## Partition model options` |
| `-g` constraints, `-ninit`, `-pers`, `-nstop` | `## Tree search parameters` |
| `--pathogen`, `--sprta`, `--tree-search` | `## Tree search for pathogen data` |
| `-bb`, `-bnni`, `-bsam`, `-wbt` | `## Ultrafast bootstrap parameters` |
| `-b`, `-bc`, `-bo` | `## Nonparametric bootstrap` |
| `-au`, `-zb`, `-z`, `-n 0` | `## Tree topology tests` |
| `-blfix`, `-alninfo`, `-wsl`, `-fconst` | `## Miscellaneous options` |

## Related skills

- `iqtree` — install, datatypes, output-file overview, troubleshooting.
- `iqtree-model-selection` — ModelFinder details behind the recipe commands.
- `iqtree-complex-models` — partition/mixture model files used in Recipes C and D.
- `iqtree-tree-search` — `-g`, `-t`, `-te`, `-blfix` semantics in depth.
- `iqtree-branch-support` — UFBoot/aLRT used by the placement and GC runs.
- `iqtree-concordance-factor` — gCF/sCF theory and the `.cf.stat` format.
- `iqtree-alisim` — full AliSim option reference for Recipe B and Recipe D.
- `iqtree-dating` — LSD2 downstream of a placed/re-optimised tree.
- `iqtree-rooting`, `iqtree-poms`, `iqtree-command-reference`.

## Resources

- `references/phylogenetic-placement.md` — Recipe A in full: `sequences.fa`/`existing_tree.nex`, `-g` placement mechanics, `-t` global re-optimisation, the from-scratch AliSim + `nw_prune` script, UShER/matOptimize/MAPLE scale advice.
- `references/goldman-cox-test.md` — Recipe B: gap-free-input requirement, the single `--alisim --num-alignments 999` command, `get_delta()` helper, the exact `.iqtree` lines to parse, ranking logic and R plotting code.
- `references/concordance-vector.md` — Recipe C: environment setup, `-S` gene trees, ASTRAL `-q`/`-t 2` annotation, the three chained IQ-TREE commands, every `.cf.stat`/`.treefile` output, the complete `concordance_vectors.csv` column dictionary and CI methodology.
- `references/Simulating-concatenated-alignments-with-multiple-trees.md` — Recipe D: two-tree AMAS workflow, the generalised bash loop, and the author's caveats on random-tree distinctness and error checking. Note the page itself redirects readers to "Topology-unlinked partition model" in the AliSim docs (Recipe D3 above).
- `references/command-reference.md` — 34-heading IQ-TREE option reference (use the greps above). Load only the section you need; never read all ~7.4k words.
- `references/index.md` — the recipes index page listing the four official recipes (20 lines; safe to read whole).
- `../iqtree-alisim/references/alisim.md` (sibling skill folder, read-only) — AliSim option table and the topology-unlinked partition model section used in Recipe D3.
