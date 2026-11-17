---
name: iqtree-branch-support
description: Assess and interpret branch support / statistical confidence on IQ-TREE phylogenies. Use when asking how confident a clade is, running or tuning ultrafast bootstrap (UFBoot, -B/-bb), standard nonparametric bootstrap (-b/-bo/-bc), single-branch tests SH-aLRT (-alrt), aBayes (-abayes), local bootstrap (-lbp), -bnni model-violation correction, or reading support values out of .treefile, .iqtree, .contree, .splits.nex, .ufboot, .suptree, .symtest.csv, .lmap.svg/.lmap.eps. Also covers likelihood mapping (-lmap), matched-pairs tests of symmetry (--symtest*), topology tests used to check supports (-z, -zb, -au, -zw, bp-RELL), the SH-aLRT >=80 / UFBoot >=95 reliability rule, and diagnosing low or suspiciously-high supports.
---

# Assessing branch support in IQ-TREE

## Overview

IQ-TREE offers three tiers of confidence measures for an internal branch: resampling-based
supports (ultrafast bootstrap or standard nonparametric bootstrap), cheap single-branch tests
(SH-aLRT, aBayes, local BP), and whole-tree topology tests (`-z`/AU). Pick the tier by dataset
size and by how much model violation you suspect, then interpret the number with the correct
probabilistic meaning — UFBoot% and bootstrap% are **not** interchangeable.

Run `iqtree -h` (or `iqtree3 -h`) for the authoritative option list; every flag below was taken
from the mirrored docs at the line numbers given in "Sources".

## Version note: `-B` vs `-bb`

The tutorials use the modern long-style flags; the Command Reference table still lists the 1.x
short names. Both refer to the same feature:

| Modern (docs' tutorial usage) | 1.x / Command-Reference table name | Meaning |
| --- | --- | --- |
| `-B 1000` | `-bb 1000` | Number of UFBoot replicates (`Tutorial.md:291`, `Command-Reference.md:492`) |
| `--sampling GENE` | `-bsam GENE` | Partition-level resampling (`Advanced-Tutorial.md:176-177`) |
| `--sampling GENESITE` | `-bsam GENESITE` | Partitions, then sites within partitions (`Advanced-Tutorial.md:183-184`) |
| `--prefix` | `-pre` | Output prefix (`Tutorial.md:152-153`) |
| `-T AUTO` | `-nt AUTO` | Thread count (`Tutorial.md:397-402`) |
| `--rate` | `-wsr` | Site rates (`Advanced-Tutorial.md:456`) |
| `--symtest-only` | (IQ-TREE 3 only) | Tests of symmetry (`assessing-assumptions.md:27`) |

Quote the 1.x equivalent as `-B` (for `-bb`) when writing commands; the docs never define `-B`
in an options table, only in example command lines.

## Choosing a support measure

| Goal | Command core | Documented minimum |
| --- | --- | --- |
| Routine branch supports (default recommendation) | `-B 1000 -alrt 1000` | 1000 replicates each (`Tutorial.md:295`, `:332`) |
| Severe model violation suspected | `-B 1000 -alrt 1000 -bnni` | ~2x slower (`Tutorial.md:305-310`) |
| Classical, conservative bootstrap | `-b 100` | >=100 replicates (`Command-Reference.md:522`, `Tutorial.md:322`) |
| Bootstrap trees only, no consensus/original-alignment analysis | `-bo 100` | (`Command-Reference.md:524`) |
| Bootstrap only, skip original-alignment run | `-bc 100` | (`Command-Reference.md:523`) |
| Extremely large data (>10,000 taxa) — single-branch tests only | `-alrt 1000 -abayes -lbp 1000` | ("recommended for extremely large data sets", `Command-Reference.md:528`) |
| Parametric aLRT instead of SH-aLRT | `-alrt 0` | (`Command-Reference.md:532`) |
| Concordance in a phylogenomic concatenation | hand off to `iqtree-concordance-factor` | FAQ says UFBoot/BS "tend to be 100%" for concatenated trees (`Frequently-Asked-Questions.md:62`) |

All of `-alrt`, `-abayes`, `-lbp` and UFBoot can be combined in one run; each branch then gets
several values separated by `/`, and the order is stated in the `.iqtree` report
(`Command-Reference.md:536`).

## Ultrafast bootstrap (UFBoot)

Enable it:

```bash
iqtree -s example.phy -m TIM2+I+G -B 1000
# for version 1.x change -B to -bb
```

How many replicates: 1000 is the documented minimum recommended number (`Tutorial.md:295`;
`-bb`: "Specify number of bootstrap replicates (>=1000)", `Command-Reference.md:492`). The
partition-model examples go up to 5000 (`Command-Reference.md:513`). No other replicate count is
recommended anywhere in the docs.

### Automatic stop rules and tuning knobs

UFBoot can finish before reaching `-B` replicates once it judges the supports converged. The docs
expose that machinery only through these options (`Command-Reference.md:492-500`):

| Option | Documented meaning / default |
| --- | --- |
| `-bcor` | Minimum correlation coefficient for the UFBoot convergence criterion. DEFAULT: 0.99 |
| `-nstep` | Iteration interval at which UFBoot convergence is checked. DEFAULT: every 100 iterations |
| `-nm` | Maximum number of iterations to stop. DEFAULT: 1000 |
| `-beps` | Small epsilon to break ties in RELL evaluation for bootstrap trees. DEFAULT: 0.5 |
| `-bnni` | Hill-climbing NNI optimisation of each bootstrap tree on its own bootstrap alignment; recommended under severe model violation; costs ~2x time but reduces overestimated supports (since IQ-TREE 1.6) |
| `-bsam` / `--sampling` | Resampling strategy for partitions: sites within partitions (default), `GENE` (resample partitions), `GENESITE` (partitions then sites) |
| `-wbt` | Write the bootstrap trees to the `.ufboot` file. DEFAULT: OFF |
| `-wbtl` | Like `-wbt` but bootstrap trees written **with branch lengths**. DEFAULT: OFF |

Not documented in this mirror: the explicit stop rule from the UFBoot2 paper (transfer-distance /
standard-deviation criteria), and any `--min-uniq-req`, `--ublo`, `--cbn`, `--bcmp`, `--bsup`,
`--bpol` style options. Those strings appear nowhere here — do not emit them.

### UFBoot with partitions

Sites are resampled *within* partitions by default; each partition's replicates are generated
separately and concatenated (`Advanced-Tutorial.md:162-171`). The same holds for the standard
bootstrap. To reduce false positives (a wrong branch getting 100%):

```bash
iqtree -s example.phy -p example.nex -B 1000 --sampling GENESITE
```

### Mapping supports onto your own tree

```bash
iqtree -sup input_tree set_of_trees      # -> .suptree  (Advanced-Tutorial.md:411)
iqtree -con mytrees                      # extended majority-rule consensus -> .contree
iqtree -con mytrees -minsup 0.5          # majority-rule consensus
iqtree -con mytrees -minsup 0.5 -bi 100  # discard first 100 trees (burn-in)
iqtree -net mytrees                      # consensus network -> .nex
```

`-minsup` DEFAULT is 0 = extended majority-rule consensus (`Command-Reference.md:619`).
`-suptag ALL` writes, per node, the IDs of trees in which that node appears
(`Command-Reference.md:622`). `-scale` (default 100, i.e. percentages) and `-prec` (default 0)
control the printed support values of `-sup` (`Command-Reference.md:623-624`).

## Standard nonparametric bootstrap

```bash
iqtree -s example.phy -m TIM2+I+G -b 100
```

`-b` runs the bootstrap *and* the analysis on the original alignment and provides a consensus
tree; output files are the same set as for UFBoot (`Tutorial.md:314-322`,
`Command-Reference.md:516-524`). It is slow — "may take weeks/months for large data sets"
(`Frequently-Asked-Questions.md:151`).

Splitting replicates across machines (documented workflow, `Frequently-Asked-Questions.md:149-199`):

```bash
# 1. N independent partial runs, distinct prefixes
iqtree -nt 4 -s input_alignment -bo 20 ... -pre boot1
# 2. merge the bootstrap tree files
cat boot*.boottrees > alltrees
# 3. consensus from the merged trees
iqtree -con -t alltrees
# 4. re-optimise branch lengths of the consensus on the original alignment
iqtree -s input_alignment -te alltrees.contree -pre alltrees.contree
# 5. alternatively map supports onto your ML tree
iqtree -sup input_alignment.treefile -t alltrees   # -> .suptree
```

Related: `-czb` collapses near-zero branches so a polytomy does not inflate short-branch
bootstrap supports (`Command-Reference.md:700`). Identical sequences are removed during analysis
and re-added at the end precisely so that their bootstrap value can legitimately be <100%
(`Frequently-Asked-Questions.md:221-225`).

## Single-branch tests

```bash
iqtree -s example.phy -m TIM2+I+G -alrt 1000 -B 1000
# for version 1.x change -B to -bb
```

| Option | Test | Documented notes |
| --- | --- | --- |
| `-alrt <n>` | SH-like approximate likelihood ratio test (SH-aLRT) | `n >= 1000` recommended (`Command-Reference.md:532`, `Tutorial.md:332`); `-alrt 0` switches to the parametric aLRT |
| `-abayes` | Approximate Bayes (aBayes) test | no replicate count needed (`Command-Reference.md:533`) |
| `-lbp <n>` | Fast local bootstrap probability | `n >= 1000` (`Command-Reference.md:534`) |
| `--alrt <n>` | SH-aLRT under the CMAPLE pathogen search | only valid with `--pathogen`/`--pathogen-force` (`Command-Reference.md:437`, `:462-468`) |
| `--sprta` | SPRTA supports, incl. external branches, in `.nex` | pathogen/CMAPLE only; also `--sprta-zero-branch`, `--sprta-other-places` (`Command-Reference.md:438-440`, `:470-476`) |

These tests are faster than any bootstrap and are recommended for extremely large datasets
(`Command-Reference.md:528`).

## Reading the numbers

Supports appear in three places:

1. `.treefile` — values on branch labels, readable by FigTree/Dendroscope/ETE (`Tutorial.md:343`).
2. `.iqtree`, section `MAXIMUM LIKELIHOOD TREE` — ASCII tree, percentages in parentheses; the
   header line states the order, e.g. `Numbers in parentheses are SH-aLRT support (%) /
   ultrafast bootstrap support (%)` (`Tutorial.md:295`, `:346-347`).
3. `.contree` (consensus with supports, branch lengths optimised on the original alignment),
   `.splits.nex` (percentage for **all** splits, including conflicting ones that never reach the
   consensus — open in SplitsTree), and `.splits` (same data, star-dot format, needs `-wsplits`)
   `.splits` (same data, star-dot format, needs `-wsplits`) (`Tutorial.md:297-299`).

### Interpretation thresholds actually given in the docs

| Statement | Source |
| --- | --- |
| "95% support correspond roughly to a probability of 95% that a clade is true" | `Frequently-Asked-Questions.md:58` |
| "For UFBoot, you should only start to rely on a branch if its support is >= 95%" | `Frequently-Asked-Questions.md:58` |
| Do **not** compare BS% with UFBoot% directly; standard bootstrap is more conservative | `Frequently-Asked-Questions.md:58` |
| "One would be more confident if a clade has its SH-aLRT >= 80% and UFboot >= 95%" | `Frequently-Asked-Questions.md:60` |
| These recommendations apply **only to single-gene trees**; for concatenated phylogenomic trees UFBoot (and even Felsenstein's BS) tend to 100% — compute concordance factors instead | `Frequently-Asked-Questions.md:62` |
| KH/SH/AU return p-values: reject a tree when p < 0.05 (marked `-`); `+` marks the 95% confidence set | `Advanced-Tutorial.md:293`, `:373-374` |
| bp-RELL and c-ELW are posterior **weights**, not p-values; they sum to 1 over the tested trees | `Advanced-Tutorial.md:295` |
| A good dataset has many fully-resolved quartets and few unresolved quartets | `assessing-assumptions.md:102` |

No numeric threshold for `-bcor`, aBayes, `-lbp`, standard bootstrap, or the `-minsup` default
beyond `0.5` = majority rule is given anywhere in this mirror.

## Topology tests as a support cross-check (brief)

Use `-z` when a branch's support looks questionable and you want to know whether the alternative
topology is statistically excluded:

```bash
iqtree -s example.phy -z example.treels -n 0 -zb 1000 -zw -au
```

| Option | Role |
| --- | --- |
| `-z <file>` | Set of trees (NEWICK) whose log-likelihoods are computed |
| `-zb <n>` | Number of RELL replicates, `n >= 1000` recommended minimum |
| `-zw` | Additionally weighted KH and weighted SH |
| `-au` | Additionally the approximately unbiased (AU) test (needs `-zb`) |
| `-n 0` | Skip full tree search; parameters estimated on an initial parsimony tree |
| `-te <tree>` | Like `-n 0` but estimate parameters on a user tree |

(`Command-Reference.md:581-608`, `Advanced-Tutorial.md:251-301`.) Without `-n 0` IQ-TREE first
reconstructs an ML tree, then scores the `-z` trees with the ML tree's estimated parameters; the
results land in the `USER TREES` section of `.iqtree`, and trees with optimized branch lengths go
to `<prefix>.treels.trees` (`Advanced-Tutorial.md:259-265`). IQ-TREE auto-drops duplicate
topologies (`Advanced-Tutorial.md:289`). Prefer AU over KH/SH: KH has no multiple-testing
correction, SH becomes too conservative with many trees (`Advanced-Tutorial.md:297-299`).

For the full treatment of `-z`, `-wsl`/CONSEL-style site likelihoods and constrained-tree
comparisons, hand off to the **iqtree-tree-search** skill.

## Testing assumptions before trusting supports

### Tests of symmetry (stationarity + homogeneity)

```bash
iqtree3 -s example.phy -p example.nex --symtest-only        # -> <prefix>.symtest.csv
iqtree3 -s example.phy -p example.nex --symtest-remove-bad  # drop SymPval < 0.05 partitions
```

| Column | Small p-value means |
| --- | --- |
| `SymPval` | Stationarity and/or homogeneity rejected for that partition |
| `MarPval` | Stationarity rejected |
| `IntPval` | Homogeneity rejected |

(`assessing-assumptions.md:22-59`.) Extra options, verbatim from `iqtree3 -h`
(`assessing-assumptions.md:72-80`): `--symtest`, `--symtest-only`, `--symtest-remove-bad`,
`--symtest-remove-good`, `--symtest-type MAR|INT`, `--symtest-pval NUMER` (default 0.05),
`--symtest-keep-zero`. After removing bad partitions, compare the "all" vs "good-only" trees with
topology tests. Note there is no `-g`-style `--symtest` shorthand in the docs; the flag is
`--symtest*` (IQ-TREE 3), while `-g` exists separately as the constraint-tree option
(`Command-Reference.md:403`).

Separately, every run starts with a composition chi-square test per sequence; a sequence is
`failed` when its composition deviates from the alignment average. Treat it as an *explorative
tool*, not an automatic removal rule — but if the tree shows an unexpected topology it can point
at the offending taxon or partition (`Frequently-Asked-Questions.md:119-139`).

### Likelihood mapping (treelikeness / phylogenetic signal)

```bash
iqtree -s example.phy -lmap 2000 -n 0 -m TEST
```

Outputs: `LIKELIHOOD MAPPING STATISTICS` in `.iqtree`, plus plots `example.phy.lmap.svg` and
`example.phy.lmap.eps` (`Command-Reference.md:120-122`, `assessing-assumptions.md:82-104`).

| Option | Meaning |
| --- | --- |
| `-lmap <n>` | Number of random quartets; `-lmap ALL` draws all unique quartets |
| `-lmclust <nex>` | NEXUS `taxset` file with 1-4 clusters (plus optional `IGNORED`) for 2-, 3- or 4-cluster mapping |
| `-n 0` | Stop after mapping, skip tree search |
| `-wql` | Write quartet log-likelihoods to `.lmap.quartetlh` |

Draw at least **25x the number of sequences** so each sequence is covered ~100 times
(`Command-Reference.md:133`). In the 7-area plot, corners 1-3 = fully resolved (informative),
rectangles 4-6 = partly informative, centre 7 = uninformative (`assessing-assumptions.md:102`).

## Diagnosing support-value problems

| Symptom | First checks (all documented) |
| --- | --- |
| Many low SH-aLRT/UFBoot values | Run likelihood mapping: high centre-area (unresolved) share means little signal, not a wrong method (`-lmap`). Then `-B 5000` for finer resolution, and consider more data. |
| High UFBoot but conflicting loci/sites | Concatenation inflates supports (`Frequently-Asked-Questions.md:62`); compute gCF/sCF → `iqtree-concordance-factor`. |
| Suspiciously uniform 100% supports in a partitioned supermatrix | Try `--sampling GENESITE`, which "may help to reduce false positives (i.e. wrong branch receiving 100% support)" (`Advanced-Tutorial.md:180`) |
| Suspected severe model violation | Add `-bnni` (`Tutorial.md:303-312`); check `--symtest*` p-values and the composition-test `failed` list; consider profile mixtures (PMSF "can effectively ameliorate long branch attraction artefacts", `Complex-Models.md:268`) → `iqtree-complex-models` |
| Long-branch-attraction-like grouping | Docs give no dedicated LBA test. Use the routes above (mixtures, PMSF, `--symtest`, AU test of the alternative topology), plus the composition test to locate the aberrant taxon (`Frequently-Asked-Questions.md:133`) |
| Support unstable between reruns | Fix `-seed` (`Command-Reference.md:74`); check tree-search adequacy with `-nstop 500 -pers 0.2` and >=10 repeats (`Command-Reference.md:405-412`) → `iqtree-tree-search` |
| Very short internal branch with high support | `-czb` to collapse near-zero branches (`Command-Reference.md:700`); inspect the second-best conflicting split in `.splits.nex` (`Tutorial.md:298`) |
| Low supports concentrated on one clade | Score the competing topologies: `-z <trees> -n 0 -zb 1000 -au`; a `+` in the 95% confidence set means undecidable, matching weak SH-aLRT/UFBoot (`Advanced-Tutorial.md:293`, `:378`) |

## Worked example 1 — single gene, both supports

Input: `example.phy` (mitochondrial DNA of 27 animals), best-fit model already known to be
`TIM2+I+G`.

```bash
iqtree -s example.phy -m TIM2+I+G -alrt 1000 -B 1000 -T AUTO
```

Relevant excerpt of `example.phy.iqtree` (`Tutorial.md:345-382`):

```
NOTE: Tree is UNROOTED although outgroup taxon 'LngfishAu' is drawn at root
Numbers in parentheses are SH-aLRT support (%) / ultrafast bootstrap support (%)
       |                  +--| (85/72)
       |                  |  |    +------------------------Crocodile
       |                  |  +----| (96.5/97)
       |               +--| (39/51)
       |               |  +---------------------------Sphenodon
                 |                |  +--| (68.3/75)
                 |                |     |  +-----Cow
```

Reading it: `Sphenodon` (39/51) and `Turtle` (85/72) fail the SH-aLRT >= 80 + UFBoot >= 95 rule,
as does the position of `Seal` (68.3/75) — reptile relationships and the seal branch are not
reliable. `(96.5/97)` for Crocodile+Bird passes both. Branches at 100/100 are solid. Because this
is a single-gene tree, the >=80/>=95 rule applies; had this been a concatenated phylogenomic
matrix, the same numbers would need concordance factors instead.

## Worked example 2 — is the weakly supported alternative excluded?

Continuing from example 1: the `Seal` branch was ambiguous, and a constrained search putting
`(Human,Seal)` as sisters is claimed to be plausible. Build unconstrained + constrained trees,
then test them.

```bash
# unconstrained
iqtree -s example.phy -m TIM2+I+G --prefix example.unconstr
# constrained: file example.constr1 contains ((Human,Seal),(Cow,Whale),Mouse);
iqtree -s example.phy -m TIM2+I+G -g example.constr1 --prefix example.constr1
# a clearly wrong constraint
iqtree -s example.phy -m TIM2+I+G -g example.constr3 --prefix example.constr3
# pool topologies and test
cat example.unconstr.treefile example.constr1.treefile example.constr3.treefile > example.treels
iqtree -s example.phy -m TIM2+I+G -z example.treels -n 0 -zb 1000 -au
```

`USER TREES` output (`Advanced-Tutorial.md:359-375`):

```
Tree      logL    deltaL  bp-RELL    p-KH     p-SH    c-ELW     p-AU
-------------------------------------------------------------------------
  1   -21152.617   0.000  0.7110 + 0.7400 + 1.0000 + 0.6954 + 0.7939 +
  2   -21156.802   4.185  0.2220 + 0.2600 + 0.5910 + 0.2288 + 0.3079 +
  4   -21339.596 186.980  0.0000 - 0.0000 - 0.0000 - 0.0000 - 0.0000 -
```

Reading it: trees 1 and 2 stay in the 95% confidence set (`+` under `p-AU`), so the seal grouping
is **undecidable** — consistent with its 68.3/75 SH-aLRT/UFBoot supports, not contradicted by
them. Tree 4 is significantly excluded (`-`, `deltaL` 187). `bp-RELL`/`c-ELW` are weights summing
to 1, not p-values, so never quote them as confidence percentages.

## Worked example 3 — supports you may not be able to trust

Partitioned phylogenomic run where UFBoot came back ~100% everywhere:

```bash
# 1. how much real signal is there?
iqtree -s example.phy -p example.nex -lmap 2000 -n 0 -m TEST
# 2. which partitions violate stationarity/homogeneity?
iqtree3 -s example.phy -p example.nex --symtest-only
# 3. re-run excluding partitions with SymPval < 0.05, with model-violation-safe UFBoot
iqtree3 -s example.phy -p example.nex --symtest-remove-bad -B 1000 -alrt 1000 -bnni --sampling GENESITE
```

Read `example.nex.symtest.csv` (`SymPval`/`MarPval`/`IntPval` per partition), the
`LIKELIHOOD MAPPING STATISTICS` block and `.lmap.svg` (many centre-quartets ⇒ weak treelikeness),
then compare the good-partitions-only tree against the all-partitions tree with
`-z ... -n 0 -zb 1000 -au`. If the concatenation still gives blanket 100%s, the supports are an
artefact of concatenation — switch to gCF/sCF (`iqtree-concordance-factor`).

## Resources

| File | What's in it | Load when |
| --- | --- | --- |
| `references/advanced-tutorial.md` | Partition UFBoot and `--sampling GENE/GENESITE`, constrained search `-g`, full topology-test workflow (`-z`, `-zb`, `-zw`, `-au`) incl. a `USER TREES` table, consensus/support mapping (`-con`, `-net`, `-sup`, `-minsup`, `-bi`), site rates, `--robust-phy` trimming | Building a multi-step command sequence for supports or topology comparisons |
| `references/assessing-assumptions.md` | Matched-pairs tests of symmetry (`.symtest.csv` columns, the `iqtree3 -h` block of `--symtest*` flags) and likelihood mapping (quartet areas, interpretation) | Deciding whether supports can be trusted at all |

Never read these files whole; pull only the needed section:

```bash
grep -n "^## " references/advanced-tutorial.md            # map sections first
grep -n "^## " references/assessing-assumptions.md
grep -n -A 25 "^## Ultrafast bootstrapping with partition model" references/advanced-tutorial.md
grep -n -A 50 "^## Tree topology tests" references/advanced-tutorial.md
grep -n -A 50 "^## Testing constrained tree" references/advanced-tutorial.md
grep -n -A 33 "^## Consensus construction and bootstrap value assignment" references/advanced-tutorial.md
grep -n -A 58 "^## Tests of symmetry" references/assessing-assumptions.md
grep -n -A 30 "^## Likelihood mapping" references/assessing-assumptions.md
```

Option tables live in the wider mirror (read-only, outside this skill folder):
`/data/zengzichao/软件技能制作/iqtree-doc/markdown/Command-Reference.md` sections
`Ultrafast bootstrap parameters` (line 486), `Nonparametric bootstrap` (516), `Single branch
tests` (526), `Tree topology tests` (581), `Constructing consensus tree` (610), `Likelihood
mapping analysis` (120); plus `Tutorial.md` lines 282-341 and
`Frequently-Asked-Questions.md` lines 56-62 and 149-199.

## Related skills

- `iqtree-tree-search` — full treatment of tree search quality, `-z` topology tests, `-wsl` site
  log-likelihoods, constraint searches; supports are only trustworthy after a converged search.
- `iqtree-model-selection` — the best-fit model (`-m MFP`/`TEST`) behind every support value;
  model misspecification is the root cause of inflated UFBoot.
- `iqtree-complex-models` — partitions, profile mixtures and PMSF when model violation (e.g. LBA
  artefacts) drives unreliable supports.
- `iqtree-concordance-factor` — gCF/sCF for phylogenomic concatenations, where the docs say
  UFBoot and standard bootstrap supports "tend to be 100%".
- `iqtree-rooting` — rootstrap support and non-reversible models for rooted-tree confidence.
- `iqtree-command-reference` — exhaustive flag lookup when a support option is not listed here.
- `iqtree` — entry point, installation, checkpointing (`-redo`), output-file overview.
