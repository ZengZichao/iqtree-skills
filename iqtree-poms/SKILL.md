---
name: iqtree-poms
description: Run IQ-TREE polymorphism-aware models (PoMo / POM) on population-level site-frequency data. Use when building or supplying an allele-count "counts file" (.cf, COUNTSFILE header with A,C,G,T counts per site), converting a polymorphic alignment or VCF into counts format, asking which PoMo substitution model to name after -m (HKY+P, GTR+P, +N virtual population size, +WB/+WH/+S sampling), how to fix the level of polymorphism Watterson's theta 4Nu with +P{EMP} or +P{0.0025}, how to get bootstrap support for PoMo trees (-b / -bb / -B), how to read PoMo branch lengths (mutations and frequency shifts per site, divide by N^2 for substitutions per site), or how to handle many identical or polimorphic sequences per population. Triggers on .cf files, cflib, great-ape example data, and species-tree inference with several chromosomes per species.
---

# IQ-TREE: Polymorphism-Aware Models (PoMo / POM)

## Overview

PoMo builds on top of DNA substitution models and describes evolution of populations through a **virtual population** of constant size `N`, so it takes allele-frequency (site-frequency-spectrum) information from multiple chromosomes per population and naturally accounts for incomplete lineage sorting. This skill covers the counts-file input, the `+P` model-string grammar (base model, virtual population size, polymorphism level, sampling method, frequencies, rate heterogeneity), bootstrapping, and branch-length interpretation.

You need **more sequences per species/population — ideally five or more chromosomes each** — otherwise there is no site-frequency information and PoMo brings nothing over a standard DNA analysis.

## Prerequisites and version check

| Item | Detail |
| --- | --- |
| Binary | The PoMo page documents the beta **1.6** line and invokes `iqtree` (not `iqtree3`). Check your build first: `iqtree -h` and scroll to the heading `POLYMORPHISM AWARE MODELS (PoMo)` — that section lists every PoMo option *and* every available PoMo model in your binary. |
| Input | Either a counts file (`.cf`) or a normal polymorphic DNA alignment. |
| Reference data | `example.cf` ships with the download; it is a subset of the great ape data analysed in Schrempf et al. (2016), with twelve great-ape populations, one to 23 individuals each (two to 56 chromosomes). |

Citation when publishing: Dominik Schrempf, Bui Quang Minh, Nicola De Maio, Arndt von Haeseler, Carolin Kosiol (2016) *Reversible polymorphism-aware phylogenetic models and their application to tree inference.* J. Theor. Biol. 407, 362-370.

## Counts files (allele counts)

Use a counts file when populations have many individuals: counting bases per position instead of listing every chromosome makes the file much smaller and faster to parse.

Structure, exactly three things:

1. One header line marking the file as a counts file and stating the number of populations (= leaves on the tree) and the number of sites, separated by white space: `COUNTSFILE  NPOP 5   NSITES N`.
2. A second header line with white-space-separated headings: `CHROM`, `POS`, then one column per sequence/population name.
3. Many lines of `A,C,G,T` counts plus their chromosome and position.

Example:

```
COUNTSFILE  NPOP 5   NSITES N
CHROM  POS  Sheep    BlackSheep  RedSheep  Wolf     RedWolf
1      1    0,0,1,0  0,0,1,0     0,0,1,0   0,0,5,0  0,0,0,1
1      2    0,0,0,1  0,0,0,1     0,0,0,1   0,0,0,5  0,0,0,1
9      8373 0,0,0,1  1,0,0,0     0,1,0,0   0,1,4,0  0,0,1,0
Y      9999 0,0,0,1  0,1,0,0     0,1,0,0   0,5,0,0  0,0,1,0
```

Notes:

- Lines before the first header line that start with `#` are treated as comments.
- Each count quadruple sums to the number of chromosomes sampled for that population at that site (e.g. `0,1,4,0` = 5 chromosomes).
- The four numbers are counts of **A, C, G, T** in that order.

### Building a counts file

| Route | Status in the local docs |
| --- | --- |
| `cflib` Python library (recommended by the docs) | Documented. See the GitHub repository of the counts-file library `cflib` (<https://github.com/pomo-dev/cflib>) for instructions; the great-ape source data live at <https://github.com/pomo-dev/data/tree/master/SystBiol2015>. |
| Your own script | Documented — the page explicitly allows "create counts files with your own scripts". |
| IQ-TREE built-in converters (`--pop-file`, `--pop-strategy`, the `POPCONV` tool) | **Not documented in these references.** Verify against your own binary with `iqtree -h` (section `POLYMORPHISM AWARE MODELS (PoMo)`) before using; do not assume option spellings. |

Same caveat for a direct VCF-to-counts path: the reference material here does not describe VCF import, so go VCF → per-population alignments → counts (via `cflib`) unless your binary advertises otherwise.

## First run

```bash
iqtree -s example.cf -m HKY+P            # single core
iqtree -nt 4 -s example.cf -m HKY+P      # multicore (OMP) version
iqtree -s example.cf -m HKY+P -pre myprefix   # custom output prefix
```

`-s` takes the alignment/counts file, `-m` takes the model. `HKY+P` means the HKY substitution model **with polymorphisms** (i.e. PoMo); `+P` is what switches PoMo on. Everything else works like a standard IQ-TREE run.

Output files (same set as the standard version):

| File | Content |
| --- | --- |
| `example.cf.iqtree` | Main self-readable report, including a textual representation of the final tree and the tree-length block (see "Branch lengths"). |
| `example.cf.treefile` | ML tree in Newick; visualise in FigTree/iTOL. |
| `example.cf.log` | Full log of the run. Send this plus the alignment when reporting bugs. |

Default prefix is the alignment file name; change it with `-pre` (documented for PoMo; `--prefix` is the modern IQ-TREE 2/3 spelling) so repeated analyses on one file do not overwrite each other.

## Model strings: the `+P` grammar

Compose the model as `<DNA-model>+P[<polymorphism>][+<sampling>][+<freq>][+G<k>]` and pass it to `-m`:

| Piece | Values | Meaning / default |
| --- | --- | --- |
| DNA base model | `HKY`, `GTR`, ... (the two the PoMo page names; full list follows the standard DNA models — confirm with `iqtree -h`) | Underlying nucleotide substitution model. Example: `-m GTR+P`. Rate parameters can be fixed, e.g. `HKY{6.0}` sets the transition/transversion ratio to 6.0. |
| `+P` | bare | PoMo with the level of polymorphism **inferred** during likelihood maximisation (default behaviour for population data). |
| `+P{EMP}` | `EMP` | Fix the level of polymorphism to the observed empirical value in the data. |
| `+P{<theta>}` | e.g. `+P{0.0025}` | Set the level of polymorphism by hand. This quantity is Watterson's theta, also called `4Nu`. |
| `+N<k>` | odd `3..19`, plus `2` and `10` | Virtual population size `N`. Default **9**, described as a good and stable value. |
| `+WB` | | Weighted **binomial** sampling — the **default**. |
| `+WH` | | Weighted **hypergeometric** sampling. |
| `+S` | | **Sampled**: randomly draw N samples with replacement; the picked sample constitutes one PoMo state given likelihood 1, all other PoMo states get likelihood 0. |
| `+F` / `+FO` | | State frequencies. Default is empirical (observed A/C/G/T frequencies in the data) = `+F`; `+FO` estimates allele frequencies together with the rate parameters during maximisation. |
| `+G4` | `+G<k>` | Discrete Gamma rate heterogeneity. As of v1.6 **only discrete Gamma** is supported for PoMo. |

Runnable examples:

```bash
iqtree -s example.cf -m GTR+P                 # GTR + PoMo
iqtree -s example.cf -m HKY+P+N19             # raise virtual population size to 19
iqtree -s example.cf -m HKY+P{EMP}            # polymorphism fixed to empirical value
iqtree -s example.cf -m HKY+P{0.0025}         # polymorphism set by hand
iqtree -s example.cf -m HKY{6.0}+P{0.0025}    # fully specified model -> tree search only
iqtree -s example.cf -m HKY+P+WB              # weighted binomial (default)
iqtree -s example.cf -m HKY+P+WH              # weighted hypergeometric
iqtree -s example.cf -m HKY+P+S               # sampled
iqtree -s example.cf -m GTR+P+FO              # optimise allele frequencies
iqtree -s example.cf -m HKY+P+G4              # 4 Gamma rate categories
```

Choosing among the sampling methods:

| Method | When |
| --- | --- |
| `+WB` (default) | Most stable in the authors' experience. Slight **overestimation of heterozygosity** is expected: monomorphic (fixed) states can be reached from polymorphic states during sampling, while polymorphic states cannot be reached from monomorphic ones (sampling involves no mutation), so heterozygosity at the leaves gets inflated to match the data. |
| `+WH` | Use when you want to avoid that overestimation effect. |
| `+S` | Simplest approximation; a single randomly drawn sample defines the state. |

### Choosing N (virtual population size)

| Data situation | Recommended `N` |
| --- | --- |
| Typical | Keep the default 9. |
| Only very few chromosomes per population (e.g. two to four) | Lower `N` to the average number of chromosomes per population. |
| Plenty of data and runtime acceptable | Increase up to a maximum of 19. |
| Allowed values | Odd values 3 to 19, plus 2 and 10. |

More drift relative to mutation as `N` grows, which changes branch lengths (below).

## Branch support

Both bootstrap schemes work as usual with PoMo:

```bash
iqtree -s example.cf -m HKY+P -b 100        # standard non-parametric bootstrap, 100 replicates
iqtree -s example.cf -m HKY+P -bb 1000      # ultrafast bootstrap (UFBoot), 1000 replicates
```

| Option | Meaning |
| --- | --- |
| `-b <n>` | Slow standard non-parametric bootstrap (recommended >= 100); performs the bootstrap plus the analysis on the original alignment and gives a consensus tree. `-bc` omits the original-alignment analysis, `-bo` does only the bootstrap. |
| `-bb <n>` | UFBoot approximation, orders of magnitude faster with relatively unbiased supports (>= 1000). Legacy spelling of the modern `-B <n>`. |

## Interpreting branch lengths

PoMo branch lengths are **not** directly comparable to DNA-model branch lengths. They measure mutations **and** frequency shifts (drift) per site; the drift share grows with `N`. To convert to substitutions per site, divide by `N^2`. The `.iqtree` report prints both scalings:

```
NOTE: The branch lengths of PoMo measure mutations and frequency shifts.
To compare PoMo branch lengths to DNA substitution models use the tree length
measured in substitutions per site.

Total tree length (sum of branch lengths)
 - measured in number of mutations and frequency shifts per site: 0.71200751
 - measured in number of substitutions per site (divided by N^2): 0.00879022
Sum of internal branch lengths
- measured in mutations and frequency shifts per site: 0.01767814 (2.48285810% of tree length)
- measured in substitutions per site: 0.01767814 (2.48285810% of tree length)
```

Rules of thumb:

- Report the **substitutions-per-site** tree length when comparing against a plain DNA analysis.
- Individual branch lengths in the Newick `.treefile` are in the polymorphism-scaled units (mutations + frequency shifts per site); apply the `/ N^2` conversion yourself if you need substitution units per branch.
- A very small internal-branch share of total tree length (as in the example, ~2.5%) means most signal sits on terminal branches.

## Runtime and scaling expectations

| Factor | Effect |
| --- | --- |
| Number of Gamma rate categories | Runtime scales **linearly** with the number of categories, because the full transition matrix cannot simply be scaled — the mutational component must be separated from genetic drift. |
| Larger `N` | More parameters/state space in the virtual population; slower. |
| Threads | `-nt <k>` (legacy) / `-T <k>` (modern) parallelises; `-nt AUTO` picks the count. |
| Fully specified model (`HKY{6.0}+P{0.0025}`) | IQ-TREE skips parameter estimation and only performs a tree search — much cheaper. |

## Identical and polymorphic sequences

- Polymorphic sites are the point of the method: within a population column the counts are split across bases (e.g. `0,1,4,0`), and the site-frequency information across those columns is what PoMo exploits. Monomorphic columns (`0,0,5,0`) are kept and contribute fixation information.
- IQ-TREE's general handling of duplicate sequences applies: by default identical sequences are removed for the analysis and added back at the end (among a group of identical sequences the first two are kept, later ones are ignored but re-added, so that bootstrap support for such a group is not forced to 100%). Force retention with `-keep-ident`. The de-duplicated alignment is written with the `.uniqueseq.phy` suffix.
- Practical consequence for PoMo inputs: collapse redundant chromosomes into **counts** rather than passing dozens of identical FASTA records — that is exactly why the counts format exists.

## Pitfalls

- Forgetting `+P` (`-m HKY` alone) silently runs a plain DNA model on the counts-derived alignment — no polymorphism modelling.
- Using fewer than ~5 chromosomes per population: too little frequency information; lower `N` accordingly or abandon PoMo.
- Trying invalid `N` values: only odd 3-19 plus 2 and 10 are offered.
- Expecting `+R` FreeRate or `+I` invariance with PoMo: only discrete Gamma (`+G<k>`) is documented as supported.
- Comparing PoMo branch lengths with DNA branch lengths without dividing by `N^2`.
- Assuming modern option spellings: the PoMo documentation uses `iqtree`, `-nt`, `-pre`, `-b`, `-bb`. Confirm with `iqtree -h`.

## Related skills

- `iqtree-model-selection` — `-m TEST`/`MFP`, DNA model catalogue, frequency/rate modifiers.
- `iqtree-complex-models` — partitioned and mixture models for multi-locus data.
- `iqtree-branch-support` — UFBoot details, convergence diagnostics, SH-aLRT.
- `iqtree-tree-search` — starting trees, constrained search, topology tests.
- `iqtree-concordance-factor` — quantifying discordance once PoMo gene trees exist.
- `iqtree-rooting` — rooting the resulting species tree (outgroup or Rootstrap).
- `iqtree-alisim` — simulating under models for benchmarking.
- `iqtree-command-reference` — full option tables.
- `iqtree` — install, formats, output overview, FAQ.

## Resources

- `references/polymorphism-aware-models.md` — the whole PoMo guide: counts-file format and example, `example.cf`, `cflib`, the `+P`/`+N`/`+WB`/`+WH`/`+S`/`+FO`/`+G4` model syntax with runnable commands, bootstrap options, and the branch-length conversion block. Load this first for any PoMo question.

Targeted greps to pull only the needed section:

```bash
cd /path/to/iqtree-poms
grep -n -A 35 "^## Counts files" references/polymorphism-aware-models.md
grep -n -A 30 "^## Virtual population size" references/polymorphism-aware-models.md
grep -n -A 26 "^## Sampling method" references/polymorphism-aware-models.md
grep -n -A 16 "^## Bootstrap branch support" references/polymorphism-aware-models.md
grep -n -A 16 "^## Interpretation of branch lengths" references/polymorphism-aware-models.md
```

## Unverified in the supplied documentation

The following were requested topics that the local PoMo reference does **not** document; confirm each against `iqtree -h` (heading `POLYMORPHISM AWARE MODELS (PoMo)`) in the installed binary before using, and treat the descriptions here as unconfirmed:

| Item asked about | Status |
| --- | --- |
| `--pop-file`, `--pop-strategy` | Not present in these references. |
| `POPCONV` conversion-tool flags | Not present; the docs point to the external Python library `cflib` instead. |
| VCF import | Not described. |
| Dedicated POM model names `POM`, `POM2`, `GTR2.0`, `GTR6.0`, `MERGE2`, `SPLIT6` | Not present. The documented route is `+P` appended to ordinary DNA model names (`HKY+P`, `GTR+P`). Note that `GTR2` in these docs is the **binary-data** model, unrelated to PoMo. |
| `--mu-rate`, `--virtual-samples` | Not present; the documented control is the `+N<k>` virtual population size in the model string. |
| Genetic-distance-based sequence sampling (`--dist-matrix`, `--sampling-method` standard vs random) | Not present. The documented "Sampling method" section refers to the likelihood-weighting options `+WB`/`+WH`/`+S` only. |
| Partitioned PoMo model selection `-m POMERGE` / `-m POPMODEL` | Not present. The documented partition-merging keywords are `TESTMERGE`, `TESTNEWMERGE`, `MF+MERGE`, `MFP+MERGE`, none specific to PoMo. |
