---
name: iqtree-complex-models
description: Fit heterogeneous IQ-TREE models for multi-gene and phylogenomic alignments (.phy/.fasta/.nex supermatrices, mitogenomes) - partitions via -q (edge-equal), -p/-spp (edge-proportional) or -Q/-sp (edge-unlinked), RAxML-style and NEXUS partition files with mixed DNA/protein/codon data; scheme merging with -m TESTMERGE, MF+MERGE, -rcluster; mixtures MIX{...}/FMIX{...}, MixtureFinder (-m MIX+MF, -qmax), profile mixtures C10-C60, LG4M/LG4X, CF4, GTRpmix, PMSF approximation (-ft, -fs) when C60 will not fit in RAM; heterotachy GHOST (+Hk); multitree MAST (+T). Use when partitions must be merged.
---

# IQ-TREE Complex Models: Partitions, Mixtures, Site-Specific Frequencies, Heterotachy, Multitrees

## Overview

Complex models relax the "one model, one tree" assumption: partition models assign different models (and optionally different branch lengths) to named site subsets, while mixture models let each site draw its model from a set without knowing the assignment in advance. This skill covers specifying both, selecting/merging partition schemes, approximating expensive profile mixtures with PMSF, modelling heterotachy with GHOST, and modelling multiple conflicting histories with MAST.

## Which heterogeneity model to reach for

| Signal in the data | Model family | Entry option |
| --- | --- | --- |
| Known loci / genes / codon positions, each needing its own model | Partition | `-p` (edge-proportional) |
| Loci may also differ in topology-driven rate trajectories (heterotachy across partitions) | Edge-unlinked partition | `-Q` |
| Unknown, unassignable shift in composition/rate along the sequence | Profile mixture | `-m LG+C20+F+G` |
| Want an unknown number of Q matrices fitted automatically | MixtureFinder | `-m MIX+MF` |
| Rate variation across sites **and** lineages, no a-priori partition | GHOST | `-m GTR+H4` |
| Conflicting gene histories (ILS, introgression, recombination) in one concatenation | MAST | `-m "GTR+FO+G+T" -te trees.nwk` |
| Expensive `C10`–`C60` run out of memory | PMSF approximation | `-m LG+C20+F+G -ft guide.treefile` |

## Partition models

Three linkages, controlled purely by which option carries the partition file:

| Option | Linkage | Comment |
| --- | --- | --- |
| `-q <file>` | Edge-equal — all partitions share one set of branch lengths | Usually unrealistic; ignores differing speeds between loci |
| `-p <file>` (`-spp` in v1.x) | Edge-proportional — shared topology, each partition gets its own rescaling rate | **Recommended default**; best-performing linkage for concatenated species-tree analyses |
| `-Q <file>` (`-sp` in v1.x) | Edge-unlinked — each partition has its own branch lengths | Most parameter-rich; overfits with many short partitions; accommodates heterotachy |

Fit all three and compare BIC to justify the choice. Bootstrapping a partition model resamples sites **within** partitions; use `--sampling GENE` to resample partitions or `--sampling GENESITE` for the two-level scheme that reduces false 100% supports.

### RAxML-style partition file

Plain text, one `DATATYPE, name = ranges` line per partition. Cannot carry per-partition models or rate types — everything inherits `-m`.

```
DNA, part1 = 1-100
DNA, part2 = 101-384
```

Non-consecutive sites and codon-position extraction (comma-separated ranges):

```
DNA, part1 = 1-100, 250-384
DNA, part2 = 101-249\3, 102-249\3
DNA, part3 = 103-249\3
```

`\3` is a modulo step: `part2` = sites 101,102,104,105,107,…,248,249 and `part3` = 103,106,…,247 — i.e. 1st/2nd/3rd codon positions.

### NEXUS partition file

A `sets` block of `charset` lines plus one `charpartition`. Uniquely allows per-partition models and rate types, per-partition source alignment files, and per-partition datatypes.

```
#nexus
begin sets;
    charset part1 = 1-100;
    charset part2 = 101-384;
    charpartition mine = HKY+G:part1, GTR+I+G:part2;
end;
```

Mixed data types across separate alignment files (site numbering restarts at 1 for every file):

```
#nexus
begin sets;
    charset part1 = dna.phy: 1-100 201-300;
    charset part2 = dna.phy: 101-200;
    charset part3 = prot.phy: 1-400;
    charset part4 = prot.phy: 401-600;
    charset part5 = codon.phy:CODON, *;
    charpartition mine = HKY:part1, GTR+G:part2, LG+G:part3, WAG+I+G:part4, GY:part5;
end;
```

Because filenames are inside the NEXUS file, `-s` may be omitted: `iqtree -p example.nex`. Taxon sets need not match across files — absent parts become missing data, and the joined matrix is written to `<prefix>.conaln`. Mark a coding partition with the `CODON` keyword (else it is read as DNA); append a code ID such as `CODON5` for non-standard genetic codes. `*` means "the whole file". Anywhere a model appears you may use the full model grammar, including 6-digit DNA codes, `+F/+FO`, `+G/+R` and `{}`-fixed parameters.

Branch-length meaning in mixed runs: codon + DNA ⇒ substitutions per **nucleotide** site; codon-only ⇒ per **codon** site (~3x longer); mixing codon, DNA and protein ⇒ per character (nucleotide or amino acid).

### Running partitioned analyses

```bash
iqtree -s example.phy -p example.partitions -m GTR+I+G   # RAxML-style: one rate type for all
iqtree -s example.phy -p example.nex                     # NEXUS: models taken from the file
iqtree -s example.phy -q example.nex                     # edge-equal
iqtree -s example.phy -Q example.nex                     # edge-unlinked
iqtree -s example.phy -p example.nex -B 1000             # UFBoot under a partition model
iqtree -s example.phy -p example.nex -B 1000 --sampling GENESITE
```

Outputs of interest: per-partition models and parameters in `.iqtree`; partition log-likelihoods with `-wpl` (`.partlh`).

## Automatic partition-scheme selection and merging

ModelFinder's greedy merge starts from the full scheme and repeatedly joins two partitions while model fit improves (PartitionFinder-equivalent). Requires a partition file via `-p`/`-spp`.

| Option | Selection only | Then tree search | Includes FreeRate |
| --- | --- | --- | --- |
| `-m TESTMERGEONLY` | yes | no | no (PartitionFinder look-alike) |
| `-m TESTMERGE` | yes | yes | no |
| `-m MF+MERGE` (alias `TESTNEWMERGEONLY`) | yes | no | yes |
| `-m MFP+MERGE` (alias `TESTNEWMERGE`) | yes | yes | yes |

```bash
iqtree -s example.phy -p example.nex -m MFP+MERGE              # best scheme, then ML tree
iqtree -s example.phy -p example.nex -m MF+MERGE               # scheme only
iqtree -s example.phy -p example.nex -m TESTMERGE              # PartitionFinder look-alike
iqtree -s example.phy -p example.nex -m TESTMERGEONLY          # PartitionFinder, no tree
iqtree -s example.phy -p example.nex -m MF+MERGE -rcluster 10  # relaxed clustering, top 10%
```

Speed-ups for large numbers of partitions: `-rcluster <pct>` (relaxed hierarchical clustering, examines only the top `%` of merging schemes), `-rclusterf <pct>` (the faster PartitionFinder2 variant, v1.6+), `-rcluster-max <n>` (absolute cap on partition pairs considered in the merging phase; default the larger of 1000 and 10x the partition count). A merge run also writes the winning scheme as `<prefix>.best_scheme.nex`, which can be fed straight back as a `-p` input on later runs — see `iqtree-rooting` for that two-step pattern. On v1.x, `...MERGE...` always searched an edge-unlinked scheme first regardless of `-spp`; v1.6+ respects the linkage.

## Mixture models

Unlike partitions, mixtures never need a site-to-class map: each site likelihood is the weighted sum over classes, and class weights are estimated.

```bash
iqtree -s example.phy -m "MIX{JC,HKY}"          # two Q-matrix components + free weights
iqtree -s example.phy -m "MIX{JC,HKY}+G4"       # x 4 Gamma rates => 8 effective components
```

Always double-quote model strings containing `{}` — otherwise the shell expands `MIX{JC,HKY}` into `MIXJC MIXHKY`. Add `-mwopt` to optimise mixture weights explicitly (auto-enabled for e.g. `LG+C20+F+G`, but **not** for `LG+C20+G`).

Frequency-only mixtures keep one rate matrix and vary state frequencies along the sequence via `FMIX{...}`, declared in a NEXUS model file and loaded with `-mdef`:

```
#nexus
begin models;
    frequency Fclass1 = 0.02549352 0.01296012 ... ;   # 20 AA freqs, order A R N D C Q E G H I L K M F P S T W Y V
    frequency Fclass2 = ...;
    frequency CF4model = FMIX{empirical,Fclass1,Fclass2,Fclass3,Fclass4};
end;
```

```bash
iqtree -s some_protein.aln -mdef mymodels.nex -m JTT+CF4model+G
```

Whole rate matrices can also be declared (`model LG4M1 = <PAML lower triangle> ...`, then `model LG4M = MIX{LG4M1,...,LG4M4}*G4;`). The `*G4` **fuses** Gamma rates with the classes so there are 4 components instead of 16 — this is how `LG4M`/`LG4X` are built internally. `frequency` and `model` commands may coexist in one file.

### Protein profile mixtures

| Model | Nature |
| --- | --- |
| `C10`, `C20`, `C30`, `C40`, `C50`, `C60` | ML variants of PhyloBayes `CAT`; assume Poisson replacement and implicitly include Gamma rates |
| `EX2`, `EX3`, `EHO`, `EX_EHO`, `UL2`, `UL3` | Structure-aware multi-matrix models (exposed/buried, extended/helix/other) |
| `LG4M`, `LG4X` | Four-matrix LG fusions (`LG4X` uses FreeRate) |
| `CF4` | Five-profile mixture |
| `G.pfam+C60`, `EAL+C60`, `ELM+C60` | Linked exchangeabilities (GTRpmix): `G.pfam` general; `ELM` for nuclear-encoded eukaryote proteins; `EAL` for eukaryote–Archaea relationships |
| `S10pfamC60`, `S28pfamC59`, `S28pfamC60`, `S10pfamC10`, `S20pfamC20`, `S30pfamC30` | GTRspmix: several exchangeability matrices, each linked to a cluster of profiles |

Combine an explicit matrix with a profile set: `LG+C20`, `LG+C20+F` (adds the empirical profile as a 21st class), `JTT+CF4+G`, `S28pfamC60+G`; replace Gamma by FreeRate with `LG+C20+R4`.

### MixtureFinder (automatic number of Q classes)

```bash
iqtree -s example.phy -m MIX+MF            # estimate the optimal Q mixture only
iqtree -s example.phy -m MIX+MFP           # select mixture, then search the tree
iqtree -s example.phy -m MIX+MF -merit AIC # rank by AIC/AICc instead of BIC
iqtree -s example.phy -m MIX+MF -lrt 0.05  # decide class count by likelihood ratio test
iqtree -s example.phy -m MIX+MF -qmax 5 -mset HKY,GTR -mrate E,I,G,I+G
```

`-qmax` caps the number of Q-mixture classes (default 10). `-mset` restricts the base model tested per added class; `-mrate` restricts the rate-heterogeneity options; `-mrate-twice 1` re-estimates the rate model after the best mixture is chosen (default 0). Prefer running this in IQ-TREE 3, whose implementation is substantially improved over IQ-TREE 2. With `-lrt`, `-merit` still governs sub-model choice within each class (default BIC).

## Site-specific frequency models (PMSF)

Posterior mean site frequencies approximate `C10`–`C60` at a fraction of the cost, largely removing long-branch-attraction artefacts. Two phases: estimate class parameters on a **guide tree** and write per-site profiles to `.sitefreq`, then analyse using those profiles instead of the mixture.

```bash
# 1. cheap guide tree
iqtree -s aln.phy -m LG+F+G -T AUTO --prefix guide

# 2. PMSF analysis (needs both the mixture spec and the guide tree)
iqtree -s aln.phy -m LG+C20+F+G -ft guide.treefile

# 3. PMSF + standard nonparametric bootstrap (feasible only because of PMSF)
iqtree -s aln.phy -m LG+C20+F+G -ft guide.treefile -b 100

# split across machines: phase 1 alone on a high-memory node
iqtree -s aln.phy -m LG+C20+F+G -ft guide.treefile -n 0     # stops after phase 1, writes .sitefreq
# then copy the .sitefreq file to a small machine
iqtree -s aln.phy -m LG+C20+F+G -fs aln.phy.sitefreq -b 100 -T 24
```

Without `-ft` the run stays a full mixture model. Phase 1 costs the same RAM as the mixture itself — `-fs` is what removes that burden from phase 2 (it also accepts your own profile file: one line per site, site ID from 1 followed by 20 state frequencies, whitespace separated, AA order `A R N D C Q E G H I L K M F P S T W Y V`). For posterior-**maximum** rather than posterior-mean profiles use `-fmax`; its exact spelling and default are listed in `../iqtree-command-reference/references/command-reference.md` under "Site-specific frequency model options".

## Heterotachy models (GHOST)

GHOST is an edge-unlinked mixture: several site classes, each with its own parameters **and** its own edge lengths on one shared topology, without requiring a pre-defined partition (v1.6+).

```bash
iqtree -s data.fst -m GTR+H4          # 4 classes, GTR rates LINKED across classes
iqtree -s data.fst -m GTR*H4          # UNLINKED: one GTR set per class
iqtree -s data.fst -m GTR+FO*H4       # unlinked + per-class optimised base frequencies
iqtree -s data.fst -m GTR+FO*H4 -wspm # write .siteprob: probability each site belongs to each class
```

`+Hk` links, `*Hk` unlinks. With `+H` only one set of empirical frequencies is inferred and shared; add `+FO` to infer frequencies separately per class. Inspect class memberships through the `-wspm` output before interpreting anything biological.

## Multitree models (MAST)

Mixtures across sites **and** trees: each input topology is a mixture component with its own branch lengths and optionally its own model, frequencies and rate model, returning ML weights per history (v2.3+, IQ-TREE 3 binary recommended).

```bash
iqtree3 -s data.fst -m "GTR+T" -te trees.nwk              # GTR linked across trees
iqtree3 -s data.fst -m "GTR+FO+G+T" -te trees.nwk         # GTR + freqs + RHAS all linked
iqtree3 -s data.fst -m "TMIX{GTR+FO+G,F81+FO+R3,HKY+FO+I}+T" -te trees.nwk   # fully unlinked
iqtree3 -s data.fst -m "GTR+FO+G+TR" -te trees.nwk        # branch-length-restricted variant
iqtree3 -s data.fst -m "GTR+FO+G+T[x,x,y]" -te trees.nwk  # trees 1 and 2 forced to equal weight
```

`trees.nwk` holds one Newick tree per line; the component count in `TMIX{...}` must equal the number of topologies. Choose linkage per element independently:

| Model string | subst | freq | RHAS |
| --- | --- | --- | --- |
| `"TMIX{GTR+FO+G,GTR+FO+G}+T"` | unlinked | unlinked | unlinked |
| `"TMIX{GTR+FO,GTR+FO}+G+T"` | unlinked | unlinked | linked |
| `"TMIX{GTR+F,GTR+F}+G+T"` | unlinked | fixed to alignment counts | linked |
| `"GTR+FO+TMIX{G,G}+T"` | linked | linked | unlinked |
| `"GTR+FO+G+T"` | linked | linked | linked |

`+TR` constrains branches that split the taxon set identically across trees to the same length. Results: `.treefile` contains all fitted trees (topologies matching the input); `.iqtree` lists per-tree parameters and weights in input order. Weights are interpretable as **tree-only** proportions solely in the fully linked case (last row above); whenever more than one substitution model is fitted, a weight reflects the joint model-and-tree class. These models are easy to over-parameterise — check branch lengths and every model parameter before interpreting weights.

## Common pitfalls

| Situation | Action |
| --- | --- |
| Shell broke `MIX{...}` / `FMIX{...}` / `+T[...]` | Quote the entire `-m` value |
| All partitions got the same `+I+G` despite wanting mixed rate types | RAxML-style file cannot express it — convert to NEXUS `charpartition` |
| Coding partition analysed under a DNA model | Tag it `charset pX = cds.phy:CODON, *;` |
| Site numbering continued past the end of the previous file in a multi-file NEXUS scheme | Each file's charsets start at 1 |
| Partition merging far too slow | `-rcluster 10` or `-rclusterf`, and cap with `-rcluster-max` |
| `LG+C60+F+G` exceeds available RAM | Use PMSF (`-ft guide.treefile`), or generate `.sitefreq` with `-n 0` on a big node and continue with `-fs` elsewhere |
| Supports of 100% on a phylogenomic concatenation | Expected artefact — compute concordance factors (`iqtree-concordance-factor`) |
| `-m MFP+MERGE` on v1.x produced an unexpected linkage | Pre-1.6 always searched edge-unlinked schemes; upgrade or interpret the second phase separately |

## Related skills

`iqtree-model-selection` (base models, ModelFinder knobs, `+G/+R/+F/+ASC`, QMaker) · `iqtree` (install, datatypes, output files, checkpointing, threads) · `iqtree-command-reference` (option defaults incl. `-q/-spp/-Q`, `-ft/-fs/-fmax`, `-wspm`, `-rcluster*`) · `iqtree-tree-search` · `iqtree-branch-support` · `iqtree-concordance-factor` · `iqtree-rooting` (linked non-reversible models across partitions via `--model-joint`) · `iqtree-alisim` (simulate under partition/mixture/GHOST models) · `iqtree-dating` (partitions and mixtures for approximate-likelihood dating) · `iqtree-recipes`.

## Resources

| File | What's in it | Load when |
| --- | --- | --- |
| `references/complex-models.md` | Authoritative specs: three partition linkages, RAxML/NEXUS partition syntax (mixed data, CODON keyword), `MIX`/`FMIX` grammar, NEXUS model files, MixtureFinder options, PMSF timings and workflow, GHOST linking rules, MAST linkage table and outputs | Writing or debugging a partition/model file, or choosing among mixture constructions |
| `references/advanced-tutorial.md` | Worked advanced tutorial: partitioned analyses (incl. mixed DNA/protein/codon), partition-scheme merging and `-rcluster`, bootstrapping partitions, constrained search `-g`, topology tests (`-z`, `-zb`, `-au`), consensus/support mapping (`-con`, `-net`, `-sup`), user-defined models, per-site rates (`--rate`, `--mlrate`), alignment trimming (`--robust-phy`) | Needing an end-to-end command sequence for a phylogenomic run or downstream tree comparisons |

Targeted extraction instead of loading a whole file:

```bash
grep -n "^## \|^### " references/complex-models.md                 # map all sections first
grep -n -A 45 "^### Partition file format" references/complex-models.md
grep -n -A 30 "^### Defining mixture models" references/complex-models.md
grep -n -A 40 "^### MixtureFinder" references/complex-models.md
grep -n -A 45 "^## Site-specific frequency models" references/complex-models.md
grep -n -A 30 "^## Heterotachy models" references/complex-models.md
grep -n -A 60 "^## Multitree models" references/complex-models.md
grep -n -A 40 "^## Choosing the right partitioning scheme" references/advanced-tutorial.md
grep -n -A 25 "^## Partitioned analysis with mixed data" references/advanced-tutorial.md
```
