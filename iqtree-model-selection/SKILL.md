---
name: iqtree-model-selection
description: Choose, test and specify IQ-TREE substitution models - run ModelFinder (-m TEST, TESTONLY, MF, MFP) with -mset, -msub, -madd, -mrate, -cmax, -merit on a .phy/.fasta/.nex alignment; pick a model by hand from the DNA (JC/HKY/TN/TIM2/TVM/GTR or a 6-digit rate code), protein (LG/WAG/JTT/mtREV/Q.pfam/NQ.yeast/GTR20), codon (MG/GY/GY2K/ECMK07/ECMrest) or binary/morphological (JC2/GTR2/MK/ORDERED) catalogues; add rate heterogeneity (+I, +G4, +R4, +I+G), frequencies (+F/+FO/+FQ) and ascertainment bias +ASC for SNP data; estimate a new amino-acid matrix with QMaker (--model-joint, --link-exchange).
---

# IQ-TREE Model Selection (ModelFinder) and Substitution Models

## Overview

IQ-TREE selects an evolutionary model either automatically (ModelFinder scores many candidate models against AIC/AICc/BIC on a starting tree) or explicitly (`-m MODEL+Freq+Rate`). This skill covers the automatic-selection options, the complete model catalogues per datatype, rate-heterogeneity and state-frequency modifiers, ascertainment-bias correction, and estimating new amino-acid exchangeability matrices from one or many alignments.

## Model-string grammar

```
-m MODEL[+FreqType][+RateType]        e.g.  GTR+F+G4,  LG+C20+F+G,  MK+FQ+ASC
```

`MODEL` is a base model name, a 6-digit DNA rate code, `TEST`/`MF`/`MFP` for automatic selection, `MIX{...}`/`FMIX{...}` for mixtures, or a **filename** holding user-defined rates. Quote any string containing `{}` so the shell does not expand it.

## Automatic selection

| Option | Behaviour |
| --- | --- |
| `-m TESTONLY` | Standard selection (jModelTest/ProtTest equivalent): no FreeRate models. With a partition file, finds a model per partition. No tree reconstruction |
| `-m TEST` | `-m TESTONLY` followed by tree inference under the winner |
| `-m MF` / `-m TESTNEWONLY` | Extended selection additionally testing FreeRate (`+R`) models. No tree reconstruction. Recommended replacement for `TESTONLY` |
| `-m MFP` / `-m TESTNEW` | `-m MF` followed by tree inference. **This is the IQ-TREE default since 1.5.4**, i.e. `iqtree -s X` == `iqtree -s X -m MFP` |
| `-m MF+MERGE`, `-m MFP+MERGE`, `-m TESTMERGE`, `-m TESTMERGEONLY` | Partition-scheme merging variants — see `iqtree-complex-models` |
| Append `+ASC` (e.g. `-m MFP+ASC`) | Restrict model testing to `+ASC` models — required for SNP/morphological data lacking constant sites |
| Append `+LM`, `+LMRY`, `+LMWS`, `+LMMK` | Additionally test Lie Markov DNA models (all / RY symmetry / WS symmetry / MK symmetry) |

Selection ranks models by BIC by default; switch with `-merit AIC` or `-merit AICc`. Scores are computed on a fixed starting tree unless `-mtree` triggers a full search per candidate (accurate, expensive — only with spare compute).

### Narrowing the search space (speed vs thoroughness)

| Option | Meaning / default |
| --- | --- |
| `-mset <list>` | Comma-separated base models to consider, e.g. `-mset WAG,LG,JTT`; also accepts a program name `raxml`, `phyml`, `mrbayes` to restrict to that tool's catalogue |
| `-msub <src>` | Restrict AA models by origin: `nuclear`, `mitochondrial`, `chloroplast`, `viral` |
| `-madd <list>` | Additionally include mixture models, e.g. `-madd LG4M,LG4X` (LG4X is excluded from `-m MF` by default because it is also a mixture model) |
| `-mfreq <list>` | Frequency types to test. Defaults: `FU,F` for protein; `,F1x4,F3x4,F` for codon |
| `-mrate <list>` | Rate types to test. Defaults: `E,I,G,I+G` for `-m TEST*`, `E,I,G,I+G,R` for `-m MF*` (`E` = equal/homogeneous) |
| `-cmax <k>` | Maximum FreeRate categories tested. Default `10`; `-cmax 15` tests up to `+R15` on long alignments. `-cmin` sets the lower bound (default 2) |
| `-merit AIC\|AICc\|BIC` | Optimality criterion used to rank candidates. Default: all three are computed and BIC decides |
| `-mtree` | Full tree search per tested model (default: fixed starting tree) |
| `-mredo` | Ignore the existing `.model.gz` checkpoint and re-test everything |
| `-mdef <file.nex>` | Load user-defined models from a NEXUS model file before selection |

Comma-separated lists must contain **no spaces**. Interrupted selection resumes automatically from the `.model` / `.model.gz` checkpoint. For the authoritative defaults of every `-m*` knob, grep the sibling command reference:

```bash
grep -n -A 15 "^## Automatic model selection" ../iqtree-command-reference/references/command-reference.md
```

### Canonical commands

```bash
iqtree -s example.phy -m MFP                 # select model then infer tree (the default)
iqtree -s example.phy -m MF                  # model selection only, no tree search
iqtree -s example.phy -m TESTONLY            # jModelTest/ProtTest-style, no FreeRate
iqtree -s example.phy -m MF -mset WAG,LG     # restrict AA matrices to save time
iqtree -s example.phy -m MF -msub nuclear    # general AA models only
iqtree -s example.phy -m MF -cmax 15         # test +R2 .. +R15
iqtree -s example.phy -m MF -mtree           # thorough: tree search per model
iqtree -s prot.phy  -m MF -madd LG4M,LG4X    # include protein mixture models
iqtree -s SNP_data.phy -m MFP+ASC            # force ascertainment-corrected candidates
iqtree -s example.phy -m MFP -merit AIC      # rank by AIC instead of BIC
```

Read the ranked table in `<prefix>.model` (all tested models) and the winner plus its parameter estimates in `<prefix>.iqtree`. Once known, hard-code the model to skip re-testing: `iqtree -s example.phy -m TIM2+I+G`.

## DNA models

Ordered by complexity; `df` = free parameters, `Code` = equality constraints on the six relative rates A-C, A-G, A-T, C-G, C-T, G-T.

| Model | df | Code | Model | df | Code |
| --- | --- | --- | --- | --- | --- |
| JC / JC69 | 0 | 000000 | TPM3 / TPM3u | 2 / 5 | 012012 |
| F81 | 3 | 000000 | TIM / TIMe | 6 / 3 | 012230 |
| K80 / K2P | 1 | 010010 | TIM2 / TIM2e | 6 / 3 | 010232 |
| HKY / HKY85 | 4 | 010010 | TIM3 / TIM3e | 6 / 3 | 012032 |
| TN / TN93 | 5 | 010020 | TVM / TVMe | 7 / 4 | 012314 |
| TNe | 2 | 010020 | SYM | 5 | 012345 |
| K81 / K3P | 2 | 012210 | GTR | 8 | 012345 |
| K81u | 5 | 012210 | TPM2 / TPM2u | 2 / 5 | 010212 |

Digit ordering is irrelevant (`101101` == `010010`); the G-T rate and all rates sharing its digit are fixed to 1.0 because rates are relative. Any restricted model may be given directly as a code: `-m 012012+G`.

Fix rates inline with curly braces, ordered by first appearance of each distinct digit:
`GTR{1.0,2.0,1.5,3.7,2.8}` (A-C, A-G, A-T, C-G, C-T; G-T = 1.0), `TIM2{4.39,5.30,12.1}` == `GTR{4.39,5.30,4.39,1.0,12.1}`, `TN{2.0,3.0}+G8{0.5}+I{0.15}`.

Lie Markov models (v1.6+, many non-reversible, consistent under concatenation unlike GTR) are named `<npar>.<variant>` with optional prefix `RY` (default), `WS`, `MK`: from `1.1` (== JC) through `12.12` (== UNREST). See the reference table for reversibility and free-frequency counts. Non-reversible variants are the entry point for root inference — see `iqtree-rooting`.

## Protein models

Empirical exchangeability matrices, grouped by intended data source (`-msub` accepts the nuclear / mitochondrial / chloroplast / viral groupings marked below):

| Region | Models |
| --- | --- |
| nuclear (general) | LG, WAG, JTT, JTTDCMut, Dayhoff, DCMut, VT, PMB, Blosum62*, Q.pfam, Q.mammal, Q.bird, Q.insect, Q.plant, Q.yeast, ELM, EAL, G.pfam |
| mitochondrial | mtREV, mtMAM, mtART, mtZOA, mtMet, mtVer, mtInv |
| chloroplast | cpREV |
| viral | HIVb, HIVw, FLU, rtREV |
| none / general | Poisson (equal rates and frequencies) |
| non-reversible | NQ.pfam, NQ.mammal, NQ.bird, NQ.insect, NQ.plant, NQ.yeast |
| parameter-rich | GTR20 (190 rates — unstable without very long alignments) |

\* Blosum62 was designed for alignment scoring, not phylogenetics; avoid it.

Profile and multi-matrix protein mixtures (`C10`–`C60`, `EX2`, `EX3`, `EHO`, `UL2`, `UL3`, `EX_EHO`, `LG4M`, `LG4X`, `CF4`, `S10pfamC60`, `S28pfamC59`, `S28pfamC60`, `S10pfamC10`, `S20pfamC20`, `S30pfamC30`) combine with an exchangeability matrix: `LG+C20`, `LG+C20+F`, `JTT+CF4+G`, `S28pfamC60+G`, and Gamma may be replaced by FreeRate (`LG+C20+R4`). Note `C10`–`C60` assume Poisson replacement and implicitly include Gamma rate heterogeneity. Fitting these is covered by `iqtree-complex-models`.

User-defined matrix: if the name matches no built-in model, IQ-TREE reads it as a **file in PAML format** (lower-triangular 190 exchangeabilities + 20 frequencies, AA order `A R N D C Q E G H I L K M F P S T W Y V`).

## Codon models

Require `-st CODON` (or `-st CODON<codeID>`); otherwise DNA models are applied.

| Model | Description |
| --- | --- |
| MG | Muse-Gaut, omega = dn/ds |
| MGK | MG + ts/tv ratio |
| MG1KTS / MGKAP2, MG1KTV / MGKAP3, MG2K / MGKAP4 | MG with ts only, tv only, or both rates (Kosiol 2007) |
| GY0K / GYKAP1 | Goldman-Yang, omega only |
| GY | GY0K + ts/tv |
| GY1KTS / GYKAP2, GY1KTV / GYKAP3, GY2K / GYKAP4 | GY with ts, tv, or both |
| ECM / ECMK07 / KOSI07, ECMrest, ECMS05 / SCHN05 | Empirical codon models — standard genetic code only |
| `<empirical>_<mechanistic>` | Combined product model, e.g. `ECMK07_GY2K` |

Parameters inline: `MG2K{1.0,0.3,0.5}` = omega 1.0, ts 0.3, tv 0.5. Unrecognised names are read as files holding a non-stop-codon exchangeability matrix, frequencies and codon list (e.g. Kosiol's `ECMrest.dat`). Genetic-code IDs follow NCBI: `CODON1` standard, `CODON2` vertebrate mito, `CODON3` yeast mito, `CODON5` invertebrate mito, `CODON11` bacterial/plastid, etc. Branch lengths count substitutions per **codon** site (~3x DNA).

## Binary and morphological models

Valid states: binary `0`/`1`; morphological `0`-`9` and `A`-`V` (states 10-31). Alphabet size auto-detected.

| Model | Description |
| --- | --- |
| JC2 | Jukes-Cantor type for binary data (equal frequencies) |
| GTR2 | General time reversible for binary data — the only one with unequal state frequencies |
| MK | Equal-rates JC-type model for morphology |
| GTRX (alias `GTR`) | Unequal rates for multistate data — heavily parameterised, verify fit |
| ORDERED | Only neighbouring states may exchange |

Common practice for morphology is `MK+FQ+ASC` (or `ORDERED+FQ+ASC` for additive ordered characters). Models with unequal rates/frequencies (`GTR2+FO`, `MK+FO`, `GTRX+FQ`, `GTRX+FO`) are only defensible when state labels are non-arbitrary (recoded amino acids, certain genomic characters), not for ordinary transformational morphological characters.

## Rate heterogeneity across sites

| Type | Meaning |
| --- | --- |
| `+I` | Proportion of invariable sites |
| `+G` | Discrete Gamma, default 4 categories; change with `+G8` |
| `+I+G` | Invariable sites + discrete Gamma |
| `+R` | FreeRate (relaxes the Gamma assumption); categories via `+R6`, default 4. Fits better than `+G` and is recommended for large datasets |
| `+I+R` | Invariable sites + FreeRate |
| `+GC` | Continuous Gamma — AliSim simulation only |

Fix parameters inline: `+I{0.2}`, `+G{0.9}`, `+I{0.2}+G{0.9}`, FreeRate `+Rk{w1,r1,...,wk,rk}`; `-a` sets the Gamma shape and `-i` the pinvar proportion (both estimated by default). `+G` uses the **mean** approximation (as RAxML/PhyML do); `-gmedian` switches to the median approximation used by TREE-PUZZLE, which makes log-likelihoods incomparable across programs. `-m MFP` tests FreeRate models, plain `-m TEST` does not.

## State frequencies

| Type | Meaning |
| --- | --- |
| `+F` | Empirical frequencies counted from the alignment. Default whenever the model has unequal frequencies |
| `+FO` | Frequencies optimised by ML (letter O, not zero) |
| `+FQ` | Equal frequencies |
| `+F1X4` | Codon: unequal nucleotide frequencies, identical across the three positions |
| `+F3X4` | Codon: unequal nucleotide frequencies per codon position |

Default codon frequency is `+F3X4` for `MG`-type, `+F` for `GY`-type, model-supplied for empirical codon models. Fix them explicitly: `GTR+F{0.1,0.2,0.3,0.4}` (must sum to 1), or supply 20 values after `+F{...}` in the standard AA order. For protein data the model's own frequencies apply unless `+F`/`+FO`/`+FQ` is appended (e.g. `WAG+G+F`).

## Ascertainment bias correction (`+ASC`)

Apply `+ASC` when the alignment contains **no constant sites** — SNP panels and most morphological matrices: `GTR+ASC`, `MK+ASC`, `MK+FQ+ASC`. It conditions the likelihood on variable sites; omitting it inflates branch lengths.

If IQ-TREE reports `Invaid use of +ASC because of ... invariant sites` or skips candidates with `Skipped since +ASC is not applicable`, the alignment still holds invariant columns — constant, partially constant (`G`/`-`/`N` mixed), or ambiguously constant (`C`/`Y`/`-`). Use the auto-written `.varsites` alignment (written from v1.5.0), or remove those columns yourself before rerunning.

## Estimating new amino-acid models (QMaker / nQMaker)

Three-step workflow. Step 1 fits per-partition reversible models and a guide tree; step 2 estimates one joint matrix; step 3 extracts it in PAML format for reuse via `-m Q.mygroup`.

```bash
# A. Concatenated alignment + NEXUS partition file
iqtree3 --seed 1 -T AUTO -s alignment.nex -p train.nex -m MFP \
        -mset LG,WAG,JTT -cmax 4 --prefix train_plant          # -> .best_model.nex + .treefile
iqtree3 -seed 1 -T AUTO -s alignment.nex \
        -p train_plant.best_model.nex -te train_plant.treefile \
        --init-model LG --model-joint GTR20+FO --prefix plant_GTR20_FO
grep -A 21 "can be used as input for IQ-TREE" plant_GTR20_FO.iqtree | tail -n20 > Q.plant

# B. Folder of individual alignments (-S gives each alignment its own tree)
iqtree3 -seed 1 -T AUTO -S train_plant -mset LG,WAG,JTT -cmax 4 -pre train_plant
iqtree3 -seed 1 -T AUTO -S train_plant.best_model.nex -te train_plant.treefile \
        --model-joint GTR20+FO --init-model LG -pre train_plant.GTR20
grep -A 21 "can be used as input for IQ-TREE" train_plant.GTR20.iqtree | tail -n20 > Q.plant
```

Swap `--model-joint GTR20+FO` for `--model-joint NONREV+FO` to estimate a **non-reversible** matrix (nQMaker, rooted trees); extract with `grep -A 22 "..." | tail -n21 > NQ.plant`. `--init-model LG` supplies the starting matrix; omit `-p` entirely when there is no partition file. Assess root confidence afterwards with `-B 1000` under the non-reversible model (see `iqtree-rooting`).

Linked exchangeabilities under profile mixtures (GTRpmix, v2.3.5+):

```bash
iqtree3 -s aln.phy -m GTR20+C60+G4 --link-exchange -te guide.treefile -me 0.99
iqtree3 -s aln.phy -mdef run.GTRPMIX.nex -m GTRPMIX+C60+G4     # reuse without re-estimating
```

`--link-exchange` links the GTR20 rates across all 60 classes (otherwise 60 matrices are estimated). `-te` fixes a guide tree and avoids a search intertwined with estimation; `-me 0.99` loosens the convergence threshold for speed; `--init-exchange LG` starts optimisation from a sensible matrix and doubles as a robustness check across starting points.

Multiple exchangeability + profile models (GTRspmix, v3.1.3+): add the `SXXpfamCYY` series to selection with `-m MFP --madd S10pfamC60+G4,S28pfamC59+G4,S28pfamC60+G4`, or load a custom model file with `-mdef model_best.nex -m SXXFYY+G4`.

## Common pitfalls

| Situation | Action |
| --- | --- |
| Model chosen by BIC differs from AIC choice | Expected; BIC is more parsimonious and is the default criterion |
| Selection takes very long | Add `-mset`, `-msub`, `-cmax 4`, drop `-mtree` |
| `LG4X` never selected by `-m MF` | It is excluded as a mixture model — add `-madd LG4X` |
| Very parameter-rich fit on short alignment (`GTR20`, `GTRX`) | Over-fitting risk; prefer simpler model or gather more sites |
| Codon run silently used a DNA model | `-st CODON` was omitted |
| Log-likelihoods disagree with TREE-PUZZLE | `-gmedian` needed to match the median approximation |
| Shell mangled `MIX{...}` or `+R{...}` | Quote the whole `-m` argument |

## Related skills

`iqtree` (install, first run, output files, datatypes) · `iqtree-complex-models` (partitions, merging, mixtures, PMSF, GHOST, MAST) · `iqtree-command-reference` (exact option syntax and defaults) · `iqtree-tree-search` (what happens after the model is fixed) · `iqtree-branch-support` · `iqtree-rooting` (non-reversible models and root tests) · `iqtree-alisim` (simulating under a chosen model) · `iqtree-poms`.

## Resources

| File | What's in it | Load when |
| --- | --- | --- |
| `references/substitution-models.md` | Full DNA/protein/codon/binary/morphological catalogues with df and 6-digit codes, Lie Markov list, frequency and rate-heterogeneity tables, `+ASC` rules, user-defined PAML-format matrices | Resolving an exact model name, parameter count, allowed modifier, or checking whether a model exists |
| `references/estimating-aa-models.md` | QMaker and nQMaker walkthroughs (single alignment, folder of alignments, non-reversible, linked exchangeabilities, GTRspmix), PAML output format, citations | Building a clade-specific AA matrix or reusing `--model-joint` / `--link-exchange` |

Large files — pull only the needed section:

```bash
grep -n "^## \|^### " references/substitution-models.md          # map all sections first
grep -n -A 30 "^### Base substitution rates" references/substitution-models.md
grep -n -A 45 "^### Lie Markov models" references/substitution-models.md
grep -n -A 45 "^### Amino-acid exchange rate matrices" references/substitution-models.md
grep -n -A 25 "^### Protein mixture models" references/substitution-models.md
grep -n -A 40 "^## Codon models" references/substitution-models.md
grep -n -A 25 "^## Binary and morphological models" references/substitution-models.md
grep -n -A 20 "^## Ascertainment bias correction" references/substitution-models.md
grep -n -A 20 "^## Rate heterogeneity across sites" references/substitution-models.md
grep -n -A 35 "^## Estimating linked exchangeabilities" references/estimating-aa-models.md
```
