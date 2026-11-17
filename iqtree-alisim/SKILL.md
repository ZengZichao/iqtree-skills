---
name: iqtree-alisim
description: Simulate DNA, amino-acid, codon, binary, morphological, or mixed sequence alignments with AliSim (iqtree3 --alisim / iqtree2 --alisim). Use for generating test MSAs from a tree plus a model, parametric bootstrap, method benchmarking, and machine-learning training data. Triggers include --alisim, simulating alignments, .phy/.fa output, GTR+I+G4/R rates, indels, partitioned or multi-tree simulations, random Yule/birth-death trees, mixture and GHOST models, pre-defined mutations, and parallel -nt or MPI runs.
---

# AliSim Sequence Simulation

## Overview

AliSim is the phylogenetic sequence simulator bundled into IQ-TREE 2 (**requires IQ-TREE 2.4.0 or later** for the full feature set described here) and also shipped standalone. It simulates biologically realistic alignments under a wide range of complex models by adaptively combining rate-matrix and probability-matrix computation: ~1.3 h and 1.3 GB RAM for one million sequences or sites, versus 2-5 h and 50-500 GB for Seq-Gen, Dawg and INDELible.

Invoke it through IQ-TREE:

```
iqtree3 --alisim <OUTPUT_PREFIX> -m <MODEL> -t <TREEFILE>
```

Default behaviour: writes `OUTPUT_PREFIX.phy` (PHYLIP), **1000 sites**, one alignment, DNA/JC if no `-m` given.

## Command skeleton

| Task | Minimal command |
| --- | --- |
| Tree + model | `iqtree3 --alisim aln -m "JC" -t tree.nwk` |
| Set length | `iqtree3 --alisim aln -m "JC" -t tree.nwk --length 5000` |
| FASTA out | `iqtree3 --alisim aln -m "JC" -t tree.nwk --out-format fasta` |
| Reproducible | add `-seed 123` to any of the above |
| Mimic a real alignment | `iqtree3 --alisim aln_mimic -s example.phy` |
| Random tree, no alignment needed | `iqtree3 --alisim aln_yh -t "RANDOM{yh/1000}"` |
| Many alignments | add `--num-alignments 3` |
| Partitioned | `iqtree3 --alisim part -q multi_genes.nex -t tree.nwk` |

## Full option reference (verbatim from the docs)

| Option | Usage and meaning |
| --- | --- |
| `--alisim <OUTPUT_ALIGNMENT>` | Activate AliSim and specify the output alignment filename. |
| `-t <TREE_FILE>` | Set the input tree file name. |
| `--seqtype <SEQUENCE_TYPE>` | Sequence type (`BIN`, `DNA`, `AA`, `CODON`, `MORPH{NUMBER_STATES}`). By default AliSim auto-detects from the model name. |
| `-m <MODEL>` | Model name (see Substitution Models / Complex Models). |
| `-st <DATA_TYPE>` | Number of states for morphology, e.g. `-st "MORPH{20}"`. |
| `-mdef <MODEL_FILE>` | NEXUS model file defining new models usable with `-m`. |
| `--fundi <TAXON_1>,...,<TAXON_N>,<RHO>` | FunDi model; last number RHO = proportion of sites randomly permuted within the sequences of the given taxa. |
| `--indel <INS>,<DEL>` | Insertion/deletion rates relative to the substitution rate. |
| `--indel-size <INS_DIS>,<DEL_DIS>` | Indel size distributions. Default `POW{1.7/100}` (Zipfian exponent 1.7, max size 100). |
| `--no-unaligned` | Do not output the unaligned-sequences file when using indels. Default writes `.unaligned.fa`. |
| `-q` / `-p` / `-Q <PARTITION>` | Edge-equal / edge-proportional / edge-unlinked (topology-unlinked) partition models. |
| `--distribution <FILE>` | File defining custom parameter lists/distributions. |
| `--branch-distribution <DISTRIBUTION>` | Draw (and override) branch lengths from a list or built-in distribution. |
| `--branch-scale <SCALE>` | Scale all branch lengths of the input tree. |
| `--length <LENGTH>` | Root sequence length. Default: 1,000. |
| `--num-alignments <NUMBER>` | Number of output datasets. Default: 1. |
| `--root-seq <ALN_FILE>,<SEQ_NAME>` | Use a sequence from an alignment as the root sequence; sets `--length` to its length. |
| `-t RANDOM{<MODEL>/<NUM_TAXA>}` | Generate a random tree (see below). |
| `-rlen <MIN_LEN> <MEAN_LEN> <MAX_LEN>` | Branch-length bounds for random trees. Default `-rlen 0.001 0.1 0.999`. |
| `-s <ALIGNMENT>` | Input alignment in PHYLIP, FASTA, NEXUS, CLUSTAL or MSF format (mimic mode). |
| `--no-copy-gaps` | Disable copying gaps from the input alignment. |
| `--site-freq <OPTION>` | `MEAN` (default), `SAMPLING`, or `MODEL` — how to mimic site frequencies under mixtures. |
| `--site-rate <OPTION>` | `MEAN` (default), `SAMPLING`, or `MODEL` — how to mimic discrete rate heterogeneity. |
| `-nt <NUM_THREADS>` | Threads for OpenMP simulation of large alignments. |
| `--openmp-alg <ALG>` | `IM` (AliSim-OpenMP-IM, default) or `EM` (AliSim-OpenMP-EM). |
| `--mem-limit <FACTOR>` | Memory limit factor for OpenMP-IM, `0 < FACTOR <= 1`. Default 0.2. |
| `--keep-seq-order` | With OpenMP-EM, keep tip preorder instead of arbitrary order. |
| `--no-merge` | With OpenMP-EM, skip concatenation; emit K sub-alignments of L/K sites. |
| `--single-output` | Output all alignments into a single file. |
| `--write-all` | Enable outputting internal sequences. |
| `-seed <NUMBER>` | Random seed. Default: PC clock — always set it for reproducibility. |
| `-gz` | Compress output (slower). Default off. |
| `--out-format <FORMAT>` | `fasta`, `phy`, or `maple`. Default `phy`. (`-af fasta` is the short alias used in the placement recipe; it is not in the option table, so prefer `--out-format`.) |
| `--mutation <FILE>` | Pre-defined mutations per node/taxon. |
| `--index-from-one` | Mutation site indexing starts at 1 instead of 0 (VGsim compatibility). |
| `--indel-size` distributions | `GEO{mean}`, `NB{mean/var}`, `POW{a/max}`, `LAV{a/max}` (values shown in the indel section). |
| `--no-export-sequence-wo-gaps` | Do not write the `<prefix>_withoutgaps.fa` gap-free copy produced by indel runs. |

### Options requested but NOT found in these references (unverified — do not use without checking `--alisim` help)

The following spellings appear nowhere in `references/alisim.md` or `references/recipe-simulating-multi-tree.md`: `--tree-height`, `--write-tree`, `--write-clocktree`, `--no-length`, `--write-ancestral`, `--haa`, `--hasr`, `--stratiestick`, `--model-only`, `--taxon-num`, `--tree-length`, `--length-template` / `-lt`, `--heterotachy`, `--install-del-rates`, `--ref-gap`, `--mutations` (the documented spelling is `--mutation`). Treat each as **unverified**; the documented equivalents are:

| Intent | Verified way to do it |
| --- | --- |
| Save the simulated tree | AliSim writes a tree file next to the alignment — recipes consume it as `sequences.treefile` (placement recipe, line 139) and `-t concat.treefile` is the usual inference output of the same name. **Partially verified**: no `--write-tree` style switch exists in these references, so confirm with `ls <PREFIX>*` after each run |
| Ultrametric / clock-like tree | Generate with a Birth-Death random tree `-t "RANDOM{bd{0.1/0.05}/1000}"` |
| Fix total tree depth | Scale every branch with `--branch-scale <SCALE>` |
| Per-partition tree lengths | Append `{tree_length}` to each partition in the `charpartition` statement, used with `-p` |
| Non-overlapping taxon sets per locus | One NEWICK per line in the tree file + `-Q` |
| Heterotachy | `[&model=HKY{2.0}*H4,lengths=...]` in the tree, or `+H4{weights}` in `-m` |
| Mutations file | `--mutation mutations.txt` |
| Ancestral/internal sequences | `--write-all` |

## Datatypes and models

| Data type | Command |
| --- | --- |
| DNA | `iqtree3 --alisim aln_dna -m "JC" -t tree.nwk` |
| Amino acid (empirical, e.g. LG) | `iqtree3 --alisim aln_aa -m "LG" -t tree.nwk` |
| Codon (dN/dS = 2.0, unequal NT freqs) | `iqtree3 --alisim aln_codon -m "MG{2.0}+F1X4{0.2/0.3/0.4/0.1}" -t tree.nwk` |
| Binary | `iqtree3 --alisim aln_bin -m "JC2" -t tree.nwk` |
| Morphological, 20 states | `iqtree3 --alisim aln_morph -m "MK" -t tree.nwk -st "MORPH{20}"` |
| Morphology w/o constant sites (ascertainment) | `iqtree3 --alisim aln_morph_asc -m "MK+ASC" -t tree.nwk -st "MORPH{20}"` |
| Non-reversible Lie Markov 12.12 (= UNREST) | `iqtree3 --alisim aln_lie -m "12.12{0.5/0.6/0.9/0.2/0.1/0.4/0.7/0.8/0.3/0.15/0.65}+F{0.1/0.2/0.4/0.3}" -t tree.nwk` |
| Non-reversible protein | `-m NONREV` |

Multifurcating (polytomy) input trees are fully supported, e.g. `(A:0.3544,(B:0.1905,C:0.1328,D:0.0898):0.05,E:0.1);`.

## Rate heterogeneity across sites

| Syntax | Meaning | Example |
| --- | --- | --- |
| `+I{p}` | proportion of invariant sites | `-m "JC+I{0.2}"` |
| `+GC{shape}` | continuous Gamma | `-m "JC+GC{0.5}"` |
| `+Gk{shape}` | k discrete Gamma categories | `-m "JC+G4{0.5}"` |
| `+Rk{w1/r1/.../wk/rk}` | Free-rate, weights then rates | `-m "JC+R3{0.5,1.5,0.2,0.7,0.3,2.0}"` → rates 1.5/0.7/2.0 with weights 0.5/0.2/0.3 |
| combinable | `+I` with `+G`/`+R` | `-m "JC+G4{0.5}+I{0.2}"` |

Omitting the braces lets AliSim draw parameters randomly (see "Specifying exact parameters").

## Insertions and deletions

```
iqtree3 --alisim alignment_indel -m "JC" -t tree.nwk --indel 0.03,0.1
```

`--indel INS,DEL` gives insertion and deletion rates relative to substitutions (here 3 insertions and 10 deletions per 100 substitutions). Extra outputs: `alignment_indel_withoutgaps.fa` (gap-free sequences; disable with `--no-export-sequence-wo-gaps`) and a `.unaligned.fa` file (disable with `--no-unaligned`).

Indel sizes default to Zipfian `POW{1.7/100}` (exponent 1.7, max 100). Override with `--indel-size`:

```
iqtree3 --alisim x -m "JC" -t tree.nwk --indel 0.1,0.05 --indel-size "GEO{5},GEO{4}"          # geometric, mean 5 / mean 4 (variance derived from mean)
iqtree3 --alisim x -m "JC" -t tree.nwk --indel 0.1,0.05 --indel-size "NB{5/20},POW{1.5/10}"   # NegBinomial mean5/var20 ; Zipfian a=1.5 max=10
iqtree3 --alisim x -m "JC" -t tree.nwk --indel 0.1,0.05 --indel-size "LAV{1.5/10},LAV{1.5/10}" # Lavalette both sides
```

With indels, `--length` is the **root sequence** length — the alignment can come out longer because gaps are inserted.

## Specifying exact parameters vs. random draws

Any model may carry fixed parameters in braces; separators are `/` inside a component and `,` between rate-category pairs:

```
iqtree3 --alisim aln_HKY -t tree.nwk -m "HKY{2.0}+F{0.2/0.3/0.1/0.4}"     # kappa=2, freqs A C G T
iqtree3 --alisim aln_HKY -t tree.nwk -m "HKY{2.0}"                        # freqs drawn randomly
```

If state frequencies are neither specified nor inferable from an input alignment, AliSim draws them from empirical distributions (Generalized-logistic, Exponential-normal, Power-log-normal, Exponential-Weibull) fitted to large empirical datasets. Sequencing error can be added with `+E{p}`:

```
iqtree3 --alisim aln_err -t tree.nwk -m "HKY{2.0}+F{0.2/0.3/0.1/0.4}+E{0.01}"
```

### Custom parameter distributions

Write one list per line (`name` then space-separated numbers) into `custom_distributions.txt`, then reference list names inside `{...}` and load the file with `--distribution`:

```
F_A 0.363799 0.313203 0.277533 0.24235 0.260252
R_A 10.363799 20.313203 10.277533 5.24235 3.26025
I_A 0.257679 0.417313 0.290922 0.301826 0.292324
```

```
iqtree3 --alisim aln_GTR_custom -t tree.nwk \
  -m "GTR{1.5/R_A/R_B/0.5/R_C}+F{Generalized_logistic/0.3/F_A/0.2}+I{F_D}+G{F_C}" \
  --distribution custom_distributions.txt
```

Frequencies are renormalized to sum to 1. Lists may replace substitution rates, state frequencies, dN/dS, transition/transversion rates and category weights for DNA, protein, codon, binary, morphological, Lie Markov, heterotachy and mixture models. Built-in names usable in place of a list: `uniform`, `Generalized_logistic`, `Exponential_normal`, `Power_log_normal`, `Exponential_Weibull`.

## Mimicking a real alignment

```
iqtree3 --alisim alignment_mimic -s example.phy
```

Internally runs IQ-TREE to infer the tree and best-fit ModelFinder model, simulates an alignment of the same length, and copies the gap pattern of `example.phy` into the output (`--no-copy-gaps` to stop). Combine with an explicit mixture/discrete-rate model and control how posterior categories are turned into per-site values:

```
iqtree3 --alisim aln_mimic -s example.phy \
  -m "MIX{GTR{2/3/4/5/6}+F{0.2/0.3/0.4/0.1},HKY{2}+F{0.3/0.2/0.1/0.4},JC}+G{0.5}" \
  --site-freq SAMPLING --site-rate MODEL
```

`--site-freq`/`--site-rate`: `MEAN` (default) uses the posterior mean, `SAMPLING` samples from the posterior, `MODEL` uses the component/model specification. Caveats: with `--length` and no indel model the output may be shorter or longer than the original and gaps are only copied over the overlapping region; with an indel model the original gaps are ignored entirely.

## Reusing parameters from an inference run

Every IQ-TREE inference run prints a ready-made AliSim command at the end of its `.iqtree` report (shown here for the `concat` prefix; this is how the concordance-factor workflow reuses model parameters):

```
tail concat.iqtree
```

```
ALISIM COMMAND
--------------
--alisim simulated_MSA -t concat.treefile -m "Q.plant+I{0.177536}+R8{0.147295,0.0935335,...}" --length 432014
```

Copy everything after `-m` (drop the `--length` part unless you want that length) to simulate under the exact ML model you just estimated. Always substitute your own values.

## Simulating along a random tree

```
iqtree3 --alisim alignment_yh -t "RANDOM{yh/1000}"                              # Yule-Harding, 1000 taxa
iqtree3 --alisim alignment_bd -t "RANDOM{bd{0.1/0.05}/1000}"                    # birth-death b=0.1 d=0.05
```

Tree generators: `yh` (Yule-Harding), `bd{birth/death}` (Birth-Death), `u` (Uniform), `cat` (Caterpillar), `bal` (Balanced). Taxon count `<NUM_TAXA>` may be a number, a list `{10/20/50}`, or `U{LOWER/UPPER}`. Branch lengths follow an exponential distribution whose minimum/mean/maximum you set with `-rlen MIN MEAN MAX` (default `0.001 0.1 0.999`). No `-m` means JC on DNA.

```
iqtree3 --alisim aln_rand -t "RANDOM{yh/1000}" --branch-distribution F_A --distribution custom_distributions.txt
```

`--branch-distribution` overrides branch lengths of random **and user-provided** trees with draws from a list or built-in distribution; `--branch-scale <SCALE>` multiplies all input branch lengths instead.

## Branch-specific models

Tag individual branches in the Newick with `[&model=...]`:

```
(A:0.1,(B:0.1,C:0.2[&model=HKY]),(D:0.3,E:0.1[&model=GTR{0.5/1.7/3.4/2.3/1.9}+F{0.2/0.3/0.4/0.1}+I{0.2}+G{0.5}]):0.2);
iqtree3 --alisim alignment_example_1 -t input_tree.nwk -m "JC"
```

The `-m` model applies to untagged branches. Two rules: a branch model is **not inherited** by descendants (tag every branch of a subtree to change the whole subtree), and rate-heterogeneity specifications on different branches operate independently.

Root sequence from chosen frequencies (rooted tree, `[&model=...,freqs=...]`):

```
(A:0.1,(B:0.1,C:0.2),(D:0.3,E:0.1):0.2):0.3[&model=GTR,freqs=0.2/0.3/0.1/0.4];
```

Per-branch heterotachy via multiple lengths:

```
(A:0.1,(B:0.1,C:0.2[&model=HKY{2.0}*H4,lengths=0.1/0.2/0.15/0.3]),(D:0.3,E:0.1):0.2);
```

## Partitioned simulations

Choose the option letter according to what partitions share — this is the single most common AliSim mistake:

| Option | Model | Requirement |
| --- | --- | --- |
| `-q` | edge-**equal**: all partitions share one tree and one set of branch lengths | plain partition file |
| `-p` | edge-**proportional**: shared topology, partition-specific tree lengths | `charset` entries carry `{tree_length}` |
| `-Q` | topology-**unlinked**: each partition has its own tree (taxon sets may differ) | multi-NEWICK tree file, one tree per line |

Partition file (`multi_genes.nex`):

```
#nexus
begin sets;
    charset gene_1 = DNA, 1-846;
    charset gene_2 = DNA, 847-1368;
    charpartition mine = HKY{2}+F{0.2/0.3/0.1/0.4}:gene_1,
        GTR{1.2/0.8/0.5/1.7/1.5}+F{0.1/0.2/0.3/0.4}+G{0.5}:gene_2,
        JC:gene_3, HKY{1.5}+I{0.2}:gene_4, K80{2.5}:gene_5;
end;
```

```
iqtree3 --alisim partition_multi_genes -q multi_genes.nex -t tree.nwk     # -> partition_multi_genes.phy
```

Edge-proportional — append the tree length per partition (`gene_1{0.26019}` rescales branches so that partition's total tree length is 0.26019; equals input-tree length × `partition_rate`, e.g. 0.26019/0.8673 = 0.3):

```
iqtree3 --alisim partition -p multi_genes.nex -t tree.nwk
```

Using `-q` here silently ignores the per-partition rates.

Topology-unlinked — `multi_trees.nwk` holds five trees (one per gene), some missing taxa:

```
iqtree3 --alisim multi_alignment -Q multi_genes.nex -t multi_trees.nwk
```

Missing taxa are padded with gaps (e.g. taxon D in `gene_2`). Using `-q` here would instead load only the **first** tree.

Mixed datatypes in one run — declare states in the `charset` lines (`DNA`, `AA`, `MORPH{6}`, `MORPH{30}`) and AliSim writes one file per datatype:

```
iqtree3 --alisim partition_mix -q example_mix.nex -t tree.nwk
# -> partition_mix_DNA.phy, partition_mix_MORPH6.phy, partition_mix_MORPH30.phy, partition_mix_AA.phy
```

## Mixture models and heterotachy (GHOST)

```
iqtree3 --alisim aln_mix_C10 -m "C10" -t tree.nwk                                  # profile mixture
iqtree3 --alisim aln_mix -m "MIX{JC,HKY{2.0}+F{0.2/0.4/0.1/0.3}}+G4" -t tree.nwk   # user-defined mixture
iqtree3 --alisim aln_new -m "C10" -mdef my_models.nex -t tree.nwk                  # model name defined in the -mdef file
```

GHOST with 4 categories needs a **multi-length tree** — each branch carries `n` slash-separated lengths:

```
(A[0.067/0.151/0.562/1.269],(B[0.001/0.078/0.319/1.724],C[0.076/0.101/0.002/1.061])[0.043/0.086/0.003/0.002],D[0.002/0.136/0.002/0.001]);
```

```
iqtree3 --alisim aln_ghost -m "GTR{2/3/4/5/6}+F{0.2/0.3/0.1/0.4}+H4{0.15/0.2/0.35/0.3}" -t ghost_tree.nwk
iqtree3 --alisim aln_ghost_unlink -m "MIX{GTR{2/3/4/5/6},GTR{3/4/5/6/7},GTR{4/5/6/7/8},GTR{5/6/7/8/9}}+F{0.2/0.3/0.1/0.4}*H4{0.15/0.2/0.35/0.3}" -t ghost_tree.nwk
```

`*H4` unlinks GTR parameters per category; add per-component `+F{...}` inside `MIX{...}` to unlink frequencies too. Omit the weights for a uniform distribution. To mimic an alignment known to evolve under GHOST: `-m "GTR+H4"`, `"GTR*H4"` (unlink rates) or `"GTR+FO*H4"` (unlink rates and frequencies) together with `-s example.phy`.

## Functional divergence (FunDi)

```
iqtree3 --alisim alignment_fundi -t tree.nwk -m "JC" --fundi A,C,0.1
```

With `RHO = 0.1` and the default 1000-site length, 100 random sites in the sequences of taxa A and C are permuted with each other (same permutation applied to both).

## Pre-defining mutations

Mutations file (`mutations.txt`) against tree `(T1:0.2,(T2:0.3,T4:0.1)Node5:0.4,T3:0.1);`:

```
Node5	C39G,T17A,G25C
T2		C25A,A5G
```

Format: taxon or internal-node name, whitespace, comma-separated `parentState + position + newState`. **Positions are 0-based** (so `C39G` forces the 40th site) for VGsim compatibility; use `--index-from-one` for 1-based input.

```
iqtree3 --alisim example_mutations -t tree_example.nwk -m "JC" --mutation mutations.txt
```

## Parallel simulation

```
iqtree3 --alisim large_alignment -t tree.nwk --length 1000000 -m "JC" -nt 4
iqtree3 --alisim large -t tree.nwk --length 1000000 -m "JC" -nt 4 --openmp-alg EM --keep-seq-order
mpirun -np 10 iqtree3-mpi --alisim many_alignment -t tree.nwk -m "JC" --num-alignments 100
mpirun -np 10 --map-by node:PE=4 --rank-by core iqtree3-mpi --alisim many_large -t tree.nwk --length 1000000 -m "JC" --num-alignments 100 -nt 4
```

- OpenMP algorithms: `IM` (default) or `EM`. `EM` writes sequences in arbitrary order (`--keep-seq-order` restores preorder at a runtime cost) and can skip merging with `--no-merge`, producing K sub-alignments of L/K sites.
- `IM` performance depends on `--mem-limit <FACTOR>` (default 0.2): small → slower, large → more RAM.
- Choose the number of MPI processes as a divisor of `--num-alignments` for best efficiency.
- **Indels require MPI** — the OpenMP implementation does not support them.

## End-to-end worked examples

### Example 1 — basic DNA simulation and sanity check

Input `tree.nwk`:

```
(A:0.3544,(B:0.1905,C:0.1328):0.0998,D:0.0898);
```

```
iqtree3 --alisim alignment -m "JC" -t tree.nwk -seed 123 --out-format fasta
```

Check `alignment.phy` / `alignment.fa`: header `4 1000` (4 sequences × 1000 sites by default), taxon names match the tree, no gaps under JC. Rerun with the same `-seed` to confirm byte-identical output; rerun **without** `-seed` twice and diff — identical files mean the clock-based seed collided, which is why the docs insist on setting it.

### Example 2 — replicate set for a parametric bootstrap

```
iqtree3 --alisim boot -m "GTR{1.2/0.8/0.5/1.7/1.5}+F{0.1/0.2/0.3/0.4}+G4{0.5}" -t tree.nwk \
  --num-alignments 100 --length 5000 --out-format phy -seed 2024
```

Check: files `boot_1.phy` … `boot_100.phy` exist, each with the right taxon count and 5000 columns; add `--single-output` to get one concatenated file instead.

### Example 3 — concatenate loci simulated under different trees

```
iqtree3 --alisim part1 -t RANDOM{yh/20} -m GTR+I+R8 --length 9000
iqtree3 --alisim part2 -t RANDOM{yh/20} -m GTR+I+R8 --length 1000
AMAS.py concat -i part1.phy part2.phy -f phylip -u phylip -t part1_part2_concat.phy -d dna
```

Check: final alignment has 20 sequences and 10000 columns; verify each intermediate file exists before concatenating (the recipe warns AMAS fails opaquely if a simulation errored). Prefer the native route for the same goal: one partition file with `{tree_length}` per charset plus `-p`, or `-Q` with a multi-tree file.

### Example 4 — feed a simulated alignment into a full inference

```
iqtree3 --alisim sim -t "RANDOM{bd{0.1/0.05}/100}" -m "HKY{2.0}+G4" --length 20000 -seed 7
iqtree3 -s sim.phy -m MFP -B 1000 -T AUTO -pre inferred
```

Check: `inferred.treefile` recovers the clades of the tree AliSim used (list the run's outputs with `ls sim*` and compare topologies, e.g. with Robinson-Foulds), and the ModelFinder-selected model is close to HKY+G — a mismatch usually means too few sites or an unsuitable branch-length scale.

## Common pitfalls

| Symptom | Cause / fix |
| --- | --- |
| Two runs gave identical alignments | Clock-derived seed collision — pass distinct `-seed` values |
| Alignment longer than `--length` | Expected with indels: `--length` is the root sequence length |
| Per-partition rates ignored | Used `-q` instead of `-p` |
| Only the first tree was used | Used `-q`/`-p` instead of `-Q` for topology-unlinked data |
| Missing taxa vanished from output | Normal — `-Q` pads absent taxa with gaps |
| No `.unaligned.fa` wanted | Add `--no-unaligned`; gap-free copy controlled by `--no-export-sequence-wo-gaps` |
| Indel run ignores threads | Indels are MPI-only; use `iqtree3-mpi` |
| Wrong positions mutated | Mutation coordinates are 0-based; add `--index-from-one` |
| Base frequencies unexpected | Unspecified `+F` is randomly drawn — always pin `+F{...}` for reproducible models |
| `--alisim` unknown option | IQ-TREE older than 2.4.0 |

## Loading only what you need

```bash
# Complete option table (last section of the manual)
grep -n -A 40 "^## All AliSim options" references/alisim.md

# Specific features
grep -n -A 30 "^## Insertion and deletion models" references/alisim.md
grep -n -A 66 "^## Specifying model parameters" references/alisim.md
grep -n -A 32 "^## Mimicking a real alignment" references/alisim.md
grep -n -A 32 "^## Simulating along a random tree" references/alisim.md
grep -n -A 52 "^## Branch-specific models" references/alisim.md
grep -n -A 115 "^## Partition models" references/alisim.md
grep -n -A 48 "^## Heterotachy GHOST model" references/alisim.md
grep -n -A 32 "^## Parallel sequence simulations" references/alisim.md

# Section index for the rest
grep -n "^## " references/alisim.md
```

## Resources

- `references/alisim.md` — the complete AliSim manual (~6.1k words): every feature section plus the authoritative "All AliSim options" table at line 692. Load the relevant section with the grep recipes above; always cross-check an option there before using it.
- `references/recipe-simulating-multi-tree.md` — recipe for alignments assembled from multiple trees: two-locus example, the AMAS concatenation workflow, and a bash loop for arbitrary numbers of loci/lengths. Note its own header points to the simpler native "Topology-unlinked partition model" (`-Q`).

## Related skills

`iqtree` (entry point and version checks), `iqtree-model-selection` and `iqtree-complex-models` (model strings reused verbatim in `-m`), `iqtree-tree-search` and `iqtree-branch-support` (analysing simulated data), `iqtree-concordance-factor` (uses the `ALISIM COMMAND` block for large-dataset speed-ups), `iqtree-dating` (simulate dated/pathogen benchmarks), `iqtree-command-reference`, `iqtree-recipes`.
