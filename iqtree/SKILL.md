---
name: iqtree
description: Entry point for IQ-TREE 2/3 maximum-likelihood phylogenetics - install or run iqtree/iqtree3 on a PHYLIP/FASTA/NEXUS/CLUSTAL/MSF alignment (.phy, .fasta, .nex), pick the datatype with -st (DNA, AA, CODON, BIN, MORPH, NT2AA), infer a first ML tree, set core counts via -T/-T AUTO, read .iqtree/.treefile/.log/.contree output, resume from a .ckp.gz checkpoint with -redo, use the web server, fix common errors (+ASC invariant sites, duplicate taxon names), and compile from source with CMake.
---

# IQ-TREE 2/3: Getting Started and Core Workflow

## Overview

IQ-TREE is a command-line maximum-likelihood (ML) phylogenetics program that jointly performs model selection (ModelFinder), tree search, branch-support assessment and complex-model fitting from a single multiple-sequence alignment. This skill covers installation, the first inference run, input formats and datatypes, multicore usage, output interpretation, checkpointing, big-data strategy, the web server, FAQ troubleshooting, and compiling from source.

## Version and option-name hazards

IQ-TREE 1.x, 2.x and 3.x differ in both binary name and option spelling. Confirm which is installed before building a command line (`iqtree -h`, `iqtree2 -h`, or `iqtree3 -h`).

| Task | IQ-TREE 1.x | IQ-TREE 2.x / 3.x |
| --- | --- | --- |
| Binary | `iqtree` (or `iqtree-omp`) | `iqtree`, compiled 3.x binary is `iqtree3` |
| CPU cores | `-nt` | `-T` |
| Ultrafast bootstrap | `-bb` | `-B` |
| Output prefix | `-pre` | `--prefix` |
| Edge-proportional partitions | `-spp` | `-p` |
| Edge-unlinked partitions | `-sp` | `-Q` |
| Partition resampling | `-bsam` | `--sampling` |
| Per-site rates | `-wsr` | `--rate` |

Prefer the modern spelling and note the legacy equivalent in a comment when the user's version is unknown.

## Installation

Pick one:

| Route | Command |
| --- | --- |
| Debian/Ubuntu | `sudo apt-get install iqtree` |
| Arch (AUR) | `yay -S iqtree` or `paru -S iqtree` |
| Arch (BioArchLinux) | `sudo pacman -S iqtree` |
| Conda/bioconda | `conda install -c bioconda iqtree` |
| Homebrew | `brew install brewsci/bio/iqtree3` |
| FreeBSD | `pkg install iqtree` |
| Manual download | <http://www.iqtree.org/#download>; extract `iqtree-X.Y.Z-OS`, copy everything in `bin/` onto `PATH` |

Run from a terminal/Command Prompt; double-clicking `iqtree.exe` does nothing. Quick sanity check on the bundled example alignment:

```bash
iqtree -s example.phy          # 7 DNA seqs x 28 sites, ModelFinder default
```

## First inference run

```bash
# Minimal: ModelFinder selects the best-fit model, then full ML tree search
iqtree -s example.phy

# Known model, 4 cores, plus UFBoot and SH-aLRT supports (the standard publishable run)
iqtree -s example.phy -m GTR+I+G -T 4 -B 1000 -alrt 1000

# Skip re-doing model selection once a best-fit model is known
iqtree -s example.phy -m TIM2+I+G -B 1000 -bnni --prefix run1
```

Key facts to carry forward:

- `-s` (alignment) is the only mandatory input for a plain run.
- Since 1.5.4 the default model behaviour **is** `-m MFP`, so `iqtree -s X` already runs ModelFinder Plus.
- The ML tree is always **unrooted**; IQ-TREE merely draws the first taxon in the alignment at the root. Rooting/outgroups are a separate concern.
- Sequence names may only contain alphanumerics, `_`, `-`, `.`, `\`, `|`; all other characters are replaced by `_`. Collisions after sanitising cause a fatal duplicate-name error — rename in the alignment, not in the tree.
- A directory passed instead of a file is loaded and concatenated automatically (v2 onwards).

## Input alignment formats

Accepted: PHYLIP, FASTA, NEXUS, CLUSTAL, MSF. PHYLIP needs the `nseq nsite` header line; FASTA does not. Unaligned raw sequences must be aligned first (MAFFT/ClustalW) — IQ-TREE never aligns.

## Datatype selection (`-st`)

IQ-TREE auto-detects the sequence type; override only when detection is wrong or when a codon/translation model is wanted.

| Goal | Command | Notes |
| --- | --- | --- |
| Nucleotide | `iqtree -s dna.phy -st DNA` | Default for ACGT alignments |
| Amino acid | `iqtree -s prot.phy -st AA` | Empirical AA matrices |
| Codon model | `iqtree -s cds.phy -st CODON` | **Always required** — otherwise DNA models are applied. Alignment length must be divisible by 3 |
| Alternative genetic code | `iqtree -s cds.phy -st CODON5` | `CODON{1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25}`; IDs follow the NCBI table |
| Binary (0/1) | `iqtree -s bin.phy -st BIN` | Alphabet size auto-detected |
| Morphological | `iqtree -s morph.phy -st MORPH` | States `0`-`9` and `A`-`Z` (for states 10-31); pair with `-m MK+ASC` when constant sites are absent |
| Translate DNA to AA | `iqtree -s cds.phy -st NT2AA` | Analysis then runs on the translation; `-st NT2AA5` for a non-standard code |

Under codon models branch lengths are substitutions per **codon** site (~3x longer than under DNA models). SNP data with no constant sites need ascertainment correction: `-m GTR+ASC`, or `-m MFP+ASC` to keep it inside model selection.

## Multicore usage

```bash
iqtree -s example.phy -m GTR+I+G -T 4        # fixed core count
iqtree -s example.phy -m GTR+I+G -T AUTO     # benchmark and pick the best count
iqtree -s example.phy -m GTR+I+G -T AUTO -ntmax 8   # cap what AUTO may allocate
```

Parallelism runs along the alignment, so efficiency collapses on short alignments — using more cores can be *slower*. Trust the `BEST NUMBER OF THREADS:` line emitted by `-T AUTO` and reuse that number for later runs on the same data. MPI/hybrid runs use `mpirun -np N iqtree-mpi -nt M ...`; the OpenMP `iqtree` binary must **not** be launched under `mpirun`.

## Output files

All outputs are prefixed by the alignment filename unless `--prefix`/`-pre` changes it.

| Suffix | Contents | When it appears |
| --- | --- | --- |
| `.iqtree` | Main self-readable report: best model, parameters, tree with supports, topology tests | Always — read this first |
| `.treefile` | ML tree in NEWICK (supports embedded if bootstrapping) | Always |
| `.log` | Full transcript of the run | Always — attach this to any bug report |
| `.model` / `.model.gz` | Log-likelihoods of every tested model; doubles as model-selection checkpoint | Model testing runs |
| `.ckp.gz` | Gzip checkpoint of the run | Written periodically |
| `.contree` | Consensus tree with assigned supports, branch lengths optimised on the real alignment | With `-B`/`-b` |
| `.splits.nex` | Support % for **all** splits, including conflicting ones excluded from the consensus (view in SplitsTree) | With `-B`/`-b` |
| `.splits` | Same as above in star-dot format | With `-wsplits` |
| `.varsites` | Variable-sites-only copy of the alignment | When `+ASC` was requested but invariant sites exist |
| `.uniqueseq.phy` | Deduplicated alignment (max two copies per identical group retained) | Identical sequences present |
| `.state` | Ancestral state reconstruction per node/site | Ancestral sequence reconstruction runs |
| `.rate` / `.mlrate` | Per-site rates (empirical Bayes / ML) | Per-site rate estimation runs |
| `.sitelh` | Per-site log-likelihoods (TREE-PUZZLE format; feeds CONSEL) | Site log-likelihood output requested |
| `.suptree` | Input tree with support values mapped on | With `-sup` |
| `.conaln` | Concatenated alignment built from a multi-file NEXUS partition scheme | Mixed-data partition runs |
| `.GTRPMIX.nex` | Optimised linked exchangeabilities | With `--link-exchange` |

## Checkpointing and resuming

Rerunning the **exact same command line** resumes from the last checkpoint and appends to the existing `.log`. A finished run refuses to rerun:

```
ERROR: Checkpoint (example.phy.ckp.gz) indicates that a previous run successfully finished
Use `-redo` option if you really want to redo the analysis and overwrite all output files.
```

```bash
iqtree -s example.phy -redo            # force a full redo, overwriting outputs
iqtree -s example.phy -m MFP -mredo    # ignore only the .model checkpoint, keep tree run
iqtree -s example.phy -cptime 60       # checkpoint interval (default 20 s)
```

Never edit or delete `.ckp.gz` (checkpoint lost) or `.model.gz` (model scores lost). Use distinct `--prefix` values for parallel analyses of the same alignment so outputs do not collide.

## Analyzing big data (>=1000 sequences or >=10000 sites)

The dedicated page is a stub in the official docs; assemble the strategy from these documented mechanisms:

- `-safe` — safe numerical mode avoiding likelihood underflow for thousands of sequences; auto-enabled above ~2000 sequences.
- `-mem 64G` — cap RAM usage.
- Checkpointing + `-T AUTO`/`-ntmax` — survive cluster wall-time limits, resume freely.
- Prefer single branch tests (`-alrt 1000`, `-abayes`, `-lbp 1000`) over `-b`/`-B` on very large trees; they are explicitly recommended above ~10,000 taxa.
- If splitting standard bootstrap across machines, run `N` jobs with `-bo` and distinct prefixes, `cat boot*.boottrees > alltrees`, then `iqtree -con -t alltrees` and `iqtree -sup input.treefile -t alltrees`.
- For expensive profile-mixture models, use the PMSF approximation (see `iqtree-complex-models`).
- Fix a guide tree with `-te` to avoid a full search intertwined with parameter estimation.

## Web server

For small/medium jobs without local installation: <http://iqtree.cibiv.univie.ac.at>. Upload a Phylip/Fasta/Nexus/Clustal/MSF alignment; Tree Inference runs defaults plus auto model detection, UFBoot and SH-aLRT; the Model Selection tab finds a model without tree reconstruction. Results zip includes `.iqtree`, `.log`, `.treefile`, `.svg`, `.pdf`, `.contree`, and `.ckp.gz` when the job was stopped. Limits: 24 h CPU and 1 GB RAM — download a stopped job and resume locally from its checkpoint.

## Compiling from source

Requirements: C++ compiler (GCC >= 4.8, Clang, MSVC or Intel), CMake >= 2.8.10, Boost, Eigen3 (required for IQ-TREE >= 1.6), OpenMP. Clang binaries run ~5-10% faster than GCC.

```bash
git clone --recursive https://github.com/iqtree/iqtree3.git
cd iqtree3 && mkdir build && cd build
cmake ..                     # -DEIGEN3_INCLUDE_DIR=<dir> .. if Eigen3 is not found
make -j                      # produces ./iqtree3
```

Variants via `-DIQTREE_FLAGS=`: `omp` (multicore), `mpi` → `iqtree-mpi`, `omp-mpi` (hybrid), `m32`/`omp-m32` (32-bit, max 2 GB RAM, no AVX — avoid for big data), `KNL` (AVX-512 Xeon Phi), `static-libcxx`. macOS: `cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..`; ARM builds need GCC 10 exactly or Clang >= 14 (macOS ARM: Clang >= 17). Library build: `cmake -DBUILD_LIB=ON ..`. Windows: use Clang with `-G "MinGW Makefiles"` and `mingw32-make -j4`.

## FAQ triage

| Symptom | Cause and fix |
| --- | --- |
| `ERROR: Invaid use of +ASC because of ... invariant sites` / `Skipped since +ASC is not applicable` | Alignment has constant, partially-constant or ambiguously-constant columns. Use the auto-written `.varsites` alignment, or add constant sites back with `-fconst` |
| Duplicate sequence names error | Illegal characters were rewritten to `_` and collided — rename in the alignment |
| Run won't restart | Finished-run checkpoint guard — add `-redo` |
| More cores = slower | Alignment too short for along-site parallelism — use `-T AUTO` or fewer cores |
| Sequences flagged `failed` by the composition chi-square test | Character composition deviates from the alignment mean. Explorative only; do **not** delete by default. Investigate per partition, or try protein profile mixtures (`C10`..`C60`) which absorb compositional heterogeneity |
| Unexpectedly high supports on a concatenated tree | UFBoot/SH-aLRT thresholds (UFBoot >= 95%, SH-aLRT >= 80%) apply to **single-gene** trees only; for phylogenomic concatenation compute concordance factors instead |
| Standard bootstrap infeasible | Use `-B` (UFBoot) with `-bnni`, or split `-bo` across machines and combine `.boottrees` |
| Crash / assertion failure | Report with description, first and last 10 lines of `.log`, the assertion text; send data privately to developers |
| `*` or `!` in protein alignment | Stop codons translate to `*`; frameshift/padding positions to `!`; both count as unknown AA |

## Routing table

| User goal | Load this skill |
| --- | ---|
| "Which model should I use?", `-m TEST`/`MFP`, restrict candidate models, rate heterogeneity `+I/+G/+R`, `+F/+FO/+FQ`, `+ASC`, estimate a new AA matrix (QMaker/nQMaker) | `iqtree-model-selection` |
| Multi-gene/supermatrix partitions, partition files, merging schemes, mixture models, `C10`-`C60`, PMSF, GHOST `+H`, MAST multitree | `iqtree-complex-models` |
| Bootstrap, SH-aLRT, aBayes, local bootstrap probability, UFBoot stopping rules, interpreting supports | `iqtree-branch-support` |
| Tree search control, starting/fixed trees, constrained search `-g`, NNI/SPR, pathogen/CMAPLE search | `iqtree-tree-search` |
| Exact syntax, defaults, or an option not listed here | `iqtree-command-reference` |
| Tip dates, molecular clock dating, LSD2, MCMCtree gradients | `iqtree-dating` |
| Simulating alignments, AliSim, benchmarks | `iqtree-alisim` |
| Gene/site concordance factors, species-tree vs gene-tree discordance, phylogenetic placement | `iqtree-concordance-factor` |
| Rooting without an outgroup, non-reversible models, root confidence | `iqtree-rooting` |
| PoMo / polymorphism-aware species trees from population counts | `iqtree-poms` |
| End-to-end worked examples combining IQ-TREE with other tools | `iqtree-recipes` |
| Topology tests (AU/KH/SH), likelihood mapping, symmetry tests, alignment trimming | `iqtree-recipes` + `iqtree-command-reference`; see also `references/advanced-tutorial.md` in `iqtree-complex-models` |

## Resources

| File | What's in it | Load when |
| --- | --- | --- |
| `references/quickstart.md` | Install commands per platform, Windows/macOS terminal walkthrough, minimal command-line examples, v1-to-v2 option renames | Setting IQ-TREE up, or needing a canonical short command |
| `references/home.md` | What IQ-TREE is, key features/citations, web server, support channels, full per-feature citation list | Answering "what is IQ-TREE", writing Methods/citations sections |
| `references/tutorial.md` | Beginner tutorial: input formats and name rules, first run, output files, checkpoint/`-redo`/`--prefix`, model selection intro, codon and binary/morph/SNP data, UFBoot vs `-b` vs SH-aLRT, `-T AUTO` | Any first-run, datatype or output-file question |
| `references/documentation-map.md` | Every official doc page with its section anchors and target file names | Deciding which upstream document holds a topic |
| `references/faq.md` | Help/bug protocol, UFBoot interpretation, gap/ambiguity tables (IUPAC + AA), mixed DNA+protein partitions, codon+DNA branch-length meaning, composition test, thread choice, split-bootstrap recipe, `+ASC` errors, identical sequences, sites vs patterns, folder concatenation | Troubleshooting an error message or an odd result |
| `references/analyzing-big-data.md` | Big-data page (stub: `>=1000` seqs or `>=10000` sites) plus the doc index footer | Confirming the big-data guidance is unlisted upstream and falling back to the mechanisms above |
| `references/web-server-tutorial.md` | Three web-server tabs, supported uploads, result-file suffix table, 24 h / 1 GB limits, resuming a stopped server job | User prefers the online server or hits a server limit |
| `references/compilation-guide.md` | Requirements, Linux/macOS/Windows builds, `IQTREE_FLAGS` variants (omp/mpi/m32/KNL/static), lib build, Xcode project, precompiled-binary notes | Building from source or a custom HPC target |
| `references/developer-guide.md` | Internal `Alignment`, substitution-model and `PhyloTree` class design | Extending IQ-TREE or explaining why a model is unsupported |

Targeted extraction instead of loading a whole file:

```bash
grep -n -A 30 "^## Utilizing multi-core CPUs" references/tutorial.md
grep -n -A 20 "^## How does IQ-TREE treat gap/missing/ambiguous characters" references/faq.md
grep -n -A 20 "^## Why does IQ-TREE complain" references/faq.md
grep -n -A 25 "^## Compiling MPI version" references/compilation-guide.md
grep -n -A 68 "^## Minimal command-line examples" references/quickstart.md
grep -n "for version 1.x" references/quickstart.md
```
