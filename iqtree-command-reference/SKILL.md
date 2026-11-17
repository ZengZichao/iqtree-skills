---
name: iqtree-command-reference
description: Assemble, validate and interpret full IQ-TREE command lines flag by flag - the complete option catalogue from Command-Reference (general options -s/-st/-T/-nt/-mem/--prefix/-seed/-o, checkpointing -redo/-cptime/-keep-ident, likelihood mapping -lmap, model selection -m TEST/MF/MFP/TESTMERGE, rate heterogeneity +I/+G/+R, partitions -q/-p/-Q, tree search -g/-n/-nstop/-pers, pathogen -pathogen, UFBoot -B/-bb/-bnni/-bcor/-beps, standard bootstrap -b/-bc/-bo, branch tests -alrt/-abayes/-lbp, ASR -asr, topology tests -z/-zb/-zw/-au, consensus -con/-minsup/-bi/-sup, RF distances -rf, random trees -r/-ru/-rcat/-rbal/-rcsg, misc -alninfo/-fconst/-me/-blmin/-blmax). Also maps IQ-TREE 1.x legacy spellings (-nt, -pre, -bb, -spp, -sp, -wsr) to their 2.x/3.x names, explains every output file (.iqtree .treefile .log .ckp.gz .contree .splits.nex .ufboot .state .rate), and ships three helper scripts (build_command.py, check_env.py, parse_results.py).
---

# IQ-TREE Command Reference

## Overview

This is the flag-level reference for IQ-TREE command lines. Use it when the question is "which option does X", "why is my command rejected", "what spelling does my version use" or "what does this output file contain" - as opposed to the thematic skills (`iqtree-model-selection`, `iqtree-tree-search`, `iqtree-branch-support`, ...) which explain the underlying methodology. Everything here is verified against the official Command Reference page; nothing is invented.

## Legacy (1.x) vs modern (2.x/3.x) spellings

The docs repeat one pattern: *"for version 1.x change X to Y"*. The two styles must never be mixed in one command:

| Concept | Modern (IQ-TREE 2/3) | Legacy (IQ-TREE 1.x) |
| --- | --- | --- |
| CPU cores / AUTO | `-T 4`, `-T AUTO` | `-nt 4`, `-nt AUTO` |
| Max cores for AUTO | `-ntmax 8` (same) | `-ntmax 8` |
| Output prefix | `--prefix run1` | `-pre run1` |
| Ultrafast bootstrap | `-B 1000` | `-bb 1000` |
| Edge-proportional partitions | `-p PART_FILE` | `-spp PART_FILE` |
| Edge-unlinked partitions | `-Q PART_FILE` | `-sp PART_FILE` |
| Write per-site rates | `--rate` | `-wsr` |
| Binary name | `iqtree2` / `iqtree3` | `iqtree` (or `iqtree-omp`, `iqtree.mpich`) |

Unchanged in both: `-s`, `-m`, `-st`, `-seed`, `-o`, `-t`, `-te`, `-g`, `-alrt`, `-q`, `-redo`, `-mem`, `-n`, `-nstop`, `-pers`, `-bnni`.

## Section map of `references/command-reference.md`

| Section | Key options |
| --- | --- |
| General options | `-s` `-st` `-T/-nt` `-ntmax` `-mem` `--prefix/-pre` `-seed` `-o` `-t` `-te` `-redo` `-keep-ident` `-v` `-quiet` |
| Checkpointing | `.ckp.gz` auto-resume, `-redo` to discard, `-cptime`, `.model.gz` (>=1.6) vs `.model` (<=1.5) |
| Likelihood mapping | `-lmap N`, `-lmclust`, `.lmap.svg/.eps`, `.lmap.quad` |
| Automatic model selection | `-m TEST/TESTONLY/MF/MFP`, `TESTMERGE/MF+MERGE`, `-mset` `-msub` `-madd` `-mdef` `-mrate` `-cmin` `-cmax` `-merit` `-mtree` |
| Specifying substitution models | `-m JC/GTR/HKY/LG/WAG/MF...`, mixture `MIX{...}`, `-fconst`, `+ASC`, `+I/+G/+R` |
| Rate heterogeneity | `+I`, `+G[k]`, `+R[k]`, `+I+G`, `--opt-gamma-inv` |
| Partition model options | `-q` `-p/-spp` `-Q/-sp`, partition file formats, `-spartition` |
| Site-specific frequency model | `+FQ`, `-mwopt` |
| Tree search parameters | `-n` `-nstop` `-pers` `-sprrad` `-allnni` `-fast` `-djc` `-g` constraint, `-ninit` `-ntop` `-nbest` |
| Tree search for pathogen data | `--pathogen` `--pathogen-force`, `-T` tuning note |
| Ultrafast bootstrap | `-B/-bb` `-bcor` `-beps` `-bnni` `-bsam` `-nm` `-nstep` `-wbt` `-wbtl` `-uwboot` |
| Nonparametric bootstrap | `-b` `-bc` `-bo` `-wbt` |
| Single branch tests | `-alrt` `-abayes` `-lbp` `-czb` |
| Ancestral reconstruction | `-asr` `-asr-min`, `.state` file |
| Tree topology tests | `-z` `-zb` `-zw` `-au` `-te`, `.iqtree` test table |
| Consensus tree | `-con` `-minsup` `-bi` `-sup` `-suptag` `-scale` |
| Robinson-Foulds distance | `-rf` `-rf_all` `-rf_adj` |
| Random trees | `-r` `-ru` `-rcat` `-rbal` `-rcsg` `-rlen` |
| Miscellaneous | `-alninfo` `-blfix` `-blmin` `-blmax` `-me` `-wpl` `-wt` `-fconst` `-prec` |

`references/lsd2-help.md` is the verbatim `lsd2 -h` output (native LSD2 flags, used through IQ-TREE's `--date-options "..."`); see the `iqtree-dating` skill for the wrapper-to-native mapping. `references/analyzing-big-data.md` is the official "Analyzing big data" page, which is still a TODO stub on the site - practical big-data advice lives in `iqtree-recipes` and the pathogen section of the command reference.

## Helper scripts (`scripts/`)

All three are dependency-free Python 3 and print doc-verified options only; each has a `--self-test` that passes.

### `scripts/build_command.py` - assemble and validate a command line

```
python3 scripts/build_command.py --alignment aln.phy --model GTR+I+G -T 4 \
    --bootstrap 1000 --alrt 1000 --prefix run1          # prints the command (dry run)
python3 scripts/build_command.py --alignment aln.phy --style legacy ...  # 1.x spellings
python3 scripts/build_command.py --alignment aln.phy --extra "-wsplits -quiet"
python3 scripts/build_command.py --alignment aln.phy --run               # execute
```

- Emits one spelling style only (`--style modern|legacy`, default modern).
- Rejects bad values before running: bad `-T`/`-mem`/`-pers` formats, `-t` together with `-te`, `--bootstrap-nni` without `--bootstrap`, missing `--alignment`.
- Warns on bootstrap < 1000, `-ntmax` without AUTO, `TESTMERGE` without a partition file, the cost of `--rate` on big alignments.
- Flags outside its map are passed verbatim through `--extra "..."` (checked against the known-option set; nothing is invented).

### `scripts/check_env.py` - find the binary, recommend threads/memory

```
python3 scripts/check_env.py                 # locate iqtree/iqtree2/iqtree3/... and report
python3 scripts/check_env.py --require       # exit 2 if no executable found
python3 scripts/check_env.py --path /opt/iqtree/bin
```

Reports build flavour (OpenMP/MPI), version, CPU count, RAM, then recommends concrete flags: `-T AUTO` vs `-T N`, `-ntmax N`, `-mem NG` (reserving ~25% headroom).

### `scripts/parse_results.py` - digest a finished run

```
python3 scripts/parse_results.py --prefix run1           # run1.iqtree, run1.log, run1.treefile ...
python3 scripts/parse_results.py --prefix run1 --list-files
python3 scripts/parse_results.py --prefix run1 --verbose
```

Prints an LLM-friendly digest: best-fit model with BIC/AICc, log-likelihood, tree length, taxa/sites/patterns, topology-test table with the 95% confidence sets, CPU/wall time, threads/RAM, checkpoint status, trailing NOTE/WARNING/ERROR lines, and a `.treefile` census. Missing files are reported as ABSENT and still exit 0; exit 1 only when no output file matches the prefix at all.

## Typical agent workflow

1. `scripts/check_env.py --require` - which binary, how many threads to give it.
2. `scripts/build_command.py ...` - assemble a validated command line with correct spellings; append site-specific flags via `--extra`.
3. Run IQ-TREE.
4. `scripts/parse_results.py --prefix ...` - read the results back out.

## Loading only what you need

```bash
# One option group (section names are the "## " headings listed above)
grep -n -A 40 "^## General options" references/command-reference.md
grep -n -A 30 "^## Ultrafast bootstrap parameters" references/command-reference.md
grep -n -A 25 "^## Tree topology tests" references/command-reference.md

# Every occurrence of one flag
grep -n -- "-alrt" references/command-reference.md
grep -n -- "-mem" references/command-reference.md

# Legacy-to-modern rename sentences scattered through the docs
grep -rn "for version 1.x" references/

# Native LSD2 flags (tab-indented in lsd2 -h output)
grep -n "^	-" references/lsd2-help.md
```

## Resources

- `references/command-reference.md` - the full official Command Reference page (all option groups with example usages). Load the relevant `##` section, not the whole file.
- `references/lsd2-help.md` - verbatim `lsd2 -h` output; every native LSD2 flag with defaults.
- `references/analyzing-big-data.md` - official page, TODO stub on the site; do not expect content.
- `scripts/build_command.py`, `scripts/check_env.py`, `scripts/parse_results.py` - command assembly, environment check, result digestion (each `--self-test` verified).

## Related skills

`iqtree` (getting started, first runs, output-file basics), `iqtree-model-selection` (the `-m` catalogue and ModelFinder), `iqtree-tree-search` (search control, topology tests, ASR in depth), `iqtree-branch-support` (UFBoot/bootstrap/branch tests in depth), `iqtree-complex-models` (partition files and mixtures in depth), `iqtree-dating` (`--date`/`--dating` front-end and lsd2 flag mapping), `iqtree-concordance-factor` (`--gcf`/`--scf`), `iqtree-alisim` (`--alisim`), `iqtree-rooting` (`-o`, `--root-test`), `iqtree-recipes` (multi-tool workflows).
