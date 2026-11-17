---
title: Introduction
source_url: https://iqtree.github.io/doc/Home.html
headings: 8
---

# Introduction

> Source: <https://iqtree.github.io/doc/Home.html>

## Table of contents

- [Why IQ-TREE?](#why-iq-tree)
- [Key features](#key-features)
- [Free web server](#free-web-server)
- [User support](#user-support)
- [Documentation](#documentation)
- [How to cite IQ-TREE?](#how-to-cite-iq-tree)
- [Development team](#development-team)
- [Credits and acknowledgements](#credits-and-acknowledgements)

_Last update: Jun 27, 2024, Contributors: Jana Trifinopoulos, Minh Bui, Lam-Tung Nguyen_

## Why IQ-TREE?

Thanks to the recent advent of next-generation sequencing techniques, the amount of phylogenomic/transcriptomic data have been rapidly accumulated. This extremely facilitates resolving many “deep phylogenetic” questions in the tree of life. At the same time it poses major computational challenges to analyze such big data, where most phylogenetic software cannot handle. Moreover, there is a need to develop more complex probabilistic models to adequately capture realistic aspects of genomic sequence evolution.

This trends motivated us to develop the IQ-TREE software with a strong emphasis on phylogenomic inference. Our goals are:

  * **Accuracy** : Proposing novel computational methods that perform better than existing approaches.
  * **Speed** : Allowing fast analysis on big data sets and utilizing high performance computing platforms.
  * **Flexibility** : Facilitating the inclusion of new (phylogenomic) models and sequence data types.
  * **Versatility** : Implementing a broad range of commonly-used maximum likelihood analyses.

IQ-TREE has been developed since 2011 and freely available at <http://www.iqtree.org/> as open-source software under the [GNU-GPL license version 2](http://www.gnu.org/licenses/licenses.en.html). It is actively maintained by the core development team (see below) and a number of collabrators.

The name IQ-TREE comes from the fact that it is the successor of [**IQ** PNNI](http://www.cibiv.at/software/iqpnni/) and [**TREE** -PUZZLE](http://www.tree-puzzle.de/) software.

## Key features

  * **Efficient search algorithm** : Fast and effective stochastic algorithm to reconstruct phylogenetic trees by maximum likelihood. IQ-TREE compares favorably to RAxML and PhyML in terms of likelihood while requiring similar amount of computing time ([Nguyen et al., 2015](https://doi.org/10.1093/molbev/msu300)).
  * **Ultrafast bootstrap** : An ultrafast bootstrap approximation (UFBoot) to assess branch supports. UFBoot is 10 to 40 times faster than RAxML rapid bootstrap and obtains less biased support values ([Minh et al., 2013](https://doi.org/10.1093/molbev/mst024); [Hoang et al., 2018](https://doi.org/10.1093/molbev/msx281)).
  * **Ultrafast model selection** : An ultrafast and automatic model selection (ModelFinder) which is 10 to 100 times faster than jModelTest and ProtTest. ModelFinder also finds best-fit partitioning scheme like PartitionFinder ([Kalyaanamoorthy et al., 2017](https://doi.org/10.1038/nmeth.4285)).
  * **Simulating sequences** : A fast sequence alignment simulator (AliSim) which is much more realistic than Seq-Gen and INDELible ([Ly-Trong et al., 2023](https://doi.org/10.1093/bioinformatics/btad540)).
  * **Big Data Analysis** : Supporting huge datasets with thousands of sequences or millions of alignment sites via [checkpointing](../../iqtree-command-reference/references/command-reference.md#checkpointing-to-resume-stopped-run), safe numerical and low memory mode. [Multicore CPUs](tutorial.md#utilizing-multi-core-cpus) and [parallel MPI system](compilation-guide.md#compiling-mpi-version) are utilized to speedup analysis.
  * **Phylogenetic testing** : Several fast branch tests like SH-aLRT and aBayes test ([Anisimova et al., 2011](https://doi.org/10.1093/sysbio/syr041)) and tree topology tests like the approximately unbiased (AU) test ([Shimodaira, 2002](https://doi.org/10.1080/10635150290069913)).

The strength of IQ-TREE is the availability of a wide variety of phylogenetic models:

  * **Common models** : All [common substitution models](../../iqtree-model-selection/references/substitution-models.md) for DNA, protein, codon, binary and morphological data with [rate heterogeneity among sites](../../iqtree-model-selection/references/substitution-models.md#rate-heterogeneity-across-sites) and [ascertainment bias correction](../../iqtree-model-selection/references/substitution-models.md#ascertainment-bias-correction) for e.g. SNP data.
  * **[Partition models](../../iqtree-complex-models/references/complex-models.md#partition-models)** : Allowing individual models for different genomic loci (e.g. genes or codon positions), mixed data types, mixed rate heterogeneity types, linked or unlinked branch lengths between partitions.
  * **Mixture models** : [fully customizable mixture models](../../iqtree-complex-models/references/complex-models.md#mixture-models) and [empirical protein mixture models](../../iqtree-model-selection/references/substitution-models.md#protein-mixture-models) and.
  * **Polymorphism-aware models** : Accounting for _incomplete lineage sorting_ to infer species tree from genome-wide population data ([Schrempf et al., 2016](https://doi.org/10.1016/j.jtbi.2016.07.042)).

## Free web server

For a quick start you can also try the IQ-TREE web server, which performs online computation using a dedicated computing cluster. It is very easy to use with as few as just 3 clicks! Try it out at

<http://iqtree.cibiv.univie.ac.at>

## User support

Please refer to the [user documentation](http://www.iqtree.org/doc/) and [frequently asked questions](http://www.iqtree.org/doc/Frequently-Asked-Questions).

If you find a bug (e.g. when IQ-TREE prints a crash message) or want to request a new feature, please post an issue on GitHub: <https://github.com/iqtree/iqtree3/issues>. For other questions and feedback, please ask in GitHub discussions: <https://github.com/iqtree/iqtree3/discussions>

## Documentation

IQ-TREE has an extensive documentation with several tutorials and manual:

  * [Getting started guide](quickstart.md): recommended for users who just downloaded IQ-TREE.
  * [Web Server Tutorial](web-server-tutorial.md): A quick starting guide for the IQ-TREE Web Server.
  * [Beginner’s tutorial](tutorial.md): recommended for users starting to use IQ-TREE.
  * [Advanced tutorial](../../iqtree-complex-models/references/advanced-tutorial.md): recommended for more experienced users who want to explore more features of IQ-TREE.
  * [Command Reference](../../iqtree-command-reference/references/command-reference.md): Comprehensive documentation of command-line options available in IQ-TREE.
  * [Substitution Models](../../iqtree-model-selection/references/substitution-models.md): All common substitution models and usages.
  * [Complex Models](../../iqtree-complex-models/references/complex-models.md): Complex models such as partition and mixture models.
  * [Polymorphism Aware Models](../../iqtree-poms/references/polymorphism-aware-models.md): Polymorphism-aware phylogenetic Models (PoMo) related documentation.
  * [Compilation guide](compilation-guide.md): for advanced users who wants to compile IQ-TREE from source code.
  * [Frequently asked questions (FAQ)](faq.md): recommended to have a look before you post a question in the [IQ-TREE group](https://groups.google.com/d/forum/iqtree).

## How to cite IQ-TREE?

> To maintain IQ-TREE, support users and secure fundings, it is important for us that you cite the following papers, whenever the corresponding features were applied for your analysis.**
>
> Example 1: We obtained branch supports with the ultrafast bootstrap (Hoang et al., 2018) implemented in the IQ-TREE 2 software (Minh et al., 2020).
>
> Example 2: We used IQ-TREE 2 (Minh et al., 2020) to infer the maximum-likelihood tree using the edge-linked partition model (Chernomor et al., 2016).

General citation for IQ-TREE 3:

  * Thomas KF Wong, Nhan Ly-Trong, Huaiyan Ren, Piyumal Demotte, Hector Banos, Andrew J Roger, Edward Susko, Chris Bielow, Nicola De Maio, Nick Goldman, Matthew W Hahn, Mario dos Reis, Le Sy Vinh, Gavin Huttley, Robert Lanfear, Bui Quang Minh (2026) IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models. _Mol. Biol. Eviol._. <https://doi.org/10.1093/molbev/msag117>

When using MixtureFinder please cite:

  * Huaiyan Ren, Thomas KF Wong, Bui Quang Minh, Robert Lanfear (2025) MixtureFinder: Estimating DNA Mixture Models for Phylogenetic Analyses. _Mol. Biol. Evol._ , 42:msae264. <https://doi.org/10.1093/molbev/msae264>

When using GTRpmix please cite:

  * Hector Banos, Thomas KF Wong, Justin Daneau, Edward Susko, Bui Quang Minh, Robert Lanfear, Matthew W Brown, Laura Eme, Andrew J Roger (2025) GTRpmix: A linked general-time reversible model for profile mixture models. _Mol. Biol. Evol._ , 92:msae174. <https://doi.org/10.1093/molbev/msae174>

When using CMAPLE please cite:

  * Nhan Ly-Trong, Chris Bielow, Nicola De Maio, Bui Quang Minh (2024) CMAPLE: Efficient phylogenetic inference in the pandemic era. _Mol. Biol. Evol._ , 41:msae134. <https://doi.org/10.1093/molbev/msae134>

When using the MAST model please cite:

  * Thomas KF Wong, Caitlin Cherryh, Allen G Rodrigo, Matthew W Hahn, Bui Quang Minh, Robert Lanfear (2024) MAST: Phylogenetic Inference with Mixtures Across Sites and Trees. _Syst. Biol._ , 73:375–391. <https://doi.org/10.1093/sysbio/syae008>

When using site/gene concordance factor please cite:

  * Yu K Mo, Robert Lanfear, Matthew W Hahn, Bui Quang Minh (2023) Updated site concordance factors minimize effects of homoplasy and taxon sampling. _Bioinformatics_ , 39:btac741. <https://doi.org/10.1093/bioinformatics/btac741>

When using AliSim please cite:

  * Nhan Ly-Trong, Giuseppe MJ Barca, Bui Quang Minh (2023) AliSim-HPC: parallel sequence simulator for phylogenetics. _Bioinformatics_ , 39:btad540. <https://doi.org/10.1093/bioinformatics/btad540>

When using QMaker for estimating protein Q matrix please cite:

  * Bui Quang Minh, Cuong Cao Dang, Le Sy Vinh, Robert Lanfear (2021) QMaker: Fast and Accurate Method to Estimate Empirical Models of Protein Evolution. _Syst. Biol._ , 70:1046–1060. <https://doi.org/10.1093/sysbio/syab010>

When using nQMaker for estimating non-reversible Q matrix please cite:

  * Cuong Cao Dang, Bui Quang Minh, Hanon McShea, Joanna Masel, Jennifer E James, Le Sy Vinh, Robert Lanfear (2022) nQMaker: Estimating Time Nonreversible Amino Acid Substitution Models. _Syst. Biol._ , 71:1110–1123. <https://doi.org/10.1093/sysbio/syac007>

When using Rootstrap please cite:

  * Suha Naser-Khdour, Bui Quang Minh, Robert Lanfear (2022) Assessing Confidence in Root Placement on Phylogenies: An Empirical Study Using Nonreversible Models for Mammals. _Syst. Biol._ , 71:959–972. <https://doi.org/10.1093/sysbio/syab067>

When using the GHOST "+H" model please cite:

  * Stephen M Crotty, Bui Quang Minh, Nigel G Bean, Barbara R Holland, Jonathan Tuke, Lars S Jermiin, Arndt Von Haeseler (2020) GHOST: Recovering Historical Signal from Heterotachously Evolved Sequence Alignments. _Syst. Biol._ , 69:249-264. <https://doi.org/10.1093/sysbio/syz051>

When using tests of symmetry please cite:

  * Suha Naser-Khdour, Bui Quang Minh, Wenqi Zhang, Eric A Stone, Robert Lanfear (2019) The Prevalence and Impact of Model Violations in Phylogenetic Analysis. _Genome Biol. Evol._ , 11:3341-3352. <https://doi.org/10.1093/gbe/evz193>

When using ModelFinder please cite:

  * Subha Kalyaanamoorthy, Bui Quang Minh, Thomas KF Wong, Arndt von Haeseler, and Lars S Jermiin (2017) ModelFinder: Fast model selection for accurate phylogenetic estimates. _Nat. Methods_ , 14:587–589. <https://doi.org/10.1038/nmeth.4285>

When using polymorphism-aware models please cite:

  * Dominik Schrempf, Bui Quang Minh, Arndt von Haeseler, Carolin Kosiol (2019) Polymorphism-Aware Species Trees with Advanced Mutation Models, Bootstrap, and Rate Heterogeneity. _Mol. Biol. Evol._ , 36:1294–1301. <https://doi.org/10.1093/molbev/msz043>

When using site-specific frequency model please cite:

  * Huai-Chun Wang, Edward Susko, Bui Quang Minh, and Andrew J. Roger (2018) Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. _Syst. Biol._ , 67:216–235. <https://doi.org/10.1093/sysbio/syx068>

When using ultrafast bootstrap (UFBoot) please cite:

  * Diep Thi Hoang, Olga Chernomor, Arndt von Haeseler, Bui Quang Minh, and Le Sy Vinh (2018) UFBoot2: Improving the ultrafast bootstrap approximation. _Mol. Biol. Evol._ , 35:518–522. <https://doi.org/10.1093/molbev/msx281>

When using partition models please cite:

  * Olga Chernomor, Arndt von Haeseler, and Bui Quang Minh (2016) Terrace aware data structure for phylogenomic inference from supermatrices. _Syst. Biol._ , 65:997-1008. <https://doi.org/10.1093/sysbio/syw037>

When using terrace analysis please cite:

  * Biczok R, Bozsoky P, Eisenmann P, Ernst J, Ribizel T, Scholz F, Trefzer A, Weber F, Hamann M, Stamatakis A. (2018) Two C++ libraries for counting trees on a phylogenetic terrace. _Bioinformatics_ , 34:3399–3401. <https://doi.org/10.1093/bioinformatics/bty384>

When using least square dating (LSD) please cite:

  * Thu-Hien To, Matthieu Jung, Samantha Lycett, Olivier Gascuel (2016) Fast dating using least-squares criteria and algorithms. _Syst. Biol._ , 65:82-97. <https://doi.org/10.1093/sysbio/syv068>

When using IQ-TREE web server please cite:

  * J. Trifinopoulos, L.-T. Nguyen, A. von Haeseler, B.Q. Minh (2016) W-IQ-TREE: a fast online phylogenetic tool for maximum likelihood analysis. _Nucleic Acids Res._ , 44:W232-W235. <https://doi.org/10.1093/nar/gkw256>

When using IQ-TREE 2 please cite:

  * B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. _Mol. Biol. Evol._ , 37:1530-1534. <https://doi.org/10.1093/molbev/msaa015>

When using IQ-TREE version 1 please cite:

  * L. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. _Mol. Biol. Evol._ , 32:268-274. <https://doi.org/10.1093/molbev/msu300>

## Development team

IQ-TREE is actively developed by:

**Bui Quang Minh** , _Team leader_ , Designs and implements software core, tree search, ultrafast bootstrap, model selection.

**Robert Lanfear** , _Co-leader_ , Co-leading various projects (incl. new model selection) and documentations, acting as a bridge with biologists.

**Thomas Wong** , _Developer_ , Mixture models, model selection.

**Nhan Ly-Trong** , _Developer_ , sequence simulator, pathogen phylogeny reconstruction.

**Piyumal Demotte** , time tree inference.

**Olga Chernomor** , _Developer_ , Implements partition models.

**Heiko A. Schmidt** , _Developer_ , Integrates TREE-PUZZLE features, user supports and documentations.

**Dominik Schrempf** , _Developer_ , Implements polymorphism-aware models (PoMo).

**Michael Woodhams** , _Developer_ , Implements Lie Markov models.

**Diep Thi Hoang** , _Developer_ , Improves ultrafast bootstrap.

**Arndt von Haeseler** , _Advisor_.

Past members:

**Lam Tung Nguyen** , _Developer_ , Implemented tree search algorithm.

**Jana Trifinopoulos** , _Developer_ , Implemented web service.

**James Barbetti** , _Developer_ , Code optimization.

## Credits and acknowledgements

Some parts of the code were taken from the following packages/libraries: [Phylogenetic likelihood library](http://www.libpll.org), [TREE-PUZZLE](http://www.tree-puzzle.de), [BIONJ](https://doi.org/10.1093/oxfordjournals.molbev.a025808), [Nexus Class Libary](https://doi.org/10.1093/bioinformatics/btg319), [Eigen library](http://eigen.tuxfamily.org/), [SPRNG library](http://www.sprng.org), [Zlib library](http://www.zlib.net), gzstream library, [vectorclass library](http://www.agner.org/optimize/), [GNU scientific library](https://www.gnu.org/software/gsl/).

IQ-TREE was funded by the [Austrian Science Fund](http://www.fwf.ac.at/) (grant no. I760-B17 from 2012-2015 and I 2508-B29 from 2016-2017), the [University of Vienna](https://www.univie.ac.at/) (Initiativkolleg I059-N from 2012-2015), the [Australian National University](https://www.anu.edu.au) (2018-onwards), [Chan-Zuckerberg Initiative](https://chanzuckerberg.com) (2020).
