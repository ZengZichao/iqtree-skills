---
title: Introduction
source_url: https://iqtree.github.io/doc/index.html
headings: 20
---

# Introduction

> Source: <https://iqtree.github.io/doc/index.html>

## Table of contents

- [Introduction](#introduction)
- [Getting Started](#getting-started)
- [Web Server Tutorial](#web-server-tutorial)
- [Beginner's Tutorial](#beginners-tutorial)
- [Advanced Tutorial](#advanced-tutorial)
- [Assessing Phylogenetic Assumptions](#assessing-phylogenetic-assumptions)
- [Concordance Factor](#concordance-factor)
- [Phylogenetic Dating](#phylogenetic-dating)
- [Rooting phylogenies](#rooting-phylogenies)
- [Simulating sequence alignments](#simulating-sequence-alignments)
- [Frequently Asked Questions](#frequently-asked-questions)
- [Recipes (How-To)](#recipes-how-to)
- [Substitution Models](#substitution-models)
- [Complex Models](#complex-models)
- [Estimating amino acid substitution models](#estimating-amino-acid-substitution-models)
- [Polymorphism-Aware Models](#polymorphism-aware-models)
- [Command Reference](#command-reference)
- [Compilation Guide](#compilation-guide)
- [Developer Guide](#developer-guide)
- [Analyzing Big Data](#analyzing-big-data)

#### [__IQ-TREE Tutorials and Manual in one PDF](../../assets/files/iqtree-doc.pdf)

## Introduction

  * [Why IQ-TREE?](home.md#why-iq-tree)
  * [Key features](home.md#key-features)
  * [Free web server](home.md#free-web-server)
  * [User support](home.md#user-support)
  * [Documentation](home.md#documentation)
  * [How to cite IQ-TREE?](home.md#how-to-cite-iq-tree)
  * [Development team](home.md#development-team)
  * [Credits and acknowledgements](home.md#credits-and-acknowledgements)

[→ Read All Introduction Here](home.md)

## Getting Started

Recommended for users who just downloaded IQ-TREE the first time.

  * [IQ-TREE web server](quickstart.md#iq-tree-web-server)
  * [Installation](quickstart.md#installation)
  * [For Windows users](quickstart.md#for-windows-users)
  * [For Mac OS X users](quickstart.md#for-mac-os-x-users)
  * [Minimal command-line examples](quickstart.md#minimal-command-line-examples)

[→ Read All Getting Started Here](quickstart.md)

## Web Server Tutorial

A quick starting guide for the IQ-TREE Web Server.

  * [Tree Inference](web-server-tutorial.md#tree-inference)
  * [Model Selection](web-server-tutorial.md#model-selection)
  * [Analysis Results](web-server-tutorial.md#analysis-results)

[→ Read All Web Server Tutorial Here](web-server-tutorial.md)

## Beginner's Tutorial

This tutorial gives a beginner's guide.

  * [Input data](tutorial.md#input-data)
  * [First example](tutorial.md#first-running-example)
  * [Model selection](tutorial.md#choosing-the-right-substitution-model)
  * [Using codon models](tutorial.md#using-codon-models)
  * [Binary, Morphological, SNPs](tutorial.md#binary-morphological-and-snp-data)
  * [Ultrafast bootstrap](tutorial.md#assessing-branch-supports-with-ultrafast-bootstrap-approximation)
  * [Reducing impact of severe model violations with UFBoot](tutorial.md#reducing-impact-of-severe-model-violations-with-ufboot)
  * [Nonparametric bootstrap](tutorial.md#assessing-branch-supports-with-standard-nonparametric-bootstrap)
  * [Single branch tests](tutorial.md#assessing-branch-supports-with-single-branch-tests)
  * [Utilizing multi-core CPUs](tutorial.md#utilizing-multi-core-cpus)

[→ Read All Beginner's Tutorial Here](tutorial.md)

## Advanced Tutorial

Recommended for experienced users to explore more features.

  * [Partitioned analysis](../../iqtree-complex-models/references/advanced-tutorial.md#partitioned-analysis-for-multi-gene-alignments)
  * [Partitioning with mixed data](../../iqtree-complex-models/references/advanced-tutorial.md#partitioned-analysis-with-mixed-data)
  * [Partition scheme selection](../../iqtree-complex-models/references/advanced-tutorial.md#choosing-the-right-partitioning-scheme)
  * [Bootstrapping partition model](../../iqtree-complex-models/references/advanced-tutorial.md#ultrafast-bootstrapping-with-partition-model)
  * [Constrained tree search](../../iqtree-complex-models/references/advanced-tutorial.md#constrained-tree-search)
  * [Tree topology tests](../../iqtree-complex-models/references/advanced-tutorial.md#tree-topology-tests)
  * [Testing constrained tree](../../iqtree-complex-models/references/advanced-tutorial.md#testing-constrained-tree)
  * [Consensus construction and bootstrap value assignment](../../iqtree-complex-models/references/advanced-tutorial.md#consensus-construction-and-bootstrap-value-assignment)
  * [User-defined models](../../iqtree-complex-models/references/advanced-tutorial.md#user-defined-substitution-models)
  * [Inferring site-specific rates](../../iqtree-complex-models/references/advanced-tutorial.md#inferring-site-specific-rates)
  * [Trimming alignment sites by likelihood](../../iqtree-complex-models/references/advanced-tutorial.md#trimming-alignment-sites-by-likelihood)

[→ Read All Advanced Tutorial Here](../../iqtree-complex-models/references/advanced-tutorial.md)

## Assessing Phylogenetic Assumptions

This guide is about evaluating the suitability of the data for phylogenetic analysis.

  * [Tests of symmetry](../../iqtree-branch-support/references/assessing-assumptions.md#tests-of-symmetry)
  * [Likelihood mapping](../../iqtree-branch-support/references/assessing-assumptions.md#likelihood-mapping)

[→ Read All Assessing Phylogenetic Assumptions Here](../../iqtree-branch-support/references/assessing-assumptions.md)

## Concordance Factor

Gene and site concordance factor computation.

  * [Inferring species tree](../../iqtree-concordance-factor/references/concordance-factor.md#inferring-species-tree)
  * [Inferring gene/locus trees](../../iqtree-concordance-factor/references/concordance-factor.md#inferring-genelocus-trees)
  * [Gene concordance factor (gCF)](../../iqtree-concordance-factor/references/concordance-factor.md#gene-concordance-factor-gcf)
  * [Site concordance factor (sCF)](../../iqtree-concordance-factor/references/concordance-factor.md#site-concordance-factor-scf)
  * [Putting it all together](../../iqtree-concordance-factor/references/concordance-factor.md#putting-it-all-together)

[→ Read All Concordance Factor Here](../../iqtree-concordance-factor/references/concordance-factor.md)

## Phylogenetic Dating

Inferring time trees

  * [Bayesian dating with MCMCtree](../../iqtree-dating/references/dating.md#bayesian-dating-with-mcmctree)
  * [IQ2MC workflow for time tree inference](../../iqtree-dating/references/dating.md#iq2mc-workflow-for-time-tree-inference)
  * [Estimating the gradients and the Hessian for MCMCtree dating](../../iqtree-dating/references/dating.md#estimating-the-gradients-and-the-hessian-for-mcmctree-dating)
  * [Using partitions and Mixture models for approximate likelihood dating](../../iqtree-dating/references/dating.md#using-partitions-and-mixture-models-for-approximate-likelihood-dating)
  * [How to run MCMCtree](../../iqtree-dating/references/dating.md#how-to-run-mcmctree)
  * [Least Square Dating (LSD2)](../../iqtree-dating/references/dating.md#least-square-dating-lsd2)
  * [Inferring time tree with tip dates](../../iqtree-dating/references/dating.md#inferring-time-tree-with-tip-dates)
  * [Calibrating tree using ancestral dates](../../iqtree-dating/references/dating.md#calibrating-tree-using-ancestral-dates)
  * [Dating an existing tree](../../iqtree-dating/references/dating.md#dating-an-existing-tree)
  * [Obtaining confidence intervals](../../iqtree-dating/references/dating.md#obtaining-confidence-intervals)
  * [Excluding outlier taxa/nodes](../../iqtree-dating/references/dating.md#excluding-outlier-taxanodes)
  * [Full list of LSD2 options](../../iqtree-dating/references/dating.md#full-list-of-lsd2-options)

[→ Read All Phylogenetic Dating Here](../../iqtree-dating/references/dating.md)

## Rooting phylogenies

Inferring rooted trees using non-reversible models.

  * [Inferring unrooted tree with outgroup](../../iqtree-rooting/references/rootstrap.md#inferring-unrooted-tree-with-outgroup)
  * [Inferring rooted trees without outgroup](../../iqtree-rooting/references/rootstrap.md#inferring-rooted-trees-without-outgroup)
  * [Testing root positions](../../iqtree-rooting/references/rootstrap.md#testing-root-positions)

[→ Read All Rooting phylogenies Here](../../iqtree-rooting/references/rootstrap.md)

## Simulating sequence alignments

Simulating sequence alignments

  * [Simulating an alignment from a tree and model](../../iqtree-alisim/references/alisim.md#simulating-an-alignment-from-a-tree-and-model)
  * [Simulating other datatypes](../../iqtree-alisim/references/alisim.md#simulating-other-datatypes)
  * [Non-reversible models](../../iqtree-alisim/references/alisim.md#non-reversible-models)
  * [Rate heterogeneity across sites](../../iqtree-alisim/references/alisim.md#rate-heterogeneity-across-sites)
  * [Customizing output alignments](../../iqtree-alisim/references/alisim.md#customizing-output-alignments)
  * [Insertion and deletion models](../../iqtree-alisim/references/alisim.md#insertion-and-deletion-models)
  * [Specifying model parameters](../../iqtree-alisim/references/alisim.md#specifying-model-parameters)
  * [Mimicking a real alignment](../../iqtree-alisim/references/alisim.md#mimicking-a-real-alignment)
  * [Simulating along a random tree](../../iqtree-alisim/references/alisim.md#simulating-along-a-random-tree)
  * [Branch-specific models](../../iqtree-alisim/references/alisim.md#branch-specific-models)
  * [Partition models](../../iqtree-alisim/references/alisim.md#partition-models)
  * [Mixture models](../../iqtree-alisim/references/alisim.md#mixture-models)
  * [Heterotachy GHOST model](../../iqtree-alisim/references/alisim.md#heterotachy-ghost-model)
  * [Functional divergence model](../../iqtree-alisim/references/alisim.md#functional-divergence-model)
  * [Pre-define mutations](../../iqtree-alisim/references/alisim.md#pre-define-mutations)
  * [Parallel sequence simulations](../../iqtree-alisim/references/alisim.md#parallel-sequence-simulations)
  * [Command reference](../../iqtree-alisim/references/alisim.md#all-alisim-options)

[→ Read All Simulating sequence alignments Here](../../iqtree-alisim/references/alisim.md)

## Frequently Asked Questions

For common questions and answers.

  * [How do I get help?](faq.md#how-do-i-get-help)
  * [How do I report bug?](faq.md#how-do-i-report-a-bug)
  * [How to interpret ultrafast bootstrap (UFBoot) supports?](faq.md#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values)
  * [How does IQ-TREE treat gap/missing/ambiguous characters?](faq.md#how-does-iq-tree-treat-gapmissingambiguous-characters)
  * [Can I mix DNA and protein data in a partitioned analysis?](faq.md#can-i-mix-dna-and-protein-data-in-a-partitioned-analysis)
  * [What is the interpretation of branch lengths when mixing codon and DNA data?](faq.md#what-is-the-interpretation-of-branch-lengths-when-mixing-codon-and-dna-data)
  * [What is the purpose of composition test?](faq.md#what-is-the-purpose-of-composition-test)
  * [What is the good number of CPU cores to use](faq.md#what-is-the-good-number-of-cpu-cores-to-use)
  * [How do I save time for standard bootstrap?](faq.md#how-do-i-save-time-for-standard-bootstrap)
  * [Why does IQ-TREE complain about the use of +ASC model?](faq.md#why-does-iq-tree-complain-about-the-use-of-asc-model)
  * [How does IQ-TREE treat identical sequences?](faq.md#how-does-iq-tree-treat-identical-sequences)

[→ Read All Frequently Asked Questions Here](faq.md)

## Recipes (How-To)

Collection of examples that combine IQ-TREE with other software.

  * [Phylogenetic placement](../../iqtree-recipes/references/phylogenetic-placement.md)
  * [Perform a Goldman-Cox test](../../iqtree-recipes/references/goldman-cox-test.md)
  * [Estimating gene, site, and quartet concordance vectors](../../iqtree-recipes/references/concordance-vector.md)
  * [Simulating concatenated alignments with multiple trees](../../iqtree-recipes/references/Simulating-concatenated-alignments-with-multiple-trees.md)

[→ Read All Recipes (How-To) Here](../../iqtree-recipes/references/index.md)

## Substitution Models

All common substitution models and usages.

  * [DNA models](../../iqtree-model-selection/references/substitution-models.md#dna-models)
  * [Protein models](../../iqtree-model-selection/references/substitution-models.md#protein-models)
  * [Codon models](../../iqtree-model-selection/references/substitution-models.md#codon-models)
  * [Binary, morphological models](../../iqtree-model-selection/references/substitution-models.md#binary-and-morphological-models)
  * [Ascertainment bias correction](../../iqtree-model-selection/references/substitution-models.md#ascertainment-bias-correction)
  * [Rate heterogeneity across sites](../../iqtree-model-selection/references/substitution-models.md#rate-heterogeneity-across-sites)

[→ Read All Substitution Models Here](../../iqtree-model-selection/references/substitution-models.md)

## Complex Models

Complex models such as partition and mixture models.

  * [Partition models](../../iqtree-complex-models/references/complex-models.md#partition-models)
  * [Mixture models](../../iqtree-complex-models/references/complex-models.md#mixture-models)
  * [Site-specific frequency models](../../iqtree-complex-models/references/complex-models.md#site-specific-frequency-models)
  * [Heterotachy models](../../iqtree-complex-models/references/complex-models.md#heterotachy-models)
  * [Multitree models](../../iqtree-complex-models/references/complex-models.md#multitree-models)

[→ Read All Complex Models Here](../../iqtree-complex-models/references/complex-models.md)

## Estimating amino acid substitution models

Estimating amino acid substitution models.

  * [Estimating a model from a single concatenated alignment](../../iqtree-model-selection/references/estimating-aa-models.md#estimating-a-model-from-a-single-concatenated-alignment)
  * [Estimating a model from a folder of alignments](../../iqtree-model-selection/references/estimating-aa-models.md#estimating-a-model-from-a-folder-of-alignments)
  * [Estimating a non-reversible model](../../iqtree-model-selection/references/estimating-aa-models.md#estimating-a-non-reversible-model)
  * [Estimating linked exchangeabilities](../../iqtree-model-selection/references/estimating-aa-models.md#estimating-linked-exchangeabilities)
  * [Estimating a multiple exchangeability and profile model](../../iqtree-model-selection/references/estimating-aa-models.md#estimating-a-multiple-exchangeability-and-profile-model)

[→ Read All Estimating amino acid substitution models Here](../../iqtree-model-selection/references/estimating-aa-models.md)

## Polymorphism-Aware Models

Improve inferences using population data.

  * [Counts files](../../iqtree-poms/references/polymorphism-aware-models.md#counts-files)
  * [First example](../../iqtree-poms/references/polymorphism-aware-models.md#first-running-example)
  * [Substitution models](../../iqtree-poms/references/polymorphism-aware-models.md#substitution-models)
  * [Virtual population size](../../iqtree-poms/references/polymorphism-aware-models.md#virtual-population-size)
  * [Sampling method](../../iqtree-poms/references/polymorphism-aware-models.md#sampling-method)
  * [Bootstrap branch support](../../iqtree-poms/references/polymorphism-aware-models.md#bootstrap-branch-support)
  * [Interpretation of branch lengths](../../iqtree-poms/references/polymorphism-aware-models.md#interpretation-of-branch-lengths)

[→ Read All Polymorphism-Aware Models Here](../../iqtree-poms/references/polymorphism-aware-models.md)

## Command Reference

Commprehensive documentation of command-line options.

  * [General options](../../iqtree-command-reference/references/command-reference.md#general-options)
  * [Checkpointing to resume stopped run](../../iqtree-command-reference/references/command-reference.md#checkpointing-to-resume-stopped-run)
  * [Likelihood mapping analysis](../../iqtree-command-reference/references/command-reference.md#likelihood-mapping-analysis)
  * [Automatic model selection](../../iqtree-command-reference/references/command-reference.md#automatic-model-selection)
  * [Specifying substitution models](../../iqtree-command-reference/references/command-reference.md#specifying-substitution-models)
  * [Rate heterogeneity](../../iqtree-command-reference/references/command-reference.md#rate-heterogeneity)
  * [Partition model options](../../iqtree-command-reference/references/command-reference.md#partition-model-options)
  * [Site-specific frequency model options](../../iqtree-command-reference/references/command-reference.md#site-specific-frequency-model-options)
  * [Tree search parameters](../../iqtree-command-reference/references/command-reference.md#tree-search-parameters)
  * [Tree search for pathogen data](../../iqtree-command-reference/references/command-reference.md#tree-search-for-pathogen-data)
  * [Ultrafast bootstrap parameters](../../iqtree-command-reference/references/command-reference.md#ultrafast-bootstrap-parameters)
  * [Nonparametric bootstrap](../../iqtree-command-reference/references/command-reference.md#nonparametric-bootstrap)
  * [Single branch tests](../../iqtree-command-reference/references/command-reference.md#single-branch-tests)
  * [Ancestral sequence reconstruction](../../iqtree-command-reference/references/command-reference.md#ancestral-sequence-reconstruction)
  * [Tree topology tests](../../iqtree-command-reference/references/command-reference.md#tree-topology-tests)
  * [Constructing consensus tree](../../iqtree-command-reference/references/command-reference.md#constructing-consensus-tree)
  * [Robinson-Foulds distance](../../iqtree-command-reference/references/command-reference.md#computing-robinson-foulds-distance)
  * [Generating random trees](../../iqtree-command-reference/references/command-reference.md#generating-random-trees)
  * [Miscellaneous options](../../iqtree-command-reference/references/command-reference.md#miscellaneous-options)

[→ Read All Command Reference Here](../../iqtree-command-reference/references/command-reference.md)

## Compilation Guide

For advanced users to compile IQ-TREE source code.

  * [General requirements](compilation-guide.md#general-requirements)
  * [Downloading source code](compilation-guide.md#downloading-source-code)
  * [Compiling under Linux](compilation-guide.md#compiling-under-linux)
  * [Compiling under Mac OS X](compilation-guide.md#compiling-under-mac-os-x)
  * [Compiling under Windows](compilation-guide.md#compiling-under-windows)
  * [Compiling 32-bit version](compilation-guide.md#compiling-32-bit-version)
  * [Compiling MPI version](compilation-guide.md#compiling-mpi-version)
  * [Compiling Xeon Phi Knights Landing version](compilation-guide.md#compiling-xeon-phi-knights-landing-version)
  * [About precompiled binaries](compilation-guide.md#about-precompiled-binaries)

[→ Read All Compilation Guide Here](compilation-guide.md)

## Developer Guide

This guide gives developers an overview of IQ-TREE software design, data structures and discusses possibility to incorporate new models into IQ-TREE.

  * [Alignment](developer-guide.md#alignment-class)
  * [Model of substitution](developer-guide.md#model-of-substitution)
  * [PhyloTree class](developer-guide.md#phylotree-class-phylogenetic-tree)

[→ Read All Developer Guide Here](developer-guide.md)

## Analyzing Big Data

Hints and strategies to analyze big alignments with >= 1000 sequences or >= 10,000 sites.

[→ Read All Analyzing Big Data Here](analyzing-big-data.md)
