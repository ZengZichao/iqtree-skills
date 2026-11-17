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

#### [__IQ-TREE Tutorials and Manual in one PDF](../assets/files/iqtree-doc.pdf)

## Introduction

  * [Why IQ-TREE?](Home.md#why-iq-tree)
  * [Key features](Home.md#key-features)
  * [Free web server](Home.md#free-web-server)
  * [User support](Home.md#user-support)
  * [Documentation](Home.md#documentation)
  * [How to cite IQ-TREE?](Home.md#how-to-cite-iq-tree)
  * [Development team](Home.md#development-team)
  * [Credits and acknowledgements](Home.md#credits-and-acknowledgements)

[→ Read All Introduction Here](Home.md)

## Getting Started

Recommended for users who just downloaded IQ-TREE the first time.

  * [IQ-TREE web server](Quickstart.md#iq-tree-web-server)
  * [Installation](Quickstart.md#installation)
  * [For Windows users](Quickstart.md#for-windows-users)
  * [For Mac OS X users](Quickstart.md#for-mac-os-x-users)
  * [Minimal command-line examples](Quickstart.md#minimal-command-line-examples)

[→ Read All Getting Started Here](Quickstart.md)

## Web Server Tutorial

A quick starting guide for the IQ-TREE Web Server.

  * [Tree Inference](Web-Server-Tutorial.md#tree-inference)
  * [Model Selection](Web-Server-Tutorial.md#model-selection)
  * [Analysis Results](Web-Server-Tutorial.md#analysis-results)

[→ Read All Web Server Tutorial Here](Web-Server-Tutorial.md)

## Beginner's Tutorial

This tutorial gives a beginner's guide.

  * [Input data](Tutorial.md#input-data)
  * [First example](Tutorial.md#first-running-example)
  * [Model selection](Tutorial.md#choosing-the-right-substitution-model)
  * [Using codon models](Tutorial.md#using-codon-models)
  * [Binary, Morphological, SNPs](Tutorial.md#binary-morphological-and-snp-data)
  * [Ultrafast bootstrap](Tutorial.md#assessing-branch-supports-with-ultrafast-bootstrap-approximation)
  * [Reducing impact of severe model violations with UFBoot](Tutorial.md#reducing-impact-of-severe-model-violations-with-ufboot)
  * [Nonparametric bootstrap](Tutorial.md#assessing-branch-supports-with--standard-nonparametric-bootstrap)
  * [Single branch tests](Tutorial.md#assessing-branch-supports-with-single-branch-tests)
  * [Utilizing multi-core CPUs](Tutorial.md#utilizing-multi-core-cpus)

[→ Read All Beginner's Tutorial Here](Tutorial.md)

## Advanced Tutorial

Recommended for experienced users to explore more features.

  * [Partitioned analysis](Advanced-Tutorial.md#partitioned-analysis-for-multi-gene-alignments)
  * [Partitioning with mixed data](Advanced-Tutorial.md#partitioned-analysis-with-mixed-data)
  * [Partition scheme selection](Advanced-Tutorial.md#choosing-the-right-partitioning-scheme)
  * [Bootstrapping partition model](Advanced-Tutorial.md#ultrafast-bootstrapping-with-partition-model)
  * [Constrained tree search](Advanced-Tutorial.md#constrained-tree-search)
  * [Tree topology tests](Advanced-Tutorial.md#tree-topology-tests)
  * [Testing constrained tree](Advanced-Tutorial.md#testing-constrained-tree)
  * [Consensus construction and bootstrap value assignment](Advanced-Tutorial.md#consensus-construction-and-bootstrap-value-assignment)
  * [User-defined models](Advanced-Tutorial.md#user-defined-substitution-models)
  * [Inferring site-specific rates](Advanced-Tutorial.md#inferring-site-specific-rates)
  * [Trimming alignment sites by likelihood](Advanced-Tutorial.md#trimming-alignment-sites-by-likelihood)

[→ Read All Advanced Tutorial Here](Advanced-Tutorial.md)

## Assessing Phylogenetic Assumptions

This guide is about evaluating the suitability of the data for phylogenetic analysis.

  * [Tests of symmetry](Assessing-Phylogenetic-Assumptions.md#tests-of-symmetry)
  * [Likelihood mapping](Assessing-Phylogenetic-Assumptions.md#likelihood-mapping)

[→ Read All Assessing Phylogenetic Assumptions Here](Assessing-Phylogenetic-Assumptions.md)

## Concordance Factor

Gene and site concordance factor computation.

  * [Inferring species tree](Concordance-Factor.md#inferring-species-tree)
  * [Inferring gene/locus trees](Concordance-Factor.md#inferring-genelocus-trees)
  * [Gene concordance factor (gCF)](Concordance-Factor.md#gene-concordance-factor-gcf)
  * [Site concordance factor (sCF)](Concordance-Factor.md#site-concordance-factor-scf)
  * [Putting it all together](Concordance-Factor.md#putting-it-all-together)

[→ Read All Concordance Factor Here](Concordance-Factor.md)

## Phylogenetic Dating

Inferring time trees

  * [Bayesian dating with MCMCtree](Dating.md#bayesian-dating-with-mcmctree)
  * [IQ2MC workflow for time tree inference](Dating.md#iq2mc-workflow-for-time-tree-inference)
  * [Estimating the gradients and the Hessian for MCMCtree dating](Dating.md#estimating-the-gradients-and-the-hessian-for-mcmctree-dating)
  * [Using partitions and Mixture models for approximate likelihood dating](Dating.md#partitions-and-mixture-models-for-approximate-likelihood-dating)
  * [How to run MCMCtree](Dating.md#how-to-run-mcmctree)
  * [Least Square Dating (LSD2)](Dating.md#least-square-dating)
  * [Inferring time tree with tip dates](Dating.md#inferring-time-tree-with-tip-dates)
  * [Calibrating tree using ancestral dates](Dating.md#calibrating-tree-using-ancestral-dates)
  * [Dating an existing tree](Dating.md#dating-an-existing-tree)
  * [Obtaining confidence intervals](Dating.md#obtaining-confidence-intervals)
  * [Excluding outlier taxa/nodes](Dating.md#excluding-outlier-taxanodes)
  * [Full list of LSD2 options](Dating.md#full-list-of-lsd2-options)

[→ Read All Phylogenetic Dating Here](Dating.md)

## Rooting phylogenies

Inferring rooted trees using non-reversible models.

  * [Inferring unrooted tree with outgroup](Rootstrap.md#inferring-unrooted-tree-with-outgroup)
  * [Inferring rooted trees without outgroup](Rootstrap.md#inferring-rooted-trees-without-outgroup)
  * [Testing root positions](Rootstrap.md#testing-root-positions)

[→ Read All Rooting phylogenies Here](Rootstrap.md)

## Simulating sequence alignments

Simulating sequence alignments

  * [Simulating an alignment from a tree and model](AliSim.md#simulating-an-alignment-from-a-tree-and-model)
  * [Simulating other datatypes](AliSim.md#simulating-other-datatypes)
  * [Non-reversible models](AliSim.md#non-reversible-models)
  * [Rate heterogeneity across sites](AliSim.md#rate-heterogeneity-across-sites)
  * [Customizing output alignments](AliSim.md#customizing-output-alignments)
  * [Insertion and deletion models](AliSim.md#insertion-and-deletion-models)
  * [Specifying model parameters](AliSim.md#specifying-model-parameters)
  * [Mimicking a real alignment](AliSim.md#mimicking-a-real-alignment)
  * [Simulating along a random tree](AliSim.md#simulating-along-a-random-tree)
  * [Branch-specific models](AliSim.md#branch-specific-models)
  * [Partition models](AliSim.md#partition-models)
  * [Mixture models](AliSim.md#mixture-models)
  * [Heterotachy GHOST model](AliSim.md#heterotachy-ghost-model)
  * [Functional divergence model](AliSim.md#functional-divergence-model)
  * [Pre-define mutations](AliSim.md#pre-define-mutations)
  * [Parallel sequence simulations](AliSim.md#parallel-sequence-simulations)
  * [Command reference](AliSim.md#command-reference)

[→ Read All Simulating sequence alignments Here](AliSim.md)

## Frequently Asked Questions

For common questions and answers.

  * [How do I get help?](Frequently-Asked-Questions.md#how-do-i-get-help)
  * [How do I report bug?](Frequently-Asked-Questions.md#how-do-i-report-bug)
  * [How to interpret ultrafast bootstrap (UFBoot) supports?](Frequently-Asked-Questions.md#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values)
  * [How does IQ-TREE treat gap/missing/ambiguous characters?](Frequently-Asked-Questions.md#how-does-iq-tree-treat-gapmissingambiguous-characters)
  * [Can I mix DNA and protein data in a partitioned analysis?](Frequently-Asked-Questions.md#can-i-mix-dna-and-protein-data-in-a-partitioned-analysis)
  * [What is the interpretation of branch lengths when mixing codon and DNA data?](Frequently-Asked-Questions.md#what-is-the-interpretation-of-branch-lengths-when-mixing-codon-and-dna-data)
  * [What is the purpose of composition test?](Frequently-Asked-Questions.md#what-is-the-purpose-of-composition-test)
  * [What is the good number of CPU cores to use](Frequently-Asked-Questions.md#what-is-the-good-number-of-cpu-cores-to-use)
  * [How do I save time for standard bootstrap?](Frequently-Asked-Questions.md#how-do-i-save-time-for-standard-bootstrap)
  * [Why does IQ-TREE complain about the use of +ASC model?](Frequently-Asked-Questions.md#why-does-iq-tree-complain-about-the-use-of-asc-model)
  * [How does IQ-TREE treat identical sequences?](Frequently-Asked-Questions.md#how-does-iq-tree-treat-identical-sequences)

[→ Read All Frequently Asked Questions Here](Frequently-Asked-Questions.md)

## Recipes (How-To)

Collection of examples that combine IQ-TREE with other software.

  * [Phylogenetic placement](recipes/phylogenetic-placement.md)
  * [Perform a Goldman-Cox test](recipes/goldman-cox-test.md)
  * [Estimating gene, site, and quartet concordance vectors](recipes/concordance-vector.md)
  * [Simulating concatenated alignments with multiple trees](recipes/Simulating-concatenated-alignments-with-multiple-trees.md)

[→ Read All Recipes (How-To) Here](recipes/index.md)

## Substitution Models

All common substitution models and usages.

  * [DNA models](Substitution-Models.md#dna-models)
  * [Protein models](Substitution-Models.md#protein-models)
  * [Codon models](Substitution-Models.md#codon-models)
  * [Binary, morphological models](Substitution-Models.md#binary-and-morphological-models)
  * [Ascertainment bias correction](Substitution-Models.md#ascertainment-bias-correction)
  * [Rate heterogeneity across sites](Substitution-Models.md#rate-heterogeneity-across-sites)

[→ Read All Substitution Models Here](Substitution-Models.md)

## Complex Models

Complex models such as partition and mixture models.

  * [Partition models](Complex-Models.md#partition-models)
  * [Mixture models](Complex-Models.md#mixture-models)
  * [Site-specific frequency models](Complex-Models.md#site-specific-frequency-models)
  * [Heterotachy models](Complex-Models.md#heterotachy-models)
  * [Multitree models](Complex-Models.md#multitree-models)

[→ Read All Complex Models Here](Complex-Models.md)

## Estimating amino acid substitution models

Estimating amino acid substitution models.

  * [Estimating a model from a single concatenated alignment](Estimating-amino-acid-substitution-models.md#estimating-a-model-from-a-single-concatenated-alignment)
  * [Estimating a model from a folder of alignments](Estimating-amino-acid-substitution-models.md#estimating-a-model-from-a-folder-of-alignments)
  * [Estimating a non-reversible model](Estimating-amino-acid-substitution-models.md#estimating-a-non-reversible-model)
  * [Estimating linked exchangeabilities](Estimating-amino-acid-substitution-models.md#estimating-linked-exchangeabilities)
  * [Estimating a multiple exchangeability and profile model](Estimating-amino-acid-substitution-models.md#estimating-a-multiple-exchangeability-and-profile-model)

[→ Read All Estimating amino acid substitution models Here](Estimating-amino-acid-substitution-models.md)

## Polymorphism-Aware Models

Improve inferences using population data.

  * [Counts files](Polymorphism-Aware-Models.md#counts-files)
  * [First example](Polymorphism-Aware-Models.md#first-running-example)
  * [Substitution models](Polymorphism-Aware-Models.md#substitution-models)
  * [Virtual population size](Polymorphism-Aware-Models.md#virtual-population-size)
  * [Sampling method](Polymorphism-Aware-Models.md#sampling-method)
  * [Bootstrap branch support](Polymorphism-Aware-Models.md#bootstrap-branch-support)
  * [Interpretation of branch lengths](Polymorphism-Aware-Models.md#interpretation-of-branch-lengths)

[→ Read All Polymorphism-Aware Models Here](Polymorphism-Aware-Models.md)

## Command Reference

Commprehensive documentation of command-line options.

  * [General options](Command-Reference.md#general-options)
  * [Checkpointing to resume stopped run](Command-Reference.md#checkpointing-to-resume-stopped-run)
  * [Likelihood mapping analysis](Command-Reference.md#likelihood-mapping-analysis)
  * [Automatic model selection](Command-Reference.md#automatic-model-selection)
  * [Specifying substitution models](Command-Reference.md#specifying-substitution-models)
  * [Rate heterogeneity](Command-Reference.md#rate-heterogeneity)
  * [Partition model options](Command-Reference.md#partition-model-options)
  * [Site-specific frequency model options](Command-Reference.md#site-specific-frequency-model-options)
  * [Tree search parameters](Command-Reference.md#tree-search-parameters)
  * [Tree search for pathogen data](Command-Reference.md#tree-search-for-pathogen-data)
  * [Ultrafast bootstrap parameters](Command-Reference.md#ultrafast-bootstrap-parameters)
  * [Nonparametric bootstrap](Command-Reference.md#nonparametric-bootstrap)
  * [Single branch tests](Command-Reference.md#single-branch-tests)
  * [Ancestral sequence reconstruction](Command-Reference.md#ancestral-sequence-reconstruction)
  * [Tree topology tests](Command-Reference.md#tree-topology-tests)
  * [Constructing consensus tree](Command-Reference.md#constructing-consensus-tree)
  * [Robinson-Foulds distance](Command-Reference.md#computing-robinson-foulds-distance)
  * [Generating random trees](Command-Reference.md#generating-random-trees)
  * [Miscellaneous options](Command-Reference.md#miscellaneous-options)

[→ Read All Command Reference Here](Command-Reference.md)

## Compilation Guide

For advanced users to compile IQ-TREE source code.

  * [General requirements](Compilation-Guide.md#general-requirements)
  * [Downloading source code](Compilation-Guide.md#downloading-source-code)
  * [Compiling under Linux](Compilation-Guide.md#compiling-under-linux)
  * [Compiling under Mac OS X](Compilation-Guide.md#compiling-under-mac-os-x)
  * [Compiling under Windows](Compilation-Guide.md#compiling-under-windows)
  * [Compiling 32-bit version](Compilation-Guide.md#compiling-32-bit-version)
  * [Compiling MPI version](Compilation-Guide.md#compiling-mpi-version)
  * [Compiling Xeon Phi Knights Landing version](Compilation-Guide.md#compiling-xeon-phi-knights-landing-version)
  * [About precompiled binaries](Compilation-Guide.md#about-precompiled-binaries)

[→ Read All Compilation Guide Here](Compilation-Guide.md)

## Developer Guide

This guide gives developers an overview of IQ-TREE software design, data structures and discusses possibility to incorporate new models into IQ-TREE.

  * [Alignment](Developer-Guide.md#alignment-class)
  * [Model of substitution](Developer-Guide.md#model-of-substitution)
  * [PhyloTree class](Developer-Guide.md#phylotree-class-phylogenetic-tree)

[→ Read All Developer Guide Here](Developer-Guide.md)

## Analyzing Big Data

Hints and strategies to analyze big alignments with >= 1000 sequences or >= 10,000 sites.

[→ Read All Analyzing Big Data Here](Analyzing-Big-Data.md)
