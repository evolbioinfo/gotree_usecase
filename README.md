# This repository contains the Gotree/Goalign use-case workflow


## Introduction
In this use case, we analyze a phylogenomic dataset inspired from [Vanderpool et al. PLoS biology, 2020](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000954&rev=1#pbio.3000954.ref039), in which the authors analyze a set of 1,730 genes in primates in different ways. They infer the species tree either from from individual gene trees using [ASTRAL III](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) or from gene concatenation using maximum likelihood. Our use case is inspired from the concatenation study, using available groups of primate orthologous proteins in [OrthoDB](https://www.orthodb.org/).

To do so, the workflow first maps the genbank identifiers of the 1,730 analyzed genes to their OrthoDB identifiers, and retrieves the orthologous groups of proteins shared in at least 90% of the 25 analyzed primates, and finally reconstructs a phylogenetic tree of these primates.

## Results

The tree can be visualized [here](https://itol.embl.de/tree/157996425332211619518043)

![Tree](images/final_tree.svg)

An archive with alignment, gene ids, and trees is downloadable on the [v1.0 release](https://github.com/evolbioinfo/gotree_usecase/releases/tag/v1.0).

## Running the workflow

- Pre-requisites: The workflow only needs [Singularity](https://sylabs.io/) and [Nextflow](https://www.nextflow.io/)
- Configuring the workflow: Change values of `executor`, `queue` and `clusterOptions` in `nextflow.config`
- Running the workflow: `nextflow run workflow.nf --itolkey <User iTOL key> --itolproject <iTOL project>`

