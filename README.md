# matched GWAS

## Overview  
The R package `matchedGWAS` is developed to conduct Genome-Wide Association Study (GWAS) with paired samples, where cases are genotyped from tumor tissues and controls come from normal tissues of the same patients. Most traditional GWAS  focus on independent case-control design. To address the limitation of current GWAS, this method is proposed to accommodate the structure of paired data to identify genetic variants or DNA segments that are associated with cancers.

Bayesian hierarchical methods for **tumor–normal matched GWAS** with multi-marker (segment/gene-level) aggregation.

This repository implements the MCMC sampler described in the dissertation:
- **Hierarchical model** with gene/segment status $$G$$ and SNP status $$H_j$$
- Per-SNP parameters: **RR** (relative risk), **AF** (allele frequency), **MR** (mutation rate)
- **Gibbs** sampler with **Metropolis–Hastings** updates for $$(RR, AF, MR)$$
- Support for **burn-in**, **thinning**, **multiple chains**, and **acceptance-rate** diagnostics


## Methods

## Data

The algorithm expects **tumor–normal matched genotype counts** in a 3×3 table per SNP, summarized as:

| Column | Description |
|:--------|:-------------|
| `segment_id` | Segment or gene identifier |
| `snp_id` | SNP identifier |
| `n00, n01, n02` | Counts of genotypes 0, 1, 2 in **normal** tissue |
| `t00, t01, t02` | Counts of genotypes 0, 1, 2 in **tumor** tissue |
| `AF`, `MR`, `RR_true` | Optional simulated or annotation columns |

Segmentation is recommended for large genes (>32 kb or >40 SNPs), consistent with the dissertation’s **segmentation step in Section 4.1.2**.

## Installation

install.packages(c("data.table", "tidyverse", "coda", "matrixStats"))
