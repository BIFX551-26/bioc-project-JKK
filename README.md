# RNA-seq analysis with RNAseq123 workflow
Kate, Jessica, Kardam

## Introduction

Goal: Demonstrate the RNAseq123 Bioconductor workflow

Why this dataset?

- Public mouse mammary cell RNA-seq counts from GSE63310
- Includes biologically distinct groups: Basal, LP, and ML
- Small, well-structured example for showing an end-to-end workflow

What the workflow covers:

- data packaging
- filtering and normalization
- quality assessment
- differential expression analysis
- visualization and interpretation

This presentation is a worked example of how raw count data move through
a standard Bioconductor RNA-seq workflow.

![](Index_files/figure-commonmark/unnamed-chunk-3-1.png)

## Experimental Data

- “A pooled shRNA screen for regulators of primary mammary stem and
  progenitor cells identifies roles for Asap1 and Prox1”

- This experiment used RNA-seq-based expression profiling in mouse
  mammary stem cell (MaSC)-enriched basal cells to look for candidate
  regulatory genes.

- Mammary stem cells are often used in RNA-seq analyses as they are
  highly proliferative and have a high capacity for differentiation.

<!-- -->

- Regulatory genes produce proteins that turn expression on or off and
  are therefore important target genes for studying cancer as well as
  other developmental processes.

- Using retroviral aided knockdown experiments and RNA-seq, the
  researchers were able to propose that the Asap1 gene is a negative
  regulator.

- The gene count file used in this analysis was created by aligning
  reads with the mouse reference genome using the align function
  (Rsubread) followed by gene-level summarization using featureCounts.

## Data Packaging

- We started with raw RNA-seq count files from GEO accession GSE63310
  and selected 9 samples representing Basal, LP, and ML mammary cell
  populations.

- Using edgeR, we combined the individual count files into a single
  DGEList object, then added sample metadata such as biological group
  and sequencing lane. We also used Mus.musculus to attach gene
  annotations, including gene symbols and chromosome information.

Data packaging converts separate raw files into one structured,
analysis-ready dataset.

## Data Pre-Processing

- Before comparing gene expression across groups, we transformed raw
  counts to CPM and log-CPM, filtered out genes with very low
  expression, and normalized the libraries using TMM normalization with
  edgeR.

![](Index_files/figure-commonmark/unnamed-chunk-18-1.png)

- We then used MDS plots to check whether samples clustered by biology
  rather than technical artifacts.

- This preprocessing step reduces noise, improves comparability across
  samples, and helps confirm that the data are ready for differential
  expression testing.

Preprocessing removes uninformative genes, corrects library-size
effects, and checks overall sample quality before modeling.
