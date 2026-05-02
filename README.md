# RNA-seq analysis with RNASeq123 workflow
Kate, Jessica, Kardam

Load needed packages:

Load Dependencies for flowchart

## Introduction

RNAseq123 Workflow:

## Experimental Data

“A pooled shRNA screen for regulators of primary mammary stem and
progenitor cells identifies roles for Asap1 and Prox1”

This experiment used RNA-seq-based expression profiling in mouse mammary
stem cell (MaSC)-enriched basal cells to look for candidate regulatory
genes.

Mammary stem cells are often used in RNA-seq analyses as they are highly
proliferative and have a high capacity for differentiation.

Regulatory genes produce proteins that turn expression on or off and are
therefore important target genes for studying cancer as well as other
developmental processes.

Using retroviral aided knockdown experiments and RNA-seq, the
researchers were able to propose that the Asap1 gene is a negative
regulator.

The gene count file used in this analysis was created by aligning reads
with the mouse reference genome using the align function (Rsubread)
followed by gene-level summarization using featureCounts.

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

## Data pre-processing

- Before comparing gene expression across groups, we transformed raw
  counts to CPM and log-CPM, filtered out genes with very low
  expression, and normalized the libraries using TMM normalization with
  edgeR.

- We then used MDS plots to check whether samples clustered by biology
  rather than technical artifacts.

- This preprocessing step reduces noise, improves comparability across
  samples, and helps confirm that the data are ready for differential
  expression testing.

Preprocessing removes uninformative genes, corrects library-size
effects, and checks overall sample quality before modeling.

## Differential Gene Analysis:

- To determine which genes are expressed at different levels between
  three cell population profiled.
- linear model fitting (assuming normally distributed data)
- Intercept act as anchor point for comparision against baseline.
- Without intercept, absolute expression levels can be determined.

## Heteroscedascity and Voom weights:

- Homoscedascity: In linear models, mean-variance relationship is
  assumed to be linear. Mean and variance changes equally. Exactly
  opposite is ‘Heteroscedascity’.

- In RNA-Seq count data, the Negative Binomial distribution assumes
  quadratic mean-variation relationship.

- Overdispersion observed due to large difference between house-keeping
  genes and highly expressed genes

- In limma, linear modelling is carried out on log-CPM values.

- “voom” function act as link to bridge limma(built for microarrays) to
  expression data. It addresses scale problem (calculates into CPM) and
  “normalize.method” argument to normalize library size and estimates
  mean-variance relationship (non-linear - overdispersion)

- “Voom weights” fixes the noise related to genes whereas
  Voomqualityweight addresses the noise generated from inter-sample
  variation.

Means (x-axis) and variances (y-axis) of each gene are plotted to show
the dependence between the two before voom is applied to the data (left
panel) and how the trend is removed after voom precision weights are
applied to the data (right panel)

## Examining Differentially Expressed Genes:

Summary of Differentially expressed gene:

Treat method - Method can be applied for stricter definition on
significance based on t-statistics. This allows user to define a log-FC
threshold.

## Examining individual DE genes from top to bottom:

- “topTreat”/ “topTable” - from toptreat or eBayes to list top DE genes.

- toptreat arranges DE genes chronologically in increasing order using
  log-FC, average log-CPM, moderated t-statistics, raw and adjusted
  p-value for each gene

- n=Inf -\> all genes

## Graphical representation of DE gene results:

- Mean-difference plots are ideally utlized to display log-FC from
  linear model fit against log-CPM values using “plotMD” function.

- Glimma package offer interactive interface for MDplot by “glMDPlot”.

- Glimma option allows brower viewing option - convenient for including
  them as linked files from an Rmarkdown.
