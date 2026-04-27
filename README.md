RNA123 Workflow for RNA-seq Analysis

Introduction:



Experimental Data:

The article describes a workflow to set up a matrix for a single gene and then applied to remaining genes in any dataset using limma packages for differential expression analysis. Covariates and Factors are two types of explanatory variables used for analysis. Covariates are used for quantitative measurements (for ex. Age, height) and factors are categorical variables(for ex. Genotype with two levels: wildtype or mutant). This workflow utilizes edgeR-limma as foundation to take gene-level counts as its input, processes and allows exploratory analysis before obtaining DE genes and gene signatures. Interactive graphics are obtained from the Glimma package for detailed exploration of data both sample and gene-level.
Examples are coded based on the assumption that expression values are genewise log-count per million measurements from an RNA-sequencing experimentation.



Analysis Method:

Overview of the analysis methods - “What” did you do?
1.Setting up the library
2.Reading & Extracting count data
3.Organizing the sample information  - Labelling the raw count matrix with metadata  for ex. Basal, ML etc
4.Organizing gene annotations - Gene information such as gene symbol, gene names, chromosome names and location, IDs etc retrieval using Organism specific packages i.e. Mus Musculus or biomaRt package for retrieving from Ensembl genome database.
5.Data transformation:
Unsupervised learning: Using MDS plot for exploratory analysis aimed at determining actual similarities and differences to further categorize the outlier and the sources of error or variation.



Conclusion:

1.Matrix and metadata are stored in two separate files which need to be patched up after loading matrix.

2. filterByexp - Function that removes the genes with low expression    values. However, it does not filter the control whose expression value could be near 0 or 0

3. calcNormfactor using TMM - A function works on logic where over & under expressed genes are ignored initially to plot a baseline. Once the baseline is plotted, the genes are cross-references for accurate comparison to avoid the compositional bias.

4. Voom weights vs Voomqualityweights - The first one fixes the noise related to gene size
