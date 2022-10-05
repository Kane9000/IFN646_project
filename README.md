# IFN646_project
Introduction

The introduction of high-throughput RNA sequencing technologies has allowed us to understand the nature of gene expression. 
These experiments would involve the quantification of a large number of short reads which are then aligned to a genome to create a count matrix of each gene or feature. 
RNA differential expression tools allow the user to quickly find genes that are highly expressed or suppressed between two conditions.

This project aims to compare tools that analyze differential expression of RNA sequencing data. 
Synthetic count matrices were generated and provided for the comparison which included two conditional groups, condition 1 and condition 2, with varying sample sizes from 3 to 9 per group and approximately ten thousand genes. 
A separate meta file was also provided that indicates whether condition 2 genes were either upregulated, downregulated or differentially expressed compared to condition 1.

Method

The differential expression tools we will investigate in this study include DESeq2, EdgeR and Voom limma. 
The tools will be evaluated based on their ability to determine whether genes were differential expressed or not as indicated in the meta file. 

To compare the ability to rank differentially expressed genes, we will generate receiver operating characteristic (ROC) graphs which plot true positive and false positive rates and evaluate each area under ROC curve (AUC) which calculates the performance of the model across all probability thresholds.
RNA differential expression tools are built in the R language (version 4.2.1) and scripts were written in RStudio (version 2022.07.1 build 554). Results from these tools were analysed in Jupyter notebook (version 5.7.8)running Python version 3.7.3). 
A logistic regression was trained using the output of the differential expression tool as the response variable and a meta data file as the predictor variable. Genes with a statistically significant fold change of p value=0.05 or smaller were considered differentially expressed.

Results

![alt text](/Picture500_500.png?raw=true)
![alt text](/Picture750_250.png?raw=true)
![alt text](/Picture1000_0.png?raw=true)
 
Area under the ROC curve (AUC). Area under the ROC curve for DESeq2 and EdgeR using synthetic count data.
Each file name indicates the number of genes that are upregulated and downregulated compared to in condition 2 compared to condition 1. 
Sample size varies from n = 3, 6 or 9 per group.
Each point represents an AUC score that has been generated from either DESeq2 (circle), EdgeR (square) or VoomLimma (Diamond).

AUC score was increased in all scenarios when the sample size was increased. DESeq2 has a slight advantage over EdgeR in most scenarios. EdgeR ranks lowest at low samples sizes and increases AUC score in higher samples sizes with similar scores or sometimes better than to DESeq2. VoomLimma has a slightly better AUC score than EdgeR at lower sample sizes and this advantage is diminished at higher sample size.

