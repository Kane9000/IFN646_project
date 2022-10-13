# IFN646_project
## Introduction

The introduction of high-throughput RNA sequencing technologies has allowed researchers to understand the complex interactions between cells and the environment. Gene expression analysis reveals the capacity of cells to react to external stimuli and allows for an in depth mechanistic studies that can piece together the internal signaling cascades and interactions between molecules. These experiments would involve the quantification of a large number of short reads which are then aligned to a genome to create a count matrix of each gene or feature. RNA differential expression tools apply computational models and allow the user to quickly find genes that are highly expressed or suppressed between conditions. Gene set enrichment analysis on the differential significantly expressed genes can be then used to identify functional gene set categories.

This project aims to compare tools that analyze differential gene expression using RNA sequencing count data. 
Synthetic count matrices were generated and provided for the comparison which included two conditional groups, condition 1 and condition 2, with varying sample sizes from 3 to 9 per group and approximately ten thousand genes. 
A separate meta file was also provided that indicated whether condition 2 genes were either upregulated, downregulated or differentially expressed compared to condition 1.

Differential expression tools have four basic functions, normalisation, dispersion estimation, differential expression followed by statistical analysis.

Normalisation

Normalisation is a critical step in differential expression analysis that is required to control the large differences in sequencing depths between samples so that accurate gene expression level comparisons can be done. 
DESeq computes a scaling factor by computing the median of the ratio.
EdgeR uses trimmed means of M values which computes a scaling factor by using a weighted average of a subset of genes that do not exhibit high average read counts or genes that have a large difference in expression.
The limma package uses quantile normalization which ensures that counts across all samples have the same distribution. Voom is a new normalization function added to the limma package which performs a LOWESS regression to estimate the mean-variance relation.
 
Statistical modeling of gene expression

Both DESeq and EdgeR use negative binomial distribution to estimate the relation between the variance and the mean in a poisson distribution model of gene expression across biological replicates.
However, edgeR and DESeq both have different methods of estimating the dispersion factor.
edgeR estimates the dispersion factor as a weighted combination of a gene-specific dispersion effect and a common dispersion effect calculated from all genes.
Alternatively, DESeq uses the mean expression of the gene and a second term that models the biological expression variability.
 
Differential expression test

EdgeR and DESeq use a variation of the fisher exact test adopted for negative binomial distribution to calculate significant changes in gene expression.
Limma implements a t-statistic with a modified standard error and degrees of freedom to calculate the P value.


## Method

The differential expression tools we will investigate in this study include DESeq2, EdgeR and Voom limma. The tools will be evaluated based on their ability to determine whether genes were differential expressed or not as indicated in the meta file. To compare the ability to rank differentially expressed genes, Precision-Recall F-measure scores were generated to calculate the accuracy of the tools.

A logistic regression was trained using the output of the differential expression tool as the response variable and a meta data file as the predictor variable. Genes with a statistically significant fold change with an adjusted p value=0.05 or smaller were considered differentially expressed. Additionally, fold change was required to match either upregulated or downregulated status as provided in the associated meta file to be considered differentially expressed. Source data, python and R scripts have been uploaded to this github (https://github.com/Kane9000/IFN646_project).

RNA differential expression tools are built in the R language (version 4.2.1) and scripts were written in RStudio (version 2022.07.1 build 554). Results from these tools were analyzed in a Jupyter notebook (version 5.7.8) running Python version 3.7.3). Packages imported include DESeq2_1.36.0, edgeR_3.38.4, limma_3.52.4, writexl_1.4.0.


## Results

Precision-Recall F-measure score

![alt text](/PR_500_500_padj.png?raw=true)
![alt text](/PR_750_250_padj.png?raw=true)
![alt text](/PR_1000_0_padj.png?raw=true)

Figure 1. Precision-Recall F-measure scores generated from voomLimma, DESeq2 and EdgeR using synthetic count data. Each plot indicates the number of genes that are upregulated and downregulated in condition 2 compared to condition 1. Sample size varies from n = 3, 6 or 9. Each point represents a F-measure score that has been generated from adjusted p values from 
either DESeq2 (circle), EdgeR (square) or VoomLimma (diamond).


The F-measure score in DESeq2 and edgeR was increased in all scenarios when the sample size was increased. VoomLimma had the best f-measure score in all scenarios except for 1 and performed better when all differentially expressed genes were upregulated as opposed to when a count matrix had a mixture of upregulated and downregulated genes. DESeq2 had an advantage over EdgeR in all scenarios, especially at the lowest sample size. Sample size did not appear to affect VoomLimma. At sample size 6 and 9 the difference between DESeq2 and EdgeR were minimal. 


Computational time

![alt text](/Time.png?raw=true)

Figure 2. Comparison of computational times of RNA differential expression tool DESeq2 (black bar), EdgeR (grey bar) and limma+voom (white bar) for all nine synthetic count matrices. 

The computational time for the differential expression tools to analyze all nine datasets were compared. Limma+voom had the quickest computational time at approximately 7 seconds, followed by EdgeR at 26 seconds and then DESeq2 at 44 seconds.


## Discussion

In this study, we compare RNA differential expression tools DESeq2, EdgeR and limma+voom using synthetic data. Log fold change and adjust p values were used to calculate a precision-recall curve to yield F-measure scores which were used to measure the performance of the algorithms. Limma+voom was found to have the highest ranking F-measure score in most of the provided data sets followed by DESeq2 and EdgeR. 

Our results show that there is an increase in F measure score as sample size increases. A similar study also using an imbalanced data set found a similar trend but compared differential expression tools with the Area Under Receiver Operating Characteristic Curve AUC(ROC) (Soneson & Delorenzi, 2013). AUC(ROC) is another common diagnostic tool that is used to calculate the performance of a binary classification model (Hanley JA & McNeil BJ 1982). However, studies have shown that precision-recall curves are more appropriate when there is an imbalance between each class in a dataset which applies to our dataset as there are approximately nine times more genes that are not differentially expressed (Saito T et al, 2015).  The difference between a precision-recall curve and the AUC(ROC) is that the calculations do not include the use of true negatives which comprises most of our dataset and instead analyze a model’s ability to predict the minority class.

This study found that there were notable differences in the time required for the differential expression tools to analyze all nine data sets. Limma-voom had the fastest computational time when analyzing the nine datasets followed by DESeq2 and EdgeR. A previous study by Law et al also found that limma-voom was one of the fastest statistical methods (Law CW, 2014).  


## References
Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. PLoS One. 2015 Mar 4;10(3):e0118432. doi: 10.1371/journal.pone.0118432. PMID: 25738806; PMCID: PMC4349800.
Law, C.W., Chen, Y., Shi, W. et al. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol 15, R29 (2014). https://doi.org/10.1186/gb-2014-15-2-r29
Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
Mark D. Robinson, Davis J. McCarthy, Gordon K. Smyth, edgeR: a Bioconductor package for differential expression analysis of digital gene expression data, Bioinformatics, Volume 26, Issue 1, 1 January 2010, Pages 139–140, https://doi.org/10.1093/bioinformatics/btp616
Soneson, C., Delorenzi, M. A comparison of methods for differential expression analysis of RNA-seq data. BMC Bioinformatics 14, 91 (2013). https://doi.org/10.1186/1471-2105-14-91
Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology. 1982;143: 29–36. pmid:7063747

