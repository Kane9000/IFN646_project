# IFN646_project
## Introduction

The introduction of high-throughput RNA sequencing technologies has allowed researchers to understand the complex interactions between cells and the environment. Gene expression analysis reveals the reaction and capacity of cells to external stimuli and allows for an in depth mechanism of action studies. These experiments would involve the quantification of a large number of short reads which are then aligned to a genome to create a count matrix of each gene or feature. RNA differential expression tools apply computational models and allow the user to quickly find genes that are highly expressed or suppressed between two conditions.

This project aims to compare tools that analyze differential gene expression using RNA sequencing count data. Synthetic count matrices were generated and provided for the comparison which included two conditional groups, condition 1 and condition 2, with varying sample sizes from 3 to 9 per group and approximately ten thousand genes. A separate meta file was also provided that indicated whether condition 2 genes were either upregulated, downregulated or differentially expressed compared to condition 1.

Discuss tools
-Normalisation
-estimate dispersion
-differential expression
-statistics

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

The F-measure score was increased in all scenarios when the sample size was increased. DESeq2 had an advantage over EdgeR in all scenarios. At sample size 6 and 9 the difference between DESeq2 and EdgeR were minimal. EdgeR second and VoomLimma ranked last and underperformed EdgeR, especially at the lowest sample size.


## Discussion

In this study, log fold change and adjust p values were used to calculate a precision-recall curve to yield F-measure scores which were used to measure the performance of RNA differential expression algorithms. DESeq2 was found to have the highest ranking F-measure score in all of the provided data sets followed by EdgeR and Voom-Limma.

A precision-recall curve is a plot of the recall on the x axis and precision on the y axis. Precision is referred to as the positive predictive value which describes how well a model is predicting a positive class. Precision is calculated by the true positive divided by the sum of the true positives and false positive (true positive / (true positive + false negative). Recall, which is known as sensitivity, is calculated by the true positive divided by the sum of the true positive and false negative (true positive / (true positive + false negative)). The F-measure score is the harmonic mean of the precision and recall metrics and is used for calculating the accuracy of a binary classification model from the precision-recall curve.

Studies have shown that precision-recall curves are more appropriate when there is an imbalance between each class in a dataset which applies to our dataset as there are approximately nine times more genes that are not differentially expressed. Area under Receiver Operating Characteristic AUC(ROC) is another common diagnostic tool that is used to calculate the performance of a binary classification model. The difference between a precision-recall curve and the AUC(ROC) is that the calculations do not include the use of true negatives which comprises most of our dataset and instead analyze a model’s ability to predict the minority class. The precision-recall plot can calculate the F-measure or F1 score which is the measure of a model’s accuracy where 1 is the highest score and 0 is the lowest.

Discuss p value adjusted…


## References
Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. PLoS One. 2015 Mar 4;10(3):e0118432. doi: 10.1371/journal.pone.0118432. PMID: 25738806; PMCID: PMC4349800.
