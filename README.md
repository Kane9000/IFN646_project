# IFN646_project
## Introduction

The introduction of high-throughput RNA sequencing technologies has allowed researchers to understand the complex interactions of cells to the environment. 
These experiments would involve the quantification of a large number of short reads which are then aligned to a genome to create a count matrix of each gene or feature. 
RNA differential expression tools apply computational models and allow the user to quickly find genes that are highly expressed or suppressed between two conditions.

This project aims to compare tools that analyze differential gene expression using RNA sequencing count data. 
Synthetic count matrices were generated and provided for the comparison which included two conditional groups, condition 1 and condition 2, with varying sample sizes from 3 to 9 per group and approximately ten thousand genes. 
A separate meta file was also provided that indicated whether condition 2 genes were either upregulated, downregulated or differentially expressed compared to condition 1.

## Method

The differential expression tools we will investigate in this study include DESeq2, EdgeR and Voom limma. 
The tools will be evaluated based on their ability to determine whether genes were differential expressed or not as indicated in the meta file. 
To compare the ability to rank differentially expressed genes, Precision-Recall F-measure scores were generated to calculate the accuracy of the tools.

A logistic regression was trained using the output of the differential expression tool as the response variable and a meta data file as the predictor variable. 
Genes with a statistically significant fold change of p value=0.05 or smaller were considered differentially expressed.
Additionally, fold change was required to match either upregulated or downregulated status as provided in the associated meta data file to be considered differentially expressed.
Source data, python and R scripts have been uploaded to this github.

RNA differential expression tools are built in the R language (version 4.2.1) and scripts were written in RStudio (version 2022.07.1 build 554). 
Results from these tools were analysed in Jupyter notebook (version 5.7.8) running Python version 3.7.3). 
Packages imported include DESeq2_1.36.0, edgeR_3.38.4, limma_3.52.4,  writexl_1.4.0.

## Results


###Precision-Recall F-measure score

![alt text](/PR_500_500.png?raw=true)
![alt text](/PR_750_250.png?raw=true)
![alt text](/PR_1000_0.png?raw=true)

Precision-Recall F-measure scores generated from voomLimma, DESeq2 and EdgeR using synthetic count data.
Each plot indicates the number of genes that are upregulated and downregulated in condition 2 compared to condition 1. 
Sample size varies from n = 3, 6 or 9 per group.
Each point represents a F-measure score that has been generated from either DESeq2 (circle), EdgeR (square) or VoomLimma (diamond).
 
The F-measure score was increased in all scenarios when the sample size was increased. 
DESeq2 has a slight advantage over EdgeR in most scenarios. 
EdgeR ranked lowest at low samples sizes and increases F-measure score in higher samples sizes with similar scores or sometimes better than to DESeq2. 
VoomLimma has a slightly better F-measure score than EdgeR at sample size 3 and this advantage was removed a sample sizes 6 or 9.

## Discussion

In this study, the precision-recall F-measure score is used to investigate the ability of a tool to correctly predict whether a gene was differentially expressed.
A precision-recall curve is a plot of the recall on the x axis and precision on the y axis.
Precision is referred to as the positive predictive value which describes how well a model is predicting a positive class. 
Precision is calculated by the true positive divided by the sum of the true positives and false positive (true positive / (true positive + false negative).
Recall, which is known as sensitivity, is calculated by the true positive divided by the sum of the true positive and false negative (true positive / (true positive + false negative)).
The F-measure score is the harmonic mean of the precision and recall metrics and is used for calculating the accuracy of binary classification model from the precision-recall curve.

Studies have shown that precision-recall curves are more appropriate when there is an imbalance between each class in a dataset which applies to our dataset as there are approximately nine times more genes that are not differentially expressed.
Area under Receiver Operating Characteristic AUC(ROC) is another common diagnostic tool that is used to calculate the performance of a binary classification model. 
The difference between a precision-recall curve and the AUC(ROC) is that the calculations do not include the use of true negatives which comprises most of our dataset and instead analyses a model’s ability to predict the minority class.
The precision-recall plot can calculate the F-measure or F1 score which is the measure of a model’s accuracy where 1 is the highest score and 0 is the lowest.

## References
Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. PLoS One. 2015 Mar 4;10(3):e0118432. doi: 10.1371/journal.pone.0118432. PMID: 25738806; PMCID: PMC4349800.
