## Gene-level differential expression analysis using DESeq2

## Setup - Loading libraries
########################################
#For excel output
#install.packages('writexl')
library(writexl)

### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)

#output_file = '3_500_500.xlsx'
output_file = '3_750_250.xlsx'
#output_file = '3_1000_0.xlsx'


data <- read.table("./data/6_750_250.tsv", header=T, row.names=1)
meta <- read.table("6_meta.txt",header=T, row.names=1)

view(meta)
view(data)

#Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

#Generate normalization counts
########################################

#generates size fator
dds <- estimateSizeFactors(dds)
#view normalization factor applied to each sample
sizeFactors(dds)
#normalised count matrix
normalized_counts <- counts(dds, normalized=TRUE)


#DE analysis: exploring the dispersion estimates and assessing model fit
########################################

## Run analysis
dds <- DESeq(dds)

#define contrast
#results of condition2 compared to condition1
contrast_oe <- c("sampletype","condition2","condition1")

#the shrunkening
results_table_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
results_table <- lfcShrink(dds, contrast=contrast_oe, res=results_table_unshrunken, type="ashr")

results_table <- data.frame(results_table)
view(results_table)

write_xlsx(results_table, output_file)


