################################################################################
###Gene-level differential expression analysis using DESeq2
################################################################################
#Input: synthetic data stored in folder "data" in current working directory
#Output: xlsx file for each synthetic data in current working directory

########################################
## Setup - Loading libraries
########################################
#For excel output
#install.packages('writexl')
#install.packages('ashr')
library(writexl)
library(DESeq2)


output_file_list <- list(
  ('3_500_500.xlsx'),
  ('3_750_250.xlsx'),
  ('3_1000_0.xlsx'),
  ('6_500_500.xlsx'),
  ('6_750_250.xlsx'),
  ('6_1000_0.xlsx'),
  ('9_500_500.xlsx'),
  ('9_750_250.xlsx'),
  ('9_1000_0.xlsx')
)

input_data_list <- list(
  "./data/3_500_500.tsv",
  "./data/3_750_250.tsv",
  "./data/3_1000_0.tsv",
  "./data/6_500_500.tsv",
  "./data/6_750_250.tsv",
  "./data/6_1000_0.tsv",
  "./data/9_500_500.tsv",
  "./data/9_750_250.tsv",
  "./data/9_1000_0.tsv"
)
input_group_list <- list(
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"))),
  (data.frame(sampletype =c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2")))
)

input_meta_list <- list(
  ("3_meta.txt"),
  ("3_meta.txt"),
  ("3_meta.txt"),
  ("6_meta.txt"),
  ("6_meta.txt"),
  ("6_meta.txt"),
  ("9_meta.txt"),
  ("9_meta.txt"),
  ("9_meta.txt")
)

func <- function(input_data, input_group, output_file){
  ####################################
  ####Import data
  ####################################
  data <- read.table(input_data, header=T, row.names=1)
  meta <- input_group
  
  #Create DESeq2Dataset object
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
  
  ####################################
  ####Normalizing the data
  ####################################
  
  #generates size fator
  dds <- estimateSizeFactors(dds)
  #view normalization factor applied to each sample
  sizeFactors(dds)
  #normalised count matrix
  normalized_counts <- counts(dds, normalized=TRUE)
  
  ########################################
  #DE analysis: exploring the dispersion estimates and assessing model fit
  ########################################
  
  ## Run analysis
  dds <- DESeq(dds)
  
  #define contrast
  #results of condition2 compared to condition1
  contrast_oe <- c("sampletype","condition2","condition1")
  
  #Shrink data
  results_table_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
  results_table <- lfcShrink(dds, contrast=contrast_oe, res=results_table_unshrunken, type="ashr")
  
  #View and write data
  results_table <- data.frame(results_table)
  #view(results_table)
  write_xlsx(results_table, output_file)
  
}

for (i in 1:9){
  func(input_data_list[[i]],input_group_list[[i]],output_file_list[[i]])
}

