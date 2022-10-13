################################################################################
## Gene-level differential expression analysis using EdgeR
################################################################################
#Input: synthetic data stored in folder "data" in current working directory
#Output: xlsx file for each synthetic data in current working directory

########################################
## Setup - Loading libraries
########################################
library(edgeR)
library(writexl)
library(tictoc)

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
input_group_list <- list(
  c("condition1", "condition1", "condition1", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2"),
  c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2")
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

func <- function(input_data, input_group, output_file){
  ####################################
  ####Import data and Preprocessing
  ####################################
  data <- read.table(input_data, header=T, row.names=1)
  dataGroups = input_group
  
  ####################################
  ###create DGEList
  ####################################
  d <- DGEList(counts=data,group=factor(dataGroups))
  
  ####################################
  ####Normalizing the data
  ####################################
  d <- calcNormFactors(d)
  
  ####################################
  ###Set up Design matrix
  ####################################
  design.mat <- model.matrix(~ 0 + d$samples$group)
  colnames(design.mat) <- levels(d$samples$group)
  
  ####################################
  ###Estimates of dispersion
  ####################################
  #estimate common and tagwise dispersion
  d1 <- estimateDisp(d, design.mat)
  
  ####################################
  ###Differential Expression
  ####################################
  # compare groups 1 and 2
  et12 <- exactTest(d1, pair=c(1,2))
  #add adjusted p value
  results <- topTags(et12, n = 'inf', adjust.method = "BH", sort.by = "none", p.value = 1)

  ####################################
  ###Output xlsx
  ####################################
  results <- data.frame(results)
  write_xlsx(results, output_file)
}
tic()
for (i in 1:9){
  func(input_data_list[[i]],input_group_list[[i]],output_file_list[[i]])
}
toc()



