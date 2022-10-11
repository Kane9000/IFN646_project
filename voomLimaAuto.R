################################################################################
###Gene-level differential expression analysis using VoomLimma
################################################################################
#Input: synthetic data stored in folder "data" in current working directory
#Output: xlsx file for each synthetic data in current working directory

########################################
## Setup - Loading libraries
########################################
library(writexl)
library(limma)
library(edgeR)


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
  #################################
  ###Input data and Preprocessing
  #################################~
  data <- read.table(input_data, header=T, row.names=1)
  D <- input_group
  
  ####################################
  ###Create DGEList object
  ####################################
  d0 <- DGEList(data)
  
  ####################################
  ####Normalizing the data
  ####################################
  d0 <- calcNormFactors(d0)
  
  #################################
  ###Voom transformation and calculation of variance weights
  #################################
  #Specify model. Voom uses variances of the model residuals (observed - fitted)
  mm <- model.matrix(~0 + D)
  #Run Voom
  y <- voom(d0, mm)
  
  #################################
  ###Fitting linear models in limma
  #################################
  #lmFit fits a linear model using weighted least squares for each gene
  fit <- lmFit(y, mm)
  fit2 <- eBayes(fit)

  ####################################
  ###Output xlsx
  ####################################
  results <- topTable(fit2, number=10000, coef=1, sort.by="none")
  write_xlsx(results, output_file)
}

for (i in 1:9){
  func(input_data_list[[i]],input_group_list[[i]],output_file_list[[i]])
}