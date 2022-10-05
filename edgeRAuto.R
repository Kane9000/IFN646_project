library(edgeR)
library(writexl)

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
  ####Import data
  ####################################
  data <- read.table(input_data, header=T, row.names=1)
  dataGroups = input_group
  
  #create DGEList
  d <- DGEList(counts=data,group=factor(dataGroups))
  
  ####################################
  ####Normalizing the data
  ####################################
  d <- calcNormFactors(d, method="TMM")
  
  
  ####################################
  ####Set up model
  ####################################
  #Design matrix
  design.mat <- model.matrix(~ 0 + d$samples$group)
  colnames(design.mat) <- levels(d$samples$group)
  
  design.mat
  ####################################
  ###GLM estimates of dispersion
  ####################################
  
  #estimate common and tagwise dispersion
  #d1 <- estimateDisp(d)
  
  d2 <- estimateGLMCommonDisp(d,design.mat)
  d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
  d2 <- estimateGLMTagwiseDisp(d2,design.mat)
  # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
  # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
  
  
  ####################################
  ###Differential Expression
  ####################################
  # tagwise tests using the exact negative binomial test
  # compare groups 1 and 2
  et12 <- exactTest(d2, pair=c(1,2))
  
  #most significant genes
  #topTags(et12)
  
  et12 <- data.frame(et12)
  write_xlsx(et12, output_file)
}
for (i in 1:9){
  func(input_data_list[[i]],input_group_list[[i]],output_file_list[[i]])
}



