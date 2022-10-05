library(edgeR)
library(writexl)
####################################
####Import data
####################################
#output_file = '3_500_500.xlsx'
#output_file = '3_750_250.xlsx'
#output_file = '3_1000_0.xlsx'
#output_file = '6_500_500.xlsx'
#output_file = '6_750_250.xlsx'
#output_file = '6_1000_0.xlsx'
#output_file = '9_500_500.xlsx'
#output_file = '9_750_250.xlsx'
output_file = '9_1000_0.xlsx'

data <- read.table("./data/9_1000_0.tsv", header=T, row.names=1)

#assigning groups
#dataGroups <- c("condition1", "condition1", "condition1", "condition2", "condition2", "condition2")
#dataGroups <- c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2")
dataGroups <- c("condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2", "condition2")

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


