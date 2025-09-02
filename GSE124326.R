rm(list=ls())
#	Illumina HiSeq 2500
library(GEOquery)
library(dplyr)
library(limma)
library(stringi)
library(stringr)
library(edgeR)
library(dplyr)

geo_GSE124326 = getGEO("GSE124326")
geo_GSE124326 = geo_GSE124326[[1]]
pheno_GSE124326 = geo_GSE124326@phenoData@data

#create diagnosis column
pheno_GSE124326 = pheno_GSE124326 %>%
  mutate(diagnosis = case_when(
    grepl("bipolar disorder diagnosis: BP1", characteristics_ch1) ~ "bipolar disorder",
    grepl("bipolar disorder diagnosis: BP2", characteristics_ch1) ~ "bipolar disorder",
    grepl("bipolar disorder diagnosis: Control", characteristics_ch1) ~ "control", 
    TRUE ~ "other"
  ))

#clean the column names
colnames(pheno_GSE124326) = stri_replace_all_fixed(str = colnames(pheno_GSE124326),
                                                      pattern = 'characteristics_ch1.1',
                                                      replacement = 'age')
colnames(pheno_GSE124326) = stri_replace_all_fixed(str = colnames(pheno_GSE124326),
                                                   pattern = 'characteristics_ch1.2',
                                                   replacement = 'gender')
colnames(pheno_GSE124326) = stri_replace_all_fixed(str = colnames(pheno_GSE124326),
                                                   pattern = 'characteristics_ch1.3',
                                                   replacement = 'lithium')
colnames(pheno_GSE124326) = stri_replace_all_fixed(str = colnames(pheno_GSE124326),
                                                   pattern = 'characteristics_ch1.4',
                                                   replacement = 'tobacco')

#clean the rows
pheno_GSE124326$age = gsub("age: ", "", pheno_GSE124326$age)
pheno_GSE124326$gender = gsub("Sex: ", "", pheno_GSE124326$gender)
pheno_GSE124326$lithium = gsub("lithium use (non-user=0, user = 1): ", "", pheno_GSE124326$`lithium use (non-user=0, user = 1):ch1`)
pheno_GSE124326$tobacco = gsub("tobacco use: ", "", pheno_GSE124326$tobacco)

#variable transformation
pheno_GSE124326$diagnosis = factor(pheno_GSE124326$diagnosis, levels = c("control", "bipolar disorder"))
pheno_GSE124326$age = as.numeric(pheno_GSE124326$age)
pheno_GSE124326$gender = factor(pheno_GSE124326$gender)
pheno_GSE124326$lithium = factor(pheno_GSE124326$lithium)
pheno_GSE124326$tobacco = factor(pheno_GSE124326$tobacco)

#create counts matrix
count_matrix = read.table(file = "GSE124326_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = '\t')
rownames(count_matrix) = count_matrix$GeneID

#get the annotation file
columns = c('GeneID', 'Symbol', 'chromosome', 'GeneType')
annotation_data = read.columns(file = 'Human.GRCh38.p13.annot.tsv', 
                               required.col = columns, 
                               sep = '\t', 
                               stringsAsFactors = F)
annotation_data$GeneID = factor(annotation_data$GeneID)
colnames(annotation_data)[2] = "Gene Symbol"

#fit the model
#DGE NO covars
#create design matrix
design = model.matrix(~diagnosis, data = pheno_GSE124326)
all(rownames(design) == rownames(pheno_GSE124326)) #TRUE

#create dge list
dge = DGEList(counts = count_matrix, annotation.columns = "GeneID")
keep = filterByExpr(dge, design) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] 

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
design = design[rownames(design) %in% rownames(dge$samples),]
all(rownames(dge$samples) == rownames(design)) #TRUE

#convert data to log2 scale and estimate mean-variance trend
v = voom(dge, design, plot = T)
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design)
fit = eBayes(fit)
GSE124326_topTable_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 3)]
colnames(SE) = c("SE", "GeneID")
GSE124326_topTable_NO_covars$GeneID = factor(GSE124326_topTable_NO_covars$GeneID)
GSE124326_topTable_NO_covars = inner_join(GSE124326_topTable_NO_covars, SE, by = "GeneID")

#annotate the table
GSE124326_topTable_NO_covars = inner_join(GSE124326_topTable_NO_covars, annotation_data[, c("Gene Symbol", "GeneID")], by = "GeneID")
GSE124326_topTable_NO_covars = GSE124326_topTable_NO_covars[, c(1, 11, 2, 10, 3:9)]
any(GSE124326_topTable_NO_covars$`Gene Symbol` == "") #FALSE
any(is.na(GSE124326_topTable_NO_covars$`Gene Symbol`)) #FALSE




#DGE WITH covars
#remove rows where covariate annotation causes problems
pheno_GSE124326 = pheno_GSE124326[rownames(pheno_GSE124326) != "GSM3529491",]
pheno_GSE124326 = pheno_GSE124326[rownames(pheno_GSE124326) != "GSM3529591",]
design1 = model.matrix(~diagnosis + age + gender + lithium, data = pheno_GSE124326)

count_matrix = count_matrix[, colnames(count_matrix) != "GSM3529491"]
count_matrix = count_matrix[, colnames(count_matrix) != "GSM3529591"]

#create dge list
dge1 = DGEList(counts = count_matrix, annotation.columns = "GeneID")
keep1 = filterByExpr(dge1, design1) #remove very low counts
dge1 = dge1[keep1, , keep.lib.sizes = F] 

#scale normalization
dge1 = calcNormFactors(dge1) #performs TTM normalization
#since the dimensions do not match, let's subset design1 and dge1 
setdiff(rownames(design1), rownames(dge1$samples))
common_samples = intersect(rownames(dge1$samples), rownames(design1))
dge1 = dge1[, common_samples]
design1 = design1[common_samples, , drop = FALSE]
all(rownames(dge1$samples) == rownames(design1)) #TRUE
design1 = design1[, c(1:3, 5:6)]

#convert the data into log2 scale
v1 = voom(dge1, design1, plot = T)
all(rownames(v1$targets) == rownames(v1$design)) #TRUE

#fit the model
fit1 = lmFit(v1, design1)
fit1 = eBayes(fit1)
GSE124326_topTable_WITH_covars = topTable(fit = fit1, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE1 = as.data.frame(sqrt(fit1$s2.post)*fit1$stdev.unscaled)
SE1$GeneID = rownames(SE1)
SE1 = SE1[, c(2, 6)]
colnames(SE1) = c("SE", "GeneID")
GSE124326_topTable_WITH_covars$GeneID = factor(GSE124326_topTable_WITH_covars$GeneID)
GSE124326_topTable_WITH_covars = inner_join(GSE124326_topTable_WITH_covars, SE1, by = "GeneID")


#annotate the table
GSE124326_topTable_WITH_covars = inner_join(GSE124326_topTable_WITH_covars, annotation_data[, c("Gene Symbol", "GeneID")], by = "GeneID")
GSE124326_topTable_WITH_covars = GSE124326_topTable_WITH_covars[, c(1, 11, 2, 10, 3:9)]
any(GSE124326_topTable_WITH_covars$`Gene Symbol` == "") #FALSE
any(is.na(GSE124326_topTable_WITH_covars$`Gene Symbol`)) #FALSE

#save the results
write.csv(GSE124326_topTable_NO_covars, "GSE124326_topTable_NO_covars.csv")
write.csv(GSE124326_topTable_WITH_covars, "GSE124326_topTable_WITH_covars.csv")


#clean the top tables
#delete rows with more than one gene mapped
GSE124326_topTable_NO_covars = GSE124326_topTable_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE124326_topTable_NO_covars = GSE124326_topTable_NO_covars %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error
#modify logFC and SE
GSE124326_topTable_NO_covars = GSE124326_topTable_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

#add study column
GSE124326_topTable_NO_covars$Study = "GSE124326"

#save the cleaned table
write_csv(GSE124326_topTable_NO_covars, "GSE124326_NO_covars.csv")


#clean the top tables
#delete rows with more than one gene mapped
GSE124326_topTable_WITH_covars = GSE124326_topTable_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE124326_topTable_WITH_covars = GSE124326_topTable_WITH_covars %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error
#modify logFC and SE
GSE124326_topTable_WITH_covars = GSE124326_topTable_WITH_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

#add study column
GSE124326_topTable_WITH_covars$Study = "GSE124326"

#save the cleaned table
write_csv(GSE124326_topTable_WITH_covars, "GSE124326_WITH_covars.csv")










