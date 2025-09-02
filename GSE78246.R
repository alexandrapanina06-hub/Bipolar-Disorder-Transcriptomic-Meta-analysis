###GSE78246 Affymetrix Human Gene 1.0 ST Array -> oligo or xps package

rm(list=ls())
library(R.utils)
library(dplyr)
library(limma)
library(GEOquery)
library(affycoretools)
library(affy)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligo")
library(oligo)
BiocManager::install("affycoretools")
library(AnnotationDbi)

BiocManager::install("huex10sttranscriptcluster.db")
library(huex10sttranscriptcluster.db)



#get the file
geo_GSE78246 = GEOquery::getGEO("GSE78246")
geo_GSE78246 = geo_GSE78246[[1]]
geo_GSE78246_pheno = geo_GSE78246@phenoData@data

#get the probes
GSE78246_probes = geo_GSE78246@featureData@data

#get raw files 
files = list.files("GSE78246_RAW/")
paths = paste0("GSE78246_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE78246_RAW")
cel_paths = paste0("GSE78246_RAW/", files)

#get raw data using oligo
rawdata_GSE78246 = read.celfiles(filenames = cel_paths)

#quality check
boxplot(rawdata_GSE78246, target = "core", main = "Raw Intensity Distribution") 
#boxplot for oligo objects requires specification of target

#data normalization
eset_GSE78246 = rma(rawdata_GSE78246)
eset_GSE78246_expression = exprs(eset_GSE78246)

#clean colnames in eset_GSE78246
colnames(eset_GSE78246_expression) = sapply(colnames(eset_GSE78246_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE78246_expression) == rownames(geo_GSE78246_pheno$geo_accession)) #TRUE

#create annotation file using huex10sttranscriptcluster.db
probe_ids = rownames(eset_GSE78246_expression)
annotation_data = select(huex10sttranscriptcluster.db, 
                          keys = probe_ids, 
                          columns = c("SYMBOL", "ENTREZID", "GENENAME"), 
                          keytype = "PROBEID")


#turn diagnosis and sex into factors
geo_GSE78246_pheno$`diagnosis:ch1` = factor(geo_GSE78246_pheno$`diagnosis:ch1`, levels = c("control", "bipolar disorder"))
geo_GSE78246_pheno$`Sex:ch1` = factor(geo_GSE78246_pheno$`Sex:ch1`)

##DE WITH COVARIATES (sex)
#create design matrix
design = model.matrix(~`diagnosis:ch1` + `Sex:ch1`, data = geo_GSE78246_pheno)
all(rownames(design) == rownames(geo_GSE78246_pheno)) #TRUE

#fit the model
fit = lmFit(eset_GSE78246_expression, design = design)
fit = eBayes(fit = fit)
GSE78246_topTable_WITH_covariates = topTable(fit = fit, coef = 2, adjust.method = 'fdr', number = Inf, confint = T)
GSE78246_topTable_WITH_covariates$ID = rownames(GSE78246_topTable_WITH_covariates)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post) * fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 4)]
colnames(SE) = c("SE", "ID")
GSE78246_topTable_WITH_covariates = inner_join(GSE78246_topTable_WITH_covariates, SE, by="ID")

#annotate
GSE78246_topTable_WITH_covariates = inner_join(GSE78246_topTable_WITH_covariates, 
                                               annotation_data[, c("SYMBOL", "PROBEID")], by = c("ID" = "PROBEID"))
GSE78246_topTable_WITH_covariates = GSE78246_topTable_WITH_covariates[, c(9, 11, 1:8, 10)]
GSE78246_topTable_WITH_covariates_filtered = GSE78246_topTable_WITH_covariates %>%
  filter(!is.na(SYMBOL))

#DE WITHOUT COVARIATES
#create new design matrix
design_no_covars = model.matrix(~`diagnosis:ch1`, data = geo_GSE78246_pheno)
all(rownames(design_no_covars) == rownames(geo_GSE78246_pheno)) #TRUE

fit1 = lmFit(eset_GSE78246_expression, design_no_covars)
fit1 = eBayes(fit1)
GSE78246_topTable_NO_covariates = topTable(fit = fit1, coef = 2, adjust.method = 'fdr', number = Inf, confint = T)
GSE78246_topTable_NO_covariates$ID = rownames(GSE78246_topTable_NO_covariates)

#add SE
SE1 = as.data.frame(sqrt(fit1$s2.post) * fit1$stdev.unscaled)
SE1$ID = rownames(SE1)
SE1 = SE1[, 2:3]
colnames(SE1) = c("SE", "ID")
GSE78246_topTable_NO_covariates = inner_join(GSE78246_topTable_NO_covariates, SE1, by="ID")

#annotate
GSE78246_topTable_NO_covariates = inner_join(GSE78246_topTable_NO_covariates,
                                             annotation_data[, c("SYMBOL", "PROBEID")], by = c("ID" = "PROBEID"))

GSE78246_topTable_NO_covariates = GSE78246_topTable_NO_covariates[, c(9, 11, 1:8, 10)]
GSE78246_topTable_NO_covariates_filtered = GSE78246_topTable_NO_covariates %>%
  filter(!is.na(SYMBOL))

#save results
write.csv(GSE78246_topTable_NO_covariates_filtered, "GSE78246_table_NO_covars.csv")
write.csv(GSE78246_topTable_WITH_covariates_filtered, "GSE78246_table_WITH_covars.csv")

#clean the tables
colnames(GSE78246_table_NO_covars)[3] = "Gene Symbol"
colnames(GSE78246_table_WITH_covars)[3] = "Gene Symbol"

#no covars 
#delete rows with more than one gene mapped
GSE78246_table_NO_covars_cleaned = GSE78246_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE78246_table_NO_covars_cleaned = GSE78246_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE78246_table_NO_covars_cleaned$Study = "GSE78246"
write_csv(GSE78246_table_NO_covars_cleaned, "GSE78246_no_covars_cleaned.csv")


#with covars 
#delete rows with more than one gene mapped
GSE78246_table_WITH_covars_cleaned = GSE78246_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE78246_table_WITH_covars_cleaned = GSE78246_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE78246_table_WITH_covars_cleaned$Study = "GSE78246"
write_csv(GSE78246_table_WITH_covars_cleaned, "GSE78246_with_covars_cleaned.csv")





