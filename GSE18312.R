rm(list = ls())
#Affymetrix Human Exon 1.0 ST Array

library(R.utils)
library(dplyr)
library(limma)
library(GEOquery)
library(affycoretools)
library(affy)
library(oligo)
library(huex10sttranscriptcluster.db)

geo_GSE18312 = getGEO("GSE18312")
geo_GSE18312 = geo_GSE18312[[1]]
pheno_GSE18312 = geo_GSE18312@phenoData@data


#clean colnames
colnames(pheno_GSE18312) = stri_replace_all_fixed(str = colnames(pheno_GSE18312),
                                                  pattern = "diagnosis:ch1",
                                                  replacement = "diagnosis")
colnames(pheno_GSE18312) = stri_replace_all_fixed(str = colnames(pheno_GSE18312),
                                                  pattern = "gender:ch1",
                                                  replacement = "gender")
colnames(pheno_GSE18312) = stri_replace_all_fixed(str = colnames(pheno_GSE18312),
                                                  pattern = "age:ch1",
                                                  replacement = "age")
colnames(pheno_GSE18312) = stri_replace_all_fixed(str = colnames(pheno_GSE18312),
                                                  pattern = "ancestry:ch1",
                                                  replacement = "ancestry")
colnames(pheno_GSE18312) = stri_replace_all_fixed(str = colnames(pheno_GSE18312),
                                                  pattern = "current tobacco use:ch1",
                                                  replacement = "tobacco")
colnames(pheno_GSE18312) = stri_replace_all_fixed(str = colnames(pheno_GSE18312),
                                                  pattern = "education:ch1",
                                                  replacement = "education")
colnames(pheno_GSE18312) = stri_replace_all_fixed(str = colnames(pheno_GSE18312),
                                                  pattern = "history of psychosis:ch1",
                                                  replacement = "psychosis")
pheno_GSE18312$diagnosis = sapply(pheno_GSE18312$diagnosis, tolower)

#filter only BD and control samples
BD_data = pheno_GSE18312 %>%
  filter(grepl("bipolar disorder", diagnosis, ignore.case = T)|
         grepl("control", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)


#get the raw files
files = list.files("GSE18312_RAW/")
paths = paste0("GSE18312_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE18312_RAW")
cel_paths = paste0("GSE18312_RAW/", files)

#filter only BD and control files
files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE18312_RAW/", files_filtered)

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, 
                           levels = c("control", "bipolar disorder"))
BD_data$age = as.numeric(BD_data$age)
BD_data$ancestry = factor(BD_data$ancestry)
BD_data$tobacco = factor(BD_data$tobacco)
BD_data$education = factor(BD_data$education)
BD_data$gender = factor(BD_data$gender)
BD_data$psychosis = factor(BD_data$psychosis)

#get the raw data using oligo package
rawdata_GSE18312 = read.celfiles(filenames = cel_paths_filtered)

#quality check
boxplot(rawdata_GSE18312, target = "core", main = "Raw Intensity Distribution")

#data normalization
eset_GSE18312 = rma(rawdata_GSE18312)
eset_GSE18312_expression = exprs(eset_GSE18312)

#clean colnames in eset_GSE78246
colnames(eset_GSE18312_expression) = sapply(colnames(eset_GSE18312_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE18312_expression) == rownames(BD_data$geo_accession)) #TRUE

#create annotation file using huex10sttranscriptcluster.db
probe_ids = rownames(eset_GSE18312_expression)
annotation_data = select(huex10sttranscriptcluster.db, 
                         keys = probe_ids, 
                         columns = c("SYMBOL", "ENTREZID", "GENENAME"), 
                         keytype = "PROBEID")
#DGE with covars
#create design matrix
design = model.matrix(~diagnosis+age+ancestry+tobacco+education+gender+psychosis, 
                      data = BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE18312_expression, design = design)
fit = eBayes(fit)
GSE18312_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE18312_topTable_WITH_covars$ID = rownames(GSE18312_topTable_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post) * fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 17)]
colnames(SE) = c("SE", "ID")
GSE18312_topTable_WITH_covars = inner_join(GSE18312_topTable_WITH_covars, SE, by="ID")

#annotate the table
GSE18312_topTable_WITH_covars = inner_join(GSE18312_topTable_WITH_covars, 
                                               annotation_data[, c("SYMBOL", "PROBEID")], by = c("ID" = "PROBEID"))
GSE18312_topTable_WITH_covars = GSE18312_topTable_WITH_covars[!is.na(GSE18312_topTable_WITH_covars$SYMBOL),]
colnames(GSE18312_topTable_WITH_covars)[11] = "Gene Symbol"
GSE18312_topTable_WITH_covars = GSE18312_topTable_WITH_covars[, c(9, 11, 1, 10, 2:8)]

#DGE NO covars
#create design matrix
design0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(design0) == rownames(BD_data)) #TRUE

#fit the model
fit0 = lmFit(eset_GSE18312_expression, design = design0)
fit0 = eBayes(fit0)
GSE18312_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE18312_topTable_NO_covars$ID = rownames(GSE18312_topTable_NO_covars)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post) * fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, c(2:3)]
colnames(SE0) = c("SE", "ID")
GSE18312_topTable_NO_covars = inner_join(GSE18312_topTable_NO_covars, SE0, by = "ID")

#annotate the table
GSE18312_topTable_NO_covars = inner_join(GSE18312_topTable_NO_covars,
                                         annotation_data[, c("SYMBOL", "PROBEID")], by = c("ID" = "PROBEID"))
GSE18312_topTable_NO_covars = GSE18312_topTable_NO_covars[!is.na(GSE18312_topTable_NO_covars$SYMBOL),]
colnames(GSE18312_topTable_NO_covars)[11] = "Gene Symbol"
GSE18312_topTable_NO_covars = GSE18312_topTable_NO_covars[, c(9,11, 1, 10, 2:8)]

#save the results
write.csv(GSE18312_topTable_NO_covars, "GSE18312_topTable_NO_covars.csv")
write.csv(GSE18312_topTable_WITH_covars, "GSE18312_topTable_WITH_covars.csv")

##clean the tables
GSE18312_topTable_NO_covars = GSE18312_topTable_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE18312_topTable_NO_covars = GSE18312_topTable_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 


GSE18312_topTable_WITH_covars = GSE18312_topTable_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE18312_topTable_WITH_covars = GSE18312_topTable_WITH_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
#add study column
GSE18312_topTable_NO_covars$Study = "GSE18312"
GSE18312_topTable_WITH_covars$Study = "GSE18312"
#save the results
write_csv(GSE18312_topTable_NO_covars, "GSE18312_NO_covars.csv")
write_csv(GSE18312_topTable_WITH_covars, "GSE18312_WITH_covars.csv")










