rm(list = ls())
#Affymetrix Human Genome U133 Plus 2.0 Array
library(GEOquery)
library(dplyr)
library(affy)
library(stringi)
library(stringr)
library(R.utils)
library(limma)


geo_GSE46449 = getGEO("GSE46449")
geo_GSE46449 = geo_GSE46449[[1]]
pheno_GSE46449 = geo_GSE46449@phenoData@data

#the dataset contains replicates, so let's filter them out (include into analysis only biological replicate 1)
BD_data = pheno_GSE46449 %>%
  filter(grepl("1$", title, ignore.case = T))
BD_data_ids = rownames(BD_data)

#clean colnames
colnames(BD_data)[35] = "diagnosis"
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = "age:ch1",
                                           replacement = "age")
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = "gender:ch1",
                                           replacement = "gender")

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control subject", "bipolar patient"))
BD_data$age = as.numeric(BD_data$age)
BD_data$gender = factor(BD_data$gender)


#get the raw files
files = list.files("GSE46449_RAW/")
paths = paste0("GSE46449_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE46449_RAW")
cel_paths = paste0("GSE46449_RAW/", files)

#filter only BD and control files
files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE46449_RAW/", files_filtered)

#read the rawdata
rawdata_GSE46449 = affy::ReadAffy(filenames = cel_paths_filtered)

#quality check
affy::boxplot(rawdata_GSE46449)
affy::hist(rawdata_GSE46449)

degradation = affy::AffyRNAdeg(rawdata_GSE46449)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL)

eset_GSE46449 = affy::expresso(afbatch = rawdata_GSE46449, 
                               bgcorrect.method = "rma",
                               pmcorrect.method = "pmonly",
                               summary.method = "medianpolish")

eset_GSE46449_expression = affy::exprs(eset_GSE46449)
colnames(eset_GSE46449_expression) = sapply(colnames(eset_GSE46449_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE46449_expression) == rownames(BD_data$geo_accession)) #TRUE

#DGE WITH covars
#create design matrix
design = model.matrix(~diagnosis+age, data = BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE46449_expression, design)
fit = eBayes(fit)
GSE46449_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE46449_topTable_WITH_covars$ID = rownames(GSE46449_topTable_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,4)]
colnames(SE) = c("SE", "ID")
GSE46449_topTable_WITH_covars = inner_join(GSE46449_topTable_WITH_covars, SE, by = "ID")

#annotate the table
annotation_data = geo_GSE46449@featureData@data
GSE46449_topTable_WITH_covars = inner_join(GSE46449_topTable_WITH_covars, annotation_data[, c("Gene Symbol", "ID")], by = "ID")

GSE46449_topTable_WITH_covars$`Gene Symbol`[GSE46449_topTable_WITH_covars$`Gene Symbol`==""] = NA
GSE46449_topTable_WITH_covars = GSE46449_topTable_WITH_covars[!is.na(GSE46449_topTable_WITH_covars$`Gene Symbol`),]

GSE46449_topTable_WITH_covars = GSE46449_topTable_WITH_covars[, c(9,11,1,10, 2:8)]

#DGE NO covars
#create design matrix
design0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(design0)==rownames(BD_data)) #TRUE

#fit the model
fit0 = lmFit(eset_GSE46449_expression, design0)
fit0 = eBayes(fit0)
GSE46449_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE46449_topTable_NO_covars$ID = rownames(GSE46449_topTable_NO_covars)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, c(2,3)]
colnames(SE0) = c("SE", "ID")
GSE46449_topTable_NO_covars = inner_join(GSE46449_topTable_NO_covars, SE0, by = "ID")

#annotate the table
GSE46449_topTable_NO_covars = inner_join(GSE46449_topTable_NO_covars, annotation_data[, c("Gene Symbol", "ID")], by = "ID")

GSE46449_topTable_NO_covars$`Gene Symbol`[GSE46449_topTable_NO_covars$`Gene Symbol`==""] = NA
GSE46449_topTable_NO_covars = GSE46449_topTable_NO_covars[!is.na(GSE46449_topTable_NO_covars$`Gene Symbol`),]
GSE46449_topTable_NO_covars = GSE46449_topTable_NO_covars[, c(9,11,1,10,2:8)]

#save the results
write.csv(GSE46449_topTable_NO_covars, "GSE46449_topTable_NO_covars.csv")
write.csv(GSE46449_topTable_WITH_covars, "GSE46449_topTable_WITH_covars.csv")

#clean the tables
GSE46449_topTable_NO_covars = GSE46449_topTable_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE46449_topTable_NO_covars = GSE46449_topTable_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE46449_topTable_WITH_covars = GSE46449_topTable_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE46449_topTable_WITH_covars = GSE46449_topTable_WITH_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

#add study column
GSE46449_topTable_NO_covars$Study = "GSE46449"
GSE46449_topTable_WITH_covars$Study = "GSE46449"

#save
write_csv(GSE46449_topTable_NO_covars, "GSE46449_NO_covars.csv")
write_csv(GSE46449_topTable_WITH_covars, "GSE46449_WITH_covars.csv")







