##GSE5392 Affymetrix U133 array
rm(list = ls())
library(R.utils)
library(dplyr)
library(affy)
library(stringi)
library(limma)
library(GEOquery)
geo_GSE5392 = GEOquery::getGEO("GSE5392")
geo_GSE5392 = geo_GSE5392[[1]]
geo_GSE5392_pheno = geo_GSE5392@phenoData@data

#get the raw data
files = list.files("GSE5392_RAW/")
paths = paste0("GSE5392_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE5392_RAW/")
cel_paths = paste0("GSE5392_RAW/", files)

#read the raw data
rawdata_GSE5392 = affy::ReadAffy(filenames = cel_paths)

#quality check
affy::boxplot(rawdata_GSE5392)
affy::hist(rawdata_GSE5392)

degradation = affy::AffyRNAdeg(rawdata_GSE5392)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL)

#normalization
eset_GSE5392 = affy::expresso(afbatch = rawdata_GSE5392, 
                                bgcorrect.method = "rma",
                                pmcorrect.method = "pmonly",
                                summary.method = "medianpolish")
eset_GSE5392_expression = affy::exprs(eset_GSE5392)

#clean eset colnames
colnames(eset_GSE5392_expression) = sapply(colnames(eset_GSE5392_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})
all(colnames(eset_GSE5392_expression) == rownames(geo_GSE5392_pheno$geo_accession)) #TRUE

#create diagnosis column
geo_GSE5392_pheno = geo_GSE5392_pheno %>%
  mutate(diagnosis = case_when(
    grepl("Healthy control", `Disease_status:ch1`) ~ "control",
    grepl("Bipolar disorder", `Disease_status:ch1`) ~ "bipolar disorder", 
    TRUE ~ "Other"
  ))

#clear colnames in phenodata
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                    pattern = 'Age (years):ch1',
                                                    replacement = 'age')

colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Age of onset (years):ch1',
                                                     replacement = 'age_onset')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Alcohol abuse (ratings scale:ch1',
                                                     replacement = 'alcohol')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Brain pH:ch1',
                                                     replacement = 'ph')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Drug abuse (ratings scale:ch1',
                                                     replacement = 'drug_abuse')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Duration of illness (years):ch1',
                                                     replacement = 'duration_illness')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Electroconvulsive therapy:ch1',
                                                     replacement = 'ECT')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Fluphenazine mg. Equivalents:ch1',
                                                     replacement = 'fluphenazine')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Gender:ch1',
                                                     replacement = 'gender')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Lithium treatment:ch1',
                                                     replacement = 'lithium')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Post mortem interval (hours):ch1',
                                                     replacement = 'pmi')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Side of brain:ch1',
                                                     replacement = 'side')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Suicide:ch1',
                                                     replacement = 'suicide')
colnames(geo_GSE5392_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE5392_pheno),
                                                     pattern = 'Valproate treatment:ch1',
                                                     replacement = 'valproate')

#variable transformation
geo_GSE5392_pheno$diagnosis = factor(geo_GSE5392_pheno$diagnosis, levels = c("control", "bipolar disorder"))
geo_GSE5392_pheno$age = as.numeric(geo_GSE5392_pheno$age)
geo_GSE5392_pheno$age_onset = as.numeric(geo_GSE5392_pheno$age_onset)
geo_GSE5392_pheno$ph = as.numeric(geo_GSE5392_pheno$ph)
geo_GSE5392_pheno$gender = factor(geo_GSE5392_pheno$gender)
geo_GSE5392_pheno$pmi = as.numeric(geo_GSE5392_pheno$pmi)
geo_GSE5392_pheno$side = factor(geo_GSE5392_pheno$side)
geo_GSE5392_pheno$suicide = factor(geo_GSE5392_pheno$suicide)
geo_GSE5392_pheno$valproate = factor(geo_GSE5392_pheno$valproate)
geo_GSE5392_pheno$lithium = factor(geo_GSE5392_pheno$lithium)

geo_GSE5392_pheno$drug_abuse = gsub("0 \\(no use\\) to 5 \\(heavy use\\)\\): ", "", geo_GSE5392_pheno$drug_abuse)
geo_GSE5392_pheno$alcohol = gsub("0 \\(no use\\) to 6 \\(heavy use\\)\\): ", "", geo_GSE5392_pheno$alcohol)
geo_GSE5392_pheno$fluphenazine = gsub("Unknown", "NA", geo_GSE5392_pheno$fluphenazine)
geo_GSE5392_pheno$alcohol = gsub("Unknown", "NA", geo_GSE5392_pheno$alcohol)
geo_GSE5392_pheno$drug_abuse = gsub("Unknown", "NA", geo_GSE5392_pheno$drug_abuse)
geo_GSE5392_pheno$ECT = gsub("Yes \\(60 treatments)", "Yes", geo_GSE5392_pheno$ECT)


geo_GSE5392_pheno$ECT = factor(geo_GSE5392_pheno$ECT)
geo_GSE5392_pheno$duration_illness = as.numeric(geo_GSE5392_pheno$duration_illness)
geo_GSE5392_pheno$fluphenazine = as.numeric(geo_GSE5392_pheno$fluphenazine)
geo_GSE5392_pheno$alcohol = factor(geo_GSE5392_pheno$alcohol)
geo_GSE5392_pheno$drug_abuse = factor(geo_GSE5392_pheno$drug_abuse)


###The dataset contains two tissues(dorsolateral prefrontal cortex - pfc, orbitofrontal cortex)
##dorsolateral pre-frontal cortex

pfc = geo_GSE5392_pheno %>%
  filter(grepl("DLPFC", title,ignore.case = T))
pfc_ids = rownames(pfc)


#read the raw data
files_pfc = files[grepl(paste(pfc_ids, collapse = "|"), files)]
cel_paths_pfc = paste0("GSE5392_RAW/", files_pfc)

rawdata_GSE5392_pfc = affy::ReadAffy(filenames = cel_paths_pfc)

#normalization
eset_GSE5392_pfc = affy::expresso(afbatch = rawdata_GSE5392_pfc, 
                              bgcorrect.method = "rma",
                              pmcorrect.method = "pmonly",
                              summary.method = "medianpolish")
eset_GSE5392_pfc_expression = affy::exprs(eset_GSE5392_pfc)

#clean eset colnames
colnames(eset_GSE5392_pfc_expression) = sapply(colnames(eset_GSE5392_pfc_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})
all(colnames(eset_GSE5392_pfc_expression) == rownames(geo_GSE5392_pheno$geo_accession)) #TRUE

##DE NO COVARS
#create design matrix
design0 = model.matrix(~diagnosis, data = pfc)
all(rownames(design0) == rownames(pfc)) #TRUE

#fit the model
fit0 = lmFit(eset_GSE5392_pfc_expression, design0)
fit0 = eBayes(fit0)
GSE5392_topTable_pfc_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE5392_topTable_pfc_NO_covars$ID = rownames(GSE5392_topTable_pfc_NO_covars)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "ID")
GSE5392_topTable_pfc_NO_covars = inner_join(GSE5392_topTable_pfc_NO_covars, SE0, by="ID")

#annotate the table
#get annotation file
GSE5392_probes = geo_GSE5392@featureData@data
GSE5392_topTable_pfc_NO_covars = inner_join(GSE5392_topTable_pfc_NO_covars, GSE5392_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE5392_topTable_pfc_NO_covars = GSE5392_topTable_pfc_NO_covars[, c(9,11, 1:8, 10)]

any(is.na(GSE5392_topTable_pfc_NO_covars$`Gene Symbol`)) #FALSE

write.csv(GSE5392_topTable_pfc_NO_covars, "GSE5392_table_pfc_NO_covars.csv")

##DE with covars
design = model.matrix(~diagnosis+age+ph+gender+pmi+side+suicide,data= pfc)
#since additional covariates cause loose of more than 50% of data
#ECT, lithium, valproate cause more than 50% of data loss


#fit the model 
fit = lmFit(eset_GSE5392_pfc_expression, design)
fit = eBayes(fit)
GSE5392_topTable_pfc_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE5392_topTable_pfc_WITH_covars$ID = rownames(GSE5392_topTable_pfc_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 9)]
colnames(SE) = c("SE", "ID")
GSE5392_topTable_pfc_WITH_covars = inner_join(GSE5392_topTable_pfc_WITH_covars, SE, by = "ID")

#annotate the table
GSE5392_topTable_pfc_WITH_covars = inner_join(GSE5392_topTable_pfc_WITH_covars, GSE5392_probes[, c("Gene Symbol", "ID")], by="ID")
GSE5392_topTable_pfc_WITH_covars = GSE5392_topTable_pfc_WITH_covars[, c(9,11, 1:8, 10)]
any(is.na(GSE5392_topTable_pfc_WITH_covars$`Gene Symbol`)) #FALSE

write.csv(GSE5392_topTable_pfc_WITH_covars, "GSE5392_table_pfc_WITH_covars.csv")





#orbitofrontal cortex - ofc
ofc = geo_GSE5392_pheno %>%
  filter(grepl("OFC", title,ignore.case = T))
ofc_ids = rownames(ofc)

#read the raw data
files_ofc = files[grepl(paste(ofc_ids, collapse = "|"), files)]
cel_paths_ofc = paste0("GSE5392_RAW/", files_ofc)

rawdata_GSE5392_ofc = affy::ReadAffy(filenames = cel_paths_ofc)

#normalization
eset_GSE5392_ofc = affy::expresso(afbatch = rawdata_GSE5392_ofc, 
                                  bgcorrect.method = "rma",
                                  pmcorrect.method = "pmonly",
                                  summary.method = "medianpolish")
eset_GSE5392_ofc_expression = affy::exprs(eset_GSE5392_ofc)

#clean eset colnames
colnames(eset_GSE5392_ofc_expression) = sapply(colnames(eset_GSE5392_ofc_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})
all(colnames(eset_GSE5392_ofc_expression) == rownames(geo_GSE5392_pheno$geo_accession)) #TRUE

##DE NO COVARS
#create design matrix
design0 = model.matrix(~diagnosis, data = ofc)
all(rownames(design0) == rownames(ofc)) #TRUE

#fit the model
fit0 = lmFit(eset_GSE5392_ofc_expression, design0)
fit0 = eBayes(fit0)
GSE5392_topTable_ofc_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE5392_topTable_ofc_NO_covars$ID = rownames(GSE5392_topTable_ofc_NO_covars)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "ID")
GSE5392_topTable_ofc_NO_covars = inner_join(GSE5392_topTable_ofc_NO_covars, SE0, by="ID")

#annotate the table
#get annotation file
GSE5392_probes = geo_GSE5392@featureData@data
GSE5392_topTable_ofc_NO_covars = inner_join(GSE5392_topTable_ofc_NO_covars, GSE5392_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE5392_topTable_ofc_NO_covars = GSE5392_topTable_ofc_NO_covars[, c(9,11, 1:8, 10)]

any(is.na(GSE5392_topTable_ofc_NO_covars$`Gene Symbol`)) #FALSE

write.csv(GSE5392_topTable_ofc_NO_covars, "GSE5392_table_ofc_NO_covars.csv")

##DE with covars
design = model.matrix(~diagnosis+age+ph+gender+pmi+side+suicide,data= ofc)
#since additional covariates cause loose of more than 50% of data
#ECT, lithium, valporate cause more than 50% of data loss


#fit the model 
fit = lmFit(eset_GSE5392_ofc_expression, design)
fit = eBayes(fit)
GSE5392_topTable_ofc_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE5392_topTable_ofc_WITH_covars$ID = rownames(GSE5392_topTable_ofc_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 9)]
colnames(SE) = c("SE", "ID")
GSE5392_topTable_ofc_WITH_covars = inner_join(GSE5392_topTable_ofc_WITH_covars, SE, by = "ID")

#annotate the table
GSE5392_topTable_ofc_WITH_covars = inner_join(GSE5392_topTable_ofc_WITH_covars, GSE5392_probes[, c("Gene Symbol", "ID")], by="ID")
GSE5392_topTable_ofc_WITH_covars = GSE5392_topTable_ofc_WITH_covars[, c(9,11, 1:8, 10)]
any(is.na(GSE5392_topTable_ofc_WITH_covars$`Gene Symbol`)) #FALSE

write.csv(GSE5392_topTable_ofc_WITH_covars, "GSE5392_table_ofc_WITH_covars.csv")

#clean the tables
#pfc no covars
#delete rows with more than one gene mapped
GSE5392_table_pfc_NO_covars_cleaned = GSE5392_table_pfc_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
pfc_NO_covars_cleaned_unique_copy = GSE5392_table_pfc_NO_covars_cleaned
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE5392_table_pfc_NO_covars_cleaned = GSE5392_table_pfc_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE5392_table_pfc_NO_covars_cleaned$Study = "GSE5392"
write.csv(GSE5392_table_pfc_NO_covars_cleaned, "GSE5392_pfc_no_covars_cleaned.csv")



#pfc with covars
#delete rows with more than one gene mapped
GSE5392_table_pfc_WITH_covars_cleaned = GSE5392_table_pfc_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
pfc_WITH_covars_cleaned_unique_copy = GSE5392_table_pfc_WITH_covars_cleaned
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE5392_table_pfc_WITH_covars_cleaned = GSE5392_table_pfc_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE5392_table_pfc_WITH_covars_cleaned$Study = "GSE5392"
write.csv(GSE5392_table_pfc_WITH_covars_cleaned, "GSE5392_pfc_with_covars_cleaned.csv")


##ofc no covars
#delete rows with more than one gene mapped
GSE5392_table_ofc_NO_covars_cleaned = GSE5392_table_ofc_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
ofc_NO_covars_cleaned_copy = GSE5392_table_ofc_NO_covars_cleaned
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE5392_table_ofc_NO_covars_cleaned = GSE5392_table_ofc_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE5392_table_ofc_NO_covars_cleaned$Study = "GSE5392"
write.csv(GSE5392_table_ofc_NO_covars_cleaned, "GSE5392_ofc_no_covars_cleaned.csv")



#ofc with covars
#delete rows with more than one gene mapped
GSE5392_table_ofc_WITH_covars_cleaned = GSE5392_table_ofc_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
ofc_WITH_covars_cleaned_copy = GSE5392_table_ofc_WITH_covars_cleaned
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE5392_table_ofc_WITH_covars_cleaned = GSE5392_table_ofc_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  )
GSE5392_table_ofc_WITH_covars_cleaned$Study = "GSE5392"

write.csv(GSE5392_table_ofc_WITH_covars_cleaned, "GSE5392_ofc_with_covars_cleaned.csv")



















