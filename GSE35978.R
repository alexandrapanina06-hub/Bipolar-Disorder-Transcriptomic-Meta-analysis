###GSE35978 Affymetrix Human gene ST1 array
rm(list=ls())
library(oligo)
library(dplyr)
library(stringi)
library(huex10sttranscriptcluster.db)
library(limma)

#get the file 
geo_GSE35978 = GEOquery::getGEO("GSE35978")
geo_GSE35978 = geo_GSE35978[[1]]
geo_GSE35978_pheno = geo_GSE35978@phenoData@data

#add grouo column to the pheno data
geo_GSE35978_pheno = geo_GSE35978_pheno %>%
  mutate(Group = case_when(
    grepl("bipolar", `disease status:ch1`) ~ "bipolar disorder",
    grepl("unaffected", `disease status:ch1`) ~ "control",
    grepl("schizophrenia", `disease status:ch1`) ~ "schizophrenia",
    TRUE ~ "Other"
  ))
#filter only bipolar disorder and controls
BD_data = geo_GSE35978_pheno %>%
  filter(grepl("bipolar disorder", Group, ignore.case = T)|
           grepl("control", Group, ignore.case = T))
BD_data_ids = rownames(BD_data)

#get the raw data files 
files = list.files("GSE35978_RAW/")
paths = paste0("GSE35978_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE35978_RAW")
cel_paths = paste0("GSE35978_RAW/", files)

#filter files with BD and controls
filtered_files = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE35978_RAW/", filtered_files)

rawdata_GSE35978 = read.celfiles(filenames = cel_paths_filtered)

#check quality
boxplot(rawdata_GSE35978, target = "core",  main = "Raw Intensity Distribution")

#normalize the data
eset_GSE35978 = rma(rawdata_GSE35978)
eset_GSE35978_expression = exprs(eset_GSE35978)

#clean colnames
colnames(eset_GSE35978_expression) = sapply(colnames(eset_GSE35978_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE35978_expression) == rownames(geo_GSE35978_pheno$geo_accession)) #TRUE

#clean colnames in BD_data
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                                    pattern = 'age:ch1',
                                                    replacement = 'age')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'ph:ch1',
                                           replacement = 'ph')

colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'tissue:ch1',
                                           replacement = 'tissue')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'Sex:ch1',
                                           replacement = 'sex')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'post-mortem interval (pmi):ch1',
                                           replacement = 'pmi')

#variable transformation
BD_data$Group = factor(BD_data$Group, levels = c("control", "bipolar disorder"))
BD_data$tissue = factor(BD_data$tissue)
BD_data$sex = factor(BD_data$sex)
BD_data$ph = as.numeric(BD_data$ph)
BD_data$age = as.numeric(BD_data$age)
BD_data$pmi = as.numeric(BD_data$pmi)

#create annotation file
probe_ids = rownames(eset_GSE35978_expression)
annotation_data = read.delim("GPL6244.an.txt", header = FALSE, sep = "\t", skip = 4, row.names = NULL)
colnames(annotation_data) = c("ID", "Gene Symbol", "Gene name", "GOTerms", "GemmaID", "NCBlids")
annotation_data = annotation_data[6:nrow(annotation_data), ]
annotation_data$`Gene Symbol`[annotation_data$`Gene Symbol`==""] = NA

#DE with covariates
#create design matrix
design = model.matrix(~Group + tissue+sex+ph+age+pmi, data = BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE35978_expression, design)
fit = eBayes(fit)
GSE35978_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE35978_topTable_WITH_covars$ID = rownames(GSE35978_topTable_WITH_covars)

#add SE to teh table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,8)]
colnames(SE) = c("SE", "ID")
GSE35978_topTable_WITH_covars = inner_join(GSE35978_topTable_WITH_covars, SE, by="ID")

#annotate the table
GSE35978_topTable_WITH_covars = inner_join(GSE35978_topTable_WITH_covars, annotation_data, by = "ID")
GSE35978_topTable_WITH_covars = GSE35978_topTable_WITH_covars[, c(9, 11, 1:8, 10)]

#remove rows with missing gene symbol
GSE35978_topTable_WITH_covars_filtered = GSE35978_topTable_WITH_covars[!is.na(GSE35978_topTable_WITH_covars$`Gene Symbol`),]

##DE without covariates
#create design matrix
design_0 = model.matrix(~Group, data = BD_data)
all(rownames(design_0) == rownames(BD_data)) #TRUE

#fit the model
fit_0 = lmFit(eset_GSE35978_expression, design_0)
fit_0 = eBayes(fit_0)
GSE35978_topTable_NO_covars = topTable(fit =fit_0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE35978_topTable_NO_covars$ID = rownames(GSE35978_topTable_NO_covars)

#add SE to the table
SE_0 = as.data.frame(sqrt(fit_0$s2.post)*fit_0$stdev.unscaled)
SE_0$ID = rownames(SE_0)
SE_0 = SE_0[, 2:3]
colnames(SE_0) = c("SE", "ID")

GSE35978_topTable_NO_covars = inner_join(GSE35978_topTable_NO_covars, SE, by = "ID")

#annotate the table
GSE35978_topTable_NO_covars = inner_join(GSE35978_topTable_NO_covars, annotation_data[, c("Gene Symbol", "ID", "Gene name")], by="ID")
GSE35978_topTable_NO_covars = GSE35978_topTable_NO_covars[, c(9,11, 1:8, 10)]

#remove rows with missing gene symbol
GSE35978_topTable_NO_covars_filtered = GSE35978_topTable_NO_covars[!is.na(GSE35978_topTable_NO_covars$`Gene Symbol`),]

#save the results
write.csv(GSE35978_topTable_NO_covars_filtered, "GSE35978_table_NO_covars.csv")
write.csv(GSE35978_topTable_WITH_covars_filtered, "GSE35978_table_WITH_covars.csv")


####since the dataset includes two tissues per person (cerebellum and parietal cortex), let's analyze them independently
BD_data_cerebellum = BD_data %>%
  filter(grepl("cerebellum", tissue, ignore.case = T))
BD_data_cerebellum_ids = rownames(BD_data_cerebellum)

#get the raw data files 
files = list.files("GSE35978_RAW/")
paths = paste0("GSE35978_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE35978_RAW")
cel_paths = paste0("GSE35978_RAW/", files)

#filter files with BD and controls
filtered_files = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths = paste0("GSE35978_RAW/", filtered_files)

files_cerebellum = filtered_files[grepl(paste(BD_data_cerebellum_ids, collapse = "|"), filtered_files)]
cel_paths_cerebellum = paste0("GSE35978_RAW/", files_cerebellum)

rawdata_GSE35978_cerebellum = read.celfiles(filenames = cel_paths_cerebellum)

#check quality
boxplot(rawdata_GSE35978_cerebellum, target = "core",  main = "Raw Intensity Distribution")

#normalize the data
eset_GSE35978_cerebellum = rma(rawdata_GSE35978_cerebellum)
eset_GSE35978_cerebellum_expression = exprs(eset_GSE35978_cerebellum)

#clean colnames
colnames(eset_GSE35978_cerebellum_expression) = sapply(colnames(eset_GSE35978_cerebellum_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE35978_cerebellum_expression) == rownames(geo_GSE35978_pheno$geo_accession)) #TRUE

#DE with covariates
#create design matrix
design = model.matrix(~Group + sex+ph+age+pmi, data = BD_data_cerebellum)
all(rownames(design) == rownames(BD_data_cerebellum)) #TRUE

#fit the model
fit = lmFit(eset_GSE35978_cerebellum_expression, design)
fit = eBayes(fit)
GSE35978_topTable_cerebellum_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE35978_topTable_cerebellum_WITH_covars$ID = rownames(GSE35978_topTable_cerebellum_WITH_covars)

#add SE to teh table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,7)]
colnames(SE) = c("SE", "ID")
GSE35978_topTable_cerebellum_WITH_covars = inner_join(GSE35978_topTable_cerebellum_WITH_covars, SE, by="ID")

#annotate the table
GSE35978_topTable_cerebellum_WITH_covars = inner_join(GSE35978_topTable_cerebellum_WITH_covars, annotation_data, by = "ID")
GSE35978_topTable_cerebellum_WITH_covars = GSE35978_topTable_cerebellum_WITH_covars[, c(9, 11, 1:8, 10)]

#remove rows with missing gene symbol
GSE35978_topTable_cerebellum_WITH_covars = GSE35978_topTable_cerebellum_WITH_covars[!is.na(GSE35978_topTable_cerebellum_WITH_covars$`Gene Symbol`),]

write.csv(GSE35978_topTable_cerebellum_WITH_covars, "GSE35978_table_cerebellum_WITH_covars.csv")

##DE without covariates
#create design matrix
design_0 = model.matrix(~Group, data = BD_data_cerebellum)
all(rownames(design_0) == rownames(BD_data_cerebellum)) #TRUE

#fit the model
fit_0 = lmFit(eset_GSE35978_cerebellum_expression, design_0)
fit_0 = eBayes(fit_0)
GSE35978_topTable_cerebellum_NO_covars = topTable(fit =fit_0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE35978_topTable_cerebellum_NO_covars$ID = rownames(GSE35978_topTable_cerebellum_NO_covars)

#add SE to the table
SE_0 = as.data.frame(sqrt(fit_0$s2.post)*fit_0$stdev.unscaled)
SE_0$ID = rownames(SE_0)
SE_0 = SE_0[, 2:3]
colnames(SE_0) = c("SE", "ID")

GSE35978_topTable_cerebellum_NO_covars = inner_join(GSE35978_topTable_cerebellum_NO_covars, SE, by = "ID")

#annotate the table
GSE35978_topTable_cerebellum_NO_covars = inner_join(GSE35978_topTable_cerebellum_NO_covars, annotation_data[, c("Gene Symbol", "ID", "Gene name")], by="ID")
GSE35978_topTable_cerebellum_NO_covars = GSE35978_topTable_cerebellum_NO_covars[, c(9,11, 1:8, 10)]

#remove rows with missing gene symbol
GSE35978_topTable_cerebellum_NO_covars = GSE35978_topTable_cerebellum_NO_covars[!is.na(GSE35978_topTable_cerebellum_NO_covars$`Gene Symbol`),]

write.csv(GSE35978_topTable_cerebellum_NO_covars, "GSE35978_table_cerebellum_NO_covars.csv")



###parietal cortex
BD_data_parietalcortex = BD_data %>%
  filter(grepl("parietal cortex", tissue, ignore.case = T))
BD_data_parietalcortex_ids = rownames(BD_data_parietalcortex)


#get the raw data files 
files = list.files("GSE35978_RAW/")
paths = paste0("GSE35978_RAW/", files)


files = list.files("GSE35978_RAW")
cel_paths = paste0("GSE35978_RAW/", files)

#filter files with BD and controls
parietalcortex_files = files[grepl(paste(BD_data_parietalcortex_ids, collapse = "|"), files)]
cel_paths_parietalcortex = paste0("GSE35978_RAW/", parietalcortex_files)

rawdata_GSE35978_parietalcortex = read.celfiles(filenames = cel_paths_parietalcortex)

#check quality
boxplot(rawdata_GSE35978_parietalcortex, target = "core",  main = "Raw Intensity Distribution")

#normalize the data
eset_GSE35978_parietalcortex = rma(rawdata_GSE35978_parietalcortex)
eset_GSE35978_parietalcortex_expression = exprs(eset_GSE35978_parietalcortex)

#clean colnames
colnames(eset_GSE35978_parietalcortex_expression) = sapply(colnames(eset_GSE35978_parietalcortex_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE35978_parietalcortex_expression) == rownames(geo_GSE35978_pheno$geo_accession)) #TRUE

#DE with covariates
#create design matrix
design = model.matrix(~Group + sex+ph+age+pmi, data = BD_data_parietalcortex)
all(rownames(design) == rownames(BD_data_parietalcortex)) #TRUE

#fit the model
fit = lmFit(eset_GSE35978_parietalcortex_expression, design)
fit = eBayes(fit)
GSE35978_topTable_parietalcortex_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE35978_topTable_parietalcortex_WITH_covars$ID = rownames(GSE35978_topTable_parietalcortex_WITH_covars)

#add SE to teh table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,7)]
colnames(SE) = c("SE", "ID")
GSE35978_topTable_parietalcortex_WITH_covars = inner_join(GSE35978_topTable_parietalcortex_WITH_covars, SE, by="ID")

#annotate the table
GSE35978_topTable_parietalcortex_WITH_covars = inner_join(GSE35978_topTable_parietalcortex_WITH_covars, annotation_data, by = "ID")
GSE35978_topTable_parietalcortex_WITH_covars = GSE35978_topTable_parietalcortex_WITH_covars[, c(9, 11, 1:8, 10)]

#remove rows with missing gene symbol
GSE35978_topTable_parietalcortex_WITH_covars = GSE35978_topTable_parietalcortex_WITH_covars[!is.na(GSE35978_topTable_parietalcortex_WITH_covars$`Gene Symbol`),]

write.csv(GSE35978_topTable_parietalcortex_WITH_covars, "GSE35978_table_parietalcortex_WITH_covars.csv")

##DE without covariates
#create design matrix
design_0 = model.matrix(~Group, data = BD_data_parietalcortex)
all(rownames(design_0) == rownames(BD_data_parietalcortex)) #TRUE

#fit the model
fit_0 = lmFit(eset_GSE35978_parietalcortex_expression, design_0)
fit_0 = eBayes(fit_0)
GSE35978_topTable_parietalcortex_NO_covars = topTable(fit =fit_0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE35978_topTable_parietalcortex_NO_covars$ID = rownames(GSE35978_topTable_parietalcortex_NO_covars)

#add SE to the table
SE_0 = as.data.frame(sqrt(fit_0$s2.post)*fit_0$stdev.unscaled)
SE_0$ID = rownames(SE_0)
SE_0 = SE_0[, 2:3]
colnames(SE_0) = c("SE", "ID")

GSE35978_topTable_parietalcortex_NO_covars = inner_join(GSE35978_topTable_parietalcortex_NO_covars, SE, by = "ID")

#annotate the table
GSE35978_topTable_parietalcortex_NO_covars = inner_join(GSE35978_topTable_parietalcortex_NO_covars, annotation_data[, c("Gene Symbol", "ID", "Gene name")], by="ID")
GSE35978_topTable_parietalcortex_NO_covars = GSE35978_topTable_parietalcortex_NO_covars[, c(9,11, 1:8, 10)]

#remove rows with missing gene symbol
GSE35978_topTable_parietalcortex_NO_covars = GSE35978_topTable_parietalcortex_NO_covars[!is.na(GSE35978_topTable_parietalcortex_NO_covars$`Gene Symbol`),]

write.csv(GSE35978_topTable_parietalcortex_NO_covars, "GSE35978_table_parietalcortex_NO_covars.csv")


#clean the tables 
#cerebelum
#no covars 
#delete rows with more than one gene mapped
GSE35978_table_cerebellum_NO_covars_cleaned = GSE35978_table_cerebellum_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE35978_table_cerebellum_NO_covars_cleaned = GSE35978_table_cerebellum_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE35978_table_cerebellum_NO_covars_cleaned$Study = "GSE35978"

write.csv(GSE35978_table_cerebellum_NO_covars_cleaned, "GSE35978_cerebellum_no_covars_cleaned.csv")



##with covars
#delete rows with more than one gene mapped
GSE35978_table_cerebellum_WITH_covars_cleaned = GSE35978_table_cerebellum_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE35978_table_cerebellum_WITH_covars_cleaned = GSE35978_table_cerebellum_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE35978_table_cerebellum_WITH_covars_cleaned$Study = "GSE35978"
write.csv(GSE35978_table_cerebellum_WITH_covars_cleaned, "GSE35978_cerebellum_with_covars_cleaned.csv")



##parietal cortex
GSE35978_table_parietalcortex_NO_covars_cleaned = GSE35978_table_parietalcortex_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE35978_table_parietalcortex_NO_covars_cleaned = GSE35978_table_parietalcortex_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE35978_table_parietalcortex_NO_covars_cleaned$Study = "GSE35978"
write.csv(GSE35978_table_parietalcortex_NO_covars_cleaned, "GSE35978_parietalcortex_no_covars_cleaned.csv")



##with covars
#delete rows with more than one gene mapped
GSE35978_table_parietalcortex_WITH_covars_cleaned = GSE35978_table_parietalcortex_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE35978_table_parietalcortex_WITH_covars_cleaned = GSE35978_table_parietalcortex_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE35978_table_parietalcortex_WITH_covars_cleaned$Study = "GSE35978"
write.csv(GSE35978_table_parietalcortex_WITH_covars_cleaned, "GSE35978_parietalcortex_with_covars_cleaned.csv")




