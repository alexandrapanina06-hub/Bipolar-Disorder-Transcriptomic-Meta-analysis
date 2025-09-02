rm(list=ls())
library(affy)

###GSE12679 Affymetrix Human U133
#get the file

geo_GSE12679 = GEOquery::getGEO("GSE12679")
geo_GSE12679 = geo_GSE12679[[1]]
geo_GSE12679_pheno = geo_GSE12679@phenoData@data

#create separate columns for each characteristics (diagnosis, cell type, age, sex, pmi)
geo_GSE12679_pheno = geo_GSE12679_pheno %>%
  mutate(diagnosis = case_when(
    grepl("bipolar", `Cell type:ch1`) ~ "bipolar disorder",
    grepl("schizophrenia", `Cell type:ch1`) ~ "schizophrenia",
    grepl("control", `Cell type:ch1`) ~ "control",
    TRUE ~ "Other"
  ))

geo_GSE12679_pheno = geo_GSE12679_pheno %>%
  mutate(cell = case_when(
    grepl("neuronal", `Cell type:ch1`) ~ "neuronal",
    grepl("endothelial", `Cell type:ch1`) ~ "endothelial",
    TRUE ~ "Other"
  ))

geo_GSE12679_pheno = geo_GSE12679_pheno %>%
  mutate(sex = case_when(
    grepl("gender:male", `Cell type:ch1`) ~ "male",
    grepl("gender:female", `Cell type:ch1`) ~ "female",
    TRUE ~ "Other"
  ))

geo_GSE12679_pheno$age = sapply(geo_GSE12679_pheno$`Cell type:ch1`, function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern=","))
  x = x[2]
  return(x)
})

geo_GSE12679_pheno$pmi = sapply(geo_GSE12679_pheno$`Cell type:ch1`, function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern=","))
  x = x[5]
  return(x)
})

geo_GSE12679_pheno$age = gsub("age:", "", geo_GSE12679_pheno$age)
geo_GSE12679_pheno$pmi = gsub("PMI:", "", geo_GSE12679_pheno$pmi)
geo_GSE12679_pheno$pmi = gsub("h", "", geo_GSE12679_pheno$pmi)

#filter only BD and controls
BD_data = geo_GSE12679_pheno %>%
  filter(grepl("bipolar disorder", diagnosis, ignore.case = T)|
           grepl("control", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)

#get supplementary files
files = list.files("GSE12679_RAW")
paths = paste0("GSE12679_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE12679_RAW")
cel_paths = paste0("GSE12679_RAW/", files)

files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE12679_RAW/", files_filtered)

#read raw data
rawdata_GSE12679 = affy::ReadAffy(filenames = cel_paths_filtered)

#quality check

affy::boxplot(rawdata_GSE12679)
affy::hist(rawdata_GSE12679)

degradation = affy::AffyRNAdeg(rawdata_GSE12679)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL)

eset_GSE12679 = affy::expresso(afbatch = rawdata_GSE12679, 
                               bgcorrect.method = "rma",
                               pmcorrect.method = "pmonly",
                               summary.method = "medianpolish")

eset_GSE12679_expression = affy::exprs(eset_GSE12679)
colnames(eset_GSE12679_expression) = sapply(colnames(eset_GSE12679_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE12679_expression) == rownames(BD_data)) #TRUE

#transform variables in BD_data
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$cell = factor(BD_data$cell)
BD_data$sex = factor(BD_data$sex)
BD_data$age = as.numeric(BD_data$age)
BD_data$pmi = as.numeric(BD_data$pmi)

#DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis + age+sex+pmi+cell, data = BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE12679_expression, design)
fit = eBayes(fit)
GSE21679_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE21679_topTable_WITH_covars$ID = rownames(GSE21679_topTable_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,7)]
colnames(SE) = c("SE", "ID")
GSE21679_topTable_WITH_covars = inner_join(GSE21679_topTable_WITH_covars, SE, by = "ID")

#annotate the table
#get annotation file
GSE21679_probes = geo_GSE12679@featureData@data
GSE21679_topTable_WITH_covars = inner_join(GSE21679_topTable_WITH_covars, GSE21679_probes[, c("Gene Symbol", "ID")], by="ID")
GSE21679_topTable_WITH_covars = GSE21679_topTable_WITH_covars[, c(9, 11, 1:8, 10)]

#remove rows with empty gene symbol
GSE21679_topTable_WITH_covars$`Gene Symbol`[GSE21679_topTable_WITH_covars$`Gene Symbol`==""]=NA
GSE21679_topTable_WITH_covars_filtered = GSE21679_topTable_WITH_covars[!is.na(GSE21679_topTable_WITH_covars$`Gene Symbol`),]

#DE NO COVARS
#create design matrix
design_0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(design_0) == rownames(BD_data)) #TRUE

#fit the model
fit_0 = lmFit(eset_GSE12679_expression, design_0)
fit_0 = eBayes(fit_0)
GSE21679_topTable_NO_covars = topTable(fit = fit_0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE21679_topTable_NO_covars$ID = rownames(GSE21679_topTable_NO_covars)

#add SE to teh table
SE_0 = as.data.frame(sqrt(fit_0$s2.post)*fit_0$stdev.unscaled)
SE_0$ID = rownames(SE_0)
SE_0 = SE_0[, 2:3]
colnames(SE_0) = c("SE", "ID")
GSE21679_topTable_NO_covars = inner_join(GSE21679_topTable_NO_covars, SE_0, by = "ID")

#annotate the table
GSE21679_topTable_NO_covars = inner_join(GSE21679_topTable_NO_covars, GSE21679_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE21679_topTable_NO_covars = GSE21679_topTable_NO_covars[, c(9, 11, 1:8, 10)]

#remove rows with empty gene symbol 
GSE21679_topTable_NO_covars$`Gene Symbol`[GSE21679_topTable_NO_covars$`Gene Symbol`==""] = NA
GSE21679_topTable_NO_covars_filtered = GSE21679_topTable_NO_covars[!is.na(GSE21679_topTable_NO_covars$`Gene Symbol`),]

#save results
write.csv(GSE21679_topTable_NO_covars_filtered, "GSE12679_table_NO_covars.csv")
write.csv(GSE21679_topTable_WITH_covars_filtered, "GSE12679_table_WITH_covars.csv")

#clean the tables
#delete rows with more than one gene mapped
GSE12679_table_NO_covars_cleaned = GSE12679_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE12679_table_NO_covars_cleaned = GSE12679_table_NO_covars_cleaned
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE12679_table_NO_covars_cleaned = GSE12679_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE12679_table_NO_covars_cleaned$Study = "GSE12679"
write.csv(GSE12679_table_NO_covars_cleaned, "GSE12679_no_covars_cleaned.csv")


#delete rows with more than one gene mapped
GSE12679_table_WITH_covars_cleaned = GSE12679_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE12679_table_WITH_covars_cleaned = GSE12679_table_WITH_covars_cleaned


#modify logFC and SE
GSE12679_table_WITH_covars_cleaned = GSE12679_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE12679_table_WITH_covars_cleaned$Study = "GSE12679"
write.csv(GSE12679_table_WITH_covars_cleaned, "GSE12679_with_covars_cleaned.csv")




