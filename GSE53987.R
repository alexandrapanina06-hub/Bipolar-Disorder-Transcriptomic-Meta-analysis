rm(list=ls())

##GSE53987 Affymetrix U133
#get the file
geo_GSE53987 = GEOquery::getGEO("GSE53987")
geo_GSE53987 = geo_GSE53987[[1]]
geo_GSE53987_pheno = geo_GSE53987@phenoData@data

#get supplementary files
files = list.files('GSE53987_RAW/') #get the files (untar)
paths = paste0('GSE53987_RAW/', files) #paths to files
sapply(paths, gunzip) #unzip

files = list.files("GSE53987_RAW")
cel_paths = paste0("GSE53987_RAW/", files)

#filter only bipolar disorder and controls
BD_data = geo_GSE53987_pheno %>%
  filter(grepl("bipolar disorder", `disease state:ch1`, ignore.case = T)|
           grepl("control", `disease state:ch1`, ignore.case = T))
BD_data_ids = rownames(BD_data)

filtered_files = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE53987_RAW/", filtered_files)

#read raw data
rawdata_GSE53987 = affy::ReadAffy(filenames = cel_paths_filtered)

#quality check
affy::boxplot(rawdata_GSE53987)
affy::hist(rawdata_GSE53987)

degradation = AffyRNAdeg(rawdata_GSE53987)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL) #slope is consistent

eset_GSE53987 = affy::expresso(afbatch = rawdata_GSE53987, 
                               bgcorrect.method = "rma",
                               pmcorrect.method = "pmonly",
                               summary.method = "medianpolish")
eset_GSE53987_expression = affy::exprs(eset_GSE53987)

colnames(eset_GSE53987_expression) = sapply(colnames(eset_GSE53987_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE53987_expression) == rownames(BD_data$geo_accession)) #TRUE

#transform data into factor and numeric type
BD_data$`disease state:ch1` = factor(BD_data$`disease state:ch1`, levels = c("control", "bipolar disorder"))

BD_data$`gender:ch1` = factor(BD_data$`gender:ch1`, levels = c("M", "F"))
BD_data$`ph:ch1` = as.numeric(BD_data$`ph:ch1`)
BD_data$`rin:ch1` = as.numeric(BD_data$`rin:ch1`)
BD_data$`pmi:ch1` = as.numeric(BD_data$`pmi:ch1`)
BD_data$`tissue:ch1` = factor(BD_data$`tissue:ch1`)
BD_data$`race:ch1` = factor(BD_data$`race:ch1`)
BD_data$`age:ch1` = as.numeric(BD_data$`age:ch1`)

##DE WITHOUT COVARS
design = model.matrix(~`disease state:ch1`, data = BD_data)

#fit the model
fit = lmFit(eset_GSE53987_expression, design)
fit = eBayes(fit)
GSE53987_topTable_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE53987_topTable_NO_covars$ID = rownames(GSE53987_topTable_NO_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "ID")
GSE53987_topTable_NO_covars = inner_join(GSE53987_topTable_NO_covars, SE, by = "ID")

#annotate the table
#get annotation file
GSE53987_probes = geo_GSE53987@featureData@data
GSE53987_topTable_NO_covars = inner_join(GSE53987_topTable_NO_covars, GSE53987_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE53987_topTable_NO_covars = GSE53987_topTable_NO_covars[, c(9, 11, 1:8, 10)]

#filter rows with missing gene symbols
GSE53987_topTable_NO_covars$`Gene Symbol`[GSE53987_topTable_NO_covars$`Gene Symbol` ==""] = NA
GSE53987_topTable_NO_covars_filtered = GSE53987_topTable_NO_covars[!is.na(GSE53987_topTable_NO_covars$`Gene Symbol`),]

## DE WITH COVARS
#make a design matrix with addition of covariates
design_addition = model.matrix(~`disease state:ch1` + `gender:ch1`+`ph:ch1`+`rin:ch1`+
                                 `pmi:ch1` + `tissue:ch1` + `race:ch1`+`age:ch1`, 
                               data = BD_data)
#fit the model
fit_addition = lmFit(eset_GSE53987_expression, design_addition)
fit_addition = eBayes(fit_addition)
GSE53987_topTable_WITH_covars_add = topTable(fit = fit_addition, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE53987_topTable_WITH_covars_add$ID = rownames(GSE53987_topTable_WITH_covars_add)

#add SE to the table
SE_add = as.data.frame(sqrt(fit_addition$s2.post)*fit_addition$stdev.unscaled)
SE_add$ID = rownames(SE_add)
SE_add = SE_add[, c(2,11)]
colnames(SE_add) = c("SE", "ID")
GSE53987_topTable_WITH_covars_add = inner_join(GSE53987_topTable_WITH_covars_add, SE, by = "ID")

#annotate the table
GSE53987_topTable_WITH_covars_add = inner_join(GSE53987_topTable_WITH_covars_add, 
                                               GSE53987_probes[, c("Gene Symbol", "ID")], by="ID" )
GSE53987_topTable_WITH_covars_add = GSE53987_topTable_WITH_covars_add[, c(9,11, 1:8, 10)]

#remove rows with missing gene symbols
GSE53987_topTable_WITH_covars_add$`Gene Symbol`[GSE53987_topTable_WITH_covars_add$`Gene Symbol`==""] = NA
GSE53987_topTable_WITH_covars_add_filtered = GSE53987_topTable_WITH_covars_add[!is.na(GSE53987_topTable_WITH_covars_add$`Gene Symbol`),]

#save the tables
write.csv(GSE53987_topTable_NO_covars_filtered, "GSE53987_table_NO_covars.csv")
write.csv(GSE53987_topTable_WITH_covars_add_filtered, "GSE53987_table_WITH_covars.csv")




###here, there are 3 tissues: hippocampus, pre-frontal cortex, and associative striatum
#DGE in hippocampus

BD_data_hippocampus = BD_data %>%
  filter(grepl("hippocampus", source_name_ch1, ignore.case = T))
BD_data_hippocampus_ids = rownames(BD_data_hippocampus)
files_hippocampus = files[grepl(paste(BD_data_hippocampus_ids, collapse = "|"), files)]
cel_paths_hippocampus = paste0("GSE53987_RAW/", files_hippocampus)

rawdata_GSE53987_hippocampus = affy::ReadAffy(filenames = cel_paths_hippocampus)
eset_GSE53987_hippocampus = affy::expresso(afbatch = rawdata_GSE53987_hippocampus, 
                               bgcorrect.method = "rma",
                               pmcorrect.method = "pmonly",
                               summary.method = "medianpolish")
eset_GSE53987_hippocampus_expression = affy::exprs(eset_GSE53987_hippocampus)

colnames(eset_GSE53987_hippocampus_expression) = sapply(colnames(eset_GSE53987_hippocampus_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE53987_hippocampus_expression) == rownames(BD_data_hippocampus$geo_accession)) #TRUE

#transform data into factor and numeric type
BD_data_hippocampus$`disease state:ch1` = factor(BD_data_hippocampus$`disease state:ch1`, levels = c("control", "bipolar disorder"))
BD_data_hippocampus$`gender:ch1` = factor(BD_data_hippocampus$`gender:ch1`)
BD_data_hippocampus$`ph:ch1` = as.numeric(BD_data_hippocampus$`ph:ch1`)
BD_data_hippocampus$`pmi:ch1` = as.numeric(BD_data_hippocampus$`pmi:ch1`)
BD_data_hippocampus$`race:ch1` = factor(BD_data_hippocampus$`race:ch1`)
BD_data_hippocampus$`age:ch1` = as.numeric(BD_data_hippocampus$`age:ch1`)

##DE WITHOUT COVARS
design = model.matrix(~`disease state:ch1`, data = BD_data_hippocampus)

#fit the model
fit = lmFit(eset_GSE53987_hippocampus_expression, design)
fit = eBayes(fit)
GSE53987_topTable_hippo_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE53987_topTable_hippo_NO_covars$ID = rownames(GSE53987_topTable_hippo_NO_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "ID")
GSE53987_topTable_hippo_NO_covars = inner_join(GSE53987_topTable_hippo_NO_covars, SE, by = "ID")

#annotate the table
#get annotation file
GSE53987_probes = geo_GSE53987@featureData@data
GSE53987_topTable_hippo_NO_covars = inner_join(GSE53987_topTable_hippo_NO_covars, GSE53987_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE53987_topTable_hippo_NO_covars = GSE53987_topTable_hippo_NO_covars[, c(9, 11, 1:8, 10)]

#filter rows with missing gene symbols
GSE53987_topTable_hippo_NO_covars$`Gene Symbol`[GSE53987_topTable_hippo_NO_covars$`Gene Symbol` ==""] = NA
GSE53987_topTable_hippo_NO_covars_filtered = GSE53987_topTable_hippo_NO_covars[!is.na(GSE53987_topTable_hippo_NO_covars$`Gene Symbol`),]

write.csv(GSE53987_topTable_hippo_NO_covars_filtered, "GSE53987_table_hippo_NO_covars.csv")

## DE WITH COVARS
#make a design matrix with addition of covariates
design_addition = model.matrix(~`disease state:ch1` + `gender:ch1`+`ph:ch1`+
                                 `pmi:ch1` + `race:ch1`+`age:ch1`, 
                               data = BD_data_hippocampus)
#fit the model
fit_addition = lmFit(eset_GSE53987_hippocampus_expression, design_addition)
fit_addition = eBayes(fit_addition)
GSE53987_topTable_hippo_WITH_covars = topTable(fit = fit_addition, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE53987_topTable_hippo_WITH_covars$ID = rownames(GSE53987_topTable_hippo_WITH_covars)

#add SE to the table
SE_add = as.data.frame(sqrt(fit_addition$s2.post)*fit_addition$stdev.unscaled)
SE_add$ID = rownames(SE_add)
SE_add = SE_add[, c(2,8)]
colnames(SE_add) = c("SE", "ID")
GSE53987_topTable_hippo_WITH_covars = inner_join(GSE53987_topTable_hippo_WITH_covars, SE, by = "ID")

#annotate the table
GSE53987_topTable_hippo_WITH_covars = inner_join(GSE53987_topTable_hippo_WITH_covars, 
                                               GSE53987_probes[, c("Gene Symbol", "ID")], by="ID" )
GSE53987_topTable_hippo_WITH_covars = GSE53987_topTable_hippo_WITH_covars[, c(9,11, 1:8, 10)]

#remove rows with missing gene symbols
GSE53987_topTable_hippo_WITH_covars$`Gene Symbol`[GSE53987_topTable_hippo_WITH_covars$`Gene Symbol`==""] = NA
GSE53987_topTable_hippo_WITH_covars_filtered = GSE53987_topTable_hippo_WITH_covars[!is.na(GSE53987_topTable_hippo_WITH_covars$`Gene Symbol`),]

write.csv(GSE53987_topTable_hippo_WITH_covars_filtered, "GSE53987_table_hippo_WITH_covars.csv")


##dorsolateral pre-frontal cortex (BA46) - pfc
BD_data_pfc = BD_data %>%
  filter(grepl("(BA46)", source_name_ch1, ignore.case = T))
BD_data_pfc_ids = rownames(BD_data_pfc)


files_pfc = files[grepl(paste(BD_data_pfc_ids, collapse = "|"), files)]
cel_paths_pfc = paste0("GSE53987_RAW/", files_pfc)

rawdata_GSE53987_pfc = affy::ReadAffy(filenames = cel_paths_pfc)
eset_GSE53987_pfc = affy::expresso(afbatch = rawdata_GSE53987_pfc, 
                                           bgcorrect.method = "rma",
                                           pmcorrect.method = "pmonly",
                                           summary.method = "medianpolish")
eset_GSE53987_pfc_expression = affy::exprs(eset_GSE53987_pfc)

colnames(eset_GSE53987_pfc_expression) = sapply(colnames(eset_GSE53987_pfc_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE53987_pfc_expression) == rownames(BD_data_pfc$geo_accession)) #TRUE

#transform data into factor and numeric type
BD_data_pfc$`disease state:ch1` = factor(BD_data_pfc$`disease state:ch1`, levels = c("control", "bipolar disorder"))
BD_data_pfc$`gender:ch1` = factor(BD_data_pfc$`gender:ch1`)
BD_data_pfc$`ph:ch1` = as.numeric(BD_data_pfc$`ph:ch1`)
BD_data_pfc$`pmi:ch1` = as.numeric(BD_data_pfc$`pmi:ch1`)
BD_data_pfc$`race:ch1` = factor(BD_data_pfc$`race:ch1`)
BD_data_pfc$`age:ch1` = as.numeric(BD_data_pfc$`age:ch1`)

##DE WITHOUT COVARS
design = model.matrix(~`disease state:ch1`, data = BD_data_pfc)

#fit the model
fit = lmFit(eset_GSE53987_pfc_expression, design)
fit = eBayes(fit)
GSE53987_topTable_pfc_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE53987_topTable_pfc_NO_covars$ID = rownames(GSE53987_topTable_pfc_NO_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "ID")
GSE53987_topTable_pfc_NO_covars = inner_join(GSE53987_topTable_pfc_NO_covars, SE, by = "ID")

#annotate the table
#get annotation file
GSE53987_probes = geo_GSE53987@featureData@data
GSE53987_topTable_pfc_NO_covars = inner_join(GSE53987_topTable_pfc_NO_covars, GSE53987_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE53987_topTable_pfc_NO_covars = GSE53987_topTable_pfc_NO_covars[, c(9, 11, 1:8, 10)]

#filter rows with missing gene symbols
GSE53987_topTable_pfc_NO_covars$`Gene Symbol`[GSE53987_topTable_pfc_NO_covars$`Gene Symbol` ==""] = NA
GSE53987_topTable_pfc_NO_covars_filtered = GSE53987_topTable_pfc_NO_covars[!is.na(GSE53987_topTable_pfc_NO_covars$`Gene Symbol`),]

write.csv(GSE53987_topTable_pfc_NO_covars_filtered, "GSE53987_table_pfc_NO_covars.csv")

## DE WITH COVARS
#make a design matrix with addition of covariates
design_addition = model.matrix(~`disease state:ch1` + `gender:ch1`+`ph:ch1`+
                                 `pmi:ch1` + `race:ch1`+`age:ch1`, 
                               data = BD_data_pfc)
#fit the model
fit_addition = lmFit(eset_GSE53987_pfc_expression, design_addition)
fit_addition = eBayes(fit_addition)
GSE53987_topTable_pfc_WITH_covars = topTable(fit = fit_addition, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE53987_topTable_pfc_WITH_covars$ID = rownames(GSE53987_topTable_pfc_WITH_covars)

#add SE to the table
SE_add = as.data.frame(sqrt(fit_addition$s2.post)*fit_addition$stdev.unscaled)
SE_add$ID = rownames(SE_add)
SE_add = SE_add[, c(2,8)]
colnames(SE_add) = c("SE", "ID")
GSE53987_topTable_pfc_WITH_covars = inner_join(GSE53987_topTable_pfc_WITH_covars, SE, by = "ID")

#annotate the table
GSE53987_topTable_pfc_WITH_covars = inner_join(GSE53987_topTable_pfc_WITH_covars, 
                                                 GSE53987_probes[, c("Gene Symbol", "ID")], by="ID" )
GSE53987_topTable_pfc_WITH_covars = GSE53987_topTable_pfc_WITH_covars[, c(9,11, 1:8, 10)]

#remove rows with missing gene symbols
GSE53987_topTable_pfc_WITH_covars$`Gene Symbol`[GSE53987_topTable_pfc_WITH_covars$`Gene Symbol`==""] = NA
GSE53987_topTable_pfc_WITH_covars_filtered = GSE53987_topTable_pfc_WITH_covars[!is.na(GSE53987_topTable_pfc_WITH_covars$`Gene Symbol`),]

write.csv(GSE53987_topTable_pfc_WITH_covars_filtered, "GSE53987_table_pfc_WITH_covars.csv")

#clean and modify the tables 
#pre-frontal cortex - pfc
#no covars 
#delete rows with more than one gene mapped
GSE53987_table_pfc_NO_covars_cleaned = GSE53987_table_pfc_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE53987_table_pfc_NO_covars_cleaned = GSE53987_table_pfc_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE53987_table_pfc_NO_covars_cleaned$Study = "GSE53987"
write.csv(GSE53987_table_pfc_NO_covars_cleaned, "GSE53987_pfc_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE53987_table_pfc_WITH_covars_cleaned = GSE53987_table_pfc_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE53987_table_pfc_WITH_covars_cleaned = GSE53987_table_pfc_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE53987_table_pfc_WITH_covars_cleaned$Study = "GSE53987"
write.csv(GSE53987_table_pfc_WITH_covars_cleaned, "GSE53987_pfc_with_covars_cleaned.csv")


#hippocampus - hippo
#no covars 
#delete rows with more than one gene mapped
GSE53987_table_hippo_NO_covars_cleaned = GSE53987_table_hippo_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE53987_table_hippo_NO_covars_cleaned = GSE53987_table_hippo_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE53987_table_hippo_NO_covars_cleaned$Study = "GSE53987"
write.csv(GSE53987_table_hippo_NO_covars_cleaned, "GSE53987_hippo_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE53987_table_hippo_WITH_covars_cleaned = GSE53987_table_hippo_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE53987_table_hippo_WITH_covars_cleaned = GSE53987_table_hippo_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE53987_table_hippo_WITH_covars_cleaned$Study = "GSE53987"
write.csv(GSE53987_table_hippo_WITH_covars_cleaned, "GSE53987_hippo_with_covars_cleaned.csv")












