rm(list=ls())
library(affy)
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
##GSE120340 Affymetrix human gene U133 array
#get the file
geo_GSE120340 = GEOquery::getGEO("GSE120340")
geo_GSE120340 = geo_GSE120340[[1]]
geo_GSE120340_pheno = geo_GSE120340@phenoData@data

geo_GSE120340_pheno = geo_GSE120340_pheno %>%
  mutate(diagnosis = case_when(
    grepl("control", `disease state:ch1`) ~ "control",
    grepl("BD", `disease state:ch1`) ~ "bipolar disorder", 
    TRUE ~ "Other"
  ))

#filter only BD and controls

BD_data = geo_GSE120340_pheno %>%
  filter(grepl("control", diagnosis, ignore.case = T)|
           grepl("bipolar disorder", diagnosis, ignore.case = T))

BD_data_ids = rownames(BD_data)

#get supplementary files
files = list.files("GSE120340_RAW")
paths = paste0("GSE120340_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE120340_RAW")
cel_paths = paste0("GSE120340_RAW/", files)

files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE120340_RAW/", files_filtered)

#read raw data
rawdata_GSE120340 = affy::ReadAffy(filenames = cel_paths_filtered)

#quality check
affy::boxplot(rawdata_GSE120340)
affy::hist(rawdata_GSE120340)

degradation = affy::AffyRNAdeg(rawdata_GSE120340)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL)

#normalization
eset_GSE120340 = affy::expresso(afbatch = rawdata_GSE120340, 
                                bgcorrect.method = "rma",
                                pmcorrect.method = "pmonly",
                                summary.method = "medianpolish")
eset_GSE120340_expression = affy::exprs(eset_GSE120340)

colnames(eset_GSE120340_expression) = sapply(colnames(eset_GSE120340_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

all(colnames(eset_GSE120340_expression) == rownames(BD_data)) #TRUE

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$`laterality:ch1` = factor(BD_data$`laterality:ch1`)
BD_data$`disease state:ch1` = factor(BD_data$`disease state:ch1`, levels = c("control", "BD(-)", "BD(+)"))

#DE NO COVARS
#create design matrix
design_0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(design_0) == rownames(BD_data)) #TRUE

#fit the model
fit_0 = lmFit(eset_GSE120340_expression, design_0)
fit_0 = eBayes(fit_0)
GSE120340_topTable_NO_covars = topTable(fit = fit_0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE120340_topTable_NO_covars$ID = rownames(GSE120340_topTable_NO_covars)

#add SE to teh table
SE_0 = as.data.frame(sqrt(fit_0$s2.post)*fit_0$stdev.unscaled)
SE_0$ID = rownames(SE_0)
SE_0 = SE_0[, 2:3]
colnames(SE_0) = c("SE", "ID")
GSE120340_topTable_NO_covars = inner_join(GSE120340_topTable_NO_covars, SE_0, by = "ID")

#annotate the table
#annotation file
GSE120340_probes = geo_GSE120340@featureData@data

probe_ids = rownames(eset_GSE120340_expression)
annotation_data = select(hgu133plus2.db, 
                         keys = probe_ids, 
                         columns = c("SYMBOL", "ENTREZID", "GENENAME"), 
                         keytype = "PROBEID")
colnames(annotation_data) = c("ID", "Gene Symbol", "ENTREZID", "Gene name")
GSE120340_topTable_NO_covars = inner_join(GSE120340_topTable_NO_covars, annotation_data[, c("Gene Symbol", "ID")], by="ID")
GSE120340_topTable_NO_covars = GSE120340_topTable_NO_covars[, c(9, 11, 1:8, 10)]

#remove rows with missing gene symbol
GSE120340_topTable_NO_covars_filtered = GSE120340_topTable_NO_covars[!is.na(GSE120340_topTable_NO_covars$`Gene Symbol`),]

#DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis + `laterality:ch1`, data = BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE120340_expression, design)
fit = eBayes(fit)
GSE120340_topTable_WITH_covars = topTable(fit=fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE120340_topTable_WITH_covars$ID = rownames(GSE120340_topTable_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,4)]
colnames(SE) = c("SE", "ID")
GSE120340_topTable_WITH_covars = inner_join(GSE120340_topTable_WITH_covars, SE, by="ID")

#annotate the table
GSE120340_topTable_WITH_covars = inner_join(GSE120340_topTable_WITH_covars, annotation_data[, c("Gene Symbol", "ID")], by="ID")
GSE120340_topTable_WITH_covars = GSE120340_topTable_WITH_covars[,c(9,11, 1:8, 10)]

#remove rows with missing gene symbol
GSE120340_topTable_WITH_covars_filtered = GSE120340_topTable_WITH_covars[!is.na(GSE120340_topTable_WITH_covars$`Gene Symbol`),]

#save the results
write.csv(GSE120340_topTable_NO_covars_filtered, file = "GSE120340_table_NO_covars.csv")
write.csv(GSE120340_topTable_WITH_covars_filtered, file = "GSE120340_table_WITH_covars.csv")

#clean the tables
#no covars 
#delete rows with more than one gene mapped
GSE120340_table_NO_covars_cleaned = GSE120340_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE120340_table_NO_covars_cleaned = GSE120340_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE120340_table_NO_covars_cleaned = GSE120340_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE120340_table_NO_covars_cleaned$Study = "GSE120340"
write_csv(GSE120340_table_NO_covars_cleaned, "GSE120340_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE120340_table_WITH_covars_cleaned = GSE120340_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE120340_table_WITH_covars_cleaned = GSE120340_table_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE120340_table_WITH_covars_cleaned = GSE120340_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE120340_table_WITH_covars_cleaned$Study = "GSE120340"
write_csv(GSE120340_table_WITH_covars_cleaned, "GSE120340_with_covars_cleaned.csv")









