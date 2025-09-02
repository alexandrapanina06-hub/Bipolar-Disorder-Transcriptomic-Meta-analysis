rm(list=ls())
##GSE92538 Affymetrix U133 array
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
geo_GSE92538 = GEOquery::getGEO("GSE92538")
geo_GSE92538 = geo_GSE92538[[1]] #[HG-U133_Plus_2]
geo_GSE92538_pheno = geo_GSE92538@phenoData@data

#add diagnosis column to phenodata
geo_GSE92538_pheno = geo_GSE92538_pheno %>%
  mutate(diagnosis = case_when(
    grepl("Control", `diagnosis:ch1`) ~ "control",
    grepl("Bipolar Disorder", `diagnosis:ch1`) ~ "bipolar disorder",
    grepl("Major Depressive Disorder", `diagnosis:ch1`) ~ "MDD",
    grepl("Schizophrenia", `diagnosis:ch1`) ~ "schizophrenia",
    TRUE ~ "other"
  ))

#filter bipolar disorder and controls
BD_data = geo_GSE92538_pheno %>%
  filter(grepl("control", diagnosis, ignore.case = T)|
           grepl("bipolar disorder", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)


#get the raw data
files = list.files("GSE92538_RAW/")
paths = paste0("GSE92538_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE92538_RAW/")
cel_paths = paste0("GSE92538_RAW/", files)

files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE92538_RAW/", files_filtered)

#read raw data
rawdata_GSE92538 = affy::ReadAffy(filenames = cel_paths_filtered)

#quakity check
affy::boxplot(rawdata_GSE92538)
affy::hist(rawdata_GSE92538)

degradation = affy::AffyRNAdeg(rawdata_GSE92538)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL)

#normalize the data
eset_GSE92538 = affy::expresso(afbatch = rawdata_GSE92538, 
                              bgcorrect.method = "rma",
                              pmcorrect.method = "pmonly",
                              summary.method = "medianpolish")
eset_GSE92538_expression = affy::exprs(eset_GSE92538)

#clean expression colnames
colnames(eset_GSE92538_expression) = sapply(colnames(eset_GSE92538_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

all(colnames(eset_GSE92538_expression) == rownames(geo_GSE92538_pheno$geo_accession)) #TRUE

#clean colnames in BD_data
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                                     pattern = 'age:ch1',
                                                     replacement = 'age')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'agonal factor:ch1',
                                           replacement = 'agonal_factor')

colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'gender:ch1',
                                           replacement = 'gender')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'post-mortem interval:ch1',
                                           replacement = 'pmi')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'race:ch1',
                                           replacement = 'race')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'suicide (1=yes):ch1',
                                           replacement = 'suicide')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'tissue ph (cerebellum):ch1',
                                           replacement = 'ph')
#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, level = c("control", "bipolar disorder"))
BD_data$age = as.numeric(BD_data$age)
BD_data$agonal_factor = factor(BD_data$agonal_factor)
BD_data$gender = factor(BD_data$gender)
BD_data$pmi = as.numeric(BD_data$pmi)
BD_data$race = factor(BD_data$race)
BD_data$suicide = factor(BD_data$suicide)
BD_data$ph = as.numeric(BD_data$ph)

##DE NO COVARS
#create design matrix
design0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(BD_data) == rownames(BD_data)) #TRUE

#fit the model
fit0 = lmFit(eset_GSE92538_expression, design0)
fit0 = eBayes(fit0)
GSE92538_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE92538_topTable_NO_covars$ID = rownames(GSE92538_topTable_NO_covars)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "ID")
GSE92538_topTable_NO_covars = inner_join(GSE92538_topTable_NO_covars, SE0, by = "ID")

#annotate the table
#get the annotation file
GSE92538_probes = geo_GSE92538@featureData@data
probe_ids = rownames(eset_GSE92538_expression)
annotation_data = select(hgu133plus2.db, keys = probe_ids,
                         columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                         keytype = "PROBEID")
GSE92538_topTable_NO_covars = inner_join(GSE92538_topTable_NO_covars, annotation_data[, c("PROBEID", "SYMBOL")], by=c("ID"="PROBEID"))
GSE92538_topTable_NO_covars = GSE92538_topTable_NO_covars[, c(9,11,1:8, 10)]

#remove rows with missing gene symbols
GSE92538_topTable_NO_covars_filtered = GSE92538_topTable_NO_covars[!is.na(GSE92538_topTable_NO_covars$SYMBOL),]

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+age+agonal_factor+gender+pmi+race+suicide+ph, data = BD_data)
#exclude ph because it causes loose of data
all(rownames(design) == rownames(BD_data[!is.na(BD_data$ph),])) #TRUE

eset_GSE92538_expression_filtered = eset_GSE92538_expression[, colnames(eset_GSE92538_expression) %in% rownames(design)]
all(colnames(eset_GSE92538_expression_filtered)==rownames(design)) #TRUE

#fit the model
fit = lmFit(eset_GSE92538_expression_filtered, design)
fit = eBayes(fit)
GSE92538_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE92538_topTable_WITH_covars$ID = rownames(GSE92538_topTable_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,11)]
colnames(SE) = c("SE", "ID")
GSE92538_topTable_WITH_covars = inner_join(GSE92538_topTable_WITH_covars, SE, by="ID")

#annotate the table
GSE92538_topTable_WITH_covars = inner_join(GSE92538_topTable_WITH_covars, annotation_data[, c("SYMBOL", "PROBEID")], by=c("ID" = "PROBEID"))
GSE92538_topTable_WITH_covars = GSE92538_topTable_WITH_covars[, c(9,11,1:8,10)]

#remove rows with missing gene symbol
GSE92538_topTable_WITH_covars_filtered = GSE92538_topTable_WITH_covars[!is.na(GSE92538_topTable_WITH_covars$SYMBOL),]

#save the results
write.csv(GSE92538_topTable_NO_covars_filtered, "GSE92538_table_NO_covars.csv")
write.csv(GSE92538_topTable_WITH_covars_filtered, "GSE92538_table_WITH_covars.csv")

#clean the tables

colnames(GSE92538_table_NO_covars)[3] = "Gene Symbol"
colnames(GSE92538_table_WITH_covars)[3] = "Gene Symbol"

#no covars 
#delete rows with more than one gene mapped
GSE92538_table_NO_covars_cleaned = GSE92538_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE92538_table_NO_covars_cleaned = GSE92538_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE92538_table_NO_covars_cleaned = GSE92538_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE92538_table_NO_covars_cleaned$Study = "GSE92538"
write_csv(GSE92538_table_NO_covars_cleaned, "GSE92538_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE92538_table_WITH_covars_cleaned = GSE92538_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE92538_table_WITH_covars_cleaned = GSE92538_table_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE92538_table_WITH_covars_cleaned = GSE92538_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE92538_table_WITH_covars_cleaned$Study = "GSE92538"
write_csv(GSE92538_table_WITH_covars_cleaned, "GSE92538_with_covars_cleaned.csv")








