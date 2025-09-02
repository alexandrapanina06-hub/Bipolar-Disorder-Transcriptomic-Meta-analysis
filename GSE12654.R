##GSE12654 Affymetrix Human Genome U95 Version 2 Array
rm(list=ls())
#get the file

geo_GSE12654 = GEOquery::getGEO("GSE12654")
geo_GSE12654 = geo_GSE12654[[1]]
geo_GSE12654_pheno = geo_GSE12654@phenoData@data

#add column with diagnosis
geo_GSE12654_pheno = geo_GSE12654_pheno %>%
  mutate(diagnosis = case_when(
    grepl("bipolar disorder",source_name_ch1) ~ "bipolar disorder",
    grepl("depression", source_name_ch1) ~ "depression",
    grepl("control", source_name_ch1) ~ "control",
    grepl("schizophrenia", source_name_ch1) ~ "schizophrenia", 
    TRUE ~ "Other"
  ))

#filter only BD and control data
BD_data = geo_GSE12654_pheno %>%
  filter(grepl("bipolar disorder", diagnosis, ignore.case = T)|
           grepl("control", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)


#get raw data files
files = list.files("GSE12654_RAW")
paths = paste0("GSE12654_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE12654_RAW")
cel_paths = paste0("GSE12654_RAW/", files)

files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE12654_RAW/", files_filtered)

#read raw data
rawdata_GSE12654 = affy::ReadAffy(filenames = cel_paths_filtered)
#Warning: found an empty line where not expected in GSE12654_RAW/GSM317653.CEL.

test_cel = tryCatch({
  affy::ReadAffy(filenames = "GSE12654_RAW/GSM317653.CEL")
}, error = function(e) { e }, warning = function(w) { w })
print(test_cel)
summary(exprs(test_cel))
missing_rows = which(is.na(exprs(test_cel)) | is.infinite(exprs(test_cel)))
print(missing_rows) #character(0)

##since dataset is small and only one cel is missing, let's proceed without removing the file


#quality check
affy::boxplot(rawdata_GSE12654)
affy::hist(rawdata_GSE12654)

#normalization
eset_GSE12654 = affy::expresso(afbatch = rawdata_GSE12654, 
                               bgcorrect.method = "rma",
                               pmcorrect.method = "pmonly",
                               summary.method = "medianpolish")
eset_GSE12654_expression = affy::exprs(eset_GSE12654)

#clear colnames
colnames(eset_GSE12654_expression) = sapply(colnames(eset_GSE12654_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE12654_expression)==rownames(BD_data)) #TRUE

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))

#DE WITHOUT COVARS
#create design matrix
design = model.matrix(~diagnosis, data = BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE12654_expression, design)
fit = eBayes(fit)
GSE12654_topTable_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE12654_topTable_NO_covars$ID = rownames(GSE12654_topTable_NO_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "ID")
GSE12654_topTable_NO_covars = inner_join(GSE12654_topTable_NO_covars, SE, by="ID")

#annotate the table
#get the annotation file
GSE12654_probes = geo_GSE12654@featureData@data
GSE12654_topTable_NO_covars = inner_join(GSE12654_topTable_NO_covars, GSE12654_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE12654_topTable_NO_covars = GSE12654_topTable_NO_covars[, c(9,11,1:8, 10)]
any(is.na(GSE12654_topTable_NO_covars$`Gene Symbol`)) #FALSE

#save the results
write.csv(GSE12654_topTable_NO_covars, "GSE12654_table_NO_covars.csv")

#clean the tables
#delete rows with more than one gene mapped
GSE12654_table_NO_covars_cleaned = GSE12654_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE12654_table_NO_covars_cleaned = GSE12654_table_NO_covars_cleaned
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE12654_table_NO_covars_cleaned = GSE12654_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE12654_table_NO_covars_cleaned$Study = "GSE12654"

write.csv(GSE12654_table_NO_covars_cleaned, "GSE12654_no_covars_cleaned.csv")



