rm(list=ls())
##GSE12649 Affymetrix human gene U133 array

#get the file
geo_GSE12649 = GEOquery::getGEO("GSE12649")
geo_GSE12649 = geo_GSE12649[[1]]
geo_GSE12649_pheno = geo_GSE12649@phenoData@data

#add diagnosis column
geo_GSE12649_pheno = geo_GSE12649_pheno %>%
  mutate(diagnosis = case_when(
    grepl("bipolar", source_name_ch1) ~ "bipolar disorder",
    grepl("schizophrenia", source_name_ch1) ~ "schizophrenia",
    grepl("control", source_name_ch1) ~ "control",
    TRUE ~ "Other"
  ))

#filter bipolar disorder and controls
BD_data = geo_GSE12649_pheno %>%
  filter(grepl("bipolar disorder", diagnosis, ignore.case = T)|
           grepl("control", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)

#get supplementary files
files = list.files("GSE12649_RAW")
paths = paste0("GSE12649_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE12649_RAW")
cel_paths = paste0("GSE12649_RAW/", files)

files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE12649_RAW/", files_filtered)

#read raw data
rawdata_GSE12649 = affy::ReadAffy(filenames = cel_paths_filtered)

#quality check
affy::boxplot(rawdata_GSE12649)
affy::hist(rawdata_GSE12649)

degradation = affy::AffyRNAdeg(rawdata_GSE12649)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL)

#normalize the data
eset_GSE12649 = affy::expresso(afbatch = rawdata_GSE12649, 
                               bgcorrect.method = "rma",
                               pmcorrect.method = "pmonly",
                               summary.method = "medianpolish")
eset_GSE12649_expression = affy::exprs(eset_GSE12649)
colnames(eset_GSE12649_expression) = sapply(colnames(eset_GSE12649_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="."))
  x = x[1]
  return(x)
})

all(colnames(eset_GSE12649_expression)==rownames(BD_data)) #TRUE

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))

##DE NO COVARS
#create design matrix
design = model.matrix(~diagnosis, data=BD_data)
all(rownames(design)==rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE12649_expression, design)
fit = eBayes(fit)
GSE12649_topTable_NO_covars = topTable(fit=fit, coef=2, adjust.method = "fdr", number = Inf, confint = T)
GSE12649_topTable_NO_covars$ID = rownames(GSE12649_topTable_NO_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "ID")
GSE12649_topTable_NO_covars = inner_join(GSE12649_topTable_NO_covars, SE, by="ID")

#annotate the table
GSE12649_probes = geo_GSE12649@featureData@data
GSE12649_topTable_NO_covars = inner_join(GSE12649_topTable_NO_covars, GSE12649_probes[, c("Gene Symbol", "ID")], by = "ID")
GSE12649_topTable_NO_covars = GSE12649_topTable_NO_covars[, c(9,11,1:8, 10)]

#remove rows with missing gene symbols
GSE12649_topTable_NO_covars$`Gene Symbol`[GSE12649_topTable_NO_covars$`Gene Symbol`==""] = NA
GSE12649_topTable_NO_covars_filtered = GSE12649_topTable_NO_covars[!is.na(GSE12649_topTable_NO_covars$`Gene Symbol`),]

#save the results
write.csv(GSE12649_topTable_NO_covars_filtered, "GSE12649_table_NO_covars.csv")

#clean the tables
GSE12649_table_NO_covars_cleaned = GSE12649_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
cleaned_copy = GSE12649_table_NO_covars_cleaned
#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE12649_table_NO_covars_cleaned = GSE12649_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE12649_table_NO_covars_cleaned$Study = "GSE12649"
write.csv(GSE12649_table_NO_covars_cleaned, "GSE12649_no_covars_cleaned.csv")








