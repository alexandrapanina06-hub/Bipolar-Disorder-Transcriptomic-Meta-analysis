##GSE87610 "Affymetrix Human Genome U219 Array

library(edgeR)
library(R.utils)
library(dplyr)
library(affy)
library(limma)
rm(list=ls())
#get the files
geo_GSE87610 = GEOquery::getGEO("GSE87610")
geo_GSE87610 = geo_GSE87610[[1]]
geo_GSE87610_pheno = geo_GSE87610@phenoData@data

##filter only bipolar disorder and controls
bd_control_filter = geo_GSE87610_pheno %>% 
  filter(grepl("Bipolar", characteristics_ch1.3, ignore.case = TRUE) | 
           grepl("Unaffected Comparison subject", characteristics_ch1.3, ignore.case = TRUE))
bd_control_filter_ids = rownames(bd_control_filter)

#get supplementary files
files = list.files('GSE87610_RAW/') #get the files (untar)
paths = paste0('GSE87610_RAW/', files) #paths to files
sapply(paths, gunzip) #unzip

files = list.files("GSE87610_RAW")
cel_paths = paste0("GSE87610_RAW/", files)

filtered_files = files[grepl(paste(bd_control_filter_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE87610_RAW/", filtered_files)

#read raw data
rawdata_GSE87610 = affy::ReadAffy(filenames = cel_paths_filtered)


degradation = AffyRNAdeg(rawdata_GSE87610)
plotAffyRNAdeg(degradation, transform="shift.scale", cols=NULL) 
#slope is consistent

##data processing
Eset_GSE87610 = affy::expresso(afbatch = rawdata_GSE87610,
                               bgcorrect.method = "rma",
                               pmcorrect.method = "pmonly",
                               summary.method = "medianpolish")
Eset_GSE87610_expression = affy::exprs(Eset_GSE87610)

colnames(Eset_GSE87610_expression) = sapply(colnames(Eset_GSE87610_expression), 
  function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

#check alignment
all(colnames(Eset_GSE87610_expression) == colnames(geo_GSE87610_pheno$geo_accession)) #TRUE

#differential expression using limma
bd_control_filter = bd_control_filter %>% mutate(Group = case_when(
  grepl("Bipolar", `genotype:ch1`) ~ "bipolar disorder", 
  grepl("Unaffected Comparison subject", `genotype:ch1`) ~ "control", 
  TRUE ~ "Other"
))

bd_control_filter$Group = factor(bd_control_filter$Group, levels = c("control", "bipolar disorder"))

#differential expression WITHOUT COVARIATES
design = model.matrix(~Group, data = bd_control_filter)

#allignment check
all(rownames(design) == rownames(bd_control_filter)) #TRUE 

#fit the model
fit = lmFit(Eset_GSE87610_expression, design = design)
fit = eBayes(fit)
GSE_87610_topTable_NO_covar = limma::topTable(fit = fit, coef=2, adjust.method = 'fdr', number = Inf, confint = T)
GSE_87610_topTable_NO_covar$ID = rownames(GSE_87610_topTable_NO_covar)
#add SE to topTable
SE = as.data.frame(sqrt(fit$s2.post) * fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c("Groupbipolar disorder", "ID")]
colnames(SE) = c("SE", "ID")
GSE_87610_topTable_NO_covar = inner_join(GSE_87610_topTable_NO_covar, SE, by="ID")

#annotate
GSE87610_probes = geo_GSE87610@featureData@data
GSE_87610_topTable_NO_covar = inner_join(GSE_87610_topTable_NO_covar, GSE87610_probes[, c("Gene Symbol", "ID")], by="ID")

#rearrange columns in topTable
GSE_87610_topTable_NO_covar = GSE_87610_topTable_NO_covar[, c(9, 11, 1:8, 10)]
any(is.na(GSE_87610_topTable_NO_covar$`Gene Symbol`)) #FALSE


#save results 
write.csv(GSE_87610_topTable_NO_covar, file = "GSE87610_table_NO_covars.csv")

#clean the tables
#no covars 
#delete rows with more than one gene mapped
GSE87610_table_NO_covars_cleaned = GSE87610_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE87610_table_NO_covars_cleaned = GSE87610_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE87610_table_NO_covars_cleaned = GSE87610_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE87610_table_NO_covars_cleaned$Study = "GSE87610"
write_csv(GSE87610_table_NO_covars_cleaned, "GSE87610_no_covars_cleaned.csv")











