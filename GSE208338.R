##GSE208338 Affymetrix Human ST1 array
rm(list=ls())
library(oligo)
library(GEOquery)
library(huex10sttranscriptcluster.db)
library(dplyr)
library(R.utils)
library(limma)

#get the file
#geo_GSE208338 = read.delim("GSE208338_series_matrix.txt", comment.char = "!", header = TRUE, sep = "\t", quote = "")

library(data.table)
geo_GSE208338 = fread("GSE208338_series_matrix.txt", 
                       skip = "ID_REF",    
                       sep = "\t", 
                       quote = "",
                       data.table = FALSE) 
# Read the metadata separately
metadata_lines = readLines("GSE208338_series_matrix.txt")

# Extract lines that contain phenotype data
pheno_lines = metadata_lines[grepl("^!", metadata_lines)]
pheno_lines = metadata_lines[grepl("^!Sample_", metadata_lines)]
# Check what metadata is available
head(pheno_lines, 20)  # Print the first 20 lines

# Split each line into columns
pheno_list = strsplit(pheno_lines, "\t")

# Convert to a structured dataframe
pheno_df = as.data.frame(do.call(rbind, pheno_list), stringsAsFactors = FALSE)

# Transpose the dataframe
pheno_long = as.data.frame(t(pheno_df), stringsAsFactors = FALSE)

# Set the first row as column names and remove it from data
colnames(pheno_long) = as.character(pheno_long[1, ])
pheno_long = pheno_long[-1, ]

#clean colnames
colnames(pheno_long) = gsub("^!", "", colnames(pheno_long))
colnames(pheno_long) = gsub("Sample_", "", colnames(pheno_long))

cols_to_rename = which(grepl("^characteristics_ch1", colnames(pheno_long)))
colnames(pheno_long)[cols_to_rename] = paste0("characteristics_ch1.", seq_along(cols_to_rename))

cols_rename = which(grepl("^data_processing", colnames(pheno_long)))
colnames(pheno_long)[cols_rename] = paste0("data_processing.", seq_along(cols_rename))

#create diagnosis column
pheno_long = pheno_long %>%
  mutate(diagnosis = case_when(
    grepl("diagnosis: BPD", characteristics_ch1.2) ~ "bipolar disorder",
    grepl("diagnosis: CTL", characteristics_ch1.2) ~ "control",
    grepl("diagnosis: MDD", characteristics_ch1.2) ~ "MDD",
    grepl("diagnosis: SCZ", characteristics_ch1.2) ~ "schizophrenia",
    TRUE ~ "other"
  ))
#filter bipolar disorder and controls
BD_data = pheno_long %>%
  filter(grepl("control", diagnosis, ignore.case = T)|
           grepl("bipolar disorder", diagnosis, ignore.case = T))
BD_data$title = gsub("DLPFC expression from sample ", "", BD_data$title)
BD_data$title = gsub('"', "", BD_data$title)
BD_data_ids = BD_data$geo_accession

#rename columns
colnames(BD_data)[13] = "age"
colnames(BD_data)[14] = "sex"
colnames(BD_data)[15] = "ph"
colnames(BD_data)[16] = "pmi"
colnames(BD_data)[19] = "suicide"
colnames(BD_data)[20] = "death"

#clean columns
BD_data$age = gsub("age: ", "", BD_data$age)
BD_data$age = gsub('"', "", BD_data$age)

BD_data$sex = gsub("Sex: ", "", BD_data$sex)
BD_data$sex = gsub('"', "", BD_data$sex)

BD_data$ph = gsub("ph: ", "", BD_data$ph)
BD_data$ph = gsub('"', "", BD_data$ph)

BD_data$pmi = gsub("pmi: ", "", BD_data$pmi)
BD_data$pmi = gsub('"', "", BD_data$pmi)

BD_data$suicide = gsub("suicide: ", "", BD_data$suicide)
BD_data$suicide = gsub('"', "", BD_data$suicide)

BD_data$death = gsub("type of death: ", "", BD_data$death)
BD_data$death = gsub('"', "", BD_data$death)

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$age = as.numeric(BD_data$age)
BD_data$sex = factor(BD_data$sex)
BD_data$ph = as.numeric(BD_data$ph)
BD_data$pmi = as.numeric(BD_data$pmi)
BD_data$suicide = factor(BD_data$suicide)
BD_data$death = factor(BD_data$death)

#get the raw data
files = list.files("GSE208338_RAW/")
paths = paste0("GSE208338_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE208338_RAW/")
cel_paths = paste0("GSE208338_RAW/", files)

BD_data_ids = gsub('\"', "", BD_data_ids)
files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE208338_RAW/", files_filtered)

rawdata_GSE208338 = read.celfiles(filenames = cel_paths_filtered)

#check quality
boxplot(rawdata_GSE208338, target = "core",  main = "Raw Intensity Distribution")

eset_GSE208338 = rma(rawdata_GSE208338)
eset_GSE208338_expression = exprs(eset_GSE208338)

#clean colnames
colnames(eset_GSE208338_expression) = sapply(colnames(eset_GSE208338_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

all(colnames(eset_GSE208338_expression)==rownames(BD_data$geo_accession)) #TRUE

#get the annotation file

annotation_data = select(huex10sttranscriptcluster.db, 
                         keys = rownames(eset_GSE208338_expression), 
                         columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                         keytype = "PROBEID")
##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+age+sex+ph+pmi+suicide+death, data = BD_data)
all(rownames(design)==rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE208338_expression, design)
fit = eBayes(fit)
GSE208338_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE208338_topTable_WITH_covars$ID = rownames(GSE208338_topTable_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 10)]
colnames(SE) = c("SE", "ID")
GSE208338_topTable_WITH_covars = inner_join(GSE208338_topTable_WITH_covars, SE, by="ID")

#annotate the table
GSE208338_topTable_WITH_covars = inner_join(GSE208338_topTable_WITH_covars, annotation_data[, c("SYMBOL", "PROBEID")], by=c("ID" = "PROBEID"))
GSE208338_topTable_WITH_covars = GSE208338_topTable_WITH_covars[, c(9,11,1:8, 10)]

#remove rows with missing gene symbol
GSE208338_topTable_WITH_covars_filtered = GSE208338_topTable_WITH_covars[!is.na(GSE208338_topTable_WITH_covars$SYMBOL), ]

##DE NO COVARS
#create design matrix
design0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(design0)==rownames(BD_data)) #TRUE

#fit the model
fit0 = lmFit(eset_GSE208338_expression, design0)
fit0 = eBayes(fit0)
GSE208338_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number=Inf, confint = T)
GSE208338_topTable_NO_covars$ID = rownames(GSE208338_topTable_NO_covars)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "ID")
GSE208338_topTable_NO_covars = inner_join(GSE208338_topTable_NO_covars, SE0, by="ID")

#annotate the table
GSE208338_topTable_NO_covars = inner_join(GSE208338_topTable_NO_covars, annotation_data[, c("SYMBOL", "PROBEID")], by=c("ID"="PROBEID"))
GSE208338_topTable_NO_covars = GSE208338_topTable_NO_covars[, c(9,11, 1:8, 10)]

#remove rows with missing gene symbol
GSE208338_topTable_NO_covars_filtered = GSE208338_topTable_NO_covars[!is.na(GSE208338_topTable_NO_covars$SYMBOL),]

#save the results
write.csv(GSE208338_topTable_NO_covars_filtered, "GSE208338_table_NO_covars.csv")
write.csv(GSE208338_topTable_WITH_covars_filtered, "GSE208338_table_WITH_covars.csv")

#clean the tables
colnames(GSE208338_table_NO_covars)[3] = "Gene Symbol"
colnames(GSE208338_table_WITH_covars)[3] = "Gene Symbol"

#no covars 
#delete rows with more than one gene mapped
GSE208338_table_NO_covars_cleaned = GSE208338_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE208338_table_NO_covars_cleaned = GSE208338_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE208338_table_NO_covars_cleaned = GSE208338_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE208338_table_NO_covars_cleaned$Study = "GSE208338"
write_csv(GSE208338_table_NO_covars_cleaned, "GSE208338_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE208338_table_WITH_covars_cleaned = GSE208338_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE208338_table_WITH_covars_cleaned = GSE208338_table_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE208338_table_WITH_covars_cleaned = GSE208338_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE208338_table_WITH_covars_cleaned$Study = "GSE208338"
write_csv(GSE208338_table_WITH_covars_cleaned, "GSE208338_with_covars_cleaned.csv")













