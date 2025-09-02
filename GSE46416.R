rm(list = ls())
#Affymetrix Human Exon 1.0 ST Array
library(oligo)

geo_GSE46416 = getGEO("GSE46416")
geo_GSE46416 = geo_GSE46416[[1]]
pheno_GSE46416 = geo_GSE46416@phenoData@data

#the dataset contains data from the same patients but in different disease state
#so, let's filter only manic state BD and controls

BD_data = pheno_GSE46416 %>%
  filter(grepl("mania", `bd phase:ch1`, ignore.case = T)|
           grepl("control", `disease status:ch1`, ignore.case = T))
BD_data_ids = rownames(BD_data)

#clean the colnames 
colnames(BD_data)[37] = "diagnosis"


#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder (BD)"))

#get the raw files
files = list.files("GSE46416_RAW/")
paths = paste0("GSE46416_RAW/", files)
sapply(paths, gunzip)

files = list.files("GSE46416_RAW")
cel_paths = paste0("GSE46416_RAW/", files)

#filter only BD and control files
files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
cel_paths_filtered = paste0("GSE46416_RAW/", files_filtered)

#get the raw data using oligo package
rawdata_GSE46416 = read.celfiles(filenames = cel_paths_filtered)

#quality check
boxplot(rawdata_GSE46416, target = "core", main = "Raw Intensity Distribution")

#normalization
eset_GSE46416 = rma(rawdata_GSE46416)
eset_GSE46416_expression = exprs(eset_GSE46416)

#clean colnames in eset_GSE78246
colnames(eset_GSE46416_expression) = sapply(colnames(eset_GSE46416_expression), function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern="_"))
  x = x[1]
  return(x)
})

#alignment check
all(colnames(eset_GSE46416_expression) == rownames(BD_data$geo_accession)) #TRUE

#create annotation file using huex10sttranscriptcluster.db
probe_ids = rownames(eset_GSE46416_expression)
annotation_data = select(huex10sttranscriptcluster.db, 
                         keys = probe_ids, 
                         columns = c("SYMBOL", "ENTREZID", "GENENAME"), 
                         keytype = "PROBEID")

#DGE NO covars
#create design matrix 
design = model.matrix(~diagnosis, data = BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(eset_GSE46416_expression, design)
fit = eBayes(fit)
GSE46416_topTable_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE46416_topTable_NO_covars$ID = rownames(GSE46416_topTable_NO_covars)

#add SE to teh table
SE = as.data.frame(sqrt(fit$s2.post) * fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,3)]
colnames(SE) = c("SE", "ID")
GSE46416_topTable_NO_covars = inner_join(GSE46416_topTable_NO_covars, SE, by="ID")

#annotate the table
GSE46416_topTable_NO_covars = inner_join(GSE46416_topTable_NO_covars, annotation_data[, c("SYMBOL", "PROBEID")], by = c("ID" = "PROBEID"))
GSE46416_topTable_NO_covars = GSE46416_topTable_NO_covars[!is.na(GSE46416_topTable_NO_covars$SYMBOL),]
GSE46416_topTable_NO_covars = GSE46416_topTable_NO_covars[, c(9, 11, 1, 10, 2:8)]
colnames(GSE46416_topTable_NO_covars)[2] = "Gene Symbol"

#save the results
write.csv(GSE46416_topTable_NO_covars, "GSE46416_topTable_NO_covars.csv")

##clean the tables
GSE46416_topTable_NO_covars = GSE46416_topTable_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE46416_topTable_NO_covars = GSE46416_topTable_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
#add study column
GSE46416_topTable_NO_covars$Study = "GSE46416"

#save
write_csv(GSE46416_topTable_NO_covars, "GSE46416_NO_covars.csv")








