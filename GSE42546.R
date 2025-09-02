rm(list=ls())
library(stringi)
library(dplyr)

##GSE42546 AB SOLiD 4 System
geo_GSE42546 = GEOquery::getGEO("GSE42546")
geo_GSE42546 = geo_GSE42546[[1]]
geo_GSE42546_pheno = geo_GSE42546@phenoData@data

#clean diagnosis column in phenodata
colnames(geo_GSE42546_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE42546_pheno),
                                                          pattern = 'diagnosis:ch1',
                                                          replacement = 'diagnosis')
#filter bipolar disorder and controls
BD_data = geo_GSE42546_pheno %>%
  filter(grepl("control", diagnosis, ignore.case = T)|
           grepl("bipolar", diagnosis, ignore.case = T))
BD_data$diagnosis = gsub("bipolar", "bipolar disorder", BD_data$diagnosis)
BD_data_ids = BD_data$title

#clean colnames in BD_data
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                                      pattern = 'age at death:ch1',
                                                      replacement = 'age')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                                      pattern = 'gender:ch1',
                                                      replacement = 'gender')
#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$gender = factor(BD_data$gender)
BD_data$age = as.numeric(BD_data$age)


#get the annotation file
columns = c('GeneID', 'Symbol', 'chromosome', 'GeneType')
annotation_data = read.columns(file = 'Human.GRCh38.p13.annot.tsv', 
                               required.col = columns, 
                               sep = '\t', 
                               stringsAsFactors = F)
#create count matrix
count_matrix = read.delim(file = "GSE42546_RawCounts.txt", check.names = FALSE)
count_matrix$gene = sapply(count_matrix$gene, function(x){
  x = unlist(stringi::stri_split_fixed(x, pattern=";"))
  x = x[2]
  return(x)
})

cleaned_count_matrix = read.delim(file = "GSE42546_CleanedRawCounts.txt", check.names = FALSE)
rownames(cleaned_count_matrix) = cleaned_count_matrix$GeneSymbol

#filter counts matrix
count_matrix_filtered = cleaned_count_matrix[, colnames(cleaned_count_matrix) %in% BD_data_ids]
count_matrix_filtered$GeneID = rownames(cleaned_count_matrix)

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis + age + gender, data = BD_data)
rownames(design) = BD_data$title
all(rownames(design) == BD_data$title) #TRUE

#create DGEList
dge = DGEList(counts = count_matrix_filtered, annotation.columns = "GeneID")
keep = filterByExpr(dge, design) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
all(rownames(dge$samples) == rownames(design)) #FALSE
#rearrange rows in design
design = design[match(rownames(dge$samples), rownames(design)), ]
all(rownames(dge$samples) == rownames(design)) #TRUE


v = voom(dge, design, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design)
fit = eBayes(fit)
GSE42546_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 5)]
colnames(SE) = c("SE", "GeneID")
GSE42546_topTable_WITH_covars = inner_join(GSE42546_topTable_WITH_covars, SE, by="GeneID")

#annotate the table
GSE42546_topTable_WITH_covars = inner_join(GSE42546_topTable_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by=c("GeneID"="Symbol"))
GSE42546_topTable_WITH_covars = GSE42546_topTable_WITH_covars[, c(11, 1, 2:10)]
any(GSE42546_topTable_WITH_covars$GeneID =="") #FALSE

colnames(GSE42546_topTable_WITH_covars)[2] = "Gene Symbol"
colnames(GSE42546_topTable_WITH_covars)[1] = "Gene ID"


##DE NO COVARS
#create design matrix
design0 = model.matrix(~diagnosis, data=BD_data)
all(rownames(design0) == rownames(BD_data)) #TRUE
rownames(design0) = BD_data$title

#create DGEList
dge0 = DGEList(counts = count_matrix_filtered, annotation.columns = "GeneID")
keep = filterByExpr(dge0, design0) #remove very low counts
dge0 = dge0[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge0 = calcNormFactors(dge0) #performs TTM normalization
all(rownames(dge0$samples) == rownames(design0)) #FALSE
#rearrange rows in design
design0 = design0[match(rownames(dge0$samples), rownames(design0)), ]
all(rownames(dge0$samples) == rownames(design0)) #TRUE

v0 = voom(dge0, design0, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v0$targets) == rownames(v0$design)) #TRUE

#fit the model
fit0 = lmFit(v0, design0)
fit0 = eBayes(fit0)
GSE42546_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE42546_topTable_NO_covars = inner_join(GSE42546_topTable_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE42546_topTable_NO_covars = inner_join(GSE42546_topTable_NO_covars, annotation_data[, c("GeneID", "Symbol")], by=c("GeneID"="Symbol"))
GSE42546_topTable_NO_covars = GSE42546_topTable_NO_covars[, c(1, 11, 2:10)]
GSE42546_topTable_NO_covars = GSE42546_topTable_NO_covars[, c(2,1,3:11)]
colnames(GSE42546_topTable_NO_covars)[1] = "Gene ID"
colnames(GSE42546_topTable_NO_covars)[2] = "Gene Symbol"
any(GSE42546_topTable_NO_covars$`Gene Symbol`=="") #FALSE

#save the results
write.csv(GSE42546_topTable_NO_covars, "GSE42546_table_NO_covars.csv")
write.csv(GSE42546_topTable_WITH_covars, "GSE42546_table_WITH_covars.csv")

#clean the tables
#no covars 
#delete rows with more than one gene mapped
GSE42546_table_NO_covars_cleaned = GSE42546_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE42546_table_NO_covars_cleaned = GSE42546_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE42546_table_NO_covars_cleaned = GSE42546_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE42546_table_NO_covars_cleaned$Study = "GSE42546"
write_csv(GSE42546_table_NO_covars_cleaned, "GSE42546_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE42546_table_WITH_covars_cleaned = GSE42546_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE42546_table_WITH_covars_cleaned = GSE42546_table_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE42546_table_WITH_covars_cleaned = GSE42546_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE42546_table_WITH_covars_cleaned$Study = "GSE42546"
write_csv(GSE42546_table_WITH_covars_cleaned, "GSE42546_with_covars_cleaned.csv")









