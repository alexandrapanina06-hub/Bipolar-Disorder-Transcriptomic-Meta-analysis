rm(list = ls())

##GSE80336 Illumina HiSeq 2000
geo_GSE80336 = GEOquery::getGEO("GSE80336")
geo_GSE80336 = geo_GSE80336[[1]]
geo_GSE80336_pheno = geo_GSE80336@phenoData@data

#create diagnosis column
geo_GSE80336_pheno = geo_GSE80336_pheno %>%
  mutate(diagnosis = case_when(
    grepl("bipolar", source_name_ch1) ~ "bipolar disorder",
    grepl("control", source_name_ch1) ~ "control", 
    TRUE ~ "other"
  ))
#clean colnames
colnames(geo_GSE80336_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE80336_pheno),
                                                      pattern = 'age (years):ch1',
                                                      replacement = 'age')
colnames(geo_GSE80336_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE80336_pheno),
                                                      pattern = 'postmortem interval (hours):ch1',
                                                      replacement = 'pmi')
colnames(geo_GSE80336_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE80336_pheno),
                                                      pattern = 'Sex:ch1',
                                                      replacement = 'sex')

#variable transformation
geo_GSE80336_pheno$diagnosis = factor(geo_GSE80336_pheno$diagnosis, levels = c("control", "bipolar disorder"))
geo_GSE80336_pheno$age = as.numeric(geo_GSE80336_pheno$age)
geo_GSE80336_pheno$sex = factor(geo_GSE80336_pheno$sex)
geo_GSE80336_pheno$pmi = as.numeric(geo_GSE80336_pheno$pmi)

#create counts matrix
count_matrix = read.table(file = "GSE80336_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = '\t')
rownames(count_matrix) = count_matrix$GeneID

#get annotation file
columns = c('GeneID', 'Symbol', 'chromosome', 'GeneType')
annotation_data = read.columns(file = 'Human.GRCh38.p13.annot.tsv', 
                               required.col = columns, 
                               sep = '\t', 
                               stringsAsFactors = F)
annotation_data$GeneID = factor(annotation_data$GeneID)
colnames(annotation_data)[2] = "Gene Symbol"

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+age+sex+pmi, data = geo_GSE80336_pheno)
all(rownames(design) == rownames(geo_GSE80336_pheno)) #TRUE

#create DGEList
dge = DGEList(counts = count_matrix, annotation.columns = "GeneID")
keep = filterByExpr(dge, design) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
all(rownames(dge$samples) == rownames(design)) #TRUE


v = voom(dge, design, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design)
fit = eBayes(fit)
GSE80336_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 6)]
colnames(SE) = c("SE", "GeneID")
GSE80336_topTable_WITH_covars$GeneID = factor(GSE80336_topTable_WITH_covars$GeneID)
GSE80336_topTable_WITH_covars = inner_join(GSE80336_topTable_WITH_covars, SE, by = "GeneID")

#annotate the table
GSE80336_topTable_WITH_covars = inner_join(GSE80336_topTable_WITH_covars, annotation_data[, c("Gene Symbol", "GeneID")], by = "GeneID")
GSE80336_topTable_WITH_covars = GSE80336_topTable_WITH_covars[, c(1,11,2:10)]
any(GSE80336_topTable_WITH_covars$`Gene Symbol`=="") #FALSE

##DE NO COVARS
#create design matrix
design0 = model.matrix(~diagnosis, data = geo_GSE80336_pheno)
all(rownames(design0) == rownames(geo_GSE80336_pheno)) #TRUE

#create DGEList
dge0 = DGEList(counts = count_matrix, annotation.columns = "GeneID")
keep0 = filterByExpr(dge0, design0) #remove very low counts
dge0 = dge0[keep0, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge0 = calcNormFactors(dge0) #performs TTM normalization
all(rownames(dge0$samples) == rownames(design0)) #TRUE


v0 = voom(dge0, design0, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v0$targets) == rownames(v0$design)) #TRUE

#fit the model
fit0 = lmFit(v0, design0)
fit0 = eBayes(fit0)
GSE80336_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE80336_topTable_NO_covars$GeneID = factor(GSE80336_topTable_NO_covars$GeneID)
GSE80336_topTable_NO_covars = inner_join(GSE80336_topTable_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE80336_topTable_NO_covars = inner_join(GSE80336_topTable_NO_covars, annotation_data[, c("Gene Symbol", "GeneID")], by="GeneID")
GSE80336_topTable_NO_covars = GSE80336_topTable_NO_covars[, c(1,11, 2:10)]
any(is.na(GSE80336_topTable_NO_covars$`Gene Symbol`)) #FALSE

#save the results
write.csv(GSE80336_topTable_NO_covars, "GSE80336_table_NO_covars.csv")
write.csv(GSE80336_topTable_WITH_covars, "GSE80336_table_WITH_covars.csv")

#clean the tables
#no covars 
#delete rows with more than one gene mapped
GSE80336_table_NO_covars_cleaned = GSE80336_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE80336_table_NO_covars_cleaned = GSE80336_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE80336_table_NO_covars_cleaned = GSE80336_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE80336_table_NO_covars_cleaned$Study = "GSE80336"
write_csv(GSE80336_table_NO_covars_cleaned, "GSE80336_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE80336_table_WITH_covars_cleaned = GSE80336_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE80336_table_WITH_covars_cleaned = GSE80336_table_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE80336_table_WITH_covars_cleaned = GSE80336_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE80336_table_WITH_covars_cleaned$Study = "GSE80336"
write_csv(GSE80336_table_WITH_covars_cleaned, "GSE80336_with_covars_cleaned.csv")



