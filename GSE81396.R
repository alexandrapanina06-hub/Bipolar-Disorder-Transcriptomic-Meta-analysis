##GSE81396 Illumina Hiseq 2000
rm(list = ls())

geo_GSE81396 = GEOquery::getGEO("GSE81396")
geo_GSE81396 = geo_GSE81396[[1]]
geo_GSE81396_pheno = geo_GSE81396@phenoData@data

#create diagnosis column
geo_GSE81396_pheno = geo_GSE81396_pheno %>%
  mutate(diagnosis = case_when(
    grepl("healthy control", `subject status:ch1`) ~ "control",
    grepl("Bipolar Disorder individual", `subject status:ch1`) ~"bipolar disorder",
    TRUE ~ "other"
  ))
#clean colnames
colnames(geo_GSE81396_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE81396_pheno),
                                                      pattern = 'tissue subtype:ch1',
                                                      replacement = 'tissue')
#variable transformation
geo_GSE81396_pheno$diagnosis = factor(geo_GSE81396_pheno$diagnosis, levels = c("control", "bipolar disorder"))
geo_GSE81396_pheno$tissue = factor(geo_GSE81396_pheno$tissue)

#get the counts matrix
count_matrix = read.table(file = "GSE81396_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = '\t')
rownames(count_matrix) = count_matrix$GeneID

#get the annotation file
columns = c('GeneID', 'Symbol', 'chromosome', 'GeneType')
annotation_data = read.columns(file = 'Human.GRCh38.p13.annot.tsv', 
                               required.col = columns, 
                               sep = '\t', 
                               stringsAsFactors = F)
##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+tissue, data=geo_GSE81396_pheno)
all(rownames(design)==rownames(geo_GSE81396_pheno)) #TRUE

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
GSE81396_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[,c(2, 4)]
colnames(SE) = c("SE", "GeneID")
GSE81396_topTable_WITH_covars$GeneID = factor(GSE81396_topTable_WITH_covars$GeneID)
GSE81396_topTable_WITH_covars = inner_join(GSE81396_topTable_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
colnames(annotation_data)[2] = "Gene Symbol"
GSE81396_topTable_WITH_covars = inner_join(GSE81396_topTable_WITH_covars, annotation_data[, c("GeneID", "Gene Symbol")], by="GeneID")
GSE81396_topTable_WITH_covars = GSE81396_topTable_WITH_covars[, c(1,11, 2:10)]
any(GSE81396_topTable_WITH_covars$`Gene Symbol`=="") #FALSE

#DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data = geo_GSE81396_pheno)
all(rownames(design0) == rownames(geo_GSE81396_pheno)) #TRUE

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
GSE81396_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE81396_topTable_NO_covars$GeneID = factor(GSE81396_topTable_NO_covars$GeneID)
GSE81396_topTable_NO_covars = inner_join(GSE81396_topTable_NO_covars, SE0, by="GeneID")

#annotate the table
GSE81396_topTable_NO_covars = inner_join(GSE81396_topTable_NO_covars, annotation_data[, c("GeneID", "Gene Symbol")], by="GeneID")
GSE81396_topTable_NO_covars = GSE81396_topTable_NO_covars[, c(1,11, 2:10)]
any(GSE81396_topTable_NO_covars$`Gene Symbol`=="") #FALSE

#save the results
write.csv(GSE81396_topTable_NO_covars, "GSE81396_table_NO_covars.csv")
write.csv(GSE81396_topTable_WITH_covars, "GSE81396_table_WITH_covars.csv")


###dataset includes 2 parts of dorsal striatum: putamen and caudate nucleus
##putamen
#filter putamen samples
putamen = geo_GSE81396_pheno %>%
  filter(grepl("putamen", tissue, ignore.case = T))
putamen_ids = rownames(putamen)

#filter counts matrix
counts_putamen = count_matrix[, colnames(count_matrix) %in% putamen_ids]
counts_putamen$GeneID = rownames(counts_putamen)

##DE NO COVARS
#create design matrix
design = model.matrix(~diagnosis, data=putamen)
all(rownames(design)==rownames(putamen)) #TRUE

#create DGEList
dge = DGEList(counts = counts_putamen, annotation.columns = "GeneID")
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
GSE81396_topTable_putamen_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[,c(2, 3)]
colnames(SE) = c("SE", "GeneID")
GSE81396_topTable_putamen_NO_covars$GeneID = factor(GSE81396_topTable_putamen_NO_covars$GeneID)
GSE81396_topTable_putamen_NO_covars = inner_join(GSE81396_topTable_putamen_NO_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
colnames(annotation_data)[2] = "Gene Symbol"
GSE81396_topTable_putamen_NO_covars = inner_join(GSE81396_topTable_putamen_NO_covars, annotation_data[, c("GeneID", "Gene Symbol")], by="GeneID")
GSE81396_topTable_putamen_NO_covars = GSE81396_topTable_putamen_NO_covars[, c(1,11, 2:10)]
any(GSE81396_topTable_putamen_NO_covars$`Gene Symbol`=="") #FALSE

#save the results 
write.csv(GSE81396_topTable_putamen_NO_covars, "GSE81396_table_putamen_NO_covars.csv")

##caudate nucleus
#filter cuadate nucleus samples
caudate = geo_GSE81396_pheno %>%
  filter(grepl("cuadate nucleus", tissue, ignore.case = T))
caudate_ids = rownames(caudate)

#filter counts matrix
counts_caudate = count_matrix[, colnames(count_matrix) %in% caudate_ids]
counts_caudate$GeneID = rownames(counts_caudate)

##DE NO COVARS
#create design matrix
design = model.matrix(~diagnosis, data=caudate)
all(rownames(design)==rownames(caudate)) #TRUE

#create DGEList
dge = DGEList(counts = counts_caudate, annotation.columns = "GeneID")
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
GSE81396_topTable_caudate_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[,c(2, 3)]
colnames(SE) = c("SE", "GeneID")
GSE81396_topTable_caudate_NO_covars$GeneID = factor(GSE81396_topTable_caudate_NO_covars$GeneID)
GSE81396_topTable_caudate_NO_covars = inner_join(GSE81396_topTable_caudate_NO_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
colnames(annotation_data)[2] = "Gene Symbol"
GSE81396_topTable_caudate_NO_covars = inner_join(GSE81396_topTable_caudate_NO_covars, annotation_data[, c("GeneID", "Gene Symbol")], by="GeneID")
GSE81396_topTable_caudate_NO_covars = GSE81396_topTable_caudate_NO_covars[, c(1,11, 2:10)]
any(GSE81396_topTable_caudate_NO_covars$`Gene Symbol`=="") #FALSE

#save the results 
write_csv(GSE81396_topTable_caudate_NO_covars, "GSE81396_table_caudate_NO_covars.csv")


#clean the tables
##caudate
#no covars 
#delete rows with more than one gene mapped
GSE81396_table_caudate_NO_covars_cleaned = GSE81396_table_caudate_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE81396_table_caudate_NO_covars_cleaned = GSE81396_table_caudate_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))


#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE81396_table_caudate_NO_covars_cleaned = GSE81396_table_caudate_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE81396_table_caudate_NO_covars_cleaned$Study = "GSE81396"
write_csv(GSE81396_table_caudate_NO_covars_cleaned, "GSE81396_caudate_no_covars_cleaned.csv")

##putamen
#no covars 
#delete rows with more than one gene mapped
GSE81396_table_putamen_NO_covars_cleaned = GSE81396_table_putamen_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE81396_table_putamen_NO_covars_cleaned = GSE81396_table_putamen_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE81396_table_putamen_NO_covars_cleaned = GSE81396_table_putamen_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE81396_table_putamen_NO_covars_cleaned$Study = "GSE81396"
write_csv(GSE81396_table_putamen_NO_covars_cleaned, "GSE81396_putamen_no_covars_cleaned.csv")













