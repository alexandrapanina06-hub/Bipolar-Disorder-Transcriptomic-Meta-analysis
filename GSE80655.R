rm(list=ls())
## GSE80655 Illumina HiSeq 2000
library(limma)
library(edgeR)
#get the file 
geo_GSE80655 = GEOquery::getGEO("GSE80655")
geo_GSE80655 = geo_GSE80655[[1]]
geo_GSE80655_pheno = geo_GSE80655@phenoData@data

#get count matrix
count_matrix = read.table(file = "GSE80655_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = '\t')
rownames(count_matrix) = count_matrix$GeneID

#add diagnosis column
geo_GSE80655_pheno = geo_GSE80655_pheno %>%
  mutate(diagnosis = case_when(
    grepl("Control", `clinical diagnosis:ch1`) ~ "control",
    grepl("Bipolar Disorder", `clinical diagnosis:ch1`) ~"bipolar disorder",
    grepl("Schizophrenia", `clinical diagnosis:ch1`) ~ "schizophrenia",
    grepl("Major Depression", `clinical diagnosis:ch1`) ~ "major depression",
    TRUE ~ "Other"
  ))

#filter bipolar disorder and controls
BD_data = geo_GSE80655_pheno %>%
  filter(grepl("control", diagnosis, ignore.case = T)|
           grepl("bipolar disorder", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)


#clean colnames in BD_data
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'age at death:ch1',
                                           replacement = 'age_death')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'brain ph:ch1',
                                           replacement = 'ph')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'brain region:ch1',
                                           replacement = 'brain_region')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'ethnicity:ch1',
                                           replacement = 'ethnicity')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'gender:ch1',
                                           replacement = 'gender')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'post-mortem interval:ch1',
                                           replacement = 'pmi')
#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$age_death = as.numeric(BD_data$age_death)
BD_data$ph = as.numeric(BD_data$ph)
BD_data$brain_region = factor(BD_data$brain_region)
BD_data$ethnicity = factor(BD_data$ethnicity)
BD_data$gender = factor(BD_data$gender)
BD_data$pmi = as.numeric(BD_data$pmi)



#get annotation file
columns = c('GeneID', 'Symbol', 'chromosome', 'GeneType')
annotation_data = read.columns(file = 'Human.GRCh38.p13.annot.tsv', 
                               required.col = columns, 
                               sep = '\t', 
                               stringsAsFactors = F)
##DE WITH COVARS

#create design matrix
design = model.matrix(~diagnosis + age_death+ph+brain_region+ethnicity+gender+pmi, data=BD_data)
all(rownames(design) == rownames(BD_data)) #TRUE

#filter countmatrix 
count_matrix_filtered = count_matrix[, colnames(count_matrix) %in% BD_data_ids]
count_matrix_filtered$GeneID = rownames(count_matrix_filtered)

#create DGElist
dge = DGEList(counts = count_matrix_filtered, annotation.columns = "GeneID")
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
GSE80655_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 13)]
colnames(SE) = c("SE", "GeneID")
GSE80655_topTable_WITH_covars = inner_join(GSE80655_topTable_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE80655_topTable_WITH_covars = inner_join(GSE80655_topTable_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_WITH_covars = GSE80655_topTable_WITH_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_WITH_covars$Symbol)) #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(design0)==rownames(BD_data)) #TRUE


#create DGElist
dge = DGEList(counts = count_matrix_filtered, annotation.columns = "GeneID")
keep0 = filterByExpr(dge, design0)
dge = dge[keep0, , keep.lib.sizes = F]

#scale normalization
dge = calcNormFactors(dge)
all(rownames(dge$samples) == rownames(design0)) #TRUE

v = voom(dge, design0, plot = T)
all(rownames(v$targets) == rownames(v$design0)) #TRUE

#fit the model
fit0 = lmFit(v, design0)
fit0 = eBayes(fit0)
GSE80655_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE80655_topTable_NO_covars = inner_join(GSE80655_topTable_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE80655_topTable_NO_covars = inner_join(GSE80655_topTable_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_NO_covars = GSE80655_topTable_NO_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_NO_covars$Symbol)) #FALSE

#save the results
write.csv(GSE80655_topTable_NO_covars, "GSE80655_table_NO_covars.csv")
write.csv(GSE80655_topTable_WITH_covars, "GSE80655_table_WITH_covars.csv")




###the dataset includes 3 different tissues from each patient
### anterior cingulate cortex, dorsolateral prefrontal cortex, and nucleus accumbens

#anterior cingulate cortex (acc)
#filter acc data
BD_data_acc = BD_data %>%
  filter(grepl("AnCg", brain_region, ignore.case = T))
BD_data_acc_ids = rownames(BD_data_acc)

#filter counts matrix
counts_acc = count_matrix[, colnames(count_matrix) %in% rownames(BD_data_acc)]
counts_acc$GeneID = rownames(counts_acc)

##DE WITH COVARS 
#create design matrix
design_acc = model.matrix(~diagnosis+age_death+ph+ethnicity+gender+pmi, data=BD_data_acc)
all(rownames(design_acc) == rownames(BD_data_acc)) #TRUE

##create DGElist
dge = DGEList(counts = counts_acc, annotation.columns = "GeneID")
keep = filterByExpr(dge, design_acc) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
all(rownames(dge$samples) == rownames(design_acc)) #TRUE

v = voom(dge, design_acc, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design_acc)
fit = eBayes(fit)
GSE80655_topTable_acc_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 11)]
colnames(SE) = c("SE", "GeneID")
GSE80655_topTable_acc_WITH_covars = inner_join(GSE80655_topTable_acc_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE80655_topTable_acc_WITH_covars = inner_join(GSE80655_topTable_acc_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_acc_WITH_covars = GSE80655_topTable_acc_WITH_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_acc_WITH_covars$Symbol)) #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data = BD_data_acc)
all(rownames(design0)==rownames(BD_data_acc)) #TRUE


#create DGElist
dge0 = DGEList(counts = counts_acc, annotation.columns = "GeneID")
keep0 = filterByExpr(dge0, design0)
dge0 = dge0[keep0, , keep.lib.sizes = F]

#scale normalization
dge0 = calcNormFactors(dge0)
all(rownames(dge0$samples) == rownames(design0)) #TRUE

v0 = voom(dge0, design0, plot = T)
all(rownames(v0$targets) == rownames(v0$design0)) #TRUE

#fit the model
fit0 = lmFit(v0, design0)
fit0 = eBayes(fit0)
GSE80655_topTable_acc_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE80655_topTable_acc_NO_covars = inner_join(GSE80655_topTable_acc_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE80655_topTable_acc_NO_covars = inner_join(GSE80655_topTable_acc_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_acc_NO_covars = GSE80655_topTable_acc_NO_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_acc_NO_covars$Symbol)) #FALSE

#save the results
write.csv(GSE80655_topTable_acc_NO_covars, "GSE80655_table_acc_NO_covars.csv")
write.csv(GSE80655_topTable_acc_WITH_covars, "GSE80655_table_acc_WITH_covars.csv")


#dorsolateral prefrontal cortex (pfc)
#filter pfc data
BD_data_pfc = BD_data %>%
  filter(grepl("DLPFC", brain_region, ignore.case = T))
BD_data_pfc_ids = rownames(BD_data_pfc)


#filter counts matrix
counts_pfc = count_matrix[, colnames(count_matrix) %in% rownames(BD_data_pfc)]
counts_pfc$GeneID = rownames(counts_pfc)


##DE WITH COVARS 
#create design matrix
design_pfc = model.matrix(~diagnosis+age_death+ph+ethnicity+gender+pmi, data=BD_data_pfc)
all(rownames(design_pfc) == rownames(BD_data_pfc)) #TRUE

##create DGElist
dge = DGEList(counts = counts_pfc, annotation.columns = "GeneID")
keep = filterByExpr(dge, design_pfc) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
all(rownames(dge$samples) == rownames(design_pfc)) #TRUE

v = voom(dge, design_pfc, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design_pfc)
fit = eBayes(fit)
GSE80655_topTable_pfc_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 11)]
colnames(SE) = c("SE", "GeneID")
GSE80655_topTable_pfc_WITH_covars = inner_join(GSE80655_topTable_pfc_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE80655_topTable_pfc_WITH_covars = inner_join(GSE80655_topTable_pfc_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_pfc_WITH_covars = GSE80655_topTable_pfc_WITH_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_pfc_WITH_covars$Symbol)) #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data = BD_data_pfc)
all(rownames(design0)==rownames(BD_data_pfc)) #TRUE


#create DGElist
dge0 = DGEList(counts = counts_pfc, annotation.columns = "GeneID")
keep0 = filterByExpr(dge0, design0)
dge0 = dge0[keep0, , keep.lib.sizes = F]

#scale normalization
dge0 = calcNormFactors(dge0)
all(rownames(dge0$samples) == rownames(design0)) #TRUE

v0 = voom(dge0, design0, plot = T)
all(rownames(v0$targets) == rownames(v0$design0)) #TRUE

#fit the model
fit0 = lmFit(v0, design0)
fit0 = eBayes(fit0)
GSE80655_topTable_pfc_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE80655_topTable_pfc_NO_covars = inner_join(GSE80655_topTable_pfc_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE80655_topTable_pfc_NO_covars = inner_join(GSE80655_topTable_pfc_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_pfc_NO_covars = GSE80655_topTable_pfc_NO_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_pfc_NO_covars$Symbol)) #FALSE

#save the results
write.csv(GSE80655_topTable_pfc_NO_covars, "GSE80655_table_pfc_NO_covars.csv")
write.csv(GSE80655_topTable_pfc_WITH_covars, "GSE80655_table_pfc_WITH_covars.csv")


#nucleus accumbens (nac)
#filter data for nac
BD_data_nac = BD_data %>% 
  filter(grepl("nAcc", brain_region, ignore.case = T))
BD_data_nac_ids = rownames(BD_data_nac)

#filter counts matrix
counts_nac = count_matrix[, colnames(count_matrix) %in% BD_data_nac_ids]
counts_nac$GeneID = rownames(counts_nac)

##DE WITH COVARS 
#create design matrix
design_nac = model.matrix(~diagnosis+age_death+ph+ethnicity+gender+pmi, data=BD_data_nac)
all(rownames(design_nac) == rownames(BD_data_nac)) #TRUE

##create DGElist
dge = DGEList(counts = counts_nac, annotation.columns = "GeneID")
keep = filterByExpr(dge, design_nac) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
all(rownames(dge$samples) == rownames(design_nac)) #TRUE

v = voom(dge, design_nac, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design_nac)
fit = eBayes(fit)
GSE80655_topTable_nac_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2, 11)]
colnames(SE) = c("SE", "GeneID")
GSE80655_topTable_nac_WITH_covars = inner_join(GSE80655_topTable_nac_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE80655_topTable_nac_WITH_covars = inner_join(GSE80655_topTable_nac_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_nac_WITH_covars = GSE80655_topTable_nac_WITH_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_nac_WITH_covars$Symbol)) #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data = BD_data_nac)
all(rownames(design0)==rownames(BD_data_nac)) #TRUE


#create DGElist
dge0 = DGEList(counts = counts_nac, annotation.columns = "GeneID")
keep0 = filterByExpr(dge0, design0)
dge0 = dge0[keep0, , keep.lib.sizes = F]

#scale normalization
dge0 = calcNormFactors(dge0)
all(rownames(dge0$samples) == rownames(design0)) #TRUE

v0 = voom(dge0, design0, plot = T)
all(rownames(v0$targets) == rownames(v0$design0)) #TRUE

#fit the model
fit0 = lmFit(v0, design0)
fit0 = eBayes(fit0)
GSE80655_topTable_nac_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE80655_topTable_nac_NO_covars = inner_join(GSE80655_topTable_nac_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE80655_topTable_nac_NO_covars = inner_join(GSE80655_topTable_nac_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE80655_topTable_nac_NO_covars = GSE80655_topTable_nac_NO_covars[, c(1,11, 2:10)]

any(is.na(GSE80655_topTable_nac_NO_covars$Symbol)) #FALSE

#save the results
write.csv(GSE80655_topTable_nac_NO_covars, "GSE80655_table_nac_NO_covars.csv")
write.csv(GSE80655_topTable_nac_WITH_covars, "GSE80655_table_nac_WITH_covars.csv")

#clean the tables

colnames(GSE80655_table_pfc_NO_covars)[3] = "Gene Symbol"
colnames(GSE80655_table_pfc_WITH_covars)[3] = "Gene Symbol"


#no covars 
#delete rows with more than one gene mapped
GSE80655_table_pfc_NO_covars_cleaned = GSE80655_table_pfc_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE80655_table_pfc_NO_covars_cleaned = GSE80655_table_pfc_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE80655_table_pfc_NO_covars_cleaned = GSE80655_table_pfc_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE80655_table_pfc_NO_covars_cleaned$Study = "GSE80655"
write_csv(GSE80655_table_pfc_NO_covars_cleaned, "GSE80655_pfc_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE80655_table_pfc_WITH_covars_cleaned = GSE80655_table_pfc_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE80655_table_pfc_WITH_covars_cleaned = GSE80655_table_pfc_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE80655_table_pfc_WITH_covars_cleaned = GSE80655_table_pfc_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE80655_table_pfc_WITH_covars_cleaned$Study = "GSE80655"
write_csv(GSE80655_table_pfc_WITH_covars_cleaned, "GSE80655_pfc_with_covars_cleaned.csv")










