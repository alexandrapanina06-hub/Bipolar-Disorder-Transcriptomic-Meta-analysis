##GSE78936 Illumina Hiseq 2000
rm(list=ls())


geo_GSE78936 = GEOquery::getGEO("GSE78936")
geo_GSE78936 = geo_GSE78936[[1]]
geo_GSE78936_pheno = geo_GSE78936@phenoData@data

#create diagnosis column
geo_GSE78936_pheno = geo_GSE78936_pheno %>%
  mutate(diagnosis = case_when(
    grepl("Bipolar Disorder", `disease stage:ch1`) ~"bipolar disorder",
    grepl("Control", `disease stage:ch1`) ~ "control",
    grepl("Schizophrenia", `disease stage:ch1`) ~"schizophrenia",
    TRUE ~ "other"
  ))

#filter bipolar disorder and controls
BD_data = geo_GSE78936_pheno %>%
  filter(grepl("bipolar disorder", diagnosis, ignore.case = T)|
           grepl("control", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)

#clean colnames
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = "brain region:ch1",
                                           replacement = "brain_region")

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$brain_region = factor(BD_data$brain_region)

#create counts matrix
count_matrix = read.table(file = "GSE78936_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = '\t')
rownames(count_matrix) = count_matrix$GeneID

#filter countmatrix 
count_matrix_filtered = count_matrix[, colnames(count_matrix) %in% BD_data_ids]
count_matrix_filtered$GeneID = rownames(count_matrix_filtered)

#get the annotation file
columns = c('GeneID', 'Symbol', 'chromosome', 'GeneType')
annotation_data = read.columns(file = 'Human.GRCh38.p13.annot.tsv', 
                               required.col = columns, 
                               sep = '\t', 
                               stringsAsFactors = F)

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+brain_region, data = BD_data)
all(rownames(design)==rownames(BD_data)) #TRUE

#create DGEList
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
GSE78936_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2,5)]
colnames(SE) = c("SE", "GeneID")
GSE78936_topTable_WITH_covars = inner_join(GSE78936_topTable_WITH_covars, SE, by="GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE78936_topTable_WITH_covars = inner_join(GSE78936_topTable_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE78936_topTable_WITH_covars = GSE78936_topTable_WITH_covars[, c(1,11, 2:10)]
colnames(GSE78936_topTable_WITH_covars)[2] = "Gene Symbol"
any(GSE78936_topTable_WITH_covars$`Gene Symbol`=="") #FALSE

##DE NO COVARS
#create design matrix
design0 = model.matrix(~diagnosis, data=BD_data)
all(rownames(design0) == rownames(BD_data)) #TRUE

#create DGEList
dge0 = DGEList(counts = count_matrix_filtered, annotation.columns = "GeneID")
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
GSE78936_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE78936_topTable_NO_covars = inner_join(GSE78936_topTable_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE78936_topTable_NO_covars = inner_join(GSE78936_topTable_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE78936_topTable_NO_covars = GSE78936_topTable_NO_covars[, c(1,11, 2:10)]
colnames(GSE78936_topTable_NO_covars)[2] = "Gene Symbol"
any(GSE78936_topTable_NO_covars$`Gene Symbol`=="") #FALSE

#save the results
write.csv(GSE78936_topTable_NO_covars, "GSE78936_table_NO_covars.csv")
write.csv(GSE78936_topTable_WITH_covars, "GSE78936_table_WITH_covars.csv")



####the dataset includes 3 tissues from diffreent brain areas from the same patients

##BA9
#filter BA9 data
BD_data_ba9 = BD_data %>%
  filter(grepl("Brodmann's Area 9", brain_region, ignore.case = T))
BD_data_ba9_ids = rownames(BD_data_ba9)

#filter counts matrix
counts_ba9 = count_matrix[, colnames(count_matrix) %in% BD_data_ba9_ids]
counts_ba9$GeneID = rownames(counts_ba9)

#variable transformation
BD_data_ba9$diagnosis = factor(BD_data_ba9$diagnosis, levels = c("control", "bipolar disorder"))

#DE NO COVARS
#create design matrix
design = model.matrix(~diagnosis, data = BD_data_ba9)
all(rownames(design)==rownames(BD_data_ba9)) #TRUE

#create DGEList
dge = DGEList(counts = counts_ba9, annotation.columns = "GeneID")
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
GSE78936_topTable_ba9_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "GeneID")
GSE78936_topTable_ba9_NO_covars = inner_join(GSE78936_topTable_ba9_NO_covars, SE, by="GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE78936_topTable_ba9_NO_covars = inner_join(GSE78936_topTable_ba9_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE78936_topTable_ba9_NO_covars = GSE78936_topTable_ba9_NO_covars[, c(1,11, 2:10)]
colnames(GSE78936_topTable_ba9_NO_covars)[2] = "Gene Symbol"
any(GSE78936_topTable_ba9_NO_covars$`Gene Symbol`=="") #FALSE

#save the results
write.csv(GSE78936_topTable_ba9_NO_covars, "GSE78936_table_ba9_NO_covars.csv")

##BA11
#filter BA11 data
BD_data_ba11 = BD_data %>%
  filter(grepl("Brodmann's Area 11", brain_region, ignore.case = T))
BD_data_ba11_ids = rownames(BD_data_ba11)

#filter counts matrix
counts_ba11 = count_matrix[, colnames(count_matrix) %in% BD_data_ba11_ids]
counts_ba11$GeneID = rownames(counts_ba11)

#variable transformation
BD_data_ba11$diagnosis = factor(BD_data_ba11$diagnosis, levels = c("control", "bipolar disorder"))

#DE NO COVARS
#create design matrix
design = model.matrix(~diagnosis, data = BD_data_ba11)
all(rownames(design)==rownames(BD_data_ba11)) #TRUE

#create DGEList
dge = DGEList(counts = counts_ba11, annotation.columns = "GeneID")
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
GSE78936_topTable_ba11_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "GeneID")
GSE78936_topTable_ba11_NO_covars = inner_join(GSE78936_topTable_ba11_NO_covars, SE, by="GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE78936_topTable_ba11_NO_covars = inner_join(GSE78936_topTable_ba11_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE78936_topTable_ba11_NO_covars = GSE78936_topTable_ba11_NO_covars[, c(1,11, 2:10)]
colnames(GSE78936_topTable_ba11_NO_covars)[2] = "Gene Symbol"
any(GSE78936_topTable_ba11_NO_covars$`Gene Symbol`=="") #FALSE

#save the results
write.csv(GSE78936_topTable_ba11_NO_covars, "GSE78936_table_ba11_NO_covars.csv")

##BA24
#filter BA24 data
BD_data_ba24 = BD_data %>%
  filter(grepl("Brodmann's Area 24", brain_region, ignore.case = T))
BD_data_ba24_ids = rownames(BD_data_ba24)

#filter counts matrix
counts_ba24 = count_matrix[, colnames(count_matrix) %in% BD_data_ba24_ids]
counts_ba24$GeneID = rownames(counts_ba24)

#variable transformation
BD_data_ba24$diagnosis = factor(BD_data_ba24$diagnosis, levels = c("control", "bipolar disorder"))

#DE NO COVARS
#create design matrix
design = model.matrix(~diagnosis, data = BD_data_ba24)
all(rownames(design)==rownames(BD_data_ba24)) #TRUE

#create DGEList
dge = DGEList(counts = counts_ba24, annotation.columns = "GeneID")
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
GSE78936_topTable_ba24_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, 2:3]
colnames(SE) = c("SE", "GeneID")
GSE78936_topTable_ba24_NO_covars = inner_join(GSE78936_topTable_ba24_NO_covars, SE, by="GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE78936_topTable_ba24_NO_covars = inner_join(GSE78936_topTable_ba24_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE78936_topTable_ba24_NO_covars = GSE78936_topTable_ba24_NO_covars[, c(1,11, 2:10)]
colnames(GSE78936_topTable_ba24_NO_covars)[2] = "Gene Symbol"
any(GSE78936_topTable_ba24_NO_covars$`Gene Symbol`=="") #FALSE

#save the results
write.csv(GSE78936_topTable_ba24_NO_covars, "GSE78936_table_ba24_NO_covars.csv")

#clean the tables
##BA9 
#no covars
#delete rows with more than one gene mapped
GSE78936_table_ba9_NO_covars_cleaned = GSE78936_table_ba9_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE78936_table_ba9_NO_covars_cleaned = GSE78936_table_ba9_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE78936_table_ba9_NO_covars_cleaned = GSE78936_table_ba9_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE78936_table_ba9_NO_covars_cleaned$Study = "GSE78936"
write_csv(GSE78936_table_ba9_NO_covars_cleaned, "GSE78936_ba9_no_covars_cleaned.csv")


##BA11 
#no covars
#delete rows with more than one gene mapped
GSE78936_table_ba11_NO_covars_cleaned = GSE78936_table_ba11_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE78936_table_ba11_NO_covars_cleaned = GSE78936_table_ba11_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE78936_table_ba11_NO_covars_cleaned = GSE78936_table_ba11_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE78936_table_ba11_NO_covars_cleaned$Study = "GSE78936"
write_csv(GSE78936_table_ba11_NO_covars_cleaned, "GSE78936_ba11_no_covars_cleaned.csv")



##BA24
#no covars
#delete rows with more than one gene mapped
GSE78936_table_ba24_NO_covars_cleaned = GSE78936_table_ba24_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE78936_table_ba24_NO_covars_cleaned = GSE78936_table_ba24_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE78936_table_ba24_NO_covars_cleaned = GSE78936_table_ba24_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE78936_table_ba24_NO_covars_cleaned$Study = "GSE78936"
write_csv(GSE78936_table_ba24_NO_covars_cleaned, "GSE78936_ba24_no_covars_cleaned.csv")





















