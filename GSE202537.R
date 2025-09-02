rm(list=ls())

##GSE202537 Illumina NextSeq 500

geo_GSE202537 = GEOquery::getGEO("GSE202537")
geo_GSE202537 = geo_GSE202537[[1]]
geo_GSE202537_pheno = geo_GSE202537@phenoData@data

#create diagnosis column
geo_GSE202537_pheno = geo_GSE202537_pheno %>%
  mutate(diagnosis = case_when(
    grepl("match control", `disease state:ch1`) ~ "control",
    grepl("psychosis_bipolar", `disease state:ch1`) ~"bipolar disorder",
    grepl("psychosis_schizophrenia", `disease state:ch1`) ~ "schizophrenia",
    TRUE ~ "other"
  ))

#filter bipolar disorder and controls
BD_data = geo_GSE202537_pheno %>%
  filter(grepl("control", diagnosis, ignore.case = T)|
           grepl("bipolar disorder", diagnosis, ignore.case = T))

BD_data_ids = rownames(BD_data)

#clean the column names in BD_data
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "age:ch1",
                                           replacement = 'age')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "age:ch1",
                                           replacement = 'age')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "bmi:ch1",
                                           replacement = 'bmi')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "gender:ch1",
                                           replacement = 'gender')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "manner of death:ch1",
                                           replacement = 'manner_death')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "ph:ch1",
                                           replacement = 'ph')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "pmi:ch1",
                                           replacement = 'pmi')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "race:ch1",
                                           replacement = 'race')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "tissue:ch1",
                                           replacement = 'tissue')
colnames(BD_data) = stri_replace_all_fixed(str=colnames(BD_data),
                                           pattern = "tissuestoragetime:ch1",
                                           replacement = 'tissue_stor_time')

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$age = as.numeric(BD_data$age)
BD_data$bmi = as.numeric(BD_data$bmi)
BD_data$manner_death = factor(BD_data$manner_death)
BD_data$ph = as.numeric(BD_data$ph)
BD_data$pmi = as.numeric(BD_data$pmi)
BD_data$race = factor(BD_data$race)
BD_data$tissue = factor(BD_data$tissue)
BD_data$tissue_stor_time = as.numeric(BD_data$tissue_stor_time)

#create counts matrix
count_matrix = read.table(file = "GSE202537_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = '\t')
rownames(count_matrix) = count_matrix$GeneID

#filter counts matrix 
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
design = model.matrix(~diagnosis+age+gender+manner_death+ph+pmi+race+tissue+tissue_stor_time+bmi, data = BD_data)
#with bmi loose 6 samples

#create DGEList
dge = DGEList(counts = count_matrix_filtered, annotation.columns = "GeneID")
keep = filterByExpr(dge, design) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
all(rownames(dge$samples) == rownames(design)) #FALSE
dge = dge[, intersect(rownames(dge$samples), rownames(design))]
all(rownames(dge$samples) == rownames(design)) #TRUE


v = voom(dge, design, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design)
fit = eBayes(fit)
GSE202537_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 15)]
colnames(SE) = c("SE", "GeneID")
GSE202537_topTable_WITH_covars = inner_join(GSE202537_topTable_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE202537_topTable_WITH_covars = inner_join(GSE202537_topTable_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_WITH_covars = GSE202537_topTable_WITH_covars[, c(1,11, 2:10)]
colnames(GSE202537_topTable_WITH_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_WITH_covars$`Gene Symbol`=="") #FALSE

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
GSE202537_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE202537_topTable_NO_covars = inner_join(GSE202537_topTable_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE202537_topTable_NO_covars = inner_join(GSE202537_topTable_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_NO_covars = GSE202537_topTable_NO_covars[, c(1,11,2:10)]
colnames(GSE202537_topTable_NO_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_NO_covars$`Gene Symbol`=="") #FALSE


#save the results
write.csv(GSE202537_topTable_NO_covars, "GSE202537_table_NO_covars.csv")
write.csv(GSE202537_topTable_WITH_covars, "GSE202537_table_WITH_covars.csv")


###putamen
#filter putamen data
putamen = BD_data %>%
  filter(grepl("Putamen", tissue, ignore.case = T))
putamen_ids = rownames(putamen)


#filter counts 
counts_putamen = count_matrix[, colnames(count_matrix) %in% putamen_ids]
counts_putamen$GeneID = rownames(counts_putamen)

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+age+gender+manner_death+ph+pmi+race+tissue_stor_time+bmi, data = putamen)
counts_putamen_filtered = counts_putamen[, colnames(counts_putamen) %in% rownames(design)]
counts_putamen_filtered$GeneID = rownames(counts_putamen_filtered)

#create DGEList
dge = DGEList(counts = counts_putamen_filtered, annotation.columns = "GeneID")
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
GSE202537_topTable_putamen_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 13)]
colnames(SE) = c("SE", "GeneID")
GSE202537_topTable_putamen_WITH_covars = inner_join(GSE202537_topTable_putamen_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE202537_topTable_putamen_WITH_covars = inner_join(GSE202537_topTable_putamen_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_putamen_WITH_covars = GSE202537_topTable_putamen_WITH_covars[, c(1,11, 2:10)]
colnames(GSE202537_topTable_putamen_WITH_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_putamen_WITH_covars$`Gene Symbol`=="") #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data=putamen)
all(rownames(design0) == rownames(putamen)) #TRUE

#create DGEList
dge0 = DGEList(counts = counts_putamen, annotation.columns = "GeneID")
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
GSE202537_topTable_putamen_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE202537_topTable_putamen_NO_covars = inner_join(GSE202537_topTable_putamen_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE202537_topTable_putamen_NO_covars = inner_join(GSE202537_topTable_putamen_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_putamen_NO_covars = GSE202537_topTable_putamen_NO_covars[, c(1,11,2:10)]
colnames(GSE202537_topTable_putamen_NO_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_putamen_NO_covars$`Gene Symbol`=="") #FALSE


#save the results
write.csv(GSE202537_topTable_putamen_NO_covars, "GSE202537_table_putamen_NO_covars.csv")
write.csv(GSE202537_topTable_putamen_WITH_covars, "GSE202537_table_putamen_WITH_covars.csv")


#caudate
#filter caudate data
caudate = BD_data %>%
  filter(grepl("Caudate", tissue, ignore.case = T))
caudate_ids = rownames(caudate)

#filter counts 
counts_caudate = count_matrix[, colnames(count_matrix) %in% caudate_ids]
counts_caudate$GeneID = rownames(counts_caudate)

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+age+gender+manner_death+ph+pmi+race+tissue_stor_time+bmi, data = caudate)
counts_caudate_filtered = counts_caudate[, colnames(counts_caudate) %in% rownames(design)]
counts_caudate_filtered$GeneID = rownames(counts_caudate_filtered)

#create DGEList
dge = DGEList(counts = counts_caudate_filtered, annotation.columns = "GeneID")
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
GSE202537_topTable_caudate_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 13)]
colnames(SE) = c("SE", "GeneID")
GSE202537_topTable_caudate_WITH_covars = inner_join(GSE202537_topTable_caudate_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE202537_topTable_caudate_WITH_covars = inner_join(GSE202537_topTable_caudate_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_caudate_WITH_covars = GSE202537_topTable_caudate_WITH_covars[, c(1,11, 2:10)]
colnames(GSE202537_topTable_caudate_WITH_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_caudate_WITH_covars$`Gene Symbol`=="") #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data=caudate)
all(rownames(design0) == rownames(caudate)) #TRUE

#create DGEList
dge0 = DGEList(counts = counts_caudate, annotation.columns = "GeneID")
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
GSE202537_topTable_caudate_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE202537_topTable_caudate_NO_covars = inner_join(GSE202537_topTable_caudate_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE202537_topTable_caudate_NO_covars = inner_join(GSE202537_topTable_caudate_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_caudate_NO_covars = GSE202537_topTable_caudate_NO_covars[, c(1,11,2:10)]
colnames(GSE202537_topTable_caudate_NO_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_caudate_NO_covars$`Gene Symbol`=="") #FALSE


#save the results
write.csv(GSE202537_topTable_caudate_NO_covars, "GSE202537_table_caudate_NO_covars.csv")
write.csv(GSE202537_topTable_caudate_WITH_covars, "GSE202537_table_caudate_WITH_covars.csv")


##nucleus accumbens, nac
#filter nac data
nac = BD_data %>%
  filter(grepl("Nac", tissue, ignore.case = T))
nac_ids = rownames(nac)

#filter counts 
counts_nac = count_matrix[, colnames(count_matrix) %in% nac_ids]
counts_nac$GeneID = rownames(counts_nac)

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+age+gender+manner_death+ph+pmi+race+tissue_stor_time+bmi, data = nac)
counts_nac_filtered = counts_nac[, colnames(counts_nac) %in% rownames(design)]
counts_nac_filtered$GeneID = rownames(counts_nac_filtered)

#create DGEList
dge = DGEList(counts = counts_nac_filtered, annotation.columns = "GeneID")
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
GSE202537_topTable_nac_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 13)]
colnames(SE) = c("SE", "GeneID")
GSE202537_topTable_nac_WITH_covars = inner_join(GSE202537_topTable_nac_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE202537_topTable_nac_WITH_covars = inner_join(GSE202537_topTable_nac_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_nac_WITH_covars = GSE202537_topTable_nac_WITH_covars[, c(1,11, 2:10)]
colnames(GSE202537_topTable_nac_WITH_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_nac_WITH_covars$`Gene Symbol`=="") #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data=nac)
all(rownames(design0) == rownames(nac)) #TRUE

#create DGEList
dge0 = DGEList(counts = counts_nac, annotation.columns = "GeneID")
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
GSE202537_topTable_nac_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE202537_topTable_nac_NO_covars = inner_join(GSE202537_topTable_nac_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE202537_topTable_nac_NO_covars = inner_join(GSE202537_topTable_nac_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_nac_NO_covars = GSE202537_topTable_nac_NO_covars[, c(1,11,2:10)]
colnames(GSE202537_topTable_nac_NO_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_nac_NO_covars$`Gene Symbol`=="") #FALSE


#save the results
write.csv(GSE202537_topTable_nac_NO_covars, "GSE202537_table_nac_NO_covars.csv")
write.csv(GSE202537_topTable_nac_WITH_covars, "GSE202537_table_nac_WITH_covars.csv")

##dorsal striatum (putamen + NAc)  - dsm
#filter nac data
dsm = BD_data %>%
  filter(grepl("Putamen", tissue, ignore.case = T)|
           grepl("Caudate", tissue, ignore.case = T))
dsm_ids = rownames(dsm)

#filter counts 
counts_dsm = count_matrix[, colnames(count_matrix) %in% dsm_ids]
counts_dsm$GeneID = rownames(counts_dsm)

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+age+gender+manner_death+ph+pmi+race+tissue_stor_time+bmi, data = dsm)
counts_dsm_filtered = counts_dsm[, colnames(counts_dsm) %in% rownames(design)]
counts_dsm_filtered$GeneID = rownames(counts_dsm_filtered)

#create DGEList
dge = DGEList(counts = counts_dsm_filtered, annotation.columns = "GeneID")
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
GSE202537_topTable_dsm_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 13)]
colnames(SE) = c("SE", "GeneID")
GSE202537_topTable_dsm_WITH_covars = inner_join(GSE202537_topTable_dsm_WITH_covars, SE, by = "GeneID")

#annotate the table
annotation_data$GeneID = factor(annotation_data$GeneID)
GSE202537_topTable_dsm_WITH_covars = inner_join(GSE202537_topTable_dsm_WITH_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_dsm_WITH_covars = GSE202537_topTable_dsm_WITH_covars[, c(1,11, 2:10)]
colnames(GSE202537_topTable_dsm_WITH_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_dsm_WITH_covars$`Gene Symbol`=="") #FALSE

##DE NO COVARS
#create design matrix 
design0 = model.matrix(~diagnosis, data=dsm)
all(rownames(design0) == rownames(dsm)) #TRUE

#create DGEList
dge0 = DGEList(counts = counts_dsm, annotation.columns = "GeneID")
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
GSE202537_topTable_dsm_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$GeneID = rownames(SE0)
SE0 = SE0[, 2:3]
colnames(SE0) = c("SE", "GeneID")
GSE202537_topTable_dsm_NO_covars = inner_join(GSE202537_topTable_dsm_NO_covars, SE0, by = "GeneID")

#annotate the table
GSE202537_topTable_dsm_NO_covars = inner_join(GSE202537_topTable_dsm_NO_covars, annotation_data[, c("GeneID", "Symbol")], by="GeneID")
GSE202537_topTable_dsm_NO_covars = GSE202537_topTable_dsm_NO_covars[, c(1,11,2:10)]
colnames(GSE202537_topTable_dsm_NO_covars)[2] = "Gene Symbol"
any(GSE202537_topTable_dsm_NO_covars$`Gene Symbol`=="") #FALSE


#save the results
write.csv(GSE202537_topTable_dsm_NO_covars, "GSE202537_table_dsm_NO_covars.csv")
write.csv(GSE202537_topTable_dsm_WITH_covars, "GSE202537_table_dsm_WITH_covars.csv")

#clean the tables
#caudate
#no covars 
#delete rows with more than one gene mapped
GSE202537_table_caudate_NO_covars_cleaned = GSE202537_table_caudate_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE202537_table_caudate_NO_covars_cleaned = GSE202537_table_caudate_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE202537_table_caudate_NO_covars_cleaned = GSE202537_table_caudate_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE202537_table_caudate_NO_covars_cleaned$Study = "GSE202537"
write_csv(GSE202537_table_caudate_NO_covars_cleaned, "GSE202537_caudate_no_covars_cleaned.csv")

##putamen
#no covars 
#delete rows with more than one gene mapped
GSE202537_table_putamen_NO_covars_cleaned = GSE202537_table_putamen_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE202537_table_putamen_NO_covars_cleaned = GSE202537_table_putamen_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE202537_table_putamen_NO_covars_cleaned = GSE202537_table_putamen_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE202537_table_putamen_NO_covars_cleaned$Study = "GSE202537"
write_csv(GSE202537_table_putamen_NO_covars_cleaned, "GSE202537_putamen_no_covars_cleaned.csv")
`

#caudate
#with covars 
#delete rows with more than one gene mapped
GSE202537_table_caudate_WITH_covars_cleaned = GSE202537_table_caudate_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE202537_table_caudate_WITH_covars_cleaned = GSE202537_table_caudate_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))


#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE202537_table_caudate_WITH_covars_cleaned = GSE202537_table_caudate_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE202537_table_caudate_WITH_covars_cleaned$Study = "GSE202537"
write_csv(GSE202537_table_caudate_WITH_covars_cleaned, "GSE202537_caudate_with_covars_cleaned.csv")

##putamen
#with covars 
#delete rows with more than one gene mapped
GSE202537_table_putamen_WITH_covars_cleaned = GSE202537_table_putamen_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE202537_table_putamen_WITH_covars_cleaned = GSE202537_table_putamen_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE202537_table_putamen_WITH_covars_cleaned = GSE202537_table_putamen_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE202537_table_putamen_WITH_covars_cleaned$Study = "GSE202537"
write_csv(GSE202537_table_putamen_WITH_covars_cleaned, "GSE202537_putamen_with_covars_cleaned.csv")








