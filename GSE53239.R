rm(list=ls())
##GSE53239

geo_GSE53239 = GEOquery::getGEO("GSE53239")
geo_GSE53239_2 = geo_GSE53239[[2]] #Illumina Genome Analyzer II 
geo_GSE53239_2_pheno = geo_GSE53239_2@phenoData@data

geo_geo_GSE53239_1 = geo_GSE53239[[1]] #Illumina HiSeq 1000
geo_geo_GSE53239_1_pheno = geo_geo_GSE53239_1@phenoData@data

#join two datasets
all(colnames(geo_geo_GSE53239_1_pheno)==colnames(geo_GSE53239_2_pheno)) #TRUE
geo_GSE53239 = bind_rows(geo_geo_GSE53239_1_pheno, geo_GSE53239_2_pheno)

#create diagnosis column
geo_GSE53239 = geo_GSE53239 %>%
  mutate(diagnosis = case_when(
    grepl("BD", `disease state:ch1`) ~ "bipolar disorder",
    grepl("Control", `disease state:ch1`) ~ "control",
    TRUE ~ "other"
  ))
#variable transformation
geo_GSE53239$diagnosis = factor(geo_GSE53239$diagnosis, levels = c("control", "bipolar disorder"))

#create counts matrix
count_matrix = read.table(file = "GSE53239_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = '\t')
rownames(count_matrix) = count_matrix$GeneID

#get annotation file
columns = c('GeneID', 'Symbol', 'chromosome', 'GeneType')
annotation_data = read.columns(file = 'Human.GRCh38.p13.annot.tsv', 
                               required.col = columns, 
                               sep = '\t', 
                               stringsAsFactors = F)
annotation_data$GeneID = factor(annotation_data$GeneID)
colnames(annotation_data)[2] = "Gene Symbol"

##DE NO COVARS
#create design matrix
design = model.matrix(~diagnosis, data = geo_GSE53239)
all(rownames(design)==rownames(geo_GSE53239)) #TRUE

#create DGEList
dge = DGEList(counts = count_matrix, annotation.columns = "GeneID")
keep = filterByExpr(dge, design) #remove very low counts
dge = dge[keep, , keep.lib.sizes = F] #filter dge based on keep and adjust library sizes

#scale normalization
dge = calcNormFactors(dge) #performs TTM normalization
all(rownames(dge$samples) == rownames(design)) #FALSE

#rearrange rows in design matrix
design = design[intersect(colnames(count_matrix), rownames(design)), , drop = FALSE]
all(rownames(dge$samples) == rownames(design)) #TRUE

v = voom(dge, design, plot = T) #convert data to log2 scale, estimates mean-variance trend
all(rownames(v$targets) == rownames(v$design)) #TRUE

#fit the model
fit = lmFit(v, design)
fit = eBayes(fit)
GSE53239_topTable_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$GeneID = rownames(SE)
SE = SE[, c(2, 3)]
colnames(SE) = c("SE", "GeneID")
GSE53239_topTable_NO_covars$GeneID = factor(GSE53239_topTable_NO_covars$GeneID)
GSE53239_topTable_NO_covars = inner_join(GSE53239_topTable_NO_covars, SE, by="GeneID")

#annotate the table
GSE53239_topTable_NO_covars = inner_join(GSE53239_topTable_NO_covars, annotation_data[, c("Gene Symbol", "GeneID")], by="GeneID")
GSE53239_topTable_NO_covars = GSE53239_topTable_NO_covars[, c(1,11,2:10)]
any(is.na(GSE53239_topTable_NO_covars$`Gene Symbol`)) #FALSE

#save the results
write_csv(GSE53239_topTable_NO_covars, "GSE53239_table_NO_covars.csv")


#clean the tables
#delete rows with more than one gene mapped
GSE53239_table_NO_covars_cleaned = GSE53239_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE53239_table_NO_covars_cleaned = GSE53239_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE53239_table_NO_covars_cleaned = GSE53239_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE53239_table_NO_covars_cleaned$Study = "GSE53239"
write_csv(GSE53239_table_NO_covars_cleaned, "GSE53239_no_covars_cleaned.csv")











