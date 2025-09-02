rm(list = ls())
#Illumina HumanHT-12 V4.0 expression beadchip
###The expression matrix doesn't have GSMXXX column, but only SAMPLEX
#let's assume, that X corresponds to the row number

geo_GSE39653 = getGEO("GSE39653")
geo_GSE39653 = geo_GSE39653[[1]]
pheno_GSE39653 = geo_GSE39653@phenoData@data

#clean the colnames
colnames(pheno_GSE39653)[32] = "diagnosis"

#filter only BD and controls
BD_data = pheno_GSE39653 %>%
  filter(grepl("bipolar disorder", diagnosis, ignore.case = T)|
           grepl("healthy control", diagnosis, ignore.case = T))
BD_data_ids = rownames(BD_data)

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("healthy control", "bipolar disorder"))

expr_data = read.delim("GSE39653_non_normalized.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Remove ProbeID column and make probeIDs rownames
probe_ids = expr_data$ID_REF
expr_data = expr_data[, -1]

# Extract expression columns (non-"Detection" columns)
expr_cols = grep("^SAMPLE", colnames(expr_data))
exprs = as.matrix(expr_data[, expr_cols])
rownames(exprs) = probe_ids


# Extract detection p-value columns
det_cols = grep("Detection", colnames(expr_data))
detection_pvals = as.matrix(expr_data[, det_cols])
rownames(detection_pvals) = probe_ids


elist_raw = list(E = exprs,
                 other = list(Detection = detection_pvals))
class(elist_raw) = "EListRaw"

#dimensions check
dim(exprs)
dim(detection_pvals)

#Apply the neqc function to calibrate the background level, normalize and transform
#the intensities from each sample
elist_norm = neqc(elist_raw)
colnames(elist_norm$E) = rownames(pheno_GSE39653)

#filter elist based on BD-data ids
elist_norm = elist_norm[,colnames(elist_norm) %in% BD_data_ids]

#DGE NO covars
#create design matrix
design = model.matrix(~diagnosis, data = BD_data)
all(rownames(design)==rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(elist_norm, design)
fit = eBayes(fit)
GSE39653_topTable_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE39653_topTable_NO_covars$ID = rownames(GSE39653_topTable_NO_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,3)]
colnames(SE) = c("SE", "ID")
GSE39653_topTable_NO_covars = inner_join(GSE39653_topTable_NO_covars, SE, by = "ID")

#annotate the table
annotation_data = geo_GSE39653@featureData@data
GSE39653_topTable_NO_covars = inner_join(GSE39653_topTable_NO_covars, annotation_data[, c("ILMN_Gene", "ID")], by = "ID")
colnames(GSE39653_topTable_NO_covars)[11] = "Gene Symbol"
GSE39653_topTable_NO_covars = GSE39653_topTable_NO_covars[!is.na(GSE39653_topTable_NO_covars$`Gene Symbol`),]
GSE39653_topTable_NO_covars = GSE39653_topTable_NO_covars[, c(9,11,1,10, 2:8)]

#save the results
write.csv(GSE39653_topTable_NO_covars, "GSE39653_topTable_NO_covars.csv")

#clean the tables
GSE39653_topTable_NO_covars = GSE39653_topTable_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE39653_topTable_NO_covars = GSE39653_topTable_NO_covars %>% 
  filter(!grepl("Mar", `Gene Symbol`))
GSE39653_topTable_NO_covars = GSE39653_topTable_NO_covars %>% 
  filter(!grepl("Dec", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE39653_topTable_NO_covars = GSE39653_topTable_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

#add study column to the table
GSE39653_topTable_NO_covars$Study = "GSE39653"

#save
write_csv(GSE39653_topTable_NO_covars, "GSE39653_NO_covars.csv")





