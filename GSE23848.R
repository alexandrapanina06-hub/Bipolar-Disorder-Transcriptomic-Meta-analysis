#Sentrix Human-6 v2 Expression BeadChip

rm(list=ls())
BiocManager :: install (c(" beadarray " , " limma " , " GEOquery " ,
" illuminaHumanv1 . db " , " illuminaHumanv2 . db " , " illuminaHumanv3 . db ",
" BeadArrayUseCases " ))
BiocManager :: install (c(" GOstats " , " GenomicRanges " , " Biostrings " ))
library(beadarray)
BiocManager::install("beadarray")

geo_GSE23848 = getGEO("GSE23848")
geo_GSE23848 = geo_GSE23848[[1]]
pheno_GSE23848 = geo_GSE23848@phenoData@data

#clean colnames
colnames(pheno_GSE23848)[48] = "diagnosis" 
colnames(pheno_GSE23848) = stri_replace_all_fixed(str = colnames(pheno_GSE23848), 
                                                  pattern = "gender:ch1", 
                                                  replacement = "gender")
colnames(pheno_GSE23848) = stri_replace_all_fixed(str = colnames(pheno_GSE23848), 
                                                  pattern = "age (in years):ch1", 
                                                  replacement = "age")
colnames(pheno_GSE23848) = stri_replace_all_fixed(str = colnames(pheno_GSE23848), 
                                                  pattern = "mood stabilizers, antipsychotic medications:ch1", 
                                                  replacement = "medication")
colnames(pheno_GSE23848) = stri_replace_all_fixed(str = colnames(pheno_GSE23848), 
                                                  pattern = "race/ ethnicity:ch1", 
                                                  replacement = "race")
#varibale transformation
pheno_GSE23848$diagnosis = factor(pheno_GSE23848$diagnosis, levels = c("Healthy control", "Bipolar depressed"))
pheno_GSE23848$age = as.numeric(pheno_GSE23848$age)
pheno_GSE23848$gender = factor(pheno_GSE23848$gender)
pheno_GSE23848$race = factor(pheno_GSE23848$race)
pheno_GSE23848$medication = factor(pheno_GSE23848$medication)


#get the raw data
expr_data = read.delim("GSE23848_non_normalized_data.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Remove ProbeID column and make probeIDs rownames
probe_ids = expr_data$ProbeID

# Extract expression columns (non-"Detection" columns)
expr_cols = grep("Detection", colnames(expr_data), invert = TRUE)
exprs = as.matrix(expr_data[, expr_cols])
rownames(exprs) = probe_ids
exprs = exprs[, c(2:36)]

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

#get the annotation file
annotation_data = geo_GSE23848@featureData@data

#DGE WITH covars
#create design matrix
design = model.matrix(~diagnosis+age+gender+race, data = pheno_GSE23848)
all(rownames(design) == rownames(pheno_GSE23848)) #TRUE

#fit the model
fit = lmFit(elist_norm, design)
fit = eBayes(fit)
GSE23848_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE23848_topTable_WITH_covars$ID = rownames(GSE23848_topTable_WITH_covars)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,10)]
colnames(SE) = c("SE", "ID")
GSE23848_topTable_WITH_covars = inner_join(GSE23848_topTable_WITH_covars, SE, by="ID")

#annotate the table
annotation_data$ID = factor(annotation_data$ID)
GSE23848_topTable_WITH_covars = inner_join(GSE23848_topTable_WITH_covars, annotation_data[, c("Symbol", "ID")], by = "ID")
colnames(GSE23848_topTable_WITH_covars)[11] = "Gene Symbol"
GSE23848_topTable_WITH_covars[GSE23848_topTable_WITH_covars$`Gene Symbol`=="",] = NA
GSE23848_topTable_WITH_covars = GSE23848_topTable_WITH_covars[!is.na(GSE23848_topTable_WITH_covars$`Gene Symbol`),]
GSE23848_topTable_WITH_covars = GSE23848_topTable_WITH_covars[, c(9, 11, 1, 10, 2:8)]

#DGE NO covars
#create design matrix
design0 = model.matrix(~diagnosis, data = pheno_GSE23848)
all(rownames(design0) == rownames(pheno_GSE23848)) #TRUE

#fit the model
fit0 = lmFit(elist_norm, design0)
fit0 = eBayes(fit0)
GSE23848_topTable_NO_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)
GSE23848_topTable_NO_covars$ID = rownames(GSE23848_topTable_NO_covars)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, c(2,3)]
colnames(SE0) = c("SE", "ID")
GSE23848_topTable_NO_covars = inner_join(GSE23848_topTable_NO_covars, SE0, by = "ID")

#annotate the table
GSE23848_topTable_NO_covars = inner_join(GSE23848_topTable_NO_covars, annotation_data[, c("Symbol", "ID")], by = "ID")
colnames(GSE23848_topTable_NO_covars)[11] = "Gene Symbol"
GSE23848_topTable_NO_covars[GSE23848_topTable_NO_covars$`Gene Symbol`=="",]=NA
GSE23848_topTable_NO_covars = GSE23848_topTable_NO_covars[!is.na(GSE23848_topTable_NO_covars$`Gene Symbol`),]
GSE23848_topTable_NO_covars = GSE23848_topTable_NO_covars[, c(9, 11, 1, 10, 2:8)]

#save the results
write.csv(GSE23848_topTable_NO_covars, "GSE23848_topTable_NO_covars.csv")
write.csv(GSE23848_topTable_WITH_covars, "GSE23848_topTable_WITH_covars.csv")

#clean the tables
GSE23848_topTable_WITH_covars = GSE23848_topTable_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE23848_topTable_WITH_covars = GSE23848_topTable_WITH_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

#clean the tables
GSE23848_topTable_NO_covars = GSE23848_topTable_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE23848_topTable_NO_covars = GSE23848_topTable_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE23848_topTable_NO_covars$Study = "GSE23848"
GSE23848_topTable_WITH_covars$Study = "GSE23848"

#save the results
write_csv(GSE23848_topTable_NO_covars, "GSE23848_NO_covars.csv")
write_csv(GSE23848_topTable_WITH_covars, "GSE23848_WITH_covars.csv")





