##GSE210064  Agilent-072363 SurePrint G3 Human GE v3 8x60K Microarray 039494
rm(list = ls())
library(GEOquery)
library(R.utils)
library(limma)
library(dplyr)
library(stringi)
library(stringr)

#get the file 
geo_GSE210064 = GEOquery::getGEO("GSE210064")
geo_GSE210064 = geo_GSE210064[[1]]
geo_GSE210064_pheno = geo_GSE210064@phenoData@data

#get raw data
files = list.files("GSE210064_RAW/")
paths = paste0("GSE210064_RAW/", files)
sapply(paths, gunzip)

rawdata_GSE210064 = read.maimages(paths, source = 'agilent', green.only = T)
summary(rawdata_GSE210064$E)
head(rawdata_GSE210064$E)

#get the annotation file
GSE210064_probes = geo_GSE210064@featureData@data

#check distribution before normalization
plot(density(rawdata_GSE210064$E[,1]))
#pre-process the data
rawdata_GSE210064 = backgroundCorrect(rawdata_GSE210064, method = 'normexp')
rawdata_GSE210064 = normalizeBetweenArrays(rawdata_GSE210064, method="quantile")
#after normalization
plot(density(rawdata_GSE210064$E[,1])) #somewhat like log2 scale

rownames(rawdata_GSE210064$E) = rawdata_GSE210064$genes$ProbeName

#filter probes
controls = rawdata_GSE210064$genes$ControlType != 0 #controlType = 0 - regular probe
hist(rowSums(rawdata_GSE210064$E > log2(10)), breaks=20, main="Distribution of Expressed Samples per Probe")
IsExpr = rowSums(rawdata_GSE210064$E > log2(10)) >= 6

#filter rawdata
rawdata_GSE210064 = rawdata_GSE210064[!controls&IsExpr,]

#clean colnames
colnames(rawdata_GSE210064$E) <-  gsub(".*/([^/_]+)_.*", "\\1", colnames(rawdata_GSE210064$E))
all(rownames(geo_GSE210064_pheno) == colnames(rawdata_GSE210064$E)) #TRUE


#create diagnosis column in phenodata
geo_GSE210064_pheno = geo_GSE210064_pheno %>%
  mutate(diagnosis = case_when(
    grepl("Control", source_name_ch1) ~ "control",
    grepl("Bipolar Disorder", source_name_ch1) ~ "bipolar disorder",
    TRUE ~ "Other"
  ))

#variable transformation in phenotype data
geo_GSE210064_pheno$diagnosis = factor(geo_GSE210064_pheno$diagnosis, levels = c("control", "bipolar disorder"))

colnames(geo_GSE210064_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE210064_pheno),
                                                    pattern = 'age:ch1',
                                                    replacement = 'age')
colnames(geo_GSE210064_pheno) = stri_replace_all_fixed(str = colnames(geo_GSE210064_pheno),
                                                       pattern = 'Sex:ch1',
                                                       replacement = 'sex')
geo_GSE210064_pheno$age = as.numeric(geo_GSE210064_pheno$age)
geo_GSE210064_pheno$sex = factor(geo_GSE210064_pheno$sex)


##DE WITH COVARS
#create design matrix 
design = model.matrix(~diagnosis+age+sex, data = geo_GSE210064_pheno)
all(rownames(design) == rownames(geo_GSE210064_pheno)) #TRUE

#fit the model
fit = lmFit(rawdata_GSE210064$E, design)
fit = eBayes(fit)
GSE210064_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,5)]
colnames(SE) = c("SE", "ID")
GSE210064_topTable_WITH_covars = inner_join(GSE210064_topTable_WITH_covars, SE, by = "ID")

#annoatte the table
GSE210064_topTable_WITH_covars = inner_join(GSE210064_topTable_WITH_covars, 
                                            GSE210064_probes[, c("GENE_SYMBOL", "NAME")], by=c("ID" = "NAME"))
GSE210064_topTable_WITH_covars = GSE210064_topTable_WITH_covars[, c(1,11, 2:10)]

#remove rows with missing genes
GSE210064_topTable_WITH_covars[GSE210064_topTable_WITH_covars$GENE_SYMBOL == "", ] = NA
GSE210064_topTable_WITH_covars_filtered = GSE210064_topTable_WITH_covars[!is.na(GSE210064_topTable_WITH_covars$GENE_SYMBOL),]


#DE WITHOUT COVARS
#create design matrix
design_0 = model.matrix(~diagnosis, data = geo_GSE210064_pheno)
all(rownames(design_0) == rownames(geo_GSE210064_pheno)) #TRUE

#fit the model
fit_0 = lmFit(rawdata_GSE210064$E, design_0)
fit_0 = eBayes(fit_0)
GSE210064_topTable_NO_covars = topTable(fit = fit_0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE_0 = as.data.frame(sqrt(fit_0$s2.post)*fit_0$stdev.unscaled)
SE_0$ID = rownames(SE_0)
SE_0 = SE_0[, 2:3]
colnames(SE_0) = c("SE", "ID")
GSE210064_topTable_NO_covars = inner_join(GSE210064_topTable_NO_covars, SE_0, by = "ID")

#annotate the table
GSE210064_topTable_NO_covars = inner_join(GSE210064_topTable_NO_covars, 
                                          GSE210064_probes[, c("GENE_SYMBOL", "NAME")], by=c("ID" = "NAME"))
GSE210064_topTable_NO_covars = GSE210064_topTable_NO_covars[, c(1, 11, 2:10)]

#remove rows with missing genes
GSE210064_topTable_NO_covars[GSE210064_topTable_NO_covars$GENE_SYMBOL=="", ] = NA
GSE210064_topTable_NO_covars_filtered = GSE210064_topTable_NO_covars[!is.na(GSE210064_topTable_NO_covars$GENE_SYMBOL),]


#save the results
write.csv(GSE210064_topTable_NO_covars_filtered, "GSE210064_bale_NOcovars.csv")
write.csv(GSE210064_topTable_WITH_covars_filtered, "GSE210064_bale_withcovars.csv")

#clean the tables

#no covars 
#delete rows with more than one gene mapped
GSE210064_table_NO_covars_cleaned = GSE210064_table_NO_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE210064_table_NO_covars_cleaned = GSE210064_table_NO_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE210064_table_NO_covars_cleaned = GSE210064_table_NO_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 

GSE210064_table_NO_covars_cleaned$Study = "GSE210064"
write_csv(GSE210064_table_NO_covars_cleaned, "GSE210064_no_covars_cleaned.csv")


#with covars
#delete rows with more than one gene mapped
GSE210064_table_WITH_covars_cleaned = GSE210064_table_WITH_covars %>% 
  filter(!grepl("///", `Gene Symbol`))
GSE210064_table_WITH_covars_cleaned = GSE210064_table_WITH_covars_cleaned %>%
  filter(!grepl("---", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE210064_table_WITH_covars_cleaned = GSE210064_table_WITH_covars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE210064_table_WITH_covars_cleaned$Study = "GSE210064"
write_csv(GSE210064_table_WITH_covars_cleaned, "GSE210064_with_covars_cleaned.csv")




