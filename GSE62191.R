rm(list=ls())
##GSE62191 Agilent-014850 Whole Human Genome Microarray 4x44K G4112F
library(limma)
#get the file 
geo_GSE62191 = GEOquery::getGEO("GSE62191")
geo_GSE62191 = geo_GSE62191[[1]]
geo_GSE62191_pheno = geo_GSE62191@phenoData@data

#filter BD and controls
BD_data = geo_GSE62191_pheno %>%
  filter(grepl("bipolar disorder", characteristics_ch1.1, ignore.case = T)|
           grepl("control", characteristics_ch1.1, ignore.case = T))
BD_data_ids = rownames(BD_data)


#get the raw data
files = list.files("GSE62191_RAW")
paths = paste0("GSE62191_RAW/", files)
sapply(paths, gunzip)

files_filtered = files[grepl(paste(BD_data_ids, collapse = "|"), files)]
paths_filtered = paste0("GSE62191_RAW/", files_filtered)

rawdata_GSE62191 = read.maimages(paths_filtered, source = 'agilent', green.only = T)

#data pre-processing
rawdata_GSE62191 = backgroundCorrect(rawdata_GSE62191, method = "normexp")
rawdata_GSE62191 = normalizeBetweenArrays(rawdata_GSE62191, method = "quantile")
plot(density(rawdata_GSE62191$E[,1]))

rownames(rawdata_GSE62191$E) = rawdata_GSE62191$genes$ProbeName

#filter probes
controls = rawdata_GSE62191$genes$ControlType != 0
hist(rowSums(rawdata_GSE62191$E > log2(10)), breaks=20, main="Distribution of Expressed Samples per Probe")
IsExpr = rowSums(rawdata_GSE62191$E > log2(10)) >= 6

#filter raw data
rawdata_GSE62191 = rawdata_GSE62191[!controls&IsExpr, ]

#clean colnames
colnames(rawdata_GSE62191$E) <-  gsub(".*/([^/_]+)_.*", "\\1", colnames(rawdata_GSE62191$E))
all(rownames(BD_data) == colnames(rawdata_GSE62191$E)) #TRUE

#create diagnosis column in BD_data
BD_data = BD_data %>%
  mutate(diagnosis = case_when(
    grepl("bipolar disorder", characteristics_ch1.1) ~ "bipolar disorder",
    grepl("control", characteristics_ch1.1) ~ "control",
    TRUE ~ 'Other'
  ))

#clean colnames
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                                       pattern = 'age:ch1',
                                                       replacement = 'age')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'age of onset:ch1',
                                           replacement = 'age_onset')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'dsm-iv:ch1',
                                           replacement = 'dsm-iv')
colnames(BD_data) = stri_replace_all_fixed(str = colnames(BD_data),
                                           pattern = 'population:ch1',
                                           replacement = 'population')

BD_data$age = gsub("yr", "", BD_data$age)
BD_data$age_onset = gsub("yr", "", BD_data$age_onset)

#variable transformation
BD_data$diagnosis = factor(BD_data$diagnosis, levels = c("control", "bipolar disorder"))
BD_data$population = factor(BD_data$population)
BD_data$age = as.numeric(BD_data$age)
BD_data$age_onset = as.numeric(BD_data$age_onset)
BD_data$`dsm-iv`= as.numeric(BD_data$`dsm-iv`)

##DE WITH COVARS
#create design matrix
design = model.matrix(~diagnosis+population+age, data = BD_data)
all(rownames(design)==rownames(BD_data)) #TRUE

#fit the model
fit = lmFit(rawdata_GSE62191$E, design)
fit = eBayes(fit)
GSE62191_topTable_WITH_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,5)]
colnames(SE) = c("SE", "ID")
GSE62191_topTable_WITH_covars = inner_join(GSE62191_topTable_WITH_covars, SE, by="ID")

#annotate the table
#get annotation file
GSE62191_probes = geo_GSE210064@featureData@data

GSE62191_topTable_WITH_covars = inner_join(GSE62191_topTable_WITH_covars, 
                                           GSE62191_probes[, c("GENE_SYMBOL", "NAME")], by = c("ID"="NAME"))
GSE62191_topTable_WITH_covars = GSE62191_topTable_WITH_covars[, c(1,11,2:10)]

#remove rows with missing gene symbols
GSE62191_topTable_WITH_covars$GENE_SYMBOL[GSE62191_topTable_WITH_covars$GENE_SYMBOL==""] = NA
GSE62191_topTable_WITH_covars_filtered = GSE62191_topTable_WITH_covars[!is.na(GSE62191_topTable_WITH_covars$GENE_SYMBOL),]

##DE NO COVARS
#create design matrix 
design_0 = model.matrix(~diagnosis, data = BD_data)
all(rownames(design_0) == rownames(BD_data)) #TRUE

#fit the model
fit_0 = lmFit(rawdata_GSE62191$E, design_0)
fit_0 = eBayes(fit_0)
GSE62191_topTable_NO_covars = topTable(fit = fit_0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE_0 = as.data.frame(sqrt(fit_0$s2.post)*fit_0$stdev.unscaled)
SE_0$ID = rownames(SE_0)
SE_0 = SE_0[, 2:3]
colnames(SE_0) = c("SE", "ID")

GSE62191_topTable_NO_covars = inner_join(GSE62191_topTable_NO_covars, SE_0, by = "ID")

#annotate the table
GSE62191_topTable_NO_covars = inner_join(GSE62191_topTable_NO_covars, 
                                         GSE62191_probes[, c("GENE_SYMBOL", "NAME")], by = c("ID"="NAME"))
GSE62191_topTable_NO_covars = GSE62191_topTable_NO_covars[, c(1,11,2:10)]

#remove rows with missing gene symbols
GSE62191_topTable_NO_covars$GENE_SYMBOL[GSE62191_topTable_NO_covars$GENE_SYMBOL==""] = NA
GSE62191_topTable_NO_covars_filtered = GSE62191_topTable_NO_covars[!is.na(GSE62191_topTable_NO_covars$GENE_SYMBOL),]


#save results
write_csv(GSE62191_topTable_NO_covars_filtered, "GSE62191_array_table_nocovars.csv")
write_csv(GSE62191_topTable_WITH_covars_filtered, "GSE62191_array_table_withcovars.csv")

#clean the tables 
colnames(GSE62191_array_table_nocovars)[2] = "Gene ID"
colnames(GSE62191_array_table_withcovars)[2] = "Gene ID"
colnames(GSE62191_array_table_nocovars)[3] = "Gene Symbol"
colnames(GSE62191_array_table_withcovars)[3] = "Gene Symbol"

#no covars 
#delete rows with more than one gene mapped
GSE62191_array_table_nocovars_cleaned = GSE62191_array_table_nocovars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE62191_array_table_nocovars_cleaned = GSE62191_array_table_nocovars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE62191_array_table_nocovars_cleaned$Study = "GSE62191"
write_csv(GSE62191_array_table_nocovars_cleaned, "GSE62191_no_covars_cleaned.csv")


#with covars 
#delete rows with more than one gene mapped
GSE62191_array_table_withcovars_cleaned = GSE62191_array_table_withcovars %>% 
  filter(!grepl("///", `Gene Symbol`))

#for genes that are mapped to more than one ID
#assign average logFC and maximal error

#modify logFC and SE
GSE62191_array_table_withcovars_cleaned = GSE62191_array_table_withcovars_cleaned %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE)
  ) 
GSE62191_array_table_withcovars_cleaned$Study = "GSE62191"
write_csv(GSE62191_array_table_withcovars_cleaned, "GSE62191_with_covars_cleaned.csv")


