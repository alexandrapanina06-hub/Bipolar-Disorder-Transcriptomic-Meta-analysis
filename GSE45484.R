#GSE45484
#Illumina HumanHT-12 V4.0 expression beadchip

library(GEOquery)
library(limma)
library(edgeR)

geo_GSE45484 = getGEO("GSE45484")
geo_GSE45484 = geo_GSE45484[[1]]
pheno_GSE45484 = geo_GSE45484@phenoData@data

#clean colnames
colnames(pheno_GSE45484) = gsub(":ch1", "", colnames(pheno_GSE45484))

#variable transformation
pheno_GSE45484$age = as.numeric(pheno_GSE45484$age)
pheno_GSE45484$sex = factor(pheno_GSE45484$sex)
pheno_GSE45484$`time point` = factor(pheno_GSE45484$`time point`)
pheno_GSE45484$`treatment group` = factor(pheno_GSE45484$`treatment group`)
pheno_GSE45484$responder = factor(pheno_GSE45484$responder, levels = c("NO", "YES"))

#get expression matrix
expr_data = read.delim("GSE45484_non-normalized.txt", header = TRUE, sep = "\t", check.names = FALSE)
expr_ids = expr_data$`GEO ID`
expr = as.matrix(expr_data[2:nrow(expr_data), 2:ncol(expr_data)])
rownames(expr) = expr_ids[2:47324]

all(colnames(expr) %in% rownames(pheno_GSE45484)) #TRUE
expr = expr[, rownames(pheno_GSE45484)]
storage.mode(expr) = "numeric"

expr = log2(expr + 0.00001)

expr = normalizeBetweenArrays(expr, method = "quantile")

#DE without covariates
#create design matrix
design = model.matrix(~responder, data = pheno_GSE45484)
all(rownames(design) == pheno_GSE45484$geo_accession) #TRUE

#fit the linear model
fit = lmFit(expr, design)
fit = eBayes(fit)
GSE45484_topTable_NO_covars = topTable(fit = fit, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE = as.data.frame(sqrt(fit$s2.post)*fit$stdev.unscaled)
SE$ID = rownames(SE)
SE = SE[, c(2,3)]
colnames(SE) = c("SE", "ID")

GSE45484_topTable_NO_covars = GSE45484_topTable_NO_covars %>% tibble::rownames_to_column("ID")
GSE45484_topTable_NO_covars = as.data.frame(GSE45484_topTable_NO_covars)
GSE45484_topTable_NO_covars = inner_join(GSE45484_topTable_NO_covars, SE, by = "ID")

#annotate the table
GSE45484_probes = geo_GSE45484@featureData@data
GSE45484_topTable_NO_covars = inner_join(GSE45484_topTable_NO_covars, GSE45484_probes[, c("ID", "Symbol")], by = "ID")
GSE45484_topTable_NO_covars = GSE45484_topTable_NO_covars[, c(1, 11, 2:10)]

#merge various IDs mapped to the same gene
GSE45484_topTable_NO_covars = GSE45484_topTable_NO_covars %>% 
  filter(!grepl("///", Symbol))
GSE45484_topTable_NO_covars$Symbol = na_if(GSE45484_topTable_NO_covars$Symbol, "")
GSE45484_topTable_NO_covars = GSE45484_topTable_NO_covars[!is.na(GSE45484_topTable_NO_covars$Symbol),]


GSE45484_topTable_NO_covars = GSE45484_topTable_NO_covars %>%
  group_by(Symbol) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE), 
    p_value = min(P.Value)
  ) 

GSE45484_topTable_NO_covars = GSE45484_topTable_NO_covars[14:nrow(GSE45484_topTable_NO_covars),]
write_csv(GSE45484_topTable_NO_covars, "GSE45484_no_covars.csv")

#with covariates
#create design matrix
design0 = model.matrix(~responder + age + sex + `time point` + `treatment group`, data = pheno_GSE45484)
all(rownames(design0) == pheno_GSE45484$geo_accession) #TRUE

#fit the linear model
fit0 = lmFit(expr, design0)
fit0 = eBayes(fit0)
GSE45484_topTable_WITH_covars = topTable(fit = fit0, coef = 2, adjust.method = "fdr", number = Inf, confint = T)

#add SE to the table
SE0 = as.data.frame(sqrt(fit0$s2.post)*fit0$stdev.unscaled)
SE0$ID = rownames(SE0)
SE0 = SE0[, c(2,7)]
colnames(SE0) = c("SE", "ID")

GSE45484_topTable_WITH_covars = GSE45484_topTable_WITH_covars %>% tibble::rownames_to_column("ID")
GSE45484_topTable_WITH_covars = as.data.frame(GSE45484_topTable_WITH_covars)
GSE45484_topTable_WITH_covars = inner_join(GSE45484_topTable_WITH_covars, SE0, by = "ID")

#annotate the table
GSE45484_probes = geo_GSE45484@featureData@data
GSE45484_topTable_WITH_covars = inner_join(GSE45484_topTable_WITH_covars, GSE45484_probes[, c("ID", "Symbol")], by = "ID")
GSE45484_topTable_WITH_covars = GSE45484_topTable_WITH_covars[, c(1, 11, 2:10)]

#merge various IDs mapped to the same gene
GSE45484_topTable_WITH_covars = GSE45484_topTable_WITH_covars %>% 
  filter(!grepl("///", Symbol))
GSE45484_topTable_WITH_covars$Symbol = na_if(GSE45484_topTable_WITH_covars$Symbol, "")
GSE45484_topTable_WITH_covars = GSE45484_topTable_WITH_covars[!is.na(GSE45484_topTable_WITH_covars$Symbol),]


GSE45484_topTable_WITH_covars = GSE45484_topTable_WITH_covars %>%
  group_by(Symbol) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE), 
    p_value = min(P.Value)
  ) 
GSE45484_topTable_WITH_covars = GSE45484_topTable_WITH_covars[14:nrow(GSE45484_topTable_WITH_covars),]
write_csv(GSE45484_topTable_WITH_covars, "GSE45484_with_covars.csv")

no_covars_sign = GSE45484_topTable_NO_covars[GSE45484_topTable_NO_covars$p_value < 0.05 & 
                                               abs(GSE45484_topTable_NO_covars$logFC) > 0.2, ]
with_covars_sign = GSE45484_topTable_WITH_covars[GSE45484_topTable_WITH_covars$p_value < 0.05 &
                                                   abs(GSE45484_topTable_WITH_covars$logFC) > 0.2,]


intersect(blood_no_covars[blood_no_covars$meta_pval < 0.05,]$Gene.Symbol, no_covars_sign$Symbol)
intersect(blood_with_covars[blood_with_covars$meta_pval < 0.05,]$Gene.Symbol, with_covars_sign$Symbol)

write_csv(no_covars_sign, "no_covars_sign.csv")
write_csv(with_covars_sign, "with_covars_sign.csv")

#Functional enrichment analysis
library(org.Hs.eg.db)
library(clusterProfiler)
library(hgu133plus2cdf)
library(AnnotationDbi)
library(purrr)

#no covars
gene_universe = unique(GSE45484_topTable_NO_covars$Symbol)
head(gene_universe) 
#gene symbols are gene names, cluster profiler needs ENTREZ ID

#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                        keys=gene_universe, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create list of dif expr genes
genes_de = no_covars_sign$Symbol
genes_de = mapIds(org.Hs.eg.db, 
                  keys=genes_de, 
                  column="ENTREZID", 
                  keytype="SYMBOL")
#enrichment for biological process (BP)
enrichment = clusterProfiler::enrichGO(gene = genes_de,
                                       OrgDb = "org.Hs.eg.db",
                                       ont = "BP",
                                       universe = genes_entrezid,
                                       minGSSize=5)
ego = as.data.frame(enrichment) #check geneRatio and Bgratio
enrichment = setReadable(enrichment, "org.Hs.eg.db", "ENTREZID")

barplot(enrichment, showCategory = 10, font.size=15)
write_csv(enrichment@result, "enrichment_nocovars.csv")

kegg = enrichKEGG(gene = genes_de, 
                  organism = "hsa",
                  pvalueCutoff = 0.05)
dotplot(kegg, showCategory = 10)

#with covars
gene_universe = unique(GSE45484_topTable_WITH_covars$Symbol)
head(gene_universe) 
#gene symbols are gene names, cluster profiler needs ENTREZ ID

#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                        keys=gene_universe, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create list of dif expr genes
genes_de = with_covars_sign$Symbol
genes_de = mapIds(org.Hs.eg.db, 
                  keys=genes_de, 
                  column="ENTREZID", 
                  keytype="SYMBOL")
#enrichment for biological process (BP)
enrichment = clusterProfiler::enrichGO(gene = genes_de,
                                       OrgDb = "org.Hs.eg.db",
                                       ont = "BP",
                                       universe = genes_entrezid,
                                       minGSSize=5)
ego = as.data.frame(enrichment) #check geneRatio and Bgratio
enrichment = setReadable(enrichment, "org.Hs.eg.db", "ENTREZID")

barplot(enrichment, showCategory = 10, font.size=15)
write_csv(enrichment@result, "enrichment_withcovars.csv")

kegg = enrichKEGG(gene = genes_de, 
                  organism = "hsa",
                  pvalueCutoff = 0.05)
dotplot(kegg, showCategory = 10)

