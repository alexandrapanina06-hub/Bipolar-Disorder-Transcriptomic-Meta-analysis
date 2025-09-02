#Enrichment analysis
library(org.Hs.eg.db)
library(clusterProfiler)
library(hgu133plus2cdf)
library(AnnotationDbi)
library(purrr)

#meta-analysis A, no covars
#create gene universe (unique genes in metadata)
gene_universe = unique(metadata_no_covars$`Gene Symbol`)
head(gene_universe) 
#gene symbols are gene names, cluster profiler needs ENTREZ ID

#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                   keys=gene_universe, 
                   column="ENTREZID", 
                   keytype="SYMBOL")
#create list of dif expr genes
genes_de = meta_analysis_no_covars[meta_analysis_no_covars$meta_pval<0.05,]$Gene.Symbol
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
write_csv(enrichment@result, "enrichment_BP.csv")

cnetplot(enrichment,
         circulr = FALSE,
         colorEdge = TRUE,
         color_category = 'firebrick',
         color_gene = 'steelblue',
         layout = 'kk',
         showCategory = 10,
         cex_gene = 0.7,
         cex_label_gene = 0.8,
         shadowtex = "gene",
         max.overlaps = 1000,
         force = 5,
         force_pull = 0.8,
         max.time = 2)

heatplot(enrichment, showCategory = 5)
edo = pairwise_termsim(enrichment)
treeplot(edo)
upsetplot(enrichment)


###Meta-analysis A with covars
gene_universe_covars = unique(metadata_with_covars$`Gene Symbol`)
#convert gene symbol to ENTREZ ID
genes_entrezid_covars = mapIds(org.Hs.eg.db, 
                        keys=gene_universe_covars, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create list of dif expr genes
genes_de_sign_cov = meta_analysis_with_covars[meta_analysis_with_covars$meta_pval<0.05,]$Gene.Symbol
genes_de_sign_cov = mapIds(org.Hs.eg.db, 
                  keys=genes_de_sign_cov, 
                  column="ENTREZID", 
                  keytype="SYMBOL")

#enrichment for biological process (BP)
enrichment_cov = clusterProfiler::enrichGO(gene = genes_de_sign_cov,
                                       OrgDb = "org.Hs.eg.db",
                                       ont = "BP",
                                       universe = genes_entrezid_covars,
                                       minGSSize=5)

enrichment_cov = setReadable(enrichment_cov, "org.Hs.eg.db", "ENTREZID")
barplot(enrichment_cov, showCategory = 10, font.size=10)
#save the results
write.csv(enrichment_cov@result, "enrichment_cov.csv")


cnetplot(enrichment_cov,
         circulr = FALSE,
         colorEdge = TRUE,
         color_category = 'firebrick',
         color_gene = 'steelblue',
         layout = 'kk',
         showCategory = 10,
         cex_gene = 0.7,
         cex_label_gene = 0.8,
         shadowtex = "gene",
         max.overlaps = 1000,
         force = 5,
         force_pull = 0.8,
         max.time = 2)
edo2 = pairwise_termsim(enrichment_cov)
treeplot(edo2)
upsetplot(enrichment_cov)

#KEGG pathway analysis 
#without covariates
kegg = enrichKEGG(gene = genes_de, 
                  organism = "hsa",
                  pvalueCutoff = 0.05)
dotplot(kegg, showCategory = 10)

#with covariates
kegg_covar = enrichKEGG(gene = genes_de_sign_cov, 
                        organism = "hsa",
                        pvalueCutoff = 0.05)
dotplot(kegg_covar, showCategory = 10)



###BA9 enrichment analysis (meta-analysis B)
#no covariates 
gene_universe = unique(ba9_metadata_no_covars$`Gene Symbol`)

#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                        keys=gene_universe, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create list of dif expr genes
genes_ba9_de = ba9_no_covars[ba9_no_covars$meta_pval<0.05,]$Gene.Symbol
genes_ba9_de = mapIds(org.Hs.eg.db, 
                  keys=genes_ba9_de, 
                  column="ENTREZID", 
                  keytype="SYMBOL")
#enrichment for biological process (BP)
enrichment_ba9 = clusterProfiler::enrichGO(gene = genes_ba9_de,
                                       OrgDb = "org.Hs.eg.db",
                                       ont = "BP",
                                       universe = genes_entrezid,
                                       minGSSize=5)

enrichment_ba9 = setReadable(enrichment_ba9, "org.Hs.eg.db", "ENTREZID")
barplot(enrichment_ba9, showCategory = 10, font.size=15)
write_csv(enrichment_ba9@result, "enrichment_BP_ba9.csv")

kegg_ba9 = enrichKEGG(gene = genes_ba9_de, 
                      organism = "hsa",
                      pvalueCutoff = 0.05)
dotplot(kegg_ba9, showCategory = 10)

cnetplot(enrichment_ba9,
         circulr = FALSE,
         colorEdge = TRUE,
         color_category = 'firebrick',
         color_gene = 'steelblue',
         layout = 'kk',
         showCategory = 10,
         cex_gene = 0.7,
         cex_label_gene = 0.8,
         shadowtex = "gene",
         max.overlaps = 1000,
         force = 5,
         force_pull = 0.8,
         max.time = 2)
treeplot(pairwise_termsim(enrichment_ba9))
upsetplot(enrichment_ba9)

#with covariates
gene_universe_covars = unique(ba9_metadata_with_covars$`Gene Symbol`)
#convert gene symbol to ENTREZ ID
genes_entrezid_covars = mapIds(org.Hs.eg.db, 
                               keys=gene_universe_covars, 
                               column="ENTREZID", 
                               keytype="SYMBOL")
#create a list of differentially expressed genes
genes_ba9_de_cov = ba9_with_covars[ba9_with_covars$meta_pval<0.05,]$Gene.Symbol
genes_ba9_de_cov = mapIds(org.Hs.eg.db, 
                      keys=genes_ba9_de_cov, 
                      column="ENTREZID", 
                      keytype="SYMBOL")
#enrichment for biological process (BP)
enrichment_ba9_cov = clusterProfiler::enrichGO(gene = genes_ba9_de_cov,
                                           OrgDb = "org.Hs.eg.db",
                                           ont = "BP",
                                           universe = genes_entrezid,
                                           minGSSize=5)

enrichment_ba9_cov = setReadable(enrichment_ba9_cov, "org.Hs.eg.db", "ENTREZID")
barplot(enrichment_ba9_cov, showCategory = 10, font.size=15)
write_csv(enrichment_ba9_cov@result, "enrichment_BP_ba9_cov.csv")

##KEGG pathway enrichment for ba9
kegg_ba9_cov = enrichKEGG(gene = genes_ba9_de_cov, 
                      organism = "hsa",
                      pvalueCutoff = 0.05)
dotplot(kegg_ba9_cov, showCategory = 10)



###Enrichment analysis for BA9 + BA46 (pfc)
#no covariates 
gene_universe = unique(pfc_metadata_nocov$`Gene Symbol`)

#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                        keys=gene_universe, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create list of dif expr genes
genes_pfc_de = pfc_no_covars[pfc_no_covars$meta_pval<0.05, ]$Gene.Symbol
genes_pfc_de = mapIds(org.Hs.eg.db, 
                      keys=genes_pfc_de, 
                      column="ENTREZID", 
                      keytype="SYMBOL")
#enrichment for biological process (BP)
enrichment_pfc = clusterProfiler::enrichGO(gene = genes_pfc_de,
                                           OrgDb = "org.Hs.eg.db",
                                           ont = "BP",
                                           universe = genes_entrezid,
                                           minGSSize=5)

enrichment_pfc = setReadable(enrichment_pfc, "org.Hs.eg.db", "ENTREZID")
barplot(enrichment_pfc, showCategory = 10, font.size=15)
write_csv(enrichment_pfc@result, "enrichment_BP_pfc.csv")
kegg_pfc = enrichKEGG(gene = genes_pfc_de, 
                      organism = "hsa",
                      pvalueCutoff = 0.05)
dotplot(kegg_pfc, showCategory = 10)

cnetplot(enrichment_pfc,
         circulr = FALSE,
         colorEdge = TRUE,
         color_category = 'firebrick',
         color_gene = 'steelblue',
         layout = 'kk',
         showCategory = 10,
         cex_gene = 0.7,
         cex_label_gene = 0.8,
         shadowtex = "gene",
         max.overlaps = 1000,
         force = 5,
         force_pull = 0.8,
         max.time = 2)
treeplot(pairwise_termsim(enrichment_pfc))
upsetplot(enrichment_pfc)

#with covariates
gene_universe = unique(pfc_metadata_withcovar$`Gene Symbol`)

#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                        keys=gene_universe, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create a list of differentially expressed genes
genes_pfc_de_cov = pfc_with_covars[pfc_with_covars$meta_pval<0.05,]$Gene.Symbol
genes_pfc_de_cov = mapIds(org.Hs.eg.db, 
                          keys=genes_pfc_de_cov, 
                          column="ENTREZID", 
                          keytype="SYMBOL")
#enrichment for biological process (BP)
enrichment_pfc_cov = clusterProfiler::enrichGO(gene = genes_pfc_de_cov,
                                               OrgDb = "org.Hs.eg.db",
                                               ont = "BP",
                                               universe = genes_entrezid,
                                               minGSSize=5)

enrichment_pfc_cov = setReadable(enrichment_pfc_cov, "org.Hs.eg.db", "ENTREZID")
barplot(enrichment_pfc_cov, showCategory = 10, font.size=15)
write_csv(enrichment_pfc_cov@result, "enrichment_BP_pfc_cov.csv")

##KEGG pathway enrichment for ba9 + ba46
kegg_pfc_cov = enrichKEGG(gene = genes_pfc_de_cov, 
                          organism = "hsa",
                          pvalueCutoff = 0.05)
dotplot(kegg_pfc_cov, showCategory = 10)

cnetplot(enrichment_pfc_cov,
         circulr = FALSE,
         colorEdge = TRUE,
         color_category = 'firebrick',
         color_gene = 'steelblue',
         layout = 'kk',
         showCategory = 10,
         cex_gene = 0.7,
         cex_label_gene = 0.8,
         shadowtex = "gene",
         max.overlaps = 1000,
         force = 5,
         force_pull = 0.8,
         max.time = 2)
treeplot(pairwise_termsim(enrichment_pfc_cov))
upsetplot(enrichment_pfc_cov)




#Network analysis (meta-analysis A, no covars, significant genes)
library('wordcloud2')
library('webshot')
library('htmlwidgets')
library(networkD3)
webshot::install_phantomjs()

subset_data = meta_analysis_no_covars_significant[, c("Gene.Symbol", "meta_LFc")]

simpleNetwork(subset_data, charge = -1000, width = 2000,
              height = 2000, fontSize = 20, opacity = 1) %>%
  saveWidget(file = 'network_A_nocov_signif.html')
webshot("network_A_nocov_signif.html", "network_A_nocov_signif.png", delay = 5, vwidth = 2000, vheight = 2000)





