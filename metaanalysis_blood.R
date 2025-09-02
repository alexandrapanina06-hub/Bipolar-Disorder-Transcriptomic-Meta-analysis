library(dplyr)
library(purrr)
BiocManager::install("metafor")
library(metafor)
library(ggrepel)
library(gridExtra)
library(grid)
library(edgeR)
library(R.utils)

#curate genotype in each table
curated_genes = function(table) {
  reference = Homo_sapiens %>% dplyr::select(V5,V3)
  table = table %>%
    left_join(reference, by = c("Gene Symbol" = "V5")) %>%  # Match old names to gene column
    mutate(`Gene Symbol` = coalesce(V3, `Gene Symbol`)) %>%  # Replace with new names where available
    dplyr::select(-V3)  # Remove the reference column
  return(table)
}

#GSE124326
GSE124326_NO_covars = curated_genes(GSE124326_NO_covars)
GSE124326_NO_covars = GSE124326_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE124326"
  ) 

GSE124326_WITH_covars = curated_genes(GSE124326_WITH_covars)
GSE124326_WITH_covars = GSE124326_WITH_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE124326"
  ) 

#GSE18312
GSE18312_NO_covars = curated_genes(GSE18312_NO_covars)
GSE18312_NO_covars = GSE18312_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE18312"
  ) 

GSE18312_WITH_covars = curated_genes(GSE18312_WITH_covars)
GSE18312_WITH_covars = GSE18312_WITH_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE18312"
  ) 

#GSE23848
GSE23848_NO_covars = curated_genes(GSE23848_NO_covars)
GSE23848_NO_covars = GSE23848_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE23848"
  ) 
GSE23848_WITH_covars = curated_genes(GSE23848_WITH_covars)
GSE23848_WITH_covars = GSE23848_WITH_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE23848"
  ) 

#GSE46416
GSE46416_NO_covars = curated_genes(GSE46416_NO_covars)
GSE46416_NO_covars = GSE46416_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE46416"
  ) 

#GSE39653
GSE39653_NO_covars = curated_genes(GSE39653_NO_covars)
GSE39653_NO_covars = GSE39653_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE39653"
  ) 

#GSE46449
GSE46449_NO_covars = curated_genes(GSE46449_NO_covars)
GSE46449_NO_covars = GSE46449_NO_covars %>%
  group_by(`Gene Symbol`) %>%  # Group by gene symbol
  summarise(
    logFC = mean(logFC, na.rm = TRUE),  
    SE = max(SE, na.rm = TRUE),
    Study = "GSE46449"
  ) 

#create meta-data for the meta-analysis
metadata_blood_nocovars = bind_rows(GSE124326_NO_covars, GSE18312_NO_covars,
                                    GSE46416_NO_covars, GSE46449_NO_covars,
                                    GSE23848_NO_covars, GSE39653_NO_covars)

metadata_blood_withcovars = bind_rows(GSE124326_WITH_covars, GSE18312_WITH_covars,
                                      GSE23848_WITH_covars, GSE46416_NO_covars,
                                      GSE46449_NO_covars, GSE39653_NO_covars)

write_csv(metadata_blood_nocovars, "metadata_blood_nocovars.csv")
write_csv(metadata_blood_withcovars, "metadata_blood_withcovars.csv")

#run meta-analysis
#Define function to perform meta-analysis

#loop through each gene in gene list
gene_meta_analysis = function(data) {
  unique_genes = unique(data$`Gene Symbol`)
  results <- data.frame(`Gene Symbol` = character(), meta_LFc = numeric(), meta_se = numeric(),
                        meta_pval = numeric(), tau2 = numeric(), I2 = numeric(),
                        H2 = numeric(), Q = numeric(), Q_p = numeric(), 
                        down = numeric(), up = numeric(),
                        stringsAsFactors = FALSE)
  
  for (x in unique_genes) {
    tmp_df_genes = data[data$`Gene Symbol` == x,] # x is iterated in a loop/function
    if (nrow(tmp_df_genes) >= 3) { #Ensure there is enough data for meta-analysis
      tmp_meta_model = rma.uni(yi = logFC, 
                               vi = SE^2, 
                               data = tmp_df_genes, 
                               method = "SJ", 
                               weighted = TRUE)
      results = bind_rows(results, data.frame(
        `Gene Symbol` = x,
        meta_LFc = tmp_meta_model$b, #meta-analyzed logFC
        meta_se = tmp_meta_model$se, #meta-analyzed SE
        meta_pval = tmp_meta_model$pval, #p-value
        tau2 = tmp_meta_model$tau2, #tau squared, between study variance
        I2 = tmp_meta_model$I2, #I^2 heterogeneity
        H2 = tmp_meta_model$H2, #H^2 heterogeneity
        Q = tmp_meta_model$QE, #Cochran's Q statistic
        Q_p = tmp_meta_model$QEp, #Q statistic p value
        down = sum(tmp_df_genes$logFC < 0), #number of studies with down regulation
        up = sum(tmp_df_genes$logFC >0) #number of studies with up regulation
      )) 
    }
  }
  return(results)
}

#NO COVARIATES
blood_no_covars = gene_meta_analysis(metadata_blood_nocovars)

#add adjusted p-value
blood_no_covars$adj_pval = p.adjust(blood_no_covars$meta_pval, method = "fdr")

#filter only significantly differentially expressed genes
blood_no_covars_sign = blood_no_covars[blood_no_covars$meta_pval < 0.05 & abs(blood_no_covars$meta_LFc) > 0.3,]

write_csv(blood_no_covars, "blood_no_covars.csv")
write_csv(blood_no_covars_sign, "blood_no_covars_sign.csv")


#WITH COVARIATES
blood_with_covars = gene_meta_analysis(metadata_blood_withcovars)

#add adjusted p-values
blood_with_covars$adj_pval = p.adjust(blood_with_covars$meta_pval, method = "fdr")

blood_with_covars_sign = blood_with_covars[blood_with_covars$meta_pval<0.05 & abs(blood_with_covars$meta_LFc)>0.3,]

write_csv(blood_with_covars, "blood_with_covars.csv")
write_csv(blood_with_covars_sign, "blood_with_covars_sign.csv")

#Venn diagrams for the overlap of significantly differentially expressed genes between two models
library("ggVennDiagram")
library(ggplot2)

venn_blood = ggVennDiagram(x = list(blood_no_covars_sign$Gene.Symbol, 
                       blood_with_covars_sign$Gene.Symbol)[1:2], label_alpha = 0,
              category.names = c("no covariates", "with covariates")) +
  ggplot2::scale_fill_gradient(low="white",high = "orange2") +
  coord_flip() +
  ggtitle("Blood") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        aspect.ratio = 0.7)

#combine with brain datasets diagram
venn_brain = ggVennDiagram(x = list(meta_analysis_no_covars_sign$Gene.Symbol, 
                                    meta_analysis_with_covars_significant$Gene.Symbol)[1:2], label_alpha = 0,
                           category.names = c("no covariates", "with covariates")) +
  ggplot2::scale_fill_gradient(low="white",high = "blue") +
  coord_flip() +
  ggtitle("Brain") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        aspect.ratio = 0.7)

library(gridExtra)
figure2 = grid.arrange(venn_brain, venn_blood)

library(ggrepel)
#build volcano plots for two models
#no covars
#create columns for expression
blood_no_covars = blood_no_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.3 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.3 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))
top_up = blood_no_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 10)
top_down = blood_no_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 10)

top_genes = bind_rows(top_up, top_down)

plot_no_covars = ggplot(blood_no_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
  geom_point(aes(color = factor(expression, 
                                levels = c("None", "Down", "Up"))), alpha=0.8, size=2) +
  scale_color_manual(values=c("grey","cyan", "red"))

plot_no_covars = plot_no_covars +
  geom_text_repel(data = top_genes,
                  aes(label = Gene.Symbol),  
                  size = 5,
                  box.padding = 0.35,
                  point.padding = 0.3,
                  segment.color = 'black')
#add p-value <0.05 line
plot_no_covars = plot_no_covars +geom_hline(yintercept = -log10(0.05), linetype="dashed",
                                                color="darkgrey", linewidth=0.8)
#add p-value <0.001 line
plot_no_covars = plot_no_covars +geom_hline(yintercept = -log10(0.001), linetype="dashed",
                                                color="darkgrey", linewidth=0.8)
#add logFC line
if (min(blood_no_covars$meta_LFc) < -0.3){
  plot_no_covars = plot_no_covars + geom_vline(xintercept = -0.3, linetype="dashed",
                                                   color="darkgrey", linewidth=0.8)
}

if (max(blood_no_covars$meta_LFc) > 0.3){
  plot_no_covars = plot_no_covars + geom_vline(xintercept = 0.3, linetype="dashed",
                                                   color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_no_covars = plot_no_covars + 
  labs(title="No covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 15)


#with covars
blood_with_covars = blood_with_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.3 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.3 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))
top_up = blood_with_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 5)
top_down = blood_with_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 5)

top_genes = bind_rows(top_up, top_down)

plot_with_covars = ggplot(blood_with_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
  geom_point(aes(color = factor(expression, 
                                levels = c("None", "Down", "Up"))), alpha=0.8, size=2) +
  scale_color_manual(values=c("grey","dodgerblue3", "red"))

plot_with_covars = plot_with_covars +
  geom_text_repel(data = top_genes,
                  aes(label = Gene.Symbol),  
                  size = 5,
                  box.padding = 0.35,
                  point.padding = 0.3,
                  segment.color = 'black')
#add p-value <0.05 line
plot_with_covars = plot_with_covars +geom_hline(yintercept = -log10(0.05), linetype="dashed",
                                            color="darkgrey", linewidth=0.8)
#add p-value <0.001 line
plot_with_covars = plot_with_covars +geom_hline(yintercept = -log10(0.001), linetype="dashed",
                                            color="darkgrey", linewidth=0.8)
#add logFC line
if (min(blood_with_covars$meta_LFc) < -0.3){
  plot_with_covars = plot_with_covars + geom_vline(xintercept = -0.3, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

if (max(blood_with_covars$meta_LFc) > 0.3){
  plot_with_covars = plot_with_covars + geom_vline(xintercept = 0.3, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_with_covars = plot_with_covars + 
  labs(title="With covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 15)

#combine two plots
#create a function to arrange both plots  
library(grid)
grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }
volcano_blood = grid_arrange_shared_legend(plot_no_covars, plot_with_covars)
ggsave("volcano_blood.pdf", plot = volcano_blood, width = 30, height = 18, units = "cm", useDingbats = FALSE)



#make FOREST PLOTS for the top-6 up and down regulated genes in both models

#create a dataframe with number of BPD and control individuals per dataset
sample_size = data.frame(study = c(
  "GSE124312", "GSE18312", "GSE46416", 
  "GSE46449", "GSE23848", "GSE39653"),
  BD = c(240, 9, 11, 28, 20, 8),
  CTRL = c(240, 8, 10, 25, 15, 24))

#add sample size to metadata
metadata_blood_nocovars = inner_join(metadata_blood_nocovars, sample_size, by = c("Study" = "study"))
metadata_blood_withcovars = inner_join(metadata_blood_withcovars, sample_size, by = c("Study" = "study"))

#NO COVARS UP
#filter the data and perform meta-analysis per gene

#NDUFS5
#filter the data and perform meta-analysis per gene
ndufs5 = metadata_blood_nocovars %>%
  filter(grepl("^NDUFS5$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
ndufs5$sample_size = paste0(ndufs5$BD, "/", ndufs5$CTRL)
#perform meta-analysis
res1 = rma.uni(yi = logFC, sei = SE, data = ndufs5, method = "SJ")
ndufs5$weight = round(weights(res1),1)

"HP"
#filter the data and perform meta-analysis per gene
hp = metadata_blood_nocovars %>%
  filter(grepl("^HP$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
hp$sample_size = paste0(hp$BD, "/", hp$CTRL)
#perform meta-analysis
res2 = rma.uni(yi = logFC, sei = SE, data = hp, method = "SJ")
hp$weight = round(weights(res2),1)


#filter the data and perform meta-analysis per gene
snora10 = metadata_blood_nocovars %>%
  filter(grepl("^SNORA10$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
snora10$sample_size = paste0(snora10$BD, "/", snora10$CTRL)
#perform meta-analysis
res3 = rma.uni(yi = logFC, sei = SE, data = snora10, method = "SJ")
snora10$weight = round(weights(res3),1)

anxa3 = metadata_blood_nocovars %>%
  filter(grepl("^ANXA3$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
anxa3$sample_size = paste0(anxa3$BD, "/", anxa3$CTRL)
#perform meta-analysis
res4 = rma.uni(yi = logFC, sei = SE, data = anxa3, method = "SJ")
anxa3$weight = round(weights(res4),1)

ifit1 = metadata_blood_nocovars %>%
  filter(grepl("^IFIT1$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
ifit1$sample_size = paste0(ifit1$BD, "/", ifit1$CTRL)
#perform meta-analysis
res5 = rma.uni(yi = logFC, sei = SE, data = ifit1, method = "SJ")
ifit1$weight = round(weights(res5),1)

s100p = metadata_blood_nocovars %>%
  filter(grepl("^S100P$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
s100p$sample_size = paste0(s100p$BD, "/", s100p$CTRL)
#perform meta-analysis
res6 = rma.uni(yi = logFC, sei = SE, data = s100p, method = "SJ")
s100p$weight = round(weights(res6),1)

#make a forest plot
png("combined_forest_top6_UP_blood_nocovars.png", width = 2200, height = 2400, res = 150)
par(mfrow = c(3, 2), mar = c(5, 10, 4, 2), oma = c(1, 1, 2, 1))  # for 6 plots in 3 rows, 2 columns

forest_hp = forest(res2, 
                             slab = hp$Study,
                             ilab = cbind(hp$sample_size, hp$weight),
                             ilab.xpos = c(-1.6, -1.2),
                             xlab = "log2 fold change", 
                             main = "HP, Haptoglobin", 
                             cex = 1)
text(-1.6, nrow(hp) + 2, "BD/CTRL", font = 2)
text(-1.2, nrow(hp) + 2, "Weight, %", font = 2)

forest_ndufs5 = forest(res1, 
                         slab = ndufs5$Study,
                         ilab = cbind(ndufs5$sample_size, ndufs5$weight),
                         ilab.xpos = c(-1.7, -1.2),
                         xlab = "log2 fold change", 
                         main = "NDUFS5, NADH:Ubiquinone Oxidoreductase Subunit", 
                         cex = 1)
text(-1.7, nrow(ndufs5) + 2, "BD/CTRL", font = 2)
text(-1.2, nrow(ndufs5) + 2, "Weight, %", font = 2)

forest_snora10 = forest(res3, 
                     slab = snora10$Study,
                     ilab = cbind(snora10$sample_size, snora10$weight),
                     ilab.xpos = c(-0.9, -0.6),
                     xlab = "log2 fold change", 
                     main = "SNORA10, Small Nucleolar RNA", 
                     cex = 1)
text(-0.9, nrow(snora10) + 2, "BD/CTRL", font = 2)
text(-0.6, nrow(snora10) + 2, "Weight, %", font = 2)

forest_anxa3 = forest(res4, 
                             slab = anxa3$Study,
                             ilab = cbind(anxa3$sample_size, anxa3$weight),
                             ilab.xpos = c(-1.5, -1.0),
                             xlab = "log2 fold change", 
                             main = "ANXA3, Annexin A3", 
                             cex = 1)
text(-1.5, nrow(anxa3) + 2, "BD/CTRL", font = 2)
text(-1.0, nrow(anxa3) + 2, "Weight, %", font = 2)

forest_ifit = forest(res5, 
                      slab = ifit1$Study,
                      ilab = cbind(ifit1$sample_size, ifit1$weight),
                      ilab.xpos = c(-2.4, -1.7),
                      xlab = "log2 fold change", 
                      main = "IFIT1, Interferon Induced Protein", 
                      cex = 1)
text(-2.4, nrow(ifit1) + 2, "BD/CTRL", font = 2)
text(-1.7, nrow(ifit1) + 2, "Weight, %", font = 2)

forest_s100p = forest(res6, 
                             slab = s100p$Study,
                             ilab = cbind(s100p$sample_size, s100p$weight),
                             ilab.xpos = c(-2.1, -1.6),
                             xlab = "log2 fold change", 
                             main = "S100P, Calcium-binding protein", 
                             cex = 1)
text(-2.1, nrow(s100p) + 2, "BD/CTRL", font = 2)
text(-1.6, nrow(s100p) + 2, "Weight, %", font = 2)

dev.off()

#TOP DOWN_REGULATED
ccz1b = metadata_blood_nocovars %>%
  filter(grepl("^CCZ1B$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
ccz1b$sample_size = paste0(ccz1b$BD, "/", ccz1b$CTRL)
#perform meta-analysis
res1 = rma.uni(yi = logFC, sei = SE, data = ccz1b, method = "SJ")
ccz1b$weight = round(weights(res1),1)

erv3 = metadata_blood_nocovars %>%
  filter(grepl("^ERV3-1$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
erv3$sample_size = paste0(erv3$BD, "/", erv3$CTRL)
#perform meta-analysis
res2 = rma.uni(yi = logFC, sei = SE, data = erv3, method = "SJ")
erv3$weight = round(weights(res2),1)

smg1 = metadata_blood_nocovars %>%
  filter(grepl("^SMG1$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
smg1$sample_size = paste0(smg1$BD, "/", smg1$CTRL)
#perform meta-analysis
res3 = rma.uni(yi = logFC, sei = SE, data = smg1, method = "SJ")
smg1$weight = round(weights(res3),1)

fam102b = metadata_blood_nocovars %>%
  filter(grepl("^FAM102B$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
fam102b$sample_size = paste0(fam102b$BD, "/", fam102b$CTRL)
#perform meta-analysis
res4 = rma.uni(yi = logFC, sei = SE, data = fam102b, method = "SJ")
fam102b$weight = round(weights(res4),1)

golga8a = metadata_blood_nocovars %>%
  filter(grepl("^GOLGA8A$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
golga8a$sample_size = paste0(golga8a$BD, "/", golga8a$CTRL)
#perform meta-analysis
res5 = rma.uni(yi = logFC, sei = SE, data = golga8a, method = "SJ")
golga8a$weight = round(weights(res5),1)

ptgdr = metadata_blood_nocovars %>%
  filter(grepl("^PTGDR$", ignore.case = T, metadata_blood_nocovars$`Gene Symbol`))
#add sample size
ptgdr$sample_size = paste0(ptgdr$BD, "/", ptgdr$CTRL)
#perform meta-analysis
res6 = rma.uni(yi = logFC, sei = SE, data = ptgdr, method = "SJ")
ptgdr$weight = round(weights(res6),1)


#make a forest plot
png("combined_forest_top6_DOWN_blood_nocovars.png", width = 2200, height = 2400, res = 150)
par(mfrow = c(3, 2), mar = c(5, 10, 4, 2), oma = c(1, 1, 2, 1))  # for 6 plots in 3 rows, 2 columns

forest_ccz1b = forest(res1, 
                   slab = ccz1b$Study,
                   ilab = cbind(ccz1b$sample_size, ccz1b$weight),
                   ilab.xpos = c(-3.1, -2.5),
                   xlab = "log2 fold change", 
                   main = "CCZ1B, Vacuolar protein", 
                   cex = 1)
text(-3.1, nrow(ccz1b) + 2, "BD/CTRL", font = 2)
text(-2.5, nrow(ccz1b) + 2, "Weight, %", font = 2)

forest_erv3 = forest(res2, 
                       slab = erv3$Study,
                       ilab = cbind(erv3$sample_size, erv3$weight),
                       ilab.xpos = c(-3.0, -2.4),
                       xlab = "log2 fold change", 
                       main = "ERV3-1, Endogenous envelope protein", 
                       cex = 1)
text(-3.0, nrow(erv3) + 2, "BD/CTRL", font = 2)
text(-2.4, nrow(erv3) + 2, "Weight, %", font = 2)

forest_smg1 = forest(res3, 
                        slab = smg1$Study,
                        ilab = cbind(smg1$sample_size, smg1$weight),
                        ilab.xpos = c(-2.7, -2.2),
                        xlab = "log2 fold change", 
                        main = "SMG1, Serine/threonine protein kinase", 
                        cex = 1)
text(-2.7, nrow(smg1) + 2, "BD/CTRL", font = 2)
text(-2.2, nrow(smg1) + 2, "Weight, %", font = 2)

forest_fam102b = forest(res4, 
                      slab = fam102b$Study,
                      ilab = cbind(fam102b$sample_size, fam102b$weight),
                      ilab.xpos = c(-1.8, -1.4),
                      xlab = "log2 fold change", 
                      main = "FAM102B, EEIG family member", 
                      cex = 1)
text(-1.8, nrow(fam102b) + 2, "BD/CTRL", font = 2)
text(-1.4, nrow(fam102b) + 2, "Weight, %", font = 2)

forest_golga8a = forest(res5, 
                     slab = golga8a$Study,
                     ilab = cbind(golga8a$sample_size, golga8a$weight),
                     ilab.xpos = c(-3.3, -2.7),
                     xlab = "log2 fold change", 
                     main = "GOLGA8A, Golgi apparatus protein", 
                     cex = 1)
text(-3.3, nrow(golga8a) + 2, "BD/CTRL", font = 2)
text(-2.7, nrow(golga8a) + 2, "Weight, %", font = 2)

forest_ptgdr = forest(res6, 
                      slab = ptgdr$Study,
                      ilab = cbind(ptgdr$sample_size, ptgdr$weight),
                      ilab.xpos = c(-3.2, -2.6),
                      xlab = "log2 fold change", 
                      main = "PTGDR, Prostaglandin D2 Receptor", 
                      cex = 1)
text(-3.2, nrow(ptgdr) + 2, "BD/CTRL", font = 2)
text(-2.6, nrow(ptgdr) + 2, "Weight, %", font = 2)

dev.off()


#ENRICHMENT ANALYSIS
library(org.Hs.eg.db)
library(clusterProfiler)
library(hgu133plus2cdf)
library(AnnotationDbi)
library(purrr)
#no covars
#create gene universe (unique genes in metadata)
gene_universe = unique(metadata_blood_nocovars$`Gene Symbol`)
head(gene_universe) 
#gene symbols are gene names, cluster profiler needs ENTREZ ID
#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                        keys=gene_universe, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create list of dif expr genes
genes_de = blood_no_covars[blood_no_covars$meta_pval<0.05,]$Gene.Symbol
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

enrichment = setReadable(enrichment, "org.Hs.eg.db", "ENTREZID")
barplot(enrichment, showCategory = 10, font.size=15)
write_csv(enrichment@result, "enrichment_blood_NO_covars.csv")

#KEGG pathway analysis
kegg = enrichKEGG(gene = genes_de, 
                  organism = "hsa",
                  pvalueCutoff = 0.05)
dotplot(kegg, showCategory = 10)

#with covariates
gene_universe = unique(metadata_blood_withcovars$`Gene Symbol`)

#convert gene symbol to ENTREZ ID
genes_entrezid = mapIds(org.Hs.eg.db, 
                        keys=gene_universe, 
                        column="ENTREZID", 
                        keytype="SYMBOL")
#create list of dif expr genes
genes_de = blood_with_covars[blood_with_covars$meta_pval<0.05,]$Gene.Symbol
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

enrichment = setReadable(enrichment, "org.Hs.eg.db", "ENTREZID")
barplot(enrichment, showCategory = 10, font.size=15)
write_csv(enrichment@result, "enrichment_blood_WITH_covars.csv")

#KEGG pathway analysis
kegg = enrichKEGG(gene = genes_de, 
                  organism = "hsa",
                  pvalueCutoff = 0.05)
dotplot(kegg, showCategory = 10)

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
library(ggtree)
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
tree = pairwise_termsim(enrichment)

upsetplot(enrichment)

#correlation between metaLFc in brain and blood samples
blood = blood_no_covars[blood_no_covars$meta_pval<0.05, c(1,2)]
brain = meta_analysis_no_covars[meta_analysis_no_covars$meta_pval<0.05,c(1,2)]
colnames(blood) = c("Gene.Symbol", "meta_logFC_blood")
colnames(brain) = c("Gene.Symbol", "meta_logFC_brain")
df = merge(blood, brain, by = "Gene.Symbol")

library(ggrepel)
ggplot(df, aes(x = meta_logFC_blood, y = meta_logFC_brain)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red3", se = FALSE, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = Gene.Symbol), size = 3, max.overlaps = 50, color = alpha("black", 0.7)) + 
  theme_minimal() +
  labs(x = "Meta logFC (Blood)", y = "Meta logFC (Brain)", 
       caption = paste0("Spearman r = ", round(cor_result$estimate, 2), 
                        ", p = ", signif(cor_result$p.value, 3)))

cor_result = cor.test(df$meta_logFC_brain, df$meta_logFC_blood, method = "spearman")












