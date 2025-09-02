rm(list=ls())
library(dplyr)
library(purrr)
BiocManager::install("metafor")
library(metafor)
library(ggrepel)
library(gridExtra)
library(grid)

#clean the tables and add Study column
#GSE5392
GSE5392_no_covars = GSE5392_pfc_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE5392_no_covars$Study = "GSE5392"
write_csv(GSE5392_no_covars, "GSE5392_no_covars.csv")

GSE5392_with_covars = GSE5392_pfc_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE5392_with_covars$Study = "GSE5392"
write_csv(GSE5392_with_covars, "GSE5392_with_covars.csv")

#GSE12649 
GSE12649_no_covars = GSE12649_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE12649_no_covars$Study = "GSE12649"
write_csv(GSE12649_no_covars, "GSE12649_no_covars.csv")

#GSE12679
GSE12679_no_covars = GSE12679_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE12679_no_covars$Study = "GSE12679"
write_csv(GSE12679_no_covars, "GSE12679_no_covars.csv")

GSE12679_with_covars = GSE12679_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE12679_with_covars$Study = "GSE12679"
write_csv(GSE12679_with_covars, "GSE12679_with_covars.csv")

#GSE35978
GSE35978_no_covars = GSE35978_parietalcortex_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE35978_no_covars$Study = "GSE35978"
write_csv(GSE35978_no_covars, "GSE35978_no_covars.csv")

GSE35978_with_covars = GSE35978_parietalcortex_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE35978_with_covars$Study = "GSE35978"
write_csv(GSE35978_with_covars, "GSE35978_with_covars.csv")

#GSE53987
GSE53987_no_covars = GSE53987_pfc_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE53987_no_covars$Study = "GSE53987"
write_csv(GSE53987_no_covars, "GSE53987_no_covars.csv")

GSE53987_with_covars = GSE53987_pfc_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE53987_with_covars$Study = "GSE53987"
write_csv(GSE53987_with_covars, "GSE53987_with_covars.csv")

#GSE62191
GSE62191_no_covars = GSE62191_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE62191_no_covars$Study = "GSE62191"
write_csv(GSE62191_no_covars, "GSE62191_no_covars.csv")

GSE62191_with_covars = GSE62191_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE62191_with_covars$Study = "GSE62191"
write_csv(GSE62191_with_covars, "GSE62191_with_covars.csv")

#GSE78246
GSE78246_no_covars = GSE78246_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE78246_no_covars$Study = "GSE78246"
write_csv(GSE78246_no_covars, "GSE78246_no_covars.csv")


GSE78246_with_covars = GSE78246_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE78246_with_covars$Study = "GSE78246"
write_csv(GSE78246_with_covars, "GSE78246_with_covars.csv")

#GSE87610
GSE87610_no_covars = GSE87610_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE87610_no_covars$Study = "GSE87610"
write_csv(GSE87610_no_covars, "GSE87610_no_covars.csv")

#GSE92538
GSE92538_no_covars = GSE92538_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE92538_no_covars$Study = "GSE92538"
write_csv(GSE92538_no_covars, "GSE92538_no_covars.csv")

GSE92538_with_covars = GSE92538_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE92538_with_covars$Study = "GSE92538"
write_csv(GSE92538_with_covars, "GSE92538_with_covars.csv")



#GSE120340
GSE120340_no_covars = GSE120340_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE120340_no_covars$Study = "GSE120340"
write_csv(GSE120340_no_covars, "GSE120340_no_covars.csv")

GSE120340_with_covars = GSE120340_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE120340_with_covars$Study = "GSE120340"
write_csv(GSE120340_with_covars, "GSE120340_with_covars.csv")

#GSE208338
GSE208338_no_covars = GSE208338_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE208338_no_covars$Study = "GSE202338"
write_csv(GSE208338_no_covars, "GSE208338_no_covars.csv")

GSE208338_with_covars = GSE208338_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE208338_with_covars$Study = "GSE202338"
write_csv(GSE208338_with_covars, "GSE208338_with_covars.csv")

#GSE210064
GSE210064_no_covars = GSE210064_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE210064_no_covars$Study = "GSE210064"
write_csv(GSE210064_no_covars, "GSE210064_no_covars.csv")

GSE210064_with_covars = GSE210064_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE210064_with_covars$Study = "GSE210064"
write_csv(GSE210064_with_covars, "GSE210064_with_covars.csv")

#GSE42546
GSE42546_no_covars = GSE42546_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE42546_no_covars$Study = "GSE42546"
write_csv(GSE42546_no_covars, "GSE42546_no_covars.csv")

GSE42546_with_covars = GSE42546_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE42546_with_covars$Study = "GSE42546"
write_csv(GSE42546_with_covars, "GSE42546_with_covars.csv")

#GSE53239
GSE53239_no_covars = GSE53239_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE53239_no_covars$Study = "GSE53239"
write_csv(GSE53239_no_covars, "GSE53239_no_covars.csv")

#GSE78936
GSE78936_no_covars = GSE78936_ba11_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE78936_no_covars$Study = "GSE78936"
write_csv(GSE78936_no_covars, "GSE78936_no_covars.csv")

#GSE80336
GSE80336_no_covars = GSE80336_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE80336_no_covars$Study = "GSE80336"
write_csv(GSE80336_no_covars, "GSE80336_no_covars.csv")

GSE80336_with_covars = GSE80336_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE80336_with_covars$Study = "GSE80336"
write_csv(GSE80336_with_covars, "GSE80336_with_covars.csv")

#GSE80655
GSE80655_no_covars = GSE80655_pfc_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE80655_no_covars$Study = "GSE80655"
write_csv(GSE80655_no_covars, "GSE80655_no_covars.csv")

GSE80655_with_covars = GSE80655_pfc_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE80655_with_covars$Study = "GSE80655"
write_csv(GSE80655_with_covars, "GSE80655_with_covars.csv")

#GSE81396
GSE81396_no_covars = GSE81396_putamen_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE81396_no_covars$Study = "GSE81396"
write_csv(GSE81396_no_covars, "GSE81396_no_covars.csv")

#GSE202537
GSE202537_no_covars = GSE202537_putamen_no_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE202537_no_covars$Study = "GSE202537"
write_csv(GSE202537_no_covars, "GSE202537_no_covars.csv")

GSE202537_with_covars = GSE202537_putamen_with_covars_curated[, c("Gene Symbol", "logFC", "SE")]
GSE202537_with_covars$Study = "GSE202537"
write_csv(GSE202537_with_covars, "GSE202537_with_covars.csv")

#create lists with data for metaanalysis
no_covars = c(GSE5392_no_covars, GSE12649_no_covars, GSE12679_no_covars,
              GSE35978_no_covars, GSE42546_no_covars, GSE53239_no_covars,
              GSE53987_no_covars, GSE62191_no_covars, GSE78246_no_covars,
              GSE78936_no_covars, GSE80336_no_covars, GSE80655_no_covars,
              GSE81396_no_covars, GSE87610_no_covars, GSE92538_no_covars,
              GSE120340_no_covars, GSE202537_no_covars, GSE208338_no_covars,
              GSE210064_no_covars)
with_covars = c(GSE5392_with_covars, GSE12649_no_covars, GSE12679_with_covars,
                GSE35978_with_covars, GSE42546_with_covars, GSE53239_no_covars,
                GSE53987_with_covars, GSE62191_with_covars, GSE78246_with_covars,
                GSE78936_no_covars, GSE80336_with_covars, GSE80655_with_covars,
                GSE81396_no_covars, GSE87610_no_covars, GSE92538_with_covars,
                GSE120340_with_covars, GSE202537_with_covars, GSE208338_with_covars,
                GSE210064_with_covars)

#merge data into one table
#without covars
metadata_no_covars = bind_rows(GSE5392_no_covars, GSE12649_no_covars, GSE12679_no_covars,
                               GSE35978_no_covars, GSE42546_no_covars, GSE53239_no_covars,
                               GSE53987_no_covars, GSE62191_no_covars, GSE78246_no_covars,
                               GSE78936_no_covars, GSE80336_no_covars, GSE80655_no_covars,
                               GSE81396_no_covars, GSE87610_no_covars, GSE92538_no_covars,
                               GSE120340_no_covars, GSE202537_no_covars, GSE208338_no_covars,
                               GSE210064_no_covars)
#with covars
metadata_with_covars = bind_rows(GSE5392_with_covars, GSE12649_no_covars, GSE12679_with_covars,
                                 GSE35978_with_covars, GSE42546_with_covars, GSE53239_no_covars,
                                 GSE53987_with_covars, GSE62191_with_covars, GSE78246_with_covars,
                                 GSE78936_no_covars, GSE80336_with_covars, GSE80655_with_covars,
                                 GSE81396_no_covars, GSE87610_no_covars, GSE92538_with_covars,
                                 GSE120340_with_covars, GSE202537_with_covars, GSE208338_with_covars,
                                 GSE210064_with_covars)
#extract unique gene symbols
genes1 = unique(metadata_no_covars$`Gene Symbol`)
genes2 = unique(metadata_with_covars$`Gene Symbol`)

gene_list = reduce(list(genes1, genes2), union) #list of unique gene symbols

#calculate number of up/down regulated genes
gene_counts_per_study = metadata_no_covars %>%
  group_by(Study) %>%
  summarise(
    Upregulated = sum(logFC > 0, na.rm = TRUE),
    Downregulated = sum(logFC < 0, na.rm = TRUE),
    Total_Genes = n() #how many genes total per cohort
  )
gene_counts_per_study = gene_counts_per_study %>%
  mutate(Upregulated_Percent = (Upregulated / Total_Genes) * 100,
         Downregulated_Percent = (Downregulated / Total_Genes) * 100)


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
    if (nrow(tmp_df_genes) >= 5) { #Ensure there is enough data for meta-analysis
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
        up = sum(tmp_df_genes$logFC >0) #number of studies withup regulation
        )) 
  }
}
return(results)
}


###META-ANALYSIS A
##Full meta-analysis no covars
meta_analysis_no_covars_copy = meta_analysis_no_covars 
meta_analysis_no_covars = gene_meta_analysis(metadata_no_covars)
#add adjusted p-values
meta_analysis_no_covars$adj_pval = p.adjust(meta_analysis_no_covars$meta_pval, method = "fdr")

meta_analysis_no_covars_significant = meta_analysis_no_covars[meta_analysis_no_covars$meta_pval < 0.05 
                                                             & abs(meta_analysis_no_covars$meta_LFc) >0.3, ]
nrow(meta_analysis_no_covars_significant[meta_analysis_no_covars_significant$meta_LFc<0,])

#Full meta-analysis with covars 
meta_analysis_with_covars_copy = meta_analysis_with_covars
meta_analysis_with_covars = gene_meta_analysis(metadata_with_covars)

#add adjusted meta p-value
meta_analysis_with_covars$adj_pval = p.adjust(meta_analysis_with_covars$meta_pval, method = "fdr")

meta_analysis_with_covars_significant = meta_analysis_with_covars[meta_analysis_with_covars$meta_pval < 0.05
                                                                  & abs(meta_analysis_with_covars$meta_LFc) >0.25,]


#save the results
write_csv(meta_analysis_no_covars_significant, "meta-analysis_no_covars_sign.csv")
write_csv(meta_analysis_with_covars_significant, "meta-analysis_with_covars_significant.csv")
write_csv(meta_analysis_no_covars, "meta_analysis_no_covars.csv")
write_csv(meta_analysis_with_covars, "meta_analysis_with_covars.csv")
write_csv(metadata_no_covars, "metadata_no_covars.csv")
write_csv(metadata_with_covars, "metadata_with_covars.csv")


####sub meta-analysis, META-ANALYSIS B
##meta-analysis BA9, dorsolateral prefrontal cortex
#list the studies
ba9_list_no_covars = bind_rows(GSE87610_no_covars, GSE12679_no_covars, GSE5392_no_covars,
                       GSE92538_no_covars, GSE80655_no_covars, GSE78936_no_covars,
                       GSE208338_no_covars)

ba9_list_with_covars = bind_rows(GSE87610_no_covars, GSE12679_with_covars, GSE5392_with_covars,
                                 GSE92538_with_covars, GSE80655_with_covars, GSE78936_no_covars,
                                 GSE208338_with_covars)
write_csv(ba9_list_no_covars, "ba9_metadata_no_covars.csv")
write_csv(ba9_list_with_covars, "ba9_metadata_with_covars.csv")

#no covars
ba9_no_covars = gene_meta_analysis(ba9_list_no_covars)

#add adjusted p-value
ba9_no_covars$adj_pval = p.adjust(ba9_no_covars$meta_pval, method = 'fdr')
ba9_no_covars_significant = ba9_no_covars[ba9_no_covars$meta_pval <0.05
                                                                  &abs(ba9_no_covars$meta_LFc)>0.2,]
write_csv(ba9_no_covars_significant, "ba9_no_covars_sign.csv")
write_csv(ba9_no_covars, "ba9_no_covars.csv")


#with covars
metaanalysis_ba9_with_covars = gene_meta_analysis(ba9_list_with_covars)
nrow(meta_analysis_with_covars[meta_analysis_with_covars$meta_pval<0.05,])
#add adjusted p-value
ba9_with_covars$adj_pval = p.adjust(ba9_with_covars$meta_pval, method = 'fdr')
ba9_with_covars_significant = ba9_with_covars[ba9_with_covars$meta_pval<0.05
                                                                        &abs(ba9_with_covars$meta_LFc)>0.2,]

write_csv(ba9_with_covars_significant, "ba9_with_covars_sign.csv")
write_csv(ba9_with_covars, "ba9_with_covars.csv")


#create volcano plots for full meta-analysis
BiocManager::install("ggplot2")
library(ggplot2)
library(ggrepel)

###VOLCANO PLOTS FOR META-ANALYSIS A
##without covariates
#create column for expression
meta_analysis_no_covars = meta_analysis_no_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.3 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.3 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))

top_up = meta_analysis_no_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 5)
top_down = meta_analysis_no_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 9)

top_genes = bind_rows(top_up, top_down)


plot_no_covars = ggplot(meta_analysis_no_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
                          geom_point(aes(color = factor(expression, 
                                                        levels = c("None", "Down", "Up"))), alpha=0.8, size=2) +
                          scale_color_manual(values=c("grey","dodgerblue3", "red"))

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
if (min(meta_analysis_no_covars$meta_LFc) < -0.3){
  plot_no_covars = plot_no_covars + geom_vline(xintercept = -0.3, linetype="dashed",
                           color="darkgrey", linewidth=0.8)
}

if (max(meta_analysis_no_covars$meta_LFc) > 0.3){
  plot_no_covars = plot_no_covars + geom_vline(xintercept = 0.3, linetype="dashed",
                           color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_no_covars = plot_no_covars +
  labs(title="No covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 20)
plot_no_covars

##with covariates
#create columns for expression
meta_analysis_with_covars = meta_analysis_with_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.3 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.3 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))
top_up = meta_analysis_with_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 5)
top_down = meta_analysis_with_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 5)

top_genes = bind_rows(top_up, top_down)

plot_with_covars = ggplot(meta_analysis_with_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
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
if (min(meta_analysis_with_covars$meta_LFc) < -0.3){
  plot_with_covars = plot_with_covars + geom_vline(xintercept = -0.3, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

if (max(meta_analysis_with_covars$meta_LFc) > 0.3){
  plot_with_covars = plot_with_covars + geom_vline(xintercept = 0.3, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_with_covars = plot_with_covars + 
  labs(title="With covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 20)
plot_with_covars

#create a function to arrange both plots   
library(gridExtra)
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
volcano = grid_arrange_shared_legend(plot_no_covars, plot_with_covars)
ggsave("volcano_brain.pdf", plot = volcano, width = 30, height = 18, units = "cm", useDingbats = FALSE)
library(grid)

###volcano plots for Meta-analysis B (BA9)
##without covariates
#create column for expression
ba9_no_covars = ba9_no_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.2 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.2 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))
top_up = ba9_no_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 5)
top_down = ba9_no_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 5)

top_genes = bind_rows(top_up, top_down)

plot_ba9_no_covars = ggplot(ba9_no_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
  geom_point(aes(color = factor(expression, 
                                levels = c("None", "Down", "Up"))), alpha=0.8, size=2) +
  scale_color_manual(values=c("grey","dodgerblue3", "red"))

plot_ba9_no_covars = plot_ba9_no_covars +
  geom_text_repel(data = top_genes,
                  aes(label = Gene.Symbol),  
                  size = 5,
                  box.padding = 0.35,
                  point.padding = 0.3,
                  segment.color = 'black')

#add p-value <0.05 line
plot_ba9_no_covars = plot_ba9_no_covars +geom_hline(yintercept = -log10(0.05), linetype="dashed",
                                            color="darkgrey", linewidth=0.8)
#add p-value <0.001 line
plot_ba9_no_covars = plot_ba9_no_covars +geom_hline(yintercept = -log10(0.001), linetype="dashed",
                                            color="darkgrey", linewidth=0.8)
#add logFC line
if (min(ba9_no_covars$meta_LFc) < -0.2){
  plot_ba9_no_covars = plot_ba9_no_covars + geom_vline(xintercept = -0.2, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

if (max(ba9_no_covars$meta_LFc) > 0.2){
  plot_ba9_no_covars = plot_ba9_no_covars + geom_vline(xintercept = 0.2, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_ba9_no_covars = plot_ba9_no_covars + 
  labs(title="No covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 15)



##with covariates
#create columns for expression
ba9_with_covars = ba9_with_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.2 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.2 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))
top_up = ba9_with_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 5)
top_down = ba9_with_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 5)

top_genes = bind_rows(top_up, top_down)

plot_ba9_with_covars = ggplot(ba9_with_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
  geom_point(aes(color = factor(expression, 
                                levels = c("None", "Down", "Up"))), alpha=0.8, size=2) +
  scale_color_manual(values=c("grey","dodgerblue3", "red"))

plot_ba9_with_covars = plot_ba9_with_covars +
  geom_text_repel(data = top_genes,
                  aes(label = Gene.Symbol),  
                  size = 5,
                  box.padding = 0.35,
                  point.padding = 0.3,
                  segment.color = 'black')

#add p-value <0.05 line
plot_ba9_with_covars = plot_ba9_with_covars +geom_hline(yintercept = -log10(0.05), linetype="dashed",
                                                    color="darkgrey", linewidth=0.8)
#add p-value <0.001 line
plot_ba9_with_covars = plot_ba9_with_covars +geom_hline(yintercept = -log10(0.001), linetype="dashed",
                                                    color="darkgrey", linewidth=0.8)
#add logFC line
if (min(ba9_with_covars$meta_LFc) < -0.2){
  plot_ba9_with_covars = plot_ba9_with_covars + geom_vline(xintercept = -0.2, linetype="dashed",
                                                       color="darkgrey", linewidth=0.8)
}

if (max(ba9_with_covars$meta_LFc) > 0.2){
  plot_ba9_with_covars = plot_ba9_with_covars + geom_vline(xintercept = 0.2, linetype="dashed",
                                                       color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_ba9_with_covars = plot_ba9_with_covars + 
  labs(title="With covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 15)

volcano_B = grid_arrange_shared_legend(plot_ba9_no_covars, plot_ba9_with_covars)


###Sub META-ANALYSIS C
##dorsolateral pre-frontal cortex: BA9, BA46
#create table with all data
pfc_no_covars = bind_rows(ba9_list_no_covars, GSE120340_no_covars, GSE12649_no_covars, 
                          GSE53987_no_covars, GSE53239_no_covars)
write.csv(pfc_no_covars, "pfc_metadata_nocov.csv")
#run meta-analysis
pfc_no_covars = gene_meta_analysis(pfc_no_covars)
#add adjusted p-values
pfc_no_covars$adj_pval = p.adjust(pfc_no_covars$meta_pval, method = "fdr")
pfc_no_covars_significant = pfc_no_covars[pfc_no_covars$meta_pval<0.05
                                                                  & abs(pfc_no_covars$meta_LFc)>0.25,]
#save the results
write_csv(pfc_no_covars_significant, "pfc_no_covars_sign.csv")
write_csv(pfc_no_covars, "pfc_no_covars.csv")


#create table with all data
pfc_with_covars = bind_rows(ba9_list_with_covars, GSE120340_with_covars, GSE12649_no_covars, 
                          GSE53987_with_covars, GSE53239_no_covars)
write.csv(pfc_with_covars, "pfc_metadata_withcovar.csv")

#run meta-analysis
pfc_with_covars = gene_meta_analysis(pfc_with_covars)
#add adjusted p-value
pfc_with_covars$adj_pval = p.adjust(pfc_with_covars$meta_pval, method = "fdr")

pfc_with_covars_significant = pfc_with_covars[pfc_with_covars$meta_pval<0.05
                                                                      & abs(pfc_with_covars$meta_LFc)>0.25,]
#save the results
write_csv(pfc_with_covars_significant, "pfc_with_covars_sign.csv")
write_csv(pfc_with_covars, "pfc_with_covars.csv")


##VOLCANO PLOT FOR SUB META_ANALYSIS C
#no covars
#create column for expression
pfc_no_covars = pfc_no_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.25 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.25 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))
top_up = pfc_no_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 5)
top_down = pfc_no_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 5)

top_genes = bind_rows(top_up, top_down)

plot_no_covars = ggplot(pfc_no_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
  geom_point(aes(color = factor(expression, 
                                levels = c("None", "Down", "Up"))), alpha=0.8, size=2) +
  scale_color_manual(values=c("grey","dodgerblue3", "red"))

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
if (min(pfc_no_covars$meta_LFc) < -0.25){
  plot_no_covars = plot_no_covars + geom_vline(xintercept = -0.25, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

if (max(pfc_no_covars$meta_LFc) > 0.25){
  plot_no_covars = plot_no_covars + geom_vline(xintercept = 0.25, linetype="dashed",
                                               color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_no_covars = plot_no_covars + 
  labs(title="No covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 15)



##with covariates
#create columns for expression
pfc_with_covars = pfc_with_covars %>%
  mutate(expression = case_when(
    meta_LFc > 0.25 & meta_pval < 0.05 ~ "Up",
    meta_LFc < -0.25 & meta_pval < 0.05 ~ "Down",
    .default = "None"
  ))
top_up = pfc_with_covars %>%
  filter(expression == "Up") %>%
  arrange(desc(meta_LFc)) %>%
  slice_head(n = 5)
top_down = pfc_with_covars %>%
  filter(expression == "Down") %>%
  arrange(meta_LFc) %>%
  slice_head(n = 5)

top_genes = bind_rows(top_up, top_down)

plot_with_covars = ggplot(pfc_with_covars, aes(x=meta_LFc, y=-log10(meta_pval))) +
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
if (min(pfc_with_covars$meta_LFc) < -0.25){
  plot_with_covars = plot_with_covars + geom_vline(xintercept = -0.25, linetype="dashed",
                                                   color="darkgrey", linewidth=0.8)
}

if (max(pfc_with_covars$meta_LFc) > 0.25){
  plot_with_covars = plot_with_covars + geom_vline(xintercept = 0.25, linetype="dashed",
                                                   color="darkgrey", linewidth=0.8)
}

#add axis names 
plot_with_covars = plot_with_covars + 
  labs(title="With covariates", x="log2 fold change", y="-log10 p-value", col="expression") +
  theme_bw(base_size = 15)


#put two plots together
volcano_C = grid_arrange_shared_legend(plot_no_covars, plot_with_covars)


##Venn diagram 
#Full meta-analysis (A)
#No covariates VS with covariates (significant genes)
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

#all tissues covars/no covars sign genes
venn_brain = ggVennDiagram(x = list(meta_analysis_no_covars_sign$Gene.Symbol, 
                                    meta_analysis_with_covars_significant$Gene.Symbol)[1:2], label_alpha = 0,
                           category.names = c("no covariates", "with covariates")) +
  ggplot2::scale_fill_gradient(low="white",high = "palegreen2") +
  coord_flip() +
  ggtitle("Brain") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        aspect.ratio = 0.7) 




#BA9 covars/no covars sign genes
ggVennDiagram(x = list(ba9_no_covars_sign$Gene.Symbol, 
                       ba9_with_covars_sign$Gene.Symbol)[1:2], label_alpha = 0,
              category.names = c("no covariates", "with covariates")) +
  ggplot2::scale_fill_gradient(low="white",high = "tan2") +
  coord_flip()

#BA9+BA46 covars/no covars sign genes
ggVennDiagram(x = list(pfc_no_covars_sign$Gene.Symbol, 
                       pfc_with_covars_sign$Gene.Symbol)[1:2], label_alpha = 0,
              category.names = c("no covariates", "with covariates")) +
  ggplot2::scale_fill_gradient(low="white",high = "steelblue2") +
  coord_flip()


###FOREST PLOTS for top 6 up genes in all tissues (meta-analysis A)
#create a dataframe with number of BPD and control individuals per dataset
sample_size = data.frame(
  study = c("GSE87610", "GSE78246", "GSE35978", 
            "GSE62191", "GSE210064", "GSE120340", 
            "GSE12679", "GSE12649", "GSE53987", 
            "GSE5392", "GSE92538", "GSE80655", 
            "GSE42546", "GSE78936", "GSE202537", 
            "GSE81396", "GSE80336", "GSE208338", 
            "GSE53239", "GSE12654"),
  BD = c(73, 9, 46, 29, 11, 10,
         5, 33, 17, 30, 12,
         23, 16, 16, 8, 4, 18,
        15, 11, 11),
  CTRL = c(72, 11, 50, 30, 11, 10,
           11, 34, 19, 31, 56, 24,
           29, 12, 35, 4, 18, 62,
           11, 15))
#add number of BD and controls
metadata_no_covars = inner_join(metadata_no_covars, sample_size, by=c("Study" = "study"))

###6 top up-regulated genes based on logFC
##no covars
#prepare data
xist = metadata_no_covars %>%
  filter(grepl("XIST", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
xist$sample_size = paste0(xist$BD, "/", xist$CTRL)
#perform meta-analysis
res = rma.uni(yi = logFC, sei = SE, data = xist, method = "SJ")
#add weights to the res
#since we have random effects model, the weights are computed
#using the formula 1/(vi+tau2), where vi = SE^2, and tau2 - estimated between-study variance
xist$weights = round(weights(res), 1)


##EFCAB3P1
efcab3p1 = metadata_no_covars %>%
  filter(grepl("EFCAB3P1", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
efcab3p1$sample_size = paste0(efcab3p1$BD, "/", efcab3p1$CTRL)

#perform meta-analysis
res1 = rma.uni(yi = logFC, sei = SE, data = efcab3p1, method = "SJ")
#add weights
efcab3p1$weights = round(weights(res1), 1)

##BAALC-AS1
baalc_as1 = metadata_no_covars %>%
  filter(grepl("BAALC-AS1", ignore.case = T, metadata_no_covars$`Gene Symbol`))

#add sample size
baalc_as1$sample_size = paste0(baalc_as1$BD, "/", baalc_as1$CTRL)

#perform meta-analysis
res3 = rma.uni(yi = logFC, sei = SE, data = baalc_as1, method = "SJ")
#add weights
baalc_as1$weights = round(weights(res3), 1)

##LOC112268133
loc112268133 = metadata_no_covars %>%
  filter(grepl("LOC112268133", ignore.case = T, metadata_no_covars$`Gene Symbol`))

#add sample size
loc112268133$sample_size = paste0(loc112268133$BD, "/", loc112268133$CTRL)
#perfrom meta-analysis
res4 = rma.uni(yi = logFC, sei = SE, data = loc112268133, method = "SJ")
#add weights
loc112268133$weights = round(weights(res4),1)

#LINC02884
linc02884 = metadata_no_covars %>%
  filter(grepl("LINC02884", ignore.case = T, metadata_no_covars$`Gene Symbol`))

#add sample size
linc02884$sample_size = paste0(linc02884$BD, "/", linc02884$CTRL)
#perfrom meta-analysis
res5 = rma.uni(yi = logFC, sei = SE, data = linc02884, method = "SJ")
#add weights
linc02884$weights = round(weights(res5),1)

#TRNN
trnn = metadata_no_covars %>%
  filter(grepl("^TRNN$", ignore.case = T, metadata_no_covars$`Gene Symbol`))

#add sample size
trnn$sample_size = paste0(trnn$BD, "/", trnn$CTRL)
#perform meta-analysis
res10 = rma.uni(yi = logFC, sei = SE, data = trnn, method = "SJ")
#add weights
trnn$weights = round(weights(res10),1)


#combine forest plots for top-6 up-regulated genes
png("combined_forest_top6_UP_A_nocovars.png", width = 2200, height = 2400, res = 150)
par(mfrow = c(3, 2), mar = c(5, 10, 4, 2), oma = c(1, 1, 2, 1))  # for 6 plots in 3 rows, 2 columns
forest_xist = forest(res, 
                     slab = xist$Study,
                     ilab = cbind(xist$sample_size, xist$weights),
                     ilab.xpos = c(-19, -14.5),
                     xlab = "log2 fold change",
                     main = "XIST, Long non-coding RNA", 
                     cex = 1)
text(-19, nrow(xist) + 2, "BD/CTRL", font = 2)
text(-14.5, nrow(xist) + 2, "Weight, %", font = 2)

forest_efcab3p1 = forest(res1, 
                         slab = efcab3p1$Study, 
                         ilab = cbind(efcab3p1$sample_size, efcab3p1$weights),
                         ilab.xpos = c(-2.2, -1.5),
                         xlab = "log2 fold change", 
                         main = "EFCAB3P1, Pseudogene", 
                         cex = 1)
text(-2.2, nrow(efcab3p1) + 2, "BD/CTRL", font = 2)
text(-1.5, nrow(efcab3p1) + 2, "Weight, %", font = 2)

forest_baalc_as1 = forest(res3, 
                          slab = baalc_as1$Study, 
                          ilab = cbind(baalc_as1$sample_size, baalc_as1$weights),
                          ilab.xpos = c(-3.4, -2.4),
                          xlab = "log2 fold change", 
                          main = "BAALC-AS1, BAALC anti-sense RNA",
                          cex = 1)
text(-3.4, nrow(baalc_as1) + 2, "BD/CTRL", font = 2)
text(-2.4, nrow(baalc_as1) + 2, "Weight, %", font = 2)

forest_loc112268133 = forest(res4, 
                             slab = loc112268133$Study,
                             ilab = cbind(loc112268133$sample_size, loc112268133$weights),
                             ilab.xpos = c(-3.0, -2.2),
                             xlab = "log2 fold change", 
                             main = "LOC112268133, Uncharacterized/Predicted gene",
                             cex = 1)

text(-3.0, nrow(loc112268133) + 2, "BD/CTRL", font = 2)
text(-2.2, nrow(loc112268133) + 2, "Weight, %", font = 2)

forest_linc02884 = forest(res5, 
                          slab = linc02884$Study,
                          ilab = cbind(linc02884$sample_size, linc02884$weights),
                          ilab.xpos = c(-4.3, -3.4),
                          xlab = "log2 fold change", 
                          main = "LINC02884, Long intergenic non-coding RNA",
                          cex = 1)

text(-4.3, nrow(linc02884) + 2, "BD/CTRL", font = 2)
text(-3.4, nrow(linc02884) + 2, "Weight, %", font = 2)

forest_trnn = forest(res10, 
                     slab = trnn$Study,
                     ilab = cbind(trnn$sample_size, trnn$weights),
                     ilab.xpos = c(-5.3, -4.0),
                     xlab = "log2 fold change", 
                     main = "TRNN, tRNA-related gene",
                     cex = 1)

text(-5.3, nrow(trnn) + 2, "BD/CTRL", font = 2)
text(-4.0, nrow(trnn) + 2, "Weight, %", font = 2)

dev.off()

##top 6 down-regulated genes based on meta-estimated logFC 
#RPS4Y1
rps4y1 = metadata_no_covars %>%
  filter(grepl("^RPS4Y1", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
rps4y1$sample_size = paste0(rps4y1$BD, "/", rps4y1$CTRL)

#perform meta-analysis
res9 = rma.uni(yi = logFC, sei = SE, data = rps4y1, method = "SJ")
#add weights
rps4y1$weights = round(weights(res9),1)


#LINC02192
linc02192 = metadata_no_covars %>%
  filter(grepl("^LINC02192", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
linc02192$sample_size = paste0(linc02192$BD, "/", linc02192$CTRL)
#perform meta-analysis
res10 = rma.uni(yi = logFC, sei = SE, data = linc02192, method = "SJ")
#add weights
linc02192$weights = round(weights(res10),1)

#SST
sst = metadata_no_covars %>%
  filter(grepl("^SST$", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
sst$sample_size = paste0(sst$BD, "/", sst$CTRL)
#perform meta-analysis
res11 = rma.uni(yi = logFC, sei = SE, data = sst, method = "SJ")
#add weights
sst$weights = round(weights(res11),1)

#LOC105372442
loc105372442 = metadata_no_covars %>%
  filter(grepl("^LOC105372442$", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
loc105372442$sample_size = paste0(loc105372442$BD, "/", loc105372442$CTRL)
#perform meta-analysis
res12 = rma.uni(yi = logFC, sei = SE, data = loc105372442, method = "SJ")
#add weights
loc105372442$weights = round(weights(res12),1)


#CCL3
ccl3 = metadata_no_covars %>%
  filter(grepl("^CCL3$", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
ccl3$sample_size = paste0(ccl3$BD, "/", ccl3$CTRL)
#perform meta-analysis
res13 = rma.uni(yi = logFC, sei = SE, data = ccl3, method = "SJ")
#add weights
ccl3$weights = round(weights(res13), 1)

#DDX3Y
ddx3y = metadata_no_covars %>%
  filter(grepl("^DDX3Y$", ignore.case = T, metadata_no_covars$`Gene Symbol`))
#add sample size
ddx3y$sample_size = paste0(ddx3y$BD, "/", ddx3y$CTRL)
#perform meta-analysis
res14 = rma.uni(yi = logFC, sei = SE, data = ddx3y, method = "SJ")
#add weights
ddx3y$weights = round(weights(res14), 1)


#combine forest plots for top-6 down-regulated genes
png("combined_forest_top6_DOWN_A_nocovars.png", width = 2200, height = 2400, res = 150)
par(mfrow = c(3, 2), mar = c(5, 10, 4, 2), oma = c(1, 1, 2, 1))  # for 6 plots in 3 rows, 2 columns

forest_rps4y1 = forest(res9, 
                       slab = rps4y1$Study, 
                       ilab = cbind(rps4y1$sample_size, rps4y1$weights),
                       ilab.xpos = c(-17, -14.2),
                       xlab = "log2 fold change", 
                       main = "RPS4Y1, Ribosome protein", 
                       cex = 1)
text(-17, nrow(rps4y1) + 2, "BD/CTRL", font = 2)
text(-14.2, nrow(rps4y1) + 2, "Weight, %", font = 2)

forest_ddx3y = forest(res14, 
                      slab = ddx3y$Study,
                      ilab = cbind(ddx3y$sample_size, ddx3y$weights),
                      ilab.xpos = c(-17.0, -14.0),
                      xlab = "log2 fold change", 
                      main = "DDX3Y, DEAD-box RNA helicase", 
                      cex = 1)
text(-17.0, nrow(ddx3y) + 2, "BD/CTRL", font = 2)
text(-14.0, nrow(ddx3y) + 2, "Weight, %", font = 2)

forest_linc02192 = forest(res10, 
                          slab = linc02192$Study,
                          ilab = cbind(linc02192$sample_size, linc02192$weights), 
                          ilab.xpos = c(-5.0, -4.2),
                          xlab = "log2 fold change", 
                          main = "LINC02192, Long intergenic non-coding RNA", 
                          cex = 1)
text(-5.0, nrow(linc02192) + 2, "BD/CTRL", font = 2)
text(-4.2, nrow(linc02192) + 2, "Weight, %", font = 2)

forest_sst = forest(res11, 
                    slab = sst$Study, 
                    ilab = cbind(sst$sample_size, sst$weights),
                    ilab.xpos = c(-7.2, -6.0),
                    xlab = "log2 fold change", 
                    main = "SST, Somatostatin", 
                    cex = 1)
text(-7.2, nrow(sst) + 2, "BD/CTRL", font = 2)
text(-6.0, nrow(sst) + 2, "Weight, %", font = 2)

forest_loc105372442 = forest(res12, 
                             slab = loc105372442$Study,
                             ilab = cbind(loc105372442$sample_size, loc105372442$weights),
                             ilab.xpos = c(-4.9, -4.0),
                             xlab = "log2 fold change", 
                             main = "LOC105372442, Uncharacterized/Predicted gene",
                             cex = 1)
text(-4.9, nrow(loc105372442) + 2, "BD/CTRL", font = 2)
text(-4.0, nrow(loc105372442) + 2, "Weight, %", font = 2)

forest_ccl3 = forest(res13, 
                     slab = ccl3$Study,
                     ilab = cbind(ccl3$sample_size, ccl3$weights),
                     ilab.xpos = c(-7.9, -6.6),
                     xlab = "log2 fold change", 
                     main = "CCL3, Small inducible cytokine", 
                     cex = 1)
text(-7.9, nrow(ccl3) + 2, "BD/CTRL", font = 2)
text(-6.6, nrow(ccl3) + 2, "Weight, %", font = 2)

dev.off()


###forest plots for the model with covariate adjustment
#add sample size to metadata
metadata_with_covars = inner_join(metadata_with_covars, sample_size, by=c("Study" = "study"))

#top up
loc102723630 = metadata_with_covars %>%
  filter(grepl("^LOC102723630$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
loc102723630$sample_size = paste0(loc102723630$BD, "/", loc102723630$CTRL)
#perform meta-analysis
res = rma.uni(yi = logFC, sei = SE, data = loc102723630, method = "SJ")
#add weights
loc102723630$weights = round(weights(res), 1)

efcab3p1 = metadata_with_covars %>%
  filter(grepl("^EFCAB3P1$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
efcab3p1$sample_size = paste0(efcab3p1$BD, "/", efcab3p1$CTRL)
#perform meta-analysis
res1 = rma.uni(yi = logFC, sei = SE, data = efcab3p1, method = "SJ")
#add weights
efcab3p1$weights = round(weights(res1), 1)

#HBG2
hbg2 = metadata_with_covars %>%
  filter(grepl("^HBG2$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
hbg2$sample_size = paste0(hbg2$BD, "/", hbg2$CTRL)
#perform meta-analysis
res2 = rma.uni(yi = logFC, sei = SE, data = hbg2, method = "SJ")
#add weights
hbg2$weights = round(weights(res2), 1)

#LOC102724046
loc102724046 = metadata_with_covars %>%
  filter(grepl("^LOC102724046$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
loc102724046$sample_size = paste0(loc102724046$BD, "/", loc102724046$CTRL)
#perform meta-analysis
res3 = rma.uni(yi = logFC, sei = SE, data = loc102724046, method = "SJ")
#add weights
loc102724046$weights = round(weights(res3), 1)

#BAALC-AS1
baalc_as1 = metadata_with_covars %>%
  filter(grepl("^BAALC-AS1$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
baalc_as1$sample_size = paste0(baalc_as1$BD, "/", baalc_as1$CTRL)
#perform meta-analysis
res4 = rma.uni(yi = logFC, sei = SE, data = baalc_as1, method = "SJ")
#add weights
baalc_as1$weights = round(weights(res4), 1)

#LOC100653133
loc100653133 = metadata_with_covars %>%
  filter(grepl("^LOC100653133$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
loc100653133$sample_size = paste0(loc100653133$BD, "/", loc100653133$CTRL)
#perform meta-analysis
res5 = rma.uni(yi = logFC, sei = SE, data = loc100653133, method = "SJ")
#add weights
loc100653133$weights = round(weights(res5), 1)

#combine plots for top6 up-regulated genes
png("combined_forest_top6_UP_A_withcovars.png", width = 2200, height = 2400, res = 150)
par(mfrow = c(3, 2), mar = c(5, 10, 4, 2), oma = c(1, 1, 2, 1))  # for 6 plots in 3 rows, 2 columns

forest_loc102723630 = forest(res, 
                             slab = loc102723630$Study,
                             ilab = cbind(loc102723630$sample_size, loc102723630$weights),
                             ilab.xpos = c(-7.5, -5.0),
                             xlab = "log2 fold change", 
                             main = "LOC102723630, Uncharacterized/Predicted gene", 
                             cex = 1)
text(-7.5, nrow(loc102723630) + 2, "BD/CTRL", font = 2)
text(-5.0, nrow(loc102723630) + 2, "Weight, %", font = 2)

forest_efcab3p1 = forest(res1, 
                         slab = efcab3p1$Study,
                         ilab = cbind(efcab3p1$sample_size, efcab3p1$weights),
                         ilab.xpos = c(-2.0, -1.0),
                         xlab = "log2 fold change", 
                         main = "EFCAB3P1, Pseudogene", 
                         cex = 1)
text(-2.0, nrow(efcab3p1) + 2, "BD/CTRL", font = 2)
text(-1.0, nrow(efcab3p1) + 2, "Weight, %", font = 2)

forest_hbg2 = forest(res2, 
                     slab = hbg2$Study,
                     ilab = cbind(hbg2$sample_size, hbg2$weights),
                     ilab.xpos = c(-6.3, -5.0),
                     xlab = "log2 fold change", 
                     main = "HBG2, Hemoglobin subunit", 
                     cex = 1)
text(-6.3, nrow(hbg2) + 2, "BD/CTRL", font = 2)
text(-5.0, nrow(hbg2) + 2, "Weight, %", font = 2)

forest_loc102724046 = forest(res3, 
                             slab = loc102724046$Study,
                             ilab = cbind(loc102724046$sample_size, loc102724046$weights),
                             ilab.xpos = c(-4.0, -3.0),
                             xlab = "log2 fold change", 
                             main = "LOC102724046, Uncharacterized/Predicted gene", 
                             cex = 1)
text(-4.0, nrow(loc102724046) + 2, "BD/CTRL", font = 2)
text(-3.0, nrow(loc102724046) + 2, "Weight, %", font = 2)

forest_baalc = forest(res4, 
                      slab = baalc_as1$Study,
                      ilab = cbind(baalc_as1$sample_size, baalc_as1$weights),
                      ilab.xpos = c(-3.0, -2.0),
                      xlab = "log2 fold change", 
                      main = "BAALC-AS1, BBALC anti-sense RNA", 
                      cex = 1)
text(-3.0, nrow(baalc_as1) + 2, "BD/CTRL", font = 2)
text(-2.0, nrow(baalc_as1) + 2, "Weight, %", font = 2)

forest_loc100653133 = forest(res5, 
                             slab = loc100653133$Study,
                             ilab = cbind(loc100653133$sample_size, loc100653133$weights),
                             ilab.xpos = c(-3.5, -2.5),
                             xlab = "log2 fold change", 
                             main = "LOC100653133, Uncharacterized/Predicted gene", 
                             cex = 1)
text(-3.5, nrow(loc100653133) + 2, "BD/CTRL", font = 2)
text(-2.5, nrow(loc100653133) + 2, "Weight, %", font = 2)

dev.off()


##top6 down-regulated genes with covars
#prepare the data and perform gene-level meta-analysis
linc02192 = metadata_with_covars %>%
  filter(grepl("^LINC02192$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
linc02192$sample_size = paste0(linc02192$BD, "/", linc02192$CTRL)
#perform meta-analysis
res6 = rma.uni(yi = logFC, sei = SE, data = linc02192, method = "SJ")
#add weights
linc02192$weights = round(weights(res6), 1)

linc01445 = metadata_with_covars %>%
  filter(grepl("^LINC01445$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
linc01445$sample_size = paste0(linc01445$BD, "/", linc01445$CTRL)
#perform meta-analysis
res7 = rma.uni(yi = logFC, sei = SE, data = linc01445, method = "SJ")
#add weights
linc01445$weights = round(weights(res7), 1)

sst = metadata_with_covars %>%
  filter(grepl("^SST$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
sst$sample_size = paste0(sst$BD, "/", sst$CTRL)
#perform meta-analysis
res8 = rma.uni(yi = logFC, sei = SE, data = sst, method = "SJ")
#add weights
sst$weights = round(weights(res8), 1)

loc107985303 = metadata_with_covars %>%
  filter(grepl("^LOC107985303$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
loc107985303$sample_size = paste0(loc107985303$BD, "/", loc107985303$CTRL)
#perform meta-analysis
res9 = rma.uni(yi = logFC, sei = SE, data = loc107985303, method = "SJ")
#add weights
loc107985303$weights = round(weights(res9), 1)

linc01978 = metadata_with_covars %>%
  filter(grepl("^LINC01978$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
linc01978$sample_size = paste0(linc01978$BD, "/", linc01978$CTRL)
#perform meta-analysis
res10 = rma.uni(yi = logFC, sei = SE, data = linc01978, method = "SJ")
#add weights
linc01978$weights = round(weights(res10), 1)

snord35a = metadata_with_covars %>%
  filter(grepl("^SNORD35A$", ignore.case = T, metadata_with_covars$`Gene Symbol`))
#add sample size
snord35a$sample_size = paste0(snord35a$BD, "/", snord35a$CTRL)
#perform meta-analysis
res11 = rma.uni(yi = logFC, sei = SE, data = snord35a, method = "SJ")
#add weights
snord35a$weights = round(weights(res11), 1)

#build combined forest plot
png("combined_forest_top6_DOWN_A_withcovars.png", width = 2200, height = 2400, res = 150)
par(mfrow = c(3, 2), mar = c(5, 10, 4, 2), oma = c(1, 1, 2, 1))  # for 6 plots in 3 rows, 2 columns

forest_sst = forest(res8, 
                             slab = sst$Study,
                             ilab = cbind(sst$sample_size, sst$weights),
                             ilab.xpos = c(-7.5, -6.0),
                             xlab = "log2 fold change", 
                             main = "SST, Somatostatin", 
                             cex = 1)
text(-7.5, nrow(sst) + 2, "BD/CTRL", font = 2)
text(-6.0, nrow(sst) + 2, "Weight, %", font = 2)

forest_linc02192 = forest(res6, 
                         slab = linc02192$Study,
                         ilab = cbind(linc02192$sample_size, linc02192$weights),
                         ilab.xpos = c(-5.0, -4.0),
                         xlab = "log2 fold change", 
                         main = "LINC02192, Long intergenic non-coding RNA", 
                         cex = 1)
text(-5.0, nrow(linc02192) + 2, "BD/CTRL", font = 2)
text(-4.0, nrow(linc02192) + 2, "Weight, %", font = 2)

forest_linc01445 = forest(res7, 
                     slab = linc01445$Study,
                     ilab = cbind(linc01445$sample_size, linc01445$weights),
                     ilab.xpos = c(-7.0, -5.5),
                     xlab = "log2 fold change", 
                     main = "LINC01445, Long intergenic non-coding RNA", 
                     cex = 1)
text(-7.0, nrow(linc01445) + 2, "BD/CTRL", font = 2)
text(-5.5, nrow(linc01445) + 2, "Weight, %", font = 2)

forest_loc107985303 = forest(res9, 
                             slab = loc107985303$Study,
                             ilab = cbind(loc107985303$sample_size, loc107985303$weights),
                             ilab.xpos = c(-3.0, -2.3),
                             xlab = "log2 fold change", 
                             main = "LOC107985303, Uncharacterized/Predicted gene", 
                             cex = 1)
text(-3.0, nrow(loc107985303) + 2, "BD/CTRL", font = 2)
text(-2.3, nrow(loc107985303) + 2, "Weight, %", font = 2)

forest_linc01978 = forest(res10, 
                      slab = linc01978$Study,
                      ilab = cbind(linc01978$sample_size, linc01978$weights),
                      ilab.xpos = c(-4.5, -3.3),
                      xlab = "log2 fold change", 
                      main = "LINC01978, Long intergenic non-coding RNA", 
                      cex = 1)
text(-4.5, nrow(linc01978) + 2, "BD/CTRL", font = 2)
text(-3.3, nrow(linc01978) + 2, "Weight, %", font = 2)

forest_snord35a = forest(res11, 
                             slab = snord35a$Study,
                             ilab = cbind(snord35a$sample_size, snord35a$weights),
                             ilab.xpos = c(-6.6, -5.3),
                             xlab = "log2 fold change", 
                             main = "SNORD35A, 	Small Nucleolar RNA", 
                             cex = 1)
text(-6.6, nrow(snord35a) + 2, "BD/CTRL", font = 2)
text(-5.3, nrow(snord35a) + 2, "Weight, %", font = 2)

dev.off()



