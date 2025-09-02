
#combine results to one xlxs with multiple lists
s1 = enrichment_BP
s2 = enrichment_cov
s3 = enrichment_BP_ba9
s4 = enrichment_BP_ba9_cov
s5 = enrichment_BP_pfc
s6 = enrichment_BP_pfc_cov
s7 = enrichment_blood_NO_covars
s8 = enrichment_blood_WITH_covars

library(openxlsx)

dfs = list(s1, s2, s3, s4, s5, s6, s7, s8)  

wb = createWorkbook()

for (i in seq_along(dfs)) {
  sheet_name <- i
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, dfs[[i]])
}

saveWorkbook(wb, "Supplementary_Table_3.xlsx", overwrite = TRUE)
