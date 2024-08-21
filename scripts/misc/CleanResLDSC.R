#!/user/bin/env Rscript
### LAST VERSION UPDATE 5 APRIL 2023 (v1.0) -
### THIS SCRIPT LOADS ALL THE LDSC RESULT FILES AND PROCESSES THEM INTO A NEAT TABLE

.libPaths("../../software/RLibrary/")

library(dplyr)
library(readr)
library(stringr)

res_files <- list.files("DATA/OUTPUT/RES")
pheno_pairs <- res_files %>% str_remove("_ldsc\\.log$")

res_list <- list()
for(i in 1:length(pheno_pairs)){
  
  res <- read_lines("H:/Projects/OVERLAP_RA-CVD/TMP/FAKE_LDSC.log", skip_empty_rows = T)
  
  h2_res <- res[res %>% str_which("Total Observed scale h2")]
  rg_res <- res[res %>% str_which("Genetic Correlation\\:")]
  rgp_res <- res[res %>% str_which("P\\:")]
  
  res_df <- data.frame(h21 = h2_res[1], h22 = h2_res[2], rg = rg_res, rgp = rgp_res)
  res_df_clean <- res_df %>% mutate(h21 = str_extract(h21, "\\:.+$"), h22 = str_extract(h22, "\\:.+$"), rg = str_extract(rg, "\\:.+$"), rgp = str_extract(rgp, "\\:.+$")) %>%
    mutate(h21 = str_remove(h21, "^\\: "), h22 = str_remove(h22, "^\\: "), rg = str_remove(rg, "^\\: "), rgp = str_remove(rgp, "^\\: "))
  res_list[[i]] <- res_df_clean
}
res_table <- bind_rows(res_list)
rownames(res_table) <- pheno_pairs

write.table(res_table, "DATA/OUTPUT/RES/RES_LDSC_CLEAN.txt", row.names = T, col.names = F, quote = F, sep = "\t")

### TO DO:
# 1.1.
### NOTES:
# 2.1.