#!/user/bin/env Rscript
### LAST VERSION UPDATE - 21 FEB 2023 (v2.1.1) - FIXED A MINOR BUG WITH RESPECT TO FORMAT BY ADDING AN `as.numeric()`; PRETTIER SOLUTIONS WARRANTED
### TAKES THE PLINK2 PCA OUTPUT AS INPUT FILTERS OUT INDIVIDUALS ANALOGOUSLY TO EIGENSOFT OUTPUTS A LIST OF IDS FOR EXCLUSION VIA PLINK

#! NOTE THAT MORE RECENT VERSIONS ARE LIKELY AVAILABLE, SEE FOR INSTANCE THE GENOTYPE IMPUTATION PROJECT

.libPaths('/home2/genetics/antobe/software/RLibrary/')
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(readr)))

pca_df <- read_delim("pca.eigenvec", show_col_types = F, col_names = c("FID", "IID", paste0("PC", 1:10)), skip = 1)    #SKIP LETS YOU CUT THE HEADER AND USE YOUR OWN
sigma <- 6    #N STANDARD DEVIATIONS USED FOR FILTRATIONS - SIX IS DEFAULT IN EIGENSOFT

mean_pci <- pca_df %>% summarise(across(3:ncol(pca_df), mean))
sd_pci <- pca_df %>% summarise(across(3:ncol(pca_df), sd))
target_Lbound <- mean_pci - sigma * sd_pci
target_Ubound <- mean_pci + sigma * sd_pci

outlier_IDs <- list()
for(i in 1:(ncol(pca_df)-2)){
  outlier_IDs[[i]] <- pca_df %>% select(1:2, PC_TARGET = i + 2) %>% mutate(FAILURE = ifelse(between(PC_TARGET, target_Lbound[i] %>% as.numeric(), target_Ubound[i] %>% as.numeric()), 0, 1)) %>% filter(FAILURE == 1) %>% select(1:2)
}

bind_rows(outlier_IDs) %>% distinct() %>% write.table("pca_outliers.txt", quote = F, row.names = F, col.names = F, sep = "\t")
