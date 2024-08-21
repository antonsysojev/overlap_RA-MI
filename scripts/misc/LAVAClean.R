### LAST VERSION UPDATE 14 JULY 2023 (v1.0)
### THIS SCRIPT PROCESSES THE RAW LAVA OUTPUT

.libPaths("H:/Programs/RLibrary/")

library(dplyr)
library(stringr)
library(foreach)
library(tidyr)
library(readr)
library(purrr)

PHENOTYPE <- "SPOS_RA"

### 1. IDENTIFY WHICH FILES TO READ

TARGET_FILES <- foreach(CHROM = 1:22, .combine = list, .multicombine = T) %do%{
  
  dir("H:/Projects/OVERLAP_RA-CVD/TEMPORARY/CACHE/") %>% 
    str_subset(str_c("LAVA_", CHROM, "_")) %>%     #EXTRACTS THE FILENAMES FOR THE TARGET CHR
    str_extract("\\d+\\-\\d+") %>%     #EXTRACTS THE LOCI RANGES
    str_split_fixed("\\-", 2) %>% as.data.frame() %>% as_tibble() %>%
    mutate(START = V1 %>% as.numeric(), STOP = V2 %>% as.numeric()) %>% select(START, STOP) %>%
    arrange(START, STOP) %>%
    group_by(START) %>% filter(STOP == max(STOP)) %>% ungroup() %>%
    mutate(CHR = CHROM, .before = everything()) %>%
    mutate(READSTRING = str_c("LAVA_", CHR, "_", START, "-", STOP, ".rds"))
  
}
TARGET_FILES <- TARGET_FILES %>% bind_rows()

### 2. READ AND PROCESS EACH INPUT
### ### 2.1. READ AND PROCESS PRIMARY INPUT

source("H:/Projects/OVERLAP_RA-CVD/SCRIPTS/MISC/PROCESS_LAVA_fun.R")

FIRST_CLEAN <- foreach(i = 1:nrow(TARGET_FILES), .combine = list, .multicombine = T) %do% {
  
  if(i %% 10 == 0) print(str_c("PROCESSING SEGMENT ", i, "..."))
  
  LAVA_RES <- readRDS(str_c("H:/Projects/OVERLAP_RA-CVD/TEMPORARY/CACHE/", TARGET_FILES[i, ]$READSTRING))
  TARGET_FILES_SEQUENCE <- seq(TARGET_FILES[i, ]$START, TARGET_FILES[i, ]$STOP)
  LAVA_RES_SUB <- LAVA_RES[TARGET_FILES_SEQUENCE]
  
  map2(LAVA_RES_SUB, TARGET_FILES_SEQUENCE, PROCESS_LAVA, PHENOTYPE = PHENOTYPE, CHROM = TARGET_FILES[i, ]$CHR)
  
}
FIRST_CLEAN_df <- FIRST_CLEAN %>% bind_rows() %>% as_tibble()

rm(FIRST_CLEAN); rm(LAVA_RES); rm(LAVA_RES_SUB); rm(TARGET_FILES); rm(CHROM); rm(i)

### ### 2.2. READ AND PROCESS SECONDARY INPUT

LAVA_RES <- readRDS("H:/Projects/OVERLAP_RA-CVD/TEMPORARY/LAVA_MISC_1.rds")
SECOND_CLEAN_1 <- lapply(LAVA_RES, PROCESS_LAVA, PHENOTYPE = PHENOTYPE) %>% bind_rows() %>% filter(CHR == 1)    #SEE NOTE

LAVA_RES <- readRDS("H:/Projects/OVERLAP_RA-CVD/TEMPORARY/LAVA_MISC_2.rds")
SECOND_CLEAN_2 <- lapply(LAVA_RES[[1]], PROCESS_LAVA, PHENOTYPE = PHENOTYPE) %>% bind_rows()
SECOND_CLEAN_3 <- lapply(LAVA_RES[2:length(LAVA_RES)], PROCESS_LAVA, PHENOTYPE = PHENOTYPE) %>% bind_rows()

SECOND_CLEAN_df <- SECOND_CLEAN_1 %>% bind_rows(SECOND_CLEAN_2) %>% bind_rows(SECOND_CLEAN_3) %>% select(CHR, LOCUS = LOC, everything()) %>% as_tibble()

### 3. COMBINE AND PROCESS TOTAL INPUT

LOCI <- read_delim("H:/Projects/OVERLAP_RA-CVD/DATA/RAW/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile", delim = " ")
LOCI_ALT <- LOCI %>% group_by(CHR) %>% mutate(IDX = row_number(CHR)) %>% ungroup() %>% mutate(ID = str_c(CHR, "_", IDX)) %>% select(ID, START, STOP)

FIRST_CLEAN_df %>% bind_rows(SECOND_CLEAN_df) %>% arrange(CHR, LOCUS) %>% distinct() %>%
  select(CHR, LOCUS, h2_RA = h2.obs_RA_DECODE, p_RA = p_RA_DECODE, h2_MI = h2.obs_MI, p_MI, rho, rho.lower, rho.upper, p, r2, r2.lower, r2.upper) %>% 
  group_by(CHR) %>% mutate(IDX = row_number(CHR), .after = IDX) %>% ungroup() %>% select(-LOCUS) %>%     #SEE NOTE
  mutate(ID = str_c(CHR, "_", IDX)) %>% inner_join(LOCI_ALT, by = "ID") %>%
  select(CHR, IDX, START, STOP, h2_RA, p_RA, h2_MI, p_MI, rho, rho.lower, rho.upper, p, r2, r2.lower, r2.upper) %>% 
  write.table(str_c("H:/Projects/OVERLAP_RA-CVD/DATA/OUTPUT/RES/", PHENOTYPE, "_LAVA_CLEAN.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

### TO DO:
# 1.1. You need to deal with the missing IDs for loci somehow. 
#       We can do it so that the function ONLY uses input if there isn't already input (i.e. length 2).
#       Then we can simply always give it an input.
# 1.2. Final chunk of code is bugged and tries to extract columns with names that are hardcoded and do not
#       take the target phenotype into account. Can be fixed by doing positional extractions and then
#       renaming accordingly but I'll leave this for now.
### NOTES:
# 2.1. There's an issue with duplicates. I deal with it by cutting only to the loci on CHR1, which was the target
#       loci on the processing run.