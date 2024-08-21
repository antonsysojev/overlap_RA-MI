### LAST VERSION UPDATE 29 MAY 2023 (v1.2) - REMOVED THE VALIDATION DATA
### THIS SCRIPT CREATES SUBFILES OF THE SRQ LINKAGE THAT CAN BE COPIED OVER TO THE LINUX SERVER FOR EASIER USE

library(dplyr)
library(haven)
setwd("H:/Projects/OVERLAP_RA-CVD/")

srq <- read_sas("K:/Reuma/RASPA 2021/01. Data Warehouse/01. Processed Data/01. SRQ/srq_basdata.sas7bdat") %>% select(pid, biobank_ID)
eira <- read_sas("K:/Reuma/RASPA 2021/01. Data Warehouse/01. Processed Data/02. EIRA/eira.sas7bdat") %>% select(pid, eira)

srq_eira_key <- srq %>% distinct() %>% mutate(COHORT_ID = as.character(biobank_ID), TYPE = "SRQ") %>% select(pid, COHORT_ID, TYPE) %>%
  bind_rows(eira %>% distinct() %>% mutate(TYPE = "EIRA") %>% select(pid, COHORT_ID = eira, TYPE)) %>% distinct() %>%
  filter(!is.na(COHORT_ID))

write.table(srq_eira_key, "DATA/RAW/KEY_pid-GWAS.txt", col.names = T, row.names = F, sep = "\t", quote = F)