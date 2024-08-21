### LAST VERSION UPDATE 5 MAY 2023 (v1.0)
### THIS SCRIPT EXTRACTS DATA FOR THE EIRA-SRQB GWAS PARTICIPANTS

setwd("H:/Projects/OVERLAP_RA-CVD/")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(haven)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(readxl)))

eira_ids <- read_sas("K:/Reuma/RASPA 2021/01. Data Warehouse/01. Processed Data/02. EIRA/eira.sas7bdat") %>% select(pid, COHORT_ID = eira) %>% distinct()
srqb_ids <- read_sas("K:/Reuma/RASPA 2021/01. Data Warehouse/01. Processed Data/01. SRQ/srq_basdata.sas7bdat") %>% mutate(COHORT_ID = as.character(biobank_ID)) %>% select(pid, COHORT_ID) %>% distinct()
eira_srqb_ids <- eira_ids %>% bind_rows(srqb_ids) %>% distinct()

eira_key <- read_tsv("DATA/RAW/Genotyped_EIRA_from_Leonid_20190823_IDconv.txt", show_col_types = F) %>% select(COHORT_ID = EIRA, GWAS = SMP.number) %>% distinct()
srqb_key <- read_sas("DATA/RAW/key_20230126.sas7bdat") %>% mutate(COHORT_ID = as.character(SRQ_id), GWAS = as.character(barcode_deCODE_b1)) %>% distinct(COHORT_ID, GWAS)
eira_srqb_key <- eira_key %>% bind_rows(srqb_key) %>% distinct()

EIRA_SRQB <- eira_srqb_ids %>% left_join(eira_srqb_key, by= "COHORT_ID", na_matches = "never") %>% arrange(pid, COHORT_ID, GWAS) %>% distinct()

### 1 - APPEND AGE

EIRA_SRQB_AGE <- EIRA_SRQB %>% 
  left_join(read_sas("K:/Reuma/RASPA 2021/01. Data Warehouse/01. Processed Data/02. EIRA/eira.sas7bdat") %>% distinct(pid, ALDER), by = "pid", relationship = "many-to-many") %>%     #SEE NOTES
  left_join(read_sas("K:/Reuma/RASPA 2021/01. Data Warehouse/01. Processed Data/01. SRQ/srq_basdata.sas7bdat") %>% distinct(pid, debutalder), by = "pid", relationship = "many-to-many")    #SEE NOTES

### 2 - APPEND SEROSTATUS

eira_rf <- read_xlsx("DATA/RAW/221215 serologi by Johan Rönnelid.xlsx") %>% select(COHORT_ID = `Sample clean ID`, RF = `IgG RF, microg/mL`) %>% mutate(RF = ifelse(RF < 25, 0, 1))
eira_acpa <- read_xlsx("DATA/RAW/221206 EIRA CCP till Helga.xlsx", col_names = paste0("V", 1:13), skip = 1) %>% mutate(ACPA = ifelse(V2 < 25, 0, 1)) %>% select(COHORT_ID = V1, ACPA)
srqb_rf_acpa <- read_sas("DATA/RAW/anton_ccp_rf.sas7bdat") %>% mutate(biobank_id = as.character(biobank_id)) %>% select(COHORT_ID = biobank_id, RF = RF_IgM, ACPA = antiCCP)

EIRA_SRQB_AGE_SERO <- EIRA_SRQB_AGE %>% 
  left_join(eira_rf %>% distinct(), by = "COHORT_ID") %>% 
  left_join(eira_acpa %>% distinct(), by = "COHORT_ID") %>% 
  left_join(srqb_rf_acpa %>% distinct(), by = "COHORT_ID", relationship = "many-to-many", suffix = c("_EIRA", "_SRQB"))     #SEE NOTES

### 3 - WRITE TABLE

write.table(EIRA_SRQB_AGE_SERO, "DATA/RAW/EIRA_SRQB_COV.txt", col.names = T, row.names = F, sep = "\t", quote = F)

### TO DO:
### NOTES:
# 2.1. The EIRA data contains multiple duplicates in `pid`, hence the need for `relationship = "many-to-many"`. Some of these are simply duplicated lines of data, 
#       some seem to be people once sampled as a case, once as a control. However, some of these are seemingly repeated sampling of individuals at later stages, which are somewhat unclear.
#       I chose to keep all of these, leaving the filtering of duplicates to a later stage.
# 2.2. Something similar is occurring with SRQb, though there's actually only one person who is duplicated. I chose to deal with this the same way: adding them onto my data and doing the selection later.
# 2.3. And again, we have multiples within the `srqb_rf_acpa` file that end up giving us warnings. I've silenced the warning here to deal with it at a later stage.
# 2.4. Note that this script can only be executed from H and not from the Linux-server, as it tries to reach folders on K that I can't copy over to the Linux-server. I recommend running this from
#      H, and moving the output files over to the Linux-server and using them that way.