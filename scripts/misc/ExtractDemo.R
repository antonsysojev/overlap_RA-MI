### LAST VERSION UPDATE 18 AUG 2023 (v2.0.1) - MINOR QOL TO UPDATE TO NEW FOLDER STRUCTURE
### THIS SCRIPT EXTRACTS COHORT DEMOGRAPHICS FOR THE RA GWAS POPULATION

.libPaths("H:/Programs/RLibrary/")
library(haven)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
setwd("H:/Projects/OVERLAP_RA-CVD/")

fam <- read.table("DATA/OUTPUT/DATA/RA.fam", header = F) %>% as_tibble() %>% select(IID = V2)
eira_srqb_ID <-read.table("TEMPORARY/tmp-5/eira-plus-others-imputed.fam", header = F) %>% as_tibble() %>% select(IID = V2) %>% mutate(TYPE = "EIRA-SRQB")
twingene_ID <- read_tsv("TEMPORARY/tmp-5/TwinGene_HRC_imputation.psam", show_col_types = F) %>% mutate(IID = paste0(`#FID`, "_", IID)) %>% select(IID) %>% mutate(TYPE = "TWINGENE")
salty_ID <- read_tsv("TEMPORARY/tmp-5/SALTY_HRC_imputation.psam", show_col_types = F) %>% mutate(IID = paste0(`#FID`, "_", IID)) %>% select(IID) %>% mutate(TYPE = "SALTY")

### 1 - APPEND TYPE ONTO THE RAW IDS

long_ID <- eira_srqb_ID %>% bind_rows(salty_ID) %>% bind_rows(twingene_ID) %>% distinct()
fam_TYPE <- fam %>% inner_join(long_ID, by = "IID")

rm(eira_srqb_ID); rm(twingene_ID); rm(salty_ID); rm(long_ID); rm(fam)

### 2 - APPEND SEX AND AGE ONTO THE RAW IDS

age <- read_tsv("TEMPORARY/tmp-5/QCOV.txt") %>% select(IID, AGE)    #SEE NOTES
sex <- read.table("DATA/OUTPUT/DATA/RA.fam", header = F) %>% as_tibble() %>% select(IID = V2, SEX = V5)

fam_TYPE_SEX_AGE <- fam_TYPE %>% left_join(age, by = "IID") %>% left_join(sex, by = "IID")    #USE LEFT_JOIN TO ALLOW MISSINGNESS

rm(fam_TYPE); rm(age); rm(sex)

### 3 - APPEND PHENOTYPE STATUS ONTO THE RAW IDS

ra <- read.table("DATA/OUTPUT/DATA/RA.fam", header = F) %>% as_tibble() %>% select(IID = V2, RA = V6)

fam_TYPE_SEX_AGE_RA <- fam_TYPE_SEX_AGE %>% inner_join(ra, by = "IID") %>% mutate(RA = ifelse(RA %in% c(-9, 1), 1, RA))

rm(fam_TYPE_SEX_AGE); rm(ra)

### 4 - APPEND SEROSTATUS ONTO THE RAW IDS

spos_cohort <- read.table("TEMPORARY/tmp-5/SPOS_RA_TargetInds.txt", header = F) %>% as_tibble() %>% mutate(SPOS_COHORT = 1) %>% select(IID = V2, SPOS_COHORT)
sneg_cohort <- read.table("TEMPORARY/tmp-5/SNEG_RA_TargetInds.txt", header = F) %>% as_tibble() %>% mutate(SNEG_COHORT = 1) %>% select(IID = V2, SNEG_COHORT)    #SEE NOTES

fam_TYPE_SEX_AGE_RA_SPOS <- fam_TYPE_SEX_AGE_RA %>% left_join(spos_cohort, by = "IID") %>% left_join(sneg_cohort, by = "IID") %>%
  mutate(SEROSTATUS = case_when(RA == 2 & SPOS_COHORT == 1 ~ "SEROPOS RA", RA == 2 & SNEG_COHORT == 1 ~ "SERONEG RA", RA == 1 ~ "CONTROL", .default = "MISSING")) %>%
  select(-SPOS_COHORT, -SNEG_COHORT)

rm(spos_cohort); rm(sneg_cohort); rm(fam_TYPE_SEX_AGE_RA)

### 5 - APPEND CASE CHARACTERISTICS

srq <- read_sas("K:/Reuma/RASPA 2021/01. Data Warehouse/01. Processed Data/01. SRQ/srq_basdata.sas7bdat") %>% distinct(pid, debutalder, disease_debut_1, inklusions_ar, inklusions_manad)

srq_clean <- srq %>% mutate(inclusion_time = str_c(inklusions_ar, "-", inklusions_manad, "-01") %>% as_date()) %>% select(-inklusions_ar, -inklusions_manad) %>%
  mutate(duration_months = round(interval(disease_debut_1, inclusion_time) / months(1), 0)) %>% select(-disease_debut_1, -inclusion_time)

eira_key <- read_tsv("DATA/RAW/Genotyped_EIRA_from_Leonid_20190823_IDconv.txt", show_col_types = F) %>% select(EIRA_id = EIRA, GWAS = SMP.number)
srqb_key <- read_sas("DATA/RAW/key_20230126.sas7bdat") %>% mutate(SRQ_id = as.character(SRQ_id), BARCODE = as.character(barcode_deCODE_b1)) %>% select(SRQ_id, GWAS = BARCODE)
KEY <- read_tsv("DATA/RAW/KEY_pid-GWAS.txt", show_col_types = F)     #SEE NOTES ABOUT THIS
KEY_FULL <- KEY %>% left_join(eira_key %>% select(COHORT_ID = EIRA_id, IID = GWAS) %>% bind_rows(srqb_key %>% select(COHORT_ID = SRQ_id, IID = GWAS)), by = "COHORT_ID") %>% distinct(pid, IID) %>% mutate(IID = str_c("FAM001_", IID))

fam_TYPE_SEX_AGE_RA_SPOS_PID <- fam_TYPE_SEX_AGE_RA_SPOS %>% left_join(KEY_FULL, by = "IID") %>% distinct()
fam_TYPE_SEX_AGE_RA_SPOS_PID_SRQ <- fam_TYPE_SEX_AGE_RA_SPOS_PID %>% left_join(srq_clean, by = "pid") %>% select(IID, pid, everything())

rm(eira_key); rm(fam_TYPE_SEX_AGE_RA_SPOS); rm(fam_TYPE_SEX_AGE_RA_SPOS_PID); rm(KEY); rm(KEY_FULL); rm(srq); rm(srq_clean); rm(srqb_key)

### 6 - WRITE TABLE

write.table(fam_TYPE_SEX_AGE_RA_SPOS_PID_SRQ, "TEMPORARY/tmp-5/Demographics_RAW.txt", quote = F, row.names = F, col.names = T, sep = "\t")

### TO DO
### NOTES
# 2.1. The file `QCOV.txt` should be in /TEMPORARY/tmp-5/ given that the below code has been ran previously
#       from within the shell at the OVERLAP projects homefolder:
#       > Rscript SCRIPTS/MISC/ExtractQCov 1
#       > mv TEMPORARY/tmp-2/QCOV.txt TEMPORARY/tmp-5/
#       Then everything should work as expected.
# 2.2. Just like for the file `QCOV.txt`, the files `SPOS_RA_TargetInds.txt` and its seronegative analogue will
#       be in the appropriate folders given that the below lines of code are ran from within the shell at 
#       the OVERLAP projects homefolder prior:
#       > Rscript SCRIPTS/MISC/ExtractSubPheno SPOS_RA
#       > Rscript SCRIPTS/MISC/ExtractSubPheno SNEG_RA
#       > mv TEMPORARY/tmp-2/SPOS_RA_TargetInds.txt TEMPORARY/tmp-5/
#       > mv TEMPORARY/tmp-5/SNEG_RA_TargetInds.txt TEMPORARY/tmp-5/
#       Then everything should be fine!
# 2.3. This key is produced by the helper script `ExtractGWASKey.R`, a copy of the script `CleanEIRA-SRQB_HELPER.R` from the MTX_GWAS project.
#       The key is pretty much identical, although it does not contain the non-Batch 1 SRQb individuals. I've only copied the script so as to not rely on files outside of this project folder.