#!/usr/bin/env Rscript
### LAST UPDATED 26 MAY 2023 (v1.2) - ADDED FUNCTIONALITY TO SKIP THE PCA PARTS TO ONLY EXTRACT AGE (FOR DEMOGRAPHICS)
### THIS SCRIPT EXTRACTS THE AGE OF ALL INDIVIDUALS WITHIN THE STR DATA

.libPaths("/home2/genetics/antobe/software/RLibrary/")
suppressWarnings(suppressMessages(library(haven)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(lubridate)))
suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(argparser)))

### 1 - SETUP

args <- arg_parser('EMPTY FN') %>% add_argument('noPCA', help = 'GIVE "1" IF YOU WISH TO SKIP THE PCA STUFF AND JUST GET THE AGES', type = 'character') %>% parse_args()
noPCA <- as.logical(as.numeric(args$noPCA))
#print(noPCA)

if(!noPCA){pca_df <- read.table("TEMPORARY/tmp-2/PCA.txt", header = F); colnames(pca_df) <- c("FID", "IID", str_c("V", seq(3, ncol(pca_df)))); pca_df <- pca_df %>% as_tibble()}
#print("HI!")
fam_df <- read_delim("DATA/OUTPUT/DATA/RA.fam", delim = " ", col_names = c("FID", "IID", "V3", "V4", "SEX", "PHENO"), show_col_types = F)
str_var <- read_sas("/home2/genetics/antobe/data/STR/variables/ja_all_individ.sas7bdat")

### 2 - EXTRACT AGE FOR EIRA-SRQB

eira_srqb_AGE <- read_tsv("DATA/RAW/EIRA_SRQB_COV.txt", show_col_types = F) %>% distinct(GWAS, ALDER, debutalder)

eira_srqb_AGE_CLEAN <- eira_srqb_AGE %>% filter(!is.na(GWAS)) %>%     #REMOVE ALL WITH MISSING GWAS ID SINCE THESE CAN NOT BE LINKED TO GENOTYPED DATA ANYWAYS
  mutate(GWAS = str_c("FAM001_", GWAS)) %>% distinct(GWAS, ALDER, debutalder) %>%     #REMOVE PID AND REMOVE DUPLICATES
  mutate(AGE = case_when(is.na(ALDER) & is.na(debutalder) ~ NA,
                         !is.na(ALDER) & is.na(debutalder) ~ ALDER,
                         is.na(ALDER) & !is.na(debutalder) ~ debutalder,
                         !is.na(ALDER) & !is.na(debutalder) ~ debutalder)) %>% distinct(GWAS, AGE) %>%     #AGGREGATE THE VARIABLES INTO ONE AGE VARIABLE AND REMOVE DUPLICATES
  filter(!is.na(AGE)) %>%     #REMOVE THOSE WITH MISSING AGE - THESE WILL NOT BE AVAILABLE REGARDLESS
  group_by(GWAS) %>% filter(!(n() > 1 & AGE == min(AGE))) %>% ungroup()    #REMOVE A COPY OF THOSE THAT HAVE TWO REPORTED AGES - SEE NOTES

fam_df_1 <- fam_df %>% left_join(eira_srqb_AGE_CLEAN, by = c("IID" = "GWAS")) %>% select(everything(), ES_AGE = AGE)

### 3 - EXTRACT AGE FOR TWINGENE

tg_df <- read_sas("/home2/genetics/antobe/data/STR/variables/ja_consent_date.sas7bdat")
tg_df_AGE_CLEAN <- tg_df %>% inner_join(str_var %>% select(bbTWIN, fodelse), by = "bbTWIN") %>%
	mutate(CONSENT_DATE = as_date(CONSENT_DATE), fodelse = as_date(fodelse)) %>%
	mutate(AGE = floor(interval(fodelse, CONSENT_DATE) / years(1))) %>%
	filter(!is.na(AGE)) %>%
	mutate(IID = paste0(bbPAIR, "_", bbTWIN)) %>%
	distinct(IID, AGE) %>% suppressWarnings()    #SEE NOTE

fam_df_2 <- fam_df_1 %>% left_join(tg_df_AGE_CLEAN, by = "IID") %>% select(everything(), TG_AGE = AGE)

### 4 - EXTRACT AGE FOR SALTY

salty_df <- read_sas("/home2/genetics/antobe/data/STR/variables/ja_salty_received.sas7bdat")
salty_df_AGE_CLEAN <- salty_df %>% mutate(bbPAIR = str_sub(bbTWIN, 1, -2)) %>%
	inner_join(str_var %>% select(bbTWIN, fodelse), by = "bbTWIN") %>% 
	mutate(received = as_date(received), fodelse = as_date(fodelse)) %>% 
	mutate(AGE = floor(interval(fodelse, received) / years(1))) %>% 
	mutate(IID = paste0(bbPAIR, "_", bbTWIN)) %>%
	distinct(IID, AGE)

fam_df_3 <- fam_df_2 %>% left_join(salty_df_AGE_CLEAN, by = "IID") %>% select(everything(), SALTY_AGE = AGE)

### 5 - COMBINE TOGETHER AND OUTPUT

QCOV <- fam_df_3 %>% mutate(AGE = case_when(!is.na(ES_AGE) ~ ES_AGE, !is.na(TG_AGE) ~ TG_AGE, !is.na(SALTY_AGE) ~ SALTY_AGE)) %>% select(FID, IID, AGE) %>% filter(AGE > 0)    #SEE NOTE ABOUT THIS
if(!noPCA){
  QCOV <- QCOV %>% inner_join(pca_df %>% select(-FID), by = "IID")
  colnames(QCOV) <- c("FID", "IID", "AGE", paste0("PC", seq(1, ncol(pca_df) - 2)))
}else{
  colnames(QCOV) <- c("FID", "IID", "AGE")
}

write.table(QCOV, "TEMPORARY/tmp-2/QCOV.txt", col.names = T, row.names = F, quote = F, sep = "\t")

### TO DO:
# 1.1. May require some future cleaning of data prior to this? We have 410 individuals from STR without an AGE as well, which isn't great...
### NOTES:
# 2.1. Note the use of `suppressWarnings()`. Five individuals within the `tg_df` data has impossible dates for the variable `CONSENT_DATE` 
#       (you can find them by going `filter(is.na(AGE))` at the same step where we currently remove the missing AGEs.
#	      Using `as_date()` on these leads to warnings being raised, but I've checked these and the current handling of them (making them NA) is FINE since we can't figure out what they were supposed to be.
# 2.2. Note the use of `case_when()`. This function is hierarchical, meaning that it does the first argument first, then the second one and then the third one.
#	      This means that it NATURALLY does what we want it to do with respect to those with multiple values.
#	      Here, there will be no overlapping values between STR-cohorts and the EIRA-SRQB individuals. However, there may be between TwinGene and SALTY.
#	      In the case of the latter, we wish to prioritize the TG_AGE, since it is the age of the participant at the time of consent, while SALTY has the age when the sample received.
# 2.3. Note that we filter out individuals with AGE < 0. We do this because one EIRA-SRQB participant has an AGE of ~-2000. This should not go into our data.
# 2.4. There are people in EIRA who are reported twice, which are the final ones to be cut here. Here, I simply chose to cut the first record of them, keeping their "oldest" age. This is 
#       simply thinking that the second one may be a change from control to case or something like that.
