#!/user/bin/env Rscript
### LAST UPDATED 8 MAY 2023 (v1.1.1) - FIXED A MINOR BUG
### THIS SCRIPT FILTERS THE GENOTYPED RA DATA INTO SUBSETS FOR ANALYSIS OF SEROPOSITIVE AND SERONEGATIVE RA

.libPaths("/home2/genetics/antobe/software/RLibrary")
suppressWarnings(suppressMessages(library(haven)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(readr)))

args <- arg_parser('EMPTY FOR NOW') %>% add_argument('PHENO', help = 'TARGET PHENOTYPE, VALID ARE `RA`, `SPOS_RA` AND `SNEG_RA`', type = 'character') %>% parse_args()
PHENO <- args$PHENO

fam <- read.table("DATA/OUTPUT/DATA/RA.fam") %>% as_tibble()
SEROSTAT <- read_tsv("DATA/RAW/EIRA_SRQB_COV.txt", show_col_types = F) %>% distinct(pid, GWAS, RF_EIRA, ACPA_EIRA, RF_SRQB, ACPA_SRQB)

if(PHENO == "RA"){

	fam_sub <- fam %>% distinct(V1, V2)	#IF TARGET PHENOTYPE IS RA THEN DO NOTHING

}else if(PHENO == "SPOS_RA"){
  
  SEROPOS <- SEROSTAT %>% filter(!is.na(GWAS)) %>% distinct(GWAS, RF_EIRA, ACPA_EIRA, RF_SRQB, ACPA_SRQB) %>%
    mutate(SEROPOS_EIRA = ifelse(RF_EIRA == 1 | ACPA_EIRA == 1, 1, NA), SEROPOS_SRQB = ifelse(RF_SRQB == 1 | ACPA_SRQB == 1, 1, NA)) %>% distinct(GWAS, SEROPOS_EIRA, SEROPOS_SRQB) %>%
    mutate(SEROPOS = ifelse(SEROPOS_EIRA == 1 | SEROPOS_SRQB == 1, 1, NA)) %>% distinct(GWAS, SEROPOS) %>%
    group_by(GWAS) %>% filter(!(n() > 1 & is.na(SEROPOS))) %>% ungroup() %>%
    mutate(GWAS = paste0("FAM001_", GWAS))
  
	fam_sub <- fam %>% left_join(SEROPOS, by = c("V2" = "GWAS")) %>% filter((V6 == 2 & is.na(SEROPOS))) %>% distinct(V1, V2)

}else if(PHENO == "SNEG_RA"){
  
  SERONEG <- SEROSTAT %>% filter(!is.na(GWAS)) %>% distinct(GWAS, RF_EIRA, ACPA_EIRA, RF_SRQB, ACPA_SRQB) %>%
    mutate(SERONEG_EIRA = ifelse(RF_EIRA == 0 & ACPA_EIRA == 0, 1, NA), SERONEG_SRQB = ifelse(RF_SRQB == 0 & ACPA_SRQB == 0, 1, NA)) %>% distinct(GWAS, SERONEG_EIRA, SERONEG_SRQB) %>%
    mutate(SERONEG = ifelse(SERONEG_EIRA == 1 | SERONEG_SRQB == 1, 1, NA)) %>% distinct(GWAS, SERONEG) %>%
    group_by(GWAS) %>% filter(!(n() > 1 & is.na(SERONEG))) %>% ungroup() %>%
    mutate(GWAS = paste0("FAM001_", GWAS))

	fam_sub <- fam %>% left_join(SERONEG, by = c("V2" = "GWAS")) %>% filter((V6 == 2 & is.na(SERONEG))) %>% distinct(V1, V2)

}

write.table(fam_sub, paste0("TEMPORARY/tmp-2/", PHENO, "_TargetInds.txt"), col.names = F, row.names = F, sep = " ", quote = F)

### TO DO:
# 1.1.
### NOTES:
# 2.1. 