#!/user/bin/env Rscript
### LAST VERSION UPDATE 30 MAY 2023 (v1.0.1) - EDITED THE OUTPUT OF THE RA GWAS FILE TO CONTAIN THE EXACT SAME COLUMNS AS THE OTHERS
### THIS SCRIPT PERFORMS SOME BASIC CLEANING OF GWAS DATA PRIOR TO MUNGING VIA LDSC

.libPaths("/home2/genetics/antobe/software/RLibrary/")
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(argparser)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(stringr)))

args <- arg_parser('EMPTY FN') %>% add_argument('GWASPATH', help = 'FILEPATH TO THE TARGET GWAS', type = 'character') %>% parse_args()
FILEPATH <- args$GWASPATH
PHENO <- FILEPATH %>% str_remove_all("_GWAS\\..+$") %>% str_extract("\\/[A-Z_]+$") %>% str_remove_all("\\/")

### 3.4.1. READ THE DATA

if(PHENO %in% c("RA", "SPOS_RA", "SNEG_RA")){gwas <- read_tsv(FILEPATH, show_col_types = F, col_types = "ccdccdddddd")}    #READ TAB-SEP WITH FIXING OF COLUMN TYPES
if(PHENO %in% c("MI", "CRP", "LDL", "SMKINIT", "SMKAMT")){gwas <- read_tsv(FILEPATH, show_col_types = F)}    #READ TAB-SEP NATURALLY
if(PHENO %in% c("SBP", "DBP", "WHR", "BMI")){gwas <- read_delim(FILEPATH, delim = " ", show_col_types = F)}    #READ SPACE-SEP
if(PHENO == "PP"){gwas <- read_delim(FILEPATH, delim = " ", skip = 1, col_names = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P", "TotalSampleSize", "N_effective"), show_col_types = F)}
if(str_detect(PHENO, "DECODE")){gwas <- read_table(FILEPATH, show_col_types = F, skip = 1, col_names = c("CHR", "POS", "SNP", "OA", "EA", "EAF", "Cohorts", "OR", "P", "Phet", "I2"), col_types = "ccccccccccc")}    #SEE NOTE

### 3.4.2. CLEAN THE DATA

if(PHENO == "MI"){N_IND <- 638176}
if(PHENO == "CRP"){N_IND <- 575531}
if(PHENO == "LDL"){N_IND <- 431167}

if(PHENO %in% c("RA", "SPOS_RA", "SNEG_RA")){gwas_clean <- gwas %>% select(SNP, A1, A2, P, BETA, N)}
if(PHENO %in% c("MI", "CRP", "LDL")){gwas_clean <- gwas %>% mutate(SNP = str_c(chromosome, ":", base_pair_location), N = N_IND) %>% select(SNP, A1 = effect_allele, A2 = other_allele, P = p_value, BETA = beta, N)}
if(PHENO %in% c("SBP", "DBP", "PP")){gwas_clean <- gwas %>% mutate(SNP = str_remove(MarkerName, "\\:SNP$"), A1 = str_to_upper(Allele1), A2 = str_to_upper(Allele2)) %>% select(SNP, A1, A2, BETA = Effect, P, N = TotalSampleSize)}
if(PHENO %in% c("WHR", "BMI")){gwas_clean <- gwas %>% mutate(SNP = str_c(CHR, ":", POS)) %>% select(SNP, A1 = Tested_Allele, A2 = Other_Allele, N, P, BETA)}
if(PHENO %in% c("SMKINIT", "SMKAMT")){gwas_clean <- gwas %>% mutate(SNP = str_c(CHROM, ":", POS)) %>% select(SNP, A1 = ALT, A2 = REF, P = PVALUE, BETA, N)}

if(str_detect(PHENO, "DECODE")){
	gwas_filt <- gwas %>% filter(CHR != "chrX") %>% mutate(start = as.numeric(POS) - 1L, end = as.numeric(POS), id = seq_along(POS))	#NEED TO LIFT THE DECODE DATA FROM 38 TO 37
	info_BED <- gwas_filt %>% select(chrom = CHR, start, end, id)
	BED <- tempfile(fileext = ".BED")
	data.table::fwrite(info_BED, BED, col.names = F, sep = " ", scipen = 50)
	lifted <- tempfile(fileext = ".BED")
	system2("/home2/genetics/antobe/software/liftOver/liftOver", c(BED, "/home2/genetics/antobe/software/liftOver/hg38ToHg19.over.chain.gz", lifted, tempfile(fileext = ".txt")))
	
	info_BED_LIFTED <- read_table(lifted, col_names = c("CHR", "START", "END", "ID"), show_col_types = F)
	gwas_lift <- gwas_filt %>% inner_join(info_BED_LIFTED %>% select(START_L = START, END_L = END, ID), by = c("id" = "ID")) %>% select(CHR, POS = END_L, OA, EA, Cohorts, OR, P) %>% distinct()
	
	if(str_detect(PHENO, "SPOS")){N_vec <- c(315450, 29104, 15878, 91814, 408565, 148812)
	}else if(str_detect(PHENO, "SNEG")){N_vec <- c(323877, 28972, 11288, 89618, 408565, 144748)
	}else{N_vec <- c(345401, 29238, 18076, 94626, 408565, 130624)}
	COHORTS <- (gwas_lift %>% mutate(Cohorts = str_replace_all(Cohorts, "\\?", "0")) %>% mutate(Cohorts = str_replace_all(Cohorts, "\\+|\\-", "1")))$Cohorts     #SEE NOTES
	gwas_lift$N <- str_split_fixed(COHORTS, "", 6) %>% apply(FUN = as.numeric, MAR = c(1, 2)) %*% N_vec
	gwas_clean <- gwas_lift %>% mutate(SNP = str_c(CHR, ":", POS) %>% str_remove("chr"), BETA = log(as.numeric(OR))) %>% select(SNP, A1 = EA, A2 = OA, BETA, P, N)
}

write.table(gwas_clean, paste0("TEMPORARY/tmp-3/", PHENO, "_pMunge.txt"), row.names = F, col.names = T, quote = F, sep = " ")

### TO DO
# 1.1.
### NOTES
# 2.1. Just for your information, I've thought carefully about the allele pairings and made sure that they are correctly interpreted by LDSC. There should be no issue here but please, do double-check it against
#	README-files for the original GWAS data downloads if you are worried.
# 2.2. For smoking, things are a bit misleading. Data comes with reports of SNPs being `REF` and `ALT`. Checking these against dbSNP (10 random SNPs) gives me an overall good agreement and indicates that the reported allele
#	frequencies refer to the presence of the `ALT` allele. README tells me that the BETA effect is also with respect to the `ALT` allele, which means it makes sense to swap their positions as I want `A1` to refer
#	to the EFFECT allele with respect to the BETA.
# 2.3. Reading the deCODE GWAS data is a bit awkward and returns an error when we do it the normal way. It runs into some issues with the final column (I believe?) but it is mostly related to how the HEADER connects to the
#	actual lines of data. I fix this by skipping the first line (the header) and simply naming the columns myself. It can probably be read regardless, but this line of code means I avoid any warnings. Also, note that
#	I set all observations to 'character'. It doesn't really matter, but this seems like a safe way to read everything into R and deal with character types later.
# 2.4. The DECODE data does not contain explicit information on SNP sample size. However, it does contain implicit information via the `Cohorts` variable, which details whether a SNP provided a positive/negative effect
#	to the meta-analysis or whether the SNP was even available across the individual cohorts. Using this information as well as secondary data on which character corresponds to which cohort (both available within
#	the data README), allows us to infer the per-SNP sample size to a closer degree than if we had just used the full N. I then find the sample size contributed by each of the six individual cohorts and use this
#	information here as the cap.
