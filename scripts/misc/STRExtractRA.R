#!/usr/bin/env Rscript
### LAST UPDATED 24 MARCH 2023 (v1.2) - UPDATED SCRIPT TO EXCLUDE ALSO THOSE WITH AN ICD-10 CODE FOR RA PREDATING INCLUSION
### THIS SCRIPT EXTRACTS ALL INDIVIDUALS WHO HAVE AN ICD-10 CODE FOR RA FROM THE STR DATA BASED ON LINKED INFORMATION FROM THE NPR

.libPaths("/home2/genetics/antobe/software/RLibrary/")
suppressWarnings(suppressMessages(library(haven, quietly = T)))
suppressWarnings(suppressMessages(library(dplyr, quietly = T)))
suppressWarnings(suppressMessages(library(stringr, quietly = T)))
suppressWarnings(suppressMessages(library(lubridate, quietly = T)))

npr_df <- read_sas("/home2/genetics/antobe/data/STR/variables/ja_all_indiv_npr.sas7bdat")

RA_ICD10 <- c("M05", "M053", "M058", "M058A", "M058B", "M058C", "M058D", "M058F", "M058G", "M058H", "M058L", "M058M", "M058N", "M058X",
	      "M059", "M059A", "M059B", "M059C", "M059CD", "M059F", "M059G", "M059H", "M059L", "M059M", "M059N", "M059X",
	      "M06", "M060", "M060A", "M060B", "M060C", "M060D", "M060F", "M060G", "M060H", "M060L", "M060M", "M060N", "M060X",
	      "M068", "M068A", "M068B", "M068C", "M068D", "M068F", "M068G", "M068H", "M068L", "M068M", "M068N", "M068X",
	      "M069", "M069A", "M069B", "M069C", "M069D", "M069F", "M069G", "M069H", "M069L", "M069M", "M069N", "M069X")

npr_df_RA <- npr_df %>% filter(str_detect(hdia, str_c(RA_ICD10, collapse = "|"))) %>% mutate(GWAS_ID = paste0(bbPAIR, "_", bbTWIN)) %>% distinct(GWAS_ID)    #IDENTIFY ALL TO BE EXCLUDED

write.table(npr_df_RA, "TEMPORARY/tmp-1/STR_wRA.txt", quote = F, col.names = F, row.names = F, sep = "\t")

### TO DO:
### NOTES:
# 2.1. I've based the ICD-10 codes on the ones written down in the manuscript (version AS1). The ones in the manuscript are based on those found on internetmedicin.se. I've then checked the ones I've written down in the
#      manuscript towards the contents of the STR NPR file (all distinct ICD-10 with an "M") and made sure none are missing unintentionally. Some appear within the STR NPR file but have been intentionally excluded (Still's and
#      Felty's, along with some rheumatoid nodules and such).
# 2.2. I currently consider all patients with at least one occurrence of the above described ICD-10 codes in `hdia` as an RA-patient. One could incorporate information
#      from the secondary diagnoses too, or alternatively require multiple occurences of the relevant ICD-10 codes, but this seemed a simply enough filtration.
