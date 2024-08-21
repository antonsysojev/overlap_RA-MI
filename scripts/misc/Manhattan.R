#!/user/bin/env bash
### LAST VERSION UPDATE - 21 FEB 2023
### THIS SCRIPT PRODUCES A MANHATTAN PLOT FOR THE INPUT GWAS DATA

.libPaths("/home2/genetics/antobe/software/RLibrary")
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(argparser)))

args <- arg_parser('EMPTY FOR NOW') %>% add_argument('GWAS', help = 'FILEPATH TO THE GWAS DATA FOR WHICH THE MANHATTAN SHOULD BE CONSTRUCTED', type = 'character') %>% parse_args()
GWAS <- args$GWAS
PHENO <- str_remove_all(GWAS, "^.+\\/") %>% str_remove_all("\\..+$")
GWAS_df <- read_tsv(GWAS, col_types = c("dcdccddddddddd"), show_col_types = F)

GWAS_CLEAN <- GWAS_df %>% inner_join(GWAS_df %>% group_by(CHR) %>% summarise(MAX_POS = as.numeric(max(POS))) %>% mutate(POS_ADD = lag(cumsum(MAX_POS), default = 0)) %>% select(CHR, POS_ADD), by = "CHR") %>%    #BE CAREFUL ABOUT POSSIBLE INTEGER OVERFLOW HERE...
	mutate(CUMUL_POS = POS + POS_ADD, COL_HELPER = if_else(CHR %% 2 == 0, 1, 0))
CHR_MED_POS <- GWAS_CLEAN %>% group_by(CHR) %>% summarise(MEDIAN_POS = round(median(CUMUL_POS), 0))    #HELPER DATA FOR X-AXIS

tiff(filename = paste0("DATA/OUTPUT/DATA/", PHENO, ".tiff"), width = 1400, height = 600, units = "px")
ggplot(GWAS_CLEAN, aes(x = CUMUL_POS, y = -log10(P), col = factor(COL_HELPER))) +
	geom_point(alpha = 0.75) +
	scale_x_continuous(breaks = CHR_MED_POS$MEDIAN_POS, labels = CHR_MED_POS$CHR) %>%
	scale_y_continuous(expand = c(0, 0), limits = c(0, -log10(5e-9))) %>%
	scale_color_manual(values = c("Grey 5", "Grey 40")) +
	xlab("") + ylab("-log10(p)") +
	theme_minimal() +
	theme(legend.position = "none",
	      panel.border = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.grid.minor.x = element_blank())
dev.off() %>% invisible()   #THE `%>% invisible()` SHUTS `dev.off()` UP ACTUALLY!

### TO DO:
# 1.1. Can we make the pixels controlable? I think 1400 x 600 is fine, but it might be nicer to make these smaller directly in R, instead of messing around in paint or something.
# 1.2. You should consider using `ggsave()` insted, which seems more stable.
### NOTES:
# 2.1.
