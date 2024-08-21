### LAST VERSION UPDATE 14 JULY 2023 (v1.0)
### THIS SCRIPT CONTAINS A FUNCTION THAT PROCESSES THE LAVA OUTPUT

PROCESS_LAVA <- function(LAVA_LIST, PHENOTYPE, CHROM = NA, LOCUS_NO = NA){
  
  LOCUS <- c(CHROM, LOCUS_NO)
  UNIV <- data.frame(h2.obs_1 = NA, h2.obs_2 = NA, p_1 = NA, p_2 = NA)
  BIVAR <- data.frame(rho = NA, rho.lower = NA, rho.upper = NA, r2 = NA, r2.lower = NA, r2.upper = NA, p = NA)
  
  if(!is.null(LAVA_LIST[[1]])) UNIV <- LAVA_LIST[[1]] %>% filter(phen %in% c(PHENOTYPE, "MI")) %>% pivot_wider(names_from = phen, values_from = c(h2.obs, p))
  if(!is.null(LAVA_LIST[[2]])) BIVAR <- LAVA_LIST[[2]] %>% filter(phen1 == PHENOTYPE) %>% select(-phen1, -phen2)
  
  if(length(LAVA_LIST) == 3) LOCUS <- c(LAVA_LIST[[3]][1], LAVA_LIST[[3]][2])
  
  if(nrow(UNIV) == 0) UNIV <- data.frame(h2.obs_1 = NA, h2.obs_2 = NA, p_1 = NA, p_2 = NA)    #AD-HOC FIX TO A BUG
  if(nrow(BIVAR) == 0) BIVAR <- data.frame(rho = NA, rho.lower = NA, rho.upper = NA, r2 = NA, r2.lower = NA, r2.upper = NA, p = NA)
  
  data.frame(CHR = LOCUS[1], LOCUS = LOCUS[2]) %>% bind_cols(UNIV) %>% bind_cols(BIVAR)
  
}

### TO DO:
### NOTES:
# 2.1. If the univariate data is not null, but there are no univariate estimates for MI or your phenotype, then
#       we get a data.frame but of dim (0, 0) which will fuck up the data.frame in the end. The ad-hoc
#       solution fixes this, but is bad code.