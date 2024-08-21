### LAST VERSION UPDATE 11 JULY 2023 (v2.1) - NOW WITH ADDED COMMENTS...
### THIS SCRIPT PERFORMS THE LAVA LOCAL GENETIC CORRELATION ANALYSIS - BY RUNNING OUTSIDE OF THE LINUX SHELL

.libPaths("H:/Programs/RLibrary/")

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(LAVA)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(doParallel)))
suppressMessages(suppressWarnings(library(foreach)))

loci <- read.loci(str_c(RAW, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
input <- process.input(input.info.file = str_c(TMP, "RA-MI_InputInfo.txt"), sample.overlap.file = str_c(TMP, "RA-MI_SampleOverlap.txt"), ref.prefix = str_c(DATA, "REFPAN"), phenos = c("RA", "SPOS_RA", "SNEG_RA", "RA_DECODE", "MI"))

SubsetToRemainingLoci <- function(CHROM, LOCI){
  LAVA_CONTENTS <- dir(CACHE) %>% str_subset(str_c("LAVA_", CHROM, "_"))
  LAVA_RANGE_DF <- LAVA_CONTENTS %>% str_extract("\\d+\\-\\d+") %>% str_split_fixed("\\-", 2) %>% as.data.frame()
  COMPLETED_LOCI <- mapply(seq, from = LAVA_RANGE_DF$V1, to = LAVA_RANGE_DF$V2) %>% unlist(use.names = F) %>% sort() %>% unique()
  
  LOCI_SUB <- LOCI %>% filter(CHR == CHROM) %>% mutate(IDX = 1:nrow(.)) %>% filter(!(IDX %in% COMPLETED_LOCI))
  LOCI_SUB
}
loci_sub <- lapply(1:22, SubsetToRemainingLoci, loci) %>% bind_rows()
### THIS APPROACH IS FINE NORMALLY, BUT IT WILL NOT TAKE INTO ACCOUNT THE FILES CREATED VIA THE `LAVA_MISC` PROCEDURE...

cl <- makeCluster(30, outfile = "H:/Projects/OVERLAP_RA-CVD/TEMPORARY/tmp-4/OUTFILE_X.log")
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("H:/Programs/RLibrary/"))

checkpoints <- ceiling(quantile(1:nrow(loci_sub), seq(0.05, 1, 0.05)))
#RES_LIST <- list()

LAVA_RES <- foreach(i = 1:nrow(loci_sub), .combine = list, .multicombine = T, .packages = c("dplyr", "LAVA", "foreach")) %dopar% {

  #SOME INFORMATION HERE SO WE KNOW THAT WE STARTED LOOPING? OR NOT...
  
  LOCUS <- process.locus(loci_sub[i, ], input)
  
  UNIV_RES <- NULL; BIVAR_RES <- NULL; LOCUS_ID <- loci_sub[i, c("CHR", "LOC")]
  if(!is.null(LOCUS)){
    UNIV_RES <- run.univ(LOCUS)
    if("MI" %in% LOCUS$phenos & length(LOCUS$phenos) > 1){
      BIVAR_RES <- run.bivar(LOCUS, target = "MI")
    }
  }
  
  #RES_LIST[[i]] <- list(UNIV_RES, BIVAR_RES, LOCUS_ID)
  
  if(i %in% checkpoints){
    cat(paste0("LAVA STATUS UPDATE: LOCUS 1-", i, " PROCESSED, OUT OF THE TARGET ", nrow(loci_sub), "... \n"))
    #saveRDS(RES_LIST, paste0("H:/Projects/OVERLAP_RA-CVD/TEMPORARY/LAVA_MISC_1-", i, ".rds"))
    #SEE NOTE ABOUT THIS PROCEDURE
    #SEE NOTE ABOUT THE BUG
  }
  
  list(UNIV_RES, BIVAR_RES, LOCUS_ID)
  
}

stopCluster(cl)

### TO DO:
### NOTES:
# 1. There *might* be a problem with the current way of saving data. 
#     Suppose for instance that 30 parallel procedures are started at the same time, with a checkpoint at iteration i = 25.
#     Suppose then that iteration i = 25 finishes first: then it will save the .rds, which currently only
#     contains something at position 25.
#     This means that there MAY be holes in the output .rds, but this will have to be fixed at a later stage...
# 2. Checkpointing steps may be poor when starting from scratch. I've only got about 500 loci left which
#     means it hits at each 30th iteration almost which, assuming similar computational time on each node,
#     would mean one full go-through of a cluster. But more often may be better when starting from all 2495.
# 3. I removed the writing, which no longer works if I run it in parallel in the above way.
#     With the current version of the script (2.0) we save an RDS of `RES_LIST` at each checkpoint, but the
#     list only contains the contents of the current iteration. 
#     I'm not sure why this is happening, but I decided to circumvent it by running it for everything that remains
#     on a given chromosome, which is large enough that it should give me a reasonable output.
#     I can increase the size of the targets later.