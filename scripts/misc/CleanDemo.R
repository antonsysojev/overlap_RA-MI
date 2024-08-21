### LAST VERSION UPDATE 27 OCTOBER 2023 (v2.1) - FIXED A BUG THAT NOW GETS RID OF DUPLICATES
### CLEANS THE DEMOGRAPHICS FOR THE TABLE 1 AND THE TABLE S3

.libPaths("H:/Programs/RLibrary/")
library(dplyr)
library(readr)
library(stringr)

RAW.df <- read_tsv("H:/Projects/OVERLAP_RA-CVD/TEMPORARY/tmp-5/Demographics_RAW.txt")

### TABLE 1

RAW.df %>% distinct(IID, TYPE, RA) %>% group_by(TYPE, RA) %>% summarise(N = n())
RAW.df %>% distinct(IID, TYPE, RA, SEX) %>% group_by(TYPE, RA, SEX) %>% summarise(N = n()) %>% ungroup() %>% group_by(TYPE, RA) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(N_PERC = N / N_TOT)
RAW.df %>% distinct(IID, TYPE, RA, AGE) %>% group_by(TYPE, RA) %>% summarise(MU = mean(AGE, na.rm = T), SIGMA = sd(AGE, na.rm = T)) %>% ungroup()
RAW.df %>% distinct(IID, TYPE, RA, SEROSTATUS) %>% group_by(TYPE, RA, SEROSTATUS) %>% summarise(N = n()) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(N_PERC = N / N_TOT)
RAW.df %>% distinct(IID, TYPE, RA, debutalder) %>% group_by(TYPE, RA) %>% summarise(MU = mean(debutalder, na.rm = T), SIGMA = sd(debutalder, na.rm = T))
RAW.df %>% distinct(IID, TYPE, RA, duration_months) %>% group_by(TYPE, RA) %>% summarise(MU = mean(duration_months, na.rm = T), SIGMA = sd(duration_months, na.rm = T))

RAW.df %>% distinct(IID, RA) %>% group_by(RA) %>% summarise(N = n())
RAW.df %>% distinct(IID, RA, SEX) %>% group_by(RA, SEX) %>% summarise(N = n()) %>% ungroup() %>% group_by(RA) %>% mutate(N_TOT = sum(N)) %>% mutate(N_PERC = N / N_TOT)
RAW.df %>% distinct(IID, RA, AGE) %>% group_by(RA) %>% summarise(MU = mean(AGE, na.rm = T), SIGMA = sd(AGE, na.rm = T))

### TABLE S3

RAW.df %>% distinct(IID, TYPE, SEROSTATUS) %>% group_by(TYPE, SEROSTATUS) %>% summarise(N = n())
RAW.df %>% distinct(IID, TYPE, RA, SEROSTATUS, SEX) %>% group_by(TYPE, RA, SEROSTATUS, SEX) %>% summarise(N = n()) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(N_PERC = N / N_TOT)
RAW.df %>% distinct(IID, TYPE, RA, SEROSTATUS, AGE) %>% group_by(TYPE, RA, SEROSTATUS) %>% summarise(MU = mean(AGE, na.rm = T), SIGMA = sd(AGE, na.rm = T))
RAW.df %>% distinct(IID, TYPE, RA, SEROSTATUS, debutalder) %>% group_by(TYPE, RA, SEROSTATUS) %>% summarise(MU = mean(debutalder, na.rm = T), SIGMA = sd(debutalder, na.rm = T))
RAW.df %>% distinct(IID, TYPE, RA, SEROSTATUS, duration_months) %>% group_by(TYPE, RA, SEROSTATUS) %>% summarise(MU = mean(duration_months, na.rm = T), SIGMA = sd(duration_months, na.rm = T))

### NOTE SURE THE BELOW IS RELEVANT ANY MORE... (27 OCTOBER 2023)
### ### SEROPOSITIVE PART

RAW_SEROPOS.df <- RAW.df %>% filter(SEROSTATUS != "SERONEG RA")
RAW_SEROPOS.df %>% group_by(TYPE, RA) %>% group_by(TYPE, RA) %>% summarise(N = n())
RAW_SEROPOS.df %>% group_by(TYPE, RA, SEX) %>% summarise(N = n()) %>% ungroup() %>% group_by(TYPE, RA) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(N_PERC = N / N_TOT)
RAW_SEROPOS.df %>% group_by(TYPE, RA) %>% summarise(MU = mean(AGE, na.rm = T), SIGMA = sd(AGE, na.rm = T)) %>% ungroup()
RAW_SEROPOS.df %>% group_by(TYPE, RA, SEROSTATUS) %>% summarise(N = n()) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(N_PERC = N / N_TOT)
RAW_SEROPOS.df %>% group_by(TYPE, RA) %>% summarise(MU = mean(debutalder, na.rm = T), SIGMA = sd(debutalder, na.rm = T))
RAW_SEROPOS.df %>% group_by(TYPE, RA) %>% summarise(MU = mean(duration_months, na.rm = T), SIGMA = sd(duration_months, na.rm = T))

### ### SERONEGATIVE PART

RAW_SERONEG.df <- RAW.df %>% filter(SEROSTATUS != "SEROPOS RA")
RAW_SERONEG.df %>% group_by(TYPE, RA) %>% group_by(TYPE, RA) %>% summarise(N = n())
RAW_SERONEG.df %>% group_by(TYPE, RA, SEX) %>% summarise(N = n()) %>% ungroup() %>% group_by(TYPE, RA) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(N_PERC = N / N_TOT)
RAW_SERONEG.df %>% group_by(TYPE, RA) %>% summarise(MU = mean(AGE, na.rm = T), SIGMA = sd(AGE, na.rm = T)) %>% ungroup()
RAW_SERONEG.df %>% group_by(TYPE, RA, SEROSTATUS) %>% summarise(N = n()) %>% mutate(N_TOT = sum(N)) %>% ungroup() %>% mutate(N_PERC = N / N_TOT)
RAW_SERONEG.df %>% group_by(TYPE, RA) %>% summarise(MU = mean(debutalder, na.rm = T), SIGMA = sd(debutalder, na.rm = T))
RAW_SERONEG.df %>% group_by(TYPE, RA) %>% summarise(MU = mean(duration_months, na.rm = T), SIGMA = sd(duration_months, na.rm = T))
