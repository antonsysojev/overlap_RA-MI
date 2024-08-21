df <- read.xlsx("C:/Users/antsys/Desktop/Nikpay, Goel and Won et al - Supplementary Material - A comprehensive guide....xlsx", sheet = 2) %>% as_tibble()
df.sub <- df %>% slice(-c(1:3)) %>% slice(-(n():(n()-3))) %>% select(1, 2, 5, 6, 9, 10)
colnames(df.sub) <- c("ID1", "ID2", "ETHNICITY", "N_TOTAL", "N_CASE", "N_CTRL")

df.sub %>% mutate(MI_PROP = str_extract(N_CASE, "\\[.+\\]") %>% str_remove("\\[") %>% str_remove("\\]") %>% str_remove("\\%") %>% as.numeric()) %>%
           mutate(MI_PROP = ifelse(is.na(MI_PROP), 0, MI_PROP)) %>%
           mutate(N_CASE_CLEAN = str_remove(N_CASE, " \\[.+") %>% as.numeric()) %>%
           mutate(N_MI = N_CASE_CLEAN * (MI_PROP / 100)) %>%
           mutate(N_CTRL_CLEAN = str_remove(N_CTRL, "\\*") %>% str_remove(",") %>% str_remove("\\^") %>% as.numeric()) %>%
           mutate(N_MI_TOT = N_MI + N_CTRL_CLEAN) %>%
           group_by(ETHNICITY) %>% summarise(N_ETHNICITY = sum(N_MI_TOT)) %>% mutate(N_TOT = sum(N_ETHNICITY)) %>% mutate(N_PERC = N_ETHNICITY / N_TOT) %>%
           filter(!(ETHNICITY %in% c("White European", "white European", "white European (Finnish)", " white european"))) %>%
           summarise(TOT = sum(N_PERC))