### LAST VERSION UPDATE 16 OCT 2023 (V1.2) - CHANGED THE FORMAT UP
### THIS SCRIPT CREATES A NEAT FOREST PLOT IN LINE WITH THIS ONE: https://www.khstats.com/blog/forest-plots/

library(tidyverse)    #AT LEAST `stringr`, `readr`, `dplyr`
library(ggplot2)     #IS IT NOT PART OF `tidyverse`? IF SO, CAN BE SKIPPED HERE
library(patchwork)

raw <- read_tsv("C:/Users/antsys/Desktop/MAIN/Projects/OVERLAP_RA-CVD/TMP/LDSC-RES.txt", show_col_types = F)

### 1 - CLEAN DATA

clean <- raw %>% select(TYPE, ESTIMATE = RG, SERROR = RG_SE, P = P) %>%
  mutate(CI_L = ESTIMATE - qnorm(1 - 0.05 / 2) * SERROR, CI_U = ESTIMATE + qnorm(1 - 0.05 / 2) * SERROR) %>% select(-SERROR) %>%
  mutate(ESTIMATE_CLEAN = round(ESTIMATE, 2) %>% format(nsmall = 2), CI_L_CLEAN = round(CI_L, 2) %>% format(nsmall = 2), CI_U_CLEAN = round(CI_U, 2) %>% format(nsmall = 2), P_CLEAN = round(P, 2) %>% format(nsmall = 2)) %>%
  mutate(RESULTS_CLEAN = str_c(ESTIMATE_CLEAN, " (", CI_L_CLEAN, "-", CI_U_CLEAN, ")"))

clean_extra <- rbind(c("TRAIT", NA, "p-value", 0, 0, "", 0, 0, "p-value", "Genetic correlation"), clean) %>% mutate(ESTIMATE = as.numeric(ESTIMATE), CI_L = as.numeric(CI_L), CI_U = as.numeric(CI_U))
clean_extra$TYPE <- c("TRAIT", "Myocardial Infarction", "C-Reactive Protein", "Body Mass Index", "Waist-to-hip-ratio", "Systolic Blood Pressure", "Diastolic Blood Pressure", "Pulse Pressure", "Low-density Lipoprotein Cholesterol", "Smoking (Cigarettes per day)", "Smoking (ever/never)")

### 2 - CONSTRUCT PLOT SECTIONS

MAIN_SECTION <- clean_extra %>% ggplot(aes(y = factor(TYPE, levels = rev(TYPE)))) + 
  geom_point(aes(x = ESTIMATE), shape = 15, size = 3) +
  geom_linerange(aes(xmin = CI_L, xmax = CI_U)) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(1, nrow(clean_extra)), xlim = c(-0.15, 0.3)) +
  labs(x = bquote(r[g]), y = "") +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

LEFT_SECTION <- clean_extra %>% ggplot(aes(y = factor(TYPE, levels = rev(TYPE)))) + 
  geom_text(aes(x = 0, label = TYPE), hjust = 0, fontface = "bold", size = 5) +
  #geom_text(aes(x = 1, label = RESULTS_CLEAN), hjust = 0, fontface = ifelse(clean_extra$RESULTS_CLEAN == "Genetic correlation", "bold", "plain"), size = 5) +
  coord_cartesian(xlim = c(0, 4)) + 
  theme_void()

RIGHT_SECTION <- clean_extra %>% ggplot() + 
  geom_text(aes(x = 0, y = factor(TYPE, levels = rev(TYPE)), label = RESULTS_CLEAN), hjust = 0, fontface = ifelse(clean_extra$RESULTS_CLEAN == "Genetic correlation", "bold", "plain"), size = 5) +
  #geom_text(aes(x = 0, y = factor(TYPE, levels = rev(TYPE)), label = P_CLEAN), hjust = 0, fontface = ifelse(clean_extra$P_CLEAN == "p-value", "bold", "plain"), size = 5) + 
  theme_void()

layout <- c(area(t = 0, l = 0, b = 30, r = 8), area(t = 1, l = 8, b = 30, r = 15), area(t = 0, l = 11, b = 30, r = 20))
#layout <- c(area(t = 0, l = 0, b = 30, r = 4), area(t = 1, l = 4, b = 30, r = 11))

FULL_FIGURE <- LEFT_SECTION + MAIN_SECTION + RIGHT_SECTION + plot_layout(design = layout) 
#FULL_FIGURE <- LEFT_SECTION + MAIN_SECTION + plot_layout(design = layout)
FULL_FIGURE
ggsave("C:/Users/antsys/Desktop/FOREST.tiff", width = 3225, height = 1075, units = "px")

### TO DO:
# 1.1. Update the input filepath to be given via `argparser` instead of as is done here...
# 1.2. Preferably, we would have some type of modified cleaning of data, that interprets the read data...§
# 1.3. Need to get the footer for the r_g in the x-axis label.
# 1.4. Range of the x-axis should be modified too, to go with what data contains...
### NOTES: