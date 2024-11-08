### 2. Fix names in xlab to be appropriate.

### LAST VERSION UPDATE 9 JUNE (V1.2) - FIXED A BUG THAT MESSED UP MY FIGURE
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

clean_extra <- rbind(c("TRAIT", 0, "p-value", "", "", "", "", "", "p-value", "Genetic correlation"), clean)

### 2 - CONSTRUCT PLOT SECTIONS

MAIN_SECTION <- clean %>% ggplot(aes(y = factor(TYPE, levels = rev(TYPE)))) + 
  geom_point(aes(x = ESTIMATE), shape = 15, size = 3) +
  geom_linerange(aes(xmin = CI_L, xmax = CI_U)) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(0, nrow(clean)) + 1, xlim = c(-0.25, 0.5)) +
  labs(x = bquote(r[g]), y = "") +
  annotate("text", x = -.1, y = 12, label = "") +
  annotate("text", x = .1, y = 12, label = "") +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

LEFT_SECTION <- clean_extra %>% ggplot(aes(y = factor(TYPE, levels = rev(TYPE)))) + 
  geom_text(aes(x = 0, label = TYPE), hjust = 0, fontface = "bold") +
  geom_text(aes(x = 1, label = RESULTS_CLEAN), hjust = 0, fontface = ifelse(clean_extra$RESULTS_CLEAN == "Genetic correlation", "bold", "plain")) +
  coord_cartesian(xlim = c(0, 4)) + 
  theme_void()

RIGHT_SECTION <- clean_extra %>% ggplot() + 
  geom_text(aes(x = 0, y = factor(TYPE, levels = rev(TYPE)), label = P_CLEAN), hjust = 0, fontface = ifelse(clean_extra$P_CLEAN == "p-value", "bold", "plain")) + 
  theme_void()

layout <- c(area(t = 0, l = 0, b = 30, r = 3), area(t = 1, l = 4, b = 30, r = 9), area(t = 0, l = 9, b = 30, r = 11))

FULL_FIGURE <- LEFT_SECTION + MAIN_SECTION + RIGHT_SECTION + plot_layout(design = layout) 
FULL_FIGURE
ggsave("C:/Users/antsys/Desktop/FOREST.tiff", width = 12, height = 4)

### TO DO:
# 1.1. Update the input filepath to be given via `argparser` instead of as is done here...
# 1.2. Preferably, we would have some type of modified cleaning of data, that interprets the read data...§
# 1.3. Need to get the footer for the r_g in the x-axis label.
# 1.4. Range of the x-axis should be modified too, to go with what data contains...
### NOTES: