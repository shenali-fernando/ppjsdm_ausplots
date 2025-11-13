#########################################################
# Figures for publication

library(ggplot2)
library(dplyr)
library(ggpubr)


# Using df_all_3_5
df <- df_all_3_5

## Within-species interactions 
within <- df %>% 
  filter(species_from == species_to)

within <- within %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA))


within <- within %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA)) %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
    group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
    group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
    group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
    TRUE ~ group))

within$class_int <- factor(within$class_int,
                           levels = c("large_large", "small_large", "small_small" ))

###### Overall within-sp interactions 
#1..coloured
within %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) %>% 
  ggplot(aes(x = alpha, y = class_int, colour = cc_to)) + 
  geom_violin(scale = "width", 
              width = 0.65,
              position = position_dodge(0.7), 
              size = 0.75) +
  geom_boxplot(position = position_dodge(0.7),
               outlier.shape = NA, 
               width = 0.35, 
               size = 0.75) +
  geom_point(aes(colour = cc_to), 
             position = position_jitterdodge(seed = 21, 
                                           jitter.width = 0.3,
                                           dodge.width = 0.7), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.35) +
  scale_colour_manual(values = c("#35B779FF", "#39568CFF"), 
                      breaks = c("Subcanopy", "Canopy"),
                      name = "Strata") + 
  stat_summary(aes(fill = cc_to), fun = mean,
               geom = "point",
               shape = 19,
               size = 3, 
               color="black",
               position = position_dodge(width=0.7, preserve = "single"), 
               show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw(base_size = 14) + 
  labs(x = "Interaction Coefficient", y = "") +
  ylab("") + 
  ggtitle("Within-species interactions") 

#2. only points coloured
withinfig <- within %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) %>% 
  ggplot(aes(x = alpha, y = class_int, groups = cc_to)) + 
  geom_vline(xintercept = 0, colour = "red") +
  annotate("rect", 
           xmin = -0.1, xmax = 0.1, 
           ymin = -Inf, ymax = Inf,
           fill = "#EBB6B3") + 
  geom_violin(colour = "black", 
              scale = "width", 
              width = 0.7,
              fill = NA,
              position = position_dodge(0.75)) +
  geom_boxplot(colour = "black", position = position_dodge(0.75),
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  geom_point(aes(colour = cc_to), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.5) +
  scale_colour_manual(values = c("#35B779FF", "#39568CFF"), 
                      breaks = c("Subcanopy", "Canopy"),
                      name = "Strata") + 
  stat_summary(aes(fill = cc_to), fun = mean,
               geom = "point",
               shape = 19,
               size = 3, 
               color="black",
               position = position_dodge(width=0.7, preserve = "single"), 
               show.legend = FALSE) +
  theme_bw() + # Theme
  theme(text = element_text(size = 18), 
        legend.text=element_text(size= 12), 
        axis.title.x = element_text(margin = margin(8, 0, 0, 0))) + 
  guides(colour = guide_legend(override.aes = list(size=3.5))) + 
  labs(x = "Interaction Coefficient", y = "") +
  ylab("") + 
  ggtitle("Within-species interactions") 

ggsave("within.png", withinfig, 
       width = 11, height = 6.75)


########### To do ANOVA must check normality of data and equal variance 
summary <- within %>% group_by(cc_from, class_int) %>% 
  summarise(mean = mean(alpha), 
            count = n())


model <- aov(alpha ~ class_int * cc_to, data = within)
summary(model)

plot(model)         # residuals vs fitted and QQ plot
shapiro.test(residuals(model))   # normality
### Data is not normal!!!! p is less than 0.05
bartlett.test(alpha ~ interaction(class_int, cc_to), data = within)  # homogeneity
### Variance is not homogenous!!!!! 

## Need to use a non-parametric adjustment to the test
library(rstatix)

# Differences among groups within each subgroup
## We can use the Kruskal–Wallis Test 
within %>%
  group_by(cc_to) %>%
  kruskal_test(alpha ~ class_int)

# Or differences between canopy levels within each group
within %>%
  group_by(class_int) %>%
  kruskal_test(alpha ~ cc_to)

#Ad-hoc parwise comparisons using Pairwise Wilcox test 
  within %>%
  group_by(cc_to) %>%
  pairwise_wilcox_test(alpha ~ class_int, p.adjust.method = "BH")
    
  within %>%
    group_by(class_int) %>%
    pairwise_wilcox_test(alpha ~ cc_to, p.adjust.method = "BH")
  
    
within_sig <- withinfig +
geom_signif(comparisons = list(c("Small ↔ Small", "Large ↔ Large"), 
                              c("Small ↔ Large", "Large ↔ Large")), 
                        map_signif_level = T, 
                        y_position = c(1.25, 1.75)) + 
  geom_signif(y_position = c(0.5, 0.75), 
              xmin = c(0.8, 2.8), 
              xmax = c(1.2, 3.2),
              annotation = c("**", "*"), 
              tip_length = 0)  


ggsave("withinsig.png", within_sig, 
       width = 14, height = 8, 
       dpi = 300)






#################################################################################################
#################################################################################################
################## BETWEEN INTS #######################

between <- df %>% 
  filter(!(species_from == species_to)) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int))


between <- between %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA)) %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
    group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
    group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
    group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
    TRUE ~ group))


#split into two parts 
bwfig <- between %>% 
  filter(group %in% c("Subcanopy.small_Subcanopy.small", 
                      "Subcanopy.large_Subcanopy.large", 
                      "Canopy.large_Canopy.large", 
                      "Canopy.small_Canopy.small")) %>%  
  mutate(group_int = case_when(
    group == "Subcanopy.small_Subcanopy.small"~ "Subcanopy small ↔ Subcanopy small", 
    group == "Subcanopy.large_Subcanopy.large" ~ "Subcanopy large ↔ Subcanopy large", 
    group == "Canopy.large_Canopy.large" ~ "Canopy large ↔ Canopy large", 
    group == "Canopy.small_Canopy.small" ~ "Canopy small ↔ Canopy small")) %>% 
  ggplot(aes(x = alpha, y = group_int)) + 
  annotate("rect", 
           xmin = -0.1, xmax = 0.1, 
           ymin = -Inf, ymax = Inf,
           fill = "#EBB6B3") + 
  geom_vline(xintercept = 0, colour = "red") + 
  geom_violin(colour = "black", 
              scale = "width", 
              width = 0.7,
              fill = NA) +
  geom_boxplot(colour = "black",
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  geom_point(aes(colour = region), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.1), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.5) + 
  stat_summary(aes(fill = group_int),
               fun = mean,
               geom = "point",
               shape = 19,
               size = 3, 
               color="black",
               position = position_dodge(width=0.7, preserve = "single"), 
               show.legend = FALSE) +
  theme_bw() + # Theme
  theme(text = element_text(size = 18), 
        legend.text=element_text(size= 12), 
        axis.title.x = element_text(margin = margin(8, 0, 0, 0))) + 
  guides(colour = guide_legend(override.aes = list(size=3.5))) + 
  labs(x = "Interaction Coefficient", y = "") +
  ylab("") + 
  ggtitle("Between-species interactions") 

bwfig

ggsave("bwfig.png", bwfig, 
       width = 14, height = 8, 
       dpi = 300)
########################################
bwfig1 <- between %>% 
  filter(group %in% c("Canopy.small_Subcanopy.small", 
                      "Canopy.small_Subcanopy.large",
                      "Canopy.large_Subcanopy.small", 
                      "Canopy.large_Subcanopy.large")) %>%   
  mutate(group_int = case_when(
    group == "Canopy.small_Subcanopy.small"~ "Canopy small ↔ Subcanopy small", 
    group == "Canopy.small_Subcanopy.large" ~ "Canopy small ↔ Subcanopy large", 
    group == "Canopy.large_Subcanopy.small" ~ "Canopy large ↔ Subcanopy small", 
    group == "Canopy.large_Subcanopy.large" ~ "Canopy large ↔ Subcanopy large")) %>% 
  ggplot(aes(x = alpha, y = group_int)) + 
  annotate("rect", 
           xmin = -0.1, xmax = 0.1, 
           ymin = -Inf, ymax = Inf,
           fill = "#EBB6B3") + 
  geom_vline(xintercept = 0, colour = "red") + 
  geom_violin(colour = "black", 
              scale = "width", 
              width = 0.7,
              fill = NA) +
  geom_boxplot(colour = "black",
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  geom_point(aes(colour = region), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.1), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.5) + 
  stat_summary(aes(fill = group_int),
               fun = mean,
               geom = "point",
               shape = 19,
               size = 3, 
               color="black",
               position = position_dodge(width=0.7, preserve = "single"), 
               show.legend = FALSE) +
  theme_bw() + # Theme
  theme(text = element_text(size = 18), 
        legend.text=element_text(size= 12), 
        axis.title.x = element_text(margin = margin(8, 0, 0, 0))) + 
  guides(colour = guide_legend(override.aes = list(size=3.5))) + 
  labs(x = "Interaction Coefficient", y = "") +
  ylab("") + 
  ggtitle("Between-species interactions") 

bwfig1


ggsave("bwfig1.png", bwfig1, 
       width = 14, height = 8, 
       dpi = 300)


model <- aov(alpha ~ group, data = between)
summary(model)

plot(model)         # residuals vs fitted and QQ plot
shapiro.test(residuals(model))   # normality
### Data is not normal!!!! p is less than 0.05
bartlett.test(alpha ~ interaction(class_int, cc_to), data = within)  # homogeneity
### Variance is not homogenous!!!!! 


## We can use the Kruskal–Wallis Test 
between %>%
  kruskal_test(alpha ~ group)

#Ad-hoc parwise comparisons using Pairwise Wilcox test 
wilcox <- between %>%
  pairwise_wilcox_test(alpha ~ group, p.adjust.method = "BH")


