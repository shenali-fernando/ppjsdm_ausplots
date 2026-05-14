#########################################################
# Figures for publication
library(forcats)
library(ggbeeswarm)
library(ggplot2)
library(dplyr)
library(ggpubr)


# Using df_all_3_5
df <- df_add3_5_t15
# 
# species_class <- species_class %>% select(Species, edge_case)
# species_class <- species_class %>% rename(species_to = species_from)
# 
# df <- left_join(df, species_class, by = "species_to")
# df <- df %>% rename(edge_species_to = edge_case)

# within <- within %>%
#   mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA)) %>%
#   mutate(group = paste0(edge_species_from, ".", class_from,  "_", edge_species_to, ".", class_to)) %>%
#   mutate(group = case_when(
#     group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small",
#     group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large",
#     group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small",
#     group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large",
#     group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small",
#     group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
#     TRUE ~ group))


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


###### Overall within-sp interactions 

within <- within %>% 
  mutate(alpha2 = ifelse(alpha < -2.25, -2.33, alpha))

withinfig <- within %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) %>% 
  ggplot(aes(x = alpha2, y = class_int, groups = cc_to)) + 
   geom_vline(xintercept = 0, colour = "red") +
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
  geom_text(aes(label = ifelse(alpha2 == -2.33, "-3.3", "")), 
            hjust = 1,
            vjust = 2.5, 
            size = 2.25) +
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

withinfig

ggsave("within.png", withinfig, 
       width = 11, height = 6.75)



#weighted mean by se or cis = ave ints
#lm weights, no predictor mod, 

### Same fig but with weighted means instead
#Adjust weighted means df 
within_wmean <- w1 %>% 
  mutate(group_1 = c("Canopy", "Canopy", "Canopy", 
                     "Subcanopy", "Subcanopy", "Subcanopy")) %>% 
  mutate(class_int = c("Large ↔ Large", "Small ↔ Large", "Small ↔ Small",
                       "Large ↔ Large", "Small ↔ Large", "Small ↔ Small")) %>% 
  mutate(group_1 = fct_recode())

withinfig <- within %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) %>% 
  ggplot(aes(x = alpha2, y = fct_rev(class_int), groups = fct_rev(cc_to))) + 
  geom_vline(xintercept = 0, colour = "red") +
  geom_violin(colour = "black", 
              scale = "width", 
              width = 0.7,
              fill = NA,
              position = position_dodge(0.75)) +
  geom_boxplot(colour = "black", position = position_dodge(0.75),
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  geom_point(aes(colour = fct_rev(cc_to)), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.5) +
  scale_colour_manual(values = c("#39568CFF", "#35B779FF"), 
                      labels = c("Canopy", "Subcanopy"),
                      limits = c("Canopy", "Subcanopy"),
                      name = "Strata") + 
  geom_text(aes(label = ifelse(alpha2 == -2.33, "-3.3", "")), 
            hjust = 1,
            vjust = -6.5, 
            size = 2.25) +
  geom_point(data = within_wmean, 
             aes(x = `w.mean`, y = class_int, groups = fct_rev(group_1)),
             position = position_dodge(0.75), 
             size = 2.75) +
  theme_bw() + # Theme
  theme(text = element_text(size = 15), 
        legend.title = element_text(size = 14),
        legend.text=element_text(size= 12), 
        axis.title.x = element_text(margin = margin(8, 0, 0, 0), 
                                    size = 12)) + 
  guides(colour = guide_legend(override.aes = list(size=3.5))) + 
  labs(x = "Interaction Coefficient", y = "") 

withinfig






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
  geom_signif(y_position = c(0.5), 
              xmin = c(2.8), 
              xmax = c(3.2),
              annotation = c("**"), 
              tip_length = 0)  

within_sig
ggsave("withinsig.png", within_sig, 
       dpi = 300)



###### Lines showing change fit 

#make median df 
median_df <- within %>% 
  group_by(cc_to, class_int) %>% 
  summarise(median = median(alpha)) %>% 
  ungroup() %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) 

mean_df <- within %>% 
  group_by(cc_to, class_int) %>% 
  summarise(mean = mean(alpha)) %>% 
  ungroup() %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) 



lines <- within %>% 
  mutate(species_site = paste0(species_from, sep = "_", site)) %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  arrange(species_site, class_int) %>% 
  ggplot() + 
  geom_hline(yintercept = 0, colour = "red") + 
  geom_beeswarm(aes(x = class_int, y = alpha), 
                alpha = 0.35) + 
  geom_line(aes(x = class_int, y = alpha, group = species_site), alpha = 0.1, na.rm = T) + 
  geom_point(data = median_df, aes(x = class_int, y = median),
             colour = "#EB7C0E", size = 2) +
  geom_point(data = mean_df, aes(x = class_int, y = mean),
             colour = "#6BF5F7", size = 2) +
  facet_wrap(~cc_to, scales = 'free_y') +
  xlab("") +
  ylab("Alpha coefficient") +
  theme_bw() + 
  theme(axis.title.y = element_text(size = 10))

lines

ggsave("lines.png", lines, 
       dpi =700)


####### Within scatterplots 

a <- within %>% 
  pivot_wider(id_cols = c(species_from, georegion, site, cc_to), 
              names_from = class_int, 
              values_from = alpha)

a <- a %>% 
  rename("Small ↔ Small" = "small_small") %>% 
  rename("Small ↔ Large" = "small_large") %>% 
  rename("Large ↔ Large" = "large_large")


### Small-small :: large-large

ggplot(data = a, 
         aes(x = `Small ↔ Small`,
             y = `Large ↔ Large`, 
             colour = cc_to)) + 
    geom_point(shape = 19, size = 2) + 
  geom_smooth(method=lm, se=FALSE) +
  theme_bw() + 
  ylim(c(-1.75, 0.5))


### small-small :: small-large 

ggplot(data = a, 
       aes(y = `Small ↔ Small`,
           x = `Small ↔ Large`, 
           colour = cc_to)) + 
  geom_point(shape = 19, size = 2) + 
  geom_smooth(method=lm, se=FALSE) +
  theme_bw()



### large-large :: small-large
ggplot(data = a, 
       aes(x = `Small ↔ Large`,
           y = `Large ↔ Large`, 
           colour = cc_to)) + 
  geom_point(shape = 19, size = 2) + 
  geom_smooth(method=lm, se=FALSE) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)+
  theme_bw() + 
  ylim(c(-1.75, 0.5))



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
           xmin = -0.05, xmax = 0.05, 
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


##Use weighted means 

bw_mean <- b1 %>% 
  filter(group %in% c("Subcanopy.small_Subcanopy.small", 
                      "Subcanopy.large_Subcanopy.large", 
                      "Canopy.large_Canopy.large", 
                      "Canopy.small_Canopy.small")) %>%  
  mutate(group_int = case_when(
    group == "Subcanopy.small_Subcanopy.small"~ "Subcanopy small ↔ Subcanopy small", 
    group == "Subcanopy.large_Subcanopy.large" ~ "Subcanopy large ↔ Subcanopy large", 
    group == "Canopy.large_Canopy.large" ~ "Canopy large ↔ Canopy large", 
    group == "Canopy.small_Canopy.small" ~ "Canopy small ↔ Canopy small"))


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
  geom_vline(xintercept = 0, colour = "red") +
  # annotate("rect", 
  #          xmin = -0.05, xmax = 0.05, 
  #          ymin = -Inf, ymax = Inf,
  #          fill = "#EBB6B3") + 
  # geom_vline(xintercept = 0, colour = "red") + 
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
             alpha = 0.45) + 
  geom_point(data = bw_mean, 
             aes(x = `w.mean`, y = group_int),
             position = position_dodge(0.75), 
             size = 2.5) +
  theme_bw(base_size = 10) + 
  scale_color_discrete(labels = c("Northeast Aus", "Southeast Aus", "Tasmania", "Western Aus"),
                       name = "Region") + 
  theme(axis.title.x = element_text(margin = margin(8, 0, 0, 0)), 
                title = element_text(size = 9), 
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10)) +  
  guides(colour = guide_legend(override.aes = list(size=3.5))) + 
  labs(x = "", y = "") + 
  ggtitle("a. Between-species, within functional group")

bwfig


ggsave("bwbwfig", bwfig, 
       dpi = 300)



########################################

bw_mean2 <- b1 %>% 
  filter(group %in% c("Canopy.small_Subcanopy.small", 
                      "Canopy.small_Subcanopy.large",
                      "Canopy.large_Subcanopy.small", 
                      "Canopy.large_Subcanopy.large")) %>%  
  mutate(group_int = case_when(
    group == "Canopy.small_Subcanopy.small"~ "Canopy small ↔ Subcanopy small", 
    group == "Canopy.small_Subcanopy.large" ~ "Canopy small ↔ Subcanopy large", 
    group == "Canopy.large_Subcanopy.small" ~ "Canopy large ↔ Subcanopy small", 
    group == "Canopy.large_Subcanopy.large" ~ "Canopy large ↔ Subcanopy large"))




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
             alpha = 0.35) + 
  geom_point(data = bw_mean2, 
             aes(x = `w.mean`, y = group_int),
             position = position_dodge(0.75), 
             size = 2.5) +
  theme_bw(base_size = 10) + # Theme
  scale_color_discrete(guide = "none") +
  theme(axis.title.x = element_text(margin = margin(8, 0, 0, 0), size = 10), 
        axis.title.y = element_text(size = 10),
        title = element_text(size = 9)) + 
  #guides(colour = guide_legend(override.aes = list(size=3.5))) + 
  labs(x = "Interaction Coefficient", y = "") + 
  ggtitle("b. Between-species, between functional group")

bwfig1


# ggsave("bwfig1.png", bwfig1, 
#        width = 14, height = 8, 
#        dpi = 300)

##add together
bw <- bwfig / bwfig1 + plot_layout(guides = "collect")
bw

ggsave("betweenfig2.png", bw, 
       width = 8,height = 8,
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


