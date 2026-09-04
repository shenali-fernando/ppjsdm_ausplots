library(forcats)
library(ggbeeswarm)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)


# Read in model specification dfs 
#or updated sp_size model 
sp_size <- read.csv("../../scripts/model_specifications/species_diameter/sp_size_df_t15_misc_updated.csv")
#and fg_size model 
fg_size <- read.csv("../../scripts/model_specifications/fg_size/fg_size_df_t15_updated.csv")


#Also need to read in the linear model predicted coefficients 
#for sp_size
#Read in the linear model post-hoc intraspecific coefficients 
sp_preds <- read.csv("scripts/analysis/lm_sp_size_intra_updated.csv")

#and for df_size
fg_preds <- read.csv("scripts/analysis/lm_fg_size_intra.csv")

#######################################################################################
############## SCATTERPLOT 
## The easiest visual comparison of the two model is a scatterplot of a-priori (fg_size)
#  v the post-hoc (linear model of the species_size) models. There will be the 6 classes 
# for each, and we can see where these line up on the one-to-one line 


#first clean
intra_fg_size <- fg_size %>% 
  filter(class_from == class_to)  %>% 
  mutate(size_int = case_when(size_int == "small small" ~ "Small ↔ Small", 
                               size_int == "small large" ~ "Small ↔ Large", 
                               size_int == "large large" ~ "Large ↔ Large"))
  

#then summarise
intra_fg_size_sum <- intra_fg_size %>% 
  group_by(size_int, class_from) %>% 
  summarise(median = median(alpha)) %>% 
  ungroup()

#cleaning to bind 
intra_fg_size_sum <- intra_fg_size_sum %>% 
  mutate(model = "a_priori")

intra_sp_preds <- intra_sp_preds %>% 
  select(x, group, predicted) %>% 
  mutate(model = "post_hoc") %>% 
  rename(median = predicted) %>% 
  rename(size_int = x) %>% 
  rename(class_from = group)

intra_sum_df <- rbind(intra_fg_size_sum, intra_sp_preds)

#oh dang, need to pivot_wider

intra_sum_df_wide <- pivot_wider(data = intra_sum_df, 
                                 names_from = model,
                                 values_from = median)

## then plot 
ggplot(intra_sum_df_wide, 
       aes(x = post_hoc, 
           y = a_priori)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 0, color = "gray70") + 
  geom_hline(yintercept = 0, color = "gray70") + 
  theme_bw() + 
  ylab("a priori (fg+size model)") + 
  xlab("post-hoc (linear model of species+size model)") + 
  xlim(c(-0.6, 0.3)) + 
  ylim(c(-0.6, 0.3))


#to add colour and shape to different types of points take away geom_smooth
ggplot(intra_sum_df_wide, 
       aes(x = post_hoc, 
           y = a_priori, 
           colour = class_from,
           shape = size_int)) + 
  geom_point(size = 4.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 0, color = "gray70") + 
  geom_hline(yintercept = 0, color = "gray70") + 
  theme_bw() + 
  ylab("a priori (fg+size model)") + 
  xlab("post-hoc (linear model of species+size model)") + 
  xlim(c(-0.6, 0.3)) + 
  ylim(c(-0.6, 0.3)) + 
  labs(
    color = "FG",
    shape = "Size Interaction"
  )


## okay, so there actually are some major differences in how the two models are working
## the sp_size model effectively doubles the interaction between L-L individuals, 
# this is particular the case for subcanopy individuals which have a 0!!! interaction 
# when in the fg_size model. 
# however, the fg_size model has a stronger clustering effect in the S-S interaction, 
# which means that the 0 effect in species_size for canopy is 0.2 in fg_size
# this is not goodd 

## so major differences we know are occurring between the two model 
# we know that the aggregation effect of the fg model is quite large 
# we need to now understand WHY there is such a difference 









##########################################################################
############### RAW DATA COMPARISON 

#Let's fix a few things 
fg_size$model <- "fg_size"
sp_size$model <- "sp_size"

# get rid of misc in sp_size 
sp_size <- sp_size %>% 
  filter(! from %in% c("Misc_large", "Misc_small")) %>% 
  filter(! to %in% c("Misc_large", "Misc_small"))

#bind 
full_fg_sp <- bind_rows(fg_size, sp_size)

#Need to neaten up some columns for plotting 
full_fg_sp <- full_fg_sp %>% 
  mutate(size_int = case_when(
    size_int %in% c("small small", "small_small") ~ "Small ↔ Small", 
    size_int %in% c("small large", "large_small", "small_large") ~ "Small ↔ Large", 
    size_int %in% c("large large", "large_large") ~ "Large ↔ Large", 
    TRUE ~ size_int))

full_fg_sp <- full_fg_sp %>% 
  mutate(size_int = as.factor(size_int)) %>% 
  mutate(size_int = fct_relevel(size_int, 
                                "Small ↔ Small",
                                "Small ↔ Large", 
                                "Large ↔ Large")) %>% 
  mutate(model = case_when(model == "fg_size" ~ "FG + Size", 
                           model == "sp_size" ~ "Species + Size"))


###########INTRA
#filter for within-group 
intra <- full_fg_sp %>% 
  filter(
    (model == "Species + Size" & (species_from == species_to)) | 
    (model == "FG + Size" & (class_from == class_to))) 

intra <- intra %>% 
  mutate(class_from = as.factor(class_from)) %>% 
  mutate(class_from = fct_relevel(class_from, "Subcanopy", "Canopy"))

#and plot 
ggplot(data = intra,
       aes(x = alpha, 
           y = size_int, 
           groups = class_from)) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  geom_violin(colour = "black", 
              scale = "width", 
              width = 0.7,
              fill = NA,
              position = position_dodge(0.75)) +
  geom_boxplot(colour = "black", position = position_dodge(0.75),
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  geom_point(aes(colour = class_from), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.25) +
  scale_colour_manual(values = c("Canopy" = "#39568CFF", 
                                 "Subcanopy" = "#35B779FF")) + 
  guides(colour = guide_legend(reverse = TRUE)) +
  facet_wrap(~model, ncol = 1) + 
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient")


############### Instead of boxplot and violin, add the predicted lm model instead 
#do some adding of cols and bind 
fg_preds$model <- "FG + Size"
sp_preds$model <- "Species + Size"

preds <- rbind(fg_preds, sp_preds)

preds <- preds %>% rename(class_from = group)

#and order everything 
preds <- preds %>% 
  mutate(class_from = as.factor(class_from)) %>% 
  mutate(class_from = fct_relevel(class_from, "Subcanopy", "Canopy")) %>% 
  mutate(x = as.factor(x)) %>% 
  mutate(x = fct_relevel(x, 
                                "Small ↔ Small",
                                "Small ↔ Large", 
                                "Large ↔ Large"))

#plot in 
fig2 <- ggplot() +  
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.2, linetype = "dashed") +
  geom_violin(data = intra,
              aes(x = alpha, 
                  y = size_int,
                  groups = class_from),
              scale = "width",
              colour = "gray50",
              width = 0.7,
              fill = NA,
              position = position_dodge(0.75)) +  
  geom_point(data = intra,
             aes(x = alpha, 
                 y = size_int,
                 groups = class_from,
                 colour = class_from),
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.08) +
  geom_point(data = preds,
             aes(y = x,
                 x = predicted,
                 group = class_from,
                 colour = class_from),
             size = 1.7,
             position = position_jitterdodge(seed = 21,
                                             jitter.width = 0.05,
                                             dodge.width = 0.75)) +
  geom_linerange(data = preds,
                aes(y = x,
                    x = predicted, 
                    xmin = conf.low,
                    xmax = conf.high,
                    colour = class_from),
                position = position_jitterdodge(seed = 21, 
                                                jitter.width = 0.05,
                                                dodge.width = 0.75),
                linewidth = 1) +
  # annotate("segment", 
  #          x = 0.45, 
  #          xend = 0.45, 
  #          y = 2.75,
  #          yend = 3.25) + 
  # annotate("segment", 
  #          x = 1.4,
  #          xend = 1.4, 
  #          y = 1,
  #          yend = 3) + 
  scale_colour_manual(values = c("Canopy" = "#440154FF", 
                                 "Subcanopy" =  "#1FA187FF"), 
                      name = "Functional Group") + 
  guides(colour = guide_legend(reverse = TRUE)) +
  facet_wrap(~model, ncol = 1) + 
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient") + 
  theme(text = element_text(size = 12), 
        legend.title = element_text(size = 11, margin = margin(0, 20, 1, 0)),
        legend.text=element_text(size= 10, margin = margin(r = 13, l = 7)), 
        legend.position = "bottom", 
        axis.title.x = element_text(margin = margin(8, 0, 0, 0), 
                                    size = 12), 
        strip.background = element_rect(fill = "gray95"))

fig2


seg_data <- data.frame(
  x = c(0.45, 1.55, 0.45, 1.5, 1.4, 0.9, 1.1), 
  xend = c(0.45, 1.55, 0.45, 1.5, 1.4, 0.9, 1.1), 
  y = c(2.75, 1, 2.75, 1, 2, 1, 2.05), 
  yend = c(3.25, 3, 3.25, 3, 3, 1.95, 3.25), 
  model = c("Species + Size", "Species + Size", "FG + Size", "FG + Size", "Species + Size", "FG + Size", "FG + Size")
)

fig2 + geom_segment(
  data = seg_data, 
  mapping = aes(x = x, y = y, xend = xend, yend = yend), 
  linewidth = 0.25
)

#save out and add segment ends in powerpoint? 

ggsave("fig2.png", 
       width = 13, 
       height = 14,
       scale = 1.5,
       units = "cm",
       dpi = 300)



#There is more clustering in the fg + size model in the L-L and S-S groupings 
#So, there is more repulsion in the species + size model 
# why? 
# there are 6 groups at most in each FG + Size models plots. 

# ** TO DO ** 
#Would be quite nice to have a way to qualtify or look qualitatively at what is 
#happening at the grouping level...



###########################################################################
##### BETWEEN GROUP COMPARISON (RAW DATA)

inter <- full_fg_sp %>% 
  filter(
    (model == "Species + Size" & (!species_from == species_to)) | 
      (model == "FG + Size" & (!class_from == class_to))) 

a <- setdiff(inter, intra) #all rows are unique, that's good

#make the int column and do some cleaning 
inter <- inter %>% 
  select(-int) %>% 
  mutate(int = paste0(class_from, " ", size_from, "_", class_to, " ", size_to)) %>% 
  mutate(int = case_when(
        int %in% c("Canopy large_Canopy small", "Canopy small_Canopy large") ~ "Canopy large_Canopy small",
        int %in% c("Canopy large_Subcanopy large", "Subcanopy large_Canopy large", "Subcanopy  large_Canopy  large") ~ "Canopy large_Subcanopy large",
        int %in% c("Canopy large_Subcanopy small", "Subcanopy small_Canopy large", "Subcanopy  small_Canopy  large") ~  "Canopy large_Subcanopy small",
        int %in% c("Canopy small_Subcanopy large", "Subcanopy large_Canopy small", "Subcanopy  large_Canopy  small") ~ "Canopy small_Subcanopy large",
        int %in% c("Canopy small_Subcanopy small", "Subcanopy small_Canopy small", "Subcanopy  small_Canopy  small") ~ "Canopy small_Subcanopy small",
        int %in% c("Subcanopy large_Subcanopy small", "Subcanopy small_Subcanopy large") ~ "Subcanopy large_Subcanopy small",
        TRUE ~ int))

#right, in the FG + Size model Canopy large_Canopy large is a within-group interaction, but 
#in the Species + Size model this can be both a within- and between- group interaction 
# therefore to understand the DIFFERENCE between the two models, we can look at where the between-
# group interactions overlap: 

#plot
sp_and_fg <- inter |> 
  filter(int %in% c("Canopy large_Subcanopy large", 
                    "Canopy large_Subcanopy small", 
                    "Canopy small_Subcanopy large",
                    "Canopy small_Subcanopy small")) |> 
  mutate(int = case_when(int == "Canopy large_Subcanopy large" ~ "Canopy large ↔ Subcanopy large", 
                         int == "Canopy large_Subcanopy small" ~ "Canopy large  ↔ Subcanopy small", 
                         int == "Canopy small_Subcanopy large" ~  "Canopy small  ↔ Subcanopy large",
                         int == "Canopy small_Subcanopy small" ~ "Canopy small  ↔ Subcanopy small", 
                         TRUE ~ int)) |> 
ggplot(
       aes(x = alpha, 
           y = int, 
           groups = model)) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  geom_violin(colour = "black", 
              scale = "width", 
              width = 0.7,
              fill = NA,
              position = position_dodge(0.75)) +
  geom_point(aes(colour = model), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.15) +
  geom_boxplot(colour = "black", position = position_dodge(0.75),
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  guides(colour = guide_legend(reverse = TRUE, 
                               title = "Model Groupings")) + 
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient") 

############
#also going to plot the between-group interactions where there is no overlap for the two model 
sp <- inter %>% 
  filter(int %in% c("Canopy large_Canopy large", 
                    "Canopy large_Canopy small", 
                    "Canopy small_Canopy small", 
                    "Subcanopy large_Subcanopy large",
                    "Subcanopy large_Subcanopy small",
                    "Subcanopy small_Subcanopy small")) %>% 
  mutate(int = case_when(int == "Canopy large_Canopy large" ~ "Canopy large ↔ Canopy large", 
                         int == "Canopy large_Canopy small" ~ "Canopy large ↔ Canopy small",
                         int == "Canopy small_Canopy small" ~ "Canopy small ↔ Canopy small", 
                         int == "Subcanopy large_Subcanopy large" ~ "Subcanopy large ↔ Subcanopy large", 
                         int == "Subcanopy large_Subcanopy small" ~ "Subcanopy large ↔ Subcanopy small", 
                         int == "Subcanopy small_Subcanopy small" ~ "Subcanopy small ↔ Subcanopy small",  
                         TRUE ~ int)) |> 
  ggplot(
    aes(x = alpha, 
        y = int)) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  geom_violin(colour = "black", 
              scale = "width", 
              width = 0.7,
              fill = NA,
              position = position_dodge(0.75)) +
  geom_boxplot(colour = "black", position = position_dodge(0.75),
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  geom_point(position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.25, 
             colour = "#00BFC4") +
  geom_point(data = preds,
             aes(y = group,
                 x = predicted,
                 group = group, 
                 colour = class_from),
             size = 1.7,
             position = position_jitterdodge(seed = 21,
                                             jitter.width = 0.05,
                                             dodge.width = 0.75)) +
  geom_linerange(data = preds,
                 aes(y = x,
                     x = predicted, 
                     xmin = conf.low,
                     xmax = conf.high),
                 position = position_jitterdodge(seed = 21, 
                                                 jitter.width = 0.05,
                                                 dodge.width = 0.75),
                 linewidth = 1) +
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient") 

 

library(patchwork)
p <- sp_and_fg / sp + plot_annotation(tag_levels = "a") + plot_layout(axis_titles = "collect")
p

##########################################
############### Instead of boxplot and violin, add the predicted lm model instead 
#do some adding of cols and bind 
sp_inter_preds <- read.csv("lm_sp_size_inter.csv")

#and for df_size
fg_inter_preds <- read.csv("lm_fg_size_inter.csv")

#clean up and join together 
fg_inter_preds$model <- "FG + Size"
sp_inter_preds$model <- "Species + Size"

preds_inter <- rbind(fg_inter_preds, sp_inter_preds)

#rename things 
preds_inter <- preds_inter |> 
  mutate(group_name = case_when(x == "Canopy.large_Subcanopy.large" ~ "Canopy large ↔ Subcanopy large", 
                           x == "Canopy.large_Subcanopy.small" ~ "Canopy large ↔ Subcanopy small", 
                           x  == "Canopy.small_Subcanopy.large" ~ "Canopy small ↔ Subcanopy large",
                           x == "Canopy.small_Subcanopy.small" ~ "Canopy small ↔ Subcanopy small",
                           x == "Canopy.large_Canopy.large" ~ "Canopy large ↔ Canopy large", 
                           x  == "Canopy.large_Canopy.small" ~ "Canopy large ↔ Canopy small",
                           x == "Canopy.small_Canopy.small" ~ "Canopy small ↔ Canopy small", 
                           x == "Subcanopy.large_Subcanopy.large" ~ "Subcanopy large ↔ Subcanopy large", 
                           x == "Subcanopy.large_Subcanopy.small" ~ "Subcanopy large ↔ Subcanopy small", 
                           x == "Subcanopy.small_Subcanopy.small" ~ "Subcanopy small ↔ Subcanopy small",  
                           x == "Subcanopy small_Canopy large" ~ "Canopy large ↔ Subcanopy small", 
                           x == "Subcanopy small_Canopy small" ~ "Canopy small ↔ Subcanopy small", 
                           x == "Subcanopy large_Canopy large" ~ "Canopy large ↔ Subcanopy large", 
                           x == "Subcanopy large_Canopy small" ~ "Canopy small ↔ Subcanopy large",
                                 TRUE ~ x)) 

a_pred <- preds_inter |> 
 filter(group_name %in% c("Canopy large ↔ Subcanopy large", 
                         "Canopy large ↔ Subcanopy small", 
                         "Canopy small ↔ Subcanopy large",
                        "Canopy small ↔ Subcanopy small"))

b_pred <- preds_inter |> 
  filter(! group_name %in% c("Canopy large ↔ Subcanopy large", 
                           "Canopy large ↔ Subcanopy small", 
                           "Canopy small ↔ Subcanopy large",
                           "Canopy small ↔ Subcanopy small"))


#similar to the previous plot we need to build two plots and join together 

  
sp_and_fg_inter <- inter |> 
  filter(int %in% c("Canopy large_Subcanopy large", 
                    "Canopy large_Subcanopy small", 
                    "Canopy small_Subcanopy large",
                    "Canopy small_Subcanopy small")) |> 
  mutate(int = case_when(int == "Canopy large_Subcanopy large" ~ "Canopy large ↔ Subcanopy large", 
                         int == "Canopy large_Subcanopy small" ~ "Canopy large ↔ Subcanopy small", 
                         int == "Canopy small_Subcanopy large" ~  "Canopy small ↔ Subcanopy large",
                         int == "Canopy small_Subcanopy small" ~ "Canopy small ↔ Subcanopy small", 
                         TRUE ~ int)) |> 
  mutate(int = as.factor(int)) |> 
  mutate(int = fct_relevel(int,  "Canopy small ↔ Subcanopy small",
                            "Canopy large ↔ Subcanopy small", 
                             "Canopy small ↔ Subcanopy large",
                           "Canopy large ↔ Subcanopy large"
                            )) |> 
  ggplot(
    aes(x = alpha, 
        y = int, 
        groups = model)) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  geom_violin(colour = "gray35", 
              scale = "width", 
              width = 0.55,
              fill = NA,
              linewidth = 0.4,
              position = position_dodge(0.75)) +
  geom_point(aes(colour =model),
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.05) +
  geom_point(data = a_pred,
             aes(y = group_name,
                 x = predicted,
                 groups = model, 
                 colour = model),
             size = 1.7,
             position = position_jitterdodge(seed = 21,
                                             jitter.width = 0.05,
                                             dodge.width = 0.75)) +
  geom_linerange(data = a_pred,
                 aes(y = group_name,
                     x = predicted, 
                     xmin = conf.low,
                     xmax = conf.high, 
                     groups = model, 
                     colour = model),
                
                 position = position_jitterdodge(seed = 21, 
                                                 jitter.width = 0.05,
                                                 dodge.width = 0.75),
                 linewidth = 1) +
  scale_colour_manual(values = c("#A973BA", "black")) + 
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        text = element_text(size = 12.5)) + 
  ylab("") + 
  xlab("Interaction coefficient") 

sp_and_fg_inter



############
#also going to plot the between-group interactions where there is no overlap for the two model 
sp_inter <- inter %>% 
  filter(int %in% c("Canopy large_Canopy large", 
                    "Canopy large_Canopy small", 
                    "Canopy small_Canopy small", 
                    "Subcanopy large_Subcanopy large",
                    "Subcanopy large_Subcanopy small",
                    "Subcanopy small_Subcanopy small")) %>% 
  mutate(int = case_when(int == "Canopy large_Canopy large" ~ "Canopy large ↔ Canopy large", 
                         int == "Canopy large_Canopy small" ~ "Canopy large ↔ Canopy small",
                         int == "Canopy small_Canopy small" ~ "Canopy small ↔ Canopy small", 
                         int == "Subcanopy large_Subcanopy large" ~ "Subcanopy large ↔ Subcanopy large", 
                         int == "Subcanopy large_Subcanopy small" ~ "Subcanopy large ↔ Subcanopy small", 
                         int == "Subcanopy small_Subcanopy small" ~ "Subcanopy small ↔ Subcanopy small",  
                         TRUE ~ int)) |> 
  mutate(int = as.factor(int)) |> 
  mutate(int = fct_relevel(int, "Subcanopy small ↔ Subcanopy small", 
                           "Subcanopy large ↔ Subcanopy small",
                           "Subcanopy large ↔ Subcanopy large",
                           "Canopy small ↔ Canopy small", 
                           "Canopy large ↔ Canopy small", 
                           "Canopy large ↔ Canopy large")) |> 
  ggplot(
    aes(x = alpha, 
        y = int)) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  geom_violin(colour = "gray35", 
              scale = "width", 
              width = 0.55,
              fill = NA,
              linewidth = 0.4,
              position = position_dodge(0.75)) +
  geom_point(position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.075, 
             colour = "black") +
  geom_point(data = b_pred,
             aes(y = group_name,
                 x = predicted),
             size = 1.7,
             position = position_jitterdodge(seed = 21,
                                             jitter.width = 0.05,
                                             dodge.width = 0.75), 
             colour = "black") +
  geom_linerange(data = b_pred,
                 aes(y = group_name,
                     x = predicted, 
                     xmin = conf.low,
                     xmax = conf.high),
                 position = position_jitterdodge(seed = 21, 
                                                 jitter.width = 0.05,
                                                 dodge.width = 0.75),
                 linewidth = 1, 
                 colour = "black") +
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        text = element_text(size = 12.5)) + 
  ylab("") + 
  xlab("Interaction coefficient") 

library(patchwork)

p_inter <- sp_and_fg_inter / sp_inter + 
  plot_annotation(tag_levels = "a") + 
  plot_layout(axis_titles = "collect", guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(0, 20, 1, 0)), 
        legend.text = element_text(margin = margin(r = 13, l = 7))) &
  labs(colour = "Model")
p_inter


ggsave("fig5.jpg", 
       dpi = 300, 
       height = 19, 
       width = 16.5, 
       units = "cm", 
       scale = 1.1)

#spilt plot design 
#additional info in splitting within-group (model-based predictions?) can use? 
# inferences based on summary of glmm - what term has the strongest effect (look at summary)
# or what has 