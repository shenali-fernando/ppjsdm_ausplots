library(forcats)
library(ggbeeswarm)
library(ggplot2)
library(dplyr)
library(ggpubr)


#Load in dfs 
sp_size <- read.csv("scripts/model_specifications/species_diameter/df_add3_5_t15.csv")

fg_size <- read.csv("scripts/model_specifications/fg_size/fg_size_df_t15_updated.csv")


#######################################################################################
############## SCATTERPLOT 
## The easiest visual comparison of the two model is a scatterplot of a-priori (fg_size)
#  v the post-hoc (linear model of the species_size) models. There will be the 6 classes 
# for each, and we can see where these line up on the one-to-one line 

#Read in the linear model post-hoc intraspecific coefficients 
intra_sp_preds <- read.csv("scripts/analysis/intra_pred_df.csv")

#Get the medians for the within-sclass interactions for the a-priori (fg_size) model
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
  ylim(c(-0.6, 0.3))


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

#bind 
full_fg_sp <- bind_rows(fg_size, sp_size)

#Need to neaten up some columns for plotting 
full_fg_sp <- full_fg_sp %>% 
  mutate(size_int = case_when(
    size_int %in% c("small small", "small_small") ~ "Small ↔ Small", 
    size_int %in% c("small large", "large_small") ~ "Small ↔ Large", 
    size_int %in% c("large large", "large_large") ~ "Large ↔ Large", 
    TRUE ~ size_int))

full_fg_sp <- full_fg_sp %>% 
  mutate(size_int = as.factor(size_int)) %>% 
  mutate(size_int = fct_relevel(size_int, 
                                "Small ↔ Small",
                                "Small ↔ Large", 
                                "Large ↔ Large")) %>% 
  mutate(model = case_when(model == "fg_size" ~ "FG + Size", 
                           model == "sp_size" ~ "Species + Size")) %>% 
  mutate(model = as.factor(model))


###########INTRA
#filter for within-group 
intra <- full_fg_sp %>% 
  filter(
    (model == "Species + Size" & (species_from == species_to)) | 
    (model == "FG + Size" & (class_from == class_to))) 

#and plot 
ggplot(data = intra,
       aes(x = alpha, 
           y = size_int, 
           groups = model)) + 
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
  geom_point(aes(colour = model), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.25) +
  guides(colour = guide_legend(reverse = TRUE, 
                               title = "Model Groupings")) + 
  facet_grid(~class_from) + 
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient")


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
inter %>% 
  filter(int %in% c("Canopy large_Subcanopy large", 
                    "Canopy large_Subcanopy small", 
                    "Canopy small_Subcanopy large",
                    "Canopy small_Subcanopy small")) %>% 
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
  geom_boxplot(colour = "black", position = position_dodge(0.75),
               outlier.shape = NA, 
               fill = NA,
               width = 0.4) +
  geom_point(aes(colour = model), 
             position = position_jitterdodge(seed = 21, 
                                             jitter.width = 0.3,
                                             dodge.width = 0.75), 
             shape = 19, 
             size = 1.75, 
             alpha = 0.25) +
  guides(colour = guide_legend(reverse = TRUE, 
                               title = "Model Groupings")) + 
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient") + 
  ggtitle("Between-group interactions with overlap of the models")

############
#also going to plot the between-group interactions where there is no overlap for the two model 
inter %>% 
  filter(int %in% c("Canopy large_Canopy large", 
                    "Canopy large_Canopy small", 
                    "Canopy small_Canopy small", 
                    "Subcanopy large_Subcanopy large",
                    "Subcanopy large_Subcanopy small",
                    "Subcanopy small_Subcanopy small")) %>% 
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
             alpha = 0.25) +
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient") + 
  ggtitle("Between-group interactions for Species + Size model only")


