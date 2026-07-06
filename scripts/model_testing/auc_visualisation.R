library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(forcats)


########### Visualisation of AUC 

## Read the auc csv made from AUC_all_models.R (this was when threshold was set to 10)
size <- read.csv("scripts/model_testing/auc_mods/auc_size_model.csv")
sp <- read.csv("scripts/model_testing/auc_mods/auc_species_model.csv")
sp_size <- read.csv("scripts/model_testing/auc_mods/auc_species_size_model.csv")
fg <- read.csv("scripts/model_testing/auc_mods/auc_fg_model.csv")
fg_size <- read.csv("scripts/model_testing/auc_mods/auc_fg_size_model.csv")

#increased the threshold to 15 
size <- read.csv("scripts/model_testing/auc_mods/auc_size_model_t15.csv")
sp <- read.csv("scripts/model_testing/auc_mods/auc_species_model_t15.csv")
sp_size <- read.csv("scripts/model_testing/auc_mods/auc_species_size_model_t15.csv")
fg <- read.csv("scripts/model_testing/auc_mods/auc_fg_model_t15.csv")
fg_size <- read.csv("scripts/model_testing/auc_mods/auc_fg_size_model_t15.csv")


### put into a single df 
fg$model <- "FG"
fg_size$model <- "FG + Size"
size$model <- "Size"
sp$model <- "Species"
sp_size$model <- "Species + Size"

auc_full <- rbind(fg, 
                  fg_size,
                  size, 
                  sp, 
                  sp_size)

#get rid of anything with auc = 1 because these are sites with only one group 
auc_full <- auc_full %>% 
  filter(! auc == 1)

#Add some more columns 
auc_full2 <- auc_full %>% 
  filter(! str_detect(type, "Unid tree")) %>% 
  mutate(size = word(type, -1)) %>% 
  mutate(size = if_else(! size %in% c('large', "small"), " ", size)) %>% 
  mutate(fg = word(type, 1)) %>% 
  mutate(fg = if_else(! fg %in% c('Subcanopy', 'Canopy'), " ", fg)) %>% 
  mutate(sp = word(type, 1, 2)) %>% 
  mutate(sp = if_else(str_detect(sp, "small|large"), " ", sp)) %>% 
  mutate(sp = if_else(type == "Misc", "Misc", sp))

#relevel for easier plotting 
auc_full2 <- auc_full2 %>% 
  mutate(model = as.factor(model)) %>% 
  mutate(model = fct_relevel(model, 
                             "Size", 
                             "FG", 
                             
                             "FG + Size", 
                             "Species", 
                             "Species + Size"))

### Then make a simple visualisation: 
ggplot(data = auc_full2,
       aes(x = auc, 
           y = model)) + 
  geom_violin(scale = "width", 
              width = 0.6, 
              linewidth = 0.35) + 
  geom_boxplot(width = 0.35, 
               linewidth = 0.3) + 
  xlab("AUC Score") + 
  ylab("") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.title.x = element_text(margin = margin(t = 12)),
        panel.grid.minor.x = element_line(colour = "gray85", linewidth = 0.15)) + 
  scale_x_continuous(breaks = seq(0.2, 0.9, by= 0.1), 
                     minor_breaks = seq(0.2, 0.9, by= 0.05)) 
  
ggsave("auc_all_mods_t15.png", 
       width = 15.5, 
       height = 11.5, 
       units = "cm", 
       dpi = 300)

#similar pattern to aic, good

### Would be also good to visualise the average plot level difference in AUC score
#because there is an increase between species + size and species models, just the variability hides this 

auc_full2 %>% 
  group_by(site, model) %>% 
  summarise(mean = mean(auc)) %>% 
ggplot(aes(x = model, 
           y = mean)) + 
  geom_point(alpha = 0.4, 
             size = 0.7) +
geom_line(aes(group = site), 
              alpha = 0.45, 
              linewidth = 0.25) +
  ylab("Mean AUC score for a site") + 
  xlab("") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        panel.grid.minor.y = element_line(colour = "gray95", linewidth = 0.1), 
        axis.title.y = element_text(margin = margin(r = 12))) + 
  scale_y_continuous(breaks = seq(0.2, 0.9, by= 0.1), 
                     minor_breaks = seq(0.2, 0.9, by= 0.05)) 

ggsave("auc_site_t15.png", 
       width = 15.5, 
       height = 11.5, 
       units = "cm", 
       dpi = 300)

#Check medians at site and overall
summary_auc <- 
  auc_full2 %>% 
  group_by(site, model) %>% 
  summarise(median = median(auc), 
            mean = mean(auc)) %>% 
  ungroup()


#get overall model level summary 
summary_overall_auc <-
  auc_full2 %>% 
  group_by( model) %>% 
  summarise(median = median(auc), 
            mean = mean(auc)) %>% 
  ungroup()
