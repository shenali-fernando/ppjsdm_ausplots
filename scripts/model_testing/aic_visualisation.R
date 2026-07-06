library(dplyr)
library(ggplot2)
library(forcats)

###### Load in aic dataframes and rbind into one 
null <- read.csv("scripts/model_testing/aic_mods/aic_null.csv")
size <- read.csv("scripts/model_testing/aic_mods/aic_size.csv")
sp <- read.csv("scripts/model_testing/aic_mods/aic_species.csv")
sp_size <- read.csv("scripts/model_testing/aic_mods/aic_sp_size.csv")
fg <- read.csv("scripts/model_testing/aic_mods/aic_fg.csv")
fg_size <- read.csv("scripts/model_testing/aic_mods/aic_fg_size.csv")


## or load in t15 dataframes 
null <- read.csv("scripts/model_testing/aic_mods/aic_null_t15.csv")
size <- read.csv("scripts/model_testing/aic_mods/aic_size_t15.csv")
sp <- read.csv("scripts/model_testing/aic_mods/aic_species_t15.csv")
sp_size <- read.csv("scripts/model_testing/aic_mods/aic_sp_size_t15.csv")
fg <- read.csv("scripts/model_testing/aic_mods/aic_fg_t15.csv")
fg_size <- read.csv("scripts/model_testing/aic_mods/aic_fg_size_t15.csv")


aic_full <- rbind(null, 
                  fg, 
                  fg_size,
                  size, 
                  sp, 
                  sp_size)

## Visualisation 
aic_full %>% 
  mutate(model = as.factor(model)) %>% 
  mutate(model = case_when(model == "null" ~ "Null", 
                           model == "fg" ~ "FG", 
                           model == "size" ~ "Size", 
                           model == "species" ~ "Species", 
                           model == "fg_size" ~ "FG + \n Size", 
                           model == "sp_size" ~ "Species + \n Size")) %>% 
  mutate(model = fct_relevel(model, "Null", "Size" ,"FG", "FG + \n Size",
                              "Species",
                              "Species + \n Size")) %>% 
  ggplot(
       aes(x = model, 
           y = aic)) + 
  geom_point(size = 0.75) + 
  geom_line(aes(group = site), 
            alpha = 0.5, 
            linewidth = 0.3) + 
  xlab("AIC") +
  ylab("") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))

ggsave("AIC_all_mods_t15_2.png", 
       dpi = 300, 
       width = 19.25, 
       height = 21.5, 
       units = "cm")


#check medians quickly 

aic_full %>% group_by(model) %>% 
summarise(median = median)

#AIC penality factor - fgxsize, halved parameters
#AIC changes due to reduction of parameters
#variability of species 