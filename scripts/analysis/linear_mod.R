library(dplyr)
library(forcats)
library(glmmTMB)
library(effects)
library(ggplot2)
library(ggeffects)
library(emmeans)
library(tidyr)

#load in dataframe of coef estimates v 
df <- read.csv("scripts/model_specifications/species_diameter/df_add3_5_t15.csv")

#or 
df <- read.csv("scripts/model_specifications/species_diameter/df_add3.5_t12.csv")

#or updated sp_size model 
sp_size <- read.csv("scripts/model_specifications/species_diameter/sp_size_df_t15_misc_updated.csv")

#for intraspecific 

intra <- sp_size %>% 
  filter(species_to == species_from)

#rename size_ints  
intra <- intra %>% 
  mutate(size_int = case_when(size_int == "small_small" ~ "Small ↔ Small", 
                               size_int == "small_large" ~ "Small ↔ Large", 
                               size_int == "large_large" ~ "Large ↔ Large")) %>% 
  mutate(size_int = as.factor(size_int)) %>% 
  mutate(size_int = fct_relevel(size_int, "Small ↔ Small", 
                              "Small ↔ Large", 
                              "Large ↔ Large"
                              )) %>% 
  mutate(cc_from = as.factor(class_from)) %>% 
  mutate(cc_from = fct_relevel(class_from, "Subcanopy",  "Canopy",))


##remove misc groups cos it dont mean anything without an fg (this is what we interested in)
intra <- intra %>% 
  filter(! species_from == "Misc")


#look at distribution of response
range(intra$alpha)
hist(intra$alpha)
hist(intra$alpha)



#use package glmtmb
mod1 <- glmmTMB(alpha ~
                1 + size_int + class_from + size_int:class_from + (1|species_from) + (1|site), 
                data = intra)

s <- summary(mod1)


#look at some diagnostics to check assumptions

model_residuals <- residuals(mod1)
hist(model_residuals)

qqnorm(model_residuals)
qqline(model_residuals, col = "red")

plot(mod, 
     type = c("p", "smooth"))



#more diagnositics (perhaps better)
library(DHARMa) 

#tests are more sensitive than the actual model
testDispersion(mod1) #tests if the simulated dispersion is equal to the observed dispersion
r <- simulateResiduals(mod1, n= 1000, plot = TRUE)
#boxplot - 6 groups? ordered? residuals of other groups? 
#some groups are really skewed - non-homogenoueous/outliers 

#horizontal, tidy up 

library(performance)

check_model(mod)

########## plot effects using easystats package universe
library(see)
library(modelbased)
library(parameters)
library(marginaleffects)


#model based estimates of the response variable for 
#different combinations of predictor values 
pred <- estimate_relation(mod1) 
plot(pred)
plot(pred, show_data=T)

#Estimate average values of the response variable at each
#factor level 
estimate_means(mod)


#conditional estimates of the regression
slopes(mod1, variable = "size")
avg_slopes(mod1, variable = "fg") #average marginal estimates

#estimated marginal effects 
estimate_slopes(mod, trend = "fg")
estimate_slopes(mod, trend = "size")

#predictions on y against values of predictions on x
p1 <- plot_predictions(mod)  +
  geom_hline(yintercept = 0, linetype = "dotted")
p1



### visualisation 

ggeffect(mod1, terms = c("size_int", "class_from")) %>% 
  plot(show_data = TRUE, jitter = TRUE)




ggpredict(mod1, terms = c("size_int", "class_from")) %>% 
  plot(show_data = TRUE, 
       jitter = TRUE, 
       dot_size = 3, line_size = 1.5,
       alpha =1,
       dot_alpha = 0.1, 
       n_rows = 1, 
       use_theme = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_manual(
    values = c( "#39568CFF", "#35B779FF"),
    breaks = c("Canopy", "Subcanopy"),
    name = "Functional Group") +
  scale_fill_manual(
    values = c( "#39568CFF", "#35B779FF"),
    breaks = c("Canopy", "Subcanopy"),
    name = "Functional Group") + 
  xlab("") +
  ylab("Alpha coefficient") + 
  coord_flip() +
  ggtitle("") + 
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(), 
        axis.text = element_text(size = 10))

ggsave("linearmod_sp_size_updated.jpg", width = 20, height = 14, units = "cm", dpi = 300)

##take these predictions and put into a df 
pred <- ggpredict(mod1, terms = c("size_int", "class_from"))
intra_pred_df <-as.data.frame(pred)
write.csv(intra_pred_df,"intra_pred_df_updatedmod.csv")
#table in supplement?? visualise 








#for interspecific 

inter <- df %>% 
  filter( ! species_to == species_from)

inter <- inter %>% 
  mutate(size_int = paste0(class_from, "_", class_to)) %>% 
  mutate(size_int = ifelse(size_int == "large_small", "small_large", size_int))
  #%>% 
  # mutate(size_int = case_when(size_int == "small_small" ~ "Small ↔ Small", 
  #                              size_int == "small_large" ~ "Small ↔ Large", 
  #                              size_int == "large_large" ~ "Large ↔ Large"))




inter <- inter %>% 
  rename(group_names = cc_int)

inter %>% count(group)

inter <- inter %>%
  mutate(group = paste0(cc_from, ".", class_from, "_", cc_to, ".", class_to)) %>% 
    mutate(group = case_when(
      group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small",
      group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large",
      group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small",
      group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large",
      group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small",
      group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
      TRUE ~ group))

inter_med <- inter %>% 
  group_by(group) %>% 
  summarise(median = median(alpha))

#Rename for nicer plotting:
inter <- inter %>%
  mutate(group_names = as.factor(group_names)) %>% 
  mutate(size_int = as.factor(size_int)) %>% 
  mutate(group = as.factor(group)) %>% 
  mutate(group = fct_relevel(group, 
                          "Canopy.large_Canopy.large", 
                           "Canopy.large_Canopy.small", 
                           "Canopy.small_Canopy.small", 
                           "Subcanopy.large_Subcanopy.large", 
                           "Subcanopy.large_Subcanopy.small", 
                           "Subcanopy.small_Subcanopy.small", 
                           "Canopy.large_Subcanopy.large", 
                           "Canopy.large_Subcanopy.small", 
                           "Canopy.small_Subcanopy.large",
                           "Canopy.small_Subcanopy.small")) %>% 
  mutate(group = case_when(
    group == "Canopy.large_Canopy.large" ~ "Canopy large ↔ Canopy large", 
    group == "Canopy.large_Canopy.small" ~ "Canopy large ↔ Canopy small", 
    group == "Canopy.small_Canopy.small" ~ "Canopy small ↔ Canopy small",
    group == "Subcanopy.large_Subcanopy.large" ~ "Subcanopy large ↔ Subcanopy large",
    group == "Subcanopy.large_Subcanopy.small" ~ "Subcanopy large ↔ Subcanopy small",
    group == "Subcanopy.small_Subcanopy.small" ~ "Subcanopy small ↔ Subcanopy small",
    group == "Canopy.large_Subcanopy.large" ~ "Canopy large ↔ Subcanopy large",
    group == "Canopy.large_Subcanopy.small" ~ "Canopy large ↔ Subcanopy small",
    group == "Canopy.small_Subcanopy.large" ~ "Canopy small ↔ Subcanopy large",
    group == "Canopy.small_Subcanopy.small" ~ "Canopy small ↔ Subcanopy small",
    TRUE ~ group)) 


mod2 <- glmmTMB(alpha ~
                  1  + size_int + group_names + size_int:group_names + (1|site), 
                data = inter)

s2 <- summary(mod2)


mod3 <- glmmTMB(alpha ~
                  1  + size_int + group_names + size_int:group_names + 
                  (1|species_from) + (`|species_to`)  + (1|site), 
                data = inter)

s3 <- summary(mod3)
### including random effect of species_from and species_to makes everything crown with zero
#makes things look quite not okay 


ggpredict(mod2, terms = c("size_int", "group_names")) %>% 
  plot(show_data = TRUE, 
       jitter = TRUE, 
       dot_size = 3,
       dot_alpha = 0.05,
       n_rows = 1, 
       use_theme = FALSE, 
       line_size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("") +
  ylab("Interaction coefficient") + 
  coord_flip() +
  ggtitle("") + 
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(), 
        axis.text = element_text(size = 12))

ggeffect(mod2, terms = c("group_names", "size_int")) %>% 
  plot(show_data = TRUE, 
       dot_alpha = 0.01) + 
  coord_flip()

##sorting doesn't work, need to grab the prediction tibble and then contrust ggplot 

## Grabbing predictions 

pred <- ggpredict(mod2, terms = c("size_int", "group_names"))
inter_pred_df <-as.data.frame(pred)

inter_pred_df <- inter_pred_df %>%
  mutate(full_group = paste0(x, "_", group)) %>% 
  separate(x, into = c("size_from", "size_to"), sep = "_") %>%
  separate(group, into = c("class_from", "class_to"), sep = "_") %>%
  mutate(
    label = paste(class_from, size_from, "_", class_to, size_to)
  )

#start constructing ggplot 

 ggplot(inter_pred_df) + 
  geom_point(aes(x = predicted, 
                 y = label)) + 
  geom_errorbar(aes(xmin = `conf.low`, 
                    xmax = `conf.high`, 
                    y = label), 
                    width = 0.01) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") + 
  theme_bw() + 
   ylab("") + xlab("Predicted interaction coefficient")
