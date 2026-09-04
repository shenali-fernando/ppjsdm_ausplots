library(dplyr)
library(forcats)
library(glmmTMB)
library(ggplot2)
library(ggeffects)
library(performance)
library(tidyr)
library(DHARMa) 
library(gtsummary)
library(emmeans)

### LINEAR MODELLING 

## We want to make some predictions of what the interaction coefficient may be based on 
# fg and size. We can do this using linear modelling 

#Load alpha coefficients 
#updated sp_size model 
sp_size <- read.csv("../../scripts/model_specifications/species_diameter/sp_size_df_t15_misc_updated.csv")
#and fg_size model 
fg_size <- read.csv("../../scripts/model_specifications/fg_size/fg_size_df_t15_updated.csv")


#for intraspecific species + size:
intra <- sp_size %>% 
  filter(species_to == species_from)

#Some renaming to make life easier 
intra <- intra |> 
  rename(FG = class_from)



 #rename size_ints levels
intra <- intra %>% 
  mutate(size_int = case_when(size_int == "small_small" ~ "Small ↔ Small", 
                              size_int %in% c("small_large") ~ "Small ↔ Large", 
                              size_int   == "large_large" ~ "Large ↔ Large")) %>% 
  mutate(size_int = as.factor(size_int)) %>% 
  mutate(size_int = fct_relevel(size_int, "Small ↔ Small", 
                              "Small ↔ Large", 
                              "Large ↔ Large"
                              )) %>% 
  mutate(FG = as.factor(FG)) %>% 
  mutate(FG = fct_relevel(FG, "Subcanopy",  "Canopy",))


##remove misc groups cos it dont mean anything without an fg (this is what we interested in)
intra <- intra %>% 
  filter(! species_from == "Misc")


#look at distribution of response
range(intra$alpha)
hist(intra$alpha)
hist(intra$alpha)



#use package glmtmb
mod1 <- glmmTMB(alpha ~
                1 + size_int + FG + size_int:FG + (1|species_from) + (1|site), 
                data = intra)

 summary(mod1)

table_sp <- tbl_regression(mod1, 
               intercept = T,  
               estimate_fun = label_style_sigfig(4)) |> 
  modify_header(label = "**Fixed Effect**") 


#look at some diagnostics to check assumptions

model_residuals <- residuals(mod1)
hist(model_residuals)

qqnorm(model_residuals)
qqline(model_residuals, col = "red")

plot(mod, 
     type = c("p", "smooth"))


#more diagnositics (perhaps better)
#tests are more sensitive than the actual model
testDispersion(mod1) #tests if the simulated dispersion is equal to the observed dispersion
r <- simulateResiduals(mod1, n= 1000, plot = TRUE)
#boxplot - 6 groups? ordered? residuals of other groups? 
#some groups are really skewed - non-homogenoueous/outliers 

#using performance package
check_model(mod1)

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

ggeffect(mod1, terms = c("size_int", "FG")) %>% 
  plot(show_data = TRUE, jitter = TRUE)



ggpredict(mod1, terms = c("size_int", "FG")) %>% 
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


#Also need to know what groups are significantly different from each other 
mod1_testpred <- test_predictions(mod1, 
                 terms = c("size_int", "FG"))

predict_response(mod1, 
                 terms = c("size_int", "FG"))

################################################################################

### also want to do a similar thing for the fg + size model

#get within as well and clean
intra_fg <- fg_size %>% 
  filter(class_to == class_from)

#Some renaming to make life easier 
intra_fg <- intra_fg |> 
  rename(FG = class_from)

intra_fg <- intra_fg |> 
  rename(`Size Interaction` = size_int)

#rename size_ints  
intra_fg <- intra_fg %>% 
  mutate(`Size Interaction` = case_when(`Size Interaction` == "small small" ~ "Small ↔ Small", 
                                        `Size Interaction` == "small large" ~ "Small ↔ Large", 
                                        `Size Interaction` == "large large" ~ "Large ↔ Large")) %>% 
  mutate(`Size Interaction`= as.factor(`Size Interaction`)) %>% 
  mutate(`Size Interaction` = fct_relevel(`Size Interaction`, "Small ↔ Small", 
                                "Small ↔ Large", 
                                "Large ↔ Large"
  )) %>% 
  mutate(FG = as.factor(FG)) %>% 
  mutate(FG = fct_relevel(FG, "Subcanopy",  "Canopy"))



#make linear mod
lm_fg <- glmmTMB(alpha ~
                  1 + size_int + FG + size_int:FG + (1|site), 
                data = intra_fg) 

table_fg <-  tbl_regression(lm_fg,
                 intercept = T, 
                 estimate_fun = label_style_sigfig(4)) |> 
  modify_header(label = "**Fixed Effect**") 

#merge the intra sp and intra fg tables together into one and convert to gt
table_intra <- tbl_merge(list(table_fg, table_sp), 
          tab_spanner = c("FG + Size Model", "Species + Size Model")) |> 
  as_gt() 
  
table_intra |> 
  cols_hide(columns = c("groupname_col_1", "groupname_col_2")) |> 
  rm_source_notes() |> 
  tab_style(style = list(cell_text(weight = "bold", style = "italic")), 
            locations = cells_column_spanners(spanners = everything())) |> 
  tab_style(style = list(cell_text(weight = "bold")), 
            locations = cells_body(columns = "label", rows = c(2, 6, 9))) |> 
  tab_style(style = list(cell_text(style = "italic")), 
            locations = cells_body(columns = "label", rows = c(1, 3, 4, 5, 7, 8, 10, 11))) |> 
  opt_horizontal_padding(scale = 1.75) |> 
  opt_vertical_padding(scale = 0.7) |> 
  opt_table_lines(extent = "default") |> 
  tab_style(
    style = cell_borders(sides = "bottom", color = "white"), 
    locations = cells_body(columns = everything(),
                           rows = c(1:11))) 

  
  #gtsave(filename = "table_intra_lm.docx")

glance(lm_fg)
sum <- summary(lm_fg)

compare_performance(lm_fg, mod1)
##check models 
check_model(lm_fg)
model_performance(lm_fg)
check_singularity(lm_fg)

## for within - group models we can compare the performance 
compare_performance(mod1, lm_fg)

### Model diagnostics 
testDispersion(lm_fg) #tests if the simulated dispersion is equal to the observed dispersion
r <- simulateResiduals(lm_fg, n= 1000, plot = TRUE)


#ggpredict and visualise
ggpredict(lm_fg, terms = c("size_int", "FG")) %>% 
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


##take these predictions and put into a df 
pred <- ggpredict(lm_fg, terms = c("size_int", "class_from"))
intra_pred_df <-as.data.frame(pred)
write.csv(intra_pred_df,"lm_fg_size_intra.csv")



#Also need to know what groups are significantly different from each other 
lmfg_testpred <- test_predictions(lm_fg, 
                 terms = c("size_int", "FG"))

predict_response(lm_fg, 
                 terms = c("size_int", "FG"))








#################################################################################################################################


#################################################################################
########################### BETWEEN - GROUP COEFS ###################################
#################################################################################
# For interspecific interactions we have more of a problem 
# We cannot easily get the random effect of species into a glmm model here due to its double sided nature 
# i.e. there will be two species/groups that will need to be added into the glmm in that case 
# this has the effect of bringing to the prediction to basically 0 as with so many random effects included it 
# can be hard to find an effect (see last bit of this script). So let's try a glm instead to least get 
# this aligned with the intraspecific coefficients 


#species + size first 

inter_sp <- sp_size %>% 
  filter( ! species_to == species_from)

inter_sp <- inter_sp %>% 
  mutate(class_int = paste0(class_from, "_", class_to)) %>% 
  mutate(size_int = ifelse(size_int == "large_small", "small_large", size_int))

inter_sp <- inter_sp |> 
  filter(! species_from == "Misc") |> 
  filter(! species_to == "Misc") 


## random of site doesnt make any differece
inter_mod_sp <- glmmTMB(alpha ~
                    1  + size_int + class_int + size_int:class_int + (1 |site), 
                    data = inter_sp)
summary(inter_mod_sp)

## there are repetition of levels here 
inter_mod1_sp <- glm(alpha ~ 
                       1 + size_int + class_int + size_int:class_int, 
                     data = inter_sp)


##### This is the final inter_sp_size model we're going with 
inter_mod_sp <- glm(alpha ~ 
                      1 + group, 
                    data = inter_sp)

mod1 <- glmmTMB(alpha ~ 1 + group + (1 |site), 
                data = inter_sp)

compare_performance(inter_mod_sp, mod1)
## visualise 
ggpredict(mod1, terms = c("group")) %>% 
  plot(show_data = TRUE, 
       jitter = TRUE, 
       dot_size = 3, line_size = 1.5,
       alpha =1,
       dot_alpha = 0.1, 
       n_rows = 1, 
       use_theme = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Alpha coefficient") + 
  coord_flip() +
  ggtitle("") + 
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(), 
        axis.text = element_text(size = 10))

##Its very hard to compare mod1 direct ggpredict figure with above
#so make to df first 
inter_sp_pred <- ggpredict(inter_mod_sp, terms =  c("group"))
inter_sp_df <-as.data.frame(inter_sp_pred)

inter_sp1_df <- inter_sp1_df |> 
  mutate(group2 = paste0(str_extract(group, "[^_]+"), ".", 
                         str_extract(x, "[^_]+"), "_", 
                         str_extract(group, "(?<=_).*"), ".", 
                         str_extract(x, "(?<=_).*"))
  )#clean names 

ggplot(data = inter_sp1_df, 
       aes(x = predicted, 
           y = group2)) + 
  geom_point() + 
  geom_linerange(aes(y = group2,
                     x = predicted, 
                     xmin = conf.low,
                     xmax = conf.high)) + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(), 
        axis.text = element_text(size = 10)) + 
  xlim(c(-1.5, 1.5))


#Output as a df so I can actually visualise it 
inter_sp_pred <- ggpredict(inter_mod_sp, terms =  c("size_int", "class_int"))
inter_sp_df <-as.data.frame(inter_sp_pred)

inter_sp_df <- inter_sp_df |> rename(interacting_groups = x)
  

inter_sp_pred <- ggpredict(inter_mod1_sp, terms =  c("size_int", "class_int"))
inter_sp1_df <-as.data.frame(inter_sp_pred)

write.csv(inter_sp_df,"lm_sp_size_inter.csv")
  
  
  

### For Species + Size, it's probably easier to put everything into a dataframe and compare


compare_mods <- data.frame(interacting_groups = levels(as.factor(inter_sp$group)), 
                           raw_median = as.numeric(NA), 
                           lm_nested = as.numeric(NA), 
                           lm_direct = as.numeric(NA), 
                           lmm_nested = as.numeric(NA), 
                           lmm_direct = as.numeric(NA)
)

a <- inter_sp |> 
  group_by(group) |> 
  summarise(median(alpha)) |> 
  ungroup()

compare_mods$raw_median <- a$`median(alpha)`






###############
#and similar with fg + size 

inter_fg <- fg_size %>% 
  filter( ! class_to == class_from)

inter_fg <- inter_fg %>% 
  mutate(class_int = paste0(class_from, "_", class_to)) %>% 
  mutate(size_int = ifelse(size_int == "large_small", "small_large", size_int))

inter_fg <- inter_fg |> 
  rename(Interaction = int)


inter_mod_fg <- glm(alpha ~
                      1  + Interaction, 
                    data = inter_fg)
summary(inter_mod_fg)


ggpredict(inter_mod_fg, terms = c("group")) %>% 
  plot(show_data = TRUE, 
       jitter = TRUE, 
       dot_size = 3, line_size = 1.5,
       alpha =1,
       dot_alpha = 0.1, 
       n_rows = 1, 
       use_theme = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Alpha coefficient") + 
  coord_flip() +
  ggtitle("") + 
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(), 
        axis.text = element_text(size = 10))

#Output as a df so I can actually visualise it 
inter_fg_pred <- ggpredict(inter_mod_fg, terms = c("group"))
inter_fg_df <-as.data.frame(inter_fg_pred)
write.csv(inter_fg_df,"lm_fg_size_inter.csv")





#Rename for nicer plotting:
inter_sp <- inter_sp %>%
 # rename(Interaction = group) |> 
  mutate(Interaction = case_when(
    Interaction == "Canopy.large_Canopy.large" ~ "Canopy large ↔ Canopy large", 
    Interaction == "Canopy.large_Canopy.small" ~ "Canopy large ↔ Canopy small", 
    Interaction == "Canopy.small_Canopy.small" ~ "Canopy small ↔ Canopy small",
    Interaction == "Subcanopy.large_Subcanopy.large" ~ "Subcanopy large ↔ Subcanopy large",
    Interaction == "Subcanopy.large_Subcanopy.small" ~ "Subcanopy large ↔ Subcanopy small",
    Interaction == "Subcanopy.small_Subcanopy.small" ~ "Subcanopy small ↔ Subcanopy small",
    Interaction == "Canopy.large_Subcanopy.large" ~ "Canopy large ↔ Subcanopy large",
    Interaction == "Canopy.large_Subcanopy.small" ~ "Canopy large ↔ Subcanopy small",
    Interaction == "Canopy.small_Subcanopy.large" ~ "Canopy small ↔ Subcanopy large",
    Interaction == "Canopy.small_Subcanopy.small" ~ "Canopy small ↔ Subcanopy small",
    TRUE ~ Interaction)) |>  
  mutate(Interaction = as.factor(Interaction)) |> 
  mutate(Interaction = fct_relevel(Interaction, 
                                   "Canopy large ↔ Subcanopy large", 
                                   "Canopy large ↔ Subcanopy small", 
                                   "Canopy small ↔ Subcanopy large",
                                   "Canopy small ↔ Subcanopy small",
                                   "Canopy large ↔ Canopy large", 
                                   "Canopy large ↔ Canopy small", 
                                   "Canopy small ↔ Canopy small", 
                                   "Subcanopy large ↔ Subcanopy large", 
                                   "Subcanopy large ↔ Subcanopy small", 
                                   "Subcanopy small ↔ Subcanopy small", 
  ))


inter_fg <- inter_fg |> 
 # rename(Interaction = int) |> 
  mutate(Interaction = case_when(Interaction == "Subcanopy large_Canopy large" ~ "Canopy large ↔ Subcanopy large",
                                 Interaction == "Subcanopy large_Canopy small" ~ "Canopy small ↔ Subcanopy large",
                                 Interaction == "Subcanopy small_Canopy large" ~ "Canopy large ↔ Subanopy small",
                                 Interaction == "Subcanopy small_Canopy small" ~ "Canopy small ↔ Subcanopy small",
                                 Interaction == "Canopy large ↔ Subanopy small" ~ "Canopy large ↔ Subcanopy small",
                                 TRUE ~ Interaction)) |> 
  mutate(Interaction = as.factor(Interaction)) |> 
  mutate(Interaction = fct_relevel(Interaction, 
                                   "Canopy large ↔ Subcanopy large", 
                                   "Canopy large ↔ Subcanopy small", 
                                   "Canopy small ↔ Subcanopy large",
                                   "Canopy small ↔ Subcanopy small"))






 ##################################################################################
 #### Make between - group table 
 
 
 tab_sp <- tbl_regression(inter_mod_sp, 
                            intercept = T,  
                            estimate_fun = label_style_sigfig(4)) |> 
   modify_header(label = "**Effect**") 
   
 
 
 tab_fg <-  tbl_regression(inter_mod_fg,
                             intercept = T, 
                             estimate_fun = label_style_sigfig(4)) |> 
 modify_header(label = "**Effect**") 
 
 #merge the intra sp and intra fg tables together into one and convert to gt
 table_inter <- tbl_merge(list(tab_fg, tab_sp), 
                          tab_spanner = c("FG + Size Model", "Species + Size Model")) |> 
   as_gt() 
 
 table_inter |> 
   rm_source_notes() |> 
   tab_style(style = list(cell_text(weight = "bold", style = "italic")), 
             locations = cells_column_spanners(spanners = everything())) |> 
   tab_style(style = list(cell_text(weight = "bold")), 
             locations = cells_body(columns = "label", rows = c(2))) |> 
   tab_style(style = list(cell_text(style = "italic")), 
             locations = cells_body(columns = "label", rows = c(1, 3:12))) |> 
   opt_horizontal_padding(scale = 1.75) |> 
   opt_vertical_padding(scale = 0.7) |> 
   opt_table_lines(extent = "default") |> 
   tab_style(
     style = cell_borders(sides = "bottom", color = "white"), 
     locations = cells_body(columns = everything(),
                            rows = c(1:12))) 

 
## Also check the performance of the models 
compare_performance(inter_mod_fg, inter_mod_sp)
