library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")
source("code/make_summary_fun.R")

####### Dominant Eucalypt == E. regnans 
reg <- data_c %>% filter(Genus_Species == "Eucalyptus regnans") %>% 
  group_by(Site_Name) %>% 
  count()
reg

### ANU101 
anu101 <- data_c %>% filter(Site_Name == "ANU101")
anu101 %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

anu101_c <- anu101 %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class))  %>% 
  filter(! new_cc == "Non-euc Intermediate")

anu101_c %>% count(new_cc)

configuration_anu101 <- Configuration(anu101_c$Ausplot_X, anu101_c$Ausplot_Y, types = anu101_c$new_cc)
plot(configuration_anu101)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
nspecies <- length(levels(configuration_anu101$types))
short_range <- matrix(8, nspecies, nspecies)

fit_anu101_2 <- ppjsdm::gibbsm(configuration = configuration_anu101, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_anu101_2 <- summary(fit_anu101_2)


ppjsdm::box_plot(fit_anu101_2,fit_anu101_2,
                 summ = list(sum_anu101, sum_anu101_2),
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)

ppjsdm::box_plot(fit_anu101_2,
                 summ = sum_anu101_2,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### ANU363 
anu363 <- data_c %>% filter(Site_Name == "ANU363")
anu363 %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

anu363_c <- anu363 %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class))

anu363_c %>% count(new_cc)


configuration_anu363 <- Configuration(anu363_c$Ausplot_X, anu363_c$Ausplot_Y, types = anu363_c$new_cc)
#configuration_anu363 <- Configuration(anu363$Ausplot_X, anu363$Ausplot_Y, types = anu363$Crown_Class)
plot(configuration_anu363)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

# Set parameters
nspecies <- length(levels(configuration_anu363$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_anu363 <- ppjsdm::gibbsm(configuration = configuration_anu363, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_anu363 <- summary(fit_anu363)


ppjsdm::box_plot(fit = fit_anu363,
                 summ = sum_anu363,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)

ppjsdm::box_plot(fit_anu363, fit_anu363_2,
                 summ = list(sum_anu363, sum_anu363_2),
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)



### ANU589 
anu589 <- data_c %>% filter(Site_Name == "ANU589")
anu589 %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

anu589_c <- anu589 %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class))  %>% 
  filter(! new_cc == "Non-euc Intermediate")

anu589_c %>% count(new_cc)


configuration_anu589 <- Configuration(anu589_c$Ausplot_X, anu589_c$Ausplot_Y, types = anu589_c$new_cc)
plot(configuration_anu589)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

# Set parameters
nspecies <- length(levels(configuration_anu589$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_anu589 <- ppjsdm::gibbsm(configuration = configuration_anu589, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_anu589 <- summary(fit_anu589)


ppjsdm::box_plot(fit = fit_anu589,
                 summ = sum_anu589,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### Ada Tree 
ada <- data_c %>% filter(Site_Name == "Ada Tree")
ada %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

ada_c <- ada %>%   
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class))  %>% 
  filter(! new_cc == "Non-euc Intermediate")

ada_c %>% count(new_cc)


configuration_ada <- Configuration(ada_c$Ausplot_X, ada_c$Ausplot_Y, types = ada_c$new_cc)
plot(configuration_ada)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

ada_c <- ada_c %>%
  group_by(Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) %>%
  ungroup() 



configuration_ada <- Configuration(ada_c$x_jitter, ada_c$y_jitter, types = ada_c$new_cc)
plot(configuration_ada)

# Set parameters
nspecies <- length(levels(configuration_ada$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_ada <- ppjsdm::gibbsm(configuration = configuration_ada, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_ada <- summary(fit_ada)

ppjsdm::box_plot(fit = fit_ada,
                 summ = sum_ada,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)




### Weeaproinah 
wee <- data_c %>% filter(Site_Name == "Weeaproinah")
wee %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

wee_c <- wee %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class)) %>% 
  filter(! new_cc == "Eucalyptus obliqua Suppressed")

wee_c %>% count(new_cc)

configuration_wee <- Configuration(wee_c$Ausplot_X, wee_c$Ausplot_Y, types = wee_c$new_cc)
plot(configuration_wee)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

# Set parameters
nspecies <- length(levels(configuration_wee$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_wee <- ppjsdm::gibbsm(configuration = configuration_wee, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_wee <- summary(fit_wee)


ppjsdm::box_plot(fit = fit_wee,
                 summ = sum_wee,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)



### Turtons 
turtons <- data_c %>% filter(Site_Name == "Turtons")
turtons %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

turtons_c <- turtons %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class)) 

turtons_c %>% count(new_cc)

configuration_turtons <- Configuration(turtons_c$Ausplot_X, turtons_c$Ausplot_Y, types = turtons_c$new_cc)
plot(configuration_turtons)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


# Set parameters
nspecies <- length(levels(configuration_turtons$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_turtons <- ppjsdm::gibbsm(configuration = configuration_turtons, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_turtons <- summary(fit_turtons)


ppjsdm::box_plot(fit = fit_turtons,
                 summ = sum_turtons,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)





### Lardner 
lardner <- data_c %>% filter(Site_Name == "Lardner")
lardner2 <- lardner %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

lardner_c <- lardner %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class)) 

lardner_c %>% count(new_cc)

configuration_lardner <- Configuration(lardner_c$Ausplot_X, lardner_c$Ausplot_Y, types = lardner_c$new_cc)
plot(configuration_lardner)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


# Set parameters
nspecies <- length(levels(configuration_lardner$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_lardner <- ppjsdm::gibbsm(configuration = configuration_lardner, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_lardner <- summary(fit_lardner)


ppjsdm::box_plot(fit = fit_lardner,
                 summ = sum_lardner,
                 coefficient = "alpha",
                 which = "all",
                 text_size = 12)



### HardyCreek 
hardy <- data_c %>% filter(Site_Name == "HardyCreek")
hardy %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

hardy_c <- hardy %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class)) %>% 
  filter(! new_cc == "Eucalyptus nitens Co-dominant")

hardy_c %>% count(new_cc)

configuration_hardy <- Configuration(hardy_c$Ausplot_X, hardy_c$Ausplot_Y, types = hardy_c$new_cc)
plot(configuration_hardy)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


# Set parameters
nspecies <- length(levels(configuration_hardy$types))
short_range <- matrix(5, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_hardy <- ppjsdm::gibbsm(configuration = configuration_hardy, 
                            window = window,
                            short_range = short_range, 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)
sum_hardy <- summary(fit_hardy)


ppjsdm::box_plot(fit = fit_hardy,
                 summ = sum_hardy,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)




### NorthStyx
norths <- data_c %>% filter(Site_Name == "NorthStyx")
norths %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

norths_c <- norths %>%
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class)) %>% 
  mutate(new_cc = if_else(new_cc == "Eucalyptus regnans Emergent", "Eucalyptus regnans Dominant", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus regnans Dominant", "Eucalyptus delegatensis Dominant", "Eucalyptus delegatensis Emergent"), "Eucalyptus Dominant", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus regnans Co-dominant", "Eucalyptus delegatensis Co-dominant"), "Eucalyptus Co-dominant", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus regnans Suppressed", "Eucalyptus delegatensis Suppressed"), "Eucalyptus Suppressed", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus regnans Intermediate", "Eucalyptus delegatensis Intermediate"), "Eucalyptus Intermediate", new_cc)) %>% 
  filter(! new_cc %in% c("Non-euc Intermediate", "Non-euc Co-dominant"))

norths_c %>% count(new_cc)

configuration_norths <- Configuration(norths_c$Ausplot_X, norths_c$Ausplot_Y, types = norths_c$new_cc)
plot(configuration_norths)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))



# Set parameters
nspecies <- length(levels(configuration_norths$types))
short_range <- matrix(7, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_norths <- ppjsdm::gibbsm(configuration = configuration_norths, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_norths <- summary(fit_norths)


ppjsdm::box_plot(fit = fit_norths,
                 summ = sum_norths,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### Weld
weld <- data_c %>% filter(Site_Name == "Weld")
weld2 <- weld %>% group_by(Genus_Species, Crown_Class, Diameter) %>% 
  count()

weld_c <- weld %>% 
  mutate(species = if_else(str_starts(Genus_Species, "Eucalyptus", negate = TRUE), "Non-euc", Genus_Species)) %>% 
  mutate(new_cc = paste(species, Crown_Class)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus obliqua Co-dominant","Eucalyptus regnans Co-dominant"), "Eucalyptus Co-dominant", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus obliqua Dominant","Eucalyptus regnans Dominant"), "Eucalyptus Dominant", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus obliqua Suppressed","Eucalyptus regnans Suppressed"), "Eucalyptus Suppressed", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus regnans Intermediate", "Eucalyptus obliqua Intermediate"), "Eucalyptus Intermediate", new_cc))

weld_c %>% count(new_cc)

configuration_weld <- Configuration(weld_c$Ausplot_X, weld_c$Ausplot_Y, types = weld_c$new_cc)
plot(configuration_weld)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

# Set parameters
nspecies <- length(levels(configuration_weld$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_weld <- ppjsdm::gibbsm(configuration = configuration_weld, 
                           window = window,
                           short_range = short_range, 
                           model = "exponential", 
                           saturation = 10, 
                           nthreads = 4, 
                           fitting_package = "glmnet",
                           dummy_distribution = "stratified",
                           min_dummy = 1, dummy_factor = 1e10, 
                           max_dummy = 1e3)
sum_weld <- summary(fit_weld)


ppjsdm::box_plot(fit = fit_weld,
                 summ = sum_weld,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)








##########################################################################

fits <- list(
  "ANU101" = fit_anu101_2, 
  "ANU363" = fit_anu363, 
  "ANU689" = fit_anu589, 
  "Ada Tree" = fit_ada, 
  "HardyCreek" = fit_hardy,
  "Weeaproinah" = fit_wee,
  "Lardner" = fit_lardner, 
  "Turtons" = fit_turtons,
  "NorthStyx" = fit_norths,
  "Weld" = fit_weld
)

sums <- list(sum_anu101_2, 
             sum_anu363, 
             sum_anu589, 
             sum_ada,
             sum_hardy,
             sum_wee,
             sum_lardner, 
             sum_turtons,
             sum_norths,
             sum_weld)


df_reg <- make_sum_df(fits = fits, 
                      summ = sums)

df_reg <- df_reg %>% select(-Potential)

write.csv(df_reg, "code/crownclass/regnans_cc.csv")
