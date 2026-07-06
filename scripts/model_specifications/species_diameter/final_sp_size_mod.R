library(ppjsdm)
library(dplyr)
library(stringr)

#load functions
source("scripts/model_specifications/species_diameter/size_funs.R")
source("functions/make_summary_fun.R")

#load data
data <- read.csv("data/ausplots/data_cleaned.csv")

#run function through sites to get the raw model specifications for the species_size model 
#we now have a misc group and also thinking about a threshold = 10

sites <- unique(data$Site_Name)

#take out sites with ci_range > 5
sites <- sites[! sites %in% c("Bird", "BirdTree", "BlackRiver", "Bruxner", 
                              "Caveside", "Flowerdale", "Lardner", "Lorne","Mackenzie",
                              "MtField", "NorthStyx", 
                              "Tinebank", "WaratahMix", "Weld")]

full_df <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                       "hi_numerical", "Potential","site")

set.seed(3)
for (i in sites){
  
  site_mod <- size_sites(site = i, #single site 
                         show_size_freq = FALSE,
                         config_only = FALSE, #if TRUE, returns config only and exits function
                         threshold = 12, 
                         short_range = 10, 
                         short_model = "exponential")
  
  working_df <- make_sum_df(fits = list(site_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(site_mod$sum))
  
  working_df <- working_df %>% mutate(site = i)                              
  
  full_df<- rbind(full_df, working_df)
}



#Check unreasonable cis - 13/48 sites with one problem species each 

full_df <- full_df %>% 
  mutate(range_ci = hi - lo)

a <- full_df %>% filter(range_ci > 5)
a %>% count(site)

threshold <- 12
##############################################################################
## Bird, just Acacia melanoyxlon large 

df <- data %>%
  filter(Site_Name == "Bird")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 == "Acacia melanoxylon large", "Misc_large", species_size2))

#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                              window = window,
                              short_range = matrix(10, nspecies, nspecies), 
                              model = "exponential",
                              saturation = 10, 
                              nthreads = 4, 
                              use_regularization = FALSE, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                           summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Bird")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

###############################################################################
## BirdTree; needs threshold = 15
df <- data %>%
   filter(Site_Name == "BirdTree")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 



d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < 15 & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < 15 & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "BirdTree")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)


############################################################################
## BlackRiver
df <- data %>%
  filter(Site_Name == "BlackRiver")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 %in% c("Acacia mucronata large", "Monotoca glauca large"), "Misc_large", species_size2))

#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "BlackRiver")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)



#################################################################################
##Lorne
df <- data %>%
  filter(Site_Name == "Lorne")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < 13)) 


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Lorne")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)


###################################################################################
## Bruxner, 

df <- data %>%
  filter(Site_Name == "Bruxner")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 == "Geissois benthamii large", "Misc_large", species_size2)) %>% 
  mutate(species_size2 = ifelse(species_size2 == "Niemeyeria whitei small", "Misc_small", species_size2))

#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Bruxner")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

#################################################################################
## Caveside 

df <- data %>%
  filter(Site_Name == "Caveside")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 == "Acacia dealbata small", "Misc_small", species_size2)) %>%
  mutate(species_size2 = ifelse(species_size2 == "Olearia argophylla small", "Misc_small", species_size2)) %>% 
  filter(! species_size2 == "Misc_large")


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                              window = window,
                              short_range = matrix(10, nspecies, nspecies), 
                              model = "exponential",
                              saturation = 10, 
                              nthreads = 4, 
                              use_regularization = FALSE, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                           summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Caveside")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

#################################################################################
## Mackenzie

df <- data %>%
  filter(Site_Name == "Mackenzie")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  filter(! species_size2 == "Misc_small")


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Mackenzie")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

#################################################################################
## Flowerdale

df <- data %>%
  filter(Site_Name == "Flowerdale")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 %in% c("Nematolepis squamea large","Acacia melanoxylon large"),
                                "Misc_large", species_size2))


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Flowerdale")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)


##################################################################################
## Lardner 

df <- data %>%
  filter(Site_Name == "Lardner")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 == "Acacia melanoxylon large", "Misc_large", species_size2)) %>%
  filter(! species_size2 == "Misc_small")


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Lardner")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

###################################################################################
## MtField

df <- data %>%
  filter(Site_Name == "MtField")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 == "Eucalyptus urnigera small", "Misc_small", species_size2))


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "MtField")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

#################################################################################
## NorthStyx 
df <- data %>%
  filter(Site_Name == "NorthStyx")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  mutate(species_size2 = ifelse(species_size2 == "Nothofagus cunninghamii large", "Misc_large", species_size2))


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "NorthStyx")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)


##################################################################################
## Tinebank 

df <- data %>%
  filter(Site_Name == "Tinebank")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  filter(! species_size2 %in% c("Misc_small", "Eucalyptus pilularis large"))


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Tinebank")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

###################################################################################
## WaratahMix

df <- data %>%
  filter(Site_Name == "WaratahMix")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  filter(! species_size2 == "Eucalyptus cypellocarpa large")


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "WaratahMix")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)

####################################################################################
## Weld

df <- data %>%
  filter(Site_Name == "Weld")

df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


d <- df %>% 
  filter(!Genus_Species == "Unidentified tree") %>% 
  group_by(Genus_Species) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size_class = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(species_size = paste0(Genus_Species, " ", size_class)) %>% 
  group_by(species_size) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>%
  mutate(species_size2 = case_when(
    observation_count < threshold & str_ends(species_size, "small") ~ "Misc_small",
    observation_count < threshold & str_ends(species_size, "large") ~ "Misc_large",
    TRUE ~ species_size)) %>% 
  group_by(species_size2) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(!(obs < threshold)) 

d <- d %>% 
  filter(! species_size2 == "Misc_large")


#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
plot(configuration)

d %>% count(species_size2)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = window,
                     short_range = matrix(10, nspecies, nspecies), 
                     model = "exponential",
                     saturation = 10, 
                     nthreads = 4, 
                     use_regularization = FALSE, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)


sum <- summary(fit)
extra_df <- make_sum_df(fits = list(fit), 
                        summ = list(sum))

extra_df <- extra_df %>% mutate(site = "Weld")

extra_df <- extra_df %>% 
  mutate(range_ci = hi - lo)

at <- extra_df %>% filter(range_ci > 5)

full_df <- rbind(full_df, extra_df)


###############################################################################

## checks on full df
full_df %>% count(site)

all(!duplicated(full_df$alpha))

full_df %>% filter(range_ci > 5)



####################
# Make some useful plotting columns
full_df2 <- full_df %>% 
  mutate(class_from = str_extract(from, "\\w+$")) %>% 
  mutate(class_to = str_extract(to, "\\w+$")) %>% 
  mutate(species_from = str_extract(from, "\\w+\\s+\\w+")) %>% 
  mutate(species_to = str_extract(to, "\\w+\\s+\\w+"))

full_df2 <- full_df2 %>% mutate(class_int = paste(class_from, sep = "_", class_to))

#make regions 
full_df2 <- full_df2 %>% 
  mutate(georegion = case_when(site %in% c("Weeaproinah", "Turtons", "Lardner") ~ "S_VIC", 
                               site %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek") ~ "N_VIC", 
                               site %in% c("Dawson", "Frankland", "Clare", "Giants") ~ "S_WA",
                               site %in% c("Carey", "Dombakup", "Warren",  "Sutton","Collins") ~ "N_WA", 
                               site %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton") ~ "QLD", 
                               site %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans") ~ "N_NSW", 
                               site %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo") ~ "S_NSW", 
                               site %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx") ~ "d_TAS", 
                               site %in% c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Flowerdale", "Dip") ~ "o_TAS")) 




#Add sig 
full_df2 <- full_df2 %>% 
  mutate(sig = ifelse(lo > 0 | hi < 0, 1 , NA))


#Get rid of Potential col
full_df2 <- full_df2 %>% 
  select(-Potential)


### save this df at this point 
write.csv(full_df2, "sp_size_df_t12_withmisc.csv")


#Oop need to fix an issue with misc 
full_df2 <- full_df2 %>% 
  mutate(class_from = case_when(class_from == "Misc_small" ~ "small", 
                                class_from == "Misc_large" ~ "large", 
                                TRUE ~ class_from)) %>% 
  mutate(class_to = case_when(class_to == "Misc_small" ~ "small", 
                              class_to == "Misc_large" ~ "large", 
                              TRUE ~ class_to)) 


full_df2 <- full_df2 %>% mutate(class_int = paste(class_from, sep = "_", class_to))


#Assign a species to canopy or subcanopy 
species_class <- read.csv("data/ausplots/species_class.csv")
species_class <- species_class %>% 
  mutate(Genus_Species = iconv(Genus_Species, to = "ASCII", sub = " "))

species_class <- species_class %>% dplyr::rename(species_from = Genus_Species)
full_df2 <- left_join(full_df2, species_class, by = "species_from")
full_df2 <- full_df2 %>% rename(cc_from = Class)

species_class <- species_class %>% dplyr::rename(species_to = species_from)
full_df2 <- left_join(full_df2, species_class, by = "species_to")
full_df2 <- full_df2 %>% rename(cc_to = Class)

#For misc put misc 

full_df2 <- full_df2 %>% 
  mutate(cc_to = ifelse(is.na(cc_to), "Misc", cc_to)) 

full_df2 <- full_df2 %>% 
  mutate(cc_from = ifelse(is.na(cc_from), "Misc", cc_from)) 


full_df2 <- full_df2 %>% 
  mutate(species_to = ifelse(is.na(species_to), "Misc", species_to)) 

full_df2 <- full_df2 %>% 
  mutate(species_from = ifelse(is.na(species_from), "Misc", species_from)) 


full_df2 <- full_df2 %>% 
  mutate(cc_int = paste0(cc_from, "_", cc_to))

#check classes 
full_df2 %>% count(cc_int)

#Put the same classes together because there is no difference between S_C and C_S 
full_df2 <- full_df2 %>% 
  mutate(cc_int = case_when(cc_int == "Canopy_Misc" ~ "Misc_Canopy", 
                            cc_int == "Canopy_Subcanopy" ~ "Subcanopy_Canopy", 
                            cc_int == "Subcanopy_Misc" ~ "Misc_Subcanopy", 
                            TRUE ~ cc_int))

##Add a functional group size column 
full_df2 <- full_df2 %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to))

full_df2 %>% count(group) 
  
full_df2 <- full_df2 %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
    group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
    group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
    group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
    group %in% c("Misc.large_Canopy.large", "Canopy.large_Misc.large") ~ "Misc.large_Canopy.large", 
    group == "Misc.large_Misc.small" ~ "Misc.small_Misc.large",
    group  %in% c("Misc.large_Canopy.small", "Canopy.small_Misc.large") ~ "Misc.large_Canopy.small",
    group %in% c("Misc.large_Subcanopy.small", "Subcanopy.small_Misc.large") ~ "Misc.large_Subcanopy.small", 
    group %in% c("Misc.large_Subcanopy.large", "Subcanopy.large_Misc.large") ~ "Misc.large_Subcanopy.large",
    group  %in% c("Misc.small_Canopy.small", "Canopy.small_Misc.small") ~ "Misc.small_Canopy.small",
    group %in% c("Misc.small_Subcanopy.small", "Subcanopy.small_Misc.small") ~ "Misc.small_Subcanopy.small", 
    group %in% c("Misc.small_Subcanopy.large", "Subcanopy.large_Misc.small") ~ "Misc.small_Subcanopy.large",
    TRUE ~ group))


#Fix two issues 

full_df2 <- full_df2 %>% 
  mutate(cc_to = ifelse(species_to == "Synoum glandulosum", "Subcanopy", cc_to)) %>% 
  mutate(cc_from = ifelse(species_from == "Synoum glandulosum", "Subcanopy", cc_from)) %>% 
  mutate(cc_to = ifelse(species_to == "Lophostemon sp", "Subcanopy", cc_to)) %>% 
  mutate(cc_from = ifelse(species_from == "Lophostemon sp", "Subcanopy", cc_from))

#write out final version 

write.csv(full_df2, "sp_size_df_t12_withmisc_clean.csv")
