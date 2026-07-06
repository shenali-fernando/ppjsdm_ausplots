library(dplyr)
library(tidyr)
library(ggplot2)
library(spatstat)
library(ppjsdm)
library(ggpubr)
library(patchwork)

##### Shockingly, the AIC and BIC are computed in the fit for the fit!!!

# Cool, i need to understand how well each of the 6 model specifications are doing compared to each other:
#Null, Sp, FG, Size, FG x Size, Sp x Size

##################################################################################
########################## NULL MODEL ###################################

#make empty df 
aic_null <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_null) <- c("site", "aic", "bic", "model")

#Load the data 
data <- read.csv("data/ausplots/data_cleaned.csv")
fg_data <- read.csv("data/ausplots/species_class.csv")

#get the fg dataframe ready
fg_data <- fg_data %>% 
  mutate(Genus_Species = iconv(Genus_Species, to = "ASCII", sub = " "))

data1 <- left_join(data, fg_data, 
                   by = "Genus_Species")

data1 <- data1 %>% 
  mutate(Class = if_else(Genus_Species == "Unidentified tree", "Unid tree",Class)) %>% 
  mutate(Class = if_else(is.na(Class), "Subcanopy", Class)) #everything not in the df is subcanopy

data1 %>% count(Class)

#set sites
site_names <- unique(data$Site_Name)


#Null Model 

#get rid of duplicate points first 
data <- data %>%
  group_by(Site_Name, Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) %>%
  ungroup() 

window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

for(i in seq_along(site_names)){
  
  site <- site_names[i]
  
  df <- data %>% filter(Site_Name == site) #filter data to site
  
  configuration <- Configuration(df$x_jitter, df$y_jitter)
  
  fit <- ppjsdm::gibbsm(configuration = configuration, #fit model
                        window = window,
                        short_range = matrix(8, 1, 1), 
                        model = "exponential", 
                        saturation = 10, 
                        nthreads = 4, 
                        fitting_package = "glmnet",
                        dummy_distribution = "stratified",
                        min_dummy = 1, dummy_factor = 1e10, 
                        max_dummy = 1e3)
  
  aic_null$site[i] <- site
  
  aic_null$aic[i] <- fit$aic 
  
  aic_null$bic[i] <- fit$bic
}

aic_null$model <- "null"
write.csv(aic_null, "aic_null_t15.csv")


################################################################################
############################## SIZE MODEL #######################################

aic_size <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_size) <- c("site", "aic", "bic", "model")


for (i in seq_along(site_names)) {

  site <- site_names[i]
  
  d <- data %>% filter(Site_Name == site)
  
  d <- d %>% 
    filter(!is.na(Diameter)) %>% 
    mutate(size = if_else(Diameter < (mean(Diameter) + 3.5), "small", "large")) %>% 
    group_by(size) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!observation_count < 10) 
  
  
  d <- d %>%
    group_by(Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
    mutate(
      is_duplicated = n() > 1, #create column of TRUE/FALSE 
      #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
      x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
      y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
    ) %>%
    ungroup() 
  
  
  configuration <- Configuration(d$x_jitter, d$y_jitter, types = d$size) #make config 
  
  # Set parameters
  window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
  nspecies <- length(levels(configuration$types))
  short_range <- matrix(10, nspecies, nspecies)
  model <- "exponential"
  
  fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                       window = window,
                       short_range = short_range, 
                       model = "exponential", 
                       saturation = 10, 
                       nthreads = 4, 
                       fitting_package = "glmnet",
                       dummy_distribution = "stratified",
                       min_dummy = 1, dummy_factor = 1e10, 
                       max_dummy = 1e3)

aic_size$site[i] <- site

aic_size$aic[i] <- fit$aic 

aic_size$bic[i] <- fit$bic
}

aic_size$model <- "size"
write.csv(aic_size, "aic_size_t15.csv")

####################################################################################
################################ SPECIES MODEL ####################################

aic_sp <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_sp) <- c("site", "aic", "bic", "model")


for (i in seq_along(site_names)) {
  
  site <- site_names[i]
  
  d <- data %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    group_by(Genus_Species) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    mutate(group = if_else(observation_count < 15, "Misc", Genus_Species)) %>% 
    group_by(group) %>% 
    mutate(obs = n()) %>% 
    ungroup() %>% 
    filter(! obs < 15)
  
  
  d <- d %>%
    group_by(Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
    mutate(
      is_duplicated = n() > 1, #create column of TRUE/FALSE 
      #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
      x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
      y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
    ) %>%
    ungroup() 
  
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$group) #make config 
  
  # Set parameters
  window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
  nspecies <- length(levels(configuration$types))
  short_range <- matrix(10, nspecies, nspecies)
  model <- "exponential"
  
  fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                       window = window,
                       short_range = short_range, 
                       model = "exponential", 
                       saturation = 10, 
                       nthreads = 4, 
                       fitting_package = "glmnet",
                       dummy_distribution = "stratified",
                       min_dummy = 1, dummy_factor = 1e10, 
                       max_dummy = 1e3)
  

  aic_sp$site[i] <- site
  
  aic_sp$aic[i] <- fit$aic 
  
  aic_sp$bic[i] <- fit$bic
  }

aic_sp$model <- "species"
write.csv(aic_sp, "aic_species_t15.csv")


#################################################################################
############################### FG MODEL ########################################

aic_fg <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_fg) <- c("site", "aic", "bic", "model")

for(i in seq_along(site_names)){
  
  site <- site_names[i]
  
  
  df <- data1 %>%
    filter(Site_Name == site)
  
  
  #jitter coordinates
  df <- df %>%
    mutate(
      is_duplicated = n() > 1, #create column of TRUE/FALSE 
      #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
      x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
      y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
    ) 
  
  
  #make config
  configuration <- ppjsdm::Configuration(df$x_jitter,
                                         df$y_jitter, 
                                         types = df$Class)
  
  #set parameters 
  window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                     y_range = c(0, 100))
  
  nspecies <- length(levels(configuration$types))
  
  set.seed(1030)  
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
  
  
  aic_fg$site[i] <- site
  
  aic_fg$aic[i] <- fit$aic 
  
  aic_fg$bic[i] <- fit$bic
}

aic_fg$model <- "fg"
write.csv(aic_fg, "aic_fg_t15.csv")


###################################################################################
########################## SPECIES X SIZE MODEL ###################################

aic_sp_size <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_sp_size) <- c("site", "aic", "bic", "model")


for (i in seq_along(site_names)) {
  
  site <- site_names[i]
  
  df <- data %>%
    filter(Site_Name == site)
  
  
  #jitter coordinates
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
    mutate(new_group = paste0(Genus_Species, sep = " ", size_class)) %>% 
    group_by(new_group) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>%
    mutate(size_class2 = if_else(observation_count < 15, 
                                 "Misc", 
                                 new_group)) %>% 
    group_by(size_class2) %>% 
    mutate(obs = n()) %>% 
    ungroup() %>% 
    filter(!(obs < 15)) 
  
  
  
  #make config
  configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$size_class2)
  
  
  #set parameters 
  window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                     y_range = c(0, 100))
  
  nspecies <- length(levels(configuration$types))
  
  set.seed(1030)  
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
  aic_sp_size$site[i] <- site
  
  aic_sp_size$aic[i] <- fit$aic 
  
  aic_sp_size$bic[i] <- fit$bic
}

aic_sp_size$model <- "sp_size"
write.csv(aic_sp_size, "aic_sp_size_t15.csv")

##################################################################################
################################### FG X SIZE MODEL #############################

aic_fg_size <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_fg_size) <- c("site", "aic", "bic", "model")


for(i in seq_along(site_names)){
  
  site <- site_names[i]

df <- data1 %>%
  filter(Site_Name == site)


#jitter coordinates
df <- df %>%
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) 


#make groups of fg x size  
d <- df %>%
  group_by(Class) %>% 
  mutate(median_diameter = ceiling(median(Diameter, na.rm = TRUE)) + 3.5)  %>% 
  mutate(size = case_when(
    Diameter < median_diameter ~ "small", 
    Diameter >= median_diameter ~ "large")) %>% 
  ungroup() %>% 
  mutate(new_group = paste0(Class, sep = " ", size)) %>% 
  group_by(new_group) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(observation_count < 15))


#make config
configuration <- ppjsdm::Configuration(d$x_jitter,
                                       d$y_jitter, 
                                       types = d$new_group)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))

set.seed(1030)  
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

aic_fg_size$site[i] <- site

aic_fg_size$aic[i] <- fit$aic 

aic_fg_size$bic[i] <- fit$bic
}

aic_fg_size$model <- "fg_size"
write.csv(aic_fg_size, "aic_fg_size_t15.csv")

