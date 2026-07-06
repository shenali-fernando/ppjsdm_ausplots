library(dplyr)
library(spatstat)
library(ppjsdm)

set.seed(20)

######### AUC - Check the AUC of all models (null, species, size, fg, species x size, fg x size)

## AUC Function

source("functions/auc_function.R")

## Load Data 
data <- read.csv("data/ausplots/data_cleaned.csv")
fg <- read.csv("data/ausplots/species_class.csv")


#Okay, so let's start with the simpler models of species, size and fg 
# The inputs for the AUC is window, fit and configuration. Therefore this is done at the site level (or each fit level)
# Basically, just need to adapt the loop code to generate the coefficients to get the AUC instead 


#Set window because it's always the same 

window_o <- owin(xrange = c(0, 100), 
                 yrange = c(0, 100))

###########################################################################################
####################### SPECIES MODEL #####################


sp_model <- function(site, 
                     threshold = 15, #threshold to exclude groups
                     short_range = 10){ #short range value 
  
  d <- data %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    group_by(Genus_Species) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    mutate(group = if_else(observation_count < threshold, "Misc", Genus_Species)) %>% 
    group_by(group) %>% 
    mutate(obs = n()) %>% 
    ungroup() %>% 
    filter(! obs < threshold)
    
  
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
  short_range <- matrix(short_range, nspecies, nspecies)
  model <- "exponential"
  
set.seed(20)

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

  
  return(list(fit = fit, configuration = configuration)) #return both fit and sum to extract 
  
}



sites <- unique(data$Site_Name) 

#make empty df to fill auc to 
auc_full <- data.frame(type = as.character(), 
                       auc = numeric(), 
                       site = as.character())

for(i in sites){ 
  
  sp_mod <- sp_model(site = i, 
                     threshold= 15, #threshold to exclude groups
                     short_range = 10) #short range always ten
  
 auc <- conditional_auc(window = window_o,
                        configuration = sp_mod$configuration,
                        fit = sp_mod$fit)
  
 auc$site <- i #add site name 
 
 auc_full <- rbind(auc_full, auc)
 
 }


#save out 

write.csv(auc_full, "auc_species_model_t15.csv")



##################################################################################
############################# SIZE MODEL #########################################


size_site <- function(site, 
                      threshold = 15, #threshold to exclude groups
                      short_range = 10){ #short range value 
  
  d <- data %>% filter(Site_Name == site)
  
  d <- d %>% 
    filter(!is.na(Diameter)) %>% 
    mutate(size = if_else(Diameter < (mean(Diameter) + 3.5), "small", "large")) %>% 
    group_by(size) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!observation_count < threshold) 
  
  
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
  short_range <- matrix(short_range, nspecies, nspecies)
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

  
  return(list(fit = fit, configuration = configuration)) #return both fit and sum to extract 
  
}

sites <- unique(data$Site_Name) 

auc_size_full <- data.frame(type = as.character(), 
                       auc = numeric(), 
                       site = as.character())

for(i in sites){ 
  
  size_mod <- size_site(site = i, 
                             threshold = 15, #threshold to exclude euc groups
                             short_range = 10) #run site_ccmodel function
  
  auc <- conditional_auc(window = window_o,
                         configuration = size_mod$configuration,
                         fit = size_mod$fit)
  
  auc$site <- i #add site name 
  
  auc_size_full <- rbind(auc_size_full, auc)
}

auc_size_full %>% count(site) #threshold not in use so no difference 

write.csv(auc_size_full, "auc_size_model_t15.csv") 


####################################################################################
########################### Functional Group Model #################################

fg <- fg %>% 
  mutate(Genus_Species = iconv(Genus_Species, to = "ASCII", sub = " "))

data1 <- left_join(data, fg, 
                   by = "Genus_Species")

data1 <- data1 %>%
  mutate(Class = if_else(Genus_Species == "Unidentified tree", "Unid tree",Class)) %>% 
  mutate(Class = if_else(is.na(Class), "Subcanopy", Class)) #everything not in the df is subcanopy

data1 %>% count(Class, Site_Name)



fg_fun <- function(site, #single site 
                    threshold = 15, 
                    short_range = 10, 
                    short_model = "exponential"){ #threshold to exclude a group 
  
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
  
  return(list(fit = fit, configuration = configuration))
}

sites <- unique(data$Site_Name) 

auc_fg <- data.frame(type = as.character(), 
                            auc = numeric(), 
                            site = as.character())

for(i in sites){ 
  
  fg_mod <- fg_fun(site = i,
                        threshold = 15, #threshold to exclude groups, not in use here 
                        short_range = 10) #run site_ccmodel function
  
  auc <- conditional_auc(window = window_o,
                         configuration = fg_mod$configuration,
                         fit = fg_mod$fit)
  
  auc$site <- i #add site name 
  
  auc_fg <- rbind(auc_fg, auc)
}

auc_fg %>% count(site)

write.csv(auc_fg, "auc_fg_model_t15.csv")
#Unid trees included as own class 






###################################################################################
################################ FG x Size Model ##################################
#it seems if there are low numbers in an FG + Size group, its way under 10 so changing threshold has minimal impact

fg_size <- function(site, #single site 
                          threshold = 15, #threshold to exclude group  
                          short_range = 10, 
                          short_model = "exponential"){  
  
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
    filter(!(observation_count < threshold)) 
  
  
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
  
  return(list(fit = fit, configuration = configuration))
}


sites <- unique(data$Site_Name) 

auc_fg_size <- data.frame(type = as.character(), 
                     auc = numeric(), 
                     site = as.character())

for(i in sites){ 
  
  fg_size_mod <- fg_size(site = i,
                   threshold = 15, #threshold to exclude euc groups
                   short_range = 10) #run site_ccmodel function
  
  auc <- conditional_auc(window = window_o,
                         configuration = fg_size_mod$configuration,
                         fit = fg_size_mod$fit)
  
  auc$site <- i #add site name 
  
  auc_fg_size <- rbind(auc_fg_size, auc)
}

auc_fg_size %>% count(site)


write.csv(auc_fg_size, "auc_fg_size_model_t15.csv")



###################################################################################
########################### SPECIES x SIZE MODEL #################################

#Most complicated moodel... 

size_species_fun <- function(site, 
                       threshold = 15, 
                       short_range = 10, 
                       short_model = "exponential"){ #threshold to exclude a group 
  
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
      mutate(size_class2 = if_else(observation_count < threshold, 
                                   "Misc", 
                                   new_group)) %>% 
      group_by(size_class2) %>% 
      mutate(obs = n()) %>% 
      ungroup() %>% 
      filter(!(obs < threshold)) 

  

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
  

  return(list(fit = fit, configuration = configuration))
}

sites <- unique(data$Site_Name)

auc_species_size <- data.frame(type = as.character(), 
                          auc = numeric(), 
                          site = as.character())

for (i in sites){
  
  species_size_mod <- size_species_fun(site = i, 
                                       threshold = 15, 
                                       short_range = 10, 
                                       short_model = "exponential")

  
  auc <- conditional_auc(window = window_o,
                         configuration = species_size_mod$configuration,
                         fit = species_size_mod$fit)
  
  auc$site <- i #add site name 
  
  auc_species_size <- rbind(auc_species_size, auc)
  
}

auc_species_size %>% count(site)

write.csv(auc_species_size, "auc_species_size_model_t15.csv")

#the problem is that i manually removed groups with a range of CI larger than 5 
# this is not possible in the above code 