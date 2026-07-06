library(ppjsdm)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

#Load Data
# data <- read.csv("data/data_cleaned.csv")
# 
# #What to be done: 
# # for each species in a site, we need to spilt each species so that there are roughly even numbers of 'small' and 'large' individuals 
# #for all individuals, not just eucs 
# 
# #This part comes under the data cleaning before fitting a model for each site
# #So, we need a function that for each site produces:
# #spilts individuals into small and large and fits and compute the model summary 
# 
# #These are the sites
# sites <- unique(data$Site_Name)


########################################################################################################
## Output fit and summary for a site when using size classes 
size_sites <- function(site, #single site =
                       show_size_freq = FALSE,
                       config_only = FALSE, #if TRUE, returns config only and exits function
                       threshold = 10, 
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
 
 
  if(show_size_freq){
    d %>% 
      ggplot(aes(x = Diameter, group = Genus_Species, fill = Genus_Species)) +
      geom_density(adjust =1.5, alpha = 0.5) + 
      theme_bw() + 
      xlim(c(0,120))
  }

  #make config
  configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$species_size2)
  plot(configuration)
  
  if(config_only){
    return(configuration)
  }
  
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

  
  sum <- summary(fit)
  sum
  return(list(fit = fit, sum = sum, config = configuration))
}


