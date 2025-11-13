library(ppjsdm)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

#Load Data
data <- read.csv("data/data_cleaned.csv")


#What to be done: 
# for each species in a site, we need to spilt each species so that there are roughly even numbers of 'small' and 'large' individuals 
#for all individuals, not just eucs 

#This part comes under the data cleaning before fitting a model for each site
#So, we need a function that for each site produces:
#spilts individuals into small and large and fits and compute the model summary 

#These are the sites
sites <- unique(data$Site_Name)


########################################################################################################
## Output fit and summary for a site when using size classes 
size_sites <- function(site, #single site 
                       group_type = c("species_size", "species"),
                       show_size_freq = TRUE,
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
  
  
  #run group_type 
  if(group_type == "species"){ #does not distinguish eucs and non-eucs
    d <- df %>% 
      mutate(new_group = Genus_Species) %>% 
      group_by(Genus_Species) %>% 
      mutate(observation_count = n()) %>% 
      ungroup() %>% 
      filter(!new_group == "Unidentified tree") %>% 
      filter(!(observation_count < threshold))
  } 
  else if(group_type == "species_size"){
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
      filter(!(observation_count < threshold)) %>% 
      filter(! new_group %in% c("Eucalyptus urnigera small"))
  } 
  else
  {return("Error: group_type not recognised; check specification")}
 
  if(show_size_freq){
    d %>% 
      ggplot(aes(x = Diameter, group = Genus_Species, fill = Genus_Species)) +
      geom_density(adjust =1.5, alpha = 0.5) + 
      theme_bw() + 
      xlim(c(0,120))
  }

  #make config
  configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
  plot(configuration)
  
  if(config_only){
    return(configuration)
  }
  
  d %>% count(new_group)
  
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
                       fitting_package = "glmnet",
                       dummy_distribution = "stratified",
                       min_dummy = 1, dummy_factor = 1e10, 
                       max_dummy = 1e3)

  
  sum <- summary(fit)
  sum
  return(list(fit = fit, sum = sum))
}



############################################################################################################
##### Run AIC for each type in a site 
size_sites_aic <- function(site, #vector of sites 
                           threshold = 10,
                           group_type = c("species_size", "species")){
  
  plotlist <- list()
  
  df <- data %>% 
    filter(Site_Name == site)
  
  #jitter coordinates t
  df <- df %>%
    mutate(
      is_duplicated = n() > 1, #create column of TRUE/FALSE 
      #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
      x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
      y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
    ) 
  
  
  #run group_type 
  if(group_type == "species"){ #does not distinguish eucs and non-eucs
    d <- df %>% 
      mutate(new_group = Genus_Species) %>% 
      group_by(Genus_Species) %>% 
      mutate(observation_count = n()) %>% 
      ungroup() %>% 
      filter(!new_group == "Unidentified tree") %>% 
      filter(!(observation_count < threshold))
  } 
  else if(group_type == "species_size"){
    d <- df %>% 
      filter(!Genus_Species == "Unidentified tree") %>% 
      group_by(Genus_Species) %>% 
      mutate(observation_count = n()) %>% 
      filter(!(observation_count < 10)) %>%  #hard-code a threshold as 10 here to make it easier 
      mutate(median_diameter = ceiling(median(Diameter))) %>% 
      ungroup() %>% 
      mutate(size_class = case_when(
        Diameter < median_diameter ~ "small", 
        Diameter >= median_diameter ~ "large"
      )) %>% 
      mutate(new_group = paste0(Genus_Species, sep = " ", size_class)) %>% 
      filter(!(observation_count < threshold))
  } 
  else
  {return("Error: group_type not recognised; check specification")}
  
  configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
  plot(configuration)
  
  #set parameters 
  window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                     y_range = c(0, 100))
  
  nspecies <- length(levels(configuration$types))
  
  
  #aic optimisation
  
  for (h in levels(configuration$types)) {
    to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
      sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe
        set.seed(1)
        fit <- ppjsdm::gibbsm(configuration[h], #create the fit
                              window = window,
                              model = df$model[i],
                              short_range = matrix(df$short[i]),
                              saturation = 10,
                              dummy_distribution = "stratified",
                              min_dummy = 1, max_dummy =1e3,
                              dummy_factor = 1e10,
                              nthreads = 4,
                              fitting_package = "glmnet")
        
        fit$aic #take the aic value from the fit
      })
    }
    
    possible_short <- seq(from = 1, to = 25, length.out = 50) #possible short_range values
    possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
    df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
    df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
    df$potentials <- df$model
    
    plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
    
    plotlist[[h]] <- plot
  }
  
  for (p in plotlist) {
    print(p)
  }
  
}

# 
# 
# size_sites_aic(site = "Turtons", 
#                group_type = "species_size")


#########################################################################################
## Optimising the fit with aic 
# This function is meant to be run after the AIC function, therefore you MUST have an idea of what parameters work 

fit_opt <- function(configuration, #config to opt 
                    short_range = c(8, 10, 12), #range of short-range interaction radii values
                    medium_range = c(30, 40, 50), 
                    long_range = c(80, 100, 120), 
                    short_model = c("exponential", "bump"), #model 
                    medium_model= c("exponential", "bump"),
                    print_df = TRUE){ #if TRUE, the whole df of aic values is returned, is FALSE only parameters for minimum AIC is returned 
  
  
  ntypes <- length(levels(configuration$types))
  
  df <- expand.grid(short = short_range, 
                    medium = medium_range,
                    long = long_range, 
                    short_model = short_model, 
                    medium_model = medium_model)
  
  to_optimize2 <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration, #create the fit
                            window = window,
                            model = df$short_model[i],
                            medium_range_model = df$medium_model[i],
                            short_range = matrix(df$short[i], ntypes, ntypes),
                            medium_range = matrix(df$medium[i], ntypes, ntypes),
                            long_range = matrix(df$long[i], ntypes, ntypes),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4,
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  
  df$aic <- to_optimize2(df)
  
  if(print_df){
  return(df)}
  else{
    min <- which.min(df$aic)
    return(df[min, ])
   }
}





