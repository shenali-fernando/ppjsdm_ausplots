library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
library(ozmaps)


#### GEOGRAPHIC REGION ANALYSIS  

#### Doing a species by crown class analysis by splitting by geographic georegion 
#There are 2 ways to aggregate the sites into georegions: 
  #1. All sites of the georegion into one long line 
  #2. Stacking sites so that there are two lines of x sites 

#Data 
# data <- read.csv("data/data_cleaned.csv")
# 
# data <- data %>% 
#   mutate(georegion = Site_Name) %>% 
#   mutate(georegion = if_else(georegion %in% c("Dawson", "Frankland", "Clare", "Giants"), 
#                              "SWA", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("Carey", "Dombakup", "Warren",  "Sutton","Collins"), 
#                              "NWA", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek"), 
#                              "NVIC", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("Weeaproinah", "Turtons", "Lardner"), 
#                              "SVIC", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"), 
#                              "SE_NSW", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans"), 
#                              "N_NSW", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton"), 
#                              "QLD", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx"), 
#                              "d_TAS", georegion)) %>% 
#   mutate(georegion = if_else(georegion %in% c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Dip", "Flowerdale"),
#                              "o_TAS", georegion)) 



##### ONE LINE: GEOREGION 

#Make georegion function
make_georegion <- function(sites, #vector of sites 
                           group_type, #options of species, crown_class, species_crown_class 
                           threshold = 10, #if a group has less than 10 individuals remove
                           short_range = 10, #short-range interaction distance (set same for all species)
                           medium_range = 50, #medium-range interaction distance to specify SITE interactions v species/group interactions
                           long_range = 100, #long-range interaction distance to specify SITE interactions v species/group interactions
                           short_model = "exponential", 
                           medium_model = "tanh"){ 
  
  n <- length(sites) #number of sites 
 
  
  #we need to bind all the data into one big georegion site before anything else 
  
  final_df <- data.frame()  # empty df to store results
  
  for (i in seq_along(sites)) {
    temp <- data %>%
     filter(Site_Name == sites[i]) %>%
      mutate(new_x = Ausplot_X + (i - 1) * 100,  # add increment
             new_y = Ausplot_Y)
    
    final_df <- rbind(final_df, temp)  # combine with previous
  }
  
  #jitter coordinates
  d <- final_df %>%
    group_by(new_x, new_y) %>% #group by coordinate columns
    mutate(
      is_duplicated = n() > 1, #create column of TRUE/FALSE 
      #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
      x_jitter = if_else(is_duplicated, new_x + runif(n(), -0.025, 0.025), new_x), #create x_jitter column
      y_jitter = if_else(is_duplicated, new_y + runif(n(), -0.025, 0.025), new_y) #create y_jitter column
    ) %>%
    ungroup() 
  
  #run group_type 
   if(group_type == "species"){ #does not distinguish eucs and non-eucs
     d <- d %>% 
       mutate(new_group = Genus_Species) %>% 
       group_by(Genus_Species) %>% 
       mutate(observation_count = n()) %>% 
       ungroup() %>% 
       filter(!new_group == "Unidentified tree") %>% 
       filter(!(observation_count < threshold))
     
     
   } else if(group_type == "crown_class"){
     d <- d %>% 
       mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
       mutate(new_group = paste(species, Crown_Class)) %>% 
       mutate(new_group = if_else(new_group == "Eucalyptus Emergent", "Eucalyptus Dominant", new_group)) %>%  #thresholds for groups 
       group_by(new_group) %>% 
       mutate(observation_count = n()) %>% 
       ungroup() %>% 
       filter(!(observation_count < threshold))
     
   } else if(group_type == "species_crown_class"){
     d <- d %>% 
       mutate(group = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
       mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
       mutate(new_group = paste(group, Crown_Class)) %>% 
       group_by(new_group) %>% 
       mutate(observation_count = n()) %>% 
       ungroup() %>% 
       filter(!(observation_count < threshold)) 
   } else
   {return("Error: group_type not recognised; check specification")}
    
  #make config
  configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
  plot(configuration)
  
  #set parameters 
  window <- ppjsdm::Rectangle_window(x_range = c(0, round(max(d$x_jitter), digits = -2)), 
                                     y_range = c(0, ceiling(max(d$y_jitter))))
  
  nspecies <- length(levels(configuration$types))
  
  #fit model 
  fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                       window = window,
                       short_range = matrix(short_range, nspecies, nspecies), 
                       medium_range = matrix(medium_range, nspecies, nspecies),
                       long_range = matrix(long_range, nspecies, nspecies),
                       model = short_model,
                       medium_range_model = medium_model,
                       saturation = 10, 
                       nthreads = 4, 
                       fitting_package = "glmnet",
                       dummy_distribution = "stratified",
                       min_dummy = 1, dummy_factor = 1e10, 
                       max_dummy = 1e3)
  
  fit$coefficients$alpha
  sum <- summary(fit)
  sum
  return(list(fit = fit, sum = sum))
}



##### TWO LINE: GEOREGION 

# 

plotlist <- list()

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

  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
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
