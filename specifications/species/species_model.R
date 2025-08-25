library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)


#Load cleaned data 
data_cleaned <- read.csv("data/data_cleaned.csv")

#Source other functions 
source("code/make_summary_fun.R")

##### Species Model
## This is the species model - it is actually very similar to the null model for eucs only
## As most sites have only one species of eucalypt, therefore I doubt there is much change from the only euc null model 
## However, have to see if it is useful 

# Make a function to create a fit and summary for each site 
sp_model <- function(site, 
                       threshold_noneuc = 10, #threshold to exclude non-euc groups
                       threshold_euc = 10, #threshold to exclude euc groups
                       short_range = 8){ #short range value 
  
  d <- data_cleaned %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
    group_by(species) %>% 
    mutate(observation_count = n()) %>% 
    filter(!(str_starts(species, "Non-euc") & observation_count < threshold_noneuc)) %>% 
    filter(!(str_starts(species,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < threshold_euc)) 
  
  
  counts <- d %>% count(species)
  print(counts) # print counts to see 
  
  d <- d %>%
    group_by(Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
    mutate(
      is_duplicated = n() > 1, #create column of TRUE/FALSE 
      #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
      x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
      y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
    ) %>%
    ungroup() 
  
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$species) #make config 
  
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
  
  sum <- summary(fit)
  
  return(list(fit = fit, sum = sum)) #return both fit and sum to extract 
  
}




sites <- unique(data_cleaned$Site_Name) 
#sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]

full_df_sp <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df_sp) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                            "hi_numerical", "Potential", "site")

for(i in sites){ 
  
  sp_mod <- sp_model(site = i, 
                         threshold_noneuc = 10, #threshold to exclude non-euc groups
                         threshold_euc = 10, #threshold to exclude euc groups
                         short_range = 8) #run site_ccmodel function
  
  working_df <- make_sum_df(fits = list(sp_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(sp_mod$sum))
  
  working_df <- working_df %>% mutate(site = i) #set new column for site name 
  
  full_df_sp <- rbind(full_df_sp, working_df) #bind into larger df 
  
}


##Okay need to check unreasonable CIs 
full_df_sp2 <- full_df_sp2 %>% 
  mutate(range_ci = hi - lo)

large_ci <- full_df_sp2 %>% 
  filter(range_ci > 5) #There's nothing!! yay


### Some visualisation 

### First, only eucs to eucs 
eucs_sp <- full_df_sp2 %>% 
  filter(! from == "Non-euc") %>% 
  filter(! to == "Non-euc")

eucs_sp %>% count(to)

### Within-species interactions 

within_eucs_sp <- eucs_sp %>% 
  filter(from == to)

ggplot(data = eucs_sp, 
       aes(x = from, 
           y = alpha, 
           colour = site)) + 
  geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#well, that's a lot to look at 

#