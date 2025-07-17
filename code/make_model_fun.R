library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

#### Crown class model function for making model csv 
site_ccmodel <- function(site, 
                       group_cd = FALSE, 
                       threshold_noneuc = 10, 
                       threshold_euc = 5, 
                       short_range = 8){
  
d <- data_cleaned %>% filter(Site_Name == site)


if(group_cd == TRUE){ #group all eucs together, group C and D together 
d <- d %>% 
  mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
  mutate(new_cc = paste(species, Crown_Class))%>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus Co-dominant", "Eucalyptus Dominant"), "Eucalyptus Co/dominant", new_cc)) } 
else{
  d <- d %>% 
    mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
    mutate(new_cc = paste(species, Crown_Class)) 
}


d <- d %>%  #thresholds for groups 
  group_by(new_cc) %>% 
  mutate(observation_count = n()) %>% 
  filter(!(str_starts(new_cc, "Non-euc") & observation_count < threshold_noneuc)) %>% 
  filter(!(str_starts(new_cc, "Eucalyptus") & observation_count < threshold_euc)) 


counts <- d %>% count(new_cc)
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


configuration<- Configuration(d$Ausplot_X, d$Ausplot_Y, types = d$new_cc)

# Set parameters
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
nspecies <- length(levels(configuration$types))
short_range <- matrix(short_range, nspecies, nspecies)
model <- "exponential"

fit<- ppjsdm::gibbsm(configuration = configuration, 
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


box <- ppjsdm::box_plot(fit = fit,
                 summ = sum,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)
print(sum)
print(box)

return(list(fit = fit, sum = sum))

}



#test function
anu101 <- site_ccmodel(site = "ANU101")
anu101$fit
anu101$sum
