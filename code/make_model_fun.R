library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

#Load cleaned data 
data_cleaned <- read.csv("data/data_cleaned.csv")


#function to extract number of individuals in a site 
find_ind <- function(site, 
                     group_cd = TRUE, 
                     threshold_noneuc = 8, 
                     threshold_euc = 8){
  
    d <- data_cleaned %>% filter(Site_Name == site)
    
    
    if(group_cd == TRUE){ #group all eucs together, group C and D together 
      d <- d %>% 
        mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
        mutate(new_cc = paste(species, Crown_Class))%>% 
        mutate(new_cc = if_else(new_cc == "Eucalyptus Emergent", "Eucalyptus Dominant", new_cc)) %>% 
        mutate(new_cc = if_else(new_cc %in% c("Eucalyptus Co-dominant", "Eucalyptus Dominant"), "Eucalyptus Co/dominant", new_cc)) %>% 
        mutate(new_cc = if_else(new_cc %in% c("Non-euc Co-dominant", "Non-euc Dominant"), "Non-euc Co/dominant", new_cc))} 
    else{
      d <- d %>% 
        mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
        mutate(new_cc = paste(species, Crown_Class)) %>% 
        mutate(new_cc = if_else(new_cc == "Eucalyptus Emergent", "Eucalyptus Dominant", new_cc)) 
    }
    
    
    d <- d %>%  #thresholds for groups 
      group_by(new_cc) %>% 
      mutate(observation_count = n()) %>% 
      filter(!(str_starts(new_cc, "Non-euc ") & observation_count < threshold_noneuc))
    
    d <- d %>% 
      filter(!(str_starts(new_cc, "Eucalyptus ") & observation_count < threshold_euc)) 
    
    
    counts <- d %>% count(new_cc)
    return(counts) # print c
    
}



spcc_ind <- function(site, 
                     threshold_noneuc = 10, 
                     threshold_euc = 10){
  
  d <- data_cleaned %>% filter(Site_Name == site)
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
    mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
    mutate(new_spcc = paste(species, Crown_Class))  
  
  
  d <- d %>%  #thresholds for groups 
    group_by(new_spcc) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!(str_starts(new_spcc, "Non-euc") & observation_count < threshold_noneuc)) %>% 
    filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < threshold_euc))
  
  counts <- d %>% count(new_spcc)
  return(counts) # print c
  
}

#Test
# weld <- find_ind(site = 'Weld', 
#                  group_cd = FALSE, 
#                  threshold_noneuc = 12, 
#                  threshold_euc = 5)



# Create vector of Site Names 
sites <- unique(data_cleaned$Site_Name)
sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]

#make empty df 
inds_full <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(inds_full) <- c("new_cc", "n", "site_name")

for(i in sites){ #do all sites at once and send into df 
  
  inds <- find_ind(site = i, 
                   group_cd = TRUE, 
                   threshold_noneuc = 8, 
                   threshold_euc = 8) #function above
  
  inds <- inds %>% mutate(site_name = i) #set new column called site_name for site 
  
  inds_full <- rbind(inds_full, inds) #bind into larger df 
}



#### Crown class model function for making model summary and fit  

site_ccmodel <- function(site, 
                       group_cd = TRUE, #should we group c and d together
                       threshold_noneuc = 12, #threshold to exclude non-euc groups
                       threshold_euc = 7, #threshold to exclude euc groups
                       short_range = 8){ #short range value 
  
d <- data_cleaned %>% filter(Site_Name == site)


if(group_cd == TRUE){ #group all eucs together, group C and D together 
d <- d %>% 
  mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
  mutate(new_cc = paste(species, Crown_Class))%>% 
  mutate(new_cc = if_else(new_cc == "Eucalyptus Emergent", "Eucalyptus Dominant", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Eucalyptus Co-dominant", "Eucalyptus Dominant"), "Eucalyptus Co/dominant", new_cc)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Non-euc Co-dominant", "Non-euc Dominant"), "Non-euc Co/dominant", new_cc))} 
else{
  d <- d %>% 
    mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
    mutate(new_cc = paste(species, Crown_Class)) %>% 
    mutate(new_cc = if_else(new_cc == "Eucalyptus Emergent", "Eucalyptus Dominant", new_cc))
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


configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_cc) #make config 

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


box <- ppjsdm::box_plot(fit = fit,
                 summ = sum,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)
print(box)
return(list(fit = fit, sum = sum)) #return both fit and sum to extract 

}



#test function
site_ccmodel(site = "BirdTree")
anu101$fit
anu101$sum


## put this into a loop to extract the df 


# Create vector of Site Names 
sites <- unique(data_cleaned$Site_Name)
sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]

#make empty df 
full_df_com <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df_com) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                         "hi_numerical", "Potential", "site")

for(i in sites){ 
  
 ccmodel <- site_ccmodel(site = i, 
                         group_cd = TRUE, #should we group c and d together
                         threshold_noneuc = 8, #threshold to exclude non-euc groups
                         threshold_euc = 8, #threshold to exclude euc groups
                         short_range = 8) #run site_ccmodel function
 
 working_df <- make_sum_df(fits = list(ccmodel$fit), #use make_summary_df function to output summary as df 
                           summ = list(ccmodel$sum))
 
 working_df <- working_df %>% mutate(site = i) #set new column for site name 
 
 full_df_com <- rbind(full_df_com, working_df) #bind into larger df 
  
  }



## make new columns class and species for easier plotting in ggplot 

full_df_com <- full_df_com %>% 
  mutate(class_from = gsub("^[^ ]+ ", "\\1", from)) %>% 
  mutate(class_to = gsub("^[^ ]+ ", "\\1", to)) %>% #same for euc or non-euc 
  mutate(species_from = gsub("^([^ ]+).*", "\\1", from)) %>% 
  mutate(species_to = gsub("^([^ ]+).*", "\\1", to))

full_df_com <- full_df_com %>% ## and move site to the front to make it easier 
  relocate(site)

full_df_com %>% count(from == 'Non-euc Suppressed')


#### To add the number of individuals into the large df 
inds_full <- inds_full %>% 
  rename(to = from) 

full_df_com <- full_df_com %>%
  left_join(inds_full, by = c("site", "from"))

full_df_com <- full_df_com %>% 
  rename(ind_from = n)

full_df_com <- full_df_com %>% 
  left_join(inds_full, by = c("site", "to")) %>% 
  rename(ind_to = n)



full_df_com <- full_df_com %>% 
  mutate(range_ci = hi - lo)


a <- full_df_com %>% filter(species_from == "Non-euc") %>% 
  filter(species_to == "Non-euc")









#Make by species 

site_spmodel <- function(site, 
                         group_cd = TRUE, #should we group c and d together
                         threshold_noneuc = 12, #threshold to exclude non-euc groups
                         threshold_euc = 7, #threshold to exclude euc groups
                         short_range = 8){ #short range value 
  
  d <- data_cleaned %>% filter(Site_Name == site)
  
  
  if(group_cd == TRUE){ #group all eucs together, group C and D together 
    d <- d %>% 
      mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia||^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
      mutate(new_cc = paste(Genus_Species, Crown_Class))%>% 
      mutate(new_cc = if_else(new_cc %in% c("Eucalyptus Co-dominant", "Eucalyptus Dominant", "Eucalyptus Emergent"), "Eucalyptus Co/dominant", new_cc)) %>% 
      mutate(new_cc = if_else(new_cc %in% c("Non-euc Co-dominant", "Non-euc Dominant"), "Non-euc Co/dominant", new_cc)) %>% 
      mutate(group = paste(Genus_Species, str_extract(new_cc, "\\S+$")))} 
  else{
    d <- d %>% 
      group = paste(Genus_Species, new_cc)
  }
  
  
  d <- d %>%  #thresholds for groups 
    group_by(group) %>% 
    mutate(observation_count = n()) %>% 
    filter(!observation_count < threshold_euc) 
  
  
  counts <- d %>% count(group)
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
  
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_cc) #make config 
  
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
  
  
  box <- ppjsdm::box_plot(fit = fit,
                          summ = sum,
                          coefficient = "alpha",
                          which = "all", 
                          text_size = 10)
  print(box)
  return(list(fit = fit, sum = sum)) #return both fit and sum to extract 
  
}

#test 
site_spmodel(site = "Weld")
