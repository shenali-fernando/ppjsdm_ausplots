library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(patchwork)
library(scales)
library(austraits)
library(forcats)

set.seed(3)


###### All species - modelling species by size in each site
# Have some pretty good functions to help us here 
source("specifications/diameter/size_funs.R")
source("code/make_summary_fun.R")

#### Very similar to other model but increasing the threshold size to 15 and less 
#Because, yes, groups/points are lost but most sites this is a loss of very uncertain groups, or they were redunant 
#in the information they were giving, or there is no change by getting rid of these groups 
#This also means less post-model processing of individual sites for a cleaner model!

data <- read.csv("data/data_cleaned.csv")

sites <- unique(data$Site_Name)

full_df <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                       "hi_numerical", "Potential","site")

set.seed(3)
for (i in sites){
  
  site_mod <- size_sites(site = i, #single site 
                         group_type = "species_size",
                         show_size_freq = FALSE,
                         config_only = FALSE, #if TRUE, returns config only and exits function
                         threshold = 16, 
                         short_range = 10, 
                         short_model = "exponential")
  
  working_df <- make_sum_df(fits = list(site_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(site_mod$sum))
  
  working_df <- working_df %>% mutate(site = i)                              
  
  full_df<- rbind(full_df, working_df)
}



#Check unreasonable cis - 6 sites with one problem species each 

full_df<- full_df %>% 
  mutate(range_ci = hi - lo)

a <- full_df %>% filter(range_ci > 5)
a %>% count(site)



## Site = Bird
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
    mutate(new_group = paste0(Genus_Species, sep = " ", size_class)) %>% 
    group_by(new_group) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!(observation_count < 16)) %>% 
    filter(! new_group == "Acacia melanoxylon large")

#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
plot(configuration)

d %>% count(new_group)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
bird_fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
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


bird_sum <- summary(bird_fit)
bird_sum

bird_df <- make_sum_df(fits = list(bird_fit), 
                       summ = list(bird_sum))

bird_df <- bird_df %>% mutate(site = "Bird")

##Caveside
caveside <-  size_sites(site = "Caveside", #single site 
                    group_type = "species_size",
                    show_size_freq = FALSE,
                    config_only = FALSE, #if TRUE, returns config only and exits function
                    threshold = 24, 
                    short_range = 10, 
                    short_model = "exponential")

caveside_df <- make_sum_df(fits = list(caveside$fit), 
                       summ = list(caveside$sum))

caveside_df <- bird_df %>% mutate(site = "Caveside")

##Flowerdale 
flower <-  size_sites(site = "Flowerdale", #single site 
                        group_type = "species_size",
                        show_size_freq = FALSE,
                        config_only = FALSE, #if TRUE, returns config only and exits function
                        threshold = 18, 
                        short_range = 10, 
                        short_model = "exponential")

flower_df <- make_sum_df(fits = list(flower$fit), 
                           summ = list(flower$sum))


flower_df <- bird_df %>% mutate(site = "Flowerdale")


##Lardner
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
  mutate(new_group = paste0(Genus_Species, sep = " ", size_class)) %>% 
  group_by(new_group) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(observation_count < 16)) %>% 
  filter(! new_group == "Acacia melanoxylon large")

#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
plot(configuration)

d %>% count(new_group)

#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))


#fit model 
lardner_fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
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


lardner_sum <- summary(lardner_fit)
lardner_sum

lardner_df <- make_sum_df(fits = list(lardner_fit), 
                       summ = list(lardner_sum))

lardner_df <- lardner_df %>% mutate(site = "Lardner")

##MtField 
mtfield <-  size_sites(site = "MtField", #single site 
                      group_type = "species_size",
                      show_size_freq = FALSE,
                      config_only = FALSE, #if TRUE, returns config only and exits function
                      threshold = 19, 
                      short_range = 10, 
                      short_model = "exponential")

mtfield_df <- make_sum_df(fits = list(mtfield$fit), 
                         summ = list(mtfield$sum))

mtfield_df <- mtfield_df %>% mutate(site = "MtField")

##WaratahMix
waratah <-  size_sites(site = "WaratahMix", #single site 
                       group_type = "species_size",
                       show_size_freq = FALSE,
                       config_only = FALSE, #if TRUE, returns config only and exits function
                       threshold = 21, 
                       short_range = 10, 
                       short_model = "exponential")

waratah_df <- make_sum_df(fits = list(waratah$fit), 
                          summ = list(waratah$sum))

waratah_df <- waratah_df %>% mutate(site = "WaratahMix")


#Put all fixed sites into a dataframe 

fixed <- rbind(bird_df, caveside_df, flower_df, lardner_df, mtfield_df, waratah_df)
fixed<- fixed %>% 
  mutate(range_ci = hi - lo)

fixed %>% filter(range_ci > 5)
fixed %>% count(site)


#Get rid of fixed sites in full_df 
full_df <- full_df %>% 
  filter(! site %in% c("Bird", "Caveside", "Flowerdale", "Lardner", "MtField", "WaratahMix"))
  
full_df %>% count(site)

#Add fixed sites
full_df <- rbind(full_df, fixed)

full_df %>% count(site)


#Okay, we have a dataframe with threshold 15 now
#Let's add the visualisation columns needed 

#Get separate classes and species 
full_df2 <- full_df %>% 
  mutate(class_from = str_extract(from, "\\w+$")) %>% 
  mutate(class_to = str_extract(to, "\\w+$")) %>% 
  mutate(species_from = str_extract(from, "\\w+\\s+\\w+")) %>% 
  mutate(species_to = str_extract(to, "\\w+\\s+\\w+"))

#Make int col
full_df2 <- full_df2 %>% mutate(class_int = paste(class_from, sep = "_", class_to))

#Make georegion col 
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

full_df2 <- full_df2 %>% 
  mutate(region = case_when(georegion %in% c("S_VIC", "N_VIC", "S_NSW") ~ "SE_AUS", 
                            georegion %in% c("S_WA", "N_WA") ~ "WA",
                            georegion %in% c("QLD", "N_NSW") ~ "N_AUS",  
                            georegion %in% c("o_TAS", "d_TAS") ~ "TAS"))

         
#Add sig 
full_df2 <- full_df2 %>% 
  mutate(sig = ifelse(lo > 0 | hi < 0, 1 , NA))


#Get rid of Potential col
full_df2 <- full_df2 %>% 
  select(-Potential)


#Assign a species to canopy or subcanopy 
species_class <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/data/species_class.csv")
species_class <- species_class %>% dplyr::rename(species_to = species_from)

full_df2 <- left_join(full_df2, species_class, by = "species_to")
full_df2 <- full_df2 %>% rename(cc_to = Class)

full_df2 <- full_df2 %>% 
  mutate(cc_int = paste0(cc_from, sep = "_", cc_to))


##Add a functional group size column 
full_df2 <- full_df2 %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
    group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
    group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
    group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
    TRUE ~ group))


###Do a comparison of raw alpha, lo and hi estimate between this mod and previous (when threshold was <=12)
#Need to left_join properly due to missing rows - it's looking similar enough 

th_15 <- full_df2 %>% 
  select(from, to, site, alpha, lo, hi) 

th_12 <- df_all_3_5 %>% 
  select(from, to, site, alpha, lo, hi) 

df <- left_join(th_15, th_12, by = c("site", "from", "to")) 

df_diff <- df %>% 
  mutate(diff_alpha = alpha.x - alpha.y) %>% 
  mutate(diff_lo = lo.x - lo.y) %>% 
  mutate(diff_hi = hi.x - hi.y) %>% 
  mutate(prop_alpha = diff_alpha/alpha.y) %>% 
  mutate(prop_lo = diff_lo/lo.y) %>% 
  mutate(prop_hi = diff_hi/hi.y)


mean(df_diff$diff_alpha, na.rm = T)
mean(df_diff$diff_lo, na.rm = T)
mean(df_diff$diff_hi, na.rm = T)


mean(df_diff$prop_alpha, na.rm = T)
median(df_diff$prop_alpha, na.rm = T)
mean(df_diff$prop_lo, na.rm = T)
mean(df_diff$prop_hi, na.rm = T)



## Save out df 
write.csv(full_df2, "df_add3.5_t15.csv")


