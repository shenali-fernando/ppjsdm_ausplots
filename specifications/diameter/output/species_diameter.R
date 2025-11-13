library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(patchwork)
library(scales)

##################
### Small v Large Diameter Trees 

#Based on analyses of diameter and where it splits the crown classes, choosing 50cm at the difference between large and small trees 

#Load cleaned dataset
data_cleaned <- read.csv("data/data_cleaned.csv")

#Source other functions 
source("code/make_summary_fun.R")

##Add small v large - based on the diameter work a small tree is <50cm, big tree is >=50cm 
data <- data_cleaned %>% 
  mutate(size = if_else(Diameter < 50, "small", "large")) %>% 
  mutate(group = paste(size, sep = " ", Genus_Species))


### Okay, now run the models - again do this is a long loop 
size_model <- function(site, 
                     threshold = 10, #threshold to exclude euc groups
                     short_range = 8){ #short range value 
  
  d <- data %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
    mutate(new_group = paste(size, sep = " ", species)) %>% 
    group_by(new_group) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!observation_count < threshold) 
  
  counts <- d %>% count(new_group)
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
  
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_group) #make config 
  
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

full_df_size <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df_size) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                          "hi_numerical", "Potential", "site")

for(i in sites){ 
  
  size_mod <- size_model(site = i, 
                     threshold = 10, #threshold to exclude euc groups
                     short_range = 8) #run site_ccmodel function
  
  working_df <- make_sum_df(fits = list(size_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(size_mod$sum))
  
  working_df <- working_df %>% mutate(site = i) #set new column for site name 
  
  full_df_size <- rbind(full_df_size, working_df) #bind into larger df 
  
}


##Okay need to check unreasonable CIs 
full_df_size2 <- full_df_size %>% 
  mutate(range_ci = hi - lo)

large_ci <- full_df_size2 %>% 
  filter(range_ci > 5) 

##################### Fixing fits of specific sites 
### There are four problem sites: MtField, Bruxner, Lorne, WaratahMix. We run these separately 
## MtField 
d <- data %>% filter(Site_Name == "MtField")

d <- d %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(new_group = paste(size, sep = " ", species)) %>% 
  group_by(new_group) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!observation_count < 10) %>% 
  filter(! new_group == "small Eucalyptus coccifera")

counts <- d %>% count(new_group)
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


configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_group) #make config 

# Set parameters
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
nspecies <- length(levels(configuration$types))
short_range <- matrix(8, nspecies, nspecies)
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

mtfield <- make_sum_df(fits = list(fit), #use make_summary_df function to output summary as df 
            summ = list(sum))

mtfield <- mtfield %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "MtField")

### Lorne 
d <- data %>% filter(Site_Name == "Lorne")

d <- d %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(new_group = paste(size, sep = " ", species)) %>% 
  group_by(new_group) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!observation_count < 21)

counts <- d %>% count(new_group)
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


configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_group) #make config 

# Set parameters
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
nspecies <- length(levels(configuration$types))
short_range <- matrix(8, nspecies, nspecies)
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

lorne <- make_sum_df(fits = list(fit), #use make_summary_df function to output summary as df 
                       summ = list(sum))

lorne <- lorne %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Lorne") 
  
  

### Bruxner 
d <- data %>% filter(Site_Name == "Bruxner")

d <- d %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(new_group = paste(size, sep = " ", species)) %>% 
  group_by(new_group) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!observation_count < 11)

counts <- d %>% count(new_group)
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


configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_group) #make config 

# Set parameters
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
nspecies <- length(levels(configuration$types))
short_range <- matrix(8, nspecies, nspecies)
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

bruxner <- make_sum_df(fits = list(fit), #use make_summary_df function to output summary as df 
                     summ = list(sum))

bruxner <- bruxner %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Bruxner") 

#Okay, I'm happy with the changes to these sites 



#Remove sites from full df 
full_df_size2 <- full_df_size2 %>% 
  filter(!site %in% c("MtField", "Bruxner", "Lorne"))

#Add sites back in 
final_sizemod <- rbind(full_df_size2, mtfield, bruxner, lorne)

#check sites 
final_sizemod %>% count(site) #looks right 



############################## Visualisation 
#Need to add a sig column 
final_sizemod <- final_sizemod %>% 
  mutate(sig = ifelse((lo > 0 | hi < 0), 1, NA ))

#make columns easy for plotting
final_sizemod <- final_sizemod %>% 
  relocate(site) %>% 
  mutate(species_from = gsub("^[^ ]+ ", "\\1", from)) %>% 
  mutate(species_to = gsub("^[^ ]+ ", "\\1", to)) %>% #same for euc or non-euc 
  mutate(size_from = gsub("^([^ ]+).*", "\\1", from)) %>% 
  mutate(size_to = gsub("^([^ ]+).*", "\\1", to))

## check species names 
final_sizemod %>% count(species_from) #okay need to fix the Syncarpia names 

final_sizemod <- final_sizemod %>% 
  mutate(species_from = if_else(species_from %in% c("Syncarpia glomulifera", 
                                                    "Syncarpia glomulifera subsp glabra", 
                                                    "Syncarpia glomulifera subsp glomulifera"), 
                                "Syncarpia glomulifera", species_from)) %>% 
  mutate(species_to = if_else(species_to %in% c("Syncarpia glomulifera", 
                                                    "Syncarpia glomulifera subsp glabra", 
                                                    "Syncarpia glomulifera subsp glomulifera"), 
                                "Syncarpia glomulifera", species_to)) 

#cool - let me save out this 
write.csv(final_sizemod, "final_sizemod.csv")




final_sizemod <- final_sizemod %>% 
  mutate(region2 = region) %>% 
  mutate(region2 = if_else(region2 %in% c("d_TAS", "o_TAS"), 
                           "TAS", region2)) %>% 
  mutate(region2 = if_else(region2 %in% c("N_NSW", "QLD"), 
                           "N_AUS", region2)) %>%
  mutate(region2 = if_else(region2 %in% c("SE_NSW", "SE_VIC"), 
                           "SE_AUS", region2)) 


#Save out df 
write.csv(final_sizemod, "final_sizemod.csv")

final_sizemod <- final_sizemod %>% 
  mutate(region = site) %>% 
  mutate(region = if_else(region %in% c("Carey", "Dombakup", "Warren", "Dawson", "Giants", "Sutton", "Frankland", "Clare", "Collins"), 
                          "WA", region)) %>% 
  mutate(region = if_else(region %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Weeaproinah", "Turtons", "Lardner"), 
                          "SE_VIC", region)) %>% 
  mutate(region = if_else(region %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"), 
                          "SE_NSW", region)) %>% 
  mutate(region = if_else(region %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans"), 
                          "N_NSW", region)) %>% 
  mutate(region = if_else(region %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton"), 
                          "QLD", region)) %>% 
  mutate(region = if_else(region %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx"), 
                          "d_TAS", region)) %>% 
  mutate(region = if_else(region %in% c("BondTier", "BlackRiver", "Weld", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Dip", "Flowerdale"), 
                          "o_TAS", region)) 



############# Visualisation 
### EUCS TO EUCS 
eucs <- final_sizemod %>% 
  filter(! species_from == "Non-euc") %>% 
  filter(! species_to == "Non-euc")



### Within-species 
### Large-Large
eucs_large <- eucs %>% 
  filter(size_to == "large") %>% 
  filter(size_from == "large") %>% 
  filter(species_from == species_to) %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region), NA))

large <- ggplot(eucs_large, 
                aes(x = alpha, 
                    y = species_from, 
                    colour = as.character(region),
                    fill = fill_col
                    )) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") + 
  geom_point(shape = 21,
             size = 4,   
             stroke = 1) + 
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#1DC199", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                      ) + 
  scale_fill_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                    na.value = "white", name = "Region") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Within-Species Large-Large") +
  theme_bw() + 
  xlim(-2, 1)
#+  guides(colour = guide_legend(override.aes = 
#                                 list(fill = colours, 
 ##                                     size = 2) ) ) 

#Facetted scatterplot by region 
#4 regions - TAS, VIC/NSW, NSW/QLD, WA

#We also want to get rid of species that are only present in one plot for now 
eucs_large <- eucs_large %>% 
  filter(!(species_from %in% c("Eucalyptus andrewsii", 
                             "Eucalyptus coccifera", 
                             "Eucalyptus dalrympleana", 
                             "Eucalyptus ovata", 
                             "Eucalyptus sieberi", 
                             "Eucalyptus subcrenulata", 
                             "Eucalyptus urnigera", 
                             "Eucalyptus viminalis", 
                             "Eucalyptus guilfoylei", 
                             "Eucalyptus radiata")))


eucs_large$region2 <- factor(eucs_large$region2, c("TAS", "SE_AUS", "N_AUS", "WA"))

ggplot(data = eucs_large) + 
  geom_point(aes(x = alpha, 
                y = species_from, 
                colour = region, 
                fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~region2, drop = TRUE) + 
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#1DC199", "#6FB1E7", "#BB9DEA", "#E48FD1"), ) + 
  scale_fill_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#6FB1E7", "#BB9DEA", "#E48FD1"),
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Large Within-Species Coefficient") + 
  ylab("Species")




### Small-Small 
eucs_small <- eucs %>% 
  filter(size_to == "small") %>% 
  filter(size_from == "small") %>% 
  filter(species_from == species_to) %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region), NA))

eucs_small <- eucs_small %>% 
  filter(!(species_from %in% c("Eucalyptus andrewsii", 
                               "Eucalyptus coccifera", 
                               "Eucalyptus dalrympleana", 
                               "Eucalyptus ovata", 
                               "Eucalyptus sieberi", 
                               "Eucalyptus subcrenulata", 
                               "Eucalyptus urnigera", 
                               "Eucalyptus viminalis", 
                               "Eucalyptus guilfoylei", 
                               "Eucalyptus radiata", 
                               "Eucalyptus jacksonii")))

eucs_small$region2 <- factor(eucs_small$region2, c("TAS", "SE_AUS", "N_AUS", "WA"))

small <- ggplot(eucs_small, 
                aes(x = alpha, 
                    y = species_from, 
                    colour = as.character(region),
                    fill = fill_col
                )) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") + 
  geom_point(shape = 21,
             size = 4, 
             stroke = 1) + 
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
  ) + 
  scale_fill_manual(values = c("#ED90A4", "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#E48FD1"), 
                    na.value = "white", name = "Region") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Within-Species Small-Small") +
  theme_bw() 
  
  
small + large
  


ggplot(data = eucs_small) + 
  geom_point(aes(x = alpha, 
                 y = species_from, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~region2, drop = TRUE) + 
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1")) + 
  scale_fill_manual(values = c("#ED90A4", "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#E48FD1"), 
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Small Within-Species Coefficient") + 
  ylab("Species")


### Facet all within-species interactions for diameter 
eucs <- final_sizemod %>% 
  filter(! species_from == "Non-euc") %>% 
  filter(! species_to == "Non-euc")

#eucs-within
eucs_within <- eucs %>% 
  filter(species_from == species_to)

eucs_within <- eucs_within %>% 
  mutate(size_int = paste0(size_from, sep = "_", size_to)) %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region), NA))


eucs_within <- eucs_within %>% 
  filter(! (species_to %in% c("Eucalyptus urnigera", 
                              "Eucalyptus sieberi", 
                              "Eucalyptus subcrenulata", 
                              "Eucalyptus andrewsii", 
                              "Eucalyptus coccifera", 
                              "Eucalyptus dalrympleana")))

ggplot(data = eucs_within) + 
  geom_point(aes(x = alpha, 
                 y = species_from, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_grid(size_int ~ region2) + 
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#1DC199", "#6FB1E7", "#BB9DEA", "#E48FD1")) + 
  scale_fill_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#1DC199", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Within-Species Coefficient") + 
  ylab("Species") + 
  xlim(-2, 1)





### Between-species
## Okay, now we need to think about between species 
eucs_bw <- eucs %>% 
  filter(!(species_from == species_to))

eucs_bw <- eucs_bw %>% 
  mutate(size_int = paste0(size_from, sep = "_", size_to)) %>% 
  mutate(species_int = paste0(species_from, sep = "_", species_to)) %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region), NA))

#18 large-large interactions, 40 small-large interactions, 23 small-small interactions 
eucs_bw %>% filter(size_int == "small_small") %>% 
  count(species_int)

#let's do each of the 4 general regions separately 
## North Aus 
n_aus <- eucs_bw %>% 
  filter(region2 == "N_AUS") 

naus <- ggplot(data = n_aus) + 
  geom_point(aes(x = alpha, 
                 y = species_int, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~size_int) + 
  scale_colour_manual(values = c("#B7AE50",  "#1DC199")) + 
  scale_fill_manual(values = c( "#B7AE50",  "#1DC199"), 
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Between-Species Coefficient") + 
  ylab("Species") 
naus


## WA  
wa <- eucs_bw %>% 
  filter(region2 == "WA") 

west_aus <- ggplot(data = wa) + 
  geom_point(aes(x = alpha, 
                 y = species_int, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~size_int) + 
  scale_colour_manual(values = c("#E48FD1")) + 
  scale_fill_manual(values = c( "#E48FD1"), 
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Between-Species Coefficient") + 
  ylab("Species") 
west_aus


## SE_AUS 
se_aus <- eucs_bw %>% 
  filter(region2 == "SE_AUS") 

se_aus <- ggplot(data = se_aus) + 
  geom_point(aes(x = alpha, 
                 y = species_int, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~size_int) + 
  scale_colour_manual(values = c( "#6FB1E7", "#BB9DEA")) + 
  scale_fill_manual(values = c( "#6FB1E7", "#BB9DEA"), 
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Between-Species Coefficient") + 
  ylab("Species") 
se_aus



## TAS 
tas <- eucs_bw %>% 
  filter(region2 == "TAS") 

tass <- ggplot(data = tas) + 
  geom_point(aes(x = alpha, 
                 y = species_int, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~size_int) + 
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#7EBA68", "#00BFC8")) + 
  scale_fill_manual(values = c("#ED90A4",  "#00BFC8"), 
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Between-Species Coefficient") + 
  ylab("Species") 
tass


##Huh, let's get rid of species and just look at the overall trends 
ggplot(data = eucs_bw) + 
  geom_point(aes(x = alpha, 
                 y = region2, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~size_int) + 
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#E48FD1")) + 
  scale_fill_manual(values = c("#ED90A4", "#B7AE50","#1DC199", "#00BFC8", "#6FB1E7", "#E48FD1"),
                    na.value = "white") + 
  theme_bw() + 
  xlab("Between-Species Coefficient") + 
  ylab("Species") 



#Summary Statistics 
d <- final_sizemod %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(!species_to == "Non-euc") %>% 
  filter(species_to == species_from)




d <- d %>% filter(size_to == "small") %>% 
  filter(size_from == "small")

sum <- d %>% 
  group_by(region, species_from) %>% 
  summarise(
    count = n(), 
    median_a = median(alpha), 
    mean_a = mean(alpha), 
    min_a = min(alpha), 
    max_a = max(alpha)
  )










### Model stats 
d <- final_sizemod %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(!species_to == "Non-euc")

d %>% count(sig)
d %>% count(region)

strong <- final_sizemod %>% 
  filter(alpha > 0.49 | alpha < -0.49)



###########
# Prediction 
auc_size <- data.frame(matrix(ncol = 4))
colnames(auc_size) <- c("site", "types", "auc_score", "n_points")

data_cleaned <- data_cleaned %>% 
  mutate(size = if_else(Diameter < 50, "small", "large")) %>% 
  mutate(group = paste(size, sep = " ", Genus_Species))

auc_working <- auc_size

auc_sizef <- function(site, 
                          threshold = 10, #threshold to exclude euc groups
                          short_range = 8){ #short range value 
  
  d <- data_cleaned %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
    mutate(new_group = paste(size, sep = " ", species)) %>% 
    group_by(new_group) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!observation_count < threshold) 
  
  counts <- d %>% count(new_group)
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
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_group) #make config 
  
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
  
  model_types <- levels(configuration$types)
  
  for(k in model_types){
    
    auc_working[k, "site"] <- site
    auc_working[k, "types"] <- k
    
    p <- plot_papangelou(fit, #to find the AUC of the train data
                         window = window, 
                         drop_type_from_configuration = TRUE, 
                         type = k, 
                         use_log = FALSE, 
                         show = k, 
                         return_papangelou = TRUE)
    
    p$v[is.na(p$v)] <- mean(p)
    p$v[p$v == -Inf] <- -1e10
    
    
    X <- subset(as.ppp(configuration, W = window), marks == k)
    if(X$n <= 0){
      auc_working[k,"auc_score"] <- NA
    } else {
      auc_working[k, "auc_score"] <- auc(X = X, covariate = as.function(p))
    }
    
    auc_working[k, "n_points"] <- X$n
  }
  
  return(auc_working)
  
}


#the sites 
sites <- unique(data_cleaned$Site_Name)
sites <- sites[!sites == "ANU363"] #only  has one group, large eucalypts
#run 
for(i in sites){ 
  
  auc_working <- auc_sizef(site = i, 
                               threshold = 10, 
                               short_range = 8)
  
  auc_size <- rbind(auc_size, auc_working) #bind into larger df 
  
}

auc_size <- data.frame(auc_size, row.names = NULL)
auc_size <- auc_size %>% 
  filter(!is.na(site))

auc_size %>% count(site)
#write.csv(auc_sizesite, "auc_site_size.csv")

#median AUC 
auc_size %>% count(types)

median(auc_size$auc_score)

a <- auc_size %>% filter(!grepl("Non-euc", types))
#OR 
a <- auc_size %>% filter(! types %in% c("large Non-euc", "small Non-euc"))

median(a$auc_score)




########### Cross validation
#As we know, the AUC method used above doesn't allow for any within-group interactions to be evaluated 
#Therefore we can use leave-n-individuals out cross validation

size_sp <- function(site, 
                      threshold = 10){ #short range value 
  d <- data_cleaned %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
    mutate(new_group = paste(size, sep = " ", species)) %>% 
    group_by(new_group) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!observation_count < threshold) 
  
  counts <- d %>% count(new_group)
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
  
}


## For all sites: 
sites <-unique(ausplot_data$Site_Name)
#sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]



     
        
  
       
 
       

  

full_auc_ind <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(full_auc_ind) <- c("iter", "types", "value", "split", "site")


a <- ind_validation1(site = "BenRidge", 
                     n_ind = 10, 
                     n_it = 20)

a$site <- "BenRidge"  #set column site to k 

full_auc_ind <- rbind(full_auc_ind, a) #bind into larger df 

full_auc_ind %>% count(site)


write.csv(full_auc_ind, "auc_inds_size_species.csv")





d <- full_auc_ind %>% 
  filter(!types %in% c("large Non-euc", "small Non-euc"))

auc_inds_sum <- full_auc_ind %>% 
  group_by(split) %>% 
  summarise(
    median_a = median(value), 
    mean_a = mean(value), 
    min_a = min(value), 
    max_a = max(value)
  )
