library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)


#Load cleaned data 
data_cleaned <- read.csv("data/data_cleaned.csv")

#Source functions 
source(code/make_summary_fun.R)


# Make a function to create a fit and summary for each site 
spcc_model <- function(site, 
                       threshold_noneuc = 10, #threshold to exclude non-euc groups
                       threshold_euc = 10, #threshold to exclude euc groups
                       short_range = 8){ #short range value 
  
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
  
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_spcc) #make config 
  
  # Set parameters
  window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
  nspecies <- length(levels(configuration$types))
  short_range <- matrix(short_range, nspecies, nspecies)
  model <- "exponential"
  
  fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                       window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                       short_range = matrix(8, nspecies, nspecies), 
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
sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]

full_df_spcc <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df_spcc) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                           "hi_numerical", "Potential", "site")

for(i in sites){ 
  
  spcc_mod <- spcc_model(site = i, 
                         threshold_noneuc = 10, #threshold to exclude non-euc groups
                         threshold_euc = 10, #threshold to exclude euc groups
                         short_range = 8) #run site_ccmodel function
  
  working_df <- make_sum_df(fits = list(spcc_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(spcc_mod$sum))
  
  working_df <- working_df %>% mutate(site = i) #set new column for site name 
  
  full_df_spcc <- rbind(full_df_spcc, working_df) #bind into larger df 
  
}

#check NA's first 
is.na(full_df_spcc$alpha)


##Okay need to check unreasonable CIs before doing anything 
full_df_spcc <- full_df_spcc %>% 
  mutate(range_ci = hi - lo)

large_ci <- full_df_spcc %>% 
  filter(range_ci > 5)

#Hm, so we have 9 interactions from 7 sites that have an interaction CI range larger than abs(5)
#We need to manually go through the 7 sites and have a look at what's happening 

#The worst offender - MTFIELD
mtfield <- data_cleaned %>% filter(Site_Name == "MtField") %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
  mutate(new_spcc = paste(species, Crown_Class))  %>% 
  group_by(new_spcc) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(str_starts(new_spcc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < 10)) %>% 
  filter(! new_spcc %in%  c("Eucalyptus coccifera Intermediate", "Eucalyptus subcrenulata Intermediate"))
mtfield %>% count(new_spcc)


mtfield <- jitter_points(mtfield)
configuration<- Configuration(mtfield$x_jitter, mtfield$y_jitter, types = mtfield$new_spcc)
plot(configuration)

nspecies <- length(levels(configuration$types))

fit_mtfield <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                     window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                     short_range = matrix(8, nspecies, nspecies), 
                     model = "exponential", 
                     saturation = 10, 
                     nthreads = 4, 
                     fitting_package = "glmnet",
                     dummy_distribution = "stratified",
                     min_dummy = 1, dummy_factor = 1e10, 
                     max_dummy = 1e3)

sum_mtfield <- summary(fit_mtfield)



### CAVESIDE
caveside <- data_cleaned %>% filter(Site_Name == "Caveside") %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
  mutate(new_spcc = paste(species, Crown_Class))  %>% 
  group_by(new_spcc) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(str_starts(new_spcc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < 10)) %>% 
  filter(!(new_spcc == "Eucalyptus obliqua Co/dominant"))
caveside %>% count(new_spcc)


caveside <- jitter_points(caveside)
configuration<- Configuration(caveside$x_jitter, caveside$y_jitter, types = caveside$new_spcc)
plot(configuration)

nspecies <- length(levels(configuration$types))

fit_caveside <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                              window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                              short_range = matrix(8, nspecies, nspecies), 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)

sum_caveside <- summary(fit_caveside)



### COLLINS
collins <- data_cleaned %>% filter(Site_Name == "Collins") %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
  mutate(new_spcc = paste(species, Crown_Class))  %>% 
  group_by(new_spcc) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(str_starts(new_spcc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < 10)) %>% 
filter(!new_spcc == "Eucalyptus diversicolor Intermediate")
collins %>% count(new_spcc)


collins <- jitter_points(collins)
configuration<- Configuration(collins$x_jitter, collins$y_jitter, types = collins$new_spcc)
plot(configuration)

nspecies <- length(levels(configuration$types))

fit_collins <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                               window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                               short_range = matrix(8, nspecies, nspecies), 
                               model = "exponential", 
                               saturation = 10, 
                               nthreads = 4, 
                               fitting_package = "glmnet",
                               dummy_distribution = "stratified",
                               min_dummy = 1, dummy_factor = 1e10, 
                               max_dummy = 1e3)

sum_collins <- summary(fit_collins)


### BRUXNER 
brux <- data_cleaned %>% filter(Site_Name == "Bruxner") %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
  mutate(new_spcc = paste(species, Crown_Class))  %>% 
  group_by(new_spcc) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(str_starts(new_spcc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < 10)) %>% 
  filter(! new_spcc %in% c("Eucalyptus microcorys Suppressed", "Eucalyptus pilularis Co/dominant"))
brux %>% count(new_spcc)

brux <- jitter_points(brux)
configuration<- Configuration(brux$x_jitter, brux$y_jitter, types = brux$new_spcc)
plot(configuration)

nspecies <- length(levels(configuration$types))

fit_brux <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                              window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                              short_range = matrix(8, nspecies, nspecies), 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)

sum_brux <- summary(fit_brux)

ppjsdm::box_plot("Grandis only" = fit_brux2, "All" = fit_brux, 
                 summ = list(sum_brux2, sum_brux), 
                 coefficient = "alpha", 
                 which = "all")

#Getting rid of both Eucalyptus microcorys Suppressed and Eucalyptus pilularis Co/dominant results in a better fit 



### A-TREE 
atree <- data_cleaned %>% filter(Site_Name == "A-Tree") %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
  mutate(new_spcc = paste(species, Crown_Class))  %>% 
  group_by(new_spcc) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(str_starts(new_spcc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < 10)) %>% 
  filter(! new_spcc == "Eucalyptus microcorys Intermediate")
atree %>% count(new_spcc)

atree <- jitter_points(atree)
configuration<- Configuration(atree$x_jitter, atree$y_jitter, types = atree$new_spcc)
plot(configuration)

nspecies <- length(levels(configuration$types))

fit_atree <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                              window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                              short_range = matrix(8, nspecies, nspecies), 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)

sum_atree <- summary(fit_atree)

### WOGWAY
wogway <- data_cleaned %>% filter(Site_Name == "WogWay") %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
  mutate(new_spcc = paste(species, Crown_Class))  %>% 
  group_by(new_spcc) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(str_starts(new_spcc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < 10)) %>% 
  filter(! new_spcc == "Eucalyptus fastigata Co/dominant")
wogway %>% count(new_spcc)

wogway <- jitter_points(wogway)
configuration<- Configuration(wogway$x_jitter, wogway$y_jitter, types = wogway$new_spcc)
plot(configuration)

nspecies <- length(levels(configuration$types))

fit_wogway <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                            window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                            short_range = matrix(8, nspecies, nspecies), 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)

sum_wogway <- summary(fit_wogway)

### WARATAHMIX
waratah <- data_cleaned %>% filter(Site_Name == "WaratahMix") %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
  mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
  mutate(new_spcc = paste(species, Crown_Class))  %>% 
  group_by(new_spcc) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!(str_starts(new_spcc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < 10))
waratah %>% count(new_spcc)

waratah <- jitter_points(waratah)
configuration<- Configuration(waratah$x_jitter, waratah$y_jitter, types = waratah$new_spcc)
plot(configuration)

nspecies <- length(levels(configuration$types))

fit_waratah <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                            window = ppjsdm::Rectangle_window(c(0, 100), c(0,100)),
                            short_range = matrix(8, nspecies, nspecies), 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)

sum_waratah <- summary(fit_waratah)

#Althought Waratah mix has a large CI for Eucalyptus cypellocarpa Co/dominant, it is significant 
#In addition, this group has other significant interactions with other groups 
#I will leave it as it, although the coefficient is large, I want to keep its interactions 




############ Fix dataframe 

## Add cleaned sites together into a df 
df <- make_sum_df(fits = list(Caveside = fit_caveside, 
                              Collins = fit_collins, 
                              Bruxner = fit_brux, 
                              MtField = fit_mtfield, 
                              'A-Tree' = fit_atree, 
                              WogWay = fit_wogway), 
                  coefficient = "alpha",
                  summ = list(sum_caveside, 
                              sum_collins, 
                              sum_brux, 
                              sum_mtfield, 
                              sum_atree, 
                              sum_wogway))

#add the range_ci to make number of columns equal in each df 
df <- df %>% 
  mutate(range_ci = hi - lo)

large_ci <- df %>% 
  filter(range_ci > 5) #solved! 

#add some plotting columns 
df <- df %>% 
  mutate(species_from = gsub("\\s+\\S+$", "", from)) %>% 
  mutate(species_to = gsub("\\s+\\S+$", "", to)) %>% 
  mutate(class_from = gsub(".*\\s+", "", from)) %>% 
  mutate(class_to = gsub(".*\\s+", "", to))

#rename col
df <- df %>% 
  rename(site = Fit)

#count sites to make sure 
df %>% count(site)

#I want to look at how many individuals are in each class 
inds_full <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(inds_full) <- c("new_spcc", "n", "site_name")

for(i in sites){ #do all sites at once and send into df 
  
  inds <- spcc_ind(site = i,  
                   threshold_noneuc = 10, 
                   threshold_euc = 10) #function above
  
  inds <- inds %>% mutate(site_name = i) #set new column called site_name for site 
  
  inds_full <- rbind(inds_full, inds) #bind into larger df 
}

full_df_spcc <- full_df_spcc %>% 
  mutate(species_from = gsub("\\s+\\S+$", "", from)) %>% 
  mutate(species_to = gsub("\\s+\\S+$", "", to)) %>% 
  mutate(class_from = gsub(".*\\s+", "", from)) %>% 
  mutate(class_to = gsub(".*\\s+", "", to))

inds_full <- inds_full %>% 
  rename(from = to) 

full_df_spcc <- full_df_spcc %>%
  left_join(inds_full, by = c("site", "from")) 

d <- full_df_spcc %>% filter(ind_from %in% c("10", "11") | ind_to %in% c("10", "11"))

#Cool, just checking if I need to increase thresholds but nah 
#So get rid of my inds cols 
full_df_spcc <- full_df_spcc %>% 
  dplyr::select(-ind_from)

### Delete the sites manually cleaned 
full_df_spcc <- full_df_spcc %>% 
  filter(! site %in% c("Caveside", "Bruxner", "Collins", "MtField", "A-Tree", "WogWay"))
full_df_spcc %>% count(site)

#Cool, now bind manually done sites to full df
final_spcc_mod <- rbind(full_df_spcc, df)

final_spcc_mod %>% count(site)


## Add in significance column 
final_spcc_mod <- final_spcc_mod %>% mutate(sig = ifelse((lo > 0 | hi < 0), 1, NA))

### Write in region and region2 
final_spcc_mod <- final_spcc_mod %>% 
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
  mutate(region = if_else(region %in% c("BondTier", "BlackRiver", "Weld", "Weld", "MtField", "ZigZag", "Supersite"), 
                          "o_TAS", region)) 


final_spcc_mod <- final_spcc_mod%>% 
  mutate(region2 = region) %>% 
  mutate(region2 = if_else(region2 %in% c("Ben_TAS", "King_TAS", "S_TAS", "N_TAS"), 
                           "TAS", region2)) %>% 
  mutate(region2 = if_else(region2 %in% c("N_NSW", "QLD"), 
                           "N_AUS", region2)) %>%
  mutate(region2 = if_else(region2 %in% c("SE_NSW", "SE_VIC"), 
                           "SE_AUS", region2)) 

#save out csv 
write.csv(final_spcc_mod, "final_spcc_mod.csv")





#######################################################################################################
##################################### VISUALISATION ###################################################

## The goal is to understand what species-size-site patterns exist
#i.e. if there is a change in the amount of repulsion and attraction between species sizes across sites which span a geographic region of Australia 

######################## BETWEEN EUCS ONLY 

eucs <- final_spcc_mod %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(! species_to == "Non-euc")

#Rename syncarpia glomulifera 
eucs <- eucs %>% 
  mutate(species_from = if_else(species_from %in% c("Syncarpia glomulifera", 
                                                    "Syncarpia glomulifera subsp glabra", 
                                                    "Syncarpia glomulifera subsp glomulifera"), 
                                "Syncarpia glomulifera", species_from)) %>% 
  mutate(species_to = if_else(species_to %in% c("Syncarpia glomulifera", 
                                                    "Syncarpia glomulifera subsp glabra", 
                                                    "Syncarpia glomulifera subsp glomulifera"), 
                                "Syncarpia glomulifera", species_to))

eucs %>% count(species_to, species_from)

## Add a class_int col 
eucs <- eucs %>% 
  mutate(class_int = paste0(class_from, sep = "_", class_to))

eucs %>% count(species_to, species_from, class_int)

#### Within-group eucalypts
eucs_within <- eucs %>%
  filter(species_from == species_to)
  
eucs_within %>% count(species_from, species_to)

#### Between-group eucalypts
eucs_bw <- eucs %>% 
  filter(! (species_from == species_to))

eucs_bw %>% count(species_from, species_to)


############# WITHIN EUCS PLOTS 

## Co/dominant 
eucs_within <- eucs_within %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region), NA))

eucs_w <- eucs_within %>% 
  filter(class_int %in% c("Co/dominant_Co/dominant", 
                          "Intermediate_Intermediate",
                          "Suppressed_Suppressed"))

ggplot(data = eucs_w) + 
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
             linetype = "dashed")  +
  facet_wrap(~class_int, drop = TRUE) + 
  scale_colour_manual(values = c("#ED90A4",  "#B7AE50", "#7EBA68",  "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1")) + 
  scale_fill_manual(values = c("#ED90A4",  "#B7AE50", "#7EBA68",  "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                    na.value = "white",
                    guide = "none") + 
  scale_y_discrete(limits = rev) + 
  theme_bw() + 
  xlab("Within Eucalypt Species Coefficient") + 
  ylab("Species") + 
  xlim(-3, 1)


ggplot(data = eucs_within) + 
  geom_boxplot(aes(x = alpha, y = species_from)) + 
  geom_point(aes(x = alpha, 
                 y = species_from, 
                 colour = region)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed")  +
  facet_wrap(~class_int, drop = TRUE) + 
  theme_bw() + 
  xlim(-2, 1)



cd <- df_spcc %>% 
  filter(class_from == "Co/dominant") %>% 
  filter(class_to == "Co/dominant") %>% 
  mutate(sig = ifelse((lo > 0 | hi < 0), "1", "0" )) 
#%>% 
#  filter(!species_from %in% c("Non-euc"))

scatterplot_cd <- ggplot() +
  geom_vline(xintercept = 0, colour = "red", linewidth = .5) +
  geom_beeswarm(data = cd, 
                aes(x = alpha, 
                    y = species_from, 
                    fill = sig, 
                    text = paste("alpha:", alpha, "<br>", 
                                 "species int:", species_from, "<br>")), 
                shape = 21, 
                size = 4,
                alpha = 0.8) +
  scale_fill_manual(values = c("white", "black")) +
  theme_bw() + 
  guides(colour = guide_legend(override.aes = 
                                 list(fill = colours, 
                                      size = 3) ) ) + 
  ylab("") + 
  ggtitle("C/D - C/D") +
  xlim(-5, 1)
scatterplot_cd

"Co/dominant_Co/dominant"
"Intermediate_Intermediate"
"Suppressed_Suppressed"

d <- eucs_w %>% filter(species_from == "Syncarpia glomulifera") %>% 
  filter(class_int == "Suppressed_Suppressed")
  
summary(d$alpha)
















## Suppressed 

sup <- df_spcc %>% 
  filter(class_from == "Suppressed") %>% 
  filter(class_to == "Suppressed") %>% 
  mutate(sig = ifelse((lo > 0 | hi < 0), "1", "0" )) 


scatterplot_sup <- ggplot() +
  geom_vline(xintercept = 0, colour = "red", linewidth = .5) +
  geom_beeswarm(data = sup, 
                aes(x = alpha, 
                    y = species_from, 
                    fill = sig, 
                    text = paste("alpha:", alpha, "<br>", 
                                 "species int:", species_from, "<br>")), 
                shape = 21, 
                size = 4,
                alpha = 0.8) +
  theme_bw() + 
  scale_fill_manual(values = c("white", "black")) +
  guides(colour = guide_legend(override.aes = 
                                 list(fill = colours, 
                                      size = 3) ) ) + 
  ylab("") + 
  ggtitle("S - S") + 
  xlim(-5, 1)
scatterplot_sup



## Suppressed 

int <- df_spcc %>% 
  filter(class_from == "Intermediate") %>% 
  filter(class_to == "Intermediate") %>% 
  mutate(sig = ifelse((lo > 0 | hi < 0), "1", "0" )) 


scatterplot_int <- ggplot() +
  geom_vline(xintercept = 0, colour = "red", linewidth = .5) +
  geom_beeswarm(data = int, 
                aes(x = alpha, 
                    y = species_from, 
                    fill = sig, 
                    text = paste("alpha:", alpha, "<br>", 
                                 "species int:", species_from, "<br>")), 
                shape = 21, 
                size = 4,
                alpha = 0.8) +
  theme_bw() + 
  scale_fill_manual(values = c("white", "black")) +
  guides(colour = guide_legend(override.aes = 
                                 list(fill = colours, 
                                      size = 3) ) ) + 
  ylab("") + 
  ggtitle("I - I") + 
  xlim(-5, 1)
scatterplot_int


scatterplot_cd + scatterplot_int + scatterplot_sup & plot_layout(nrow = 1)



within_ints <- df_spcc %>%  
  filter(class_from == class_to) %>% 
  mutate(sig = ifelse((lo > 0 | hi < 0), "1", "0" ))  %>% 
  mutate(int_type = ifelse((class_from == "Co/domiannt" & class_to == "Co/dominant"), "Co/dominant", class_from))


ggplot() +
  geom_vline(xintercept = 0, colour = "red", linewidth = .5) +
  geom_beeswarm(data = within_ints, 
                aes(x = alpha, 
                    y = species_from, 
                    fill = sig, 
                    text = paste("alpha:", alpha, "<br>", 
                                 "species int:", species_from, "<br>")), 
                shape = 21, 
                size = 4,
                alpha = 0.8, 
                position = position_dodge2(width = 0.75)) +
  scale_fill_manual(values = c("white", "black")) +
  theme_bw() + 
  facet_wrap(~class_from) +
  guides(colour = guide_legend(override.aes = 
                                 list(fill = colours, 
                                      size = 3) ) ) + 
  ylab("Species") + 
  ggtitle("Within-group associations") + 
  xlim(-5, 1)
  





###########################################################################################

### Prediction 

# Need to understand performance of the model 

#make empty df 
auc_df_spcc <- data.frame(matrix(ncol = 4))
colnames(auc_df_spcc) <- c("site", "types", "auc_score", "n_points")

auc_working <- auc_df_spcc

auc_site_spcc <- function(site, 
                                   threshold_noneuc = 10, #threshold to exclude non-euc groups
                                   threshold_euc = 10, #threshold to exclude euc groups
                                   short_range = 8){ #short range value 
  
  d <- data_cleaned %>% filter(Site_Name == site)
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia")), "Non-euc", Genus_Species)) %>% 
    mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
    mutate(new_spcc = paste(species, Crown_Class))  
  
  
  d <- d %>%  #thresholds for groups 
    group_by(new_spcc) %>% 
    mutate(observation_count = n()) %>% 
    filter(!(str_starts(new_spcc, "Non-euc") & observation_count < threshold_noneuc)) %>% 
    filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia")) & observation_count < threshold_euc)) 
  
  
  counts <- d %>% count(new_spcc)
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
  
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$new_spcc) #make config 
  
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
sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]

#run 
for(i in sites){ 
  
  auc_working <- auc_site_spcc(site = i, 
                          threshold_noneuc = 8, 
                          threshold_euc = 10, 
                          short_range = 8)
  
  auc_df_spcc <- rbind(auc_df_spcc, auc_working) #bind into larger df 
  
}


auc_df_spcc2 <- data.frame(auc_df_spcc, row.names = NULL)
auc_df_spcc2 <- auc_df_spcc2 %>% 
  filter(!is.na(site))

auc_df_spcc2 %>% count(site)

#only care about euc scores so 
auc_df_spcc2 <- auc_df_spcc2 %>% 
  filter(str_starts(types, "E|C"))

auc_df_spcc2 %>% count(types)

auc_df_spcc2 <- auc_df_spcc2 %>% 
  group_by(types) %>% 
  mutate(types_count = n())
  
auc_df_spcc2 <- unique(auc_df_spcc2)

auc_df_spcc2 %>% count(types_count)

# Make another column 
auc_df_spcc2 <- auc_df_spcc2 %>% 
  mutate(class = gsub(".*\\s+", "", types))


#Visualisation 
spcc <- ggplot(data = auc_df_spcc2, 
       aes(x = class, 
           y = auc_score)) + 
  geom_boxplot() + 
  geom_point() + 
  ylim(0.4, 0.9) + 
  ggtitle("SpCC")


library(patchwork)
spcc + cc








### Model stats 
d <- final_spcc_mod %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(!species_to == "Non-euc")

d %>% count(sig)
d <- d%>% count(species_from)

strong <- final_spcc_mod %>% 
  filter(alpha > 0.49 | alpha < -0.49)




#median AUC 
auc_df_spcc2 %>% count(types)

median(auc_df_spcc2$auc_score)

a <- auc_df_spcc2 %>% filter(types %in% c("Non-euc Co/dominant", "Non-euc Intermediate", "Non-euc Suppressed"))

median(a$auc_score)




########### Cross validation
#As we know, the AUC method used above doesn't allow for any within-group interactions to be evaluated 
#Therefore we can use leave-n-individuals out cross validation

spcc <- function(site, 
                      threshold = 10){ #short range value 
  d <- data_cleaned %>% filter(Site_Name == site)
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia")), "Non-euc", Genus_Species)) %>% 
    mutate(Crown_Class = if_else(Crown_Class %in% c("Dominant", "Emergent", "Co-dominant"), "Co/dominant", Crown_Class)) %>% 
    mutate(new_spcc = paste(species, Crown_Class))  
  
  
  d <- d %>%  #thresholds for groups 
    group_by(new_spcc) %>% 
    mutate(observation_count = n()) %>% 
    filter(!(str_starts(new_spcc, "Non-euc") & observation_count < threshold)) %>% 
    filter(!(str_starts(new_spcc,  regex("^Eucalyptus|^Corymbia")) & observation_count < threshold)) 
  
  
  counts <- d %>% count(new_spcc)
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


#Nice, it's working 

## For all sites: 
sites <- unique(data_cleaned$Site_Name) 
sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]



full_auc_ind <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(full_auc_ind) <- c("iter", "types", "value", "split", "site")


a <- ind_validation3(site = "MtField" , 
                     n_ind = 10, 
                     n_it = 20)

a$site <- "MtField"  #set column site to k 

full_auc_ind <- rbind(full_auc_ind, a) #bind into larger df 


full_auc_ind %>% count(site)



write.csv(full_auc_ind, "auc_inds_spcc.csv")

d <- full_auc_ind %>% 
  filter(!types %in% c("Non-euc Co/dominant", "Non-euc Intermediate", "Non-euc Suppressed"))

auc_inds_sum <- full_auc_ind %>% 
  group_by(split) %>% 
  summarise(
    median_a = median(value), 
    mean_a = mean(value), 
    min_a = min(value), 
    max_a = max(value)
  )
