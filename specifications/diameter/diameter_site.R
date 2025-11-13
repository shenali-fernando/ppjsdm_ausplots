library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(patchwork)
library(scales)

##################
### Small v Large Diameter Trees based on site only

#Based on analyses of diameter and where it splits the crown classes, choosing 50cm at the difference between large and small trees 

#Load cleaned dataset
ausplot_data<- read.csv("data/data_cleaned.csv")

#Source other functions 
source("code/make_summary_fun.R")

data <- data_cleaned %>% 
  mutate(size = if_else(Diameter < 50, "small", "large")) %>% 
  mutate(group = paste(size, sep = " ", Genus_Species))


### Okay, now run the models - again do this is a long loop 
size_site <- function(site, 
                       threshold = 10, #threshold to exclude groups
                       short_range = 8){ #short range value 
  
  d <- data %>% filter(Site_Name == site)
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", "Eucalypt")) %>% 
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

full_site_size <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_site_size) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                            "hi_numerical", "Potential", "site")

for(i in sites){ 
  
  size_site_mod <- size_site(site = i, 
                         threshold = 10, #threshold to exclude euc groups
                         short_range = 8) #run site_ccmodel function
  
  working_df <- make_sum_df(fits = list(size_site_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(size_site_mod$sum))
  
  working_df <- working_df %>% mutate(site = i) #set new column for site name 
  
  full_site_size <- rbind(full_site_size, working_df) #bind into larger df 
  
}


#Check unreasonable CIs
full_site_size2 <- full_site_size %>% 
  mutate(range_ci = hi - lo)

large_ci <- full_site_size2 %>% 
  filter(range_ci > 5) #none, yay! 


### Adding columns 
#Need to add a sig column 
full_site_size2 <- full_site_size2 %>% 
  mutate(sig = ifelse((lo > 0 | hi < 0), 1, NA ))

#make columns easy for plotting
full_site_size2 <- full_site_size2 %>% 
  relocate(site) %>% 
  mutate(species_from = gsub("^[^ ]+ ", "\\1", from)) %>% 
  mutate(species_to = gsub("^[^ ]+ ", "\\1", to)) %>% #same for euc or non-euc 
  mutate(size_from = gsub("^([^ ]+).*", "\\1", from)) %>% 
  mutate(size_to = gsub("^([^ ]+).*", "\\1", to))



full_site_size2 <- full_site_size2 %>% 
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

full_site_size2 <- full_site_size2 %>% 
  mutate(region2 = region) %>% 
  mutate(region2 = if_else(region2 %in% c("d_TAS", "o_TAS"), 
                           "TAS", region2)) %>% 
  mutate(region2 = if_else(region2 %in% c("N_NSW", "QLD"), 
                           "N_AUS", region2)) %>%
  mutate(region2 = if_else(region2 %in% c("SE_NSW", "SE_VIC"), 
                           "SE_AUS", region2)) 


#save it out 
write.csv(full_site_size2, "final_site_size.csv")


##############################################################################################
############## VISUALISATION

### EUCS TO EUCS 
eucs <- full_site_size2 %>% 
  filter(! species_from == "Non-euc") %>% 
  filter(! species_to == "Non-euc")

eucs <- eucs %>% 
  mutate(class_int = paste0(size_from, sep = "_", size_to)) %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region), NA))




ggplot(data = eucs) + 
  geom_point(aes(x = alpha, 
                 y = class_int, 
                 colour = region, 
                 fill = fill_col), 
             size = 3, 
             shape = 21, 
             stroke = 2, 
             position = position_jitter(width = 0.25, height = 0.1)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  facet_wrap(~region2, ncol = 4) +
  scale_colour_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#1DC199", "#6FB1E7", "#BB9DEA", "#E48FD1"), ) + 
  scale_fill_manual(values = c("#ED90A4", "#DC9E70", "#B7AE50", "#1DC199", "#6FB1E7", "#BB9DEA", "#E48FD1"),
                    na.value = "white", guide = "none") + 
  theme_bw() + 
  xlab("Within-Species Coefficient") + 
  ylab("Species")



### Summary Stats
d <- full_site_size2 %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(!species_to == "Non-euc")


d <- d %>% filter(size_to == "large") %>% 
  filter(size_from == "large")

sum <- d %>% 
  group_by(region) %>% 
  summarise(
    count = n(), 
    median_a = median(alpha), 
    mean_a = mean(alpha), 
    min_a = min(alpha), 
    max_a = max(alpha)
  )




### Model stats 
d <- full_site_size2 %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(!species_to == "Non-euc")

full_site_size2 %>% count(sig)
full_site_size2 %>% count(region)

strong <- d %>% 
  filter(alpha > 0.49 | alpha < -0.49)



###########
# Prediction 
auc_sizesite <- data.frame(matrix(ncol = 4))
colnames(auc_sizesite) <- c("site", "types", "auc_score", "n_points")

auc_working <- auc_sizesite

auc_site_size <- function(site, 
                          threshold = 10, #threshold to exclude euc groups
                          short_range = 8){ #short range value 
  
  d <- data %>% filter(Site_Name == site)
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", "Eucalypt")) %>% 
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
  
  auc_working <- auc_site_size(site = i, 
                               threshold = 10, 
                               short_range = 8)
  
  auc_sizesite <- rbind(auc_sizesite, auc_working) #bind into larger df 
  
}

auc_sizesite <- data.frame(auc_sizesite, row.names = NULL)
auc_sizesite <- auc_sizesite %>% 
  filter(!is.na(site))

auc_sizesite %>% count(site)
write.csv(auc_sizesite, "auc_site_size.csv")

#median AUC 
auc_sizesite %>% count(types)

median(auc_sizesite$auc_score)

a <- auc_sizesite %>% filter(types %in% c("large Eucalypt", "small Eucalypt"))

median(a$auc_score)




########### Cross validation
#As we know, the AUC method used above doesn't allow for any within-group interactions to be evaluated 
#Therefore we can use leave-n-individuals out cross validation

site_size <- function(site, 
                      threshold = 10){ #short range value 
  d <- data_cleaned %>% 
    mutate(size = if_else(Diameter < 50, "small", "large")) %>% 
    mutate(group = paste(size, sep = " ", Genus_Species)) %>% 
    filter(Site_Name == site)
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", "Eucalypt")) %>% 
    mutate(new_group = paste0(size, sep = " ", species)) %>% 
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


#Nice, it's working 

## For all sites: 
a <- unique(ausplot_data$Site_Name) 
#sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]

                                                           

full_auc_ind <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(full_auc_ind) <- c("iter", "types", "value", "split", "site")


a <- ind_validation1(site = "ANU363", 
                     n_ind = 10, 
                     n_it = 20)
  
  a$site <- "MtMaurice" #set column site to k 
  
  full_auc_ind <- rbind(full_auc_ind, a) #bind into larger df 
  
  full_auc_ind %>% count(site)

write.csv(full_auc_ind, "auc_inds_size_site.csv")

d <- full_auc_ind %>% 
  filter(!types %in% c("large Non-euc", "small Non-euc"))

auc_inds_sum <- d %>% 
  group_by(split) %>% 
  summarise(
    median_a = median(value), 
    mean_a = mean(value), 
    min_a = min(value), 
    max_a = max(value)
  )

