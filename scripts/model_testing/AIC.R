library(dplyr)
library(tidyr)
library(ggplot2)
library(spatstat)
library(ppjsdm)
library(ggpubr)
library(patchwork)

##### Shockingly, the AIC and BIC are computed in the fit for the fit!!!!!! 
## this would've been so useful... 

# Cool, i need to understand how well each of the 5 model specifications are doing compared to each other:
#Null, Sp, CC, Size, SpCC, SpSize

#make empty df 
aic_null <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_null) <- c("Site_Name", "AIC", "BIC", "model")

#Load the data 
data <- read.csv("data/data_cleaned.csv")

site_names <- unique(data$Site_Name)


#Null Model 

#get rid of duplicate points first 
data <- data %>%
  group_by(Site_Name, Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) %>%
  ungroup() 

window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

for(i in seq_along(site_names)){
  
  site <- site_names[i]
  
  df <- data %>% filter(Site_Name == site) #filter data to site
  
  configuration <- Configuration(df$x_jitter, df$y_jitter)
  
  fit <- ppjsdm::gibbsm(configuration = configuration, #fit model
                        window = window,
                        short_range = matrix(8, 1, 1), 
                        model = "exponential", 
                        saturation = 10, 
                        nthreads = 4, 
                        fitting_package = "glmnet",
                        dummy_distribution = "stratified",
                        min_dummy = 1, dummy_factor = 1e10, 
                        max_dummy = 1e3)
  
  aic_null$Site_Name[i] <- site
  
  aic_null$AIC[i] <- fit$aic 
  
  aic_null$BIC[i] <- fit$bic
}

aic_null$model <- "null"


#Crown Class Model 
aic_cc <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_cc) <- c("Site_Name", "AIC", "BIC", "model")

site_names <- site_names[!site_names %in% c("Dip", "Bird", "Flowerdale")]

for (i in seq_along(site_names)) {
  site <- site_names[i]
  
  d <- data %>% filter(Site_Name == site) %>% 
  mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
    mutate(new_cc = paste(species, Crown_Class))%>% 
    mutate(new_cc = if_else(new_cc == "Eucalyptus Emergent", "Eucalyptus Dominant", new_cc)) %>% 
    mutate(new_cc = if_else(new_cc %in% c("Eucalyptus Co-dominant", "Eucalyptus Dominant"), "Eucalyptus Co/dominant", new_cc)) %>% 
    mutate(new_cc = if_else(new_cc %in% c("Non-euc Co-dominant", "Non-euc Dominant"), "Non-euc Co/dominant", new_cc))
  
  d <- d %>%  #thresholds for groups 
    group_by(new_cc) %>% 
    mutate(observation_count = n()) %>% 
    filter(!observation_count < 8)
  
  
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
  
  fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
                       window = window,
                       short_range = matrix(8, nspecies, nspecies), 
                       model = "exponential", 
                       saturation = 10, 
                       nthreads = 4, 
                       fitting_package = "glmnet",
                       dummy_distribution = "stratified",
                       min_dummy = 1, dummy_factor = 1e10, 
                       max_dummy = 1e3)
  aic_cc$Site_Name[i] <- site
  
  aic_cc$AIC[i] <- fit$aic 
  
  aic_cc$BIC[i] <- fit$bic
}

aic_cc$model <- "cc"



#Size Model 
aic_di <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_di) <- c("Site_Name", "AIC", "BIC", "model")

site_names <- unique(data$Site_Name)

data <- data %>% 
  mutate(size = if_else(Diameter < 50, "small", "large")) %>% 
  mutate(group = paste(size, sep = " ", Genus_Species))

for (i in seq_along(site_names)) {
  site <- site_names[i]
  
d <- data %>% filter(Site_Name == site)

d <- d %>% 
  mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", "Eucalypt")) %>% 
  mutate(new_group = paste(size, sep = " ", species)) %>% 
  group_by(new_group) %>% 
  mutate(observation_count = n()) %>% 
  ungroup() %>% 
  filter(!observation_count < 10) 

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

aic_di$Site_Name[i] <- site

aic_di$AIC[i] <- fit$aic 

aic_di$BIC[i] <- fit$bic
}

aic_di$model <- "diam"




#Species by Crown Class Model 
aic_spcc <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_spcc) <- c("Site_Name", "AIC", "BIC", "model")

site_names <- site_names[!site_names %in% c("Dip", "Bird", "Flowerdale")]

for (i in seq_along(site_names)) {
  site <- site_names[i]
  
  d <- data %>% filter(Site_Name == site) %>% 
      mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia||^Syncarpia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
      mutate(new_cc = paste(Genus_Species, Crown_Class))%>% 
      mutate(new_cc = if_else(new_cc %in% c("Eucalyptus Co-dominant", "Eucalyptus Dominant", "Eucalyptus Emergent"), "Eucalyptus Co/dominant", new_cc)) %>% 
      mutate(new_cc = if_else(new_cc %in% c("Non-euc Co-dominant", "Non-euc Dominant"), "Non-euc Co/dominant", new_cc)) %>% 
      mutate(group = paste(Genus_Species, str_extract(new_cc, "\\S+$")))

  d <- d %>%  #thresholds for groups 
    group_by(group) %>% 
    mutate(observation_count = n()) %>% 
    filter(!observation_count < 10) 
  
  
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

  aic_spcc$Site_Name[i] <- site
  
  aic_spcc$AIC[i] <- fit$aic 
  
  aic_spcc$BIC[i] <- fit$bic
  }

aic_spcc$model <- "spcc"



#Species by Size Model 
aic_spdi <- as.data.frame(matrix(ncol = 4, nrow = 48))
colnames(aic_spdi) <- c("Site_Name", "AIC", "BIC", "model")

site_names <- unique(data$Site_Name)

data <- data %>% 
  mutate(size = if_else(Diameter < 50, "small", "large")) %>% 
  mutate(group = paste(size, sep = " ", Genus_Species))

for (i in seq_along(site_names)) {
  site <- site_names[i]
  d <- data %>% filter(Site_Name == site)
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
    mutate(new_group = paste(size, sep = " ", species)) %>% 
    group_by(new_group) %>% 
    mutate(observation_count = n()) %>% 
    ungroup() %>% 
    filter(!observation_count < 10) 
  
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
  aic_spdi$Site_Name[i] <- site
  
  aic_spdi$AIC[i] <- fit$aic 
  
  aic_spdi$BIC[i] <- fit$bic
}

aic_spdi$model <- "spdiam"


################################################################################
##### Join all aic values together into a long dataframe 

aic_mods <- rbind(aic_cc, aic_di, aic_spdi, aic_null, aic_spcc)
aic_mods <- na.omit(aic_mods)
aic_mods %>% count(model) #everything looks correct 
aic_mods %>% count(Site_Name) 

#save out 
write.csv(aic_mods, "aic_mods1.csv")

#Add region 
aic_mods <- aic_mods %>% 
mutate(region = Site_Name) %>% 
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



### Make a quick plot 
aic <- ggplot(data = aic_mods, 
       aes(x = model, 
           y = AIC)) + 
  #geom_line(aes(group = region)) + 
  geom_boxplot() + 
  geom_point(aes(colour = Site_Name)) + 
  theme_bw() + 
  ylab("") + 
  ggtitle("AIC") + 
  ylim(c( -18400, 15))

bic <- ggplot(data = aic_mods, 
              aes(x = model, 
                  y = BIC)) + 
  #geom_line(aes(group = region)) + 
  geom_boxplot() +
  geom_point(aes(colour = Site_Name)) + 
  theme_bw() +
  ylab("") + 
  rremove("y.text") + 
  rremove("y.ticks") +
  ggtitle("BIC") + 
  ylim(c( -18400, 15))


aic + bic + plot_layout(guides = "collect")
