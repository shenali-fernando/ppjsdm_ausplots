library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)


#Load cleaned data 
data_cleaned <- read.csv("data/data_cleaned.csv")


# Make a function to create a fit and summary for each site 
spcc_model <- function(site, 
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



full_df_spcc2 <- full_df_spcc %>% 
  mutate(species_from = gsub("\\s+\\S+$", "", from)) %>% 
  mutate(species_to = gsub("\\s+\\S+$", "", to)) %>% 
  mutate(class_from = gsub(".*\\s+", "", from)) %>% 
  mutate(class_to = gsub(".*\\s+", "", to))



##Okay need to check unreasonable CIs 
full_df_spcc2 <- full_df_spcc2 %>% 
  mutate(range_ci = hi - lo)

large_ci <- full_df_spcc2 %>% 
  filter(range_ci > 5)




find_ind_spcc <- function(site,  
                     threshold_noneuc = 8, 
                     threshold_euc = 8){
  
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

  return(counts) # print c
  
}


#Check if number of individuals is the problem 

inds_full <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(inds_full) <- c("new_spcc", "n", "site_name")

for(i in sites){ #do all sites at once and send into df 
  
  inds <- find_ind_spcc(site = i, 
                   threshold_noneuc = 8, 
                   threshold_euc = 8) #function above
  
  inds <- inds %>% mutate(site_name = i) #set new column called site_name for site 
  
  inds_full <- rbind(inds_full, inds) #bind into larger df 
}


inds_full <- inds_full %>% 
  rename(from = to) 

full_df_spcc2 <- full_df_spcc2 %>%
  left_join(inds_full, by = c("site", "from")) 

full_df_spcc2 <- full_df_spcc2 %>% 
  rename(ind_from = n)





large_ci <- full_df_spcc2 %>% 
  filter(range_ci > 5)



df_spcc <- full_df_spcc2 %>% 
  filter(range_ci < 5)

df_spcc %>% count(species_from)

colours <- hcl.colors(7, "Viridis")

## Co/dominant 

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
