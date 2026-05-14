library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)

set.seed(1)

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
                       threshold = 10, #threshold to exclude groups
                       short_range = 10){ #short range value 
  
  d <- data_cleaned %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    group_by(Genus_Species) %>% 
    mutate(observation_count = n()) %>% 
    filter(!(observation_count < threshold)) 
  
  
  counts <- d %>% count(Genus_Species)
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
  
  
  configuration<- Configuration(d$x_jitter, d$y_jitter, types = d$Genus_Species) #make config 
  
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
                         threshold= 15, #threshold to exclude groups
                         short_range = 14) #run site_ccmodel function
  
  working_df <- make_sum_df(fits = list(sp_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(sp_mod$sum))
  
  working_df <- working_df %>% mutate(site = i) #set new column for site name 
  
  full_df_sp <- rbind(full_df_sp, working_df) #bind into larger df 
  
}


##Okay need to check unreasonable CIs 
full_df_sp2 <- full_df_sp %>% 
  mutate(range_ci = hi - lo)

large_ci <- full_df_sp2 %>% 
  filter(range_ci > 5) #There's nothing!! yay


################################################################################
### VISUALISATION 
full_df_sp2 %>% count(new_from)

# For simplicity grouping all the Syncarpai glomulifera subsp together 
full_df_sp2 <- full_df_sp %>% 
  mutate(new_from = if_else(from %in% c("Syncarpia glomulifera", 
                                        "Syncarpia glomulifera subsp glabra", 
                                        "Syncarpia glomulifera subsp glomulifera"), 
                            "Syncarpia glomulifera", from)) %>% 
  mutate(new_to = if_else(to %in% c("Syncarpia glomulifera", 
                                        "Syncarpia glomulifera subsp glabra", 
                                        "Syncarpia glomulifera subsp glomulifera"), 
                            "Syncarpia glomulifera", to))

#add geogrpahic region 
full_df_sp2 <- full_df_sp2 %>% 
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
  mutate(region = if_else(region %in% c("Flowerdale", "Dip", "Mackenzie", "Caveside"), 
                          "N_TAS", region)) %>% 
  mutate(region = if_else(region %in% c("Bird", "Supersite", "NorthStyx", "Weld", "MtField", "ZigZag"), 
                          "S_TAS", region)) %>% 
  mutate(region = if_else(region %in% c("BlackRiver", "BondTier"), 
                          "King_TAS", region)) %>% 
  mutate(region = if_else(region %in% c("BenRidge", "MtMaurice"), 
                          "Ben_TAS", region))


full_df_sp2 <- full_df_sp2%>% 
  mutate(region2 = region) %>% 
  mutate(region2 = if_else(region2 %in% c("N_NSW", "QLD"),
                           "N_AUS", region2)) %>%
  mutate(region2 = if_else(region2 %in% c("SE_NSW", "SE_VIC","Ben_TAS", "King_TAS", "N_TAS", "S_TAS"), 
                           "SE_AUS", region2)) 

#add sig colours 
full_df_sp2 <- full_df_sp2 %>% 
  mutate(sig = ifelse((lo > 0 | hi < 0), 1, NA )) %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region2), NA))


## Save out 
write.csv(full_df_sp2, "final_allsp_mod.csv")


#### EUCS 
eucs_sp <- full_df_sp2 %>% 
  filter(! from == "Non-euc") %>% 
  filter(! to == "Non-euc")


## Within-species Interactions 

sp_within <- full_df_sp2 %>% 
  filter(new_from == new_to)


sp <- ggplot(data = sp_within, 
       aes(x = alpha, 
           y = new_from, 
           colour = region2, 
           fill = fill_col)) + 
  geom_point(shape = 21, size = 3.5) +
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") + 
  scale_colour_manual(values = c("#BB9DEA", "#D3A263",  "#00BFC8"), 
                      name = "Geographic Region", 
                      guide = "none") + 
  scale_fill_manual(values = c("#BB9DEA", "#D3A263",  "#00BFC8"), 
                    na.value = "white", 
                    label = c("Northeastern Aus", "Southeastern Aus", "Western Aus", "Insignificant"),
                    name = "Region") + 
  ylab("Species") + xlab("Within-species alpha coefficient") +
  theme_bw()


ggsave("sp_within.png", sp, 
       dpi = 300, height = 8)
#Cool, but i want to only look at species with more than 1 site present 
#let's count 
eucs_within %>% count(new_from)

#I only want to think about interactions species with more than 1 site 
eucs_within2 <- eucs_within %>% 
  group_by(new_from) %>%
  filter(n() > 1) %>%   # keep only species that appear more than once
  ungroup()

eucs_within2 %>% count(new_from)


#Let's make the plot 
ggplot(data = eucs_within2, 
       aes(x = alpha, 
           y = new_from, 
           colour = region, 
           fill = fill_col)) + 
  geom_point(shape = 21, size = 3.5) +
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") + 
  scale_colour_manual(values = c("#ED90A4", "#D3A263",  "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                      name = "Geographic Region", 
                      guide = "none") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263",  "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                    na.value = "white", 
                    name = "Region") + 
  ylab("Species") + xlab("Within-species alpha coefficient") +
  theme_bw()


### I want to look at E. delegatensis and E. obliqua up close before moving on 
### E. delegatenesis 
del <- eucs_within %>% 
  filter(new_from == "Eucalyptus delegatensis") %>% 
  mutate(region = factor(region, levels = c("Ben_TAS", "N_TAS", "S_TAS"))) %>%
  mutate(site = fct_reorder(site, as.numeric(region)))

ggplot(del, aes(x = alpha, 
                  y = site, 
                  colour = region, 
                  fill = fill_col)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, 
                    xmax = hi),
                size = 1.5) + 
  geom_point(shape = 21,
             size = 3) + 
  scale_colour_manual(values = c("#DC9E70", "#1DC199", "#BB9DEA"), 
                      name = "Region") + 
  scale_fill_manual(values = c("#DC9E70", "#1DC199", "#BB9DEA"), 
                    na.value = "white", 
                    guide = "none") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("E. delegatenesis Within-Species") +
  theme_bw()  


### E. obliqua 
ob <- eucs_within %>% 
  filter(new_from == "Eucalyptus obliqua") %>%
  mutate(region = factor(region, levels = c("King_TAS", "N_TAS", "S_TAS", "SE_NSW"))) %>%
  mutate(site = fct_reorder(site, as.numeric(region)))


ggplot(ob, aes(x = alpha, 
                y = site, 
                colour = region, 
                fill = fill_col)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, 
                    xmax = hi),
                size = 1.5) + 
  geom_point(shape = 21,
             size = 3) + 
  scale_colour_manual(values = c("#ED90A4", "#B7AE50", "#1DC199", "#BB9DEA"), 
                      name = "Region") + 
  scale_fill_manual(values = c("#ED90A4", "#B7AE50", "#1DC199", "#BB9DEA"), 
                    na.value = "white", 
                    guide = "none") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("E. obliqua Within-Species") +
  theme_bw()  



## Between-species Interactions 
#There's even less going on here, not many interactions that are significant 
#Probably because there is not many sites that have good amounts of more than one species 
eucs_bw <- eucs_sp %>% 
  filter(!(new_from == new_to)) %>% 
  mutate(int = paste0(new_from, sep = "_", new_to))

eucs_bw %>% count(int)

ggplot(data = eucs_bw, 
       aes(x = alpha, 
           y = int, 
           colour = region, 
           fill = fill_col)) + 
  geom_point(shape = 21, size = 3.5) +
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") + 
  scale_colour_manual(values = c("#ED90A4", "#D3A263",  "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                      name = "Geographic Region", 
                      guide = "none") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263",  "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                    na.value = "white", 
                    name = "Region")+
  ylab("Interaction") + xlab("Alpha coefficient") +
  theme_bw()


### NON-EUCS TO EUCS 
# Need to check what's happening here 

noneucs_sp <- full_df_sp2 %>% 
  filter(to == "Non-euc" | from == "Non-euc")


#Doesn't seem to be much of a pattern
noneucs_bw <- noneucs_sp %>% 
  filter(!(to == from))

ggplot(data = noneucs_bw, 
       aes(x = alpha, 
           y = region, 
           colour = as.factor(sig))) + 
  geom_point() +
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") + 
  scale_color_manual(values = "black",
                     na.value = "gray75",
                     name = "Significance") +
  theme_bw() + 
  ggtitle("Non-euc to Euc")



### Non-eucs to non-eucs 
noneucs_within <- noneucs_sp %>% 
  filter(to == from)

noneucs_within %>% count(site) #ANU363 only has E. regnans (6 non-eucs)

noneucs_within <- noneucs_within %>% #order data 
  arrange(region, alpha) %>% 
  mutate(site = factor(site, levels = site)) 

#There's a lot of clustering between non-eucs, no geographic spilt 
ggplot(data = noneucs_within, aes(x = alpha, 
                y = site, 
                colour = region, 
                fill = fill_col)) + 
  geom_vline(xintercept = 0, 
             colour = "red", 
             linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, 
                    xmax = hi),
                size = 0.75) + 
  geom_point(shape = 21,
             size = 3) + 
  scale_colour_manual(values = c("#ED90A4", "#D3A263",  "#B7AE50", "#7EBA68", "#1DC199", "#00BFC8", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                      name = "Geographic Region") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263",  "#B7AE50", "#7EBA68", "#1DC199", "#6FB1E7", "#BB9DEA", "#E48FD1"), 
                    na.value = "white", 
                    guide = "none") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "") + 
  ggtitle("Non-euc to Non-euc") +
  theme_bw()


### Summary Statistics 
d <- final_sp_mod %>% 
  filter(!new_from == "Non-euc") %>% 
  filter(!new_to == "Non-euc") 

d <- d %>%  
  filter(new_from == new_to)

d <- d %>% 
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
  mutate(region = if_else(region %in% c("BondTier", "BlackRiver", "Weld", "Weld", "MtField","Bird", "ZigZag", "Flowerdale", "Dip", "Supersite"), 
                          "o_TAS", region)) 

sum <- d %>% 
  group_by(region, new_from) %>% 
  summarise(
    count = n(), 
    median_a = median(alpha), 
    mean_a = mean(alpha), 
    min_a = min(alpha), 
    max_a = max(alpha)
  )





###########################################################################################

### Prediction 

# Need to understand performance of the model 

#make empty df 
auc_df_sp <- data.frame(matrix(ncol = 4))
colnames(auc_df_sp) <- c("site", "types", "auc_score", "n_points")

auc_working <- auc_df_sp

auc_species <- function(site, 
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
sites <- sites[!sites %in% c("ANU363")] #ANU363 only has 1 group (E. regnans)

#run 
for(i in sites){ 
  
  auc_working <- auc_species(site = i, 
                               threshold_noneuc = 8, 
                               threshold_euc = 10, 
                               short_range = 8)
  
  auc_df_sp <- rbind(auc_df_sp, auc_working) #bind into larger df 
  
}


auc_df_sp2 <- data.frame(auc_df_sp, row.names = NULL)
auc_df_sp2 <- auc_df_sp2 %>% 
  filter(!is.na(site))

auc_df_sp2 %>% count(site)

write.csv(auc_df_sp2, "auc_species.csv")

median(auc_df_sp2$auc_score)

#only care about euc scores so 
e <- auc_df_sp2 %>% 
  filter(str_starts(types, "E|C|S"))

e %>% count(types)

median(e$auc_score)


#Visualisation 
sp <- ggplot(data = auc_df_sp2, 
               aes(x = class, 
                   y = auc_score)) + 
  geom_boxplot() + 
  geom_point() + 
  ylim(0.4, 0.9) + 
  ggtitle("SpCC")

sups









### Model stats 
d <- final_sp_mod %>% 
  filter(!new_from == "Non-euc") %>% 
  filter(!new_to == "Non-euc")

d %>% count(sig)
d <- final_sp_mod %>% count(region)

strong <- final_sp_mod %>% 
  filter(alpha > 0.49 | alpha < -0.49)



###########
# Prediction

########### Cross validation
#As we know, the AUC method used above doesn't allow for any within-group interactions to be evaluated 
#Therefore we can use leave-n-individuals out cross validation

sp <- function(site, 
               threshold = 10){ #short range value 
  d <- data_cleaned %>% filter(Site_Name == site)
  
  
  d <- d %>% 
    mutate(species = if_else(!str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia|^Syncarpia")), "Non-euc", Genus_Species)) %>% 
    group_by(species) %>% 
    mutate(observation_count = n()) %>% 
    filter(!(str_starts(species, "Non-euc") & observation_count < 10)) %>% 
    filter(!(str_starts(species,  regex("^Eucalyptus|^Corymbia|^Syncarpia")) & observation_count < threshold)) 
  
  
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
  
}



## For all sites: 
a <- unique(ausplot_data$Site_Name) 
#sites <- sites[!sites %in% c("Dip", "Bird", "Flowerdale")]

                 

full_auc_ind <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(full_auc_ind) <- c("iter", "types", "value", "split", "site")


a <- ind_validation2(site = "BlackRiver", 
                     n_ind = 10, 
                     n_it = 20)

a$site <-"BlackRiver" #set column site to k 

full_auc_ind <- rbind(full_auc_ind, a) #bind into larger df 

full_auc_ind %>% count(site)



write.csv(full_auc_ind, "auc_inds_species.csv")

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
