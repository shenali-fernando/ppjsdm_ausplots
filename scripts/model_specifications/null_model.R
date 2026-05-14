library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

####################### NULL MODEL ##########################
# In the case of the ppjsdm model, the null model is where individual do not have an identity.
# This means that the individuals do not have a species, size, or any type of label - all individuals are in the SAME, ONE group. 
# i.e. We do not set a type when setting the configuration 
# This allows us to understand if there is more or less clustering or replusion occurring compared to the Poisson point process where all points are independent of each other 
# Furthermore we can understand which sites have more replusion or clustering than others 


# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")

### Null model - all individuals in the same group 'default', i.e. no species or crown classes

#get rid of duplicate points first 
data_c <- data_c %>%
  group_by(Site_Name, Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) %>%
  ungroup() 

site_names <- c(unique(data_c$Site_Name)) #vector of site names
null_alpha <- data.frame(matrix(ncol = 4, nrow = length(site_names))) #make dataframe to put things in
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
colnames(null_alpha) <- c("site", "estimate", "lo", "hi") #site name, alpha coef, lo CI value, high CI value

for(i in seq_along(site_names)){
  
  site <- site_names[i]
  
  df <- data_c %>% filter(Site_Name == site) #filter data to site
  
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
  sum <- summary(fit)
  
  null_alpha$site[i] <- site 
  null_alpha$estimate[i] <- fit$coefficients$alpha[[1]] #extract values from summary
  
  null_alpha$lo[i] <- sum$lo$alpha[[1]]
  null_alpha$hi[i] <- sum$hi$alpha[[1]]
  
  print(null_alpha)
  
}




# manually adding in what dominant eucs in each site
null_alpha <- null_alpha %>% 
  mutate(dom_euc = site) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Carey", "Dombakup", "Warren", "Dawson", "Giants", "Sutton", "Frankland", "Clare", "Collins"), 
                           "E. diversicolor", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Lardner", "NorthStyx", "Turtons", "Weeaproinah", "Weld"), 
                           "E. regnans", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Newline", "WaratahMix", "WogWay", "Goodenia"), 
                           "E. fastigata", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank"), 
                           "E. pilularis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Bruxner", "Osullivans", "Baldy", "Herberton", "Koombooloomba", "Lamb Range"), 
                           "E. grandis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Flowerdale", "Dip", "Bird", "Supersite", "ZigZag", "BlackRiver", "BondTier", "Candelo"), 
                           "E. obliqua", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("BenRidge", "Caveside", "MtMaurice", "Mackenzie", "MtField"), 
                           "E. delegatensis", dom_euc))
  
#making the order right
data_ordered <- null_alpha %>%
  arrange(dom_euc, estimate) %>% 
  mutate(site = factor(site, levels = site)) 

## making box-plot
ggplot(data_ordered, aes(x = estimate, y = site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "Site",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("All Species") +
  theme_bw()





##########################################################################
################ Can look at a null model with ONLY EUCS 

#Same code; but create new column in original data of Genus
euc_only <- data_c %>% mutate(Genus = Genus_Species) %>% 
  mutate(Genus = sub(" .*", "", Genus)) %>% 
  filter(Genus %in% c("Eucalyptus", "Corymbia"))

#checks 
euc_only %>% count(Genus)
euc_only %>% count(Genus_Species)

# get rid of duplicate points
euc_only <- euc_only %>%
  group_by(Site_Name, Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) %>%
  ungroup() 


site_names <- c(unique(euc_only$Site_Name)) #vector of site names
null_alpha_onlyeuc <- data.frame(matrix(ncol = 4, nrow = length(site_names))) #make dataframe to put things in
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
colnames(null_alpha_onlyeuc) <- c("site", "estimate", "lo", "hi") #site name, alpha coef, lo CI value, high CI value

for(i in seq_along(site_names)){
  
  site <- site_names[i]
  
  df <- euc_only %>% filter(Site_Name == site) #filter data to site
  
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
  sum <- summary(fit)
  
  null_alpha_onlyeuc$site[i] <- site 
  null_alpha_onlyeuc$estimate[i] <- fit$coefficients$alpha[[1]] #extract values from summary
  
  null_alpha_onlyeuc$lo[i] <- sum$lo$alpha[[1]]
  null_alpha_onlyeuc$hi[i] <- sum$hi$alpha[[1]]
  
  print(null_alpha_onlyeuc)
  
}


# manually adding in what dominant eucs in each site
null_alpha_onlyeuc <- null_alpha_onlyeuc %>% 
  mutate(dom_euc = site) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Carey", "Dombakup", "Warren", "Dawson", "Giants", "Sutton", "Frankland", "Clare", "Collins"), 
                           "E. diversicolor", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Lardner", "NorthStyx", "Turtons", "Weeaproinah", "Weld"), 
                           "E. regnans", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Newline", "WaratahMix", "WogWay", "Goodenia"), 
                           "E. fastigata", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank"), 
                           "E. pilularis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Bruxner", "Osullivans", "Baldy", "Herberton", "Koombooloomba", "Lamb Range"), 
                           "E. grandis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Flowerdale", "Dip", "Bird", "Supersite", "ZigZag", "BlackRiver", "BondTier", "Candelo"), 
                           "E. obliqua", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("BenRidge", "Caveside", "MtMaurice", "Mackenzie", "MtField"), 
                           "E. delegatensis", dom_euc))

#making the order right
data_ordered_euc <- null_alpha_onlyeuc %>%
  arrange(dom_euc, estimate) %>% 
  mutate(site = factor(site, levels = unique(site))) 



## making box-plot
ggplot(data_ordered_euc, aes(x = estimate, y = site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "Site",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Only Eucalypts") +
  theme_bw()



#Read back in with added MAT and MAP columns 
df <- read.csv("C:/Users/shena/Desktop/ausplots/null_euconly.csv")

#MAP - bin these things 
df <- df %>% mutate(bin_MAP = cut_width(MAP, width=250, boundary=0)) 

#add another factor for plotting
df$bin_MAP <- factor(df$bin_MAP, 
                     levels = c("[500, 750]", "[750,1e+03]", "(1e+03,1.25e+03]", "(1.25e+03,1.5e+03]", "(1.5e+03,1.75e+03]", "(1.75e+03,2e+03]"))
ggplot(data = df, 
       aes(x = bin_MAP, y = estimate)) + 
  geom_boxplot(position = position_nudge(x = -0.5)) +
  geom_point(aes(colour = dom_euc), 
             position = position_nudge(x = -0.5),
             size = 1.75) + 
  scale_x_discrete(drop = FALSE, labels = c(750, 1000, 1250, 1500, 1750, 2000)) +
  scale_y_reverse() +
  theme_bw()


#MAT
df <- df %>% mutate(bin_MAT = cut_width(MAT, width=1, boundary=0)) 

#add another factor for plotting
df$bin_MAT <- factor(df$bin_MAT, 
                     levels = c("[4, 6]", "[6,8]", "(8,10]", "(10,12]", "(12,14]", "(14,16]", "(16,18]", "(18,20]", "(20,22]") )


ggplot(data = df, 
       aes(x = bin_MAT, y = estimate)) + 
  geom_boxplot(position = position_nudge(x = -0.5)) +
  geom_point(aes(colour = dom_euc), 
             position = position_nudge(x = -0.5),
             size = 1.75) + 
  scale_x_discrete(drop = FALSE, labels = c(6:22)) +
  scale_y_reverse() +
  theme_bw()
