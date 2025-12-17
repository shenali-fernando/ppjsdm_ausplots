library(ggplot2)
library(dplyr)
library(ppjsdm)



# Using this results df 
df <- df_add3_5_t12

################ WEIGHTED MEANS 
### Compute a weighted mean for each size + function grouping for eahc site 

df <- df %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
    group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
    group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
    group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
    TRUE ~ group))

##### Within-species interactions 
within <- df %>% 
  filter(species_from == species_to) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int))


between <- df %>% 
  filter(!(species_from == species_to)) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int))


#### Three kinds of weighted means we can compute: 
#1. Mean weighted by the standard error of the coefficient estimate
#2. Mean weighted by if coefficient is significant or not 
#3. Mean weighted by the span of the confidence interval 

## 1. Standard error weighting 
w1 <- within %>% 
  group_by(group) %>% 
  summarise(w.mean = weighted.mean(alpha, se)) %>% 
  ungroup()

b1 <- between %>% 
  group_by(group) %>% 
  summarise(w.mean = weighted.mean(alpha, se)) %>% 
  ungroup()
  
all1 <- df %>%  
  group_by(group, site) %>% 
  summarise(w.mean = weighted.mean(alpha, se)) %>% 
  ungroup()

# lm(data = within, 
#    alpha ~ group, #This is not working because all groups need to be intercept and due to current data structure is not working
#    weights = se) #can change data structure, but weighted.mean fun does the same thing





#2. Mean weighted by if coefficient is significant or not 
#DONT USE!!!
df2 <- df %>% 
  mutate(sig = ifelse(is.na(sig) == T, 1, 2)) #make sig col 1, 2 

#within
 w2 <- df2 %>% 
  filter(species_from == species_to) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int)) %>% 
  group_by(group) %>% 
  summarise(w.mean = weighted.mean(alpha, sig)) %>% 
  ungroup()

#between 
b2 <- df2 %>% 
  filter(!(species_from == species_to)) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int)) %>% 
  group_by(group) %>% 
  summarise(w.mean = weighted.mean(alpha, sig)) %>% 
  ungroup()

#all
all2 <- df2 %>%  
  group_by(group, site) %>% 
  summarise(w.mean = weighted.mean(alpha, sig)) %>% 
  ungroup()



#3. Mean weighted by the span of the confidence interval 

#within
w3 <- df %>% 
  filter(species_from == species_to) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int)) %>% 
  group_by(group) %>% 
  summarise(w.mean = weighted.mean(alpha, range_ci)) %>% 
  ungroup()

#between 
b3 <- df %>% 
  filter(!(species_from == species_to)) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int)) %>% 
  group_by(group) %>% 
  summarise(w.mean = weighted.mean(alpha, range_ci)) %>% 
  ungroup()

#all
all3 <- df %>%  
  group_by(group, site) %>% 
  summarise(w.mean = weighted.mean(alpha, range_ci)) %>% 
  ungroup()


#### Comparison of methods 

# For within only 
within_wmean <- cbind(w1, w2$w.mean, w3$w.mean)
colnames(within_wmean) <- c("group", "se", "sig", "ci_range")


# For between only 
bw_wmean <- cbind(b1, b2$w.mean, b3$w.mean)
colnames(bw_wmean) <- c("group", "se", "sig", "ci_range")

# For all data (between site)
all_wmean <- cbind(all1, all2$w.mean, all3$w.mean)
colnames(all_wmean) <- c("group", "site", "se", "sig", "ci_range")
