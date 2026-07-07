library(ggplot2)
library(dplyr)
library(ppjsdm)



# Using this results df 
df <- df_add3_5_t15

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


#### Two kinds of weighted means we can compute: 
#1. Mean weighted by the standard error of the coefficient estimate
#3. Mean weighted by the span of the confidence interval 

## 1. Standard error weighting 
w1 <- within %>% 
  group_by(group) %>% 
  summarise(w.mean =  (alpha, se), 
            w.mean.inv.se = weighted.mean(alpha, (1/se)), 
            mean = mean(alpha)) %>%
  ungroup()

b1 <- between %>% 
  group_by(group) %>% 
  summarise(w.mean = weighted.mean(alpha, se), 
            w.mean.inv.se = weighted.mean(alpha, (1/se)), 
            mean = mean(alpha)) %>% 
  ungroup()

b1 <- between %>% 
  group_by(group) %>% 
  summarise(median = median(alpha)) %>% 
  ungroup()

  
all1 <- df %>%  
  group_by(group, site) %>% 
  summarise(w.mean = weighted.mean(alpha, se)) %>% 
  ungroup()

# lm(data = within, 
#    alpha ~ group, #This is not working because all groups need to be intercept and due to current data structure is not working
#    weights = se) #can change data structure, but weighted.mean fun does the same thing


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
# within_wmean <- cbind(w1, w3$w.mean)
# colnames(within_wmean) <- c("group", "se",  "ci_range")
# 
# 
# # For between only 
# bw_wmean <- cbind(b1, b3$w.mean)
# colnames(bw_wmean) <- c("group", "se",  "ci_range")
# 
# # For all data (between site)
# all_wmean <- cbind(all1,  all3$w.mean)
# colnames(all_wmean) <- c("group", "site", "se", "ci_range")


#Calculate weighted mean manually as small-small answer is weird

#the formual is sum(x*w) / sum(w)
w1 <- within %>% 
  group_by(group) %>% 
  summarise(w.mean = sum(alpha*se)/ sum(se)) %>%
  ungroup()

#just filter to subcanopy small 

d <- df %>% 
  filter(group == "Subcanopy.small_Subcanopy.small")

up <- sum(d$alpha*d$se)
lo <- sum(d$se)
wmean <- up/lo #-0.067 compared to -0.087







###Sig testing for w1
w1 <- w1 %>% 
  mutate(class_int = c("large_large", "large_small", "small_small", "large_large", "large_small", "small_small"))

w1 <- w1 %>% 
  mutate(fg = c("Canopy","Canopy", "Canopy" ,"Subcanopy", "Subcanopy", "Subcanopy"))

model <- aov(w.mean ~ class_int * fg, data = w1)
summary(model)

plot(model)         # residuals vs fitted and QQ plot

shapiro.test(residuals(model))   # normality

