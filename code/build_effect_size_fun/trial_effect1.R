library(ppjsdm)
library(dplyr)

############## EFFECT THRESHOLDS 
# It is important to understand what coefficient size corresponds to a negligible,
#small, moderate and large interaction. 

# If we think that for that alpha, an individuals located at the same location as 
#a focal individual has a 10-20% effect on the likelihood of the focal individual being located there.
# i.e. exp(alpha) = 1.2 which is equivalent to alpha = 0.18. 
#This corresponds to 0.095 < alpha < 0.18, as you wrote for a small interaction
#Then it would follow that an alpha <0.095 is neligible, 0.18 < alpha < 0.5 is moderate 
#However this ignores distance 

#To include distance: 
#If we think that for that alpha, an individuals from the typical distance from a focal individual 
# has a 10-20% effect on the likelihood of the focal individual being located there. 
#In this case, exp(0.5 alpha) = 1.2, which corresponds to alpha = 2 * 0.18. 
#That means we are doubling the effect thresholds 

#This works only when saturation N = 1,2 
#N is proportional is alpha, an increase in N is an increase in alpha 

#. For a scarce species (think eucalypts) there might be only one or two individuals close by. 
#So an alpha = 0.18 will indeed only increase the conditional intensity by 20%, as expected, or at most exp(2 alpha) = +40%. 
#In contrast, for an abundant species (understorey) there might always be 20 individuals close by. 
#So in practice, the effect on a focal individual of all other individuals will be exp(alpha * 20) > 36
#which is clearly not small. In this case, even a small value of alpha = 0.18 has a massive effect on the 
#probability of occurrence of a given individual, even though each individual does not contribute much. 


#Cool, so we either use incredibly specific language or reduce saturation parameter 

#OR, we can test this using conditional intensities to understand what is happening 

#F(alpha, i, j):

# Compute A- = Papangelou conditional intensity computed at the locations of species i,
#conditional on the entire configuration except species j, fixing the short-distance interaction coefficient to alpha.

# Compute A+ = Papangelou conditional intensity computed at the locations of species i, 
#conditional on the entire configuration including species j, fixing the short-distance interaction coefficient to alpha.

# Return Mean(A+ / A-).

#This function corresponds to the average effect of all individuals of species j on the conditional intensity of an
#individual of species i, when the interaction coefficient is set to alpha. 

#The function is increasing in alpha when alpha > 0 and decreasing if alpha < 0. 
#Solve F(alpha, i, j) = (1.1, 1.2, 1.5) to generate the thresholds you are looking for, specific to the dataset you are analysing. 

#Crucially, this method does more than provide rough guidelines. 
#Instead, it tells you for the dataset you are analysing, and given your model parameters (potential shape, distance, 
#saturation parameter, etc), what the alpha thresholds are.


#### Take Bruxtner as an example

#Load Data
data <- read.csv("data/data_cleaned.csv")

df <- data %>%
  filter(Site_Name == "WaratahMix")


#jitter coordinates
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
  filter(!(observation_count < 16))

d %>% count(new_group)

#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
plot(configuration)


#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))

#fit model 
fit<- ppjsdm::gibbsm(configuration = configuration, #do the fit 
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


sum <- summary(fit)
sum

df <- make_sum_df(list(fit), list(sum), coefficient = "alpha")

#get beta0
b0 <- as.numeric(fit$coefficients$beta0)

#In WaratahMix, the difference in the effect for interactions x_y v y_x is max 0.129, min 0.00177, mean 0.036 

test3 <- effect_size(configuration, 
                     fit)

ggplot(data = test3, 
       aes(x = alpha, 
           y = mean_effect)) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray60") + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_point(size = 2) + 
  scale_y_continuous(breaks = seq(0, 2.65, 0.25)) +
  theme_classic() + 
  ylab("Change in likelihood of occurrence (effect)") + 
  xlab("Alpha coefficient")


#test baselines
plot_papangelou(fit, 
                window = window, 
                configuration = configuration, 
                type = "Allocasuarina torulosa small", 
                use_log = TRUE, 
                drop_type_from_configuration = TRUE, 
                show = "Allocasuarina torulosa small")



pap <- compute_papangelou(fit, 
                   x = c(0:100), #for 1x1m pixels 
                   y = c(0:100),
                   configuration = configuration, 
                   type = "Allocasuarina torulosa small", 
                   use_log = TRUE)

#compute papangelou at the locations of allocasuarina tortulosa small 
loc <- d %>% filter(new_group == "Allocasuarina torulosa small") %>% 
  select(x_jitter, y_jitter)


pap <- compute_papangelou(fit, 
                          x = loc$x_jitter,
                          y = loc$y_jitter,
                          configuration = configuration, 
                          type = "Allocasuarina torulosa small", 
                          use_log = TRUE)

loc$condint <- pap


#let's set alpha to 0.10 
#compute papangelou conditional intensity for species i conditional on config except species j fixing the short-distance coefs to alpha

#We need to build the configuration without species j (two understorey species)
# species i = "Ceratopetalum apetalum small"
# species j = "Archontophoenix cunninghamiana small"

full_configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
plot(configuration)

noj <- d %>% filter(! new_group == "Archontophoenix cunninghamiana small")
j_configuration <- ppjsdm::Configuration(noj$x_jitter, noj$y_jitter, types = noj$new_group)
plot(j_configuration)

#get a df for the locations of species i 
locs <- d %>% filter(new_group == "Ceratopetalum apetalum small") %>% 
  select(x_jitter, y_jitter)

#Compute A- (without species j)
pap1 <- compute_papangelou(configuration = j_configuration, 
                           x = locs$x_jitter, 
                           y = locs$y_jitter, 
                           type = "Ceratopetalum apetalum small", 
                           mark = 0, 
                           beta0 = b0[-2],
                           alpha = matrix(0.1, nspecies-1, nspecies-1),
                           short_range = matrix(10, nspecies-1, nspecies-1), 
                           saturation = 10, 
                           model = "exponential"
                           )

#add to locs df 
locs$withoutj <- pap1

#Compute A+ (including species j)
pap2 <- compute_papangelou(configuration = configuration, 
                           x = locs$x_jitter, 
                           y = locs$y_jitter, 
                           type = "Ceratopetalum apetalum small", 
                           mark = 0, 
                           beta0 = b0, 
                           alpha = matrix(0.1, nspecies, nspecies),
                           short_range = matrix(10, nspecies, nspecies), 
                           saturation = 10, 
                           model = "exponential"
)


locs$withj <- pap2

#compute mean difference 
mean(locs$withj / locs$withoutj) #numerator is the change 
#1.314192

#In the Bruxner plot, there is a 30% change in the conditional likelihood of predicting a individual of species i
#correctly when individuals of species j are included when alpha is set to 0.1




#Okay, but we need to check this for all values of alpha 
locs # has the locations of species i 

alpha_list <- list() #empty list
alpha_vect <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5) #alphas we want to check 

for(a in seq_along(alpha_vect)){
#Compute A- (without species j)
  alpha <- alpha_vect[a]
withoutj <- compute_papangelou(configuration = j_configuration, 
                           x = locs$x_jitter, 
                           y = locs$y_jitter, 
                           type = "Ceratopetalum apetalum small", 
                           mark = 0, 
                           drop_type_from_configuration = FALSE,
                           alpha = matrix(alpha, nspecies-1, nspecies-1),
                           short_range = matrix(10, nspecies-1, nspecies-1), 
                           saturation = 10, 
                           model = "exponential"
)

alpha_list[[a]] <- withoutj

}

df_without <- data.frame(matrix(unlist(alpha_list), nrow = length(alpha_list[[1]]), byrow = FALSE)) 

colnames(df_without) <- paste0(alpha_vect, "_withoutj")

#Compute A+ (including species j)
alpha_list <- list()
for(a in seq_along(alpha_vect)){
alpha <- alpha_vect[a]
withj <- compute_papangelou(configuration = configuration, 
                           x = locs$x_jitter, 
                           y = locs$y_jitter, 
                           type = "Ceratopetalum apetalum small", 
                           mark = 0, 
                           alpha = matrix(alpha, nspecies, nspecies),
                           short_range = matrix(10, nspecies, nspecies), 
                           saturation = 10, 
                           model = "exponential"
)

alpha_list[[a]] <- withj 

} 

df_with <- data.frame(matrix(unlist(alpha_list), nrow = length(alpha_list[[1]]), byrow = FALSE)) 

colnames(df_with) <- paste0(alpha_vect, "_withj")

#add into locs df 

locs2 <- cbind(locs, df_with, df_without)

#Cool, now need to get mean change for each alpha 

mean(locs2$`0.05_withj`/locs2$`0.05_withoutj`)
mean(locs2$`0.1_withj`/locs2$`0.1_withoutj`)
mean(locs2$`0.2_withj`/locs2$`0.2_withoutj`)
mean(locs2$`0.3_withj`/locs2$`0.3_withoutj`)
mean(locs2$`0.4_withj`/locs2$`0.4_withoutj`)
mean(locs2$`0.5_withj`/locs2$`0.5_withoutj`)


#Cool, but this is for two understorey species, what we need to check now is for two large groups 
#Need to grab another site 
#Site = WaratahMix 


df <- data %>%
  filter(Site_Name == "WaratahMix")


#jitter coordinates
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
  filter(!(observation_count < 14))

d %>% count(new_group)



#make config
configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
plot(configuration)


#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(configuration$types))

#fit model 
fit_w <- ppjsdm::gibbsm(configuration = configuration, #do the fit 
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


sum_w <- summary(fit_w)
sum_w

#test baselines
plot_papangelou(fit_w, 
                window = window, 
                configuration = configuration, 
                type = "Eucalyptus obliqua large", 
                use_log = TRUE, 
                drop_type_from_configuration = TRUE, 
                show = "Eucalyptus obliqua large")




#We need to build the configuration without species j (two understorey species)
# species i = "Eucalyptus obliqua large"
# species j = "Eucalyptus fastigata large"

full_configuration <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
plot(configuration)

nojw <- d %>% filter(! new_group == "Eucalyptus obliqua small")
nojw_configuration <- ppjsdm::Configuration(nojw$x_jitter, nojw$y_jitter, types = nojw$new_group)
plot(nojw_configuration)

#get a df for the locations of species i 
locs_w <- d %>% filter(new_group == "Eucalyptus obliqua large") %>% 
  select(x_jitter, y_jitter)


#Run through the loop

alpha_list <- list() #empty list
alpha_vect <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5) #alphas we want to check 

for(a in seq_along(alpha_vect)){
  #Compute A- (without species j)
  alpha <- alpha_vect[a]
  withoutj <- compute_papangelou(configuration = nojw_configuration, 
                                 x = c(1:100), 
                                 y = c(1:100), 
                                 type = "Eucalyptus obliqua large", 
                                 mark = 1, 
                                 drop_type_from_configuration = FALSE,
                                 alpha = matrix(alpha, nspecies-1, nspecies-1),
                                 short_range = matrix(10, nspecies-1, nspecies-1), 
                                 saturation = 10, 
                                 model = "exponential"
  )
  
  alpha_list[[a]] <- withoutj
  
}

dfw_without <- data.frame(matrix(unlist(alpha_list), nrow = length(alpha_list[[1]]), byrow = FALSE)) 

colnames(dfw_without) <- paste0(alpha_vect, "_withoutj")

#Compute A+ (including species j)
alpha_list <- list()
for(a in seq_along(alpha_vect)){
  alpha <- alpha_vect[a]
  withj <- compute_papangelou(configuration = full_configuration, 
                              x = locs_w$x_jitter, 
                              y = locs_w$y_jitter, 
                              type = "Eucalyptus obliqua large", 
                              mark = 1, 
                              alpha = matrix(alpha, nspecies, nspecies),
                              short_range = matrix(10, nspecies, nspecies), 
                              saturation = 10, 
                              model = "exponential"
  )
  
  alpha_list[[a]] <- withj 
  
} 

dfw_with <- data.frame(matrix(unlist(alpha_list), nrow = length(alpha_list[[1]]), byrow = FALSE)) 

colnames(dfw_with) <- paste0(alpha_vect, "_withj")

#add into locs df 

locs2_w <- cbind(locs_w, dfw_with, dfw_without)


mean(locs2_w$`0.05_withj`/locs2_w$`0.05_withoutj`)
mean(locs2_w$`0.1_withj`/locs2_w$`0.1_withoutj`)
mean(locs2_w$`0.2_withj`/locs2_w$`0.2_withoutj`)
mean(locs2_w$`0.3_withj`/locs2_w$`0.3_withoutj`)
mean(locs2_w$`0.4_withj`/locs2_w$`0.4_withoutj`)
mean(locs2_w$`0.5_withj`/locs2_w$`0.5_withoutj`)





mean(locs2$`0.05_withj`/locs2$`0.05_withoutj`)
mean(locs2$`0.1_withj`/locs2$`0.1_withoutj`)
mean(locs2$`0.2_withj`/locs2$`0.2_withoutj`)
mean(locs2$`0.3_withj`/locs2$`0.3_withoutj`)
mean(locs2$`0.4_withj`/locs2$`0.4_withoutj`)
mean(locs2$`0.5_withj`/locs2$`0.5_withoutj`)




##############################################################################
### New version: Use fitted alpha
# Instead of prescribing alphas to the computation of the papangelou function 
# we can use the fitted values instead to understand if the effect is relatively high or low 
# We have already fitted the alpha at this point therefore, it makes more sense to understand 
# the change in othere interactions when a species is dropped 

#Using the same examples as above: 
#Waratah

# species i = "Ceratopetalum apetalum small"
# species j = "Archontophoenix cunninghamiana small"

window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

#get a df for the locations of species i 
locs <- d %>% filter(new_group == "Eucalyptus obliqua small") %>% 
  select(x_jitter, y_jitter)


#Make A+ : where both species i and j are present
#make config
full_config <- ppjsdm::Configuration(d$x_jitter, d$y_jitter, types = d$new_group)
plot(full_config)

nspecies <- length(levels(full_config$types))

#fit model 
full_fit <- ppjsdm::gibbsm(configuration = full_config, #do the fit 
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

coef(full_fit)

a <- compute_papangelou(fit = full_fit,
                        configuration = full_config,
                        x = coords$Var1, 
                        y = coords$Var2, 
                        type = "Eucalyptus obliqua small", 
                        drop_type_from_configuration = T, 
                        nthreads = 3)

locs$a <- a
#Now need to compute papangelou for species i without species j 

#make a new config
all_types <- unique(full_config$types)
drop_types <- setdiff(all_types, "Eucalyptus obliqua large")
drop_config <- full_config[drop_types]
plot(drop_config)

# d2 <- d %>% filter(! new_group ==  "Eucalyptus obliqua large")
# config2 <- ppjsdm::Configuration(d2$x_jitter, d2$y_jitter, d2$new_group)


short <- unname(full_fit$coefficients$short_range[[1]])
short <- short[-1, -1]


dropped_fit <-  ppjsdm::gibbsm(configuration = drop_config, #do the fit 
                              window = window,
                              short_range = short, 
                              model = "exponential",
                              fitting_package = "glmnet",
                              saturation = 10, 
                              nthreads = 4, 
                              use_regularization = F,
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)

drop_fit <- ppjsdm::gibbsm(configuration = drop_config, #do the fit 
                           window = full_fit$window,
                           short_range = short, 
                           model = full_fit$parameters$model[[1]],
                           fitting_package = full_fit$fit_algorithm,
                           saturation = full_fit$parameters$saturation, 
                           nthreads = full_fit$nthreads, 
                           use_regularization = full_fit$used_regularization,
                           dummy_distribution = full_fit$dummy_distribution,
                           min_dummy = 1, dummy_factor = 1e10, 
                           max_dummy = 1e3)
coef(dropped_fit)


drop_a <- compute_papangelou(fit = dropped_fit,
                             configuration = drop_config,
                             x = coords$Var1, 
                             y = coords$Var2, 
                             type = "Eucalyptus obliqua small", 
                             drop_type_from_configuration = F, 
                             nthreads = 3)
mean(a/drop_a)

x <- 1:100
y <- 1:100
coords <- expand.grid(x, y)
pap <- data.frame(x = coords$Var1, 
                  y = coords$Var2,  
                  a = a, 
                  drop_a = drop_a)



drop_type <- plot_papangelou(fit = drop_fit,
                             type = "Eucalyptus obliqua small", 
                             show = "Eucalyptus obliqua small", 
                             drop_type_from_configuration = T, 
                             use_log = T,
                             nthreads = 3, 
                             limits = c(-3.5, -6),
                             legend_title = "Dropped")

no_drop <- plot_papangelou(fit = drop_fit,
                           type = "Eucalyptus obliqua small", 
                           show = "Eucalyptus obliqua small", 
                           drop_type_from_configuration = F, 
                           use_log = TRUE, 
                           nthreads = 3, 
                           limits = c(-3.5, -6),
                           legend_title = "Not dropped")

drop_type + no_drop + plot_layout(guides = "collect")

#Add to df 
locs$a <- 
locs$no_j <- drop_a

locs$diff <- locs$a - locs$no_a
mean(locs$a/locs$no_a)


drop <- plot_papangelou(fit = drop_fit, 
                type = "Ceratopetalum apetalum small", 
                drop_type_from_configuration = F, 
                show = "Ceratopetalum apetalum small", 
                use_log = TRUE)


full <- plot_papangelou(fit = full_fit, 
                 type = "Ceratopetalum apetalum small", 
                 drop_type_from_configuration = F, 
                 show = "Ceratopetalum apetalum small", 
                 use_log = TRUE)

library(patchwork)
drop + full + plot_layout(guides = "collect")
