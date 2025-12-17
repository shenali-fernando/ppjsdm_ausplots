library(ppjsdm)
library(dplyr)


#We want to understand for what values of alpha is there a 10%, 20%, 50% increase in the likelihood 

#Trying this we conditional intensity: 

#Taking two sites with very different spatial patterns, species and abudances: 
#Bruxner: Many clumped understorey species 
#WaratahMix: Relatively many widely-spaced large, canopy individuals 


#Load already cleaned/formatted data
data <- read.csv("brux_waratah.csv")

brux <- data %>% filter(Site_Name == "Bruxner")

#make config
full_config <- ppjsdm::Configuration(brux$x_jitter, brux$y_jitter, types = brux$new_group)
plot(full_config)


#set parameters 
window <- ppjsdm::Rectangle_window(x_range = c(0, 100), 
                                   y_range = c(0, 100))

nspecies <- length(levels(full_config$types))

#fit model 
fit_brux <- ppjsdm::gibbsm(configuration = full_config, #do the fit 
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


#get beta0
b0 <- as.numeric(fit_brux$coefficients$beta0)

#Choose two understorey species to look at: 
# species i = "Ceratopetalum apetalum small"
# species j = "Archontophoenix cunninghamiana small"

#No j config
noj <- brux %>% filter(! new_group == "Archontophoenix cunninghamiana small")
j_config <- ppjsdm::Configuration(noj$x_jitter, noj$y_jitter, types = noj$new_group)
plot(j_config)

#get a df for the locations of species i 
locs <- brux %>% filter(new_group == "Ceratopetalum apetalum small") %>% 
  select(x_jitter, y_jitter)


alpha_list <- list() #empty list
#alpha_vect <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5) #alphas we want to check 
#We want to check the sequence of alphas so: 
alpha_vect <- seq(from = 0, to = 1, length.out = 100)



#Compute A- (without species j)
for(a in seq_along(alpha_vect)){
  alpha <- alpha_vect[a]
  withoutj <- compute_papangelou(configuration = j_config, 
                                 x = locs$x_jitter, 
                                 y = locs$y_jitter, 
                                 type = "Ceratopetalum apetalum small", 
                                 mark = 0, 
                                 drop_type_from_configuration = FALSE,
                                 beta0 = b0[-2],
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
  withj <- compute_papangelou(configuration = full_config, 
                              x = locs$x_jitter, 
                              y = locs$y_jitter, 
                              type = "Ceratopetalum apetalum small", 
                              mark = 0, 
                              drop_type_from_configuration = FALSE,
                              beta0 = b0,
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





#Do the same thing for the site Waratah Mix which has many widely-spaced canopy species 

wara <- data %>% filter(Site_Name == "WaratahMix")


#make config
full_config_wara <- ppjsdm::Configuration(wara$x_jitter, wara$y_jitter, types = wara$new_group)
plot(full_config_wara)


#set parameter
nspecies <- length(levels(full_config_wara$types))

#fit model 
fit_wara <- ppjsdm::gibbsm(configuration = full_config_wara, #do the fit 
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


#get beta0
bw0 <- as.numeric(fit_wara$coefficients$beta0)


#We need to build the configuration without species j (two canopy species)
# species i = "Eucalyptus obliqua large"
# species j = "Eucalyptus fastigata large"

#Build config without species j 
nojw <- wara %>% filter(! new_group == "Eucalyptus fastigata large")
nojw_config <- ppjsdm::Configuration(nojw$x_jitter, nojw$y_jitter, types = nojw$new_group)
plot(nojw_config)


#get a df for the locations of species i 
locs_wara <- wara %>% filter(new_group == "Eucalyptus obliqua large") %>% 
  select(x_jitter, y_jitter)


#Run through the loop

alpha_list <- list() #empty list
alpha_vect <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5) #alphas we want to check 

#Compute A- (without species j)
for(a in seq_along(alpha_vect)){
  alpha <- alpha_vect[a]
  withoutj <- compute_papangelou(configuration = nojw_config, 
                                 x = locs_wara$x_jitter, 
                                 y = locs_wara$y_jitter, 
                                 type = "Eucalyptus obliqua large", 
                                 mark = 1, 
                                 drop_type_from_configuration = FALSE,
                                 beta0 = bw0[-3],
                                 alpha = matrix(alpha, nspecies-1, nspecies-1),
                                 short_range = matrix(10, nspecies-1, nspecies-1), 
                                 saturation = 10, 
                                 model = "exponential"
  )
  
  alpha_list[[a]] <- withoutj
  
}

dfwara_without <- data.frame(matrix(unlist(alpha_list), nrow = length(alpha_list[[1]]), byrow = FALSE)) 

colnames(dfwara_without) <- paste0(alpha_vect, "_withoutj")

#Compute A+ (including species j)
alpha_list <- list()
for(a in seq_along(alpha_vect)){
  alpha <- alpha_vect[a]
  withj <- compute_papangelou(configuration = full_config_wara, 
                              x = locs_wara$x_jitter, 
                              y = locs_wara$y_jitter, 
                              type = "Eucalyptus obliqua large", 
                              mark = 1, 
                              drop_type_from_configuration = FALSE,
                              alpha = matrix(alpha, nspecies, nspecies),
                              beta0 = bw0,
                              short_range = matrix(10, nspecies, nspecies), 
                              saturation = 10, 
                              model = "exponential"
  )
  
  alpha_list[[a]] <- withj 
  
} 

dfwara_with <- data.frame(matrix(unlist(alpha_list), nrow = length(alpha_list[[1]]), byrow = FALSE)) 

colnames(dfwara_with) <- paste0(alpha_vect, "_withj")

#add into locs df 

locs2_wara <- cbind(locs_wara, dfwara_with, dfwara_without)


### Waratah Mix Thresholds for E. fastigata large and E. obliqua large 
mean(locs2_wara$`0.05_withj`/locs2_wara$`0.05_withoutj`)
mean(locs2_wara$`0.1_withj`/locs2_wara$`0.1_withoutj`)
mean(locs2_wara$`0.2_withj`/locs2_wara$`0.2_withoutj`)
mean(locs2_wara$`0.3_withj`/locs2_wara$`0.3_withoutj`)
mean(locs2_wara$`0.4_withj`/locs2_wara$`0.4_withoutj`)
mean(locs2_wara$`0.5_withj`/locs2_wara$`0.5_withoutj`)



#In the WaratahMix plot: 
#there is a 5% increase in the conditional likelihood of predicting a individual of species i
#correctly when individuals of species j are included when alpha is set to 0.05?



#Compare to understorey species in Bruxner 
mean(locs2$`0.05_withj`/locs2$`0.05_withoutj`)
mean(locs2$`0.1_withj`/locs2$`0.1_withoutj`)
mean(locs2$`0.2_withj`/locs2$`0.2_withoutj`)
mean(locs2$`0.3_withj`/locs2$`0.3_withoutj`)
mean(locs2$`0.4_withj`/locs2$`0.4_withoutj`)
mean(locs2$`0.5_withj`/locs2$`0.5_withoutj`)
