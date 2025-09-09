library(ppjsdm)
library(dplyr)
library(spatstat)
library(ggplot2)
library(terra)
library(rasterVis)
library(raster)
library(maptools)
library(spatstat.explore)

#################################
### For diameter by site only (what i'm changing each time in data creation)
ind_validation1 <- function(site, 
                           n_ind = 10, #number of individuals to leave out 
                           n_it = 30) #number of iterations, default to 30
                          {
  
  ### Make empty dataframes 
  auc_train_ind <- data.frame(matrix(ncol = 3, nrow = 0))
  auc_test_ind <- data.frame(matrix(ncol = 3, nrow =0))
  colnames(auc_train_ind) <- c("iter", "types", "value")
  colnames(auc_test_ind) <- c("iter", "types", "value")
  
  auc_train_ind_working <- auc_train_ind
  auc_test_ind_working <- auc_test_ind
  
  data <- spcc(site,  #function to get the site data into new groups 
                      threshold = 10)
  
  #get the types from the model
  types <- data %>% 
       pull(new_group) %>% 
    unique() %>% 
    sort()
  
  ## If there is less than 15 individuals in a group, we do not run the AUC for that group
  model_types <- data %>% 
    group_by(new_group) %>% 
    mutate(count = n()) %>% 
    ungroup() %>% 
    filter(count > 15) %>% ### !!! If there is less than 15 individuals in a group, the AUC is not done
    pull(new_group) %>% 
    unique() %>% 
    sort()
  
  nspecies <- length(types)
  window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
  
  #Training loop 
  for (i in 1:n_it) {
    
    set.seed(i)
    print(paste0("Starting n = ", i))
    
    
    for (j in model_types){
      
      
      #### Select 20 individuals from mark j for iteration i
      df_test <- data %>% filter(new_group == j) 
      df_test <- df_test[sample(nrow(df_test), n_ind, replace = FALSE), ]
      df_train <- anti_join(data, df_test)
      # semi_join(df, leave_out_df) #check if the anti_join has worked 
      
      
      #### Make test and train configurations from sampled datasets 
      configuration_train <- Configuration(df_train$x_jitter, df_train$y_jitter, df_train$new_group) 
      configuration_test <- Configuration(df_test$x_jitter, df_test$y_jitter, df_test$new_group)
      
      ### Fit the training model
      fit_train <- ppjsdm::gibbsm(configuration_train, 
                                  window = window, 
                                  short_range = matrix(8, nspecies, nspecies), 
                                  model = "exponential",
                                  fitting_package = "glmnet",
                                  saturation = 10, 
                                  dummy_distribution = "stratified",
                                  min_dummy = 1, dummy_factor = 1e10, 
                                  max_dummy = 1e3) 
      
      print(paste0("*** Starting model = ", j))
      
      
      
      auc_train_ind_working[1, "iter"] <- i
      auc_train_ind_working[1, "types"] <- j
      
      p_train <- plot_papangelou(fit_train, #predict to training model fit 
                                 drop_type_from_configuration = TRUE, 
                                 type = j, 
                                 use_log = FALSE, 
                                 show = j, 
                                 return_papangelou = TRUE)
      
      p_train$v[is.na(p_train$v)] <- mean(p_train)
      p_train$v[p_train$v == -Inf] <- -1e10
      
      # AUC values for training fit 
      X <- subset(as.ppp(configuration_train, W = window), marks == j) 
      auc_train_ind_working[1, "value"] <- auc(X, covariate = as.function(p_train)) 
      
      auc_train_ind <- rbind(auc_train_ind, auc_train_ind_working)
      
      
      #Testing/Predicting 
      auc_test_ind_working[1, "iter"] <- i
      auc_test_ind_working[1, "types"] <- j
      
      #Using fit_train to predict to configuration_test
      p <- plot_papangelou(fit_train, 
                           configuration = configuration_train,
                           type = j,
                           drop_type_from_configuration = FALSE,
                           use_log = FALSE, 
                           return_papangelou = TRUE)
      
      X <- subset(as.ppp(configuration_test, W = window), marks == j)
      X <- ppp(x = X$x,
               y = X$y,
               window = X$window)
      auc_test_ind_working[1, "value"] <- auc(X = X, covariate = as.function(p)) #testing auc
      
      
      auc_test_ind <- rbind(auc_test_ind, auc_test_ind_working)
      
    } # closing the j loop
    
  }# closing i loop
  
  

  auc_train_ind$split <- "train"
  auc_test_ind$split <- "test"
  
  rbind(auc_train_ind, auc_test_ind)
  
} #closing function 




####################################
#### For species model 
ind_validation2 <- function(site, 
                            n_ind = 10, #number of individuals to leave out 
                            n_it = 30) #number of iterations, default to 30
{
  
  ### Make empty dataframes 
  auc_train_ind <- data.frame(matrix(ncol = 3, nrow = 0))
  auc_test_ind <- data.frame(matrix(ncol = 3, nrow =0))
  colnames(auc_train_ind) <- c("iter", "types", "value")
  colnames(auc_test_ind) <- c("iter", "types", "value")
  
  auc_train_ind_working <- auc_train_ind
  auc_test_ind_working <- auc_test_ind
  
  data <- sp(site,  #function to get the site data into new groups 
                    threshold = 10)
  
  #get the types from the model
  types <- data %>% 
    pull(species) %>% 
    unique() %>% 
    sort()
  
  ## If there is less than 15 individuals in a group, we do not run the AUC for that group
  model_types <- data %>% 
    group_by(species) %>% 
    mutate(count = n()) %>% 
    ungroup() %>% 
    filter(count > 15) %>% ### !!! If there is less than 15 individuals in a group, the AUC is not done
    pull(species) %>% 
    unique() %>% 
    sort()
  
  nspecies <- length(types)
  window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
  
  #Training loop 
  for (i in 1:n_it) {
    
    set.seed(i)
    print(paste0("Starting n = ", i))
    
    
    for (j in model_types){
      
      
      #### Select 20 individuals from mark j for iteration i
      df_test <- data %>% filter(species == j) 
      df_test <- df_test[sample(nrow(df_test), n_ind, replace = FALSE), ]
      df_train <- anti_join(data, df_test)
      # semi_join(df, leave_out_df) #check if the anti_join has worked 
      
      
      #### Make test and train configurations from sampled datasets 
      configuration_train <- Configuration(df_train$x_jitter, df_train$y_jitter, df_train$species) 
      configuration_test <- Configuration(df_test$x_jitter, df_test$y_jitter, df_test$species)
      
      ### Fit the training model
      fit_train <- ppjsdm::gibbsm(configuration_train, 
                                  window = window, 
                                  short_range = matrix(8, nspecies, nspecies), 
                                  model = "exponential",
                                  fitting_package = "glmnet",
                                  saturation = 10, 
                                  dummy_distribution = "stratified",
                                  min_dummy = 1, dummy_factor = 1e10, 
                                  max_dummy = 1e3) 
      
      print(paste0("*** Starting model = ", j))
      
      
      
      auc_train_ind_working[1, "iter"] <- i
      auc_train_ind_working[1, "types"] <- j
      
      p_train <- plot_papangelou(fit_train, #predict to training model fit 
                                 drop_type_from_configuration = TRUE, 
                                 type = j, 
                                 use_log = FALSE, 
                                 show = j, 
                                 return_papangelou = TRUE)
      
      p_train$v[is.na(p_train$v)] <- mean(p_train)
      p_train$v[p_train$v == -Inf] <- -1e10
      
      # AUC values for training fit 
      X <- subset(as.ppp(configuration_train, W = window), marks == j) 
      auc_train_ind_working[1, "value"] <- auc(X, covariate = as.function(p_train)) 
      
      auc_train_ind <- rbind(auc_train_ind, auc_train_ind_working)
      
      
      #Testing/Predicting 
      auc_test_ind_working[1, "iter"] <- i
      auc_test_ind_working[1, "types"] <- j
      
      #Using fit_train to predict to configuration_test
      p <- plot_papangelou(fit_train, 
                           configuration = configuration_train,
                           type = j,
                           drop_type_from_configuration = FALSE,
                           use_log = FALSE, 
                           return_papangelou = TRUE)
      
      X <- subset(as.ppp(configuration_test, W = window), marks == j)
      X <- ppp(x = X$x,
               y = X$y,
               window = X$window)
      auc_test_ind_working[1, "value"] <- auc(X = X, covariate = as.function(p)) #testing auc
      
      
      auc_test_ind <- rbind(auc_test_ind, auc_test_ind_working)
      
    } # closing the j loop
    
  }# closing i loop
  
  
  
  auc_train_ind$split <- "train"
  auc_test_ind$split <- "test"
  
  rbind(auc_train_ind, auc_test_ind)
  
} #closing function 



#############################################################
#### For species by diameter model 
ind_validation3 <- function(site, 
                            n_ind = 10, #number of individuals to leave out 
                            n_it = 30) #number of iterations, default to 30
{
  
  ### Make empty dataframes 
  auc_train_ind <- data.frame(matrix(ncol = 3, nrow = 0))
  auc_test_ind <- data.frame(matrix(ncol = 3, nrow =0))
  colnames(auc_train_ind) <- c("iter", "types", "value")
  colnames(auc_test_ind) <- c("iter", "types", "value")
  
  auc_train_ind_working <- auc_train_ind
  auc_test_ind_working <- auc_test_ind
  
  data <- spcc(site,  #function to get the site data into new groups 
             threshold = 10)
  
  #get the types from the model
  types <- data %>% 
    pull(new_spcc) %>% 
    unique() %>% 
    sort()
  
  ## If there is less than 15 individuals in a group, we do not run the AUC for that group
  model_types <- data %>% 
    group_by(new_spcc) %>% 
    mutate(count = n()) %>% 
    ungroup() %>% 
    filter(count > 15) %>% ### !!! If there is less than 15 individuals in a group, the AUC is not done
    pull(new_spcc) %>% 
    unique() %>% 
    sort()
  
  nspecies <- length(types)
  window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
  
  #Training loop 
  for (i in 1:n_it) {
    
    set.seed(i)
    print(paste0("Starting n = ", i))
    
    
    for (j in model_types){
      
      
      #### Select 20 individuals from mark j for iteration i
      df_test <- data %>% filter(new_spcc == j) 
      df_test <- df_test[sample(nrow(df_test), n_ind, replace = FALSE), ]
      df_train <- anti_join(data, df_test)
      # semi_join(df, leave_out_df) #check if the anti_join has worked 
      
      
      #### Make test and train configurations from sampled datasets 
      configuration_train <- Configuration(df_train$x_jitter, df_train$y_jitter, df_train$new_spcc) 
      configuration_test <- Configuration(df_test$x_jitter, df_test$y_jitter, df_test$new_spcc)
      
      ### Fit the training model
      fit_train <- ppjsdm::gibbsm(configuration_train, 
                                  window = window, 
                                  short_range = matrix(8, nspecies, nspecies), 
                                  model = "exponential",
                                  fitting_package = "glmnet",
                                  saturation = 10, 
                                  dummy_distribution = "stratified",
                                  min_dummy = 1, dummy_factor = 1e10, 
                                  max_dummy = 1e3) 
      
      print(paste0("*** Starting model = ", j))
      
      
      
      auc_train_ind_working[1, "iter"] <- i
      auc_train_ind_working[1, "types"] <- j
      
      p_train <- plot_papangelou(fit_train, #predict to training model fit 
                                 drop_type_from_configuration = TRUE, 
                                 type = j, 
                                 use_log = FALSE, 
                                 show = j, 
                                 return_papangelou = TRUE)
      
      p_train$v[is.na(p_train$v)] <- mean(p_train)
      p_train$v[p_train$v == -Inf] <- -1e10
      
      # AUC values for training fit 
      X <- subset(as.ppp(configuration_train, W = window), marks == j) 
      auc_train_ind_working[1, "value"] <- auc(X, covariate = as.function(p_train)) 
      
      auc_train_ind <- rbind(auc_train_ind, auc_train_ind_working)
      
      
      #Testing/Predicting 
      auc_test_ind_working[1, "iter"] <- i
      auc_test_ind_working[1, "types"] <- j
      
      #Using fit_train to predict to configuration_test
      p <- plot_papangelou(fit_train, 
                           configuration = configuration_train,
                           type = j,
                           drop_type_from_configuration = FALSE,
                           use_log = FALSE, 
                           return_papangelou = TRUE)
      
      X <- subset(as.ppp(configuration_test, W = window), marks == j)
      X <- ppp(x = X$x,
               y = X$y,
               window = X$window)
      auc_test_ind_working[1, "value"] <- auc(X = X, covariate = as.function(p)) #testing auc
      
      
      auc_test_ind <- rbind(auc_test_ind, auc_test_ind_working)
      
    } # closing the j loop
    
  }# closing i loop
  
  
  
  auc_train_ind$split <- "train"
  auc_test_ind$split <- "test"
  
  rbind(auc_train_ind, auc_test_ind)
  
} #closing function 





