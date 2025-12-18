library(ppjsdm)
library(dplyr)


#Understand what effect a coefficient has for a fitted ppjsdm model

effect_size <- function(configuration, #full ppjsdm::Configuration for the model
                        fit #full fitted ppjsdm model
                        ){
set.seed(12345) #set seed for reproducibility of glm

full_config <- configuration
all_types <- levels(unique(full_config$types))

estimates <- fit$coefficients$alpha[[1]]

#add other columns to df
df <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(df) <- c("from", "to", "alpha", "mean_effect", "quan_0", "quan_25", "quan_50", "quan_75", "quan_100")


for(type_i in all_types){ #open first loop 
  
  types <- setdiff(all_types, type_i)
  
  #get the locations of species i 
  locs <- configuration[type_i]
  
  #compute full papangelou for type i 
  full <- compute_papangelou(fit, 
                             x = locs$x,
                             y = locs$y,
                             configuration = configuration, 
                             type = type_i, 
                             drop_type_from_configuration = FALSE,
                             use_log = TRUE)

  for(type_j in types){ #open second loop 
  
  #construct dropped type configuration
  drop_types <- setdiff(all_types, type_j) 
  drop_config <- full_config[drop_types] #type_i needed in this config
  #plot(drop_config)
  
  #get papangelou intensities for dropped fit 
  drop <- compute_papangelou(fit = fit,
                             configuration = drop_config,
                             x = locs$x,
                             y = locs$y,
                             type = type_i, 
                             drop_type_from_configuration = FALSE, 
                             nthreads = 3)
  
  
#compute the summary stats we want to output in the final df - put into tmp
tmp <- data.frame(
    from = type_j, 
    to = type_i, 
    alpha = estimates[type_i, type_j],  
    mean_effect = mean(full/drop),
    quan_0 = quantile(full/drop)[1],
    quan_25 = quantile(full/drop)[2],
    quan_50 = quantile(full/drop)[3],
    quan_75 = quantile(full/drop)[4],
    quan_100 = quantile(full/drop)[5])

df <- rbind(df, tmp)

} #close second loop 
  
} #close first loop 


for(type_k in all_types){ #open third loop

  locs_k <- configuration[type_k] #extract locs for third loop
  
  full <- compute_papangelou(fit, 
                             x = locs_k$x,
                             y = locs_k$y,
                             configuration = configuration, 
                             type = type_k, 
                             drop_type_from_configuration = FALSE,
                             use_log = TRUE, 
                             nthreads = 3)
  
  drop <- compute_papangelou(fit = fit,
                             configuration = configuration,
                             x = locs_k$x,
                             y = locs_k$y,
                             type = type_k, 
                             drop_type_from_configuration = TRUE, 
                             nthreads = 3, 
                             use_log = TRUE)
  

tmp <- data.frame(
    from = type_k, 
    to = type_k, 
    alpha = estimates[type_k, type_k], 
    mean_effect = mean(full/drop),
    quan_0 = quantile(full/drop)[1],
    quan_25 = quantile(full/drop)[2],
    quan_50 = quantile(full/drop)[3],
    quan_75 = quantile(full/drop)[4],
    quan_100 = quantile(full/drop)[5])

df <- rbind(df, tmp)
  
  } #close third loop 

#get rid of duplicated rows (i.e. lower half of matrix)
df <- df %>%
  group_by(alpha) %>%
  slice_max(order_by = mean_effect, n = 1) %>% #keep the row with the max mean effect calculated
  ungroup()


  return(df)

} #close function



#Test Function for site: ZigZag
e <- effect_size(zigzag_config,
                 zigzag_fit)
