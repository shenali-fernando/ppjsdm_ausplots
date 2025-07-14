library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")

################## WA Sites 
####### Dominant Eucalypt == E. diversicolor/jacksonii (+ 1 S TAS site with some amount)

### Carey
carey <- data_c %>% filter(Site_Name == "Carey")
carey %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

carey_c <- carey %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc == "Trymalium odoratissimum Co-dominant")
  
carey_c %>% count(species_cc)


configuration_carey <- Configuration(carey_c$Ausplot_X, carey_c$Ausplot_Y, types = carey_c$species_cc)
plot(configuration_carey)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_carey$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_carey[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_carey$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_carey <- ppjsdm::gibbsm(configuration = configuration_carey, 
                               window = window,
                               short_range = short_range, 
                               model = "exponential", 
                               saturation = 10, 
                               nthreads = 4, 
                               fitting_package = "glmnet",
                               dummy_distribution = "stratified",
                               min_dummy = 1, dummy_factor = 1e10, 
                               max_dummy = 1e3)
sum_carey <- summary(fit_carey)


ppjsdm::box_plot(fit = fit_carey,
                 summ = sum_carey,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)

### Collins 
collins <- data_c %>% filter(Site_Name == "Collins")
collins %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

collins_c <- collins %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) 

collins_c %>% count(species_cc)

configuration_collins <- Configuration(collins_c$Ausplot_X, collins_c$Ausplot_Y, types = collins_c$species_cc)
plot(configuration_collins)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_collins$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_collins[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_collins$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_collins <- ppjsdm::gibbsm(configuration = configuration_collins, 
                               window = window,
                               short_range = short_range, 
                               model = "exponential", 
                               saturation = 10, 
                               nthreads = 4, 
                               fitting_package = "glmnet",
                               dummy_distribution = "stratified",
                               min_dummy = 1, dummy_factor = 1e10, 
                               max_dummy = 1e3)
sum_collins <- summary(fit_collins)


ppjsdm::box_plot(fit = fit_collins,
                 summ = sum_collins,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)



### Clare 
clare <- data_c %>% filter(Site_Name == "Clare")
clare %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

clare_c <- clare %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc == "Allocasuarina decussata Dominant", "Allocasuarina decussata Suppressed", species_cc)) %>% #Dominant is most likely misid
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus diversicolor Dominant", 
                                                "Eucalyptus diversicolor Co-dominant"), #too few ind
                              "Eucalyptus diversicolor C/D", species_cc)) %>% 
  filter(! species_cc %in% c("Eucalyptus guilfoylei Co-dominant", "Eucalyptus guilfoylei Intermediate")) %>% #much too few ind
  filter(! species_cc == "Eucalyptus jacksonii Intermediate") #too few, no interactions    

clare_c %>% count(species_cc)


configuration_clare <- Configuration(clare_c$Ausplot_X, clare_c$Ausplot_Y, types = clare_c$species_cc)
plot(configuration_clare)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_clare$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_clare[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_clare$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_clare <- ppjsdm::gibbsm(configuration = configuration_clare, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_clare <- summary(fit_clare)


ppjsdm::box_plot(fit = fit_clare,
                 summ = sum_clare,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)



### Dawson
dawson <- data_c %>% filter(Site_Name == "Dawson")
dawson %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

dawson_c <- dawson %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc %in% c("Corymbia calophylla Suppressed", "Trymalium odoratissimum Suppressed")) #too few, too far away to group into misc

dawson_c %>% count(species_cc)


configuration_dawson <- Configuration(dawson_c$Ausplot_X, dawson_c$Ausplot_Y, types = dawson_c$species_cc)
plot(configuration_dawson)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_dawson$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_dawson[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_dawson$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_dawson <- ppjsdm::gibbsm(configuration = configuration_dawson, 
                            window = window,
                            short_range = short_range, 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)
sum_dawson <- summary(fit_dawson)


ppjsdm::box_plot(fit = fit_dawson,
                 summ = sum_dawson,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)



### Dombakup 
dom <- data_c %>% filter(Site_Name == "Dombakup")
dom %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

dom_c <- dom %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc %in% c("Corymbia calophylla Suppressed", "Agonis flexuosa Suppressed")) %>%  #too few to group into misc
  mutate(species_cc = if_else(species_cc == "Allocasuarina decussata Co-dominant", "Allocasuarina decussata Suppressed", species_cc)) %>%  #looks like misid of cc
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus diversicolor Co-dominant", 
                                                "Eucalyptus diversicolor Dominant", 
                                                "Eucalyptus diversicolor Intermediate"), 
                              "Eucalyptus diversicolor Large", species_cc)) 


dom_c %>% count(species_cc)
 
#df <- dom %>% filter(Genus_Species == "Allocasuarina decussata") %>% group_by(Crown_Class, Diameter) %>% count()

configuration_dom <- Configuration(dom_c$Ausplot_X, dom_c$Ausplot_Y, types = dom_c$species_cc)
plot(configuration_dom)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_dom$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_dom[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_dom$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_dom <- ppjsdm::gibbsm(configuration = configuration_dom, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_dom <- summary(fit_dom)


ppjsdm::box_plot(fit = fit_dom,
                 summ = sum_dom,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)




### Giants 
giants <- data_c %>% filter(Site_Name == "Giants")
giants %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

giants_c <- giants %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc %in% c("Acacia pentadenia Suppressed", "Eucalyptus diversicolor Suppressed")) %>%  #too few to group into misc
  mutate(species_cc = if_else(species_cc == "Allocasuarina decussata Intermediate", "Allocasuarina decussata Suppressed", species_cc)) %>%  #dunno what else to do with it
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus diversicolor Co-dominant", 
                                                "Eucalyptus diversicolor Dominant", 
                                                "Eucalyptus diversicolor Intermediate"), 
                              "Eucalyptus diversicolor Large", species_cc)) 


giants_c %>% count(species_cc)

#df <- giants %>% filter(Genus_Species == "Allocasuarina decussata") %>% group_by(Crown_Class, Diameter) %>% count()

configuration_giants <- Configuration(giants_c$Ausplot_X, giants_c$Ausplot_Y, types = giants_c$species_cc)
plot(configuration_giants)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_giants$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_giants[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_giants$types))
short_range <- matrix(9, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_giants <- ppjsdm::gibbsm(configuration = configuration_giants, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_giants <- summary(fit_giants)


ppjsdm::box_plot(fit = fit_giants,
                 summ = sum_giants,
                 coefficient = "alpha",
                 which = "all", 
                 involving = "Eucalyptus diversicolor Large", 
                 how = "one",
                 text_size = 12)



### Frankland 
frankland <- data_c %>% filter(Site_Name == "Frankland")
frankland %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

frankland_c <- frankland %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc == c("Trymalium odoratissimum Suppressed")) %>%  #too few to group into misc
  filter(! species_cc %in% c("Eucalyptus guilfoylei Intermediate", "Eucalyptus diversicolor Intermediate")) #too far away, too few
  


frankland_c %>% count(species_cc)

#df <- frankland %>% filter(Genus_Species == "Allocasuarina decussata") %>% group_by(Crown_Class, Diameter) %>% count()

configuration_frankland <- Configuration(frankland_c$Ausplot_X, frankland_c$Ausplot_Y, types = frankland_c$species_cc)
plot(configuration_frankland)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_frankland$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_frankland[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_frankland$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_frankland <- ppjsdm::gibbsm(configuration = configuration_frankland, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_frankland <- summary(fit_frankland)


ppjsdm::box_plot(fit = fit_frankland,
                 summ = sum_frankland,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### Sutton 
sutton <- data_c %>% filter(Site_Name == "Sutton")
sutton %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

sutton_c <- sutton %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) 

sutton_c %>% count(species_cc)

#df <- sutton %>% filter(Genus_Species == "Allocasuarina decussata") %>% group_by(Crown_Class, Diameter) %>% count()

configuration_sutton <- Configuration(sutton_c$Ausplot_X, sutton_c$Ausplot_Y, types = sutton_c$species_cc)
plot(configuration_sutton)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_sutton$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_sutton[h], #create the fit
                            window = window, 
                            model = df$model[i],
                            short_range = matrix(df$short[i]),
                            saturation = 10,
                            dummy_distribution = "stratified",
                            min_dummy = 1, max_dummy =1e3,
                            dummy_factor = 1e10,
                            nthreads = 4, 
                            fitting_package = "glmnet")
      
      fit$aic #take the aic value from the fit
    })
  }
  
  possible_short <- seq(from = 1, to = 15, length.out = 50) #possible short_range values
  possible_model <- c("square_exponential", "exponential", "square_bump", "bump") #possible models
  df <- expand.grid(short = possible_short, model = possible_model) #creating dataframe with possible values
  df$aic <- to_optimize(df) #optimisation, run the created dataframe through the function
  df$potentials <- df$model
  
  plot <- ggplot(df) + geom_point(aes(x = short, y = aic, colour = potentials)) + ggtitle(h)
  
  plotlist[[h]] <- plot
}

for (p in plotlist) {
  print(p)
}

# Set parameters
nspecies <- length(levels(configuration_sutton$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_sutton <- ppjsdm::gibbsm(configuration = configuration_sutton, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_sutton <- summary(fit_sutton)


ppjsdm::box_plot(fit = fit_sutton,
                 summ = sum_sutton,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)

### Warren
warren <- data_c %>% filter(Site_Name == "Warren")
warren %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

warren_c <- warren %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc == "Agonis flexuosa Suppressed")

warren_c %>% count(species_cc)

#df <- warren %>% filter(Genus_Species == "Allocasuarina decussata") %>% group_by(Crown_Class, Diameter) %>% count()

configuration_warren <- Configuration(warren_c$Ausplot_X, warren_c$Ausplot_Y, types = warren_c$species_cc)
plot(configuration_warren)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()


# Set parameters
nspecies <- length(levels(configuration_warren$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_warren <- ppjsdm::gibbsm(configuration = configuration_warren, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_warren <- summary(fit_warren)


ppjsdm::box_plot(fit = fit_warren,
                 summ = sum_warren,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)

