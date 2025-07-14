library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")


####### Dominant Eucalypt == E. fastigata 
fast <- data_c %>% filter(Genus_Species == "Eucalyptus fastigata") %>% 
  group_by(Site_Name) %>% 
  count()
fast




### Goodenia 
goodenia <- data_c %>% filter(Site_Name == "Goodenia")
goodenia %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

goodenia_c <- goodenia %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc %in% c("Doryphora sassafras Suppressed", "Elaeocarpus holopetalus Suppressed", "Pittosporum undulatum Suppressed")) %>%  #too few, too far away 
  filter(! species_cc == "Eucalyptus cypellocarpa Suppressed") %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus cypellocarpa Co-dominant", 
                                                "Eucalyptus cypellocarpa Dominant"), "Eucalyptus cypellocarpa C/D", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus sieberi Dominant",
                                                "Eucalyptus sieberi Co-dominant", 
                                                "Eucalyptus sieberi Intermediate"), "Eucalyptus sieberi Large", species_cc)) 

goodenia_c %>% count(species_cc)


configuration_goodenia <- Configuration(goodenia_c$Ausplot_X, goodenia_c$Ausplot_Y, types = goodenia_c$species_cc)
plot(configuration_goodenia)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_goodenia$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_goodenia[h], #create the fit
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
nspecies <- length(levels(configuration_goodenia$types))
short_range <- matrix(9, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_goodenia <- ppjsdm::gibbsm(configuration = configuration_goodenia, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_goodenia <- summary(fit_goodenia)


ppjsdm::box_plot(fit = fit_goodenia,
                 summ = sum_goodenia,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)


ppjsdm::box_plot(fit = fit_goodenia,
                 summ = sum_goodenia,
                 coefficient = "alpha",
                 involving = "Eucalyptus sieberi Suppressed", 
                 how = "one",
                 which = "all", 
                 text_size = 10)



### Candelo 
candelo <- data_c %>% filter(Site_Name == "Candelo")
candelo %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

candelo_c <- candelo %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus cypellocarpa Co-dominant", 
                                                "Eucalyptus cypellocarpa Dominant", 
                                                "Eucalyptus cypellocarpa Intermediate"), "Eucalyptus cypellocarpa large", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Acacia melanoxylon Suppressed",
                                                "Bedfordia arborescens Suppressed", 
                                                "Eucalyptus sieberi Suppressed", 
                                                "Pittosporum undulatum Suppressed"), "Misc Suppressed", species_cc)) %>% 
  filter(! species_cc %in% c("Eucalyptus fastigata Co-dominant", "Eucalyptus fastigata Suppressed"))

candelo_c %>% count(species_cc)


configuration_candelo <- Configuration(candelo_c$Ausplot_X, candelo_c$Ausplot_Y, types = candelo_c$species_cc)
plot(configuration_candelo)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_candelo$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_candelo[h], #create the fit
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
nspecies <- length(levels(configuration_candelo$types))
short_range <- matrix(9, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_candelo <- ppjsdm::gibbsm(configuration = configuration_candelo, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_candelo <- summary(fit_candelo)


ppjsdm::box_plot(fit = fit_candelo,
                 summ = sum_candelo,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)


### Newline 
newline <- data_c %>% filter(Site_Name == "Newline")
newline %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

newline_c <- newline %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus obliqua Co-dominant", 
                                                "Eucalyptus obliqua Dominant", 
                                                "Eucalyptus obliqua Intermediate"), "Eucalyptus obliqua large", species_cc)) %>% 
  filter(! species_cc %in% c("Acacia dealbata Suppressed", "Persoonia silvatica Suppressed"))

newline_c %>% count(species_cc)


configuration_newline <- Configuration(newline_c$Ausplot_X, newline_c$Ausplot_Y, types = newline_c$species_cc)
plot(configuration_newline)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_newline$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_newline[h], #create the fit
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
nspecies <- length(levels(configuration_newline$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_newline <- ppjsdm::gibbsm(configuration = configuration_newline, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_newline <- summary(fit_newline)


ppjsdm::box_plot(fit = fit_newline,
                 summ = sum_newline,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)




### WaratahMix
waratah <- data_c %>% filter(Site_Name == "WaratahMix")
waratah %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

waratah_c <- waratah %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus cypellocarpa Co-dominant", 
                                                "Eucalyptus cypellocarpa Dominant", 
                                                "Eucalyptus cypellocarpa Intermediate"), "Eucalyptus cypellocarpa large", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus fastigata Co-dominant", 
                                                "Eucalyptus fastigata Dominant", 
                                                "Eucalyptus fastigata Intermediate"), "Eucalyptus fastigata large", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus obliqua Co-dominant", 
                                                "Eucalyptus obliqua Dominant", 
                                                "Eucalyptus obliqua Intermediate"), "Eucalyptus obliqua large", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus radiata Co-dominant", 
                                                "Eucalyptus radiata Dominant", 
                                                "Eucalyptus radiata Intermediate"), "Eucalyptus radiata large", species_cc)) %>%
  filter(! species_cc %in% c("Eucalyptus ovata Suppressed", "Eucalyptus viminalis Intermediate", "Eucalyptus viminalis Co-dominant"))

waratah_c %>% count(species_cc)


configuration_waratah <- Configuration(waratah_c$Ausplot_X, waratah_c$Ausplot_Y, types = waratah_c$species_cc)
plot(configuration_waratah)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_waratah$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_waratah[h], #create the fit
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
nspecies <- length(levels(configuration_waratah$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_waratah <- ppjsdm::gibbsm(configuration = configuration_waratah, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_waratah <- summary(fit_waratah)
fit_waratah$type_names

ppjsdm::box_plot(fit = fit_waratah,
                 summ = sum_waratah,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)


### WogWay
wogway <- data_c %>% filter(Site_Name == "WogWay")
wogway %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

wogway_c <- wogway %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc == "Acacia dealbata Co-dominant") %>% 
  mutate(species_cc = if_else(species_cc %in% c("Lomatia myricoides Suppressed", 
                                                "Prostanthera lasianthos Suppressed", 
                                                "Acacia dealbata Suppressed"), 
                              "Misc Suppressed", species_cc)) %>% 
 mutate(species_cc = if_else(species_cc %in% c("Eucalyptus obliqua Co-dominant", 
                                                "Eucalyptus obliqua Dominant"), "Eucalyptus obliqua C/D", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus fastigata Co-dominant", 
                                                "Eucalyptus fastigata Dominant"), "Eucalyptus fastigata C/D", species_cc))


wogway_c %>% count(species_cc)

configuration_wogway <- Configuration(wogway_c$Ausplot_X, wogway_c$Ausplot_Y, types = wogway_c$species_cc)
plot(configuration_wogway)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_wogway$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_wogway[h], #create the fit
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
nspecies <- length(levels(configuration_wogway$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_wogway <- ppjsdm::gibbsm(configuration = configuration_wogway, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_wogway <- summary(fit_wogway)


ppjsdm::box_plot(fit = fit_wogway,
                 summ = sum_wogway,
                 coefficient = "alpha",
                 which = "all", 
                 involving = "Eucalyptus obliqua Suppressed", 
                 how = "one",
                 text_size = 10)

