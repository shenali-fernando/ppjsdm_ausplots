library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")


####### Dominant Eucalypt == E. delegantensis (+ 1 S TAS site with some amount)

### BenRidge
benridge <- data_c %>% filter(Site_Name == "BenRidge")
benridge %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

benridge_c <- benridge %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc == "Acacia dealbata Co-dominant", "Acacia dealbata Intermediate", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Persoonia muelleri Suppressed", "Tasmannia lanceolata Suppressed", "Pittosporum bicolor Suppressed"), "Misc Suppressed", species_cc))

benridge_c %>% count(species_cc)


configuration_benridge <- Configuration(benridge_c$Ausplot_X, benridge_c$Ausplot_Y, types = benridge_c$species_cc)
plot(configuration_benridge)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_benridge$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_benridge[h], #create the fit
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
nspecies <- length(levels(configuration_benridge$types))
short_range <- matrix(5, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_benridge <- ppjsdm::gibbsm(configuration = configuration_benridge, 
                            window = window,
                            short_range = short_range, 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)
sum_benridge <- summary(fit_benridge)


ppjsdm::box_plot(fit = fit_benridge,
                 summ = sum_benridge,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### Caveside
caveside <- data_c %>% filter(Site_Name == "Caveside")
caveside %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

caveside_c <- caveside %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc %in% c("Atherosperma moschatum Suppressed", "Pomaderris apetala Co-dominant", "Pittosporum bicolor Suppressed", "Phyllocladus aspleniifolius Suppressed")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Acacia dealbata Intermediate", "Acacia dealbata Suppressed"), "Acacia dealbata I/S", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus obliqua Co-dominant", "Eucalyptus obliqua Dominant"), "Eucalyptus obliqua C/D", species_cc))
 
caveside_c %>% count(species_cc)


configuration_caveside <- Configuration(caveside_c$Ausplot_X, caveside_c$Ausplot_Y, types = caveside_c$species_cc)
plot(configuration_benridge)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_caveside$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_caveside[h], #create the fit
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
nspecies <- length(levels(configuration_caveside$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_caveside <- ppjsdm::gibbsm(configuration = configuration_caveside, 
                               window = window,
                               short_range = short_range, 
                               model = "exponential", 
                               saturation = 10, 
                               nthreads = 4, 
                               fitting_package = "glmnet",
                               dummy_distribution = "stratified",
                               min_dummy = 1, dummy_factor = 1e10, 
                               max_dummy = 1e3)
sum_caveside <- summary(fit_caveside)


ppjsdm::box_plot(fit = fit_caveside,
                 summ = sum_caveside,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 9)


### McKenzie 
mcken <- data_c %>% filter(Site_Name == "Mackenzie")
mcken %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

mcken_c <- mcken %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc %in% c("Atherosperma moschatum Suppressed", "Eucalyptus dalrympleana Suppressed", "Monotoca glauca Suppressed")) %>% 
  filter(! species_cc %in% c("Eucalyptus dalrympleana Dominant", "Olearia argophylla Suppressed", "Phyllocladus aspleniifolius Suppressed"))

mcken_c %>% count(species_cc)


configuration_mcken <- Configuration(mcken_c$Ausplot_X, mcken_c$Ausplot_Y, types = mcken_c$species_cc)
plot(configuration_mcken)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_mcken$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_mcken[h], #create the fit
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
nspecies <- length(levels(configuration_mcken$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_mcken <- ppjsdm::gibbsm(configuration = configuration_mcken, 
                               window = window,
                               short_range = short_range, 
                               model = "exponential", 
                               saturation = 10, 
                               nthreads = 4, 
                               fitting_package = "glmnet",
                               dummy_distribution = "stratified",
                               min_dummy = 1, dummy_factor = 1e10, 
                               max_dummy = 1e3)
sum_mcken <- summary(fit_mcken)


ppjsdm::box_plot(fit = fit_mcken,
                 summ = sum_mcken,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 9)




### MtMaurice
mtmau <- data_c %>% filter(Site_Name == "MtMaurice")
mtmau %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

mtmau_c <- mtmau %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc == "Eucalyptus delegatensis ", "Eucalyptus delegatensis Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus dalrympleana Co-dominant", 
                                                "Eucalyptus dalrympleana Dominant",
                                                "Eucalyptus dalrympleana Intermediate"), "Eucalyptus dalrympleana Large", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Monotoca glauca Suppressed", 
                                                "Nothofagus cunninghamii Suppressed", 
                                                "Persoonia muelleri Suppressed",
                                                "Phyllocladus aspleniifolius Suppressed",
                                                "Pittosporum bicolor Suppressed",
                                                "Tasmannia lanceolata Suppressed"), "Misc Suppressed", species_cc))
mtmau_c %>% count(species_cc)


configuration_mtmau <- Configuration(mtmau_c$Ausplot_X, mtmau_c$Ausplot_Y, types = mtmau_c$species_cc)
plot(configuration_mtmau)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_mtmau$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_mtmau[h], #create the fit
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
nspecies <- length(levels(configuration_mtmau$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_mtmau <- ppjsdm::gibbsm(configuration = configuration_mtmau, 
                            window = window,
                            short_range = short_range, 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)
sum_mtmau <- summary(fit_mtmau)


ppjsdm::box_plot(fit = fit_mtmau,
                 summ = sum_mtmau,
                 coefficient = "alpha",
                 text_size = 9)


### MtField 
mtfield <- data_c %>% filter(Site_Name == "MtField")
df <- mtfield %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

mtfield_c <- mtfield %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc == "Pittosporum bicolor ", "Pittosporum bicolor Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc == "Leptospermum lanigerum ", "Leptospermum lanigerum Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Hakea lissosperma Suppressed",
                                                "Coprosma nitida Suppressed",
                                                "Nothofagus cunninghamii Suppressed",
                                                "Ozothamnus antennaria Suppressed"), "Misc Suppressed", species_cc)) %>% 
 filter(! species_cc  %in% c("Eucalyptus urnigera Co-dominant", "Eucalyptus urnigera Intermediate")) %>% 
  mutate(species_cc = if_else(species_cc == "Eucalyptus coccifera ", "Eucalyptus coccifera Co-dominant", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc == "Eucalyptus delegatensis ", "Eucalyptus delegatensis Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus subcrenulata Co-dominant",
                                                "Eucalyptus subcrenulata Intermediate"), "Eucalyptus subcrenulata C/I", species_cc))
  
mtfield_c %>% count(species_cc)

# df <- mtfield %>% filter(Genus_Species == "Eucalyptus delegatensis") %>% 
#   group_by(Diameter, Crown_Class) %>% count()

configuration_mtfield <- Configuration(mtfield_c$Ausplot_X, mtfield_c$Ausplot_Y, types = mtfield_c$species_cc)
plot(configuration_mtfield)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_mtfield$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_mtfield[h], #create the fit
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
nspecies <- length(levels(configuration_mtfield$types))
short_range <- matrix(10, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_mtfield <- ppjsdm::gibbsm(configuration = configuration_mtfield, 
                            window = window,
                            short_range = short_range, 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)
sum_mtfield  <- summary(fit_mtfield )


ppjsdm::box_plot(fit = fit_mtfield ,
                 summ = sum_mtfield ,
                 coefficient = "alpha",
                 involving = "Misc Suppressed", 
                 how = "one",
                 text_size = 9)



### NorthStyx (not most dominant)
norths <- data_c %>% filter(Site_Name == "NorthStyx")
norths %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

norths_c <- norths %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Eucalyptus delegatensis Dominant", "Eucalyptus delegatensis Emergent" ), "Eucalyptus delegatensis Dominant", Species_CC)) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Eucalyptus regnans Dominant", "Eucalyptus regnans Emergent" ), "Eucalyptus regnans Dominant", Species_CC)) %>% 
  filter(! Species_CC %in% c("Acacia dealbata Co-dominant", "Acacia dealbata Intermediate", "Acacia melanoxylon Intermediate")) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Acacia dealbata Suppressed",
                                                "Acacia melanoxylon Suppressed", 
                                                "Olearia argophylla Suppressed", 
                                                "Pittosporum bicolor Suppressed"), "Misc Suppressed", Species_CC)) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Eucalyptus delegatensis Co-dominant",
                                                "Eucalyptus delegatensis Dominant"), 
                              "Eucalyptus delegatensis C/D", Species_CC))

norths_c %>% count(Species_CC)

configuration_norths <- Configuration(norths_c$Ausplot_X, norths_c$Ausplot_Y, types = norths_c$Species_CC)
plot(configuration_norths)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_norths$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_norths[h], #create the fit
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
nspecies <- length(levels(configuration_norths$types))
short_range <- matrix(5, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_norths <- ppjsdm::gibbsm(configuration = configuration_norths, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_norths <- summary(fit_norths)


ppjsdm::box_plot(fit = fit_norths,
                 summ = sum_norths,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12, 
                 xmin = -10, 
                 xmax = 10)
