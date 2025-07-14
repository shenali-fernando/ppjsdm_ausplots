library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")


####### Dominant Eucalypt == E. obliqua (+ 3 SE VIC sites with some amount)

### BlackRiver 
black <- data_c %>% filter(Site_Name == "BlackRiver")
black %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

black_c <- black %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc == "Monotoca glauca Dominant", "Monotoca glauca Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc == "Eucalyptus obliqua ", "Eucalyptus obliqua Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Acacia melanoxylon Suppressed", "Acacia mucronata Suppressed"), "Acacia", species_cc))

black_c %>% count(species_cc)


configuration_black <- Configuration(black_c$Ausplot_X, black_c$Ausplot_Y, types = black_c$species_cc)
plot(configuration_black)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_black$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_black[h], #create the fit
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
nspecies <- length(levels(configuration_black$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_black <- ppjsdm::gibbsm(configuration = configuration_black, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_black <- summary(fit_black)


ppjsdm::box_plot(fit = fit_black,
                 summ = sum_black,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### BondTier 
bond <- data_c %>% filter(Site_Name == "BondTier")
bond %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

bond_c <- bond %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc == "Eucalyptus obliqua ", "Eucalyptus obliqua Suppressed", species_cc)) %>% 
  filter(! species_cc %in% c("Acacia melanoxylon Intermediate", "Leptospermum scoparium Suppressed", "Monotoca glauca Suppressed", "Nothofagus cunninghamii Intermediate"))

bond_c %>% count(species_cc)


configuration_bond <- Configuration(bond_c$Ausplot_X, bond_c$Ausplot_Y, types = bond_c$species_cc)
plot(configuration_bond)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_bond$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_bond[h], #create the fit
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
nspecies <- length(levels(configuration_bond$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_bond <- ppjsdm::gibbsm(configuration = configuration_bond, 
                            window = window,
                            short_range = short_range, 
                            model = "exponential", 
                            saturation = 10, 
                            nthreads = 4, 
                            fitting_package = "glmnet",
                            dummy_distribution = "stratified",
                            min_dummy = 1, dummy_factor = 1e10, 
                            max_dummy = 1e3)
sum_bond <- summary(fit_bond)


ppjsdm::box_plot(fit = fit_bond,
                 summ = sum_bond,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)



### Supersite 
#Supersite is a problem: It is 1.6ha and has a cleared area, therefore we can subset it to 1ha to get rid of the cleared area 
supersite <- data_c %>% filter(Site_Name == "Supersite") 
configuration<- ppjsdm::Configuration(supersite$Ausplot_X, supersite$Ausplot_Y, supersite$Genus_Species)
plot(configuration)
#1ha of supersite that was not cleared
supersite_c <- supersite %>% filter(Ausplot_X <= 100)
configuration<- ppjsdm::Configuration(supersite$Ausplot_X, supersite$Ausplot_Y, supersite$Genus_Species)
plot(configuration)


supersite_c %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

supersite_c <- supersite_c %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus obliqua Emergent", "Eucalyptus obliqua Dominant"), "Eucalyptus obliqua Dominant", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Pomaderris apetala Dominant", "Pomaderris apetala Intermediate"), "Pomaderris apetala Intermediate", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc == "Nothofagus cunninghamii ", "Nothofagus cunninghamii Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Nothofagus cunninghamii Co-dominant", "Nothofagus cunninghamii Intermediate"), "Nothofagus cunninghamii C/I", species_cc)) %>%
  mutate(species_cc = if_else(species_cc %in% c("Phyllocladus aspleniifolius Intermediate", "Phyllocladus aspleniifolius Suppressed"), "Phyllocladus aspleniifolius I/S", species_cc)) %>% 
  filter(! species_cc == "Pittosporum bicolor Intermediate") %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucryphia lucida Suppressed", "Meleleuca ericifolia Suppressed", "Acacia melanoxylon Suppressed", "Monotoca glauca Suppressed"), "Misc Suppressed", species_cc))

  
df <- supersite_c %>% count(species_cc)


configuration_supersite <- Configuration(supersite_c$Ausplot_X, supersite_c$Ausplot_Y, types = supersite_c$species_cc)
plot(configuration_supersite)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_supersite$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_supersite[h], #create the fit
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
nspecies <- length(levels(configuration_supersite$types))
short_range <- matrix(10, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_supersite <- ppjsdm::gibbsm(configuration = configuration_supersite, 
                           window = window,
                           short_range = short_range, 
                           model = "exponential", 
                           saturation = 10, 
                           nthreads = 4, 
                           fitting_package = "glmnet",
                           dummy_distribution = "stratified",
                           min_dummy = 1, dummy_factor = 1e10, 
                           max_dummy = 1e3)
sum_supersite <- summary(fit_supersite)


ppjsdm::box_plot(fit = fit_supersite,
                 summ = sum_supersite,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)




### ZigZag
zigzag <- data_c %>% filter(Site_Name == "ZigZag")
zigzag %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

zigzag_c <- zigzag %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus obliqua Emergent", "Eucalyptus obliqua Dominant"), "Eucalyptus obliqua Dominant", species_cc)) %>% 
  filter(! species_cc %in% c("Acacia dealbata Suppressed", "Acacia verticillata Suppressed", "Pittosporum bicolor Suppressed", "Pomaderris apetala Intermediate", "Pomaderris apetala Co-dominant")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus regnans Co-dominant", "Eucalyptus regnans Intermediate"), "Eucalyptus regnans C/I", species_cc))

zigzag_c %>% count(species_cc)


configuration_zigzag <- Configuration(zigzag_c$Ausplot_X, zigzag_c$Ausplot_Y, types = zigzag_c$species_cc)
plot(configuration_zigzag)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_zigzag$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_zigzag[h], #create the fit
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
nspecies <- length(levels(configuration_zigzag$types))
short_range <- matrix(6, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_zigzag <- ppjsdm::gibbsm(configuration = configuration_zigzag, 
                                window = window,
                                short_range = short_range, 
                                model = "exponential", 
                                saturation = 10, 
                                nthreads = 4, 
                                fitting_package = "glmnet",
                                dummy_distribution = "stratified",
                                min_dummy = 1, dummy_factor = 1e10, 
                                max_dummy = 1e3)
sum_zigzag <- summary(fit_zigzag)


ppjsdm::box_plot(fit = fit_zigzag,
                 summ = sum_zigzag,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)




### Candelo (NSW)
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





############## Not dominant but present 

#### Weld
weld <- data_c %>% filter(Site_Name == "Weld")
weld %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

weld_c <- weld %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Anopterus glandulosus Suppressed",
                                                "Atherosperma moschatum Suppressed", 
                                                "Pittosporum bicolor Suppressed", 
                                                "Leptospermum spp Suppressed"), 
                              "Misc Suppressed", Species_CC)) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Eucalyptus obliqua Co-dominant","Eucalyptus obliqua Dominant"), "Eucalyptus obliqua Co-dominant", Species_CC)) %>% 
  filter(! Species_CC == "Nothofagus cunninghamii Co-dominant") %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Eucalyptus regnans Co-dominant", "Eucalyptus regnans Intermediate"), "Eucalyptus regnans C/I", Species_CC))

weld_c %>% count(Species_CC)

configuration_weld <- Configuration(weld_c$Ausplot_X, weld_c$Ausplot_Y, types = weld_c$Species_CC)
plot(configuration_weld)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_weld$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_weld[h], #create the fit
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
nspecies <- length(levels(configuration_weld$types))
short_range <- matrix(10, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_weld <- ppjsdm::gibbsm(configuration = configuration_weld, 
                           window = window,
                           short_range = short_range, 
                           model = "exponential", 
                           saturation = 10, 
                           nthreads = 4, 
                           fitting_package = "glmnet",
                           dummy_distribution = "stratified",
                           min_dummy = 1, dummy_factor = 1e10, 
                           max_dummy = 1e3)
sum_weld <- summary(fit_weld)


ppjsdm::box_plot(fit = fit_weld,
                 summ = sum_weld,
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
short_range <- matrix(4, nspecies, nspecies)
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


ppjsdm::box_plot(fit = fit_waratah,
                 summ = sum_waratah,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)
