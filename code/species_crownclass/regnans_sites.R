library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")


####### Dominant Eucalypt == E. regnans 
reg <- data_c %>% filter(Genus_Species == "Eucalyptus regnans") %>% 
  group_by(Site_Name) %>% 
  count()
reg



### ANU101 
anu101 <- data_c %>% filter(Site_Name == "ANU101")
anu101 %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

anu101_c <- anu101 %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc == "Pomaderris aspera Intermediate") %>% 
  filter(! species_cc == "Unidentified tree Suppressed") %>% 
 mutate(species_cc = if_else(species_cc_count <5, "Misc Suppressed", species_cc)) %>% 
   filter(! species_cc == "Acacia dealbata Intermediate")

anu101_c %>% count(species_cc)


configuration_anu101 <- Configuration(anu101_c$Ausplot_X, anu101_c$Ausplot_Y, types = anu101_c$species_cc)
plot(configuration_anu101)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_anu101$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_anu101[h], #create the fit
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
nspecies <- length(levels(configuration_anu101$types))
short_range <- matrix(9, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_anu101 <- ppjsdm::gibbsm(configuration = configuration_anu101, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_anu101 <- summary(fit_anu101)


ppjsdm::box_plot(fit = fit_anu101,
                 summ = sum_anu101,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### ANU363 
anu363 <- data_c %>% filter(Site_Name == "ANU363")
anu363 %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

anu363_c <- anu363 %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc_count <6, "Misc Suppressed", species_cc))

anu363_c %>% count(species_cc)


configuration_anu363 <- Configuration(anu363_c$Ausplot_X, anu363_c$Ausplot_Y, types = anu363_c$species_cc)
plot(configuration_anu363)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_anu363$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_anu363[h], #create the fit
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
nspecies <- length(levels(configuration_anu363$types))
short_range <- matrix(9, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_anu363 <- ppjsdm::gibbsm(configuration = configuration_anu363, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_anu363 <- summary(fit_anu363)


ppjsdm::box_plot(fit = fit_anu363,
                 summ = sum_anu363,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### ANU589 
anu589 <- data_c %>% filter(Site_Name == "ANU589")
anu589 %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

anu589_c <- anu589 %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! species_cc == "Tasmannia lanceolata Suppressed")

anu589_c %>% count(species_cc)


configuration_anu589 <- Configuration(anu589_c$Ausplot_X, anu589_c$Ausplot_Y, types = anu589_c$species_cc)
plot(configuration_anu589)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_anu589$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_anu589[h], #create the fit
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
nspecies <- length(levels(configuration_anu589$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_anu589 <- ppjsdm::gibbsm(configuration = configuration_anu589, 
                             window = window,
                             short_range = short_range, 
                             model = "exponential", 
                             saturation = 10, 
                             nthreads = 4, 
                             fitting_package = "glmnet",
                             dummy_distribution = "stratified",
                             min_dummy = 1, dummy_factor = 1e10, 
                             max_dummy = 1e3)
sum_anu589 <- summary(fit_anu589)


ppjsdm::box_plot(fit = fit_anu589,
                 summ = sum_anu589,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)


### Ada Tree 
ada <- data_c %>% filter(Site_Name == "Ada Tree")
ada %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

ada_c <- ada %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC_count <10, "Misc Suppressed", Species_CC))

ada_c %>% count(Species_CC)


configuration_ada <- Configuration(ada_c$Ausplot_X, ada_c$Ausplot_Y, types = ada_c$Species_CC)
plot(configuration_ada)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

ada_c <- ada_c %>%
  group_by(Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) %>%
  ungroup() 



configuration_ada <- Configuration(ada_c$x_jitter, ada_c$y_jitter, types = ada_c$Species_CC)
plot(configuration_ada)

plotlist <- list()

for (h in levels(configuration_ada$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_ada[h], #create the fit
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
nspecies <- length(levels(configuration_ada$types))
short_range <- matrix(4.5, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_ada <- ppjsdm::gibbsm(configuration = configuration_ada, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_ada <- summary(fit_ada)

ppjsdm::box_plot(fit = fit_ada,
                 summ = sum_ada,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)




### Weeaproinah 
wee <- data_c %>% filter(Site_Name == "Weeaproinah")
wee %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

wee_c <- wee %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC =="Nematolepis squamea subsp squamea Suppressed", "Nematolepis squamea Suppressed", Species_CC))

wee_c %>% count(Species_CC)

wee_c <- wee_c %>% 
  mutate(Species_CC = if_else(Species_CC_count <10, "Misc", Species_CC)) %>% 
  filter(!Species_CC == "Misc")

configuration_wee <- Configuration(wee_c$Ausplot_X, wee_c$Ausplot_Y, types = wee_c$Species_CC)
plot(configuration_wee)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_wee$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_wee[h], #create the fit
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
nspecies <- length(levels(configuration_wee$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_wee <- ppjsdm::gibbsm(configuration = configuration_wee, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_wee <- summary(fit_wee)


ppjsdm::box_plot(fit = fit_wee,
                 summ = sum_wee,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)



### Turtons 
turtons <- data_c %>% filter(Site_Name == "Turtons")
turtons %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

turtons_c <- turtons %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC_count <=4, "Misc Suppressed", Species_CC))

turtons_c %>% count(Species_CC)

configuration_turtons <- Configuration(turtons_c$Ausplot_X, turtons_c$Ausplot_Y, types = turtons_c$Species_CC)
plot(configuration_turtons)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_turtons$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_turtons[h], #create the fit
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
nspecies <- length(levels(configuration_turtons$types))
short_range <- matrix(4, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_turtons <- ppjsdm::gibbsm(configuration = configuration_turtons, 
                          window = window,
                          short_range = short_range, 
                          model = "exponential", 
                          saturation = 10, 
                          nthreads = 4, 
                          fitting_package = "glmnet",
                          dummy_distribution = "stratified",
                          min_dummy = 1, dummy_factor = 1e10, 
                          max_dummy = 1e3)
sum_turtons <- summary(fit_turtons)


ppjsdm::box_plot(fit = fit_turtons,
                 summ = sum_turtons,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)






### Lardner 
lardner <- data_c %>% filter(Site_Name == "Lardner")
lardner2 <- lardner %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

lardner_c <- lardner %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Coprosma quadrifida Suppressed", "Pittosporum bicolor Suppressed", "Pomaderris aspera Suppressed"), "Misc", Species_CC)) 


lardner_c %>% count(Species_CC)

configuration_lardner <- Configuration(lardner_c$Ausplot_X, lardner_c$Ausplot_Y, types = lardner_c$Species_CC)
plot(configuration_lardner)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_lardner$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_lardner[h], #create the fit
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
nspecies <- length(levels(configuration_lardner$types))
short_range <- matrix(4, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_lardner <- ppjsdm::gibbsm(configuration = configuration_lardner, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_lardner <- summary(fit_lardner)


ppjsdm::box_plot(fit = fit_lardner,
                 summ = sum_lardner,
                 coefficient = "alpha",
                 which = "all", 
                 involving = c("Acacia melanoxylon Suppressed", 
                               "Eucalyptus regnans Co-dominant",
                               "Eucalyptus regnans Dominant",                   
                               "Eucalyptus regnans Intermediate",                
                               "Eucalyptus regnans Suppressed",                    
                               "Hedycarya angustifolia Suppressed",
                               "Nematolepis squamea subsp squamea Suppressed", 
                               "Olearia argophylla Suppressed"),
                 text_size = 12)

plot(configuration[c("Eucalyptus regnans Co-dominant",
                     "Eucalyptus regnans Intermediate",                
                     "Eucalyptus regnans Suppressed")])



### HardyCreek 
hardy <- data_c %>% filter(Site_Name == "HardyCreek")
hardy %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

hardy_c <- hardy %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  filter(! Species_CC == "Eucalyptus nitens Co-dominant") %>% 
  mutate(Species_CC = if_else(Species_CC_count <=4, "Misc Suppressed", Species_CC))

hardy_c %>% count(Species_CC)

configuration_hardy <- Configuration(hardy_c$Ausplot_X, hardy_c$Ausplot_Y, types = hardy_c$Species_CC)
plot(configuration_hardy)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_hardy$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_hardy[h], #create the fit
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
nspecies <- length(levels(configuration_hardy$types))
short_range <- matrix(5, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_hardy <- ppjsdm::gibbsm(configuration = configuration_hardy, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_hardy <- summary(fit_hardy)


ppjsdm::box_plot(fit = fit_hardy,
                 summ = sum_hardy,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)




### NorthStyx
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
  mutate(Species_CC = if_else(Species_CC %in% c("Acacia dealbata Suppressed", "Acacia melanoxylon Suppressed", "Olearia argophylla Suppressed", "Pittosporum bicolor Suppressed"), "Misc Suppressed", Species_CC))

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
                 text_size = 12)


### Weld
weld <- data_c %>% filter(Site_Name == "Weld")
weld2 <- weld %>% group_by(Genus_Species, Crown_Class, Diameter) %>% 
  count() %>% 
  filter(Genus_Species == "Eucalyptus regnans")

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
  mutate(Species_CC = if_else(Species_CC %in% c("Eucalyptus regnans Co-dominant", "Eucalyptus regnans Intermediate"), "Eucalyptus regnans C/I", Species_CC)) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Eucalyptus obliqua Co-dominant", "Eucalyptus obliqua Intermediate"), 
                              "Eucalyptus obliqua C/I", Species_CC))

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



### Make df of coefficients from all fits 
#These are the fit names: 
#fit_anu101, fit_anu363, fit_anu589, fit_ada, fit_wee, 
#fit_turtons, fit_lardner, fit_hardy, fit_norths, fit_weld

estimates <- fit_weld$coefficients$alpha[[1]]
unique_names <- colnames(estimates)
df <- as.data.frame(expand.grid(from = rownames(estimates), to = colnames(estimates)))
#insert coefficient values
df$value <- sapply(seq_len(nrow(df)), function(i) {
  val <- estimates[df$from[i], df$to[i]]
  if(length(val) == 0) {
    val <- estimates[df$to[i], df$from[i]]
  }
  val
})

#insert low CI value 
df$lo <- sapply(seq_len(nrow(df)), function(i) { # Get the lower-bound of the CIs
  val <- sum_weld$lo$alpha[[1]][df$from[i], df$to[i]]
  if(length(val) == 0) {
    val <- sum$lo$alpha[[1]][df$to[i], df$from[i]]
  }
  val
})

#insert high CI value 
df$hi <- sapply(seq_len(nrow(df)), function(i) { # Get the lower-bound of the CIs
  val <- sum_weld$hi$alpha[[1]][df$from[i], df$to[i]]
  if(length(val) == 0) {
    val <- sum$hi$alpha[[1]][df$to[i], df$from[i]]
  }
  val
})

df <- df %>% mutate(Site = "Weld")

large_df <- rbind(large_df, df)
unique(large_df$Site)

#Ok, have the large df with alpha coef, lo/hi CI for each site for each group interaction
large_df <- large_df %>% 
  mutate(Genus_from = from) %>% 
  mutate(Genus = sub(" .*", "", Genus))