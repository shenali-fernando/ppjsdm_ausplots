library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")


####### Dominant Eucalypt == E. pilularis 
pil <- data_c %>% filter(Genus_Species == "Eucalyptus pilularis") %>% 
  group_by(Site_Name) %>% 
  count()
pil


### A-Tree
atree <- data_c %>% filter(Site_Name == "A-Tree")
atree %>% group_by(Genus_Species, Crown_Class) %>% 
  count()

atree_c <- atree %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc == "Allocasuarina torulosa ", "Allocasuarina torulosa Suppressed", species_cc)) %>% 
  filter(! species_cc %in% c("Unidentified tree Suppressed", "Trochocarpa laurina Suppressed", "Lophostemon sp. Suppressed")) %>% #too few 
 filter(! species_cc %in% c("Eucalyptus microcorys Co-dominant", "Eucalyptus saligna Co-dominant", "Eucalyptus saligna Intermediate", "Eucalyptus saligna Suppressed")) #too few, too far

atree_c %>% count(species_cc)


configuration_atree <- Configuration(atree_c$Ausplot_X, atree_c$Ausplot_Y, types = atree_c$species_cc)
plot(configuration_atree)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

plotlist <- list()

for (h in levels(configuration_atree$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_atree[h], #create the fit
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
nspecies <- length(levels(configuration_atree$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_atree<- ppjsdm::gibbsm(configuration = configuration_atree, 
                               window = window,
                               short_range = short_range, 
                               model = "exponential", 
                               saturation = 10, 
                               nthreads = 4, 
                               fitting_package = "glmnet",
                               dummy_distribution = "stratified",
                               min_dummy = 1, dummy_factor = 1e10, 
                               max_dummy = 1e3)
sum_atree <- summary(fit_atree)


ppjsdm::box_plot(fit = fit_atree,
                 summ = sum_atree,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10, 
                 xmin = -6, 
                 xmax = 6)


### BirdTree 
birdtree <- data_c %>% filter(Site_Name == "BirdTree")
b <- birdtree %>% group_by(Genus_Species, Crown_Class) %>% 
  count()
b

birdtree_c <- birdtree %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Archihodomyrtus beckleri Suppressed", "Caldcluvia paniculosa Suppressed", "Cissus hypoglauca Suppressed", "Clerodendrum floribundum Suppressed",
                                                "Eucalyptus pilularis Suppressed", "Neolitsea dealbata Suppressed", "Syzygium oleosum Suppressed", "Unidentified tree Suppressed", "Eupomatia laurina Suppressed"), 
                                                "Misc Suppressed", species_cc)) %>% 
  filter(! species_cc %in% c("Eucalyptus microcorys Co-dominant", "Eucalyptus microcorys Intermediate", "Eucalyptus microcorys Dominant", "Lophostemon sp. Co-dominant", "Syncarpia glomulifera subsp glabra Intermediate")) #too few, too far

birdtree_c %>% count(species_cc)


configuration_birdtree <- Configuration(birdtree_c$Ausplot_X, birdtree_c$Ausplot_Y, types = birdtree_c$species_cc)
plot(configuration_birdtree)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))

birdtree_c <- birdtree_c %>%
  group_by(Ausplot_X, Ausplot_Y) %>% #group by coordinate columns
  mutate(
    is_duplicated = n() > 1, #create column of TRUE/FALSE 
    #new_column_name = if_else(condition, true, false): so condition=column name, if true=fill with, if false=fill with
    x_jitter = if_else(is_duplicated, Ausplot_X + runif(n(), -0.025, 0.025), Ausplot_X), #create x_jitter column
    y_jitter = if_else(is_duplicated, Ausplot_Y + runif(n(), -0.025, 0.025), Ausplot_Y) #create y_jitter column
  ) %>%
  ungroup() 

configuration_birdtree <- Configuration(birdtree_c$x_jitter, birdtree_c$y_jitter, types = birdtree_c$species_cc)
plot(configuration_birdtree)

plotlist <- list()

for (h in levels(configuration_birdtree$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_birdtree[h], #create the fit
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
nspecies <- length(levels(configuration_birdtree$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_birdtree<- ppjsdm::gibbsm(configuration = configuration_birdtree, 
                           window = window,
                           short_range = short_range, 
                           model = "exponential", 
                           saturation = 10, 
                           nthreads = 4, 
                           fitting_package = "glmnet",
                           dummy_distribution = "stratified",
                           min_dummy = 1, dummy_factor = 1e10, 
                           max_dummy = 1e3)
sum_birdtree <- summary(fit_birdtree)


ppjsdm::box_plot(fit = fit_birdtree,
                 summ = sum_birdtree,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)




### BlackBull 
blackbull <- data_c %>% filter(Site_Name == "BlackBull")
b <- blackbull %>% group_by(Genus_Species, Crown_Class) %>% 
  count()
b

blackbull_c <- blackbull %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Acacia maidenii Suppressed", "Callitris sp. Suppressed", "Elaeocarpus reticulatus Suppressed", "Unidentified tree Suppressed",
                                                "Neolitsea dealbata Suppressed", "Persoonia conjuncta Suppressed", "Schizomeria ovata Suppressed"), 
                              "Misc Suppressed", species_cc)) %>% 
  filter(! species_cc %in% c("Allocasuarina torulosa Co-dominant", "Allocasuarina torulosa Intermediate", "Syncarpia glomulifera subsp glomulifera Co-dominant")) %>% #too few, too far
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus andrewsii Suppressed", "Eucalyptus saligna Suppressed"), 
                              "Euc Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus andrewsii Co-dominant", "Eucalyptus andrewsii Dominant", "Eucalyptus andrewsii Intermediate", 
                                                "Eucalyptus saligna Co-dominant", "Eucalyptus saligna Intermediate"), 
                              "Euc Large", species_cc))

blackbull_c %>% count(species_cc)


configuration_blackbull <- Configuration(blackbull_c$Ausplot_X, blackbull_c$Ausplot_Y, types = blackbull_c$species_cc)
plot(configuration_blackbull)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_blackbull$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_blackbull[h], #create the fit
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
nspecies <- length(levels(configuration_blackbull$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_blackbull <- ppjsdm::gibbsm(configuration = configuration_blackbull, 
                              window = window,
                              short_range = short_range, 
                              model = "exponential", 
                              saturation = 10, 
                              nthreads = 4, 
                              fitting_package = "glmnet",
                              dummy_distribution = "stratified",
                              min_dummy = 1, dummy_factor = 1e10, 
                              max_dummy = 1e3)
sum_blackbull <- summary(fit_blackbull)


ppjsdm::box_plot(fit = fit_blackbull,
                 summ = sum_blackbull,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)


ppjsdm::box_plot(fit = fit_blackbull,
                 summ = sum_blackbull,
                 coefficient = "alpha",
                 which = "all", 
                 involving = "Euc Large", 
                 how = "one",
                 text_size = 10)


### Bruxner - 145 unidentified suppressed
bruxner <- data_c %>% filter(Site_Name == "Bruxner")
b <- bruxner %>% group_by(Genus_Species, Crown_Class) %>% 
  count()
b

bruxner_c <- bruxner %>% 
  group_by(Genus_Species, Crown_Class) %>% 
  mutate(species_cc_count = n()) %>% 
  ungroup() %>%
  mutate(species_cc = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Acacia maidenii Suppressed", "Alphitonia excelsa Suppressed", "Callicoma serratifolia Suppressed", "Cryptocarya bidwillii Suppressed", 
                                                "Cryptocarya rigida Suppressed"), 
                              "Misc Suppressed", species_cc)) %>% 
  filter(! species_cc == c("Allocasuarina torulosa Dominant")) %>% #too few, too far
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus andrewsii Suppressed", "Eucalyptus saligna Suppressed"), 
                              "Euc Suppressed", species_cc)) %>% 
  mutate(species_cc = if_else(species_cc %in% c("Eucalyptus andrewsii Co-dominant", "Eucalyptus andrewsii Dominant", "Eucalyptus andrewsii Intermediate", 
                                                "Eucalyptus saligna Co-dominant", "Eucalyptus saligna Intermediate"), 
                              "Euc Large", species_cc))

bruxner_c %>% count(species_cc)


configuration_blackbull <- Configuration(blackbull_c$Ausplot_X, blackbull_c$Ausplot_Y, types = blackbull_c$species_cc)
plot(configuration_blackbull)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration_blackbull$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration_blackbull[h], #create the fit
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
nspecies <- length(levels(configuration_blackbull$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_blackbull <- ppjsdm::gibbsm(configuration = configuration_blackbull, 
                                window = window,
                                short_range = short_range, 
                                model = "exponential", 
                                saturation = 10, 
                                nthreads = 4, 
                                fitting_package = "glmnet",
                                dummy_distribution = "stratified",
                                min_dummy = 1, dummy_factor = 1e10, 
                                max_dummy = 1e3)
sum_blackbull <- summary(fit_blackbull)


ppjsdm::box_plot(fit = fit_blackbull,
                 summ = sum_blackbull,
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 10)


### Lorne 



### MinesRd 



### Tinebank 