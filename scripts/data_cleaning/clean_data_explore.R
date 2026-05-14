library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)


#General data cleaning
ausplots <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/data/ausplots.csv", 
                     stringsAsFactors = FALSE, 
                     allowEscapes = TRUE, 
                     flush = TRUE)

# Remove individuals with NA values for coordinates
data <- ausplots[!is.na(ausplots$Ausplot_X) & !is.na(ausplots$Ausplot_Y), ]

# Filter only trees which are alive depending on what we're doing
data <- data[data$Tree_Condition == "Alive", ]

#We only want the first stem - get rid of Tree IDs that end in b, c, d,e, f  
data_c <- data %>% filter(!str_detect(Tree_ID, "(a|b|c|d|e|f)$"))
data_c %>% count(Tree_ID) #check if its done       

#save out dataframe for easy access 
write.csv(data_c, "C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/data/data_cleaned.csv")

# Some rows have a crown class value! 
cc <- data %>% filter(!Crown_Class == "")
cc %>% count(Site_Name)
#Sites without a crown class = Flowerdale, Dip, Bird


#Supersite is a problem: It is 1.6ha and has a cleared area, therefore we can subset it to 1ha to get rid of the cleared area 
#Look at the fit when we use all the supersite data
supersite_full <- data_cleaned %>% filter(Site_Name == "Supersite") 
supersite_full <- supersite_full %>% 
  mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
  mutate(new_cc = paste(species, Crown_Class)) %>% 
  group_by(new_cc) %>% 
  mutate(observation_count = n()) %>% 
  filter(!(str_starts(new_cc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_cc, "Eucalyptus") & observation_count < 5))
configuration<- ppjsdm::Configuration(supersite_full$Ausplot_X, supersite_full$Ausplot_Y, supersite_full$new_cc)
plot(configuration)
window <- ppjsdm::Rectangle_window(c(0, 150), c(0,100))
nspecies <- length(levels(configuration$types))
supersite_full_fit <- ppjsdm::gibbsm(configuration = configuration, 
                                     window = window,
                                     short_range = matrix(8, nspecies, nspecies), 
                                     model = "exponential", 
                                     saturation = 10, 
                                     nthreads = 4, 
                                     fitting_package = "glm",
                                     dummy_distribution = "stratified",
                                     min_dummy = 1, dummy_factor = 1e10, 
                                     max_dummy = 1e3)
supersite_full_sum <- summary(supersite_full_fit)


#1ha of supersite that was not cleared - look at the fit of the 1ha not cleared 
supersite <- supersite_full %>% filter(Ausplot_X <= 100)
supersite <- supersite %>% 
  mutate(species = if_else(str_starts(Genus_Species, regex("^Eucalyptus|^Corymbia"), negate = TRUE), "Non-euc", "Eucalyptus")) %>% 
  mutate(new_cc = paste(species, Crown_Class)) %>% 
  group_by(new_cc) %>% 
  mutate(observation_count = n()) %>% 
  filter(!(str_starts(new_cc, "Non-euc") & observation_count < 10)) %>% 
  filter(!(str_starts(new_cc, "Eucalyptus") & observation_count < 5))
configuration<- ppjsdm::Configuration(supersite$Ausplot_X, supersite$Ausplot_Y, supersite$new_cc)
plot(configuration)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))
supersite_fit <- ppjsdm::gibbsm(configuration = configuration, 
                                     window = window,
                                     short_range = matrix(8, nspecies, nspecies), 
                                     model = "exponential", 
                                     saturation = 10, 
                                     nthreads = 4, 
                                     fitting_package = "glm",
                                     dummy_distribution = "stratified",
                                     min_dummy = 1, dummy_factor = 1e10, 
                                     max_dummy = 1e3)
supersite_sum <- summary(supersite_fit)

#Compare the two fites in the box plot - very similar, uncertainity perhaps changed but not by much

ppjsdm::box_plot("Full" = supersite_full_fit, "Sub" = supersite_fit,
                 summ = list(supersite_full_sum, supersite_sum),
                 coefficient = "alpha",
                 which = "all", 
                 text_size = 12)



### Due to not adding much and having an area that has ongoing disturbance - we subset to 1ha
data_cleaned <- data_cleaned %>% 
  filter(Ausplot_X <= 100)

#save out again
write.csv(data_cleaned, "C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/data/data_cleaned.csv")



### Mapping 
library(ozmaps)
library(ggplot2)
ggplot() +
  geom_sf(data = ozmap_data("country"), fill = "gray90", color = "white") +
  geom_point(data = data_cleaned, aes(x = Longitude, y = Latitude, colour = Site_Name),size = 3) + 
  xlim(144, 150) + 
  ylim(-45, -40)

#Check the species by crown class divisions of some plots
supersite2 <- supersite %>% group_by(Genus_Species, Crown_Class) %>% 
  count()
supersite2

df <- caveside %>% filter(Genus_Species == "Nothofagus cunninghamii") %>% 
  group_by(Crown_Class, Diameter) %>% 
  count()

#Looking at different sites
tinebank <- data %>% filter(Site_Name == "Tinebank")
tinebank2 <- tinebank %>% group_by(Genus_Species, Crown_Class) %>% 
  count()


#Bruxner
a <- data %>% filter(Site_Name == "Tinebank") %>% group_by(Genus_Species) %>% 
  count()

e <- data_c %>% filter(Site_Name == "BondTier")
bond_d <- bond %>% filter(Crown_Class == "Dominant") 
bond_c <- bond %>% filter(Crown_Class == "Co-dominant") 
bond_i <- bond %>% filter(Crown_Class == "Intermediate") 
bond_s <- bond %>% filter(Crown_Class == "Suppressed") 
bond_s %>% count(Diameter)
par(mfrow = c(2, 2))
d <- hist(bond_d$Diameter)
c <- hist(bond_c$Diameter)
i <- hist(bond_i$Diameter)
s <- hist(bond_s$Diameter)


dip <- data_c %>% filter(Site_Name == "Dip")
hist(dip$Diameter, xlim = c(0, 100))

#Weeaproinah site - has a good subdivision of species into crown classes 
wee <- data %>% filter(Site_Name == "Weeaproinah")
wee2 <- wee %>% group_by(Genus_Species, Crown_Class) %>% 
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

configuration <- Configuration(wee_c$Ausplot_X, wee_c$Ausplot_Y, types = wee_c$Species_CC)
plot(configuration)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))


plotlist <- list()

for (h in levels(configuration$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration[h], #create the fit
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
nspecies <- length(levels(configuration$types))
short_range <- matrix(8, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_wee <- ppjsdm::gibbsm(configuration = configuration, 
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


ppjsdm::chord_diagram_plot(fit_wee,
                           summ = sum_wee,
                           coefficient = "alpha1",
                           only_statistically_significant = FALSE, # Only show stat. significant?
                           include_self = TRUE, 
                           outward_facing_names = TRUE)




### Ada Tree
ada <- data %>% filter(Site_Name == "Ada Tree")
ada2 <- ada %>% group_by(Genus_Species, Crown_Class) %>% 
  count()


ada_c <- ada %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC_count <10, "Misc", Species_CC))

ada_c %>% count(Species_CC)


configuration <- Configuration(ada_c$Ausplot_X, ada_c$Ausplot_Y, types = ada_c$Species_CC)
plot(configuration)
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


configuration <- Configuration(ada_c$x_jitter, ada_c$y_jitter, types = ada_c$Species_CC)
plot(configuration)

plotlist <- list()

for (h in levels(configuration$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration[h], #create the fit
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
nspecies <- length(levels(configuration$types))
short_range <- matrix(4, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_ada <- ppjsdm::gibbsm(configuration = configuration, 
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




#### Lardner
lardner <- data %>% filter(Site_Name == "Lardner")
lardner2 <- lardner %>% group_by(Genus_Species, Crown_Class) %>% 
  count()


lardner_c <- lardner %>% group_by(Genus_Species, Crown_Class) %>% 
  mutate(Species_CC_count = n()) %>% 
  ungroup() %>%
  mutate(Species_CC = paste(Genus_Species, Crown_Class, sep = " ")) %>% 
  mutate(Species_CC = if_else(Species_CC %in% c("Coprosma quadrifida Suppressed", "Pittosporum bicolor Suppressed", "Pomaderris aspera Suppressed"), "Misc", Species_CC)) 
  

lardner_c %>% count(Species_CC)


configuration <- Configuration(lardner_c$Ausplot_X, lardner_c$Ausplot_Y, types = lardner_c$Species_CC)
plot(configuration)
window <- ppjsdm::Rectangle_window(c(0, 100), c(0,100))



plotlist <- list()

for (h in levels(configuration$types)) {
  to_optimize <- function(df) { #defines a new function that optimises each row in the dataframe df
    sapply(seq_len(nrow(df)), function(i) { #loops from 1 to number of rows in df, i in function is the index to iterate through the rows of the dataframe 
      set.seed(1)
      fit <- ppjsdm::gibbsm(configuration[h], #create the fit
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
nspecies <- length(levels(configuration$types))
short_range <- matrix(4, nspecies, nspecies)
model <- "exponential"
saturation <- 10
nthreads <- 4

fit_lardner <- ppjsdm::gibbsm(configuration = configuration, 
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



