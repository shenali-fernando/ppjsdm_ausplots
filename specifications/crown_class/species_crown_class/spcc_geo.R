library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
library(ozmaps)


#### GEOGRAPHIC REGION ANALYSIS  
data <- read.csv("data/data_cleaned.csv")

#### Doing a species by crown class analysis by splitting by geographic georegion 
# Use the function from make_georegion_fun.R 

#We can change the medium and long range distances to put more or less weight on the between-plot interactions 

#Nice, got a function. But need to be able to get a well-parameterised model when we have all distances specified
svic <- make_georegion(sites = c("Weeaproinah", "Turtons", "Lardner"), #vector of sites 
                      group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                      threshold = 10, #if a group has less than 10 individuals remove
                      short_range = 8, #short-range interaction distance (set same for all species)
                      medium_range = 40,
                      long_range = 100,
                      short_model = "exponential", 
                      medium_model = "tanh")
svic$sum


nvic <- make_georegion(sites = c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek"), #vector of sites 
                       group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                       threshold = 10, #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "exponential", 
                       medium_model = "tanh")
nvic$sum


swa <- make_georegion(sites = c("Dawson", "Frankland", "Clare", "Giants"), #vector of sites 
                      group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                      threshold = 10, #if a group has less than 10 individuals remove
                      short_range = 8, #short-range interaction distance (set same for all species)
                      medium_range = 40,
                      long_range = 100,
                      short_model = "exponential", 
                      medium_model = "tanh")
swa$sum


nwa <-  make_georegion(sites = c("Carey", "Dombakup", "Warren",  "Sutton","Collins"), #vector of sites 
                       group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                       threshold = 10, #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "exponential", 
                       medium_model = "tanh")
nwa$sum


qld <- make_georegion(sites = c("Baldy", "Koombooloomba", "Lamb Range", "Herberton"), #vector of sites 
                      group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                      threshold = 10, #if a group has less than 10 individuals remove
                      short_range = 8, #short-range interaction distance (set same for all species)
                      medium_range = 40,
                      long_range = 100,
                      short_model = "exponential", 
                      medium_model = "tanh")
qld$sum


nnsw <- make_georegion(sites =  c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans"), #vector of sites 
                      group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                      threshold = 10, #if a group has less than 10 individuals remove
                      short_range = 8, #short-range interaction distance (set same for all species)
                      medium_range = 40,
                      long_range = 100,
                      short_model = "exponential", 
                      medium_model = "tanh")

nnsw$sum

snsw <- make_georegion(sites = c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"), #vector of sites 
                      group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                      threshold = 10, #if a group has less than 10 individuals remove
                      short_range = 8, #short-range interaction distance (set same for all species)
                      medium_range = 40,
                      long_range = 100,
                      short_model = "exponential", 
                      medium_model = "tanh")
snsw$sum



dtas <- make_georegion(sites = c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx"), #vector of sites 
                       group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                       threshold = 12, #if a group has less than 10 individuals remove
                       short_range = 10, #short-range interaction distance (set same for all species)
                       medium_range = 0,
                       long_range =0,
                       short_model = "bump", 
                       medium_model = "bump")
dtas$sum



otas <- make_georegion(sites = c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite"), #vector of sites 
                       group_type = "species_crown_class", #options of species, crown_class, species_crown_class 
                       threshold = 12, #if a group has less than 10 individuals remove
                       short_range = 10, #short-range interaction distance (set same for all species)
                       medium_range = 25,
                       long_range = 75,
                       short_model = "bump", 
                       medium_model = "bump")
otas$sum






### Visualisation 

nvic_bp <- ppjsdm::box_plot(fit = nvic$fit, 
                            summ = nvic$sum, 
                            coefficient = "alpha", 
                            which = "within", 
                            title = "North VIC Within Alpha", 
                            text_size = 10)


svic_bp <- ppjsdm::box_plot(fit = svic$fit, 
                            summ = svic$sum, 
                            coefficient = "alpha", 
                            which = "within", 
                            title = "South VIC Within Alpha",
                            text_size = 10)

svic_bp + nvic_bp 

swa_bp <- ppjsdm::box_plot(fit = swa$fit, 
                          summ = swa$sum, 
                          coefficient = "alpha", 
                          which = "within", 
                          title = "South WA Within Alpha",
                          text_size = 10)

nwa_bp <- ppjsdm::box_plot(fit = nwa$fit, 
                              summ = nwa$sum, 
                              coefficient = "alpha", 
                              which = "within", 
                              title = "North WA Within Alpha",
                           text_size = 10)

swa_bp + nwa_bp 


qld_bp <- ppjsdm::box_plot(fit = qld$fit, 
                         summ = qld$sum, 
                         coefficient = "alpha", 
                         which = "within", 
                         title = "QLD Within Alpha",
                         text_size = 10, 
                         xmin = -2.5, 
                         xmax = 1)

nnsw_bp  <- ppjsdm::box_plot(fit = nnsw$fit, 
                          summ = nnsw$sum, 
                          coefficient = "alpha", 
                          which = "between", 
                          title = "North NSW Within Alpha",
                          involving = c("Eucalyptus grandis Co/dominant", "Eucalyptus pilularis Co/dominant"),
                          how = "only",
                          text_size = 10,
                          xmin = -2.5, 
                          xmax = 1)

                          
qld_bp + nnsw_bp                          
                          
snsw_bp <- ppjsdm::box_plot(fit = snsw$fit, 
                         summ = snsw$sum, 
                         coefficient = "alpha", 
                         which = "within", 
                         title = "South NSW Within Alpha", 
                         text_size =10)

dtas_bp <- ppjsdm::box_plot(fit = dtas$fit, 
                         summ = dtas$sum, 
                         coefficient = "alpha", 
                         which = "within", 
                         title = "d TAS Within Alpha",
                          text_size = 10)

otas_bp <- ppjsdm::box_plot(fit = otas$fit, 
                        summ = otas$sum, 
                        coefficient = "alpha", 
                        which = "within", 
                        title = "obliqua TAS Within Alpha",
                        text_size = 10, 
                        xmin = -4)

snsw_bp + dtas_bp + otas_bp + plot_layout(nrow = 1)
