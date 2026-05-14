library(stringr)
library(spatstat)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
library(ppjsdm)
library(dplyr)




#### GEOGRAPHIC REGION ANALYSIS FOR SPECIES 

source("code/make_georegion_fun.R") #source the function needed 

data <- read.csv("data/data_cleaned.csv")

#Run for sites 
svic <- make_georegion(sites = c("Weeaproinah", "Turtons", "Lardner"), #vector of sites 
                       group_type = "species", #options of species, crown_class, species 
                       threshold = , #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "exponential", 
                       medium_model = "tanh")
svic$sum


nvic <- make_georegion(sites = c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek"), #vector of sites 
                       group_type = "species", #options of species, crown_class, species 
                       threshold = 16, #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "exponential", 
                       medium_model = "tanh")
nvic$sum



swa <- make_georegion(sites = c("Dawson", "Frankland", "Clare", "Giants"), #vector of sites 
                      group_type = "species", #options of species, crown_class, species 
                      threshold = 13, #if a group has less than 10 individuals remove
                      short_range = 8, #short-range interaction distance (set same for all species)
                      medium_range = 40,
                      long_range = 100,
                      short_model = "exponential", 
                      medium_model = "tanh")
swa$sum


nwa <-  make_georegion(sites = c("Carey", "Dombakup", "Warren",  "Sutton","Collins"), #vector of sites 
                       group_type = "species", #options of species, crown_class, species 
                       threshold = 10, #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "exponential", 
                       medium_model = "tanh")
nwa$sum


qld <- make_georegion(sites = c("Baldy", "Koombooloomba", "Lamb Range", "Herberton"), #vector of sites 
                      group_type = "species", #options of species, crown_class, species 
                      threshold = 10, #if a group has less than 10 individuals remove
                      short_range = 10, #short-range interaction distance (set same for all species)
                      medium_range = 40,
                      long_range = 100,
                      short_model = "exponential", 
                      medium_model = "bump")
qld$sum


nnsw <- make_georegion(sites =  c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans"), #vector of sites 
                       group_type = "species", #options of species, crown_class, species 
                       threshold = 15, #if a group has less than 10 individuals remove
                       short_range = 10, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "bump", 
                       medium_model = "bump")

nnsw$sum



snsw <- make_georegion(sites = c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"), #vector of sites 
                       group_type = "species", #options of species, crown_class, species 
                       threshold = 10, #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "exponential", 
                       medium_model = "bump")
snsw$sum



dtas <- make_georegion(sites = c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx"), #vector of sites 
                       group_type = "species", #options of species, crown_class, species 
                       threshold = 12, #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range =100,
                       short_model = "bump", 
                       medium_model = "bump")
dtas$sum
dtas$fit$aic


otas <- make_georegion(sites = c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite"), #vector of sites 
                       group_type = "species", #options of species, crown_class, species 
                       threshold = 12, #if a group has less than 10 individuals remove
                       short_range = 8, #short-range interaction distance (set same for all species)
                       medium_range = 40,
                       long_range = 100,
                       short_model = "bump", 
                       medium_model = "bump")
otas$sum
otas$fit$aic



### Visualisation 

nvic_bp <- ppjsdm::box_plot(fit = nvic$fit, 
                            summ = nvic$sum, 
                            coefficient = "alpha", 
                            which = "all", 
                            title = "North VIC Alpha", 
                            text_size = 10)


svic_bp <- ppjsdm::box_plot(fit = svic$fit, 
                            summ = svic$sum, 
                            coefficient = "alpha", 
                            which = "all", 
                            title = "South VIC Alpha",
                            text_size = 10)

svic_bp + nvic_bp 

swa_bp <- ppjsdm::box_plot(fit = swa$fit, 
                           summ = swa$sum, 
                           coefficient = "alpha", 
                           which = "all", 
                           title = "South WA Alpha",
                           text_size = 10)

nwa_bp <- ppjsdm::box_plot(fit = nwa$fit, 
                           summ = nwa$sum, 
                           coefficient = "alpha", 
                           which = "all", 
                           title = "North WA Alpha",
                           text_size = 10)

swa_bp + nwa_bp 


qld_bp <- ppjsdm::box_plot(fit = qld$fit, 
                           summ = qld$sum, 
                           coefficient = "alpha", 
                           which = "within", 
                           title = "QLD Within Alpha",
                           text_size = 10)

nnsw_bp  <- ppjsdm::box_plot(fit = nnsw$fit, 
                             summ = nnsw$sum, 
                             coefficient = "alpha", 
                             which = "within", 
                             title = "North NSW Within Alpha",
                             how = "only",
                             text_size = 10)


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


#Chuck everything into a df
df <- make_summary_df(fits = list(s_vic = svic$fit, 
                                  n_vic = nvic$fit, 
                                  qld = qld$fit, 
                                  n_nsw = nnsw$fit, 
                                  s_nsw = snsw$fit, 
                                  s_wa = swa$fit, 
                                  n_wa = nwa$fit, 
                                  d_tas = dtas$fit, 
                                  o_tas = otas$fit), 
                summ = list(svic$sum, nvic$sum, qld$sum, nnsw$sum, snsw$sum, swa$sum, nwa$sum, dtas$sum, otas$sum), 
                coefficient = "alpha", 
                only_statistically_significant = TRUE, 
                which = "all", 
                full_names = NULL, 
                compute_confidence_intervals = TRUE, 
                classes = NULL, 
                involving = NULL, 
                how = "one")

#write.csv(df, "species_geo_df.csv")

#get only intra-specific eucalypt interactions 
df1 <- df %>% 
  filter(from == to) %>% 
  filter(str_starts(from, "(Eucalyptus|Corymbia|Syncarpia)"))

df1 <- df1 %>%
  arrange(-alpha) %>%
  mutate(from = factor(from, levels = unique(from)))

ggplot(data = df1, 
       aes(x = alpha, 
           y = from, 
           colour = Fit)) + 
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = lo, 
                    xmax = hi), 
                size = 1) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()


box_plot()