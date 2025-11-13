library(terra)
library(sp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(car)

#Load covariate data 
covars <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/data/site_covariates.csv")

#need to rename some cols 
covars <- covars %>% 
  rename(MAT = MAT..C.) %>% 
  rename(MAP = MAP..mm.) %>% 
  rename(elev = elevation..m.) %>% 
  rename(FPI = FPI..3PG.)

#Need to understand the correlation between MAT, MAP, and productivity index (FPI) 
#So, create a scatterplot matrix 
pairs(~MAT + MAP + elev + FPI, data = covars)

#Don't actually think elevation is strongly correlated with the other covariates; so
pairs(~MAT + MAP + FPI, data = covars)

#can add some variance envelopes 
scatterplotMatrix(~MAT + MAP + FPI, data = covars,
                  diagonal = FALSE,             # Remove kernel density estimates
                  regLine = FALSE,      # Linear regression line width
                  smooth = list(col.smooth = "red",   # Non-parametric mean color
                                col.spread = "blue")) # Non-parametric variance color


## Run some correlations one by one 
#MAT v MAP 
ggplot(data = covars, aes(x = MAT, y = MAP)) +
  geom_point() + # Show dots
  geom_text(aes(label = Site_Name),
    nudge_x = 0.25, nudge_y = 0.75, 
    check_overlap = F) + 
  theme_bw()

#MAT v FPI 
ggplot(data = covars, aes(x = MAT, y = FPI)) +
  geom_point() + # Show dots
  geom_text(aes(label = Site_Name),
            nudge_x = 0.25, nudge_y = 0.25, 
            check_overlap = F) + 
  geom_vline(xintercept = 11, colour = "red", linetype = "dashed", linewidth = 1) + 
  theme_bw()

#MAP v FPI
ggplot(data = covars, aes(x = MAP, y = FPI)) +
  geom_point() + # Show dots
  geom_text(aes(label = Site_Name),
            nudge_x = 0.25, nudge_y = 0.25, 
            check_overlap = F) + 
  theme_bw()

##Correlation test 

cor <- stats::cor(covars[, 4:7])
corrplot::corrplot(cor, method = "number", order = 'hclust') #not a straight correlation between MAT and FPI anyways




e <- final_ccmod %>% 
  filter(species_from == "Eucalyptus") %>% 
  filter(species_to == "Eucalyptus") %>% 
  filter(class_from == "Suppressed") %>% 
  filter(class_to == "Suppressed")

e <- e %>% dplyr::select(alpha, site)

site_covars <- site_covars %>% 
left_join(e, by = "site")


site_covars <- site_covars %>% 
  rename(cc_s =alpha)


e2 <- final_spcc_mod %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(class_from == "Co/dominant") %>% 
  filter(class_to == "Co/dominant")

sum <- e2 %>% 
  group_by(site) %>% 
  summarise(spcc_cd = median(alpha))
  

pairs(~`MAT (C)` + `FPI (3PG)`+ NDVI, data = site_covars)

site_covars$NDVI <- iconv(site_covars$NDVI, from = "", to = "UTF-8", sub = "")
site_covars$NDVI <- gsub("[^0-9.-]", "", site_covars$NDVI)  # strip non-numeric chars
site_covars$NDVI <- as.numeric(site_covars$NDVI)


ggplot(data = site_covars, 
       aes(x = cc_cd, y = `MAP (mm)`, 
           color = region)) +
  geom_point(size = 3) + 
  theme_bw()

site_covars <- site_covars %>% 
  mutate(region = site) %>% 
  mutate(region = if_else(region %in% c("Carey", "Dombakup", "Warren", "Dawson", "Giants", "Sutton", "Frankland", "Clare", "Collins"), 
                          "WA", region)) %>% 
  mutate(region = if_else(region %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Weeaproinah", "Turtons", "Lardner"), 
                          "SE_VIC", region)) %>% 
  mutate(region = if_else(region %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"), 
                          "SE_NSW", region)) %>% 
  mutate(region = if_else(region %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans"), 
                          "N_NSW", region)) %>% 
  mutate(region = if_else(region %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton"), 
                          "QLD", region)) %>% 
  mutate(region = if_else(region %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx"), 
                          "d_TAS", region)) %>% 
  mutate(region = if_else(region %in% c("BondTier", "BlackRiver", "Weld", "Weld", "MtField", "ZigZag", "Supersite"), 
                          "o_TAS", region)) 


r <- rast("data/climateEngine_download.NDVI.tif")
plot(r)


