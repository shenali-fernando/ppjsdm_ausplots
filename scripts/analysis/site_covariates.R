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
# covars <- covars %>% 
#   rename(MAT = MAT..C.) %>% 
#   rename(MAP = MAP..mm.) %>% 
#   rename(elev = elevation..m.) %>% 
#   rename(FPI = FPI..3PG.)


#add regions 
site_covariates <- site_covariates %>% 
  mutate(georegion = case_when(Site_Name %in% c("Weeaproinah", "Turtons", "Lardner") ~ "S_VIC", 
                               Site_Name %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek") ~ "N_VIC", 
                               Site_Name %in% c("Dawson", "Frankland", "Clare", "Giants") ~ "S_WA",
                               Site_Name %in% c("Carey", "Dombakup", "Warren",  "Sutton","Collins") ~ "N_WA", 
                               Site_Name %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton") ~ "QLD", 
                               Site_Name %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans") ~ "N_NSW", 
                               Site_Name %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo") ~ "S_NSW", 
                               Site_Name %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx") ~ "d_TAS", 
                               Site_Name %in% c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Flowerdale", "Dip") ~ "o_TAS")) %>% 
  mutate(region = case_when(georegion %in% c("S_VIC", "N_VIC", "S_NSW") ~ "SE_AUS", 
                            georegion %in% c("S_WA", "N_WA") ~ "WA",
                            georegion %in% c("QLD", "N_NSW") ~ "N_AUS",  
                            georegion %in% c("o_TAS", "d_TAS") ~ "TAS"))



a <- site_covariates %>% 
  group_by(region) %>% 
  summarise(mean_map = mean(MAP), 
            max_map = max(MAP), 
            min_map = min(MAP), 
            mean_mat = mean(MAT), 
            max_mat = max(MAT), 
            min_mat = min(dMAT))


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


################################################################################
##### Association between weighted site x group means and site covariates ######

df <- df_add3_5_t15 #use this results 

#add a which_species col
df <- df %>% 
  mutate(which_species = ifelse(species_from == species_to, "within", "between"))

#check groups are not doubled up 
df %>% group_by(which_species) %>% count(group)

#compute a weighted mean for each which_species x group for each site 
w_means <- df %>% 
  group_by(which_species, group, site) %>% 
  summarise(w.mean = weighted.mean(alpha, se), 
            median = median(alpha)) %>% 
  ungroup()

##start with Within-species interactions 
within_wmeans <- w_means %>% 
  filter(which_species == "within")

site_covariates <- site_covariates %>% 
  rename(site = Site_Name)

#merge everything 
sum_within <- merge(x = within_wmeans, y = site_covariates, by = "site" , all.x = TRUE)

#for plotting make some new columns 
sum_within2 <- sum_within %>%
  mutate(fg = str_extract(group, "^\\w+")) %>% 
  mutate(class_to = str_extract(group, "\\w+$")) %>% 
  mutate(class_from = str_extract(group, "(?<=\\.)\\w+(?=_)"))

sum_within2 <- sum_within2 %>% 
  mutate(class_int = paste0(class_from, sep = "_", class_to)) %>% 
  mutate(class_int = ifelse(class_int == "small_large", "large_small", class_int))



#Visualisation of association between weighted mean within-species int and site covairates 
ggplot(data = sum_within2, 
       aes(y = w.mean, 
           x = MAT,
           colour = fg
       )) + 
  facet_grid(~class_int)+
  geom_smooth(method = "glm") +
  geom_point()+ 
  theme_bw()



#Correlations: 
pairs(~w.mean + MAT + Latitude + MAP + GPP + yearlymaxNDVI + PET, data = sum_within2)


#I think there isn't a convincing pattern...