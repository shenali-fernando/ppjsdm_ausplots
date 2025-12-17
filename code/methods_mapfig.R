library(dplyr)
library(ggforce)
library(forcats)
library(patchwork)
library(OpenStreetMap)
library(osmdata)
library(ozmaps)


#Load data 
data <- read.csv("data/data_cleaned.csv")

basemap <- ggplot() +
  geom_sf(data = ozmap(), fill = "gray100") +
  geom_point(data = data, aes(x = Longitude, y = Latitude, colour = Site_Name),
             size = 2.5) + 
  xlim(c(110, 156)) + 
  theme_classic() + 
  theme(legend.position = "none")


## Add regions to data 
data2 <- data %>% 
  mutate(georegion = case_when(Site_Name %in% c("Weeaproinah", "Turtons", "Lardner") ~ "S_VIC", 
                               Site_Name %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek") ~ "N_VIC", 
                               Site_Name %in% c("Dawson", "Frankland", "Clare", "Giants") ~ "S_WA",
                               Site_Name %in% c("Carey", "Dombakup", "Warren",  "Sutton","Collins") ~ "N_WA", 
                               Site_Name %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton") ~ "QLD", 
                               Site_Name %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans") ~ "N_NSW", 
                               Site_Name %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo") ~ "S_NSW", 
                               Site_Name %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx") ~ "d_TAS", 
                               Site_Name %in% c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Flowerdale", "Dip") ~ "o_TAS")) %>%  
  mutate(region = case_when(georegion %in% c("S_VIC", "N_VIC") ~ "VIC", 
                            georegion %in% c("S_NSW") ~ "S_NSW",
                            georegion %in% c("S_WA", "N_WA") ~ "WA", 
                            georegion %in% c("QLD") ~ "QLD",
                            georegion %in% c("N_NSW") ~ "N_NSW", 
                            georegion %in% c("o_TAS", "d_TAS") ~ "TAS"))

data2 %>% count(region)


basemap <- ggplot() +
  geom_sf(data = ozmap(), fill = "gray100") + 
  geom_point(data = data2, 
             aes(x = Longitude, y = Latitude, colour = Site_Name),
             size = 2.5) + 
  geom_mark_rect(data = data2, 
                 aes(x = Longitude, y = Latitude, group = region), 
                 expand = unit(3, "mm"), 
                 fill = NA, 
                 colour = "black", 
                 linewidth = 0.35) +
  xlim(c(110, 156)) + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.line = element_blank())
basemap


wa <- data2 %>% filter(region == "WA")  

ggplot() + 
  geom_sf(data = ozmap(), fill = "gray100") + 
  geom_point(data = wa, 
             aes(x = Longitude, y = Latitude, colour = Site_Name),
             size = 2.5) + 
  xlim(c(114, 124)) + 
  ylim(c(-36, -28)) + 
  theme_classic()




map <- openmap(upperLeft = c(-1, 111), 
               lowerRight = c(-52, 158),
               type = "esri")
plot(map)
class(map)
