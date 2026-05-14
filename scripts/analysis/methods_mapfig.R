library(dplyr)
library(ggforce)
library(forcats)
library(patchwork)
library(OpenStreetMap)
library(ggmap)
library(ozmaps)
library(ggpubr)


#Load data 
data <- read.csv("data/data_cleaned.csv")

basemap <- ggplot() +
  geom_sf(data = ozmap(), fill = "gray100") +
  geom_point(data = data, aes(x = Longitude, y = Latitude, colour = Site_Name),
             size = 2.5) + 
  xlim(c(110, 156)) + 
  theme_classic() + 
  theme(legend.position = "none")



## For ggmap need to have an registered api 
register_google(key = "AIzaSyDQDBJp4_qEcXsVK1VJG8h8PPx0R9pnWBE", 
                write = TRUE) #will charge everytime, be careful 

map <- get_googlemap(center = c(lon = 134.35, lat = -25.60), 
                     zoom = 4, 
                     maptype = "satellite")



## Get only long and lat, add regions to data 
data2 <- data %>% 
  dplyr::select(Site_Name, Latitude, Longitude) %>% 
  distinct() %>% 
  mutate(georegion = case_when(Site_Name %in% c("Weeaproinah", "Turtons", "Lardner", "ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek") ~ "VIC", 
                               Site_Name %in% c("Dawson", "Frankland", "Clare", "Giants") ~ "S_WA",
                               Site_Name %in% c("Carey", "Dombakup", "Warren",  "Sutton","Collins") ~ "N_WA", 
                               Site_Name %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton") ~ "QLD", 
                               Site_Name %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans") ~ "N_NSW", 
                               Site_Name %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo") ~ "S_NSW", 
                               Site_Name %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx") ~ "d_TAS", 
                               Site_Name %in% c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Flowerdale", "Dip") ~ "o_TAS")) %>%  
  mutate(region = case_when(georegion %in% c("VIC", "S_NSW", "o_TAS", "d_TAS") ~ "SE_AUS", 
                            georegion %in% c("S_WA", "N_WA") ~ "WA", 
                            georegion %in% c("QLD") ~ "QLD", 
                            georegion %in% c("N_NSW") ~ "N_NSW"))

data2 %>% count(region)

### Using ggmap - make a basemap

basemap <- ggmap(map) + 
  geom_point(data = data2, 
             aes(x = Longitude, y = Latitude, colour = georegion), 
             size = 2, 
             position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c("#ED90A4", "#D8A06A", "#ABB150", "#62BE79", "#00C1B2",
                                "#48B8DE", "#ACA2EC", "#E190D6")) + 
  xlim(c(112, 155)) + 
  ylim(c(-44, -10)) + 
  geom_mark_rect(data = data2, 
                 aes(x = Longitude, y = Latitude, group = region), 
                 expand = unit(3, "mm"), 
                 fill = NA, 
                 colour = "white", 
                 linewidth = 0.45) +
  annotate("text", x = 113.75, y = -34.5, label = "c", colour = "white") + 
  annotate("text", x = 141.5, y = -40, label = "b", colour = "white") + 
  annotate("text", x = 150.5, y = -31, label = "a", colour = "white") + 
  annotate("text", x = 143.5, y = -17.25, label = "a", colour = "white") +
theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none")
  
basemap

# SE AUS
d <- data2 %>% filter(georegion == "VIC")
vic_map <- get_googlemap(center = c(lon = 144.5, lat = -38), 
                        zoom = 8, 
                        maptype = "satellite")
vic <- ggmap(vic_map) + 
  geom_point(data = d, 
             aes(x = Longitude, y = Latitude, colour = georegion),
             size = 2, 
             position = position_dodge(width = 0.1)) + 
  ggtitle("b. Southeastern Australia") +
  scale_color_manual(values = "#E190D6", labels = "VIC", name = "") +
  ylim(c(-39, -37.25)) + 
  xlim(c(143.2, 146))
vic

d <- data2 %>% filter(georegion %in%  c("d_TAS", "o_TAS"))
tas_map <- get_googlemap(center = c(lon = 146.5, lat = -42), 
                         zoom = 7, 
                         maptype = "satellite")
tas <- ggmap(tas_map) + 
  geom_point(data = d, 
             aes(x = Longitude, y = Latitude, colour = georegion),
             size = 2, 
             position = position_dodge(width = 0.1)) + 
  scale_colour_manual(values = c("#ED90A4", "#62BE79"), labels = c("E.delegatensis TAS", "E.obliqua TAS"), name = "") +
  xlim(c(144.5, 148.5)) + 
  ylim(c(-43.5, -40.5))
tas


d <- data2 %>% filter(georegion == "S_NSW")
snsw_map <- get_googlemap(center = c(lon = 149.5, lat = -37), 
                         zoom = 9, 
                         maptype = "satellite")
snsw <- ggmap(snsw_map) + 
  geom_point(data = d, 
             aes(x = Longitude, y = Latitude, colour = georegion),
             size = 2, 
             position = position_dodge(width = 0.1)) + 
  scale_color_manual(values = "#48B8DE", labels = "South NSW", name = "") 
snsw



# WA 
d <- data2 %>% filter(region == "WA")
wa_map <- get_googlemap(center = c(lon = 116.2, lat = -34.5), 
                          zoom = 8, 
                          maptype = "satellite")
wa <- ggmap(wa_map) + 
  geom_point(data = d, 
             aes(x = Longitude, y = Latitude, colour = georegion),
             size = 2, 
             position = position_dodge(width = 0.1)) + 
  scale_color_manual(values = c("#ABB150", "#ACA2EC"),labels = c("E.diversicolor WA", "E.jacksonii WA"), name = "") + 
  ggtitle("c. Western Australia") + 
  theme(legend.position = "bottom", 
        legend.direction = "vertical")
wa



# NE AUS 
d <- data2 %>% filter(georegion == "QLD")
qld_map <-  get_googlemap(center = c(lon = 145.6, lat = -17.5), 
                      zoom = 9, 
                      maptype = "satellite")
qld <- ggmap(qld_map) + 
  geom_point(data = d, 
             aes(x = Longitude, y = Latitude, colour = georegion), 
             size = 2.5, 
             position = position_dodge(width = 0.1)) + 
  scale_color_manual(values = "#00C1B2", labels = "QLD", name = "") +
  ggtitle("a. Northeastern Australia")
qld


d <- data2 %>% filter(georegion == "N_NSW")
nnsw_map <- get_googlemap(center = c(lon = 152.2, lat = -32), 
                          zoom = 7, 
                          maptype = "satellite")
nnsw <- ggmap(nnsw_map) + 
  geom_point(data = d, 
             aes(x = Longitude, y = Latitude, colour = georegion), 
             size = 2, 
             position = position_dodge(width = 0.1)) + 
  scale_color_manual(values = "#D8A06A", labels = "North NSW", name = "") +
  ylim(c(-33.5, -29.5))
nnsw


# chuck everything together 

bottom_row <- vic + snsw + tas + plot_layout(guides = 'collect') &
  theme(legend.margin = margin(0, 0, 0, 0), 
        legend.spacing.x = unit(0, "cm"),
        legend.spacing.y = unit(0, "cm"))

right <- qld / nnsw +  plot_layout(guides = 'collect') &
  theme(legend.margin = margin(0, 0, 0, 0), 
        legend.spacing.x = unit(0, "cm"),
        legend.spacing.y = unit(0, "cm"))

middle_row <- wa + basemap + right + 
  plot_layout(widths = c(1, 2, 1))

final_map <- middle_row / bottom_row + 
  plot_layout(heights = c(2, 1))

final_map

ggsave("map.png", final_map, 
       height = 14, width = 22, dpi = 300)







### Using ozmap 
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
