library(terra)
library(sp)
library(dplyr)
library(tidyr)
library(ggplot2)


### Dataframe to supply coords and index to 
coords <- data_cleaned %>% 
  select(Site_Name, Longitude, Latitude) %>% 
  distinct() #extract site, long and lat 

co <- data.frame(lon = coords$Longitude,
                 lat = coords$Latitude) #make a df with just long and lat 

point <- vect(co, geom=c("lon", "lat"), crs =("+init=epsg:4283")) #make a vector 



## TAS
tas <- rast("C:/Users/shena/Desktop/ausplots/Forest Productivity Index/tas/2015_001.tif")
plot(tas)
tas

value_tas <- terra::extract(tas, point)

coords$forest_productivity_index_tas <- value_tas$`2015_001`


## WA 

wa <-  rast("C:/Users/shena/Desktop/ausplots/Forest Productivity Index/wa/2015_001.tif")
plot(wa)
wa

value_wa <- terra::extract(wa, point)

coords$forest_productivity_index_wa <- value_wa$`2015_001`


## QLD 
qld <-  rast("C:/Users/shena/Desktop/ausplots/Forest Productivity Index/qld/2015_001.tif")
plot(qld)
qld

value_qld <- terra::extract(qld, point)

coords$forest_productivity_index_qld <- value_qld$`2015_001`

## East Vic
eastvic <-  rast("C:/Users/shena/Desktop/ausplots/Forest Productivity Index/eastvic/2015_001.tif")
plot(eastvic)
eastvic

value_eastvic <- terra::extract(eastvic, point)

coords$forest_productivity_index_eastvic <- value_eastvic$`2015_001`


## West Vic 
westvic <-  rast("C:/Users/shena/Desktop/ausplots/Forest Productivity Index/westvic/2015_001.tif")
plot(westvic)
westvic

value_westvic <- terra::extract(westvic, point)

coords$forest_productivity_index_westvic <- value_westvic$`2015_001`


## North NSW 
northnsw <-  rast("C:/Users/shena/Desktop/ausplots/Forest Productivity Index/northnsw/2015_001.tif")
plot(northnsw)
northnsw

#check for a point, eg. Tinebank 
lon <- 152.5267
lat <- -31.2086
df <- data.frame(lon = lon,
                 lat = lat) #make a df with just long and lat 
 
point <- vect(df, geom=c("lon", "lat"), crs =("+init=epsg:4283")) 

value_northnsw <- terra::extract(northnsw, point) #lovely, that's right 

point <- vect(co, geom=c("lon", "lat"), crs =("+init=epsg:4283")) 

value_northnsw <- terra::extract(northnsw, point)

coords$forest_productivity_index_northnsw <- value_northnsw$`2015_001`



## South NSW 
southnsw <-  rast("C:/Users/shena/Desktop/ausplots/Forest Productivity Index/southnsw/2015_001.tif")
plot(southnsw)
southnsw

value_southnsw<- terra::extract(southnsw, point)

coords$forest_productivity_index_southnsw <- value_southnsw$`2015_001`


### Edit df 

coords2 <- coords %>% 
  rowwise() %>% 
  mutate(index = sum(c_across(4:10), na.rm = TRUE))


index_2015 <- coords2 %>% 
  select(Site_Name, Longitude, Latitude, index)

write.csv(index_2015, "forest_productivity_2015.csv")         





############################################################################
########## Plot Productivity v C/D Within group interactions 

cd <- final_ccmod %>% 
  filter(class_from == "Co/dominant") %>% 
  filter(class_to == "Co/dominant") %>% 
  filter(species_from == "Eucalyptus") %>% 
  filter(species_to == "Eucalyptus")


index_2015 <- forest_productivity_2015 %>% 
  filter(! Site_Name %in% c("Dip", "Bird", "Flowerdale")) %>% 
  rename(site = Site_Name)


df_cd <- cd %>% 
  left_join(index_2015 %>% select(site, index), by = "site")

df_cd <- df_cd %>% 
  mutate(sig = ifelse(lo > 0 | hi < 0, "1", "0"))

df_cd <- df_cd %>% 
  mutate(region = site) %>% 
  mutate(region = if_else(region %in% c("Carey", "Dombakup", "Warren", "Dawson", "Giants", "Sutton", "Frankland", "Clare", "Collins"), 
                           "WA", region)) %>% 
  mutate(region = if_else(region %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Lardner", "Weeaproinah", "Turtons"),
                           "SE Highlands VIC", region)) %>% 
  mutate(region = if_else(region %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"), 
                           "S NSW", region)) %>% 
  mutate(region = if_else(region %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans"), 
                           "N NSW", region)) %>% 
  mutate(region = if_else(region %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton"), 
                           "FN QLD", region)) %>% 
  mutate(region = if_else(region %in% c("Bird", "Supersite","NorthStyx", "ZigZag", "MtField", "Weld"), 
                           "S.Ranges TAS", region)) %>% 
  mutate(region = if_else(region %in% c("Caveside",  "Mackenzie", "BlackRiver", "BondTier"), 
                           "King-N.Slopes TAS", region)) %>% 
  mutate(region = if_else(region %in% c("BenRidge", "MtMaurice"), 
                          "Ben Lomond TAS", region)) 

df_cd %>% count(region)

noregion <- ggplot(data = df_cd, 
       aes(x = index, 
           y = alpha,
           fill = sig)) + 
  geom_hline(yintercept = 0, colour = "red", linewidth = .75, linetype = 3) + 
  geom_point(shape = 21, size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), size = 0.5, colour = "gray40") + 
  scale_fill_manual(values = c("gray80", "black"), 
                    name = "", 
                    labels = c("Not Sig", "Sig")) + 
  xlab("Forest Productivity Index") + 
  ylab("Within C/D Coefficient") + 
  theme_bw()
  


df_cd <- df_cd %>% 
  mutate(sig = ifelse(lo > 0 | hi < 0, "1", NA))

df_cd <- df_cd %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(region), NA))

colours <-  hcl.colors(n = 8, "Set 2")

region <- ggplot(data = df_cd, 
                 aes(x = index, 
                     y = alpha,
                     colour = region,
                     fill = fill_col)) + 
  geom_hline(yintercept = 0, colour = "red", linewidth = .75, linetype = 3) + 
  geom_point(shape = 21, size = 3, stroke = 1.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), size = 0.5, colour = "gray40") + 
  scale_colour_manual(values = colours, 
                      name = "Site Geographic Region") +
  scale_fill_manual(values = colours,
                    na.value = "gray100",
                    name = "", 
                    labels = c("Not Significant", "Significant", "", "", "", "", "", "")) + 
  guides(colour = guide_legend(override.aes = #change the figure legend 
                                 list(fill = colours)), 
         fill = guide_legend(override.aes = list(shape = c( 
           1, 19,
           NA, NA, NA, NA, NA, NA)
         ))) + 
  xlab("Forest Productivity Index") + 
  ylab("Within C/D Coefficient") + 
  theme_bw()

region






xlibrary(lme4)

# Fit the model
model <- lm(alpha ~ index, data = df_cd)

summary(model)

#check assumptions
par(mfrow = c(2, 2))
plot(model)



