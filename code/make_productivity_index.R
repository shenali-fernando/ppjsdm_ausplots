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

