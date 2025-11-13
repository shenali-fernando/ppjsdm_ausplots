library(terra)
library(sp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(car)
library(ozmaps)

#Load covariate data 
data <- read.csv("data/data_cleaned.csv")
covars <- read.csv("data/site_covariates.csv")

ggplot() +
  geom_sf(data = ozmap_data("country"), fill = "gray90", color = "white") +
  geom_point(data = data_cleaned, aes(x = Longitude, y = Latitude, colour = Site_Name),size = 3) + 
   xlim(140, 150) + 
    ylim(-20, -15)


### Okay, so downloading each ndvi data separately as sections due to data processing limitations
#Initialise  df 
sites <- read.csv("data/data_cleaned.csv") 

sites <- sites %>% dplyr::select(Site_Name, Latitude, Longitude) %>% 
  unique()

rownames(sites) <- NULL


#make df of coordinates
co <- data.frame(lon = sites$Longitude,
                 lat = sites$Latitude) #make a df with just long and lat 

point <- vect(co, geom=c("lon", "lat"))

######## MAXIMUM YEARLY NDVI (2015)
#TAS 
tas1 <- rast("data/NDVI/ndvi_tas1.NDVI.tif")
plot(tas1) #check if no data
crs(tas1) #check crs

value_tas1 <- terra::extract(tas1, point)

sites$a <- value_tas1$ndvi_tas1.NDVI


tas3 <- rast("data/NDVI/ndvi_tas3.NDVI.tif")
plot(tas3)
crs(tas3)

value_tas3 <- terra::extract(tas3, point)

sites$b <- value_tas3$ndvi_tas3.NDVI


tas4 <- rast("data/NDVI/ndvi_tas4.NDVI.tif")
plot(tas4)
crs(tas4)

value_tas4 <- terra::extract(tas4, point)

sites$c <- value_tas4$ndvi_tas4.NDVI


tas5 <- rast("data/NDVI/ndvi_tas5.NDVI.tif")
plot(tas5)
crs(tas5)

value_tas5 <- terra::extract(tas5, point)

sites$d <- value_tas5$ndvi_tas5.NDVI

#VIC

vic1 <- rast("data/NDVI_Max/ndvi_vic1.NDVI.tif")
plot(vic1)
crs(vic1)

value_vic1 <- terra::extract(vic1, point)

sites$e <- value_vic1$ndvi_vic1.NDVI

vic3 <- rast("data/NDVI_Max/ndvi_vic3.NDVI.tif")
plot(vic3)
crs(vic3)

value_vic3 <- terra::extract(vic3, point)

sites$g <- value_vic3$ndvi_vic3.NDVI

#WA

wa1 <- rast("data/NDVI_Max/ndvi_wa1.NDVI.tif")
plot(wa1)
crs(wa1)

value_wa1 <- terra::extract(wa1, point)

sites$h <- value_wa1$ndvi_wa1.NDVI


wa2 <- rast("data/NDVI/ndvi_wa2.NDVI.tif")
plot(wa2)
crs(wa2)

value_wa2 <- terra::extract(wa2, point)

sites$i <- value_wa2$ndvi_wa2.NDVI


#SE NSW 
snsw1 <- rast("data/NDVI/ndvi_snsw1.NDVI.tif")
plot(snsw1)

value_snsw1 <- terra::extract(snsw1, point)
sites$o <- value_snsw1$ndvi_snsw1.NDVI


#N NSW 
nnsw1 <- rast("data/NDVI/ndvi_nnsw1.NDVI.tif")
plot(nnsw1)

value_nnsw1 <- terra::extract(nnsw1, point)
sites$l <- value_nnsw1$ndvi_nnsw1.NDVI

nnsw2 <- rast("data/NDVI/ndvi_nnsw2.NDVI.tif")
plot(nnsw2)

value_nnsw2 <- terra::extract(nnsw2, point)
sites$m <- value_nnsw2$ndvi_nnsw2.NDVI

nnsw3 <- rast("data/NDVI/ndvi_nnsw3.NDVI.tif")
plot(nnsw3)

value_nnsw3 <- terra::extract(nnsw3, point)
sites$n <- value_nnsw3$ndvi_nnsw3.NDVI

#QLD
qld1 <- rast("data/NDVI/ndvi_qld1.NDVI.tif")
plot(qld1)
crs(qld1)

value_qld1 <- terra::extract(qld1, point)
sites$k <- value_qld1$ndvi_qld1.NDVI


qld2 <- rast("data/NDVI/ndvi_qld2.NDVI.tif")
plot(qld2)
crs(qld2)

value_qld2 <- terra::extract(qld2, point)

sites$k <- value_qld2$ndvi_qld2.NDVI

#Combine to make yearly max NDVI
sites <- sites %>% 
  rowwise() %>% 
  mutate(yearlymaxNDVI = mean(c_across(a:o), na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(-(a:o))

#Cool, add into final site covariates df 
site_covariates <- site_covariates %>% 
  dplyr::select(-NDVI)

site_covariates1 <- merge(x = site_covariates, y = sites[, c("Site_Name", "yearlymaxNDVI")], by = "Site_Name", all.x = TRUE)

#save out 
write.csv(site_covariates1, "site_covariates.csv")




###############################################################################################################################################################
# Now, need YEARLY MEAN NDVI 

#get sites df ready again 

#make df of coordinates
co <- data.frame(lon = sites$Longitude,
                 lat = sites$Latitude) #make a df with just long and lat 

point <- vect(co, geom=c("lon", "lat"))

#S NSW 
snsw1 <- rast("data/Mean_ndvi/mndvi_snsw1.NDVI.tif")
plot(snsw1)

value <- terra::extract(snsw1, point)
sites$a <- value$mndvi_snsw1.NDVI


#N NSW 
nnsw1 <- rast("data/Mean_ndvi/mndvi_nnsw1.NDVI.tif")
plot(nnsw1)
value <- terra::extract(nnsw1, point)
sites$b <- value$mndvi_nnsw1.NDVI

nnsw2 <- rast("data/Mean_ndvi/mndvi_nnsw2.NDVI.tif")
plot(nnsw2)
value <- terra::extract(nnsw2, point)
sites$c <- value$mndvi_nnsw2.NDVI

nnsw4 <- rast("data/Mean_ndvi/mndvi_nnsw4.NDVI.tif")
plot(nnsw4)
value <- terra::extract(nnsw4, point)
sites$d <- value$mndvi_nnsw4.NDVI

nnsw5 <- rast("data/Mean_ndvi/mndvi_nnsw5.NDVI.tif")
plot(nnsw5)
value <- terra::extract(nnsw5, point)
sites$e <- value$mndvi_nnsw5.NDVI

#QLD
qld <- rast("data/Mean_ndvi/mndvi_qld.NDVI.tif")
plot(qld)
value <- terra::extract(qld, point)
sites$f <- value$mndvi_qld.NDVI

#VIC
vic1 <- rast("data/Mean_ndvi/mndvi_vic1.NDVI.tif")
plot(vic1)
value <- terra::extract(vic1, point)
sites$g <- value$mndvi_vic1.NDVI

vic2 <- rast("data/Mean_ndvi/mndvi_vic2.NDVI.tif")
plot(vic2)
value <- terra::extract(vic2, point)
sites$h <- value$mndvi_vic2.NDVI

vic3 <- rast("data/Mean_ndvi/mndvi_vic3.NDVI.tif")
plot(vic3)
value <- terra::extract(vic3, point)
sites$i <- value$mndvi_vic3.NDVI


#TAS
tas1 <- rast("data/NDVI_Mean/mndvi_tas1.NDVI.tif")
plot(tas1)
value <- terra::extract(tas1, point)
sites$j <- value$mndvi_tas1.NDVI

tas2 <- rast("data/NDVI_Mean/mndvi_tas2.NDVI.tif")
plot(tas2)
value <- terra::extract(tas2, point)
sites$k <- value$mndvi_tas2.NDVI

tas3 <- rast("data/NDVI_Mean/mndvi_tas3.NDVI.tif")
plot(tas3)
value <- terra::extract(tas3, point)
sites$l <- value$mndvi_tas3.NDVI

tas4 <- rast("data/NDVI_Mean/mndvi_tas4.NDVI.tif")
plot(tas4)
value <- terra::extract(tas4, point)
sites$m <- value$mndvi_tas4.NDVI

tas5 <- rast("data/Mean_ndvi/mndvi_tas5.NDVI.tif")
plot(tas5)
value <- terra::extract(tas5, point)
sites$n <- value$mndvi_tas5.NDVI

#WA 
wa1 <- rast("data/Mean_ndvi/mndvi_wa1.NDVI.tif")
plot(wa1)
value <- terra::extract(wa1, point)
sites$o <- value$mndvi_wa1.NDVI

wa2 <- rast("data/Mean_ndvi/mndvi_wa2.NDVI.tif")
plot(wa2)
value <- terra::extract(wa2, point)
sites$p <- value$mndvi_wa2.NDVI


#Merge

#Combine to make yearly max NDVI
sites <- sites %>% 
  rowwise() %>% 
  mutate(yearlymeanNDVI = mean(c_across(a:p), na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(-(a:p))

#Cool, add into final site covariates df 
site_covariates1 <- merge(x = site_covariates, y = sites[, c("Site_Name", "yearlymeanNDVI")], by = "Site_Name", all.x = TRUE)

#save out 
write.csv(site_covariates1, "site_covariates.csv")

