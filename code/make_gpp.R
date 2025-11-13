library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)

## Extract site covariates for gross primary production from climate engine 

#Load covariate data 
data <- read.csv("data/data_cleaned.csv")

sites <- data %>% dplyr::select(Site_Name, Latitude, Longitude) %>% 
  unique()
rownames(sites) <- NULL


#make df of coordinates
co <- data.frame(lon = sites$Longitude,
                 lat = sites$Latitude) #make a df with just long and lat 

point <- vect(co, geom=c("lon", "lat"))


#WA
wa1 <- rast("data/GPP/gpp_wa1.GPP.tif")
plot(wa1)
crs(wa1)
value <- terra::extract(wa1, point)
sites$a <- value$gpp_wa1.GPP

wa2 <- rast("data/GPP/gpp_wa2.GPP.tif")
plot(wa2)
crs(wa2)
value <- terra::extract(wa2, point)
sites$b <- value$gpp_wa2.GPP

wa3 <- rast("data/GPP/gpp_wa3.GPP.tif")
plot(wa3)
crs(wa3)
value <- terra::extract(wa3, point)
sites$c <- value$gpp_wa3.GPP


#TAS
tas1 <- rast("data/GPP/gpp_tas1.GPP.tif")
plot(tas1)
value <- terra::extract(tas1, point)
sites$d <- value$gpp_tas1.GPP


tas3 <- rast("data/GPP/gpp_tas3.GPP.tif")
plot(tas3)
value <- terra::extract(tas3, point)
sites$f <- value$gpp_tas3.GPP

tas2 <- rast("data/GPP/gpp_tas2.GPP.tif")
plot(tas2)
value <- terra::extract(tas2, point)
sites$e <- value$gpp_tas2.GPP

tas4 <- rast("data/GPP/gpp_tas4.GPP.tif")
plot(tas4)
value <- terra::extract(tas4, point)
sites$g <- value$gpp_tas4.GPP

tas5 <- rast("data/GPP/gpp_tas5.GPP.tif")
plot(tas5)
value <- terra::extract(tas5, point)
sites$h <- value$gpp_tas5.GPP

tas6 <- rast("data/GPP/gpp_tas6.GPP.tif")
plot(tas6)
value <- terra::extract(tas6, point)
sites$i <- value$gpp_tas6.GPP

#QLD
qld1 <- rast("data/GPP/gpp_qld1.GPP.tif")
plot(qld1)
value <- terra::extract(qld1, point)
sites$j <- value$gpp_qld1.GPP

#NNSW 
nnsw1 <- rast("data/GPP/gpp_nnsw1.GPP.tif")
plot(nnsw1)
value <- terra::extract(nnsw1, point)
sites$k <- value$gpp_nnsw1.GPP

nnsw2 <- rast("data/GPP/gpp_nnsw2.GPP.tif")
plot(nnsw2)
value <- terra::extract(nnsw2, point)
sites$l <- value$gpp_nnsw2.GPP

nnsw3 <- rast("data/GPP/gpp_nnsw3.GPP.tif")
plot(nnsw3)
value <- terra::extract(nnsw3, point)
sites$m <- value$gpp_nnsw3.GPP


#SE NSW 
sensw1 <- rast("data/GPP/gpp_sensw1.GPP.tif")
plot(sensw1)
value <- terra::extract(sensw1, point)
sites$n <- value$gpp_sensw1.GPP


#vic 
vic1 <- rast("data/GPP/gpp_vic1.GPP.tif")
plot(vic1)
value <- terra::extract(vic1, point)
sites$o <- value$gpp_vic1.GPP

vic2 <- rast("data/GPP/gpp_vic2.GPP.tif")
plot(vic2)
value <- terra::extract(vic2, point)
sites$p <- value$gpp_vic2.GPP



## Combine everything together 
sites <- sites %>% 
  rowwise() %>% 
  mutate(GPP = mean(c_across(a:p), na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(-(a:p))

#Cool, add into final site covariates df 
site_covariates1 <- merge(x = site_covariates, y = sites[, c("Site_Name", "GPP")], by = "Site_Name", all.x = TRUE)

#save out 
write.csv(site_covariates1, "site_covariates.csv")
