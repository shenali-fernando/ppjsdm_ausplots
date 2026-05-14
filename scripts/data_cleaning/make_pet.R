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
wa1 <- rast("data/PET/pet_wa1.PET.tif")
plot(wa1)
value <- terra::extract(wa1, point)
sites$a <- value$pet_wa1.PET

wa2 <- rast("data/PET/pet_wa2.PET.tif")
plot(wa2)
value <- terra::extract(wa2, point)
sites$b <- value$pet_wa2.PET

#VIC 
vic1 <- rast("data/PET/pet_vic1.PET.tif")
plot(vic1)
value <- terra::extract(vic1, point)
sites$c <- value$pet_vic1.PET

vic2 <- rast("data/PET/pet_vic2.PET.tif")
plot(vic2)
value <- terra::extract(vic2, point)
sites$d <- value$pet_vic2.PET

vic3 <- rast("data/PET/pet_vic3.PET.tif")
plot(vic3)
value <- terra::extract(vic3, point)
sites$e <- value$pet_vic3.PET


# SE NSW 
sensw1 <- rast("data/PET/pet_sensw1.PET.tif")
plot(sensw1)
value <- terra::extract(sensw1, point)
sites$f <- value$pet_sensw1.PET


# N NSW 
nnsw1 <- rast("data/PET/pet_nnsw1.PET.tif")
plot(nnsw1)
value <- terra::extract(nnsw1, point)
sites$g <- value$pet_nnsw1.PET

nnsw2 <- rast("data/PET/pet_nnsw2.PET.tif")
plot(nnsw2)
value <- terra::extract(nnsw2, point)
sites$h <- value$pet_nnsw2.PET

nnsw3 <- rast("data/PET/pet_nnsw3.PET.tif")
plot(nnsw3)
value <- terra::extract(nnsw3, point)
sites$i <- value$pet_nnsw3.PET


#QLD 
qld1 <- rast("data/PET/pet_qld1.PET.tif")
plot(qld1)
value <- terra::extract(qld1, point)
sites$j <- value$pet_qld1.PET


#TAS 
tas1 <- rast("data/PET/pet_tas1.PET.tif")
plot(tas1)
value <- terra::extract(tas1, point)
sites$k <- value$pet_tas1.PET

tas2 <- rast("data/PET/pet_tas2.PET.tif")
plot(tas2)
value <- terra::extract(tas2, point)
sites$l <- value$pet_tas2.PET

tas3 <- rast("data/PET/pet_tas3.PET.tif")
plot(tas3)
value <- terra::extract(tas3, point)
sites$m <- value$pet_tas3.PET

tas4 <- rast("data/PET/pet_tas4.PET.tif")
plot(tas4)
value <- terra::extract(tas4, point)
sites$n <- value$pet_tas4.PET

tas5 <- rast("data/PET/pet_tas5.PET.tif")
plot(tas5)
value <- terra::extract(tas5, point)
sites$o <- value$pet_tas5.PET



## Combine everything together 
sites <- sites %>% 
  rowwise() %>% 
  mutate(PET = mean(c_across(a:o), na.rm = TRUE)) %>% 
  ungroup()


#Cool, add into final site covariates df 
site_covariates1 <- merge(x = site_covariates, y = sites[, c("Site_Name", "PET")], by = "Site_Name", all.x = TRUE)

#save out 
write.csv(site_covariates1, "site_covariates.csv")
