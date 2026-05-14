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
wa1 <- rast("data/ETa/ETa_wa1.ET.tif")
plot(wa1)
value <- terra::extract(wa1, point)
sites$a <- value$ETa_wa1.ET

wa2 <- rast("data/ETa/ETa_wa2.ET.tif")
plot(wa2)
value <- terra::extract(wa2, point)
sites$b <- value$ETa_wa2.ET

wa3 <- rast("data/ETa/ETa_wa3.ET.tif")
plot(wa3)
value <- terra::extract(wa3, point)
sites$c <- value$ETa_wa3.ET

#VIC
vic1 <- rast("data/ETa/ETa_vic1.ET.tif")
plot(vic1)
value <- terra::extract(vic1, point)
sites$d <- value$ETa_vic1.ET

vic2 <- rast("data/ETa/ETa_vic2.ET.tif")
plot(vic2)
value <- terra::extract(vic2, point)
sites$e <- value$ETa_vic2.ET

#SE NSW
sensw1 <- rast("data/ETa/ETa_sensw1.ET.tif")
plot(sensw1)
value <- terra::extract(sensw1, point)
sites$f <- value$ETa_sensw1.ET

#TAS 
tas1 <- rast("data/ETa/ETa_tas1.ET.tif")
plot(tas1)
value <- terra::extract(tas1, point)
sites$g <- value$ETa_tas1.ET

tas2 <- rast("data/ETa/ETa_tas2.ET.tif")
plot(tas2)
value <- terra::extract(tas2, point)
sites$h <- value$ETa_tas2.ET

tas3 <- rast("data/ETa/ETa_tas3.ET.tif")
plot(tas3)
value <- terra::extract(tas3, point)
sites$i <- value$ETa_tas3.ET

tas4 <- rast("data/ETa/ETa_tas4.ET.tif")
plot(tas4)
value <- terra::extract(tas4, point)
sites$j <- value$ETa_tas4.ET

tas5 <- rast("data/ETa/ETa_tas5.ET.tif")
plot(tas5)
value <- terra::extract(tas5, point)
sites$k <- value$ETa_tas5.ET

tas6 <- rast("data/ETa/ETa_tas6.ET.tif")
plot(tas6)
value <- terra::extract(tas6, point)
sites$l <- value$ETa_tas6.ET

#QLD 
qld1 <- rast("data/ETa/ETa_qld1.ET.tif")
plot(qld1)
value <- terra::extract(qld1, point)
sites$m <- value$ETa_qld1.ET

qld2 <- rast("data/ETa/ETa_qld2.ET.tif")
plot(qld2)
value <- terra::extract(qld2, point)
sites$n <- value$ETa_qld2.ET

#N NSW
nnsw1 <- rast("data/ETa/ETa_nnsw1.ET.tif")
plot(nnsw1)
value <- terra::extract(nnsw1, point)
sites$o <- value$ETa_nnsw1.ET

nnsw2 <- rast("data/ETa/ETa_nnsw2.ET.tif")
plot(nnsw2)
value <- terra::extract(nnsw2, point)
sites$p <- value$ETa_nnsw2.ET

nnsw3 <- rast("data/ETa/ETa_nnsw3.ET.tif")
plot(nnsw3)
value <- terra::extract(nnsw3, point)
sites$aa <- value$ETa_nnsw3.ET

nnsw4 <- rast("data/ETa/ETa_nnsw4.ET.tif")
plot(nnsw4)
value <- terra::extract(nnsw4, point)
sites$bb <- value$ETa_nnsw4.ET



## Combine everything together 
sites <- sites %>% 
  rowwise() %>% 
  mutate(ETa = mean(c_across(a:bb), na.rm = TRUE)) %>% 
  ungroup()

#Cool, add into final site covariates df 
site_covariates1 <- merge(x = site_covariates, y = sites[, c("Site_Name", "ETa")], by = "Site_Name", all.x = TRUE)

#save out 
write.csv(site_covariates1, "site_covariates.csv")
