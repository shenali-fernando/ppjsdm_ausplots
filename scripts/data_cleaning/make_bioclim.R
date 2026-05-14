library(dplyr)
library(terra)
library(geodata)
library(corrplot)

#Load covariate data 
data <- read.csv("data/data_cleaned.csv")

sites <- data %>% dplyr::select(Site_Name, Latitude, Longitude) %>% 
  unique()
rownames(sites) <- NULL


#make df of coordinates
co <- data.frame(lon = sites$Longitude,
                 lat = sites$Latitude) #make a df with just long and lat 

sites_vect <- vect(co, geom=c("lon", "lat")) #make spatvector

#Download historic worldclim data (1970–2000 baseline)
bioclim <- worldclim_global(var = "bio", res = 2.5, path = tempdir())

vals <- terra::extract(bioclim, sites_vect) #extract climate vals 

vals <- cbind(sites$Site_Name, vals)

vals <- vals %>% dplyr::select(-ID)
colnames(vals) <- c("Site_Name", 
                      "AnnualMeanTemp",
                      "MeanDiurnalRange",
                      "Isothermality",
                      "TempSeasonality",
                      "MaxTempWarmestMonth",
                      "MinTempColdestMonth",
                      "TempAnnualRange",
                      "MeanTempWettestQuarter",
                      "MeanTempDriestQuarter",
                      "MeanTempWarmestQuarter",
                      "MeanTempColdestQuarter",
                      "AnnualPrecipitation",
                      "PrecipWettestMonth",
                      "PrecipDriestMonth",
                      "PrecipSeasonality",
                      "PrecipWettestQuarter",
                      "PrecipDriestQuarter",
                      "PrecipWarmestQuarter",
                      "PrecipColdestQuarter"
                    )

#Save this out 
write.csv(vals, "site_bioclim.csv")



#Correlation analysis to see which ones can be kicked out 
df <- site_bioclim1[1:5]
pairs(df)

#Selection of uncorrelated vars based on González‐Orozco, C. E., 
#Thornhill, A. H., Knerr, N., Laffan, S., & Miller, J. T. (2014). 
#Biogeographical regions and phytogeography of the eucalypts. Diversity and 
#Distributions, 20(1), 46-58.
site_bioclim1 <- vals %>% 
  dplyr::select(AnnualMeanTemp, TempSeasonality, AnnualPrecipitation,
                PrecipSeasonality, PrecipColdestQuarter)
