library(tidyr)
library(dplyr)
library(ausplotsR)

####### Make internal indices for sites 

data <- read.csv("data/data_cleaned.csv")




#df to ouput 

#Stand Basal Area Function 
data1 <- data %>% 
    mutate(basal_area = pi * ((Diameter*0.01)/2)^2) %>% 
  group_by(Site_Name) %>% 
  summarise(sum(basal_area, na.rm = TRUE))

data1 <- data1 %>% dplyr::select(Site_Name, `sum(basal_area, na.rm = TRUE)`)
#Put into covariates   
site_covariates1 <- left_join(site_covariates, data1, by="Site_Name")
site_covariates1 <- site_covariates1 %>% 
  rename(stand_basal_area = `sum(basal_area, na.rm = TRUE)`)

#save out 
write.csv(site_covariates1, "site_covariates.csv")



