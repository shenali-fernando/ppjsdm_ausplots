library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

# Load cleaned data 
data_c <- read.csv("data/data_cleaned.csv")


####### Dominant Eucalypt == E. grandis 
grand <- data_c %>% filter(Genus_Species == "Eucalyptus grandis") %>% 
  group_by(Site_Name) %>% 
  count()
grand


### Baldy 



### Bruxner (NSW)




### Herberon 



### Koombooloomba 



### Lamb Range 



### Osullivans (NSW)