library(terra)
library(sp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(car)

#Load covariate data 
covars <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/data/site_covariates.csv")

#need to rename some cols 
covars <- covars %>% 
  rename(MAP = MAP..mm.) %>% 
  rename(elev = elevation..m.) %>% 
  rename(FPI = FPI..3PG.)

#Need to understand the correlation between MAT, MAP, and productivity index (FPI) 
#So, create a scatterplot matrix 
pairs(~MAT + MAP + elev + FPI, data = covars)

#Don't actually think elevation is strongly correlated with the other covariates; so
pairs(~MAT + MAP + FPI, data = covars)

#can add some variance envelopes 
scatterplotMatrix(~MAT + MAP + FPI, data = covars,
                  diagonal = FALSE,             # Remove kernel density estimates
                  regLine = FALSE,      # Linear regression line width
                  smooth = list(col.smooth = "red",   # Non-parametric mean color
                                col.spread = "blue")) # Non-parametric variance color


## Run some correlations one by one 
#MAT v MAP 
ggplot(data = covars, aes(x = MAT, y = MAP)) +
  geom_point() + # Show dots
  geom_text(aes(label = Site_Name),
    nudge_x = 0.25, nudge_y = 0.75, 
    check_overlap = F) + 
  theme_bw()

#MAT v FPI 
ggplot(data = covars, aes(x = MAT, y = FPI)) +
  geom_point() + # Show dots
  geom_text(aes(label = Site_Name),
            nudge_x = 0.25, nudge_y = 0.25, 
            check_overlap = F) + 
  geom_vline(xintercept = 11, colour = "red", linetype = "dashed", linewidth = 1) + 
  theme_bw()

#MAP v FPI
ggplot(data = covars, aes(x = MAP, y = FPI)) +
  geom_point() + # Show dots
  geom_text(aes(label = Site_Name),
            nudge_x = 0.25, nudge_y = 0.25, 
            check_overlap = F) + 
  theme_bw()

##Correlation test 

cor <- stats::cor(covars[, 4:7])
corrplot::corrplot(cor, method = "number", order = 'hclust') #not a straight corrleation betwen MAT and FPI anyways



