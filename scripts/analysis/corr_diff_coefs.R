library(tidyr)
library(dplyr)
library(ggplot2)
library(corrplot)

# Using df_all_3_5
df <- df_add3_5_t15


## Within-species interactions 
within <- df %>% 
  filter(species_from == species_to)

#Make wider df 
wider_df <- within %>% 
  mutate(species_site = paste0(species_from, sep = "_", site)) %>% 
  select(species_site, alpha, class_int, cc_to, georegion, region) %>% 
  pivot_wider(names_from = class_int, 
              values_from = alpha)

#######################################################################################################

############################ CORRELATION ##############################

#Canopy 
can <- wider_df %>% 
  filter(cc_to == "Canopy") %>% 
  select(species_site, small_small, small_large, large_large)

can2 <- na.omit(can) #get rid of na rows 

c <- cor(can2[, 2:4])
cor(can2$large_large, can2$small_large)
plot(can2$small_large, can2$large_large)
plot(can2$large_large, can2$small_large)


plot(can2$large_large, can2$small_small)
plot(can2$small_small, can2$large_large)

plot(can2$small_small, can2$small_large)
plot(can2$small_large, can2$small_small)

corrplot(c, 
         method = "number", 
         type = "full",
         col= colorRampPalette(c("blue","white", "red"))(25),
         title = "Canopy")

#Subcanopy 
sub <- wider_df %>% 
  filter(cc_to == "Subcanopy") %>% 
  select(species_site, small_small, small_large, large_large)

sub2 <- na.omit(sub) #get rid of na rows 

s <- cor(sub2[, 2:4])
corrplot(s, 
         method = "number", 
         type = "full",
         col= colorRampPalette(c("blue","white", "red"))(25),
         title = "Subcanopy")


plot(sub2$small_large, sub2$large_large)

####################################################################################

################### Absolute change ##################### 
#abschange = final - initial 

###Canopy 
can <- can %>% 
  mutate(abschange_ss_sl = abs(small_small - small_large)) %>% 
  mutate(abschange_ss_ll = abs(large_large - small_small)) %>% 
  mutate(abschange_sl_ll = abs(large_large - small_large))

#Small-small - small-large
mean(can$abschange_ss_sl, na.rm=T)
quantile(can$abschange_ss_sl, na.rm=T)

#Small-small - large-large
mean(can$abschange_ss_ll, na.rm=T)
quantile(can$abschange_ss_ll, na.rm=T)

#Small-large - large-large
mean(can$abschange_sl_ll, na.rm=T)
quantile(can$abschange_sl_ll, na.rm=T)



###Subcanopy
sub <- sub %>% 
  mutate(abschange_ss_sl = abs(small_small - small_large)) %>% 
  mutate(abschange_ss_ll = abs(large_large - small_small)) %>% 
  mutate(abschange_sl_ll = abs(large_large - small_large))

#Small-small - small-large
mean(sub$abschange_ss_sl, na.rm=T)
quantile(sub$abschange_ss_sl, na.rm = T)

#Small-small - large-large
mean(sub$abschange_ss_ll, na.rm=T)
quantile(sub$abschange_ss_ll, na.rm=T)

#Small-large - large-large
mean(sub$abschange_sl_ll, na.rm=T)
quantile(sub$abschange_sl_ll, na.rm=T)




################## Relative change ###################
#relchange = abschange / abs(inital)

#Canopy
can <- can %>% 
  mutate(relchange_ss_sl = abschange_ss_sl / abs(small_small)) %>% 
  mutate(relchange_ss_ll = abschange_ss_ll / abs(small_small)) %>% 
  mutate(relchange_sl_ll = abschange_sl_ll / abs(small_large))

#Small-small - small-large
mean(can$relchange_ss_sl, na.rm=T)
quantile(can$relchange_ss_sl, na.rm=T)

#Small-small - large-large
mean(can$relchange_ss_ll, na.rm=T)
quantile(can$relchange_ss_ll, na.rm=T)

#Small-large - large-large
mean(can$relchange_sl_ll, na.rm=T)
quantile(can$relchange_sl_ll, na.rm=T)

         
#Subcanopy 
sub <- sub %>% 
  mutate(relchange_ss_sl = abschange_ss_sl / abs(small_small)) %>% 
  mutate(relchange_ss_ll = abschange_ss_ll / abs(small_small)) %>% 
  mutate(relchange_sl_ll = abschange_sl_ll / abs(small_large))
         
         
#Small-small - small-large
mean(sub$relchange_ss_sl, na.rm=T)
quantile(sub$relchange_ss_sl, na.rm = T)

#Small-small - large-large
mean(sub$relchange_ss_ll, na.rm=T)
quantile(sub$relchange_ss_ll, na.rm=T)

#Small-large - large-large
mean(sub$relchange_sl_ll, na.rm=T)
quantile(sub$relchange_sl_ll, na.rm=T)


############### Using raw difference 

#Canopy 
can <- can %>% 
  mutate(diff_ss_sl = small_small - small_large) %>% 
  mutate(diff_ss_ll = large_large - small_small) %>% 
  mutate(diff_sl_ll = large_large - small_large)

#Small-small - small-large
mean(can$diff_ss_sl, na.rm=T)
quantile(can$diff_ss_sl, na.rm=T)

#Small-small - large-large
mean(can$diff_ss_ll, na.rm=T)
quantile(can$diff_ss_ll, na.rm=T)

#Small-large - large-large
mean(can$diff_sl_ll, na.rm=T)
quantile(can$diff_sl_ll, na.rm=T)



can <- can %>% 
  mutate(reldiff_ss_sl = diff_ss_sl / small_small) %>% 
  mutate(reldiff_ss_ll = diff_ss_ll / small_small) %>% 
  mutate(reldiff_sl_ll = diff_sl_ll / small_large)

#Small-small - small-large
mean(can$reldiff_ss_sl, na.rm=T)
quantile(can$reldiff_ss_sl, na.rm=T)

#Small-small - large-large
mean(can$reldiff_ss_ll, na.rm=T)
quantile(can$reldiff_ss_ll, na.rm=T)

#Small-large - large-large
mean(can$reldiff_sl_ll, na.rm=T)
quantile(can$reldiff_sl_ll, na.rm=T)


#Subcanopy
sub <- sub %>% 
  mutate(diff_ss_sl = small_small - small_large) %>% 
  mutate(diff_ss_ll = large_large - small_small) %>% 
  mutate(diff_sl_ll = large_large - small_large)

#Small-small - small-large
mean(sub$diff_ss_sl, na.rm=T)
quantile(sub$diff_ss_sl, na.rm=T)

#Small-small - large-large
mean(sub$diff_ss_ll, na.rm=T)
quantile(sub$diff_ss_ll, na.rm=T)

#Small-large - large-large
mean(sub$diff_sl_ll, na.rm=T)
quantile(sub$diff_sl_ll, na.rm=T)



sub <- sub %>% 
  mutate(reldiff_ss_sl = diff_ss_sl / small_small) %>% 
  mutate(reldiff_ss_ll = diff_ss_ll / small_small) %>% 
  mutate(reldiff_sl_ll = diff_sl_ll / small_large)

#Small-small - small-large
mean(sub$reldiff_ss_sl, na.rm=T)
quantile(sub$reldiff_ss_sl, na.rm=T)

#Small-small - large-large
mean(sub$reldiff_ss_ll, na.rm=T)
quantile(sub$reldiff_ss_ll, na.rm=T)

#Small-large - large-large
mean(sub$reldiff_sl_ll, na.rm=T)
quantile(sub$reldiff_sl_ll, na.rm=T)




         
#################################################################################

############# Correlation of change to site covariates ##########################
         
#get site covariates 
site_covariates <- read_csv("data/site_covariates.csv")


covars <- site_covariates %>% 
  rename(site = Site_Name) %>% 
  select(site, Longitude, MAP, MAT, GPP, AI, 
         FPI, stand_basal_area)

#add covars into change df 
change2 <- left_join(change, covars, by = "site")

pchange2 <- left_join(pchange, covars, by = "site")

#plot some simple scatterplots to see if theres a relationship 
ggplot(data = change2,
       aes(x = MAP, y = small_large)) + 
  geom_point(aes(colour = region)) +
  geom_smooth(method = "loess") +
  facet_wrap(~cc_from) +
  theme_bw()



##need to understand whats going on 
change_c <- change2 %>% filter(cc_from == "Canopy")
change_s <- change2 %>% filter(cc_from == "Subcanopy")


mod <- lm(data = change_s, 
   D_SS ~ FPI)
summary(mod)
plot(mod)


## Correlation test 
cor.test(~D_SS + MAP, data = change_s)



#No pattern #Noyechange_s#No pattern #NoyearlymaxNDVI pattern 



############################################################################

#Calculate percentage change between class interactions for FGs 
within <- df %>% 
  filter(species_from == species_to)  

within <- within %>% 
  select(species_from, alpha, site, class_int, region, georegion, cc_from)


within_change <- within %>% 
  pivot_wider(names_from = class_int, 
              values_from = alpha)


pchange <- within_change %>% 
  mutate(S = ((small_large - small_small)/small_small) * 100) %>% 
  mutate(D_SS = ((large_large - small_small)/small_small) * 100) %>% 
  mutate(D_SL = ((large_large - small_large)/small_large) * 100)


pchange_c <- pchange %>% filter(cc_from == "Canopy")
pchange_s <- pchange %>% filter(cc_from == "Subcanopy")

median(pchange_c$D_SS, na.rm = TRUE)
median(pchange_c$D_SL, na.rm = TRUE)
median(pchange_c$S, na.rm = TRUE)

median(pchange_s$D_SS, na.rm = TRUE)
median(pchange_s$D_SL, na.rm = TRUE)
median(pchange_s$S, na.rm = TRUE)
