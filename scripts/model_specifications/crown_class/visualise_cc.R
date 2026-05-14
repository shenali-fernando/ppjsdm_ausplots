library(ppjsdm)
library(tidyr)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

### Load csvs of coefficients 
df_pil <- read.csv("code/crownclass/piluaris_cc.csv")
df_reg <- read.csv("code/crownclass/regnans_cc.csv")
df_grandis <- read.csv("code/crownclass/grangis_cc.csv")
df_diver <- read.csv("code/crownclass/westernaus_cc.csv")
df_ob <- read.csv("code/crownclass/obliqua_cc.csv")
df_fast <- read.csv("code/crownclass/fastigata_cc.csv")
df_del <- read.csv("code/crownclass/delegatensis_cc.csv")


#rbind all dfs into a single one
df_cc <- rbind(df_pil, df_reg, df_grandis, df_diver, df_ob, df_fast, df_del)

df_cc %>% count(Fit) #all Sites are present 


df_cc <- df_cc %>% rename(Site = Fit)  #rename fits col

#need to create some new columns that give what class is used: this will make subsetting easier
df_cc <- df_cc %>% 
  mutate(class_from = gsub("^[^ ]+ ", "\\1", from)) %>% 
  mutate(class_to = gsub("^[^ ]+ ", "\\1", to)) %>% #same for euc or non-euc 
  mutate(species_from = gsub("^([^ ]+).*", "\\1", from)) %>% 
  mutate(species_to = gsub("^([^ ]+).*", "\\1", to))

df_cc <- df_crownclass %>% 
  mutate(class_to = if_else(class_to == "regnans Suppressed", "Suppressed", class_to)) %>% 
  mutate(class_to = if_else(class_to == "regnans Intermediate", "Intermediate", class_to)) %>% 
  mutate(class_to = if_else(class_to == "regnans Co-dominant", "Co-dominant", class_to)) %>% 
  mutate(class_to = if_else(class_to == "regnans Dominant", "Dominant", class_to)) %>% 
  mutate(class_to = if_else(class_to == "Dominannt", "Dominant", class_to)) %>% 
  mutate(class_to = if_else(class_to == "Codominant", "Co-dominant", class_to)) %>% 
  mutate(class_from = if_else(class_from == "regnans Suppressed", "Suppressed", class_from)) %>% 
  mutate(class_from = if_else(class_from == "regnans Intermediate", "Intermediate", class_from)) %>% 
  mutate(class_from = if_else(class_from == "regnans Co-dominant", "Co-dominant", class_from)) %>% 
  mutate(class_from = if_else(class_from == "regnans Dominant", "Dominant", class_from))

df_cc %>% count(class_to)
df_cc %>% count(class_from)
## save the df at this point

df_crownclass <- write.csv(df_cc, "output/df_crownclass.csv")

df_cc <- read.csv("output/df_crownclass.csv")

df_cc %>% count(Site)

auc_df_final <- auc_df_final %>% 
  mutate(dom_euc = site) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Carey", "Dombakup", "Warren", "Dawson", "Giants", "Sutton", "Frankland", "Clare", "Collins"), 
                           "E. diversicolor", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Lardner", "NorthStyx", "Turtons", "Weeaproinah", "Weld"), 
                           "E. regnans", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Newline", "WaratahMix", "WogWay", "Goodenia"), 
                           "E. fastigata", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank"), 
                           "E. pilularis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Bruxner", "Osullivans", "Baldy", "Koombooloomba", "Lamb Range", "Herberton"), 
                           "E. grandis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Flowerdale", "Dip", "Bird", "Supersite", "ZigZag", "BlackRiver", "BondTier", "Candelo"), 
                           "E. obliqua", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("BenRidge", "Caveside", "MtMaurice", "Mackenzie", "MtField"), 
                           "E. delegatensis", dom_euc))


df_crownclass <- df_crownclass %>% 
  mutate(range = hi - lo) %>% 
  filter(range >5)

df_crownclass %>% group_by(range) %>% count(Site)


#filter dataframe 
df_euc_cc <- df_cc %>% 
  filter(species_from == "Eucalyptus") %>% 
  filter(species_to == "Eucalyptus") %>% 
  arrange(dom_euc, alpha) %>% 
  mutate(Site = factor(Site, levels = unique(Site))) 


d <- df_euc_cc %>% 
  filter(class_int == "Dominant_Dominant")%>% 
  arrange(dom_euc, alpha) %>% 
  mutate(Site = factor(Site, levels = unique(Site))) 

c <- df_euc_cc %>% 
  filter(class_int == "Co-dominant_Co-dominant")%>% 
  arrange(dom_euc, alpha) %>% 
  mutate(Site = factor(Site, levels = unique(Site))) 

i <- df_euc_cc %>% 
  filter(class_int == "Intermediate_Intermediate")%>% 
 # filter(! Site == "BirdTree") %>% 
  arrange(dom_euc, alpha) %>% 
  mutate(Site = factor(Site, levels = unique(Site))) 

s <- df_euc_cc %>% 
  filter(class_int == "Suppressed_Suppressed")%>% 
  arrange(dom_euc, alpha) %>% 
  mutate(Site = factor(Site, levels = unique(Site))) 


## making box-plots for each 
## Dominant 
dominant <- ggplot(d, aes(x = alpha, y = Site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "Site",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Dominant") +
  theme_bw()

## Codominant
codominant <- ggplot(c, aes(x = alpha, y = Site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Codominant") +
  theme_bw()

#bruxner, lorne, lambrange - dom uncertain, but eng
# c is no strong association? lumping = small effects 

#scatterplot - c to c/d, d to c/d, c and d, 
#with intervals
#plot with site numbers


### TASSIE DELEGATENSIS IS DIFFERent 

## Intermediate

intermediate <- 
ggplot(i, aes(x = alpha, y = Site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Intermediate") +
  theme_bw()


## Suppressed
suppressed <- ggplot(s, aes(x = alpha, y = site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Suppressed") +
  theme_bw()





library(ggpubr)

ggarrange(dominant, codominant, intermediate, suppressed, 
          ncol = 2, nrow = 2, 
          common.legend = TRUE, 
          legend = "right")



##### BETWEEN GROUPS 


### S v C/D
s_cd <- full_df_com %>% 
  filter(species_from == "Eucalyptus") %>% 
  filter(species_to == "Eucalyptus") %>% 
  mutate(class_int = paste(class_from, class_to)) %>% 
  filter(class_int %in% c("Suppressed Co/dominant", "Co/dominant Suppressed"))

ggplot(s_cd, aes(x = alpha, y = site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Suppressed") +
  theme_bw()


### I v C/D
i_cd <- full_df_com %>% 
  filter(species_from == "Eucalyptus") %>% 
  filter(species_to == "Eucalyptus") %>% 
  mutate(class_int = paste(class_from, class_to)) %>% 
  filter(class_int %in% c("Intermediate Co/dominant", "Co/dominant Intermediate"))

ggplot(i_cd, aes(x = alpha, y = site), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  ggtitle("Suppressed") +
  theme_bw()



### Non-euc S v C/D
nons_cd <- full_df_com %>% 
  mutate(class_int = paste(from, to)) %>% 
  filter(class_int %in% c("Non-euc Suppressed Eucalyptus Co/dominant", "Eucalyptus Co/dominant Non-euc Suppressed"))

ggplot(nons_cd, aes(x = alpha, y = site, colour = dom_euc), size = 1.5) + 
  geom_boxplot() + 
  geom_errorbar(aes(xmin = lo, xmax = hi)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Alpha Coefficient (+/-95% CI)",
       y = "",
       color = "Dominant Eucalyptus Species") + 
  scale_y_discrete(limits = rev) + #more making order right
  theme_bw()





### 
s_c <- euc_cc %>% 
  filter(class_int == "Suppressed_Co-dominant")%>% 
  arrange(dom_euc, alpha) %>% 
  mutate(Site = factor(Site, levels = unique(Site))) 
match



###############################################################
################### Ind v Coef ################################

inds_full %>% filter(to == "Non-euc Dominant") %>% 
  filter(n <=15)
#non-eucs are set to less than 12 right now, but only 4 sites have 15 and below so might look into this to make it easier to explain 
#want to very simple rule for non-eucs 

#add a column for range of the ci 
full_df_com <- full_df_com %>% 
  mutate(range_ci = hi - lo) 


large_ci <- full_df_com %>% 
  filter(range_ci > 4)

full_df_com %>% 
  
ggplot(data = large_ci) + 
  geom_point(aes(x = ind_to, y = range_ci, colour = class_to)) + 
  ylim(c(0, 75))





###### SCATTERPLOTS

### C v D 
new_full_df <- new_full_df %>% 
  mutate(class_int = paste(class_from, class_to))

c_d <- new_full_df %>% 
  filter(species_from == "Eucalyptus") %>% 
  filter(species_to == "Eucalyptus") %>% 
  filter(class_int %in% c("Dominant Dominant", "Co-dominant Co-dominant"))

new_c_d <- c_d %>% 
  select(site, class_from, alpha, lo, hi)

new_c <- new_c_d %>% 
  filter(class_from == "Co-dominant") %>% 
  rename(c = alpha) %>% 
  rename(c_lo = lo) %>% 
  rename(c_hi = hi) %>% 
  select(-class_from)

new_d <- new_c_d %>% 
  filter(class_from == "Dominant") %>% 
  rename(d = alpha) %>% 
  rename(d_lo = lo) %>% 
  rename(d_hi = hi) %>% 
  select(-class_from)

c_d <- full_join(new_c, new_d, by = 'site')

new_c %>% count(site)

cvd <- ggplot(data = c_d) + 
  geom_point(aes(x = c, y = d)) + 
  geom_errorbarh(aes(xmin = c_lo, xmax = c_hi, y = d), 
                 linewidth = 0.1, colour = "gray50") +  
  geom_errorbar(aes(ymin = d_lo, ymax = d_hi, x = c),
                linewidth = 0.1, colour = "gray50") + 
  geom_vline(xintercept = 0, linewidth = 0.75) + 
  geom_hline(yintercept = 0, linewidth = 0.75) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.75) + 
  theme_bw() + 
  xlab("Codominant") + ylab("Dominant") +
  ylim(-5, 1) + 
  xlim(-5, 1)


cvd


### D v C/D

full_df_com <- full_df_com %>% 
  mutate(class_int = paste(class_from, class_to))

cd <- full_df_com %>% 
  filter(species_from == "Eucalyptus") %>% 
  filter(species_to == "Eucalyptus") %>% 
  filter(class_int == "Co/dominant Co/dominant")

cd <- cd %>% 
  select(site, class_from, alpha, lo, hi)

cd <- cd %>% 
  rename(cd = alpha) %>% 
  rename(cd_lo = lo) %>% 
  rename(cd_hi = hi) %>% 
  select(-class_from)


cd_com <- full_join(cd, new_d, by = 'site')


cd_com_d <- ggplot(data = cd_com) + 
  geom_point(aes(x = d, y = cd)) + 
  geom_errorbarh(aes(xmin = d_lo, xmax = d_hi, y = cd), 
                 linewidth = 0.1, colour = "gray50") +  
  geom_errorbar(aes(ymin = cd_lo, ymax = cd_hi, x = d),
                linewidth = 0.1, colour = "gray50") + 
  geom_vline(xintercept = 0, linewidth = 0.75) + 
  geom_hline(yintercept = 0, linewidth = 0.75) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.75) + 
  theme_bw() + 
  xlab("Dominant") + ylab("Co/Dominant")

cd_com_d



### C v C/D

cd_com_c <- full_join(cd, new_c, by = 'site')

cd_com_c <- ggplot(data = cd_com_c) + 
  geom_point(aes(x = c, y = cd)) + 
  geom_errorbarh(aes(xmin = c_lo, xmax = c_hi, y = cd), 
                 linewidth = 0.1, colour = "gray50") +  
  geom_errorbar(aes(ymin = cd_lo, ymax = cd_hi, x = c),
                linewidth = 0.1, colour = "gray50") + g
  geom_vline(xintercept = 0, linewidth = 0.75) + 
  geom_hline(yintercept = 0, linewidth = 0.75) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.75) + 
  xlab("Codominant") + ylab("Co/Dominant") +
  theme_bw() 

cd_com_c










## Compare WITHIN Coef by site 

full_cd_coom <- full_join(cd_com_c, new_d, by = 'site')

final_sizemod <- final_sizemod %>% 
  mutate(dom_euc = site) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Carey", "Dombakup", "Warren", "Dawson", "Giants", "Sutton", "Frankland", "Clare", "Collins"), 
                           "E. diversicolor", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Lardner", "NorthStyx", "Turtons", "Weeaproinah", "Weld"), 
                           "E. regnans", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Newline", "WaratahMix", "WogWay", "Goodenia"), 
                           "E. fastigata", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank"), 
                           "E. pilularis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Bruxner", "Osullivans", "Baldy", "Koombooloomba", "Lamb Range", "Herberton"), 
                           "E. grandis", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("Flowerdale", "Dip", "Bird", "Supersite", "ZigZag", "BlackRiver", "BondTier", "Candelo"), 
                           "E. obliqua", dom_euc)) %>% 
  mutate(dom_euc = if_else(dom_euc %in% c("BenRidge", "Caveside", "MtMaurice", "Mackenzie", "MtField"), 
                           "E. delegatensis", dom_euc))

  
  ggplot(full_cd_coom) + 
    geom_boxplot(aes(x = cd, y = site), colour = "orange") + 
    geom_errorbar(aes(y = site, xmin = cd_lo, xmax = cd_hi), colour = "orange") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
    geom_boxplot(aes(x = c, y = site),  colour = "#56AFD6",) + 
    geom_errorbar(aes(y = site, xmin = c_lo, xmax = c_hi, x = c),  colour = "#56AFD6") + 
    geom_boxplot(aes(x = d, y = site),  colour = "darkblue") + 
    geom_errorbar(aes(y = site, xmin = d_lo, xmax = d_hi, x = d), colour = "darkblue") + 
    xlim(-10, 5)


  codom2
  