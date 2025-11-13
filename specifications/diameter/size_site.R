library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(patchwork)
library(scales)
library(austraits)
library(forcats)



###### ALL SPECIES - MODELLING INTERACTIONS BETWEEN SPECIES_SIZE GROUPS 
# Have some pretty good functions to help us here 
source("size_funs.R")

# It is important to have well-parameterised models... 
# From previous work on the dataset, I know the ballpark of where the parameters should be 
# Let's run the fit_opt function for all sites 

data <- read.csv("data/data_cleaned.csv")

sites <- unique(data$Site_Name)

full_df <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                           "hi_numerical", "Potential","site")
for (i in sites){
  
  site_mod <- size_sites(site = i, #single site 
               group_type = "species_size",
               show_size_freq = FALSE,
               config_only = FALSE, #if TRUE, returns config only and exits function
               threshold = 13, 
               short_range = 10, 
               short_model = "exponential")
  
  working_df <- make_sum_df(fits = list(site_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(site_mod$sum))
  
  working_df <- working_df %>% mutate(site = i)                              
                            
  full_df<- rbind(full_df, working_df)
}


full_df<- full_df %>% 
  mutate(range_ci = hi - lo)

a <- full_df %>% filter(range_ci > 5)
a %>% count(site)

#By site cleaning 
birdtree <-  size_sites(site = "BirdTree", #single site 
                       group_type = "species_size",
                       show_size_freq = FALSE,
                       config_only = FALSE, #if TRUE, returns config only and exits function
                       threshold = 14, 
                       short_range = 10, 
                       short_model = "exponential")

birdtree_df <- make_sum_df(fits = list(fit), 
            summ = list(sum))

birdtree_df <- birdtree_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "BirdTree")


birdtree_df %>% 
  filter(range_ci > 5)


bird <-  size_sites(site = "Bird", #single site 
                        group_type = "species_size",
                        show_size_freq = FALSE,
                        config_only = TRUE, #if TRUE, returns config only and exits function
                        threshold = 13, 
                        short_range = 10, 
                        short_model = "exponential")

bird_df <- make_sum_df(fits = list(fit), 
                           summ = list(sum))
bird_df <- bird_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Bird")

bird_df %>% 
  filter(range_ci > 5)


brux <- size_sites(site = "Bruxner", #single site 
                   group_type = "species_size",
                   show_size_freq = FALSE,
                   config_only = FALSE, #if TRUE, returns config only and exits function
                   threshold = 14, 
                   short_range = 10, 
                   short_model = "exponential")
brux_df <- make_sum_df(fits = list(brux$fit), 
                       summ = list(brux$sum))

brux_df <- brux_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Bruxner")

brux_df %>% 
  filter(range_ci > 5)



cave <- size_sites(site = "Caveside", #single site 
                   group_type = "species_size",
                   show_size_freq = FALSE,
                   config_only = TRUE, #if TRUE, returns config only and exits function
                   threshold = 16, 
                   short_range = 10, 
                   short_model = "exponential")
cave_df <- make_sum_df(fits = list(fit), 
                       summ = list(sum))
cave_df <- cave_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Caveside")

cave_df %>% 
  filter(range_ci > 5)


black <- size_sites(site = "BlackRiver", #single site 
                   group_type = "species_size",
                   show_size_freq = FALSE,
                   config_only = FALSE, #if TRUE, returns config only and exits function
                   threshold = 13, 
                   short_range = 12, 
                   short_model = "exponential")
black_df <- make_sum_df(fits = list(fit), 
                       summ = list(sum))
black_df <- black_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "BlackRiver")

black_df %>% 
  filter(range_ci > 5)



lardner <- size_sites(site = "Lardner", #single site 
                     group_type = "species_size",
                     show_size_freq = FALSE,
                     config_only = FALSE, #if TRUE, returns config only and exits function
                     threshold = 14, 
                     short_range = 10, 
                     short_model = "exponential")
lardner_df <- make_sum_df(fits = list(fit), 
                         summ = list(sum))
lardner_df <- lardner_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Lardner")

lardner_df %>% 
  filter(range_ci > 5)



mtf <- size_sites(site = "MtField", #single site 
                   group_type = "species_size",
                   show_size_freq = FALSE,
                   config_only = FALSE, #if TRUE, returns config only and exits function
                   threshold = 14, 
                   short_range = 10, 
                   short_model = "exponential")
mtf_df <- make_sum_df(fits = list(fit), 
                       summ = list(sum))
mtf_df <- mtf_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "MtField")

mtf_df %>% 
  filter(range_ci > 5)


north <- size_sites(site = "NorthStyx", #single site 
                  group_type = "species_size",
                  show_size_freq = FALSE,
                  config_only = FALSE, #if TRUE, returns config only and exits function
                  threshold = 16, 
                  short_range = 10, 
                  short_model = "exponential")
north_df <- make_sum_df(fits = list(fit), 
                      summ = list(sum))
north_df <- north_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "NorthStyx")

north_df %>% 
  filter(range_ci > 5)


waratah <- size_sites(site = "WaratahMix", #single site 
                    group_type = "species_size",
                    show_size_freq = FALSE,
                    config_only = FALSE, #if TRUE, returns config only and exits function
                    threshold = 22, 
                    short_range = 10, 
                    short_model = "exponential")
waratah_df <- make_sum_df(fits = list(waratah$fit), 
                        summ = list(waratah$sum))
waratah_df <- waratah_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "WaratahMix")

waratah_df %>% 
  filter(range_ci > 5)


tine <- size_sites(site = "Tinebank", #single site 
                      group_type = "species_size",
                      show_size_freq = FALSE,
                      config_only = FALSE, #if TRUE, returns config only and exits function
                      threshold = 14, 
                      short_range = 10, 
                      short_model = "exponential")
tine_df <- make_sum_df(fits = list(tine$fit), 
                          summ = list(tine$sum))
tine_df <- tine_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Tinebank")

tine_df %>% 
  filter(range_ci > 5)


flower <- size_sites(site = "Flowerdale", #single site 
                   group_type = "species_size",
                   show_size_freq = FALSE,
                   config_only = FALSE, #if TRUE, returns config only and exits function
                   threshold = 18, 
                   short_range = 10, 
                   short_model = "exponential")
flower_df <- make_sum_df(fits = list(flower$fit), 
                       summ = list(flower$sum))
flower_df <- flower_df %>% 
  mutate(range_ci = hi - lo) %>% 
  mutate(site = "Flowerdale")

flower_df %>% 
  filter(range_ci > 5)


###### Save out df 
full_df <- full_df %>%  filter(!site %in% c("BirdTree","Bird","Bruxner", "BlackRiver", "Caveside",
                                            "Flowerdale", "Lardner", "MtField", "NorthStyx", "Tinebank", "WaratahMix"))

full_df <- rbind(full_df, birdtree_df, bird_df, brux_df, black_df, cave_df, flower_df, lardner_df, mtf_df, north_df, tine_df, waratah_df)
full_df %>% count(site)

full_df %>% filter(range_ci > 5)

write.csv(full_df, "size_site_df_add3.5.csv")


### Okay! Have a dataframe, now we can visualise stuff 
#Before that let's make a few new columns 

full_df2 <- full_df %>% 
  mutate(class_from = str_extract(from, "\\w+$")) %>% 
  mutate(class_to = str_extract(to, "\\w+$")) %>% 
  mutate(species_from = str_extract(from, "\\w+\\s+\\w+")) %>% 
  mutate(species_to = str_extract(to, "\\w+\\s+\\w+"))

## Okay, now can subset the df and visualise different things 
#We want to get an understanding of the range of interactions going on 

full_df2 <- full_df2 %>% mutate(class_int = paste(class_from, sep = "_", class_to))

full_df2 <- full_df2 %>% 
  mutate(georegion = case_when(site %in% c("Weeaproinah", "Turtons", "Lardner") ~ "S_VIC", 
                               site %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek") ~ "N_VIC", 
                               site %in% c("Dawson", "Frankland", "Clare", "Giants") ~ "S_WA",
                               site %in% c("Carey", "Dombakup", "Warren",  "Sutton","Collins") ~ "N_WA", 
                               site %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton") ~ "QLD", 
                               site %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans") ~ "N_NSW", 
                               site %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo") ~ "S_NSW", 
                               site %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx") ~ "d_TAS", 
                               site %in% c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Flowerdale", "Dip") ~ "o_TAS")) 




full_df2 <- full_df2 %>% 
  mutate(region = case_when(georegion %in% c("S_VIC", "N_VIC", "S_NSW") ~ "SE_AUS", 
                            georegion %in% c("S_WA", "N_WA") ~ "WA",
                            georegion %in% c("QLD", "N_NSW") ~ "N_AUS",  
                            georegion %in% c("o_TAS", "d_TAS") ~ "TAS")

)

#save
write.csv(full_df2, "size_site_df_add3.5.csv")

within %>% group_by(class_int) %>% summarise(median = median(alpha), 
                                             mean = mean(alpha), 
                                             min = min(alpha), 
                                             max = max(alpha))

between %>% group_by(class_int) %>% summarise(median = median(alpha), 
                                              mean = mean(alpha), 
                                              min = min(alpha), 
                                              max = max(alpha))



## Within-species interactions 
within <- df %>% 
  filter(species_from == species_to)

within <- within %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA))




within <- within %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA)) %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
    group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
    group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
    group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
    TRUE ~ group))

within$class_int <- factor(within$class_int,
                           levels = c("small_small", "small_large", "large_large"))


within %>% 
  mutate(class_int = case_when(
    class_int == "small_small" ~ "Small ↔ Small", 
    class_int == "small_large" ~ "Small ↔ Large", 
    class_int == "large_large" ~ "Large ↔ Large")) %>% 
ggplot(aes(x = alpha, y = class_int, colour = cc_to)) + 
  geom_violin(aes(colour = cc_to), 
              width = 0.85, 
              position = position_dodge(0.7)) +
  geom_boxplot(aes(colour = cc_to), 
               position = position_dodge(0.7),
               outlier.shape = NA, 
               width = 0.3) +
  geom_point(aes(group = cc_to, 
                 colour = cc_to), 
             position=position_jitterdodge(seed = 21, 
                                           jitter.width = 0.3,
                                           dodge.width = 0.7), 
             shape = 19, 
             size = 1.5, 
             alpha = 0.35) +
  stat_summary(aes(fill = cc_to), fun = mean,
               geom = "point",
               shape = 19,
               size = 3, 
               color="black",
               position = position_dodge(width=0.7, preserve = "single"), 
               show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw(base_size = 13) + 
  labs(x = "Interaction Coefficient", y = "") +
  ggtitle("Within-species interactions") +
  scale_colour_manual(values = c("orange", "orchid4"), 
                      name = "Strata")
  ylab("") 



  unique(within$species_from)


naus <- within %>% 
  filter(region == "N_AUS") %>% 
ggplot(aes(x = alpha, 
           y = group, 
           colour = georegion, fill = fill_col)) + 
 geom_point(shape = 21, size = 3) + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#E78ECD","#28BBD7"), 
                     guide = "none") +
  scale_fill_manual(values = c("#E78ECD","#28BBD7"), 
                    na.value = "white") +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  ylab("N AUS")
  

se_aus <-  within %>% 
  filter(region %in% c("TAS", "SE_AUS")) %>% 
  group_by(species_from, class_int) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(obs > 2) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col
            )) + 
  geom_point(shape = 21, size = 3) + 
  scale_color_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                     guide = "none") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                    na.value = "white") +
  facet_wrap(~class_int) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()  + 
  ylab("SE AUS")




wa <- within %>% 
  filter(region == "WA") %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col)) + 
  geom_point(shape = 21, size = 3) + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#C0AB52",  "#4FBF85"), guide = "none") + 
  scale_fill_manual(values = c("#C0AB52",  "#4FBF85"), 
                    na.value = "white") +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  ylab("WA")

naus + wa + se_aus + plot_layout(nrow = 3, guides = "collect")




## Between-species interactions 

between <- full_df2 %>% 
  filter(!(species_from == species_to)) %>% 
  mutate(class_int = if_else(class_int == "large_small", "small_large", class_int))

between <- between %>% 
  #mutate(sig = ifelse((lo > 0 | hi <0), as.character(1), NA)) %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA))

between <- between %>% 
  mutate(species_int = paste(species_from, sep = "_", species_to))

naus <- between %>% 
  filter(region == "N_AUS") %>% 
  ggplot(aes(x = alpha, 
             y = species_int, 
             colour = georegion, fill = fill_col)) + 
  geom_point(shape = 21, size = 3) + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#E78ECD","#28BBD7"), 
                     guide = "none") +
  scale_fill_manual(values = c("#E78ECD","#28BBD7"), 
                    na.value = "white") +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  ylab("N AUS")


se_aus <-  between %>% 
  filter(region %in% c("TAS", "SE_AUS")) %>% 
  group_by(species_from, class_int) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(obs > 2) %>% 
  ggplot(aes(x = alpha, 
             y = species_int, 
             colour = georegion, 
             fill = fill_col
  )) + 
  geom_point(shape = 21, size = 3) + 
  scale_color_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                     guide = "none") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                    na.value = "white") +
  facet_wrap(~class_int) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()  + 
  ylab("SE AUS")



wa <- between %>% 
  filter(region == "WA") %>% 
  ggplot(aes(x = alpha, 
             y = species_int, 
             colour = georegion, 
             fill = fill_col)) + 
  geom_point(shape = 21, size = 3) + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#C0AB52",  "#4FBF85"), guide = "none") + 
  scale_fill_manual(values = c("#C0AB52",  "#4FBF85"), 
                    na.value = "white") +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  ylab("WA")

naus + wa + se_aus + plot_layout(nrow = 3, guides = "collect")



## Ratio 
# Naus 
full_df2 %>% 
  mutate(int = ifelse(species_from == species_to, "within", "between")) %>% 
  filter(region %in% c("N_AUS", "QLD")) %>% 
  group_by(int) %>% 
  count(sig == 1)
#12 sites

#se aus 
full_df2 %>% 
  mutate(int = ifelse(species_from == species_to, "within", "between")) %>% 
  filter(region %in% c("SE_AUS", "TAS")) %>% 
  group_by(int) %>% 
  count(sig == 1)
#27

  
#Assign if species is canopy or understorey 
#load 

species_class <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/data/species_class.csv")
species_class <- species_class %>% dplyr::rename(species_from = species_to)



df <- left_join(df, species_class, by = "species_from")
df <- df %>% rename(cc_from = Class)

df <- df %>% 
  mutate(cc_int = paste0(cc_from, sep = "_", cc_to))



#make visualisation 

between <- df %>% 
  filter(!(species_from == species_to))

between <- between %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) 

between %>% 
  group_by(group) %>% 
  summarise(mean = mean(alpha), 
            median = median(alpha))

ggplot(data = between, 
       aes(x = alpha, 
           y = group)) + 
  geom_boxplot() +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() 



between %>% 
  filter(group == "Subcanopy.small_Subcanopy.small")
  

ggplot(data = between, 
       aes(x = alpha, 
           y = georegion, 
           colour = as.character(sig))) + 
  geom_point(shape = 20, size = 4, alpha = 0.8) +
  geom_vline(xintercept = 0, colour = "red") + 
  facet_wrap(~factor(group, 
                     levels = c(
                       "Subcanopy.small_Subcanopy.small", 
                       "Subcanopy.small_Subcanopy.large", 
                       "Subcanopy.small_Canopy.small", 
                       "Subcanopy.small_Canopy.large", 
                       "Subcanopy.large_Subcanopy.small", 
                       "Subcanopy.large_Subcanopy.large", 
                       "Subcanopy.large_Canopy.small", 
                       "Subcanopy.large_Canopy.large", 
                       "Canopy.small_Subcanopy.small", 
                       "Canopy.small_Subcanopy.large", 
                       "Canopy.small_Canopy.small", 
                       "Canopy.small_Canopy.large", 
                       "Canopy.large_Subcanopy.small", 
                       "Canopy.large_Subcanopy.large", 
                       "Canopy.large_Canopy.small", 
                       "Canopy.large_Canopy.large")), 
             nrow = 4,  
             scales = "free_y") + 
  scale_colour_manual(
    name = "Significance",  
    values = c("darkgoldenrod3", "gray100"),
    labels = c("Sig", "Not Sig")) +
  theme_bw() 




between <- between %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA)) %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
  group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
  group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
  group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
  TRUE ~ group))

ggplot(data = between,
       aes(x = alpha, y = group)) + 
  geom_boxplot()+ 
 geom_point(aes(colour = georegion), alpha = 0.4) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() 
  


within <- within %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion), NA)) %>% 
  mutate(group = paste0(cc_from, ".", class_from,  "_", cc_to, ".", class_to)) %>% 
  mutate(group = case_when(
    group %in% c("Canopy.large_Canopy.small", "Canopy.small_Canopy.large") ~ "Canopy.large_Canopy.small", 
    group %in% c("Canopy.large_Subcanopy.large", "Subcanopy.large_Canopy.large") ~ "Canopy.large_Subcanopy.large", 
    group %in% c("Canopy.large_Subcanopy.small", "Subcanopy.small_Canopy.large") ~  "Canopy.large_Subcanopy.small", 
    group %in% c("Canopy.small_Subcanopy.large", "Subcanopy.large_Canopy.small") ~ "Canopy.small_Subcanopy.large", 
    group %in% c("Canopy.small_Subcanopy.small", "Subcanopy.small_Canopy.small") ~ "Canopy.small_Subcanopy.small", 
    group %in% c("Subcanopy.large_Subcanopy.small", "Subcanopy.small_Subcanopy.large") ~ "Subcanopy.large_Subcanopy.small",
    TRUE ~ group))

ggplot(data = between,
       aes(x = alpha, y = group)) + 
  geom_boxplot()+ 
  geom_point(aes(colour = georegion), alpha = 0.4) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() 



between %>% filter(group %in% c(
  "Canopy.large_Canopy.large", "Canopy.large_Subcanopy.large", 
  "Canopy.large_Subcanopy.small", "Canopy.small_Canopy.small", 
  "Canopy.small_Subcanopy.large", "Subcanopy.large_Subcanopy.large"
)) %>% 
  ggplot(aes(x = alpha, 
             y = from, 
             fill = fill_col,
             colour = georegion)) + 
  geom_point(shape = 21, size = 3) + 
  facet_wrap(~group, nrow = 1) + 
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()



## Plot only those interactions with significance 
between %>% mutate(species_int = paste0(species_from, sep = "_", species_to)) %>% 
  filter(group %in% c("Canopy.large_Subcanopy.large", 
                      "Canopy.large_Subcanopy.small", 
                     "Canopy.small_Subcanopy.large", 
                     "Canopy.small_Subcanopy.small")) %>% 
  filter(species_int %in% 
           c("Eucalyptus grandis_Archontophoenix cunninghamiana",
              "Eucalyptus grandis_Corymbia intermedia",
              "Eucalyptus jacksonii_Allocasuarina decussata",
               "Eucalyptus obliqua_Bedfordia salicina", 
             "Eucalyptus pilularis_Allocasuarina torulosa",
              "Eucalyptus regnans_Acacia melanoxylon",
              "Leptospermum lanigerum_Eucalyptus delegatensis",
              "Leptospermum scoparium_Eucalyptus ovata",
              "Nematolepis squamea_Eucalyptus obliqua",
              "Olearia argophylla_Eucalyptus regnans",
              "Phyllocladus aspleniifolius_Eucalyptus delegatensis",
              "Pomaderris apetala_Eucalyptus obliqua",
              "Pomaderris aspera_Eucalyptus regnans",
               "Schizomeria ovata_Eucalyptus pilularis",
              "Tasmannia lanceolata_Eucalyptus delegatensis",
               "Zieria arborescens_Eucalyptus obliqua")) %>% 
  ggplot(aes(x = alpha, 
        y = species_int, 
        colour = georegion)) + 
  geom_point(shape = 19, size = 3) + 
  facet_wrap(~group, nrow = 1) + 
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()
  


between %>% mutate(species_int = paste0(species_from, sep = "_", species_to)) %>% 
  filter(group %in% c("Canopy.large_Canopy.large", 
                      "Canopy.large_Canopy.small", 
                      "Canopy.small_Canopy.small")) %>% 
  filter(sig == 1) %>% 
  ggplot(aes(x = alpha, 
             y = species_int, 
             colour = georegion)) + 
  geom_point(shape = 19, size = 3) + 
  facet_wrap(~group, nrow = 1) + 
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()




ggplot(data = between, aes(x = alpha, 
           y = group, 
           colour = georegion)) + 
  geom_point(shape = 19, size = 3) + 
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()



#Site covariates 
site_covars <- read.csv("data/site_covariates.csv")

summs <- df %>% 
  mutate(class_int = ifelse(class_int == "large_small", "small_large", class_int)) %>% 
  group_by(site, class_int, cc_to, region) %>% 
  summarise(median = median(alpha))

summs <- summs %>% rename(Site_Name = site)

summs <- merge(x = summs, y = site_covars, by = "Site_Name", all.x = TRUE)

summs$class_int <- factor(summs$class_int,
                           levels = c("small_small", "small_large", "large_large"))


#Latitude 
ggplot(data = summs, 
       aes(y = median, 
           x = AI 
          )) + 
  facet_grid(cc_to~class_int)+
  geom_smooth(method = "lm") +
  geom_point()+ 
  theme_bw()
  



#traits 
library(austraits)
library(APCalign)
traits <- c("plant_height",
            "bark_thickness_index",
            "leaf_area",
            "leaf_thickness")

austraits <- load_austraits(version = "6.0.0", path = "data/austraits")
species <- unique(full_df2$species_from)


#new_names <- create_taxonomic_update_lookup(species) #check names are current 

#check traits for the taxa and traits we are interested in 

austraits <- load_austraits(version = "6.0.0", path = "data/austraits")

t <- austraits %>% 
  extract_trait(trait_names = traits) %>% 
  extract_taxa(taxon_name = species)

df_traits <- as.data.frame(t$traits)

summary_traits1 <- df_traits %>% 
  group_by(taxon_name, trait_name) %>% 
  summarise(minimum = min(as.numeric(value), na.rm = TRUE), 
            median = median(as.numeric(value), na.rm = TRUE), 
            mean = mean(as.numeric(value), na.rm = TRUE), 
            maximum = max(as.numeric(value), na.rm = TRUE)) %>% 
  ungroup()

summary_traits1 <- summary_traits1 %>% 
  mutate(taxon_name = gsub("\\s(subsp.*)$", "", taxon_name))

#add class
species_class <- species_class %>%  rename(taxon_name = Species)
summary_traits1 <- left_join(summary_traits1, species_class, by = "taxon_name")

#for within large
large <- full_df2 %>% 
  filter(class_from == "large") %>% 
  filter(class_to == "large") %>% 
  filter(species_from == species_to) %>% 
  group_by(georegion, species_from) %>% 
  summarise(median = median(alpha))


large <- large %>% rename(taxon_name = species_from)

joined_large <- large %>%
  left_join(summary_traits1, by = "taxon_name") %>% 
  rename(median_alpha = median.x)

joined_large %>% 
  filter(trait_name %in% c("bark_thickness_index", "leaf_area", "leaf_thickness", "plant_height")) %>%  
ggplot(aes(x = median_alpha, 
           y = maximum, 
           shape = Class)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~trait_name, scales = "free_y") + 
  ylab("maximum trait value") + 
  xlab("Within-species large-large alpha coefficient") + theme_bw()

joined_large %>%
  filter(trait_name %in% c("bark_thickness_index", "leaf_area", "leaf_thickness", "plant_height")) %>% 
  ggplot(aes(x = median_alpha, 
             y = median.y, shape = Class, colour = Class)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~trait_name, scales = "free_y") + 
  ylab("median trait value") + 
  xlab("Within-species large-large alpha coefficient") + theme_bw()






small <- full_df2 %>% 
  filter(class_from == "small") %>% 
  filter(class_to == "small") %>% 
  filter(species_from == species_to) %>% 
  group_by(georegion, species_from) %>% 
  summarise(median = median(alpha))


small <- small %>% rename(taxon_name = species_from)


joined_small <- small %>%
  left_join(summary_traits1, by = "taxon_name") %>% 
  rename(median_alpha = median.x)

joined_small %>% 
  filter(! is.na(trait_name)) %>% 
  ggplot(aes(x = median_alpha, 
             y = minimum)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~trait_name, scales = "free_y") + 
  ylab("median trait value") + 
  xlab("Within-species small-small alpha coefficient") + theme_bw()


joined_small %>% 
  filter(trait_name %in% c("bark_thickness_index", "leaf_area", "leaf_thickness", "plant_height")) %>%  
  filter(! is.na(trait_name)) %>% 
  ggplot(aes(x = median_alpha, 
             y = median.y, 
             shape = Class, colour = Class)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~trait_name, scales = "free_y") + 
  ylab("median trait value") + 
  xlab("Within-species small-small alpha coefficient") + theme_bw()



#Correlation? 

#a trend but very small correlation 