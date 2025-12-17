##### Map Figure 
library(ozmaps)
library(ggplot2)
library(dplyr)
library(ggforce)
library(forcats)
library(ggbeeswarm)
library(patchwork)

#Load data 
data <- read.csv("data/data_cleaned.csv")

basemap <- ggplot() +
  geom_sf(data = ozmap(), fill = "gray100") +
  geom_point(data = data, aes(x = Longitude, y = Latitude, colour = Site_Name),
             size = 2.5) + 
  xlim(c(110, 156)) + 
  theme_classic() + 
  theme(legend.position = "none")



basemap2 <- basemap + 
  geom_circle(data = circle, 
              aes(x0 = Longitude, y0 = Latitude, r = radius), 
              color = "black",
               inherit.aes = FALSE) +
  geom_text(data = circle,
            aes(x = Longitude, y = Latitude, label = Site_Name),
            position = position_nudge(x = 2.5),
            size = 4,
            inherit.aes = FALSE)

  
ggsave("basemap.png", basemap2, 
       width = 12, height = 8, 
       dpi = 300)

  # + geom_segment(data = lines, 
  #              aes(x = x, y = y, 
  #                  xend = xend, yend = yend), 
  #              arrow = arrow(length = unit(0.2, "cm")), 
  #              type = "closed",
  #              color = "black",
  #              inherit.aes = FALSE)


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




se_aus <- within %>% 
  filter(region %in% c("TAS", "SE_AUS")) %>%  
  group_by(species_from, class_int) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(obs > 1) %>% 
  #filter(!species_to == "Pittosporum bicolor") %>%
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
 # mutate(species_from = factor(species_from,
                             # levels = c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                             #            "Eucalyptus delegatensis", "Nothofagus cunninghamii", 
                             #            "Atherosperma moschatum", "Acacia dealbata", "Acacia melanoxylon",
                             #            "Leptospermum lanigerum", "Phyllocladus aspleniifolius", 
                             #            "Pomaderris apetala", "Olearia argophylla", "Nematolepis squamea"))) %>% 
  mutate(species_from = fct_rev(species_from)) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col)) + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 3, orientation = "y") +
  scale_color_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                     guide = "none") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                    na.value = "white", 
                    na.translate = FALSE,
                    name = "") +
  facet_wrap(~class_int) +
  scale_x_continuous(breaks = c(-1.5, -1, -.5, 0, .5, 1)) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()  +
  ylab("") + 
  coord_cartesian(xlim = c(-1.5, 1)) +
  ggtitle("b. Southeastern Australia")


se_aus



se_height <- df_traits %>% filter(taxon_name %in% c("Acacia dealbata", "Acacia melanoxylon", 
                                       "Atherosperma moschatum", "Eucalyptus delegatensis", 
                                       "Eucalyptus fastigata", "Eucalyptus obliqua", 
                                       "Eucalyptus regnans", "Leptospermum lanigerum", 
                                       "Nematolepis squamea", "Nothofagus cunninghamii",
                                       "Olearia argophylla", "Phyllocladus aspleniifolius",
                                       "Pomaderris apetala")) %>% 
  filter(!observation_id %in% c("02287", "1414", "00241", "00784")) %>% #weird outliers
  mutate(taxon_name = factor(taxon_name,
                             levels = c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                                        "Eucalyptus delegatensis", "Nothofagus cunninghamii", 
                                        "Atherosperma moschatum", "Acacia dealbata", "Acacia melanoxylon",
                                        "Leptospermum lanigerum", "Phyllocladus aspleniifolius", 
                                        "Pomaderris apetala", "Olearia argophylla", "Nematolepis squamea"))) %>% 
  mutate(taxon_name = fct_rev(taxon_name)) %>% 
  ggplot(aes(x = as.numeric(value), y = taxon_name)) + 
  geom_boxplot() +
  theme_bw()
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


se <- se_aus + se_height + 
  plot_layout(widths = c(2, 1))

se

ggsave("fig1.png", se, 
       height = 7, width = 12, dpi = 300)


wa_aus <- within %>% 
  filter(region == "WA") %>% 
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  # mutate(species_from = factor(species_from, levels = c("Eucalyptus diversicolor", 
  #                                                   "Eucalyptus jacksonii","Eucalyptus guilfoylei",  "Corymbia calophylla", 
  #                                                   "Acacia melanoxylon", "Allocasuarina decussata", 
  #                                                   "Trymalium odoratissimum"))) %>% 
  mutate(species_from = fct_rev(species_from)) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col)) + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 3, orientation = "y") + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#C0AB52",  "#4FBF85"), guide = "none") + 
  scale_fill_manual(values = c("#C0AB52",  "#4FBF85"), 
                    na.value = "white",
                    na.translate = FALSE,
                    name = "") +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  ylab("") + 
  xlab("Interaction coefficient") + 
  ggtitle("c. Western Australia")
wa_aus

wa_height <- df_traits %>% filter(taxon_name %in% c("Acacia melanoxylon", "Allocasuarina decussata", 
                                                    "Corymbia calophylla", "Eucalyptus diversicolor", 
                                                    "Eucalyptus guilfoylei", "Eucalyptus jacksonii", 
                                                    "Trymalium odoratissimum")) %>% 
  filter(!observation_id %in% c("00784", "05117", "05047", "05132")) %>% #weird outliers
  mutate(taxon_name = factor(taxon_name, levels = c("Eucalyptus diversicolor", 
                                                    "Eucalyptus jacksonii","Eucalyptus guilfoylei",  "Corymbia calophylla", 
                                                    "Acacia melanoxylon", "Allocasuarina decussata", 
                                                    "Trymalium odoratissimum"))) %>% 
  mutate(taxon_name = fct_rev(taxon_name)) %>%  
  ggplot(aes(x = as.numeric(value), y = taxon_name)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + 
  xlab("Maximum height (cm)")

wa_height 

wa <- wa_aus + wa_height+ 
  plot_layout(widths = c(2, 1))


ggsave("wa.png", wa, 
       height = 7, width = 12, dpi = 300)



ne_aus <- within %>% 
  filter(region == "N_AUS") %>% 
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  #filter(! species_from %in% c("Niemeyeria whitei", "Lophostemon sp", "Geissois benthamii", "Cryptocarya microneura", "Acmena smithii")) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
#  mutate(species_from = ifelse(species_from == "Caldcluvia paniculosa", "Ackama paniculosa", species_from)) %>% 
  # mutate(species_from = factor(species_from, levels = c("Eucalyptus grandis", "Eucalyptus pilularis",
  #                                                   "Eucalyptus microcorys","Ackama paniculosa", "Corymbia intermedia", "Cryptocarya glaucescens",
  #                                                   "Schizomeria ovata","Ceratopetalum apetalum",   "Archontophoenix cunninghamiana", 
  #                                                   "Acacia melanoxylon", "Polyscias elegans", "Sloanea langii", "Syncarpia glomulifera",
  #                                                   "Allocasuarina torulosa", "Cryptocarya rigida", "Synoum glandulosum"))) %>% 
  mutate(species_from = fct_rev(species_from)) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, fill = fill_col)) + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 3, orientation = "y") +
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#E78ECD","#28BBD7"), 
                     guide = "none") +
  scale_fill_manual(values = c("#E78ECD","#28BBD7"), 
                    na.value = "white", 
                    name = "Region", 
                    na.translate = FALSE) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  coord_cartesian(xlim = c(-1.5, 1)) +
  ylab("") + 
  ggtitle("a. Northeastern Australia")

ne_aus

ne_height <- df_traits %>% filter(taxon_name %in% c("Acacia melanoxylon", "Allocasuarina torulosa", 
                                                    "Corymbia intermedia", "Eucalyptus grandis", 
                                                    "Eucalyptus pilularis", "Eucalyptus microcorys", 
                                                 "Archontophoenix cunninghamiana", "Ackama paniculosa",
                                                 "Ceratopetalum apetalum", "Cryptocarya glaucescens",
                                                 "Cryptocarya rigida", "Polyscias elegans",
                                                 "Schizomeria ovata", "Sloanea langii", "Syncarpia glomulifera",
                                                 "Synoum glandulosum")) %>%   
  filter(!observation_id %in% c("05215", "476")) %>% #weird outliers
 mutate(taxon_name = factor(taxon_name, levels = c("Eucalyptus grandis", "Eucalyptus pilularis",
                            "Eucalyptus microcorys","Ackama paniculosa", "Corymbia intermedia", "Cryptocarya glaucescens",
                            "Schizomeria ovata","Ceratopetalum apetalum",   "Archontophoenix cunninghamiana", 
                            "Acacia melanoxylon", "Polyscias elegans", "Sloanea langii", "Syncarpia glomulifera",
                            "Allocasuarina torulosa", "Cryptocarya rigida", "Synoum glandulosum"))) %>% 
  mutate(taxon_name = fct_rev(taxon_name)) %>%  
  ggplot(aes(x = as.numeric(value), y = taxon_name)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


ne_height 

ne <- ne_aus + ne_height + 
  plot_layout(widths = c(2, 1))

ggsave("ne.png", ne, 
       height = 7, width = 12, dpi = 300)


#### Stick all together 

f <- ne + se + wa + plot_layout(nrow =3, ncol = 1)
f

ggsave("fig1.png", ne, 
       height = 7, width = 12, dpi = 300)

fig <- basemap2 + f + plot_layout(ncol = 1, heights = c(0.62, 2))

ggsave("mapfig2.png", fig, 
       height = 18.5, width = 12, dpi = 300)


#### Get some max heights 
#traits 
library(austraits)
library(APCalign)
austraits <- load_austraits(version = "6.0.0", path = "data/austraits")

#Let's check names first 
species <- unique(df$species_from)

new_names <- create_taxonomic_update_lookup(species)
new_names$accepted_name

trait <- "plant_height"

t <- austraits %>% 
  extract_trait(trait_names = trait) %>% 
  extract_taxa(taxon_name = new_names$accepted_name)

df_height <- as.data.frame(t$traits)

df_traits  <- df_height %>% 
  filter(value_type == "maximum") %>%  #only want max height 
  filter(!taxon_name %in%  c("Acacia dealbata subsp. subalpina", "Acacia mucronata subsp. dependens", "Atherosperma moschatum subsp. integrifolium", 
                             "Eucalyptus radiata subsp. robertsonii")) 
# %>% 
#   mutate(value = as.numeric(value)) %>% 
#   group_by(taxon_name, dataset_id, basis_of_record) %>% 
#   summarise(mean = mean(value), 
#             median = median(value), 
#             max = max(value), 
#             min = min(value))
  




# # to use ggmap need a google API key 
# register_google(key = "AIzaSyAI6lTAvGDB_ejEBsbCFlmmCMwjXY28jG4")
# now requires to put in billing details to use 
# map <- get_googlemap("Australia", zoom = 4, maptype = "terrain")
# 
# # Plot it
# ggmap(map) + 
#   theme_void() + 
#   ggtitle("terrain") + 
#   theme(
#     plot.title = element_text(colour = "orange"), 
#     panel.border = element_rect(colour = "grey", fill=NA, size=2)
#   )