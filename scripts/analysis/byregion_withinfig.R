##### Map Figure 
library(ozmaps)
library(ggplot2)
library(dplyr)
library(ggforce)
library(forcats)
library(ggbeeswarm)
library(patchwork)
library(austraits)
library(APCalign)
library(ggpubr)

austraits <- load_austraits(version = "6.0.0", path = "data/austraits")
species <- unique(full_df2$species_from)

#Load data 
data <- read.csv("data/data_cleaned.csv")

#Load results 
df <- df_add3_5_t15

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



######### Another way to do heights 
height <- data_cleaned %>% filter(is.na(Height) == F) ### 24% of trees in ausplots have height measurement

isheight <- data_cleaned %>% 
  group_by(Genus_Species, Site_Name) %>% 
  summarise(na_prop = mean(is.na(Height))) #prop of na for species and site



height <- height %>%
  mutate(georegion = Site_Name) %>%
  mutate(georegion = if_else(georegion %in% c("Dawson", "Frankland", "Clare", "Giants", "Carey", "Dombakup", "Warren",  "Sutton","Collins"),
                             "WA", georegion)) %>%
  mutate(georegion = if_else(georegion %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek", "Weeaproinah", "Turtons", "Lardner", "Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo", "BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx", "BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Dip", "Flowerdale"),
                             "SE_AUS", georegion)) %>%
  mutate(georegion = if_else(georegion %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans", "Baldy", "Koombooloomba", "Lamb Range", "Herberton"),
                             "NE_AUS", georegion))

height <- height %>% 
  group_by(Genus_Species, georegion) %>% 
  summarise(min = min(Height),
            mean = mean(Height), 
            median = median (Height), 
            max = max(Height), 
            q25 = quantile(Height, 0.25),
            q75 = quantile(Height, 0.75))

ggplot(data = height, aes(x = max, y = Genus_Species)) + 
  geom_point() + 
  facet_wrap(~georegion, scales = "free_y")


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
  filter(!species_to %in% c("Pittosporum bicolor", "Hedycarya angustifolia")) %>%
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  mutate(species_from = factor(species_from,
                             levels = c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                                        "Eucalyptus delegatensis", "Eucalyptus cypellocarpa", "Acacia melanoxylon",
                                        "Acacia dealbata", "Nothofagus cunninghamii", "Pomaderris apetala",
                                        "Leptospermum lanigerum",  "Nematolepis squamea",  "Atherosperma moschatum",
                                        "Phyllocladus aspleniifolius","Pomaderris aspera", "Olearia argophylla",  "Monotoca glauca"))) %>%
  mutate(species_from = fct_rev(species_from)) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col)) + 
  geom_vline(xintercept = 0, colour = "red") + 
  geom_hline(yintercept = 10.5, colour = "black", linetype = "dashed") + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 2, orientation = "y") +
  scale_color_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                     guide = "none") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                    na.value = "white", 
                    na.translate = FALSE,
                    labels = c("E.delegatensis TAS", "North VIC", "E.obliqua TAS", "South NSW", "South Vic"),
                    name = "") +
  facet_wrap(~class_int) +
  scale_x_continuous(breaks = c(-1.5, -1, -.5, 0, .5, 1)) +
  theme_bw(base_size = 8)  +
  ylab("") + 
  coord_cartesian(xlim = c(-1.5, 1)) +
  ggtitle("b. Southeastern Australia") + 
  theme(legend.key.height = unit(.5, "cm"))


se_aus

#austraits <- load_austraits(version = "6.0.0", path = "data/austraits")
# traits <- c("plant_height")
# 
# species <- c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
#              "Eucalyptus delegatensis", "Eucalyptus cypellocarpa", "Nothofagus cunninghamii",
#              "Atherosperma moschatum", "Acacia dealbata", "Acacia melanoxylon",
#              "Leptospermum lanigerum", "Phyllocladus aspleniifolius", "Pomaderris aspera",
#              "Pomaderris apetala", "Olearia argophylla", "Nematolepis squamea", "Monotoca glauca")
# t <- austraits %>% 
#   extract_trait(trait_names = traits) %>% 
#   extract_taxa(taxon_name = species)
# 
# df_traits <- as.data.frame(t$traits) %>% 
#   filter(basis_of_record == "literature") %>% 
#    filter(value_type == "maximum") %>% 
#   group_by(taxon_name) %>% 
#   summarise(median = median(as.numeric(value)), 
#             q3 = quantile(as.numeric(value), 0.75), 
#             max = max(as.numeric(value))) %>% 
#   ungroup()
# 
# 
# se_height <- df_traits %>% filter(taxon_name %in% c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
#                                                     "Eucalyptus delegatensis", "Eucalyptus cypellocarpa", "Nothofagus cunninghamii",
#                                                     "Atherosperma moschatum", "Acacia dealbata", "Acacia melanoxylon",
#                                                     "Leptospermum lanigerum", "Phyllocladus aspleniifolius", "Pomaderris aspera",
#                                                     "Pomaderris apetala", "Olearia argophylla", "Nematolepis squamea", "Monotoca glauca")) %>% 
#   mutate(taxon_name = factor(taxon_name,
#                              levels = c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
#                                         "Eucalyptus delegatensis", "Eucalyptus cypellocarpa", "Nothofagus cunninghamii",
#                                         "Atherosperma moschatum", "Acacia dealbata", "Acacia melanoxylon",  "Phyllocladus aspleniifolius",
#                                         "Leptospermum lanigerum", "Pomaderris aspera",
#                                         "Pomaderris apetala", "Olearia argophylla", "Nematolepis squamea", "Monotoca glauca"))) %>% 
#   mutate(taxon_name = fct_rev(taxon_name)) %>% 
#   ggplot(aes(y = taxon_name)) + 
#   geom_hline(yintercept = 11.5, colour = "black", linetype = "dashed") + 
#   geom_boxplot(aes(xmin = median, xlower = median, xmiddle = median, xupper = q3, xmax = max), width = 0.6, stat = "identity") +
#   theme_bw(base_size = 8) + 
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank()) 


se_height <- height %>% 
  filter(georegion == "SE_AUS") %>% 
  filter(Genus_Species %in% c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                              "Eucalyptus delegatensis", "Eucalyptus cypellocarpa", "Acacia melanoxylon",
                              "Acacia dealbata", "Nothofagus cunninghamii",  "Pomaderris apetala",
                              "Leptospermum lanigerum",  "Nematolepis squamea",  "Atherosperma moschatum",
                              "Phyllocladus aspleniifolius","Pomaderris aspera", "Olearia argophylla",  "Monotoca glauca")) %>% 
  mutate(taxon_name = factor(Genus_Species,
                             levels = c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                                        "Eucalyptus delegatensis", "Eucalyptus cypellocarpa", "Acacia melanoxylon",
                                        "Acacia dealbata", "Nothofagus cunninghamii", "Pomaderris apetala",
                                        "Leptospermum lanigerum",  "Nematolepis squamea",  "Atherosperma moschatum",
                                        "Phyllocladus aspleniifolius","Pomaderris aspera", "Olearia argophylla",  "Monotoca glauca"))) %>%
                             mutate(taxon_name = fct_rev(taxon_name)) %>%
  ggplot(aes(y = taxon_name)) +
  geom_hline(yintercept = 11.5, colour = "black", linetype = "dashed") + 
  geom_boxplot(aes(xmin = q75, xlower = q75, xmiddle = q75, xupper = q75, xmax = max), 
               width = 0.6, 
               stat = "identity",
               position = position_dodge(width = .75)) +
  theme_bw(base_size = 8) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

se_height

se <- se_aus + se_height + 
  plot_layout(nrow = 1, widths = c(3, 1),
              guides = "collect")

se

ggsave("se_within.png", se, 
       height = 3.5, width = 8, dpi = 300)

### WA 
wa_aus <- within %>% 
  filter(region == "WA") %>% 
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  mutate(species_from = factor(species_from, levels = c("Eucalyptus diversicolor",
                                                    "Eucalyptus jacksonii","Eucalyptus guilfoylei",  "Corymbia calophylla",
                                                    "Allocasuarina decussata","Acacia melanoxylon", 
                                                    "Trymalium odoratissimum"))) %>%
  mutate(species_from = fct_rev(species_from)) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col)) +
  geom_vline(xintercept = 0, colour = "red") + 
  geom_hline(yintercept = 3.5, colour = "black", linetype = "dashed") + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 2, orientation = "y") + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#C0AB52",  "#4FBF85"), guide = "none") + 
  scale_fill_manual(values = c("#C0AB52",  "#4FBF85"), 
                    na.value = "white",
                    labels = c("North WA", "South WA"),
                    na.translate = FALSE,
                    name = "") +
  theme_bw(base_size = 8) + 
  ylab("") + 
  xlab("Interaction coefficient") + 
  ggtitle("c. Western Australia") + 
  theme(legend.key.height = unit(.5, "cm"))
wa_aus

wa_height <- height %>% 
  filter(georegion == "WA") %>% 
  filter(Genus_Species %in% c("Eucalyptus diversicolor",
                              "Eucalyptus jacksonii","Eucalyptus guilfoylei",  "Corymbia calophylla",
                              "Allocasuarina decussata", "Acacia melanoxylon",
                              "Trymalium odoratissimum")) %>% 
  mutate(taxon_name = factor(Genus_Species,
                             levels = c("Eucalyptus diversicolor",
                                        "Eucalyptus jacksonii","Eucalyptus guilfoylei",  "Corymbia calophylla",
                                        "Allocasuarina decussata","Acacia melanoxylon", 
                                        "Trymalium odoratissimum"))) %>%
  mutate(taxon_name = fct_rev(taxon_name)) %>%
  ggplot(aes(y = taxon_name)) +
  geom_hline(yintercept = 3.5, colour = "black", linetype = "dashed") + 
  geom_boxplot(aes(xmin = q75, xlower = q75, xmiddle = q75, xupper = q75, xmax = max), 
               width = 0.6, 
               stat = "identity") +
  theme_bw(base_size = 8) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + 
  xlab("Maximum height (m)")
wa_height

#austraits <- load_austraits(version = "6.0.0", path = "data/austraits")
# traits <- c("plant_height")
# 
# species <- c("Acacia melanoxylon", "Allocasuarina decussata", 
#              "Corymbia calophylla", "Eucalyptus diversicolor", 
#              "Eucalyptus guilfoylei", "Eucalyptus jacksonii", 
#              "Trymalium odoratissimum")
# t <- austraits %>% 
#   extract_trait(trait_names = traits) %>% 
#   extract_taxa(taxon_name = species)
# 
# df_traits <- as.data.frame(t$traits) %>% 
#  # filter(!basis_of_record == "literature") %>% 
#   filter(value_type == "maximum") %>% 
#   group_by(taxon_name) %>% 
#   summarise(median = median(as.numeric(value)), 
#             q3 = quantile(as.numeric(value), 0.75), 
#             max = max(as.numeric(value))) %>% 
#   ungroup()
# 
# 

# 
# wa_height <- df_traits %>% filter(taxon_name %in% c("Acacia melanoxylon", "Allocasuarina decussata", 
#                                                     "Corymbia calophylla", "Eucalyptus diversicolor", 
#                                                     "Eucalyptus guilfoylei", "Eucalyptus jacksonii", 
#                                                     "Trymalium odoratissimum")) %>% 
#   mutate(taxon_name = factor(taxon_name, levels = c("Eucalyptus diversicolor", 
#                                                     "Eucalyptus jacksonii","Eucalyptus guilfoylei",  "Corymbia calophylla", 
#                                                     "Acacia melanoxylon", "Allocasuarina decussata", 
#                                                     "Trymalium odoratissimum"))) %>% 
#   mutate(taxon_name = fct_rev(taxon_name)) %>%  
#   ggplot(aes(y = taxon_name)) + 
#   geom_hline(yintercept = 3.5, colour = "black", linetype = "dashed") + 
#   geom_boxplot(aes(xmin = median, xlower = median, xmiddle = median, xupper = q3, xmax = max), width = 0.6, stat = "identity") +
#   theme_bw(base_size = 8) +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank()) + 
#   xlab("Maximum height (m)")


wa <- wa_aus + wa_height+ 
  plot_layout(nrow = 1, widths = c(3, 1), guides = "collect")
wa

ggsave("wa.png", wa, 
       height = 5, width = 8, dpi = 300)

### NE AUS
ne_aus <- within %>% 
  filter(region == "N_AUS") %>% 
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  filter(! species_from %in% c("Lophostemon sp", "Geissois benthamii", "Cryptocarya microneura", "Acmena smithii", "Synoum glandulosum")) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  mutate(species_from = ifelse(species_from == "Caldcluvia paniculosa", "Ackama paniculosa", species_from)) %>% 
  mutate(species_from = factor(species_from, levels = c( "Eucalyptus pilularis","Eucalyptus grandis",
                                                        "Eucalyptus microcorys","Syncarpia glomulifera", "Corymbia intermedia", "Allocasuarina torulosa", 
                                                        "Schizomeria ovata","Ceratopetalum apetalum",  "Acacia melanoxylon", "Archontophoenix cunninghamiana", "Ackama paniculosa", 
                                                        "Cryptocarya glaucescens", "Cryptocarya rigida", "Polyscias elegans", "Sloanea langii"))) %>%
  mutate(species_from = fct_rev(species_from)) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, fill = fill_col)) + 
  geom_vline(xintercept = 0, colour = "red") + 
  geom_hline(yintercept = 10.5, colour = "black", linetype = "dashed") + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 2, orientation = "y") +
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#E78ECD","#28BBD7"), 
                     guide = "none") +
  scale_fill_manual(values = c("#E78ECD","#28BBD7"), 
                    na.value = "white", 
                    labels = c("North NSW", "QLD"),
                    name = "Region", 
                    na.translate = FALSE) +
  theme_bw(base_size = 8) + 
  coord_cartesian(xlim = c(-1.5, 1)) +
  ylab("") + 
  ggtitle("a. Northeastern Australia") + 
  theme(legend.key.height = unit(.5, "cm"))

ne_aus

ne_height <- height %>% 
  filter(georegion == "NE_AUS") %>% 
  filter(Genus_Species %in% c( "Eucalyptus pilularis","Eucalyptus grandis",
                              "Eucalyptus microcorys","Syncarpia glomulifera", "Corymbia intermedia", "Allocasuarina torulosa", 
                              "Schizomeria ovata","Ceratopetalum apetalum",  "Acacia melanoxylon", "Archontophoenix cunninghamiana", "Ackama paniculosa", 
                              "Cryptocarya glaucescens", "Cryptocarya rigida", "Polyscias elegans", "Sloanea langii")) %>% 
  mutate(taxon_name = factor(Genus_Species,
                             levels = c("Eucalyptus pilularis", "Eucalyptus grandis", 
                                        "Eucalyptus microcorys","Syncarpia glomulifera", "Corymbia intermedia", "Allocasuarina torulosa", 
                                        "Schizomeria ovata","Ceratopetalum apetalum",  "Acacia melanoxylon", "Archontophoenix cunninghamiana", "Ackama paniculosa", 
                                        "Cryptocarya glaucescens", "Cryptocarya rigida", "Polyscias elegans", "Sloanea langii"))) %>%
  mutate(taxon_name = fct_rev(taxon_name)) %>%
  ggplot(aes(y = taxon_name)) +
  geom_hline(yintercept = 9.5, colour = "black", linetype = "dashed") + 
  geom_boxplot(aes(xmin = q75, xlower = q75, xmiddle = q75, xupper = q75, xmax = max), 
               width = 0.6, 
               stat = "identity") +
  theme_bw(base_size = 8) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 
ne_height

#austraits <- load_austraits(version = "6.0.0", path = "data/austraits")
# traits <- c("plant_height")
# species <- c("Acacia melanoxylon", "Allocasuarina torulosa", 
#              "Corymbia intermedia", "Eucalyptus grandis", 
#              "Eucalyptus pilularis", "Eucalyptus microcorys", 
#              "Archontophoenix cunninghamiana", "Ackama paniculosa",
#              "Ceratopetalum apetalum", "Cryptocarya glaucescens",
#              "Cryptocarya rigida", "Polyscias elegans",
#              "Schizomeria ovata", "Sloanea langii", "Syncarpia glomulifera")
# t <- austraits %>% 
#   extract_trait(trait_names = traits) %>% 
#   extract_taxa(taxon_name = species)
# 
# df_traits <- as.data.frame(t$traits) %>% 
#    filter(!basis_of_record == "literature") %>% 
#   filter(value_type == "maximum") %>% 
#   group_by(taxon_name) %>% 
#   summarise(median = median(as.numeric(value)), 
#             q3 = quantile(as.numeric(value), 0.75), 
#             max = max(as.numeric(value))) %>% 
#   ungroup()
# 
# ne_height <- df_traits %>% filter(taxon_name %in% c("Acacia melanoxylon", "Allocasuarina torulosa", 
#                                                     "Corymbia intermedia", "Eucalyptus grandis", 
#                                                     "Eucalyptus pilularis", "Eucalyptus microcorys", 
#                                                  "Archontophoenix cunninghamiana", "Ackama paniculosa",
#                                                  "Ceratopetalum apetalum", "Cryptocarya glaucescens",
#                                                  "Cryptocarya rigida", "Polyscias elegans",
#                                                  "Schizomeria ovata", "Sloanea langii", "Syncarpia glomulifera")) %>% 
#  mutate(taxon_name = factor(taxon_name, levels = c("Eucalyptus grandis", "Eucalyptus pilularis",
#                                                    "Eucalyptus microcorys","Ackama paniculosa", "Corymbia intermedia", "Cryptocarya glaucescens",
#                                                    "Schizomeria ovata","Ceratopetalum apetalum",   "Archontophoenix cunninghamiana",
#                                                    "Acacia melanoxylon", "Polyscias elegans", "Sloanea langii", "Syncarpia glomulifera",
#                                                    "Allocasuarina torulosa", "Cryptocarya rigida"))) %>% 
#   mutate(taxon_name = fct_rev(taxon_name)) %>%  
#   ggplot(aes(y = taxon_name)) + 
#   geom_hline(yintercept = 11.5, colour = "black", linetype = "dashed") + 
#   geom_boxplot(aes(xmin = median, xlower = median, xmiddle = median, xupper = q3, xmax = max), width = 0.6, stat = "identity") +
#   theme_bw(base_size = 8) + 
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank())
# 


ne <- ne_aus + ne_height + 
  plot_layout(nrow = 1, ncol = 2, widths = c(3, 1), guides = "collect")
ne
ggsave("ne_within.png", ne, 
       height = 5, width = 8, dpi = 300)

#### Stick all together 

f <- ne/se/wa + plot_layout(guides = "collect")

ggsave("fig4.png", f, 
       height = 11, width = 8, dpi = 300)

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
  




#Tabulate species 
dd <- data %>%  mutate(georegion = case_when(Site_Name %in% c("Weeaproinah", "Turtons", "Lardner") ~ "S_VIC", 
                                             Site_Name %in% c("ANU101", "ANU363", "ANU589", "Ada Tree", "HardyCreek") ~ "N_VIC", 
                                             Site_Name %in% c("Dawson", "Frankland", "Clare", "Giants") ~ "S_WA",
                                             Site_Name %in% c("Carey", "Dombakup", "Warren",  "Sutton","Collins") ~ "N_WA", 
                                             Site_Name %in% c("Baldy", "Koombooloomba", "Lamb Range", "Herberton") ~ "QLD", 
                                             Site_Name %in% c("MinesRd", "A-Tree", "BirdTree", "BlackBull", "Lorne", "Tinebank", "Bruxner", "Osullivans") ~ "N_NSW", 
                                             Site_Name %in% c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo") ~ "S_NSW", 
                                             Site_Name %in% c("BenRidge", "Caveside", "Mackenzie", "MtField", "MtMaurice", "NorthStyx") ~ "d_TAS", 
                                             Site_Name %in% c("BondTier", "BlackRiver", "Weld", "MtField", "ZigZag", "Supersite", "Bird", "Flowerdale", "Dip") ~ "o_TAS")) %>% 
  mutate(region = case_when(georegion %in% c("S_VIC", "N_VIC", "S_NSW") ~ "SE_AUS", 
                            georegion %in% c("S_WA", "N_WA") ~ "WA",
                            georegion %in% c("QLD", "N_NSW") ~ "N_AUS",  
                            georegion %in% c("o_TAS", "d_TAS") ~ "TAS"))


dd <- dd %>% 
  group_by(region) %>% 
  count(Genus_Species)






############ Summary figure 

within <- df_add3_5_t15 %>% 
  filter(species_from == species_to) 

within <- within %>% 
  mutate(class_int = ifelse(class_int == "large_small", "small_large", class_int))

summary <- within %>% 
  group_by(class_int, region, species_from) %>% 
  summarise(median = median(alpha),
            mean = mean(alpha), 
            max = max(alpha), 
            min = min(alpha))

  

summary %>% filter(region %in% c("SE_AUS", "TAS")) %>% 
ggplot() + 
  geom_point(aes(x = median, y = species_from, colour = region), size = 2) + 
  facet_grid(~class_int, scales = "free_y") + 
  geom_vline(xintercept = 0, colour = "red", linetype = "dotted") + 
  geom_vline(xintercept = -0.1, colour = "red", linetype = "dotted") + 
  geom_vline(xintercept = 0.1, colour = "red", linetype = "dotted") + 
  theme_bw() + 
  xlim(-.75, 1)


summary %>% filter(species_from == "Eucalyptus delegatensis")


within <- within %>% 
  mutate(sig = ifelse(is.na(sig), 0, sig))


sigs <- within %>% 
  count(region, sig)
  
  
  