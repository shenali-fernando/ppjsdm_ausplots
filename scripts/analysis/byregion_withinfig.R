##### Fig 3: Species + Size Model Only By Region and Species

library(ggplot2)
library(dplyr)
library(ggforce)
library(forcats)
library(ggbeeswarm)
library(patchwork)
library(ggpubr)

# Read in model specification dfs 
#or updated sp_size model 
sp_size <- read.csv("scripts/model_specifications/species_diameter/sp_size_df_t15_misc_updated.csv")

#austraits <- load_austraits(version = "6.0.0", path = "data/austraits")

#Load data 
data <- read.csv("data/ausplots/data_cleaned.csv")


### Thinking about this figure, we're trying to get the maximum height of the trees (as a proxy for size) and the alpha 


######### HEIGHT 
height <- data%>% filter(is.na(Height) == F) ### 24% of trees in ausplots have height measurement

isheight <- data %>% 
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

#and plot
ggplot(data = height, aes(x = max, y = Genus_Species)) + 
  geom_point() + 
  facet_wrap(~georegion, scales = "free_y")

################################################
## Within-species interactions 
within <- sp_size %>% 
  filter(species_from == species_to)

within <- within |> 
  filter(! from %in% c("Misc_large", "Misc_small")) %>% 
  filter(! to %in% c("Misc_large", "Misc_small"))

#need to adjust georegion to align with map fig 1 
within <- within |> 
  mutate(georegion2 = case_when(georegion %in% c("N_VIC", "S_VIC") ~ "Victoria", 
                                georegion == "S_NSW" ~ "Southeastern NSW",
                               georegion %in% c("o_TAS", "d_TAS") ~ "Tasmania",
                               georegion %in% c("S_WA", "N_WA") ~ "Southwestern WA",
                               georegion == "N_NSW" ~ "Northeastern NSW", 
                               georegion == "QLD" ~  "Northeastern QLD",
                              TRUE ~ georegion)) |> 
  mutate(region = case_when(georegion %in% c("o_TAS", "d_TAS", "N_VIC", "S_VIC", "S_NSW") ~ "Southeastern Australia", 
                            georegion %in% c("N_NSW", "QLD") ~ "Northeastern Australia", 
                            georegion %in% c("N_WA", "S_WA") ~ "Southwestern Australia", 
                            TRUE ~ georegion))

#create column for plotting 
within <- within %>% 
  mutate(fill_col = ifelse(sig == 1, as.character(georegion2), NA))

#fix names and levelling of size int 
within <- within |> 
  mutate(size_int = case_when(
    size_int %in% c("small small", "small_small") ~ "Small ↔ Small", 
    size_int %in% c("small large", "large_small", "small_large") ~ "Small ↔ Large", 
    size_int %in% c("large large", "large_large") ~ "Large ↔ Large", 
    TRUE ~ size_int)) |> 
  mutate(size_int = as.factor(size_int)) |>  
  mutate(size_int = fct_relevel(size_int, 
                                "Large ↔ Large", 
                                "Small ↔ Large",   
                                "Small ↔ Small"))

## remove those that dont have at least one obs in each of the size ints

############ For the figure, we need to do each panel separately as facetting won't work: 
#### SE AUS


se_aus <- within |> 
  filter(region == "Southeastern Australia") |> 
  mutate(size_int = fct_relevel(size_int, 
                                "Small ↔ Small", "Small ↔ Large",  "Large ↔ Large", )) |> 
  filter(species_from %in% c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                              "Eucalyptus delegatensis", "Eucalyptus viminalis", "Eucalyptus radiata", "Eucalyptus cypellocarpa",
                              "Acacia melanoxylon", "Acacia dealbata", "Nothofagus cunninghamii", "Pomaderris apetala",
                              "Leptospermum lanigerum",  "Nematolepis squamea",  "Atherosperma moschatum",
                              "Phyllocladus aspleniifolius", "Olearia argophylla")) |>
  mutate(species_from = factor(species_from,
                             levels = c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                              "Eucalyptus delegatensis","Eucalyptus viminalis", "Eucalyptus radiata",  "Eucalyptus cypellocarpa",
                              "Acacia melanoxylon", "Acacia dealbata", "Nothofagus cunninghamii", "Pomaderris apetala",
                              "Leptospermum lanigerum",  "Nematolepis squamea",  "Atherosperma moschatum",
                              "Phyllocladus aspleniifolius", "Olearia argophylla"))) |>
  mutate(species_from = fct_rev(species_from)) |>
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion2, 
             fill = fill_col)) + 
  geom_vline(xintercept = 0, 
             colour = "black", 
             linewidth = 0.2,
             linetype = "dotted") + 
   geom_hline(yintercept = 9.5, 
              colour = "black", 
              linetype = "dashed") + 
  geom_quasirandom(shape = 21, 
                   stroke = 1, 
                   width = 0.25, 
                   alpha = 0.7,
                   size = 2, 
                   orientation = "y") +
  scale_color_manual(values = c("#18DAC7FF", "#FEAA33FF", "#E03F08FF"), 
                     guide = "none") + 
  scale_fill_manual(values = c("#18DAC7FF","#FEAA33FF", "#E03F08FF"), 
                    na.value = "white", 
                    na.translate = FALSE,
                    name = "") +
  facet_wrap(~size_int) +
  scale_x_continuous(breaks = c(-1.5, -1, -.5, 0, .5, 1)) +
  theme_bw()  +
  ylab("") + 
  xlab("") + 
  coord_cartesian(xlim = c(-1.5, 1)) +
  ggtitle("b. Southeastern Australia") 



se_aus


se_height <- height %>% 
  filter(georegion == "SE_AUS") %>% 
  filter(Genus_Species %in%  c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                               "Eucalyptus delegatensis","Eucalyptus viminalis", "Eucalyptus radiata", "Eucalyptus cypellocarpa",
                               "Acacia melanoxylon", "Acacia dealbata", "Nothofagus cunninghamii", "Pomaderris apetala",
                               "Leptospermum lanigerum",  "Nematolepis squamea",  "Atherosperma moschatum",
                               "Phyllocladus aspleniifolius", "Olearia argophylla")) %>% 
  mutate(taxon_name = factor(Genus_Species,
                             levels =  c("Eucalyptus regnans", "Eucalyptus obliqua", "Eucalyptus fastigata",
                                         "Eucalyptus delegatensis","Eucalyptus viminalis", "Eucalyptus radiata", "Eucalyptus cypellocarpa",
                                         "Acacia melanoxylon", "Acacia dealbata", "Nothofagus cunninghamii", "Pomaderris apetala",
                                         "Leptospermum lanigerum",  "Nematolepis squamea",  "Atherosperma moschatum",
                                         "Phyllocladus aspleniifolius", "Olearia argophylla"))) %>%
                             mutate(taxon_name = fct_rev(taxon_name)) %>%
  ggplot(aes(y = taxon_name)) +
  geom_hline(yintercept = 9.5, colour = "black", linetype = "dotdash") + 
  geom_boxplot(aes(xmin = q75, xlower = q75, xmiddle = q75, xupper = q75, xmax = max), 
               width = 0.6, 
               stat = "identity",
               position = position_dodge(width = .75)) +
  coord_cartesian(xlim = c(13, 80)) +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

se_height




### WA 
wa_aus <- within |> 
  filter(region == "Southwestern Australia") |> 
  mutate(size_int = fct_relevel(size_int, 
                                "Small ↔ Small", "Small ↔ Large",  "Large ↔ Large", )) |> 
  mutate(species_from = factor(species_from, levels = c("Eucalyptus diversicolor",
                                                    "Eucalyptus jacksonii","Eucalyptus guilfoylei",  "Corymbia calophylla",
                                                    "Allocasuarina decussata","Acacia melanoxylon",
                                                    "Trymalium odoratissimum"))) |>
  mutate(species_from = fct_rev(species_from)) |>
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion2, 
             fill = fill_col)) +
  geom_vline(xintercept = 0, 
             colour = "black", 
             linewidth = 0.2,
             linetype = "dotted") + 
  geom_hline(yintercept = 3.5, colour = "black", linetype = "dashed") + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 2, orientation = "y") + 
  facet_wrap(~size_int) +
  scale_color_manual(values = c("#7A0403FF"), 
                     guide = "none") + 
  scale_fill_manual(values = c("#7A0403FF"), 
                    na.value = "white",
                    na.translate = FALSE,
                    name = "") +
  theme_bw() + 
  coord_cartesian(xlim = c(-1.5, 1)) +
  scale_x_continuous(breaks = c(-1.5, -1, -.5, 0, .5, 1)) +
  ylab("") + 
  xlab("Interaction coefficient") + 
  ggtitle("c. Western Australia")
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
  coord_cartesian(xlim = c(13, 80)) +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  xlab("Maximum height (m)")
wa_height



### NE AUS
ne_aus <- within %>% 
  filter(region == "Northeastern Australia") |> 
  mutate(size_int = fct_relevel(size_int, 
                                "Small ↔ Small", "Small ↔ Large",  "Large ↔ Large", )) |> 
  filter(species_from %in% c("Eucalyptus pilularis","Eucalyptus grandis", "Eucalyptus microcorys",
           "Syncarpia glomulifera","Corymbia intermedia", "Allocasuarina torulosa", "Caldcluvia paniculosa",
           "Schizomeria ovata","Ceratopetalum apetalum",  "Acacia melanoxylon", "Archontophoenix cunninghamiana",
           "Cryptocarya glaucescens", "Cryptocarya rigida", "Cetatopetalum apetalum")) |>
   mutate(species_from = factor(species_from, levels = c( "Eucalyptus pilularis","Eucalyptus grandis",
                                                         "Eucalyptus microcorys","Syncarpia glomulifera",
                                                        "Corymbia intermedia", "Allocasuarina torulosa",  "Caldcluvia paniculosa",
                                                         "Schizomeria ovata","Ceratopetalum apetalum",  "Acacia melanoxylon", "Archontophoenix cunninghamiana",
                                                         "Cryptocarya glaucescens", "Cryptocarya rigida", "Cetatopetalum apetalum"))) %>%
   mutate(species_from = fct_rev(species_from)) %>%
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion2,
             fill = fill_col)) + 
  geom_vline(xintercept = 0, 
             colour = "black", 
             linewidth = 0.2,
             linetype = "dotted") + 
  geom_hline(yintercept = 8.5, colour = "black", linetype = "dashed") + 
  geom_quasirandom(shape = 21, stroke = 1, width = 0.25, alpha = 0.7, size = 2, orientation = "y") +
  facet_wrap(~size_int) +
  scale_color_manual(values = c("#30123BFF", "#3E9BFEFF" ), 
                     guide = "none") +
  scale_fill_manual(values = c("#30123BFF", "#3E9BFEFF"), 
                    na.value = "white", 
                    name = "Region", 
                    na.translate = FALSE) +
  theme_bw() + 
  coord_cartesian(xlim = c(-1.5, 1)) +
  scale_x_continuous(breaks = c(-1.5, -1, -.5, 0, .5, 1)) +
  ylab("") + 
  xlab("") + 
  ggtitle("a. Northeastern Australia")

ne_aus

ne_height <- height %>% 
  filter(georegion == "NE_AUS") %>% 
  filter(Genus_Species %in% c("Eucalyptus pilularis","Eucalyptus grandis", "Eucalyptus microcorys",
                              "Syncarpia glomulifera","Corymbia intermedia", "Allocasuarina torulosa", "Caldcluvia paniculosa",
                              "Schizomeria ovata","Ceratopetalum apetalum",  "Acacia melanoxylon", "Archontophoenix cunninghamiana",
                              "Cryptocarya glaucescens", "Cryptocarya rigida", "Cetatopetalum apetalum")) %>% 
  mutate(taxon_name = factor(Genus_Species,
                             levels = c("Eucalyptus pilularis","Eucalyptus grandis", "Eucalyptus microcorys",
                                        "Syncarpia glomulifera","Corymbia intermedia", "Allocasuarina torulosa","Caldcluvia paniculosa",
                                        "Schizomeria ovata","Ceratopetalum apetalum",  "Acacia melanoxylon", "Archontophoenix cunninghamiana",
                                        "Cryptocarya glaucescens", "Cryptocarya rigida", "Cetatopetalum apetalum"))) %>%
  mutate(taxon_name = fct_rev(taxon_name)) %>%
  ggplot(aes(y = taxon_name)) +
  geom_hline(yintercept = 9.5, colour = "black", linetype = "dashed") + 
  geom_boxplot(aes(xmin = q75, xlower = q75, xmiddle = q75, xupper = q75, xmax = max), 
               width = 0.6, 
               stat = "identity") +
  coord_cartesian(xlim = c(13, 80)) +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 
ne_height


#### Stick all together 
## NOTE!!!!!! DO NOT NEST PATCHWORKS, CAUSES BIG ISSUES 


fig3 <- ne_aus + ne_height + se_aus + se_height + wa_aus + wa_height + 
  plot_layout(nrow = 3, 
              ncol= 2, 
              widths = c(3, 1), 
              guides = "collect", 
              axes = "collect_x") & 
  theme(legend.position = "bottom", 
        legend.spacing = unit(0, "pt")) &
  guides(fill = guide_legend(nrow = 2, ncol = 3))


fig3

#And save out
ggsave("fig3.png", 
       width = 16, 
       height = 23, 
       scale = 1.5,
       units = "cm", dpi = 300)










#### Get some max heights using austraits
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
  
  
  