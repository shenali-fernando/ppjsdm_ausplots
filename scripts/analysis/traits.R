library(austraits)
library(APCalign)
library(dplyr)
library(stringr)
library(ggplot2)

austraits <- load_austraits(version = "7.0.0", path = "data/austraits", update = TRUE)


s <- summarise_database(austraits, "trait_name") 

lookup_trait(austraits, "height") %>% head()

traits <- "plant_height"

traits <- c("leaf_length", 
            "plant_height",
            "bark_thickness",
            "leaf_mass_per_area", 
            "leaf_N_per_dry_mass",
            "leaf_water_content_per_fresh_mass")

species <- c("Eucalyptus regnans", 
             "Eucalyptus grandis",          
             "Corymbia intermedia", 
             "Eucalyptus pilularis",
             "Cryptocarya glaucescens", 
             "Schizomeria ovata", 
             "Ceratopetalum apetalum", 
             "Archontophoenix cunninghamiana", 
             "Acacia melanoxylon", 
             "Polyscias elegans", 
             "Sloanea langii", 
             "Ackama paniculosa",
             "Allocasuarina torulosa", 
             "Allocasuarina decussata", 
             "Cryptocarya rigida", 
             "Eucalyptus cypellocarpa", 
             "Nothofagus cunninghamii", 
             "Acacia dealbata", 
             "Phyllocladus aspleniifolius", 
             "Leptospermum lanigerum", 
             "Pomaderris aspera", 
             "Pomaderris apetala", 
             "Olearia argophylla",
             "Monotoca glauca", 
             "Trymalium odoratissimum",
             "Syncarpia glomulifera subsp. glabra", 
             "Eucalyptus microcorys", 
             "Eucalyptus radiata", 
             "Eucalyptus viminalis", 
             "Eucalyptus fastigata", 
             "Eucalyptus obliqua", 
             "Eucalyptus diversicolor", 
             "Eucalyptus guilfoylei", 
             "Eucalyptus jacksonii",
             "Corymbia calophylla", 
             "Eucalyptus delegatensis", 
             "Eucalyptus ovata")





new_names <- create_taxonomic_update_lookup(species) #check names are current 



#check traits for the taxa and traits we are interested in 
# t <- austraits %>% 
#   extract_dataset(dataset_id = "Vesk")
# 
# a <- str_extract(austraits$traits$dataset_id, "Vesk")
# 
# df_traits <- as.data.frame(t$traits)
# 
# df_traits %>% count(dataset_id)

traits <- "plant_height"
species <- c("Eucalyptus urnigera")

t <- austraits %>% 
  extract_taxa(taxon_name = species) %>% 
  extract_trait(trait_names = traits)

df_traits <- as.data.frame(t$traits)

a <- data_cleaned %>% filter(Genus_Species == "Leptospermum lanigerum") %>% count(Site_Name)

#I want to make a summary table where for each species, for each trait there is the maximum, minimum, median and mean value 
#I will get resprouting v fire-killed data for literature (depends on location)

df_traits %>% count(trait_name) 

df_traits <- df_traits %>%
  mutate(value = as.numeric(value))

summary_traits <- df_traits %>% 
  group_by(taxon_name, trait_name) %>% 
  summarise(minimum = min(value, na.rm = TRUE), 
            median = median(value, na.rm = TRUE), 
            mean = mean(value, na.rm = TRUE), 
            maximum = max(value, na.rm = TRUE)) %>% 
  ungroup()

#write.csv(summary_traits, "summary_traits.csv")

#add intra-specific coefficients 

#let's do sp model first 
sp_site <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/specifications/species/output/final_sp_mod.csv")
sp_geo <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/specifications/species/output/species_geo_df.csv")

sp_site1 <- sp_site %>% 
  filter(! from == "Non-euc") %>% 
  filter(! to == "Non-euc") %>% 
  filter(from == to) %>% 
  dplyr::select(alpha, from) %>% 
  rename(taxon_name = from)

#Multiple values per species so will get a median for each species 
sp_site2 <- sp_site1 %>% 
  group_by(taxon_name) %>% 
  summarise(median = median(alpha)) %>% 
  ungroup()

summary_traits <- summary_traits %>% 
  left_join(sp_site2, by = "taxon_name")

#adjust df
summary_traits <- summary_traits %>% 
  rename(sp_site_median = median.y)


sp_geo1 <- sp_geo %>% 
  filter(from == to) %>% 
  filter(str_starts(from, "(Eucalyptus|Corymbia|Syncarpia)")) %>% 
  dplyr::select(alpha, from, Fit) %>% 
  rename(taxon_name = from)

summary_traits <- summary_traits %>% 
  left_join(sp_geo1, by = "taxon_name")

#Some NAs so have to deal with them 
summary_traits <- summary_traits %>% 
  filter(taxon_name %in% species) 
#spcc 
spcc_site <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/specifications/species_crown_class/output/final_spcc_mod.csv")

#some tidying
summary_traits <- summary_traits %>% 
  rename(region = Fit)

#somethings up with the bark thickness so letes remove for now 
summary_trait <- summary_traits %>% 
  filter(!trait_name == "bark_thickness_index")

#Let's get some correlations up 
 


#lets add within c-c coef 
spcc_site <- read.csv("C:/Users/shena/Desktop/ausplots/ppjsdm_ausplots/specifications/species_crown_class/output/final_spcc_mod.csv")

spcc_site1 <- spcc_site %>% 
  filter(!species_from == "Non-euc") %>% 
  filter(! species_to == "Non-euc") %>% 
  filter(species_from == species_to) %>% 
  filter(class_from == "Co/dominant") %>% 
  filter(class_to == "Co/dominant") 

spcc_site2 <- spcc_site1 %>% 
  dplyr::select(alpha, species_from)
#Multiple values per species so will get a median for each species 
spcc_site2 <- spcc_site2 %>% 
  group_by(species_from) %>% 
  summarise(median = median(alpha)) %>% 
  ungroup()

spcc_site2 <- spcc_site2 %>% rename(taxon_name = species_from)

summary_traits <- summary_traits %>% 
  left_join(spcc_site2, by = "taxon_name")

summary_trait <- summary_traits %>% 
  filter(!trait_name == "bark_thickness_index")

summary_traits <- summary_traits %>% rename(spcc_site_cd = median)


ggplot(data = summary_trait, 
       aes(x = sp_geo, 
           y = spcc_site_cd, 
           colour = region)) + 
  geom_point() + 
  geom_vline(xintercept = 0, colour = "red") +
  facet_wrap(~trait_name, scales = "free_y") + 
  ylab("trait value") + 
  xlab("Within-C/D alpha coefficient")
