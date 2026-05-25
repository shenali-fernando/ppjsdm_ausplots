library(dplyr)


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


d <- data_cleaned %>% 
  group_by(Genus_Species, Site_Name) %>% 
  summarise(max = max(Height)) %>% 
  filter(! is.na(max))

d %>% count(Genus_Species)

d2 <- d %>% 
  filter(!is.na(Crown_Class)) %>% 
  filter(Genus_Species %in% c(
    "Synoum glandulosum")) 

d_long <- pivot_wider(d2, names_from = "Crown_Class", values_from = "n" )

d_long[is.na(d_long)] <- 0 

write.csv(d_long, "count_cc_species.csv")
