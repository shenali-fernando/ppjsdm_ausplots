library(stringr)
library(tidyr)
library(dplyr)

### FG + SIZE 
fg_size <- read.csv("scripts/model_specifications/fg_size/fg_size_df_t15_updated.csv")



### SPECIES + SIZE 
sp_size <- read.csv("scripts/model_specifications/species_diameter/sp_size_df_t15_misc_updated.csv")



#add sig cols
fg_size <- fg_size %>% mutate(sig = ifelse((lo > 0 | hi < 0), 1, NA))
sp_size <- sp_size %>% mutate(sig = ifelse((lo > 0 | hi < 0), 1, NA))


#median change in coefficients 
# need to make long dfs first 
#filters
f <- fg_size |> 
  filter(class_from == class_to) |> 
  mutate(class_site = paste0(class_from, sep = "_", site)) |>  
  select(class_site, alpha, size_int) |> 
  pivot_wider(
    id_cols = class_site,
    names_from = size_int,
    values_from = alpha
  ) 

sub_f <- f |> filter(str_starts(class_site, "Subcanopy"))
median(abs(sub_f$`small small` - sub_f$`large large`), na.rm = TRUE)
median(abs(sub_f$`small large` - sub_f$`large large`), na.rm = TRUE)

can_f <- f |> filter(str_starts(class_site, "Canopy"))
median(abs(can_f$`small small` - can_f$`large large`), na.rm = TRUE)
median(abs(can_f$`small large` - can_f$`large large`), na.rm = TRUE)


s <- sp_size |> 
  filter(species_from == species_to) |> 
  mutate(species_site = paste0(species_from, sep = "_", site)) |>  
  select(species_site, alpha, size_int, class_from) |> 
  pivot_wider(
    id_cols = c(species_site, class_from),
    names_from = size_int,
    values_from = alpha
  ) 

sub_s <- s |> filter(class_from == "Subcanopy")
  median(abs(sub_s$small_small - sub_s$large_large), na.rm = TRUE)
  median(abs(sub_s$small_large - sub_s$large_large), na.rm = TRUE)
  
can_s <- s |> filter(class_from == "Canopy")
  median(can_s$small_small - can_s$large_large, na.rm = TRUE)
  median(can_s$small_large - can_s$large_large, na.rm = TRUE)
  