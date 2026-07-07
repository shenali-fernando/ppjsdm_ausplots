library(terra)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ozmaps)
library(tidyterra)
library(ggplot2)
library(cowplot)

#Get data - from Australia's state for the forests report
aus <- rast("C:/Users/shena/Desktop/ausplots/basemap_aus_data/aus_for23.tif")
att <- read.csv("C:/Users/shena/Desktop/ausplots/basemap_aus_data/Aus_For23_Attributes.csv")

#So, the tif file is the shape, crs, etc 
#I don't have FOR, I want FOR_TYPE - we can use the attribute table to get this 

#Get the column we want and rename in the attribute table
att1 <- att %>% 
  select( VALUE.., FOR_TYPE)

att1 <- att1 %>% rename(VALUE = VALUE..)


############### let's get rid of everything except tall eucalypt forest 
att1 <- att1 %>% 
  mutate(FOR_TYPE = case_when(
    FOR_TYPE %in% c("Acacia", "Callitris", "Casuarina", "Melaleuca", "Eucalypt Tall Woodland", "Eucalypt Mallee Woodland", 
                    "Eucalypt Low Woodland", "Eucalypt Mallee Open", "Eucalypt Medium Woodland", "Hardwood plantation",
                    "Softwood plantation", "Mixed species plantation", "Other native forest", "Other forest - unallocated type",
                    "Mangrove", "Eucalypt Medium Closed", "Eucalypt Medium Open", "Eucalypt Low Closed", "Eucalypt Low Open", "Non-forest", 
                    "Rainforest", "Other forest", "Forest", "Non-forest", "Commercial plantation", "Native forest") ~ "Other",
    FOR_TYPE %in% c("Eucalypt Tall Closed", "Eucalypt Tall Open") ~ "Tall Eucalypt Forest",
    TRUE ~ FOR_TYPE
  ))

levels(as.factor(att1$FOR_TYPE))

#Get the FOR from the tif file 
lev <- levels(aus)[[1]]
class(lev)

#Merge the two dataframes 
lev2 <- merge(
  lev,
  att1[, c("VALUE", "FOR_TYPE")],
  by = "VALUE",
  all.x = TRUE
)

#put new raster levels in
levels(aus) <- lev2

#set FOR TYPE as a category
activeCat(aus) <- "FOR_TYPE"

#Then plot 
plot(aus)

#Kay, add state lines and change colours 
#Cool, need a border of Aus 
oz_states <- ozmap_data("states") 
oz_vect <- vect(oz_states)

#The crs of aus is in something odd, let's change it and align with oz_vect
crs(aus) #check crs
ext(aus)
aus_proj <- project(aus, "EPSG:4326")

#attempt to save this lol 
writeRaster(aus_proj, "aus_proj.tif")

#make sure oz is same crs
country_v <- project(oz_vect, crs(aus_proj))


### Let's use ggplot to make the main plot properly 
e <- ext(112, 154.5, -44.5, -10)
main <- crop(aus_proj, e)
outline <- crop(country_v, e)

## plot using terra
plot(aus_proj, 
     col = c("white", "gray25"),
     axes = TRUE,
     smooth = TRUE, 
     legend = FALSE)
lines(country_v, col = "black", lwd = 0.25)

## plot using ggplot to have more flexibility 
main_g <- ggplot() + 
  geom_spatraster(data = aus_proj,  
                  show.legend = FALSE) + 
  geom_spatvector(data = country_v, 
                  lwd = 0.15, 
                  fill = NA) +
  geom_rect(data = boxes, 
            aes(xmin = xmin, 
                xmax = xmax, 
                ymin = ymin, 
                ymax = ymax), 
            fill = NA, 
            lwd = 0.15,
            colour = "gray20") +
  annotate(geom = "text", 
           x = c(116, 143, 149.5, 151.5, 154.25, 147.25),               
           y = c(-35.75, -39, -42, -36.75, -30.75, -17.75),               
           label = c("a.", "b.", "c.", "d.", "e.", "f."), 
           color = "black", 
           size = 2.75) +
  scale_fill_manual(values = c("white", "black"), 
                    na.value ="white") + 
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.05), 
        axis.text = element_text(size = 5), 
        axis.ticks = element_line(linewidth = 0.05), 
        aspect.ratio = NULL, 
        plot.margin = unit(c(0.1, 0, 0.1, 0), "cm")) + 
  ylab("") + 
  xlab("") + 
  xlim(c(112, 155)) +
  ylim(c(-44.5, -10)) +
  coord_sf(expand = FALSE)

main_g

#make boxes 
boxes <- data.frame(
  xmin = c(114.5, 145, 152, 149, 143.25, 144.5),
  xmax = c(117, 146.5, 153.5, 150.5, 146.25, 148.5), 
  ymin = c(-35.25, -18.5, -32.75, -37.5, -39, -43.7), 
  ymax = c(-33.5, -17, -29.5, -36, -37.25, -40.5)
)

ggsave("main_aus.jpg", 
       dpi = 300, 
       width = 13.1, 
       height = 9.8, 
       units = "cm")

#############################################################################
#Kay, using the map with coloured tall eucalypt forest now need to put the plots on 
### Let's grab the coordinates of the plots 
## We can make a teraa and ggplot verisons of the plot

ausplots <- read.csv("data/ausplots/data_cleaned.csv")
ausplots <- ausplots %>% 
  select(Site_Name, Latitude, Longitude) %>% 
  distinct(Site_Name, .keep_all = T)



#Tasmania
tas <- ausplots %>% 
  filter(Longitude > 144.5 & Longitude < 148.5) %>% 
  filter(Latitude > -44 & Latitude < -40.5)

tas_points <- vect(tas, geom = c("Longitude", "Latitude"), crs = crs(aus_proj))

e_tas <- ext(144.5, 148.5, -44, -40.5)
au_tas <- crop(aus_proj, e_tas)
outline_tas <- crop(country_v, e_tas)

tas_gg <- ggplot() + 
  geom_spatraster(data = au_tas, 
                  show.legend = F) + 
  geom_spatvector(data = outline_tas, 
                  lwd = 0.5, 
                  fill = NA) +
  geom_spatvector(data = tas_points, 
                  col = "blue", 
                  cex = 6) +
  scale_fill_manual(values = c("white", "gray50"), 
                    na.value ="white") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme_void()

tas_gg

ggsave("tas.jpg", 
       dpi = 300)

### Victoria
vic <- ausplots %>% 
  filter(Longitude > 142.5 & Longitude < 147.5) %>% 
  filter(Latitude > -39 & Latitude < -37)

vic_points <- vect(vic, geom = c("Longitude", "Latitude"), crs = crs(aus_proj))

e_vic <- ext(142.5, 147.5, -39, -37)
au_vic <- crop(aus_proj, e_vic)
outline_vic <- crop(country_v, e_vic)

vic_gg <- ggplot() + 
  geom_spatraster(data = au_vic, 
                  show.legend = F) + 
  geom_spatvector(data = outline_vic, 
                  lwd = 0.5, 
                  fill = NA) +
  geom_spatvector(data = vic_points, 
                  col = "blue", 
                  cex = 5.5) +
  scale_fill_manual(values = c("white", "gray60"), 
                    na.value ="white") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void()

vic_gg

ggsave("vic.jpg", dpi = 300)


#SE-NSW 
sensw <- ausplots %>% 
  filter(Longitude > 148 & Longitude < 150.2) %>% 
  filter(Latitude > -38.25 & Latitude < -36)

sensw_points <- vect(sensw, geom = c("Longitude", "Latitude"), crs = crs(aus_proj))
e_sensw <- ext(148, 150.2, -38.25, -36)
au_sensw <- crop(aus_proj, e_sensw)
outline_sensw <- crop(country_v, e_sensw)

sensw_gg <- ggplot() + 
  geom_spatraster(data = au_sensw, 
                  show.legend = F) + 
  geom_spatvector(data = outline_sensw, 
                  lwd = 0.5, 
                  fill = NA) +
  geom_spatvector(data = sensw_points, 
                  col = "blue", 
                  cex = 3.2) +
  scale_fill_manual(values = c("white", "gray60"), 
                    na.value ="white") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void()

sensw_gg

ggsave("sensw.jpg", dpi = 300)

#Far northern NSW
nnnsw <- ausplots %>% 
  filter(Longitude > 149.5 & Longitude < 154) %>% 
  filter(Latitude > -33.5 & Latitude < -28)

nnsw_points <- vect(nnnsw, geom = c("Longitude", "Latitude"), crs = crs(aus_proj))
e_nnsw <- ext(149.5, 154, -33.5, -28)
au_nnsw <- crop(aus_proj, e_nnsw)
outline_nnsw <- crop(country_v, e_nnsw)

nnsw_gg <- ggplot() + 
  geom_spatraster(data = au_nnsw, 
                  show.legend = F) + 
  geom_spatvector(data = outline_nnsw, 
                  lwd = 0.5, 
                  fill = NA) +
  geom_spatvector(data = nnsw_points, 
                  col = "blue", 
                  cex = 4) +
  scale_fill_manual(values = c("white", "gray50"), 
                    na.value ="white") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void()

nnsw_gg

ggsave("nnsw.jpg", dpi = 300)


#QLD (wet tropics)
qld <- ausplots %>% 
  filter(Longitude > 144 & Longitude < 147) %>% 
  filter(Latitude > -19 & Latitude < -15.5)

qld_points <- vect(qld, geom = c("Longitude", "Latitude"), crs = crs(aus_proj))

e_qld <- ext(144, 147, -19, -15.5)
au_qld <- crop(aus_proj, e_qld)
outline_qld <- crop(country_v, e_qld)


qld_gg <- ggplot() + 
  geom_spatraster(data = au_qld, 
                  show.legend = F) + 
  geom_spatvector(data = outline_qld, 
                  lwd = 0.5, 
                  fill = NA) +
  geom_spatvector(data = qld_points, 
                  col = "blue", 
                  cex = 3) +
  scale_fill_manual(values = c("white", "gray50"), 
                    na.value ="white") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void()

qld_gg

ggsave("qld.jpg", dpi = 300)

#Western Australia
wa <- ausplots %>% 
  filter(Longitude > 114.5 & Longitude < 118) %>% 
  filter(Latitude > -35.5 & Latitude < -33.5)

wa_points <- vect(wa, geom = c("Longitude", "Latitude"), crs = crs(aus_proj))
e_wa <- ext(114.5, 118, -35.5, -33.5)
au_wa <- crop(aus_proj, e_wa)
outline_wa <- crop(country_v, e_wa)


wa_gg <- ggplot() + 
  geom_spatraster(data = au_wa, 
                  show.legend = F) + 
  geom_spatvector(data = outline_wa, 
                  lwd = 0.5, 
                  fill = NA) +
  geom_spatvector(data = wa_points, 
                  col = "blue", 
                  cex = 5) +
  scale_fill_manual(values = c("white", "gray60"), 
                    na.value ="white") +
  theme_void()

wa_gg

ggsave("wa.jpg", dpi = 300)
