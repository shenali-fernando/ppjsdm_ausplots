##### Map Figure 
library(ozmaps)
library(ggmap)
library(ggplot2)
library(dplyr)
library(ggforce)


#Load data 
data <- read.csv("data/data_cleaned.csv")

basemap <- ggplot() +
  geom_sf(data = ozmap(), fill = "gray100") +
  geom_point(data = data, aes(x = Longitude, y = Latitude, colour = Site_Name),
             size = 2.5) + 
  xlim(c(110, 156)) + 
  theme_classic() + 
  theme(legend.position = "none")


circle <- data.frame(
  Longitude = c(116, 146.5, 144.5, 149.5, 152.5, 145.5),  # center x
  Latitude = c(-34.5, -42, -38, -37, -31, -17.5),   # center y
  radius = c(1.5, 2.15, 2.1, 1.5, 2, 1.5),        # radius in degrees
  Site_Name = c("b", "c", "c", "c", "a", "a"))

# lines <- data.frame(
#   x = circle$Longitude,
#   y = circle$Latitude,
#   xend = c(110, 160, 160, 160, 160, 160),  # arbitrary x for second pane
#   yend = c(-34.5, -40, -40, -40, -20, -20))


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










within$class_int <- 


se_aus <- within %>% 
  filter(region %in% c("TAS", "SE_AUS")) %>%  
  group_by(species_from, class_int) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(obs > 2) %>% 
  filter(!species_to == "Pittosporum bicolor") %>%
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col)) + 
  geom_point(shape = 21, size = 3) + 
  scale_color_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                     guide = "none") + 
  scale_fill_manual(values = c("#ED90A4", "#D3A263", "#99B657",  "#00BDCE", "#94A9EC"), 
                    na.value = "white", 
                    name = "") +
  facet_wrap(~class_int) +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw()  + 
  ylab("") + 
  ggtitle("c. Southeastern Australia")


se_aus


wa <- within %>% 
  filter(region == "WA") %>% 
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, 
             fill = fill_col)) + 
  geom_point(shape = 21, size = 3) + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#C0AB52",  "#4FBF85"), guide = "none") + 
  scale_fill_manual(values = c("#C0AB52",  "#4FBF85"), 
                    na.value = "white", 
                    name = "") +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  ylab("") + 
  ggtitle("b. Western Australia")



naus <- within %>% 
  filter(region == "N_AUS") %>% 
  mutate(class_int = case_when(
    class_int == "large_large" ~ "Large ↔ Large",
    class_int == "small_large" ~ "Small ↔ Large",
    class_int == "small_small" ~ "Small ↔ Small"
  )) %>% 
  mutate(class_int = factor(class_int, levels = c("Small ↔ Small", "Small ↔ Large", "Large ↔ Large"))) %>% 
  ggplot(aes(x = alpha, 
             y = species_from, 
             colour = georegion, fill = fill_col)) + 
  geom_point(shape = 21, size = 3) + 
  facet_wrap(~class_int) +
  scale_color_manual(values = c("#E78ECD","#28BBD7"), 
                     guide = "none") +
  scale_fill_manual(values = c("#E78ECD","#28BBD7"), 
                    na.value = "white", 
                    name = "Region") +
  geom_vline(xintercept = 0, colour = "red") + 
  theme_bw() + 
  ylab("") + 
  ggtitle("a. Northeastern Australia")


f <- naus + wa + se_aus + plot_layout(nrow = 3, guides = "collect")

ggsave("fig1.png", f, 
       height = 11, width = 10, dpi = 300)

fig <- basemap2 + f + plot_layout(ncol = 1, heights = c(0.62, 2))

ggsave("figs.png", fig, 
       height = 18.5, width = 12, dpi = 300)

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