library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)

#Load cleaned data 
data_cleaned <- read.csv("data/data_cleaned.csv")

eucs <- data_cleaned %>% 
  mutate(species = if_else(str_starts(Genus_Species, 
                                       regex("^Eucalyptus|^Corymbia"), 
                                       negate = TRUE), 
                            "Non-euc", 
                            "Eucalyptus")) %>% 
  filter(species == "Eucalyptus") 

eucs %>% count(Crown_Class)

summary(eucs$Diameter)


com_eucs <- eucs %>% 
  mutate(new_cc = if_else(Crown_Class == "Emergent", 
                          "Co/dominant", 
                          Crown_Class)) %>% 
  mutate(new_cc = if_else(new_cc %in% c("Co-dominant", "Dominant", "Co/dominant"), 
                          "Co/dominant", 
                          Crown_Class)) %>% 
  filter(!new_cc == "")
  


com_eucs %>% count(new_cc)


ggplot(com_eucs) + 
  geom_point(aes(x = new_cc, y = Diameter)) 




### Histogram 

## make diameter bins first 
com_eucs <- com_eucs %>%
  mutate(diameter_bins = cut_number(Diameter, n = 30)) %>% 
  group_by(diameter_bins) %>% 
  mutate(count = n()) %>% 
  ungroup()

ss <- com_eucs %>% filter(new_cc == "Suppressed") %>% 
  count(diameter_bins)

ggplot(com_eucs, aes(fill=new_cc, y=count, x=diameter_bins)) + 
  geom_bar(position="stack", stat="identity")


ggplot(com_eucs, aes(x = Diameter, fill = new_cc)) +            
  geom_histogram( bins = 50, 
                 position = "stack") + 
  scale_fill_manual(values = c("black", "gray30", "gray60"))


library(ggpubr)
library(cowplot)

diam <- gghistogram(com_eucs, 
            x = "Diameter", 
            y = "count", 
            fill = "new_cc",
            bins = 50, 
            rug = TRUE,
            add_density = FALSE,
            palette = c("#E7B800", "#00AFBB", "#E85856")) + 
  xlim(0, 200)

diam


dense <- ggdensity(
  com_eucs,
  x = "Diameter",
  color= "new_cc", 
  y = "count",
  bins = 50, 
  size = 1,
  palette = c("#E7B800", "#00AFBB", "#E85856"),
  alpha = 0
) +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.axis")+
  rremove("ylab") +
  rremove("y.text") +
  rremove("y.ticks") +
  rremove("legend") + 
  rremove("labs") +
  xlim(0, 200)

dense


aligned_diam <- align_plots(diam, dense, align="hv", axis="none")
diam_plot <- ggdraw(aligned_diam[[1]]) + draw_plot(aligned_diam[[2]])
diam_plot




####### E. obliqua is the eucalypt in the sites with no crown class 
ob <- com_eucs %>% 
  filter(Genus_Species == "Eucalyptus obliqua")

diam_ob <- gghistogram(ob,
                    x = "Diameter", 
                    y = "count", 
                    fill = "new_cc",
                    bins = 50, 
                    rug = TRUE,
                    add_density = FALSE,
                    palette = c("#E7B800", "#00AFBB", "#E85856"))

diam_ob


dense_ob <- ggdensity(
  ob,
  x = "Diameter",
  color= "new_cc", 
  y = "count",
  bins = 50, 
  size = 1,
  palette = c("#E7B800", "#00AFBB", "#E85856"),
  alpha = 0
) +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.axis")+
  rremove("ylab") +
  rremove("y.text") +
  rremove("y.ticks") +
  rremove("legend") 

dense_ob


aligned_ob <- align_plots(diam_ob, dense_ob, align="hv", axis="none")
diam_ob <- ggdraw(aligned_ob[[1]]) + draw_plot(aligned_ob[[2]])
diam_ob


### 




##### Violin plot 

ggplot(com_eucs, 
            aes(x = Diameter, y = new_cc)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1) + 
  ylab("")

ggplot(ob, 
       aes(x = Diameter, y = new_cc)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1) + 
  ylab("")


#similar median between all eucs and only obliqua 


##################################################################

## Sites with no crown class : Bird Track, Dip, Flowerdale





###### ORDINAL REGRESSION
#i.e. what is the probability of an individual of diameter x to be in each crown class 
#But the site and species specificity of this is important 

library(MASS)

data <- com_eucs %>% 
  select(Diameter, new_cc)

data$new_cc <- factor(data$new_cc, 
                      levels = c( "Suppressed",  "Intermediate", "Co/dominant"),
                      ordered = TRUE)



model <- polr(new_cc ~ log10(Diameter), data = data, method = "logistic")
summary(model)


library(ordinal)

model <- clm(new_cc ~ log10(Diameter), data = data, link = "logit")
summary(model)

#Models are weird, 
#boosted regression tree, random forest model, 
#is it even worth it? 