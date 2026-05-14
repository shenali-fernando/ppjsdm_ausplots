library(ppjsdm)
library(dplyr)
library(ggplot2)
library(tidyr)

##Source the functions 
source("specifications/diameter/size_funs.R")
source("code/make_summary_fun.R")
source("code/build_effect_size_fun/effect_function.R")

#We need to supply fit and configuration for effect function to work
data <- read.csv("data/data_cleaned.csv")
sites <- unique(data$Site_Name)
sites_d <- sites[!sites %in% c('WaratahMix', 'Bird', 'Lardner', 'Caveside', 'Flowerdale', 'MtField')]

eff_final <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(eff_final) <- c("from", "to", "alpha", "mean_effect", "quan_0", "quan_25", "quan_50", "quan_75", "quan_100")

for(m in sites_d){
  
  print(paste0("Starting site = ", m))
  
  c <- size_sites(site = m, #single site 
                  group_type = "species_size",
                  show_size_freq = FALSE,
                  config_only = FALSE, #if TRUE, returns config only and exits function
                  threshold = 16, 
                  short_range = 10, 
                  short_model = "exponential")
  fit <- waratah$fit
  sum <- c$sum
  configuration <- waratah$config
  show(plot(configuration))
  
  e <- effect_size2(configuration, fit) #CHECK WHAT EFFECT_SIZE FUN YOU WANT
  
  e <- e %>% 
    mutate(site = "WaratahMix")
  
  
  eff_final <- rbind(eff_final, e)
}




ggplot(data = eff_final, 
       aes(x = alpha, 
           y = log(mean_effect))) +
  geom_point()


write.csv(effects, "effects_final_fixed.csv")

eff_final <- effects_computeonloc

##Add some cols for visualisation 
eff_final2 <- eff_final %>%
  mutate(class_from = str_extract(from, "\\w+$")) %>% 
  mutate(class_to = str_extract(to, "\\w+$")) 

eff_final2 <- eff_final2 %>% 
  mutate(class_int = paste0(class_from, sep = "_", class_to)) %>% 
  mutate(class_int = ifelse(class_int == "small_large", "large_small", class_int))

eff_final2 <- eff_final2 %>% 
  mutate(species_from = str_extract(from, "\\w+\\s+\\w+")) %>% 
  mutate(species_to = str_extract(to, "\\w+\\s+\\w+")) %>% 
  mutate(which_species = ifelse(species_from == species_to, "within", "between")) %>% 
  mutate(group = paste0(which_species, "_", "species", "_", class_int))



#### Visualisation 


ggplot(data = eff_final2, 
       aes(x = alpha, 
           y = log(mean_effect))) +
  facet_grid(class_int ~ which_species) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray60") + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_point(size = 2, alpha = .75) + 
  xlim(c(-.7, .7)) +
  ylim(c(-5, 5)) +
  #theme_classic() + 
  ylab("Effect") + 
  xlab("Alpha coefficient")


plot_papangelou(fit, 
                type = "Eucalyptus obliqua small", 
                show = "Eucalyptus obliqua small", 
                use_log = T)



d <- df_add3_5_t15 %>% 
  mutate(int = paste0(from, sep = "_", to)) %>% 
  select(int, site, se, sig, lo, hi, range_ci)

eff_final2 <- eff_final2 %>% mutate(int = paste0(from, sep = "_", to))
eff_final30 <- left_join(eff_final20, d, by = c("int", "site"))

full_eff <- eff_final30 %>% 
  rowwise() %>% 
  mutate(se = mean(c_across(c(20, 25)), na.rm = TRUE), 
         sig = mean(c_across(c(21, 26)), na.rm = TRUE),
         lo = mean(c_across(c(22, 27)), na.rm = TRUE), 
         hi = mean(c_across(c(23, 28)), na.rm = TRUE), 
         range_ci = mean(c_across(c(24, 29)), na.rm = TRUE)
  ) %>% 
  ungroup()

#get rid of extra cols 
full_eff <- full_eff %>% 
    dplyr::select(-(20:29))



write.csv(full_eff, "full_final_effectsize.csv")


within <- full_eff %>% filter(group %in% c("within_species_small_small", "within_species_large_large"))
bw <- full_eff %>% filter(!group %in% c("within_species_small_small", "within_species_large_large"))


w <- ggplot(data = within, 
       aes(x = alpha, 
           y = log(mean_effect), 
           fill = as.character(sig))) +
  facet_grid(class_int ~ which_species, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray60") + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_point(shape = 21, size = 3.5, alpha = .5) + 
  scale_fill_manual(values = c("#00BDCE", "white"),
                    labels = c("Y", "N"),
                    name = "Significant") +
  scale_x_continuous(breaks = seq(-3.3, 0.75, by = 0.5)) +
 # xlim(c(-.2, .2)) +
  #ylim(c(-1.25, 1.25)) +
  theme_bw() + 
  ylab("Log effect size") + 
  xlab("Alpha coefficient") + 
  ggtitle("Within species size")

w


ggsave("between_effects.png", b, dpi = 300)


b <- ggplot(data = bw, 
       aes(x = alpha, 
           y = log(mean_effect), 
           fill = factor(sig))) +
  facet_wrap(~group, scales = "free") +
  geom_point(shape = 21, size = 3, alpha = .7, colour = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray60") + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  scale_fill_manual(values = c("#00BDCE", "white"), 
                    labels = c("Y", "N"),
                    name = "Significant") +
  # xlim(c(-.2, .2)) +
  #ylim(c(-1.25, 1.25)) +
  theme_bw() + 
  ylab("Log effect size") + 
  xlab("Alpha coefficient") + 
  ggtitle("Between species size")
b



# meaningful estimate effect size = likelihood of finding inds 
# effect size is positively correlated with alpha coefs 
# though with much more noise in within type comparisons
# odd situations in sig coefs that have opposite effect size to coef 
#bw group furrowing = effect = alpha due to 1:1 but we think many are not meaningful 
#the actual important effects are varied 
#hard to be general 
#sig as a filter? 
#a coef of > 0.1 is substantial (lose 1 sig point)
#inferences are not symmetrical - strong neg alpha that are highly uncertain and so low effect due to sparsly distributed individuals 
#type 1 error and type 2 
#a more conservative approach is a alpha of 0.1 is not important effect 
# a more liberal approach is an alpha of 0.05 is not important 





               
  

### Check interactions where positive alpha but negative effect 

odd_effs <- full_eff %>% 
  filter(alpha > 0) %>% 
  filter(mean_effect < 1)
 #23 odd effects over 15 sites 
