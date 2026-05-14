library(ppjsdm)
library(dplyr)
library(stringr)
library(spatstat)
library(ggplot2)
library(patchwork)
library(scales)
library(austraits)
library(forcats)


data <- read.csv("data/data_cleaned.csv")

sites <- unique(data$Site_Name)

full_df <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(full_df) <- c("from", "to", "alpha", "lo", "hi", "lo_numerical",
                       "hi_numerical", "Potential","site")
for (i in sites){
  
  site_mod <- size_sites(site = i, #single site 
                         group_type = "species_size",
                         show_size_freq = FALSE,
                         config_only = FALSE, #if TRUE, returns config only and exits function
                         threshold = 13, 
                         short_range = 10, 
                         short_model = "exponential")
  
  working_df <- make_sum_df(fits = list(site_mod$fit), #use make_summary_df function to output summary as df 
                            summ = list(site_mod$sum))
  
  working_df <- working_df %>% mutate(site = i)                              
  
  full_df<- rbind(full_df, working_df)
}



full_df<- full_df %>% 
  mutate(range_ci = hi - lo)

a <- full_df %>% filter(range_ci > 5)
a %>% count(site)


noreg <- full_df


reg <- reg %>% 
  rename(reg_alpha = alpha) %>% 
  mutate(int = paste(from, sep = "_", to))

noreg <- noreg %>% 
  rename(noreg_alpha = alpha) %>% 
  mutate(int = paste(from, sep = "_", to))


check <- left_join(x = reg, y = noreg, by = c("site", "int"))
check2 <- check %>% select(site, int, noreg_alpha, reg_alpha)


ggplot(data = check2, 
       aes(x = noreg_alpha, 
           y = reg_alpha)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(intercept = 0, colour = "red") + 
  geom_point(size = 1.5, alpha = 0.5) + 
  ylab("Regularised alpha") + 
  xlab("Non-regularised alpha") + 
  xlim(-5, 1) + 
  ylim(-5, 1) 


#intra v inter 