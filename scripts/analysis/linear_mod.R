library(dplyr)
library(lme4)
library(glmmTMB)

#for intraspecific 

intra <- df_add3_5_t15 %>% 
  filter(species_to == species_from)

#rename class_ints because too long to read 
intra <- intra %>% 
  mutate(class_int = case_when(class_int == "large_large" ~ "L-L", 
                               class_int == "small_large" ~ "S-L", 
                               class_int == "small_small" ~ "S-S", 
                               is.na(class_int) ~ NA)
  )

intra <- intra %>% 
  rename(sint = class_int)

#look at distribution of response
range(intra$alpha)
hist(sqrt(intra$alpha))

#fit model 
mod <- lmer(alpha ~
  1 + sint + cc_from + cc_from:sint + (1|species_from) + (1|site), 
  data = intra)


summary(mod) 

plot(mod)

#look at some diagnostics to check assumptions
qqnorm(model_residuals)
qqline(model_residuals, col = "red")

model_residuals <- residuals(mod)
hist(model_residuals)



#Transform data: 
mod1 <- lmer(sqrt(alpha) ~
              1 + sint + cc_from + cc_from:sint + (1|species_from) + (1|region), 
            data = intra
)

summary(mod1)
