library(ppjsdm)
library(spatstat)


################################################################################
########### CONDITIONAL AREA UNDER THE ROC CURVE (AUC) SCORE ###################

#This function uses the conditional Papangelou intensity to compute an AUC score 
#for each type in the fit of ppjsdm 

#window <- owin(c(0, 400), c(0, 400)) #example window for 400x400m square plot 


conditional_auc <- function(window, #supply window, must be an object of class owin
                            fit, 
                            configuration){
  
  df <- data.frame(type = character(), 
                   auc = numeric())
  
  for(focal in levels(types(configuration))) { #for each type in the configuration 

    conditional_intensity <- suppressWarnings(
      
      #documentation under plot_papangelou.default
      ppjsdm::plot_papangelou(fit, #conditional the papangelou intensity on all information included in the fit (covariates, locations, estimated alpha and beta intensity)
                              type = focal,
                              drop_type_from_configuration = TRUE, #the location of the type being predicted is NOT included in the configuration, therefore the intensity is not conditioned on locations of the type 
                              return_papangelou = TRUE, #return as im object 
                              grid_steps = c(100, 100), #think of as pixel area the intensity is calculated for 
                              nthreads = 4)) 
    
    focal_configuration <- suppressWarnings(ppp(x = configuration[focal]$x, #make configuration into a ppp object to use the auc function from spatstat
                                                y = configuration[focal]$y, 
                                                window = window))
    
    df_working <- data.frame(
      type = focal, 
      auc = auc(focal_configuration, conditional_intensity) #compute auc score into a df
    )
    
    df <- rbind(df, df_working)
  }
  
  return(df)
}

#run function 
# conditional_auc(window, 
#                 configuration,
#                 fit)
# 


#As I mentioned in our meeting, conditional AUC is not a perfect method. 
#It especially underestimates when there is strong clustering effects. This is
#because we set `drop_type_from_configuration = TRUE`, meaning that all points of 
#the focal type are dropped from configuration so that these points have to predicted 
#from interspecific interactions. Often the strongest interaction is that within a 
#species, but using this method this interaction is dropped leading to an
#underestimation of AUC. 

