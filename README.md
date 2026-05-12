# ppjsdm_ausplots
This repo contains the analysis on the TERN AusPlots Forest Monitoring Network (arge tree survey 2012-2015) data using the saturated parwise interaction Gibbs point process (ppjsdm) model. The dataset can be freely accessed through the TERN website data portal [here](https://portal.tern.org.au/metadata/TERN/0e503109-2fb6-4182-969f-2d570abdbabd). 

The file structure is so that code for functions and exploratory analysis are found in the 'code' folder. The code, vignettes, and outputs for different model specifications (i.e. changing how individuals are grouped) are gound in the 'specification' folder. Specifications included: 
  1. Crown class: individuals in a site are grouped into only crown class and when individuals are grouped into species by crown class
  2. Species specifications are when individuals are grouped simply as their species
  3. Diameter specifications are when individuals are grouped into diameter or size classes (small or large by thresholding the dbh measurements) and when grouped into a species by size class (eg. Eucalyptus obliqua       large).

Crown class and diameter specifications were included as we know from previous work that a measure of size is important in understanding the spatial structure of trees in forests. Crown class was ultimately not chossen as the measure of size for this analysis as 3 of the 48 sites did not have crown class recorded. A simple species model was run first to understand how species varied without a size measurement. The diameter specification including both size and species identity was chosen as the best model as it gave most resolution into drivers of spatial structure. 
