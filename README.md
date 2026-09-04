
# PPJSDM on AusPlots Dataset

### Project Description

This repo contains the analysis on the TERN AusPlots data using the saturated parwise interaction Gibbs point process (ppjsdm) model. The AusPlot dataset is a set of 48 1-ha plots in tall eucalypt forests that have been fully mapped for all trees larger than 10cm DBH. We have a number of expectations based on the large literature of spatial studies in tropical, subtropical, temperate deciduous and coniferous forests for the way forests should be spatially structured. Yet, there is little examination of forets in temperate broadleaf forests. The aim of this study was to understand the spatial structure of temperate broadleaf forests by way of modelling the interactions between trees in the AusPlots dataset using size, species and functional type (subcanopy or canopy) identity.

### File Structure

-   /data: Data inputs of the project, including various site level covariates

    -   /ausplots: The TERN AusPlots data, including cleaned data and list of functional groups of species in the data

    -   /metadata: PDFs giving metadata (data collection method, location, etc)

-   /functions: Functions written for specific use in modelling the data

-   /scripts:

    -   /data_cleaning: Scripts to clean data and extract various covariate data

    -   /model_specifications: As ppjsdm can have the individuals grouped in various ways (we refer to this as specifications), we looked at various model specifications to see if patterns in interaction varied, or if a simpler model was better than our expected best model that involved both species and size. Below the grouping of individuals for each model specification is described.

        -   /null_model.R: All individuals in a site are the same species (i.e. species does not matter)

        -   /species: Individuals in a site are grouped as their species

        -   /diameter.R: Individuals in a site are grouped into small or large by diameter

        -   /crown_class: Individuals in a site are grouped as their crown class

        -   /species_diameter: Individuals in a site are grouped into species and size class of small or large (e.g. Acacia melanoxylon small)

        -   /species_crown_class: Individuals in a site are grouped into species and crown class (e.g. Eucalyptus diversicolor intermediate)

    -   /model_testing: Various tests of model fit and effect including AIC and AUC tests. These are mostly found in the supplementary materials of the paper. 

    -   /analysis: Scripts on the analysis of the interaction coefficients outputted by the model including comparing different models of the FG + Size and Species + Size (compare_fg_sp.R), linear modelling (linear_mod.R), creating and exploring climate covariates, and                         linear modelling climate covariates with interaction coefficients (clim_fg_size_mod.qmd, clim_sp_size_mod.qmd). This folder also includes the effects csvs from the linear modelling. All scripts in this folder were used to make final figures and numbers                       for the paper. 

-   /outputs: Outputs from model specifications, subsequent analysis, and testing
