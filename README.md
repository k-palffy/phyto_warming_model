## **phyto_warming_model**
#### *A repository for simulated phytoplankton community dynamics using a numerical model written in R*

This repository contains R codes for a model describing the dependence of algal growth as a function of temperature and nutrient availability (Thomas et al. 2017) combined with a multispecies, multi-nutrient model (Roelke & Spatharis, 2015). The model code and the simulation script in this repository were used in the following [study](https://aslopubs.onlinelibrary.wiley.com/doi/10.1002/lno.12548):

Pálffy, K. & E. Smeti, 2024. Combined effect of warming, nutrients, and species pool size on the seasonal variability of phytoplankton composition: A modeling perspective. Limnology and Oceanography 69: 1056–1069.

List of files:

**phyto_model_functions.R** - R code of model function definitions  
monod_temp(): determines species-specific growth rate as a function of temperature and nutrient (N,P) concentrations
com_growth(): determines the temporal change in species abundances and nutrient concentrations

**phyto_simulations.R** - R script for all simulations reported in the study.  
Note: Running the nested loops is computationally demanding, so the original script was executed in a cloud environment. For testing, it’s recommended to run a single simulation first.

**species_parameters_N.csv** - data table containing species-level model parameter values for nitrogen-dependent growth (species in rows, parameters in columns)

**species_parameters_P.csv** - data table containing species-level model parameter values for phosphorus-dependent growth (species in rows, parameters in columns)

**temperature_scenario.csv** - daily temperature values over a year, representing a temperature scenario (described in the study in detail)
