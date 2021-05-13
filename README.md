# Prey encounters and spatial memory influence Area Restricted Search behaviour in a marine central place forager
R codes to reproduce the analysis presented in Iorio-Merlo et al. xxxx published in xxx. add link
The scripts provided allow to repeat the analysis presented in the publication starting from the the data available at dryad....

## R and packages
The code was developed using R version 4.0.2.  
The pakages used are: 
[adehabitatHR](https://rdocumentation.org/packages/adehabitatHR/versions/0.4.19),
[flextable](https://rdocumentation.org/packages/flextable/versions/0.6.5), 
[geosphere](https://rdocumentation.org/packages/geosphere/versions/1.5-10), 
[ggplot2](https://rdocumentation.org/packages/ggplot2/versions/3.3.3), 
[ggspatial](https://rdocumentation.org/packages/ggspatial/versions/1.1.5), 
[ggsn](https://rdocumentation.org/packages/ggsn/versions/0.5.0), 
[lubridate](https://rdocumentation.org/packages/lubridate/versions/1.7.10),
[mefa](https://rdocumentation.org/packages/mefa/versions/3.2-7), 
[momentuHMM](https://rdocumentation.org/packages/momentuHMM/versions/1.5.2),
[pracma](https://rdocumentation.org/packages/pracma/versions/1.9.9),
[raster](https://rdocumentation.org/packages/raster/versions/3.4-10),
[REdaS](https://rdocumentation.org/packages/REdaS/versions/0.9.3),
[rgdal](https://rdocumentation.org/packages/rgdal/versions/1.5-23), 
[roll](https://rdocumentation.org/packages/roll/versions/1.1.6),
[sp](https://rdocumentation.org/packages/sp/versions/1.4-5), 
[sf](https://rdocumentation.org/packages/sf/versions/0.9-8), 
[tidyverse](https://rdocumentation.org/packages/tidyverse/versions/1.3.1). 
 

## Scripts
In the folder R there are 12 scripts which should be run in the order they are numbered. Scripts 1-10 are used to prepare the data, run the analysis and generate outputs. All figures and tables are created using script 11.

* **01. Model 1 - Data preparation and dive batches**  
This code is used for the initial filtering of the dive data and the grouping of dives in 5 dive batches. It filters dives that are within 2km of the haul-out sites. Then it filters for foraging trips (i.e. round-trips to the same haul-out site longer than 12 hours) and removes trips that occurred in the first week post-tagging.

* **02. Model 1 - Dive batches state classification**  
This script is used to run the first HMM model using momentuHMM package. The analysis is run using these steps: HMM data preparation (step length and turning angle), selection of model initial parameters, model running. 

* **03. Repeatability - Overlap between consecutive months**  
This script uses the kernelUD() function from the adehabitatHR package to calculate the utilizations distribution (UD) of harbour seals while they are in an ARS state (output of model 1) for two consecutive months. We then estimate the overlap (Bhattacharaya affinity) of the distribution between the two months using the overlaphr(). we also obtain a null distribution of affinity index comparing the overlap between different individuals.

* **04. Prey encounters - At-sea accelerometer data extraction**  
This short script is used to split the large accelerometer data into smaller datasets, one for each seal trip. 

* **05. Prey encounters - PrCA axis standard deviation threshold**  
This script is used to calculate an acceleration threshold to be used in the identification of Prey catch attempts in script 07. Prey encounters - Processing accelerometer data. The code takes the standard deviation of the dynamic acceleration over a moving window of 1.5 seconds and uses a kmean cluster analysis to define a "high" and "low" threshold. See Viviant et al. 2009, Cox et al. 2018

* **06. Prey encounters - Accelerometer processing functions**  
This scripts reports a series of functions used to extract swimming effort, prey catch attempts (2 methodologies, Viviant et al.2010 and Brasseur et al. 2012) from the accelerometer.

* **07. Prey encounters - Processing accelerometer data**  
This scripts uses the functions described in script 6 and other functions to calculate dive and accelerometers parameters for each dive in the dataset.

* **08. Spatial memory - Creating memory raster**  
This script uses the output from model 1 to create raster maps of areas where the seals have displayed ARS behaviour. It uses the mid-point of each dive batch to assign dive batches to a 1x1 km grid. For each grid cell it then calculates the proportion of dive batches classifed as ARS. Cells with no dive batches are assigned a mean value.

* **09. Model 2 - Data preparation and dive batches**  
This script is very similar to script 1. It filters the data and groups the dives in 5-dive batches. For each batch the accelerometer and dive parameters are a mean of the values of the 5 dives. 

* **10. Model 2 - Dive batches state classification and covariate influence**  
This script follows the same steps of script 2 to run the HMM model. Additional steps are the import of the raster maps generated in script 8 and a covariates selection step, where the inclusion of either or both covariates is assessed.

* **11. Figures and tables**  
This script reports the codes used to create all the figures and tables presented in the paper and the supplementary materials. The codes are provided in the same order as the one in which they are presented in the paper. Most figures are plotted using outputs from previous steps in the analysis. To simplify plotting Figure 2 and S2, data subsets have been created and are provided in the "Data for figures" folder in the repository. To plot Figure 4 I have used parts of codes from the plot() function in the [momentuHMM package](https://github.com/bmcclintock/momentuHMM), which requires functions specified in the Figures - momentuHMM specific functions script.

* **Figures - momentuHMM specific functions**  
This scripts contains copies of some functions from the [momentuHMM package](https://github.com/bmcclintock/momentuHMM) which are required to estimate the stationary and transition probabilities, and their confidence intervals. These functions have been copied and not modified.

### License
GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007 (GNU GPLv3)

### References
Cox, S. L., et al. (2018). "Processing of acceleration and dive data on‐board satellite relay tags to investigate diving and foraging behaviour in free‐ranging marine predators." Methods in Ecology and Evolution 9(1): 64-77. [doi.org/10.1111/2041-210X.12845](https://doi.org/10.1111/2041-210X.12845)

Brasseur, S., et al. (2012). Habitat Preferences of Harbour Seals in the Dutch Coastal Area: Analysis and Estimate of Effects of Offshore Wind Farms (Report No. OWEZ R 252 T1 20120130 C043-10), IMARES - Wageningen UR, Noordzeewind: 58.  

McClintock, B. T. and T. Michelot (2018). "momentuHMM : R package for generalized hidden Markov models of animal movement." Methods in Ecology and Evolution 9(6): 1518-1530. [doi.org/10.1111/2041-210X.12995](https://doi.org/10.1111/2041-210X.12995)

Viviant, M., et al. (2010). "Prey capture attempts can be detected in Steller sea lions and other marine predators using accelerometers." Polar Biology 33(5): 713-719. [doi.org/10.1007/s00300-009-0750-y](https://doi.org/10.1007/s00300-009-0750-y)