# Priority site identification for invasive species management
Using spatial conservation prioritisation techniques to identify areas most sensitive to the impacts of invasive alien species. Ultimately, priority sites are those sites that are determined to be both _sensitive_ and _susceptible_ to the occurrence and establishment of invasive alien species.

Clarke, D.A., Clarke, R.H. and McGeoch, M.A. (2025), How to Identify Priority Sites for Invasive Alien Species Policy and Management. Divers Distrib, 31: e13970. https://doi.org/10.1111/ddi.13970

## Note:
This code was originally written a couple of years agao. As such, I am in the process of updating all the code, where relevant. This is especially the case regarding spatial analyses and the move from `sp` and `raster` to `sf` and `terra`.

Also, given the amount and spatial extent of the data used, these analyses cannot be fully replicated in a short period of time. Additionally, some of the data used in the analyses will not be provided. This includes the IUCN Red List data, including assessment and spatial information, and the Key Biodiversity Area (KBA) shapefile. If you wish to access this data, you will need to make a request from the data sources. 

### Required R packages
The R package `bossMaps`, which I use to create bias files for use in Maxent species distribution modelling, is no longer available via CRAN. As such, the archived version will need to be installed, e.g., 

`packageurl <- "[https://cran.r-project.org/src/contrib/Archive/bossMaps/bossMaps_0.1.0.tar.gz](https://cran.r-project.org/src/contrib/Archive/bossMaps/bossMaps_0.1.0.tar.gz)"` 

`install.packages(packageurl, repos=NULL, type="source")`
