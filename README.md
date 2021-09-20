# Priority site identification for invasive species management [in progress]
Using spatial conservation prioritisation techniques to identify areas most sensitive to the impacts of invasive alien species. Ultimately, priority sites are those sites that are determined to be both _sensitive_ and _susceptible_ to the occurrence and establishment of invasive alien species.

## Note:
Given the amount and spatial extent of the data used, these analyses cannot be replicated in a short period of time. Additionally, some of the data used in the analyses will not be provided. This includes the IUCN Red List data, including assessment and spatial information, and the Key Biodiversity Area (KBA) shapefile. If you wish to access this data, you will need to make a request from the data sources. 

### Required R packages
Many R packages are required to carry out all analyses. These can be found in various R scripts, particularly `1.project_setup.R` and `results.R`. However, the R package `bossMaps` is no longer available via CRAN. As such, the archived version will need to be installed, e.g., 

`packageurl <- "[https://cran.r-project.org/src/contrib/Archive/bossMaps/bossMaps_0.1.0.tar.gz](https://cran.r-project.org/src/contrib/Archive/bossMaps/bossMaps_0.1.0.tar.gz)"` 
`install.packages(packageurl, repos=NULL, type="source")`
