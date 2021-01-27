# Riverine_Hypoxia
Contains data and scripts used for the Global Extent of Riverine Hypoxia project

**This repository contains 3 subfolders:**
1. The "data" folder contains raw data used in the analyses. Within this folder, there are multiple files containing data:

    GRDO_FinalSumStats_2020_06_09.csv  
    	- received from J. Blaszczak on 9 July 2020  
    	- contains a few fixes to lat/lon coordinates for sites within the continental U.S., and includes new data sets from Canada that were missing from the previous summary file.  
  
2. The "R" folder contains .R files that store any functions used in the main analysis scripts.  

3. The "output" folder contains any derived datasets or figures that were generated during the analysis.  

**The root folder contains the primary analysis files:**  

"Hypoxia_SpatialAnalyses.R" is a script that:  
    - brings in the raw data    
    - projects to a uniform coordinate reference system  
    - spatially joins sites within the United States to a unique NHDV2 COMID, which is then linked to the StreamCat covariate database (https://www.epa.gov/national-aquatic-resource-surveys/streamcat)  
    - spatially joins sites around the globe to a unique HydroATLAS flowline (https://www.hydrosheds.org/page/hydroatlas)  
    - estimates stream channel slope for each site
    
**Note that to run Hypoxia_SpatialAnalyses.R, you will need to download HydroATLAS data and store locally.**  
    - Download folder "RiverATLAS_Data_V10_shp.zip" [https://figshare.com/articles/HydroATLAS_version_1_0/9890531] and store in /data/spatial subfolder.  

