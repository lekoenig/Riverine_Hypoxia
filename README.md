# Riverine_Hypoxia
Contains data and scripts used for the Global Extent of Riverine Hypoxia project

**This repository contains 3 subfolders:**
1. The "data" folder contains raw data used in the analyses. Within this folder, there are multiple files containing data:

    GRDO_FinalSumStats files
        - Received from J. Blaszczak on date given in file name
        - Contain updates and fixes to lat/lon coordinates for many sites within the continental U.S.
        - Current script uses 2020_03_04 data file

    Compiled_SumStats_2019_10_08.csv
        - Received from J.  Blaszczak 8 October 2019
        - Contains updated spreadsheet with coordinate reference systems corrected/checked  

    Compiled_SumStats_2019_07_15.csv  		  
        - Received from J. Blaszczak 15 July 2019   
        - Contains updated spreadsheet with the coordinate reference systems indicated in the "Coord_Units" column  
    
    Compiled_SumStats_2019_07_11.csv		  
        - Received from J. Blaszczak 12 July 2019
        - Update contains correct lat/lon coordinates for Argentina sites  
    
    Compiled_SumStats_2019_04_26.csv  		    
        - Downloaded from GoogleDrive on 5 June 2019. Contains global oxygen sites database

2. The "R" folder contains .R files that store any functions used in the main analysis scripts.  

3. The "output" folder contains any derived datasets or figures that were generated during the analysis.  

**The root folder contains the primary analysis files:**  

"Hypoxia_SpatialAnalyses.R" is a script that:  
    - brings in the raw data    
    - projects to a uniform coordinate reference system  
    - spatially joins sites within the United States to a unique NHDV2 COMID, which is then linked to the StreamCat covariate database (https://www.epa.gov/national-aquatic-resource-surveys/streamcat)  
    - spatially joins sites around the globe to a unique HydroATLAS flowline (https://www.hydrosheds.org/page/hydroatlas)  
    - estimates stream channel slope for each site  

