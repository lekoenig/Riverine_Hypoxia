# Riverine_Hypoxia
Contains data and scripts used for the Global Extent of Riverine Hypoxia project

This repository contains 2 subfolders:
1. The "data" folder contains raw data used in the analyses. Within this folder, there are multiple files containing data:

    Compiled_SumStats_2019_07_15.csv  		  
        - Received from J. Blaszczak 15 July 2019   
        - Contains updated spreadsheet with the coordinate reference systems indicated in the "Coord_Units" column  
    
    Compiled_SumStats_2019_07_11.csv		  
        - Received from J. Blaszczak 12 July 2019  
        - Update contains correct lat/lon coordinates for Argentina sites  
    
    Compiled_SumStats_2019_04_26.csv  		    
        - Downloaded from GoogleDrive on 5 June 2019. Contains global oxygen sites database  
  
    USGS_PowellCenter_Stats.csv   			    
        - Received from J. Blaszczak 14 June 2019   
        - Contains spreadsheet with the summarized oxygen data  
        - n_time is the number of measurements per site   
        - tdiff_sec is the number of seconds between each measurement  
        - DOmgL_X5, etc. represent the different quantiles of the DO distribution per site  
        - Hyp_pr_sub2, etc. are estimates of frequency and are calculated as the number of DO measurements below the mg/L thresholds 2,3,4, and 5 divided by the total number of DO measurements for that site  
        
    Powell_NHD_sum.csv   				            
        - Received from J. Blaszczak 14 June 2019   
        - Contains spreadsheet with summary information of metabolism estimates from each site as well as information about slope, elevation, and catchment area identified from NHDPlusV2 using the COMID of each site  

2. The "R" folder contains .R files containing functions used in the main analysis script.  

The root folder the primary analysis files. "Hypoxia_spatial.R" is a script that brings in the raw data, projects the data into a single, uniform coordinate reference system (WGS84), and creates and initial plot of the data sites.  
