## The objective of this script is to link sites from the global hypoxia database (GRDO) to U.S. National Hydrography Dataset (NHD) flowlines and to global HydroSHEDS flowlines and associated HydroATLAS covariates
## Last updated 6 July 2020
## LE Koenig

# Install packages if needed
#if (!require("dplyr")) install.packages("dplyr")
#if (!require("sf")) install.packages("sf")
#if (!require("mapview")) install.packages("mapview")
##devtools::install_github("USGS-R/nhdplusTools")                                      
#if (!require("nhdplusTools")) install.packages("nhdplusTools")
#if (!require("purrr")) install.packages("purrr")
#if (!require("foreach")) install.packages("foreach")
#if (!require("doSNOW")) install.packages("doSNOW")

# Load packages
library(dplyr)           # general data cleaning/aggregating
library(sf)              # used for geospatial analyses
#library(mapview)         # plot spatial objects
library(nhdplusTools)    # USGS OWI package for interfacing with NHDV2 and NHDHR
library(purrr)           # work with functions and vectors

# Access functions:
source("./R/Analysis_Functions.R")


## ============================================================ ##
##                          Import data                         ##
## ============================================================ ##

## Import river hypoxia summary stats (compiled by J. Blaszczak)
hypoxia.dat <- read.csv("./data/GRDO_FinalSumStats_2020_06_09.csv",header=T,stringsAsFactors=FALSE)
hypoxia.dat.merged <- read.csv("./output/data_processed/GRDO_merged_20200619.csv",header=T,stringsAsFactors = FALSE)

## ============================================================ ##
##                Join all spatial covariate data               ##
## ============================================================ ##

## 1. Joined NHD data:  
GRDO_NHD <- readRDS("./output/data_processed/GRDO_JoinNHD_20200623.rds")

# Select columns and rename:
GRDO_NHD <- GRDO_NHD %>% select(.,-c(X,V1,.id,SiteID,DB_Source,Lat_WGS84,Lon_WGS84,flag)) %>%
            rename(.,NHD_COMID = COMID, NHD_GNIS_NAME = GNIS_NAME, NHD_REACHCODE = REACHCODE,
                     NHD_FTYPE = FTYPE,NHD_FCODE = FCODE,NHD_STREAMCALC = STREAMCALC,
                     NHD_STREAMORDE = STREAMORDE,NHD_AREASQKM = AREASQKM,NHD_TOTDASQKM = TOTDASQKM,
                     NHD_SLOPE = SLOPE,NHD_SLOPELENKM = SLOPELENKM,NHD_QE_MA = QE_MA,
                     NHD_VE_MA = VE_MA,NHD_VPUID = VPUID,NHD_near_dist_m = near_dist_m,
                     NHD_CatPctFull_PopDens = CatPctFull_PopDens,NHD_WsPctFull_PopDens = WsPctFull_PopDens,
                     NHD_PopDen2010Cat = PopDen2010Cat,NHD_PopDen2010Ws = PopDen2010Ws,
                     NHD_CatPctFull_RoadDens = CatPctFull_RoadDens,NHD_WsPctFull_RoadDens = WsPctFull_RoadDens,
                     NHD_RdDensCat = RdDensCat,NHD_RdDensWs = RdDensWs,NHD_CatPctFull_ImpSurf = CatPctFull_ImpSurf,
                     NHD_WsPctFull_ImpSurf = WsPctFull_ImpSurf,NHD_PctImp2011Cat = PctImp2011Cat,NHD_PctImp2011Ws = PctImp2011Ws,
                     NHD_CatPctFull_NLCD2011 = CatPctFull_NLCD2011,NHD_WsPctFull_NLCD2011 = WsPctFull_NLCD2011,
                     NHD_PctWdWet2011Cat = PctWdWet2011Cat,NHD_PctHbWet2011Cat = PctHbWet2011Cat,
                     NHD_PctWdWet2011Ws = PctWdWet2011Ws,NHD_PctHbWet2011Ws = PctHbWet2011Ws)
GRDO_NHD$NHD_SLOPE <- ifelse(GRDO_NHD$NHD_SLOPE == -9998,NA,GRDO_NHD$NHD_SLOPE)

# Replace NHD columns with NA if distance to nearest NHD flowline exceeds 250 m: 
GRDO_NHD2 <- GRDO_NHD %>% mutate_at(vars("NHD_COMID":"NHD_PctHbWet2011Ws"),funs(ifelse(NHD_near_dist_m > 250, NA, .)))


## 2. Joined HydroATLAS data:
GRDO_HydroATLAS <- readRDS("./output/data_processed/GRDO_JoinHydroATLAS_20200629.rds")
dim(GRDO_HydroATLAS) # note: does not include 486 sites in Hawaii

# Replace HydroATLAS columns with NA if distance to nearest RiverATLAS flowline exceeds 250 m: 
GRDO_HydroATLAS2 <- GRDO_HydroATLAS %>% mutate_at(vars("HYRIV_ID":"urb_pc_use"),funs(ifelse(HydroATLAS_near_dist_m > 250, NA, .)))


## 3. Joined slope:
GRDO_slope <- readRDS("./output/data_processed/GRDO_slope_20200630.rds")

# Join data:
GRDO_join <- left_join(hypoxia.dat.merged[,c("DB_ID","X","V1",".id","SiteID","DB_Source","Lat_WGS84","Lon_WGS84","flag")],GRDO_slope,by="DB_ID")
GRDO_join <- left_join(GRDO_join,select(GRDO_HydroATLAS2,-c(X,V1,.id,SiteID,DB_Source,Lat_WGS84,Lon_WGS84,flag)),by="DB_ID")
GRDO_join <- left_join(GRDO_join,GRDO_NHD2,by="DB_ID")



## ============================================================ ##
##                     Create a metadata file                   ##
## ============================================================ ##

if (!require("dataMeta")) install.packages("dataMeta")
library(dataMeta)

# Create a folder to export the spatial covariates data:
dir.create(paste("./output/data_processed/GRDO_SpatialCovariates_",format(Sys.Date(),"%Y%m%d"),sep=""))

# Define variable descriptions:
var_desc <- c("Unique GRDO site identifier", 
              "X", 
              "V1", 
              ".id",
              "SiteID",
              "Original data source",
              "Site latitude (datum = WGS84)",
              "Site longitude (datum = WGS84)", 
              "Flag indicating whether site falls over an ocean or the Great Lakes",
              "Streambed slope calculated as the elevation difference divided by the distance of the reach (unitless)",
              "Joined HydroATLAS flowline ID",
              "Geodesic distance to nearest HydroATLAS flowline feature (meters)",
              "Local catchment area, HydroATLAS (km^2)",
              "Total upstream area, HydroATLAS (km^2)",
              "Flowline length, HydroATLAS (km)",
              "Strahler stream order, HydroATLAS",
              "Annual average discharge at reach pour point based on 1971-2000 average discharge values, HydroATLAS (m^3/s)",
              "Annual average land surface runoff in reach catchment based on 1971-2000 average runoff values, HydroATLAS (mm)",
              "River surface area along flowline reach, HydroATLAS (ha)",
              "Average elevation in local catchment, HydroATLAS (meters a.s.l)",
              "Minimum elevation in local catchment, HydroATLAS (meters a.s.l)",
              "Maximum elevation in local catchment, HydroATLAS (meters a.s.l)",
              "Average terrain slope for local catchment, HydroATLAS (degrees x 10)",
              "Average terrain slope for upstream watershed, HydroATLAS (degrees x 10)",
              "Average streambed slope along flowline segment, HydroATLAS (decimeters/km)",
              "Annual average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "Annual average air temperature for upstream watershed, HydroATLAS (degrees Celsius x 10)",
              "January monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "February monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "March monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "April monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "May monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "June monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "July monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "August monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "September monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "October monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "November monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "December monthly average air temperature for local catchment, HydroATLAS (degrees Celsius x 10)",
              "Annual average precipitation for local catchment, HydroATLAS (millimeters)",
              "Annual average precipitation for upstream watershed, HydroATLAS (millimeters)",
              "January monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "February monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "March monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "April monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "May monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "June monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "July monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "August monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "September monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "October monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "November monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "December monthly average precipitation for local catchment, HydroATLAS (millimeters)",
              "Average global aridity index for local catchment, HydroATLAS (index value x 100)",
              "Average global aridity index for upstream watershed, HydroATLAS (index value x 100)",
              "Dominant land cover class that represents the spatial majority of the local catchment, HydroATLAS (unitless)",
              "Wetland cover across all wetland classes within the local catchment (excluding lakes, reservoirs, and rivers), HydroATLAS (percent cover)",
              "Wetland cover across all wetland classes within the upstream watershed (excluding lakes, reservoirs, and rivers), HydroATLAS (percent cover)",
              "Forest cover within the local catchment, HydroATLAS (percent cover)",
              "Forest cover within the upstream watershed, HydroATLAS (percent cover)",
              "Cropland extent within the local catchment, HydroATLAS (percent cover)",
              "Cropland extent within the upstream watershed, HydroATLAS (percent cover)",
              "Average population density within the local catchment, HydroATLAS (people/km^2)",
              "Average population density within the upstream watershed, HydroATLAS (people/km^2)",
              "Upstream average road density within the local catchment, HydroATLAS (meters/km^2)",
              "Upstream average road density within the upstream watershed, HydroATLAS (meters/km^2)",
              "Average human development index (HDI) for the local catchment, HydroATLAS (index value x 1000)",
              "Urban extent within the local catchment, HydroATLAS (percent cover)",
              "Urban extent within the upstream watershed, HydroATLAS (percent cover)",
              #"Calculated vector processing unit (VPU) identifier",
              "Joined NHDPlusV2 flowline COMID",
              "NHDPlusV2 feature name within Geographic Names Information System",
              "NHDPlusV2 Reach Code assigned to feature. First 8 digits = subbasin code, next 6 digits = randomly assigned",
              "NHDPlusV2 feature type",
              "NHDPlusV2 feature code",
              "NHDPlusV2 stream calculator",
              "NHDPlusV2 modified Strahler Stream Order",
              "NHDPlusV2 local catchment area (km^2)",
              "NHDPlusV2 total upstream area (km^2)",
              "NHDPlusV2 flowline slope based on smoothed elevations (meters/meters)",
              "NHDPlusV2 NHDFlowline feature length used to compute slope (km)",
              "NHDPlusV2 gage-adjusted mean annual streamflow (cfs)",
              "NHDPlusV2 gage-adjusted mean annual velocity (ft/sec)",
              "NHDPlusV2 vector processing unit (VPU) Identifier",
              "Geodesic distance to nearest NHDPlusV2 flowline feature (meters)",
              "StreamCat percent of the local catchment that is covered by the pop dens layer",
              "StreamCat percent of the upstream watershed that is covered by the pop dens layer",
              "StreamCat mean population density within the local catchment (people/km^2)",
              "StreamCat mean population density within the upstream watershed (people/km^2)",
              "StreamCat percent of the local catchment that is covered by the road dens layer",
              "StreamCat percent of the upstream watershed that is covered by the road dens layer",
              "StreamCat density of roads (2010 Census Tiger Lines) within the local catchment (km/km^2)",
              "StreamCat density of roads (2010 Census Tiger Lines) within the upstream watershed (km/km^2)",
              "StreamCat percent of the local catchment that is covered by the impervious surfaces layer",
              "StreamCat percent of the upstream watershed that is covered by the impervious surfaces layer",
              "StreamCat mean imperviousness of anthropogenic surfaces (NLCD 2011) within the local catchment",
              "StreamCat mean imperviousness of anthropogenic surfaces (NLCD 2011) within the upstream watershed",
              "StreamCat percent of the local catchment that is covered by the NLCD 2011 land cover layer",
              "StreamCat percent of the upstream watershed that is covered by the NLCD 2011 land cover layer",
              "StreamCat percent of local catchment area classified as woody wetland land cover (NLCD 2011 class 90)",
              "StreamCat percent of local catchment area classified as herbaceous wetland land cover (NLCD 2011 class 95)",
              "StreamCat percent of total watershed area classified as woody wetland land cover (NLCD 2011 class 90)",
              "StreamCat percent of total watershed area classified as herbaceous wetland land cover (NLCD 2011 class 95)")

# Define variable types:              
var_type <- c(0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# Define categorical variables (variables where type == 1):
opt_descrip2 <- c(rep("",7),
                  "EU: EuropeanUnion","WQP: U.S. Water Quality Portal","PC: DB3","AF: DB4","AR: DB5","AU: DB6","CA: DB7","GEMS: DB8","GWB: DB9",
                  "HK: DB10","IR: DB11","NEON: National Ecological Observatory","SI: DB13","SP: NSF StreamPULSE project",
                  rep("",4),
                  NA,"ocean","GreatLakes",
                  rep("",2),
                  "Class 4: Tree cover, needle-leaved, evergreen","Class 2: Tree cover, broadleaved, decidious, closed",
                  "Class 14: Sparse herbaceous or sparse shrub cover","Class 12: Shrub cover, closed-open, deciduous (with or without sparse tree layer)",
                  "Class 16: Cultivated and managed areas",NA,"Class 6: Cultivated and managed areas","Class 22: Artificial surfaces and associated areas",
                  "Class 13: Herbaceous cover, closed-open","Class 9: Mosaic: tree cover/other natural vegetation",
                  "Class 11: Shrub cover, closed-open, evergreen (with or without sparse tree layer)",
                  "Class 21: Snow and ice (natural and artificial)","Class 20: Water bodies (natural and artificial)",
                  "Class 15: Regularly flooded shrub and/or herbaceous cover","Class 18: Mosaic: cropland/shrub and/or herbaceous cover",
                  "Class 1: Tree cover, broadleaved, evergreen","Class 7: Tree cover, regularly flooded, fresh",
                  "Class 3: Tree cover, regularly flooded, fresh","Class 17: Mosaic: cropland/tree cover/other natural vegetation",
                  "Class 8: Tree cover, regularly flooded, saline, (daily variation)","Class 19: Bare areas","Class 5: Tree cover, needle-leaved, deciduous",
                  "Class 23: No data",
                  rep("",12),
                  NA,"Fcode_55800: ArtificialPath","Fcode_46006: StreamRiver-Perennial","Fcode_46003: StreamRiver-Intermittent",
                  "Fcode_33400: Connector","Fcode_33600: CanalDitch","Fcode_46000: StreamRiver","Fcode_46007: StreamRiver-Ephemeral",
                  "Fcode_42807: Pipeline-general case-underground","Fcode_42801: Pipeline-aquaduct-surface","Fcode_33601: CanalDitch-aquaduct",
                  "Fcode_42809: Pipeline-penstock",
                  #"Fcode_42803: Pipeline-aquaduct-underground",
                  NA,"Ftype_ArtificialPath","Ftype_StreamRiver","Ftype_Connector","Ftype_CanalDitch","Ftype_Pipeline",
                  rep("",72)
                  )

# Create data dictionary:
linker <- build_linker(GRDO_join, variable_description = var_desc, variable_type = var_type)
data_dictionary <- build_dict(my.data = GRDO_join, linker = linker, option_description = opt_descrip2, 
                   prompt_varopts = FALSE)
data_dictionary <- data_dictionary %>% select(.,-"variable options") %>% rename(.,"variable options" = "notes")

# Add data set attributes and join attributes with data set:
data_description = "This data set contains geospatial covariate values for the Global River Dissolved Oxygen database (Blaszczak et al.)"
GRDO_join2 <- incorporate_attr(my.data = GRDO_join, data.dictionary = data_dictionary, main_string = data_description)

GRDO_join_ls <- list(names=attributes(GRDO_join2)$names,
                     class = attributes(GRDO_join2)$class,
                     main = attributes(GRDO_join2)$main,
                     dictionary = attributes(GRDO_join2)$dictionary,
                     last_edit_date = attributes(GRDO_join2)$last_edit_date,
                     author = attributes(GRDO_join2)$author,
                     data = GRDO_join)

# Export joined data:
saveRDS(GRDO_join_ls,paste("./output/data_processed/GRDO_SpatialCovariates_",format(Sys.Date(),"%Y%m%d"),"/GRDO_JoinSpatialCovariates_",format(Sys.Date(),"%Y%m%d"),".rds",sep=""))
write.csv(GRDO_join_ls$data,paste("./output/data_processed/GRDO_SpatialCovariates_200706/GRDO_JoinSpatialCovariates_",format(Sys.Date(),"%Y%m%d"),".rds",sep=""))

# Export metadata as a text file:
write.csv(GRDO_join_ls$dictionary,paste("./output/data_processed/GRDO_SpatialCovariates_",format(Sys.Date(),"%Y%m%d"),"/GRDO_JoinSpatialCovariates_DataDictionary_",format(Sys.Date(),"%Y%m%d"),".csv",sep=""),row.names = FALSE)


# Start writing to an output file
sink(paste("./output/data_processed/GRDO_SpatialCovariates_",format(Sys.Date(),"%Y%m%d"),"/Metadata_GRDO_JoinSpatialCovariates_",format(Sys.Date(),"%Y%m%d"),".txt",sep=""))

cat("=============================\n")
cat("Data set title\n")
cat("=============================\n")
GRDO_join_ls$main

cat("\n")
cat("\n")

cat("=============================\n")
cat("Author\n")
cat("=============================\n")
GRDO_join_ls$author

cat("\n")
cat("\n")

cat("=============================\n")
cat("Last edit date\n")
cat("=============================\n")
GRDO_join_ls$last_edit_date

cat("\n")
cat("\n")

cat("=============================\n")
cat("Column names\n")
cat("=============================\n")
GRDO_join_ls$names

# Stop writing to the file
sink()


