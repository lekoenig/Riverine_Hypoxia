## The objective of this script is to link sites from the global hypoxia database (GRDO) to U.S. National Hydrography Dataset (NHD) flowlines and to global HydroSHEDS flowlines and associated HydroATLAS covariates
## Last updated 24 April 2020
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
library(foreach)         # parallelize R functions
library(doSNOW)          # parallelize R functions

# Access functions:
source("./R/Analysis_Functions.R")

# Save session info as metadata:
writeLines(capture.output(sessionInfo()), paste("./output/data_processed/Metadata",format(Sys.Date(),"%Y%m%d"),"_Rsession.txt",sep=""))
writeLines(capture.output(sf_extSoftVersion()), paste("./output/data_processed/Metadata_",format(Sys.Date(),"%Y%m%d"),"_SpatialVersion.txt",sep=""))



## ============================================================ ##
##                          Import data                         ##
## ============================================================ ##

## Import river hypoxia summary stats (compiled by J. Blaszczak)
hypoxia.dat <- read.csv("./data/GRDO_FinalSumStats_2020_03_04.csv",header=T,stringsAsFactors=FALSE)


## ============================================================ ##
##                Prep data for spatial analysis                ##
## ============================================================ ##

## Check for duplicated database entries:
which(duplicated(hypoxia.dat$SiteID))
which(duplicated(hypoxia.dat$DB_ID))

## What different CRS are represented?
unique(hypoxia.dat$Coord_Units) 

## Transform all data points to a common geodetic datum:

  ## Subset sites based on CRS of original data and reproject to WGS84:
      # NAD83 sites:
      hypoxia.dat.NAD83sub.projectWGS84 <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="NAD83"),] %>%
        # convert to spatial object:
        st_as_sf(coords=c("Longitude","Latitude"),crs=4269) %>%
        # project to WGS84: 
        st_transform(.,4326) 
      
      # WGS84 sites:
      hypoxia.dat.WGS84sub <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="WGS1984"),] %>%
        # convert to spatial object:
        st_as_sf(coords=c("Longitude","Latitude"),crs=4326) 

      # NAD27 sites:
      # note that I've noticed differences between machines in st_transform results from NAD27 only. If no shift is observed, try lwgeom::st_transform_proj:
      hypoxia.dat.NAD27sub.projectWGS84 <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="NAD27"),] %>%
        # convert to spatial object:
        st_as_sf(coords=c("Longitude","Latitude"),crs=4267) %>%
        # project to WGS84: 
        st_transform(.,4326)
        #lwgeom::st_transform_proj(.,crs=4326)      

      # WGS72 sites:
      hypoxia.dat.WGS72sub.projectWGS84 <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="WGS72"),] %>%
        # convert to spatial object:
        st_as_sf(coords=c("Longitude","Latitude"),crs=4322) %>%
        # project to WGS84: 
        st_transform(.,4326) 

      # ETRS89 sites:
      hypoxia.dat.etrs89sub.projectWGS84 <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="etrs89"),] %>%
        # convert to spatial object:
        st_as_sf(coords=c("Longitude","Latitude"),crs=4258) %>%
        # project to WGS84: 
        st_transform(.,4326) 
      
## Merge data subsets back into one dataset (all CRS = WGS84)
hypoxia.dat.merged <- rbind(hypoxia.dat.etrs89sub.projectWGS84,
                      hypoxia.dat.NAD27sub.projectWGS84,
                      hypoxia.dat.NAD83sub.projectWGS84,
                      hypoxia.dat.WGS72sub.projectWGS84,
                      hypoxia.dat.WGS84sub) 
hypoxia.dat.merged <- hypoxia.dat.merged %>% mutate("Lat_WGS84" = st_coordinates(.)[,2],
                                                    "Lon_WGS84" = st_coordinates(.)[,1])

## Check that crs reflects WGS84 and plot data: 
st_crs(hypoxia.dat.merged)
#mapview(hypoxia.dat.merged)

## Compare number of sites in merged dataset with original GRDO database:
length(hypoxia.dat$DB_ID) == length(hypoxia.dat.merged$DB_ID)


## ============================================================ ##
##                       Flag coastal sites                     ##
## ============================================================ ##  

## Find sites that are plotting over an ocean:  

  # Download the global "oceans" data from Natural Earth Data: 
  download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip", 
                destfile = './data/spatial/ne_10m_ocean.zip')
  
  # Unzip the file and read in the shapefile:
  unzip(zipfile = "./data/spatial/ne_10m_ocean.zip", exdir = "./data/spatial/ne_ocean_10m")

  # Load ocean data layer and transform crs of hypoxia data locations (using Mollweide equal area projection):
  oceans <- st_read(dsn= "./data/spatial/ne_ocean_10m/",layer="ne_10m_ocean") %>% 
            st_transform(.,crs=54009)
  pts <- hypoxia.dat.merged %>% 
         st_transform(.,crs=54009)
  
  # Identify which points fall over the ocean:
  ocean.intersect <- st_intersects(oceans,pts,sparse=FALSE)
  ocean.sites <- hypoxia.dat.merged[which(ocean.intersect=="TRUE"),] %>% st_as_sf(.)

  # Add a flag to GRDO merged dataset:
  hypoxia.dat.merged$flag <- NA
  hypoxia.dat.merged$flag[which(ocean.intersect=="TRUE")] <- "ocean"
  # How many sites intersect the ocean?
  length(which(hypoxia.dat.merged$flag=="ocean"))
  
  # remove unzipped ocean shp files to save storage space:
  unlink(x = "./data/spatial/ne_ocean_10m",recursive=T)
  
## The code chunk below is another approach for flagging sites that fall over an ocean.
## Note that this approach is perhaps more accurate for analyzing individual sites, but takes much longer to run, and a test of 10,000 sites yielded identical results as above.
  
  # Download the global "land" data from Natural Earth Data: 
  #download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip", 
  #              destfile = './data/spatial/ne_10m_land.zip')
    
  # Unzip the file and read in the shapefile:
  #unzip(zipfile = "./data/spatial/ne_10m_land.zip", exdir = "./data/spatial/ne_land_10m")
  #land <- st_read(dsn="./data/spatial/ne_land_10m/",layer="ne_10m_land")
    
  # Identify which points fall over land (using intersects.land function):
  #hypoxia.dat.merged <- hypoxia.dat.merged %>% mutate(land_flag = split(.,1:nrow(hypoxia.dat.merged)) %>% purrr::map(intersects.land) %>% unlist(.) %>% as.character(.))
  #hypoxia.dat.merged$flag <- ifelse(hypoxia.dat.merged$land_flag=="TRUE",NA,"ocean")
  

## ============================================================ ##
##                       Save merged dataset                    ##
## ============================================================ ## 

## Export merged dataset:
  
  # Write ESRI shapefile and save to output folder (optional):
  #hypoxia.dat.merged[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84")] %>%
  #                  st_write(.,"./output/data_spatial/GRDO_merged_20200317.shp")

  # Export csv:
  #hypoxia.dat.export <- hypoxia.dat.merged[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag")] %>%
  #                      st_drop_geometry(.)
  #write.csv(hypoxia.dat.export,"./output/data_processed/GRDO_merged_20200319.csv",row.names = FALSE)


## ============================================================ ##
##            U.S. sites: Join NHDPlus and StreamCat            ##
## ============================================================ ## 
  
## Subset U.S. sites:
    
  # Define bounding box for continental United States (spatial extent of EPA StreamCat: https://www.epa.gov/national-aquatic-resource-surveys/streamcat)
    
    # Load U.S. VPU Hydroregions:
    nhd.region <- read_sf(dsn="./data/spatial/",layer="USA_HydroRegions_VPU02")
    nhd.region.Lower48 <- nhd.region[-which(is.na(nhd.region$VPUID)|nhd.region$VPUID=="20"|nhd.region$VPUID=="21"|nhd.region$VPUID=="22"),] %>% st_transform(.,5070)
    nhd.region.PR <- nhd.region[which(nhd.region$VPUID=="21"),] 
    nhd.region.HI <- nhd.region[which(nhd.region$VPUID=="20"),] 
    
    # Extract bounding box extents:
    usa.bbox <- st_bbox(st_transform(nhd.region.Lower48,4326))
    HI.bbox <- st_bbox(st_transform(nhd.region.HI,4326))
    PR.bbox <- st_bbox(st_transform(nhd.region.PR,4326))
  
    # Subset contiguous USA sites:
    hypoxia.dat.merged.usa <- hypoxia.dat.merged[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag")] %>%
                              filter(.,Lat_WGS84 > usa.bbox$ymin & Lat_WGS84 < usa.bbox$ymax & Lon_WGS84 > usa.bbox$xmin & Lon_WGS84 < usa.bbox$xmax) %>%
                              # Project to albers equal area conic [https://www.epa.gov/waterdata/spatial-data-waters]:
                              st_transform(.,crs=5070)
    
    # Subset Hawaii sites:
    hypoxia.dat.merged.HI <- hypoxia.dat.merged[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag")] %>%
                             filter(.,Lat_WGS84 > HI.bbox$ymin & Lat_WGS84 < HI.bbox$ymax & Lon_WGS84 > HI.bbox$xmin & Lon_WGS84 < HI.bbox$xmax) %>%
                             # Project to Hawaii albers equal area conic [https://www.epa.gov/waterdata/spatial-data-waters]:
                             #st_transform(.,crs=102007)
                             st_transform(.,crs="+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-157 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
    
    # Subset Puerto Rico sites:
    hypoxia.dat.merged.PR <- hypoxia.dat.merged[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag")] %>%
                             filter(.,Lat_WGS84 > PR.bbox$ymin & Lat_WGS84 < PR.bbox$ymax & Lon_WGS84 > PR.bbox$xmin & Lon_WGS84 < PR.bbox$xmax) %>%
                             # Project to albers equal area conic, adjsuted for Puerto Rico [https://www.epa.gov/waterdata/spatial-data-waters]:
                             st_transform(.,crs = "+proj=aea +lat_1=8 +lat_2=18 +lat_0=3 +lon_0=-66 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
    
    # Plot usa hypoxia subset and check projection:
    #mapview(hypoxia.dat.merged.usa) 
    st_crs(hypoxia.dat.merged.usa)
    
    # How many GRDO sites fall within the USA?
    length(hypoxia.dat.merged.usa$DB_ID)
    
    # Set preferred projected coordinate system for continental U.S. and HI/PR:
    CRS.def <- st_crs(hypoxia.dat.merged.usa)
    CRS.HI <- st_crs(hypoxia.dat.merged.HI)
    CRS.PR <- st_crs(hypoxia.dat.merged.PR)
    
## Join NHDPlus_V2 **code below requires that the NHD dataset is stored locally**
  
  # Load National Hydrography Dataset (NHDPlus_V2) 
  
    # National-scale data [download here: https://www.epa.gov/waterdata/nhdplus-national-data]
    # un-comment the three lines below if flowline data has not been previously saved:
    #read.path <- "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb"
    #output.file <- "./output/data_spatial/nhdplus_flowlines_omitCoastlines.rds"
    #Read_NHDPlus(read.path,output.file)
    
    # Additional data (for HI and PR) [download here: https://www.epa.gov/waterdata/nhdplus-national-data]
    # un-comment the three lines below if flowline data has not been previously saved:
    #read.path <- "/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_HI_PR_VI_PI.gdb"
    #output.file <- "./output/data_spatial/nhdplus_flowlines_HI_PR_omitCoastlines.rds"
    #Read_NHDPlus(read.path,output.file)
    
    # **If flowline has previously been saved**:
    flowline <- readRDS("./output/data_spatial/nhdplus_flowlines_omitCoastlines.rds")
    flowline_HI_PR <- readRDS("./output/data_spatial/nhdplus_flowlines_HI_PR_omitCoastlines.rds")
    
    # Continental USA: create consistent column names 
    names(flowline)[1:(length(names(flowline))-1)] <- toupper(names(flowline)[1:(length(names(flowline))-1)])

    # HI/PR: create consistent column names
    names(flowline_HI_PR)[1:(length(names(flowline_HI_PR))-1)] <- toupper(names(flowline_HI_PR)[1:(length(names(flowline_HI_PR))-1)])
    flowline_HI <- filter(flowline_HI_PR,VPUID=="20")
    flowline_PR <- filter(flowline_HI_PR,VPUID=="21")
    rm(flowline_HI_PR)
    
    
  # Join GRDO locations with NHDPlus_V2 flowlines
      
    # Match locations to appropriate VPU hydroregion
    if(st_crs(hypoxia.dat.merged.usa)$proj4string == sf::st_crs(CRS.def)$proj4string){
      hypoxia.dat.merged.usa <- hypoxia.dat.merged.usa %>% 
        # identify VPU region:
        mutate(nhd_vpu_intsct = st_join(.,nhd.region.Lower48)$VPUID,
               nhd_vpu = nhd.region.Lower48$VPUID[st_nearest_feature(.,nhd.region.Lower48)])
    } else {paste("warning: projections for point data and hydroregion data do not match")}
   
    hypoxia.dat.merged.HI$nhd_vpu <- "20"
    hypoxia.dat.merged.PR$nhd_vpu <- "21"
    
    
    # Find the geodesic distance between each GRDO location and the closest NHDPlusV2 flowline :
      # *Note that this process will take a long time to run across the entire GRDO database (there are 4 functions in Analysis_Functions.R to choose from for joining GRDO points to NHDPlusV2 flowlines (vary by helper function, lat/lon vs. projected data, etc.)
    
    # Continental U.S.: join points to NHDPlus_V2:
    hypoxia.dat.merged.usa <- hypoxia.dat.merged.usa %>% select(-nhd_vpu_intsct)
    
      # define function to monitor progress of foreach:
      progress <- function(n) cat(sprintf("row %d is complete\n", n))
      opts <- list(progress=progress)
      # register parallel backend:
      cl <- parallel::makeCluster(4,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points to NHD flowlines (minimize geodesic distance):
      GRDO_JoinNHD_ls <- foreach(i=1:length(hypoxia.dat.merged.usa$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- link_comid4(flowline.data = flowline,points = hypoxia.dat.merged.usa[i,],vpu = hypoxia.dat.merged.usa$nhd_vpu[i])
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
      # combine list into a data frame:
      GRDO_JoinNHD <- do.call("rbind",GRDO_JoinNHD_ls)
      rm(GRDO_JoinNHD_ls)

    # Hawaii: join points to NHDPlus_V2:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points to NHD flowlines (minimize geodesic distance):
      GRDO_JoinNHD_HI_ls <- foreach(i=1:length(hypoxia.dat.merged.HI$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- link_comid4(flowline.data = flowline_HI,points = hypoxia.dat.merged.HI[i,],vpu = hypoxia.dat.merged.HI$nhd_vpu[i])
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
      # combine list into a data frame:
      GRDO_JoinNHD_HI <- do.call("rbind",GRDO_JoinNHD_HI_ls)
      rm(GRDO_JoinNHD_HI_ls)
      
    # Puerto Rico: join points to NHDPlus_V2:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points to NHD flowlines (minimize geodesic distance):
      GRDO_JoinNHD_PR_ls <- foreach(i=1:length(hypoxia.dat.merged.PR$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- link_comid4(flowline.data = flowline_PR,points = hypoxia.dat.merged.PR[i,],vpu = hypoxia.dat.merged.PR$nhd_vpu[i])
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
      # combine list into a data frame:
      GRDO_JoinNHD_PR <- do.call("rbind",GRDO_JoinNHD_PR_ls)
      rm(GRDO_JoinNHD_PR_ls)
    
      # Combine joined NHD data sets:
      GRDO_JoinNHD_all <- rbind(GRDO_JoinNHD,GRDO_JoinNHD_HI,GRDO_JoinNHD_PR)
      
      
## Join EPA StreamCat (continental United States only)
      
    # 1. Population density
    # Bring in the data (downloaded separate .zip folders by HydroRegion) and merge: 
      folder <- paste("./data/StreamCat/PopDens")
      zips <- list.files(folder,full.names=TRUE)
      StreamCat_PopDens <- do.call("rbind", lapply(zips,function(x) {
        df <- read.csv(unz(x, unzip(x, list=TRUE)$Name), header = TRUE,sep = ",")}))
      names(StreamCat_PopDens)[names(StreamCat_PopDens) == 'CatPctFull'] <- 'CatPctFull_PopDens'
      names(StreamCat_PopDens)[names(StreamCat_PopDens) == 'WsPctFull'] <- 'WsPctFull_PopDens'
      
    # 2. Road density
      folder <- paste("./data/StreamCat/RoadDens")
      zips <- list.files(folder,full.names=TRUE)
      StreamCat_RoadDens <- do.call("rbind", lapply(zips,function(x) {
        df <- read.csv(unz(x, unzip(x, list=TRUE)$Name), header = TRUE,sep = ",")}))
      names(StreamCat_RoadDens)[names(StreamCat_RoadDens) == 'CatPctFull'] <- 'CatPctFull_RoadDens'
      names(StreamCat_RoadDens)[names(StreamCat_RoadDens) == 'WsPctFull'] <- 'WsPctFull_RoadDens'
      
    # 3. Impervious surfaces
      folder <- paste("./data/StreamCat/Impervious_Surfaces")
      zips <- list.files(folder,full.names=TRUE)
      StreamCat_ImpSurf <- do.call("rbind", lapply(zips,function(x) {
        df <- read.csv(unz(x, unzip(x, list=TRUE)$Name), header = TRUE,sep = ",")}))
      names(StreamCat_ImpSurf)[names(StreamCat_ImpSurf) == 'CatPctFull'] <- 'CatPctFull_ImpSurf'
      names(StreamCat_ImpSurf)[names(StreamCat_ImpSurf) == 'WsPctFull'] <- 'WsPctFull_ImpSurf'
      
    # 4. NLCD 2011 [https://www.mrlc.gov/data/nlcd-2011-land-cover-conus-0]
      folder <- paste("./data/StreamCat/NLCD_2011")
      zips <- list.files(folder,full.names=TRUE)
      StreamCat_NLCD <- do.call("rbind", lapply(zips,function(x) {
        df <- read.csv(unz(x, unzip(x, list=TRUE)$Name), header = TRUE,sep = ",")}))
      names(StreamCat_NLCD)[names(StreamCat_NLCD) == 'CatPctFull'] <- 'CatPctFull_NLCD2011'
      names(StreamCat_NLCD)[names(StreamCat_NLCD) == 'WsPctFull'] <- 'WsPctFull_NLCD2011'
      
    # Join StreamCat covariates and merge with hypoxia data:
      StreamCat_covariates <- left_join(StreamCat_PopDens, StreamCat_RoadDens, by="COMID") %>%
                              left_join(., StreamCat_ImpSurf, by="COMID") %>%
                              left_join(.,StreamCat_NLCD,by="COMID") %>% 
                              select(COMID,CatPctFull_PopDens,WsPctFull_PopDens,PopDen2010Cat,PopDen2010Ws,
                                     CatPctFull_RoadDens,WsPctFull_RoadDens,RdDensCat,RdDensWs,
                                     CatPctFull_ImpSurf,WsPctFull_ImpSurf,PctImp2011Cat,PctImp2011Ws,
                                     CatPctFull_NLCD2011,WsPctFull_NLCD2011,PctWdWet2011Cat,PctHbWet2011Cat,PctWdWet2011Ws,PctHbWet2011Ws)

      GRDO_JoinNHD_all2 <- GRDO_JoinNHD_all %>%
                           left_join(.,StreamCat_covariates,by="COMID")
      
      # save data frame:
      saveRDS(GRDO_JoinNHD_all2,paste("./output/data_processed/GRDO_JoinNHD_",format(Sys.Date(),"%Y%m%d"),".rds",sep=""))

      # remove NHD flowlines from R environment:
      rm(flowline)
      rm(flowline_HI)
      rm(flowline_PR)
          
        
## ============================================================ ##
##                 Global sites: Join HydroATLAS                ##
## ============================================================ ## 
      
## Load RiverATLAS data (downloaded from https://www.hydrosheds.org/page/hydroatlas on 30 December 2019)
  
  # Unzip folder that contains the river flowlines:      
  #unzip("/Users/LKoenig/Documents/SpatialData/HydroATLAS/RiverATLAS_Data_v10_shp.zip",exdir="/Users/LKoenig/Documents/SpatialData/HydroATLAS/")
  unzip("./data/spatial/RiverATLAS_Data_v10_shp.zip",exdir="./data/spatial/")
  
  # Create folder to house regional joined data:
  dir.create("./output/data_processed/Regional_HydroATLAS_joins")
  
  # Define bounds for each of the shapefiles given in the HydroATLAS download:
  # [Note that Hawaii sites (n = 486) are not represented in HydroATLAS]
  shp_data <- data.frame(file = c("rvrAtlas.na","rvrAtlas.af","rvrAtlas.ar","rvrAtlas.as","rvrAtlas.au",
                                  "rvrAtlas.eu","rvrAtlas.gr","rvrAtlas.sa.north","rvrAtlas.sa.south","rvrAtlas.si"),
                         xmin = c(-137.935417,-17.99375,-179.99792,57.639583,95.10417,-24.47292,-72.66458,-91.65625,-80.7791667,59.01667),
                         ymin = c(5.510417,-34.79375,51.22708,1.272917,-54.73958,12.60208,59.80208,-20.45417,-55.8770833,45.60208),
                         xmax = c(-52.664583,54.41667,-61.11875,150.468750,179.99792,69.51667,-12.30208,-45.85625,-34.7937500,179.99792),
                         ymax = c(62.672917,37.33125,83.20625,55.900000,20.79583,81.78958,83.58958,14.87708,-0.5895833,81.18958),
                         intersect = NA)
  
  # For each site in the global GRDO database, identify which shp file intersects the site:
  continent_int <- data.frame(DB_ID = NA, number_HydroATLAS_shp_intersects = NA)
  for(i in 1:length(hypoxia.dat.merged$DB_ID)){
    location <- hypoxia.dat.merged[i,]
    for(j in 1:length(shp_data$file)){
      int <- unlist(suppressMessages(st_intersects(location,create_bounding_box(table = shp_data, group_name = shp_data$file[j],crs.code = 4326),sparse = FALSE)))
      shp_data$intersect[j] <- int
    }
    continent_int[i,"DB_ID"] <- location$DB_ID
    continent_int[i,"number_HydroATLAS_shp_intersects"] <- length(which(shp_data$intersect=="TRUE"))
    print(i)
  }
  # How many sites overlap shp regions?
  length(which(continent_int$number_HydroATLAS_shp_intersects>1))
  # Export continent_int:
  write.csv(continent_int,paste("./output/data_processed/Regional_HydroATLAS_joins/GRDO_shp_intersections_",format(Sys.Date(),"%y%m%d"),".csv",sep=""),row.names = FALSE)
  
  
## Join GRDO sites to RiverATLAS shp files:
  
  # 1. North America and Caribbean
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.na <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.na")
  
    # Load region shp file from HydroATLAS:
    rvrAtlas.na <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                           layer = "RiverATLAS_v10_na")
  
    # Join GRDO sites to HydroATLAS flowlines:
    progress <- function(n) cat(sprintf("row %d is complete\n", n))
    opts <- list(progress=progress)
      # register parallel backend:
      cl <- parallel::makeCluster(4,outfile="")
      doSNOW::registerDoSNOW(cl)
      
      # join points to NHD using nhdTools helper function:
      GRDO_JoinHydroATLAS_na_ls <- foreach(i=1:length(hypoxia.dat.merged.na$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.na[i,],rvrAtlas.data = rvrAtlas.na)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_na <- do.call("rbind",GRDO_JoinHydroATLAS_na_ls)
    rm(GRDO_JoinHydroATLAS_na_ls)
   
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_na,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_na.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.na)
   
  # 2. Africa
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.af <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.af")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.af <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                         layer = "RiverATLAS_v10_af")
  
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points to NHD using nhdTools helper function:
      GRDO_JoinHydroATLAS_af_ls <- foreach(i=1:length(hypoxia.dat.merged.af$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.af[i,],rvrAtlas.data = rvrAtlas.af)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_af <- do.call("rbind",GRDO_JoinHydroATLAS_af_ls)
    rm(GRDO_JoinHydroATLAS_af_ls)
      
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_af,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_af.rds")  
      
    # remove region shp file:
    rm(rvrAtlas.af)
    
  # 3. North American Arctic
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.ar <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.ar")
      
    # Load region shp file from HydroATLAS:
    rvrAtlas.ar <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                           layer = "RiverATLAS_v10_ar")
    
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points:
      GRDO_JoinHydroATLAS_ar_ls <- foreach(i=1:length(hypoxia.dat.merged.ar$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.ar[i,],rvrAtlas.data = rvrAtlas.ar)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_ar <- do.call("rbind",GRDO_JoinHydroATLAS_ar_ls)
    rm(GRDO_JoinHydroATLAS_ar_ls)
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_ar,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_ar.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.ar)
    
  # 4. Central and Southeast Asia
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.as <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.as")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.as <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                           layer = "RiverATLAS_v10_as")
    
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points:
      GRDO_JoinHydroATLAS_as_ls <- foreach(i=1:length(hypoxia.dat.merged.as$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.as[i,],rvrAtlas.data = rvrAtlas.as)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_as <- do.call("rbind",GRDO_JoinHydroATLAS_as_ls)
    rm(GRDO_JoinHydroATLAS_as_ls)
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_as,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_as.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.as)
    
  # 5. Australia and Oceania
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.au <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.au")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.au <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                           layer = "RiverATLAS_v10_au")
    
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points:
      GRDO_JoinHydroATLAS_au_ls <- foreach(i=1:length(hypoxia.dat.merged.au$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.au[i,],rvrAtlas.data = rvrAtlas.au)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_au <- do.call("rbind",GRDO_JoinHydroATLAS_au_ls)
    rm(GRDO_JoinHydroATLAS_au_ls)
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_au,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_au.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.au)
    
  # 6. Europe and Middle East
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.eu <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.eu")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.eu <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                           layer = "RiverATLAS_v10_eu")
    
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points:
      GRDO_JoinHydroATLAS_eu_ls <- foreach(i=1:length(hypoxia.dat.merged.eu$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.eu[i,],rvrAtlas.data = rvrAtlas.eu)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_eu <- do.call("rbind",GRDO_JoinHydroATLAS_eu_ls)
    rm(GRDO_JoinHydroATLAS_eu_ls)
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_eu,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_eu.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.eu)
    
  # 7. Greenland
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.gr <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.gr")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.gr <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                           layer = "RiverATLAS_v10_gr")
    
    # Join GRDO sites to HydroATLAS flowlines:
    GRDO_JoinHydroATLAS_gr <- foreach(i=1:length(hypoxia.dat.merged.gr$DB_ID),.combine = rbind) %do% {
      dat <- Join.RiverAtlas(location = hypoxia.dat.merged.gr[i,],rvrAtlas.data = rvrAtlas.gr)
      return(dat)
    }
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_gr,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_gr.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.gr)
    
  # 8. South America (north)
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.sa.north <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.sa.north")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.sa.north <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                                 layer = "RiverATLAS_v10_sa_north")
    
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points:
      GRDO_JoinHydroATLAS_sa_north_ls <- foreach(i=1:length(hypoxia.dat.merged.sa.north$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.sa.north[i,],rvrAtlas.data = rvrAtlas.sa.north)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_sa_north <- do.call("rbind",GRDO_JoinHydroATLAS_sa_north_ls)
    rm(GRDO_JoinHydroATLAS_sa_north_ls)
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_sa_north,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_sa_north.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.sa.north)
    
  # 9. South America (south)
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.sa.south <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.sa.south")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.sa.south <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                                 layer = "RiverATLAS_v10_sa_south")
    
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points:
      GRDO_JoinHydroATLAS_sa_south_ls <- foreach(i=1:length(hypoxia.dat.merged.sa.south$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.sa.south[i,],rvrAtlas.data = rvrAtlas.sa.south)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_sa_south <- do.call("rbind",GRDO_JoinHydroATLAS_sa_south_ls)
    rm(GRDO_JoinHydroATLAS_sa_south_ls)
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_sa_south,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_sa_south.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.sa.south)
    
  # 10. Siberia
    # Subset GRDO sites that fall within region shp file:
    hypoxia.dat.merged.si <- filter_data_HydroATLAS_shp(GRDO_data = hypoxia.dat.merged,bbox_table = shp_data,filename="rvrAtlas.si")
    
    # Load region shp file from HydroATLAS:
    rvrAtlas.si <- read_sf(dsn = "./data/spatial/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/",
                           layer = "RiverATLAS_v10_si")
  
    # Join GRDO sites to HydroATLAS flowlines:
      # register parallel backend:
      cl <- parallel::makeCluster(3,outfile="")
      doSNOW::registerDoSNOW(cl)
      # join points:
      GRDO_JoinHydroATLAS_si_ls <- foreach(i=1:length(hypoxia.dat.merged.si$DB_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
        dat <- Join.RiverAtlas(location = hypoxia.dat.merged.si[i,],rvrAtlas.data = rvrAtlas.si)
        return(dat)
      }
      # close parallel cores:
      parallel::stopCluster(cl)
      
    # combine list into a data frame:
    GRDO_JoinHydroATLAS_si <- do.call("rbind",GRDO_JoinHydroATLAS_si_ls)
    rm(GRDO_JoinHydroATLAS_si_ls)
    
    # Save join data:
    saveRDS(GRDO_JoinHydroATLAS_si,"./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_si.rds")  
    
    # remove region shp file:
    rm(rvrAtlas.si)
  
## Join regional GRDO-RiverATLAS files:
    
    # Load regional files:
    dat.na <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_na.rds")
    dat.af <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_af.rds")
    dat.ar <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_ar.rds")
    dat.as <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_as.rds")
    dat.au <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_au.rds")
    dat.eu <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_eu.rds")
    dat.gr <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_gr.rds")
    dat.sa.north <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_sa_north.rds")
    dat.sa.south <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_sa_south.rds")
    dat.si <- readRDS("./output/data_processed/Regional_HydroATLAS_joins/GRDO_JoinHydroATLAS_si.rds")
    
    # Combine regional files:
    GRDO_HydroATLAS_data <- bind_rows(dat.na,dat.af,dat.ar,dat.as,dat.au,dat.eu,
                                      dat.gr,dat.sa.north,dat.sa.south,dat.si)
    
    # Filter DB_ID's that are not duplicated across regions:
    GRDO_HydroATLAS_no_dup <- GRDO_HydroATLAS_data %>%
                              group_by(DB_ID) %>%
                              filter(n()==1)
    
    # Filter DB_ID's that are duplicated across regions:
    GRDO_HydroATLAS_dup_all <- GRDO_HydroATLAS_data[which(GRDO_HydroATLAS_data$DB_ID %in% GRDO_HydroATLAS_data$DB_ID[c(which(duplicated(GRDO_HydroATLAS_data$DB_ID)))]),]
    GRDO_HydroATLAS_dup_cut <- GRDO_HydroATLAS_dup_all %>% group_by(DB_ID) %>% slice(which.min(HydroATLAS_near_dist_m)) 
    
    # Re-combine data with HydroATLAS covariates:
    GRDO_HydroATLAS_data2 <- bind_rows(GRDO_HydroATLAS_no_dup,GRDO_HydroATLAS_dup_cut)
    
    # Choose covariate columns that are of interest for GRDO database:
    GRDO_HydroATLAS_data_export <- GRDO_HydroATLAS_data2[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag",
                              "HYRIV_ID","CATCH_SKM","UPLAND_SKM","LENGTH_KM","ORD_STRA","dis_m3_pyr","run_mm_cyr","ria_ha_csu","ele_mt_cav","ele_mt_cmn",
                              "ele_mt_cmx","slp_dg_cav","slp_dg_uav","sgr_dk_rav","tmp_dc_cyr","tmp_dc_uyr","tmp_dc_c01","tmp_dc_c02","tmp_dc_c03","tmp_dc_c04",
                              "tmp_dc_c05","tmp_dc_c06","tmp_dc_c07","tmp_dc_c08","tmp_dc_c09","tmp_dc_c10","tmp_dc_c11","tmp_dc_c12","pre_mm_cyr","pre_mm_uyr",
                              "pre_mm_c01","pre_mm_c02","pre_mm_c03","pre_mm_c04","pre_mm_c05","pre_mm_c06","pre_mm_c07","pre_mm_c08","pre_mm_c09","pre_mm_c10",
                              "pre_mm_c11","pre_mm_c12","ari_ix_cav","ari_ix_uav","glc_cl_cmj","wet_pc_cg2","wet_pc_ug2","for_pc_cse","for_pc_use","crp_pc_cse",
                              "crp_pc_use","ppd_pk_cav","ppd_pk_uav","rdd_mk_cav","rdd_mk_uav","hdi_ix_cav","urb_pc_cse","urb_pc_use")]
    
    # Check that all DB_ID's from combined data are included in export data frame:
    which(GRDO_HydroATLAS_data_export$DB_ID %in% unique(GRDO_HydroATLAS_data$DB_ID) == FALSE)
    length(GRDO_HydroATLAS_data_export$DB_ID)
    
    # Export data:
    saveRDS(GRDO_HydroATLAS_data_export,paste("./output/data_processed/GRDO_JoinHydroATLAS_",format(Sys.Date(),"%Y%m%d"),".rds",sep=""))

    # Remove unzipped RiverATLAS folder to save space:
    unlink("./data/spatial/RiverATLAS_Data_v10_shp/", recursive = TRUE)
    
  
## ============================================================ ##
##              Global sites: Calculate stream slope            ##
## ============================================================ ## 
    
library(elevatr)       # interface with global Terrain Tiles (https://registry.opendata.aws/terrain-tiles/)
library(raster)        # work with rasters in R
library(geosphere)     # geospatial analysis in R
#library(rgdal)        # geospatial analysis in R

## Download DEM's and estimate streambed slope for global GRDO sites:
  # register parallel backend:
  cl <- parallel::makeCluster(4,outfile="")
  doSNOW::registerDoSNOW(cl)  
  GRDO_slope_DEM_ls <- foreach(i=1:length(hypoxia.dat.merged$DB_ID),.packages = c("dplyr","sf","tibble","geosphere","elevatr","raster","geosphere"),
                                 .options.snow = opts,.errorhandling = 'remove') %dopar% {
                                   slope.est <- Est_streambed_slope_DEM(location = hypoxia.dat.merged[i,],slope.dist = 500,prj_dat = st_crs(hypoxia.dat.merged)$proj4string)
                                   dat <- st_drop_geometry(hypoxia.dat.merged[i,"DB_ID"])
                                   dat$slope_calc <- slope.est
                                   return(dat)
                                 }
  # close parallel cores:
  parallel::stopCluster(cl)
    
  # combine list into a data frame:
  GRDO_slope_DEM <- do.call("rbind",GRDO_slope_DEM_ls)
  rm(GRDO_slope_DEM_ls)
    
  # save slope data:
  saveRDS(GRDO_slope_DEM,paste("./output/data_processed/GRDO_slope_",format(Sys.Date(),"%Y%m%d"),".rds",sep=""))
  rm(GRDO_slope_DEM)
  
  
  
  
  
  