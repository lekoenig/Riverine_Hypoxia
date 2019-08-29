## Spatial analyses -- river hypoxia dataset
## Created 4 June 2019
## Last updated 19 August 2019

## The objective of this script is to link U.S sites to NHD reaches with a unique COMID

# Load packages:
library(dplyr)         # general data cleaning/aggregating
library(sf)            # used for geospatial analyses
library(mapview)       # plot spatial objects
#devtools::install_github("USGS-R/nhdplusTools")  # install Blodgett nhdplusTools ()
library(nhdplusTools)  # USGS OWI package for interfacing with NHDV2 and NHDHR
library(nngeo)         # find nearest neighbor flowline relative to spatial points

source("./R/Analysis_Functions.R")


## ----------------- Load datasets ----------------- ##
  
  ## Read in results of ArcGIS spatial join (for comparison with methods below):
  gis.compare <- read.csv("./output/Hypoxia_dat_usa_20190812_nhd.csv",header=TRUE,stringsAsFactors = FALSE)

  ## Read in USA subset data from 01_Prep_spatial_data.R:
  dat.usa <- readRDS("./output/Hypoxia_dat_USA_20190812.rds")
  
  ## Load National Hydrography Data (NHDV2)
  
    # 1. Access data separately by VPU hydroregion
    #nhd.region <- read_sf(dsn="/Users/LKoenig/Dropbox/UCONN/Hypoxia_HeterotrophyWorkshop/data/spatial",
    #                      layer="USA_HydroRegions_VPU02")
    nhd.region <- read_sf(dsn="C:/Users/lak17007/Documents/UCONN/StreamPULSE/Hypoxia_HeterotrophyWorkshop/data/spatial",
                          layer="USA_HydroRegions_VPU02")
    nhd.region <- nhd.region[-which(is.na(nhd.region$VPUID)|nhd.region$VPUID=="20"|nhd.region$VPUID=="21"|nhd.region$VPUID=="22"),] %>%
                    st_transform(5070)
    #nhd.path <- "/Users/LKoenig/Desktop/SpatialData/NHDPlusNationalData/Indiv_VPU"

    # 2. Access seamless national database - https://www.epa.gov/waterdata/nhdplus-national-data 
    #    warning: this dataset is very memory-intensive
      
    # Indicate where NHD national data is stored (nhdplusTools):
    nhdplusTools::nhdplus_path("/Users/LKoenig/Desktop/SpatialData/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb")  
      
    # Stage the national data so that it is easier to access (nhdplusTools). Note that this takes *a lot of memory* and so once this line has been run, can also call the saved RDS file:
    staged_data <- nhdplusTools::stage_national_data(include = "flowline",
                                                      output_path = "./output/spatial")  
      
    # Read in flowline data and filter out coastlines:
    flowline <- readRDS(staged_data$flowline)
    names(flowline)[1:10]
    flowline <- flowline[-which(flowline$FTYPE=="Coastline"),]
    saveRDS(flowline,file="./output/spatial/nhdplus_flowlines_omitCoastlines.rds")
    rm(staged_data)
    
    # If flowline has previously been saved:
    flowline <- readRDS("./output/spatial/nhdplus_flowlines_omitCoastlines.rds")
  
    
## ----------------- Join hypoxia locations with NHD flowlines ----------------- ##
     
  ## Match locations to VPU hydroregion and add X,Y coordinates in assigned coordinate reference system:
  CRS.def <- 5070   #https://spatialreference.org/ref/sr-org/epsg-5070/

  if(st_crs(dat.usa)$proj4string == st_crs(CRS.def)$proj4string){
    sites.proj <- dat.usa[c(1:5000),] %>% 
                  # subset by extent of nhd.region polygon:
                  .[nhd.region,] %>% 
                  # identify VPU region and add coordinates in correct projection:
                  mutate(vpu = st_join(.,nhd.region)$VPUID,
                         X = st_coordinates(.)[,1],
                         Y = st_coordinates(.)[,2])
  } else {
    sites.proj <- dat.usa[c(1:5000),] %>% 
                  # subset by extent of nhd.region polygon:
                  .[nhd.region,] %>% 
                  st_transform(CRS.def) %>% 
                  # identify VPU region and add coordinates in correct projection:
                  mutate(vpu = st_join(.,nhd.region)$VPUID,
                         X = st_coordinates(.)[,1],
                         Y = st_coordinates(.)[,2])
  }
  
    # How many sites were excluded because they do not fall within extent of the NHDV2 VPU dataset?
    length(dat.usa$db_ID) - length(sites.proj$db_ID)
  
  
  ## Spatial join hypoxia site locations with NHDV2 flowlines
  
    # Approach 1 - using flowline data by vpu hydroregion: 
    # note that this function will take a long time to run (outputs row number as it runs to monitor progress)
    
    nhd.dat.ls <- list()
    for(i in 1:length(sites.proj$db_ID)){
      pts <- as.data.frame(sites.proj)[i,c("db_ID","X","Y")]
      vpu <- as.character(as.data.frame(sites.proj)[i,"vpu"])
      dat <- try(link_comid2(points = pts,CRS = CRS.def,vpu=vpu,maxDist = 250))
      nhd.dat.ls[[i]] <- dat
      print(i)
    }
    
    nhd.dat <- do.call("rbind",nhd.dat.ls)
    saveRDS(nhd.dat,"./output/Hypoxia_dat_USA_JoinNHD_20190812.rds")
    
    # compare output with output from ArcGIS:
    all.equal(gis.compare$X,nhd.dat$db_ID)
    all.equal(gis.compare$COMID,nhd.dat$COMID)
    
    
    ## Approach 2 - using seamless national flowline data and R package nhdplusTools (https://github.com/USGS-R/nhdplusTools/tree/master): 
    ## note that this approach links sites to *nearest* NHDV2 COMID, without regard to a max distance
    ## discover_nhdplus_id function does not require local versions of the NHD dataset
    ## seems like discover_nhdplus_id will ignore cases where stream order!=stream calc, as well as divergences
    
    #nhd.dat.ls2 <- list()
    #for (i in seq_along(sites.proj$db_ID)){
    #  pt <- try(comid_from_point(sites.proj$db_ID[i],sites.proj$Lat[i],sites.proj$Lon[i],CRS = 4326))
    #  nhd.dat.ls2[[i]] <- pt
    #  print(i)
    #}
    
    ## The message column indicates how many COMIDs were identified
    ## - for sites with multiple COMID's, we've selected the first one as a first pass  
    #nhd.dat2 <- do.call("rbind",nhd.dat.ls2)
    #nhd.dat2[grep("Error",nhd.dat2$db_ID),] <- NA
    
    ## The USGS server seems to sporadically time out (and throw an error). My hack for now is to incorporate try() into the comid_from_point function above, and then go back and fill in locations with NA:
    #ID <- sites.proj$db_ID[which(is.na(nhd.dat2$db_ID))]
    
    #if (length(ID)>0){
    #  sites.sub <- sites.proj[sites.proj$db_ID %in% ID,]
      
    #  nhd.dat.ls2b <- list()
    #  for (i in seq_along(sites.sub$db_ID)){
    #    pt <- try(comid_from_point(sites.sub$db_ID[i],sites.sub$Lat[i],sites.sub$Lon[i],CRS = 4326))
    #    nhd.dat.ls2b[[i]] <- pt
    #    print(i)
    #  }

    #  nhd.dat.sub <- do.call("rbind",nhd.dat.ls2b)
    #  nhd.dat.sub[grep("Error",nhd.dat.sub$db_ID),] <- NA
    #  nhd.dat2 <- rbind(nhd.dat2[which(!is.na(nhd.dat2$db_ID)),],nhd.dat.sub)  
    #} 
    
   
    ## Approach 3 - use the hydrolinks package - batch process (https://github.com/lawinslow/hydrolinks): 
    ## the hydrolinks package downloads and caches NHD data and therefore doesn't call local versions of the dataset
    ## note that I modified the link_to_flowlines function a bit and now omit ftype == coastlines

    #nhd.dat3 <- link_to_flowlines2(ids = sites.proj$db_ID,lats=sites.proj$Lat,lons=sites.proj$Lon,latlonCRS=4326,buffer=250,dataset="nhdplusv2")
    #nhd.dat3 <- left_join(sites.proj,nhd.dat3,by=c("db_ID"="MATCH_ID"))
    #st_geometry(nhd.dat3) <- NULL
    
    #hydrolinks::cache_clear()
    
    ## compare output with output from ArcGIS:
    #all.equal(gis.compare$X,nhd.dat3$db_ID)
    #all.equal(gis.compare$COMID,nhd.dat3$comid)
    
    #length(which(gis.compare$COMID != nhd.dat3$comid))
    #look <- which(gis.compare$COMID != nhd.dat3$comid)
    #gis.compare[look,c("X","COMID")]
    #nhd.dat4[look,c("db_ID","comid")]
  

      
## ----------------- Join StreamCat covariate data ----------------- ##
    
  # note: can manually reference NHD VAA and StreamCat information for a given COMID using the EPA WATERS Watershed Characterization Service (https://watersgeo.epa.gov/watershedreport/?comid=4795168) by replacing comid=X at the end of the URL.
  
  # The following data sets were downloaded from EPA StreamCAT by HydroRegion and are stored in the data sub-folder:
  # 1. Population density
  #    variable name = PopDen2010Cat = mean population density (people/km^2) within local catchment
  #    variable name = PopDen2010Ws = mean population density (people/km^2) within upstream watershed
  
  # Bring in the data (downloaded separate .zip folders by HydroRegion) and merge: 
  folder <- paste("./data/StreamCat")
  zips <- list.files(folder,full.names=TRUE)
  StreamCat_PopDens <- do.call("rbind", lapply(zips,function(x) {
                                        df <- read.csv(unz(x, unzip(x, list=TRUE)$Name), header = TRUE,sep = ",")}))
  
  # Join StreamCat covariates to hypoxia data:  
  hypoxia.dat2.usa.nhd <- nhd.dat %>%
                          st_drop_geometry() %>%
                          left_join(.,StreamCat_PopDens,by="COMID")
  
  
  