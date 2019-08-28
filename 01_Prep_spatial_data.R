## Spatial analyses -- river hypoxia dataset
## Created 4 June 2019
## Last updated 19 August 2019

## The objective of this script is to import the global hypoxia dataset for formatting and subsetting for further spatial analysis

# Load packages:
library(dplyr)         # general data cleaning/aggregating
library(sf)            # used for geospatial analyses
library(mapview)       # plot spatial objects
#devtools::install_github("USGS-R/nhdplusTools")  # install Blodgett nhdplusTools ()
library(nhdplusTools)  # USGS OWI package for interfacing with NHDV2 and NHDHR
library(USAboundaries)


## --------------- Import data --------------- ##

  ## Import river hypoxia summary stats (downloaded from: https://drive.google.com/drive/folders/1YL12S62g3kp02JTjq3mBIzFI7t8FsIOc)
  hypoxia.dat <- read.csv("./data/Compiled_SumStats_2019_07_15.csv",header=T,stringsAsFactors=FALSE)


## ------------ Prep data for GIS ------------ ##
  
  ## Subset unique site IDs?
  #hypoxia.dat.sp <- hypoxia.dat[-which(duplicated(hypoxia.dat$SiteID)),c(".id","SiteID","Latitude","Longitude")]
  
  ## What CRS are represented?
  unique(hypoxia.dat$Coord_Units) # NAD83, WGS84, NAD27, WGS1984,etrs89,WGS 1984, UNKWN/NA/OTHER/OLHI/Unknown
  
  ## Consolidate WGS84 notation:
  hypoxia.dat$Coord_Units[which(hypoxia.dat$Coord_Units=="WGS1984"|hypoxia.dat$Coord_Units=="WGS 1984")] <- "WGS84"

  ## Subset sites based on CRS of original data and reproject to WGS84:
    # NAD83 sites:
    hypoxia.dat.NAD83sub.projectWGS84 <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="NAD83"),] %>%
          # convert to spatial object
          st_as_sf(coords=c("Longitude","Latitude"),crs=4269) %>%
          # project to WGS84 
          st_transform(4326)
  
    # WGS84 sites:
    hypoxia.dat.WGS84sub <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="WGS84"),] %>%
          # It looks like some sites have some odd lat/lon (check NEON sites). Filter these out for now
          filter(Latitude > -90 & Latitude < 90 & Longitude > -180 & Longitude < 180 & !is.na(Latitude) & !is.na(Longitude)) %>%
          # convert to spatial object
          st_as_sf(coords=c("Longitude","Latitude"),crs=4326) 
      # Check lat/lon for NEON sites as they're beyond the bounds of WGS84
      # Check location for NWIS siteID 20210 (database object no. 27853) - currently plotting off the coast of Antarctica?
      # Check location for Africa siteID 256 (database object no. 54305) - currently plotting off the coast of Africa?
    
    # NAD27 sites:
    hypoxia.dat.NAD27sub.projectWGS84 <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="NAD27"),] %>%
          # convert to spatial object
          st_as_sf(coords=c("Longitude","Latitude"),crs=4267) %>%
          # project to WGS84 
          st_transform(4326)
      # Check location for NWIS siteID 21359 (database object no. 3770) - currently plotting off the coast of Africa?

    # ETRS89 sites:
    hypoxia.dat.etrs89sub.projectWGS84 <- hypoxia.dat[which(hypoxia.dat$Coord_Units=="etrs89"),] %>%
      # convert to spatial object
      st_as_sf(coords=c("Longitude","Latitude"),crs=4258) %>%
      # project to WGS84 
      st_transform(4326)

  ## merge data subsets back into one dataset (all CRS = WGS84)
  hypoxia.dat2 <- rbind(hypoxia.dat.etrs89sub.projectWGS84,hypoxia.dat.NAD27sub.projectWGS84,hypoxia.dat.NAD83sub.projectWGS84,hypoxia.dat.WGS84sub) 
    
  ## check that crs reflects WGS84 and plot data (note that for now this is ignoring sites where coord_units was equal to UNKWN/NA/OTHER/OLHI/Unknown (approx. 15% of the dataset):
  st_crs(hypoxia.dat2)
  mapview(hypoxia.dat2)
    
  ## write ESRI shapefile and save data to output folder (optional):
  #st_write(hypoxia.dat2, "./output/spatial/Hypoxia_dat_20190715.shp")
  saveRDS(hypoxia.dat2, "./output/Hypoxia_dat_20190812.rds")

    
## ------------ Subset U.S. sites and link to NHD reaches ------------ ##
  
  # 1. Clip hypoxia dataset to USA extent
  
  # Define bounding box for continental United States (limited by spatial extent of NHDV2, which does not include AK, combined with StreamCat, which does not include HI or PR):
  #usa.extent <- st_read("./data/spatial/usa_extent_albers.shp",crs=102003)
  usa.extent <- USAboundaries::us_boundaries(type="state", resolution = "low") %>% 
                      filter(!state_abbr %in% c("PR", "AK", "HI")) %>%
                      st_combine(.) %>%
                      # project data to albers equal area:
                      st_transform(crs = 5070)
                
                 
  # Clip hypoxia dataset to spatial extent of continental United States:
  hypoxia.dat2.usa <- hypoxia.dat2 %>%
                      mutate(Lat = st_coordinates(.)[,2],
                             Lon = st_coordinates(.)[,1]) %>%
                      # project data to albers equal area:
                      st_transform(crs=5070) %>%
                      # subset by usa extent:
                      .[usa.extent,] 
    
  # Plot usa hypoxia subset:
  mapview(hypoxia.dat2.usa) 
  st_crs(hypoxia.dat2.usa)
    
  # rename column X to represent unique database ID:
  names(hypoxia.dat2.usa)[names(hypoxia.dat2.usa) == 'X'] <- 'db_ID'
  
  # save data to output folder (optional):
  saveRDS(hypoxia.dat2.usa, "./output/Hypoxia_dat_USA_20190812.rds")            # save as RDS file
  st_write(hypoxia.dat2.usa, "./output/spatial/Hypoxia_dat_usa_20190812.shp")   # save as ESRI shapefile
  
  
  
  