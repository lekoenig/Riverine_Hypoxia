## Spatial analyses -- river hypoxia dataset
## Created 4 June 2019
## Last updated 11 September 2019

## The objective of this script is to estimate local, stream channel slope for every site in the global hypoxia database

# Load packages:
library(dplyr)         # general data cleaning/aggregating
library(sp)            # used for geospatial analyses
library(elevatr)       # interface with global Terrain Tiles (https://registry.opendata.aws/terrain-tiles/)
library(ggplot2)       # plotting
library(geosphere)
library(raster)
library(readr)
library(rgdal)
library(rworldmap)     # visualize global data (either by country code or by 0.5-degree grid cells)

source("./R/Analysis_Functions.R")


## ----------------- Load datasets ----------------- ##
  
  # Read in USA dataset for now (will expand to global extent):
  dat.usa <- readRDS("./output/Hypoxia_dat_USA_JoinNHD_01_20190908.rds")

  # Transform coordinates to WGS84 and create new columns for geographic coordinates and slope:
  hypoxia.dat <- dat.usa %>% st_transform(.,crs=4326) %>%
             mutate(lat = st_coordinates(.)[,2],
                    lon = st_coordinates(.)[,1],
                    slope = NA)
  
  # distance over which to estimate slope (difference in elevation between two locations along river channel)
  dist <- 500
  
  
## ----------------- Loop through sites to calculate river slope ----------------- ##
  
  test <- dat.usa[c(1:10),] %>% st_drop_geometry()
  
  prj_dat <- st_crs(hypoxia.dat)$proj4string
  
  # Steps:
  # 1. Obtain the coordinates
  # 2. Retrieve the coordinates of the points located within a defined radius of the site
  # 3. Download raster DEM data for the site
  # 4. Extract the elevation of the site as well as the points along the defined radius
  # 5. Find the minum elevation within the radius that coorespond to the downstream end of the site reach
  # 6. Calculate the slope between the site and the downstream end of the site reach
  
  for(i in 1:nrow(hypoxia.dat)){
    
    #Obtain the coordinates
    coords_site <-  data.frame(x=hypoxia.dat$lon[i], y=hypoxia.dat$lat[i])
    
    #Retrieve the coordinates of the points located within a defined range around the site
    coords_corners <- destPoint(coords_site[1,], seq(0, 350, by=10), dist)
    coords_site <- coords_site %>% add_row( x=coords_corners[,1], y=coords_corners[,2])
    
    #Download a raster of DEM of that area
    dem_point <- get_elev_raster(coords_site, prj = prj_dat,z = 12, src = "aws")
    
    coords_site <-   tibble::add_column(coords_site, z=extract(dem_point, coords_site) )
    
    # Calculate the slope between the site and the downstream site
    dataset$slope[i] <- abs(min(coords_site$z[-1])-coords_site$z[1])/range
    
    print(paste("we've done", i, "slope is", round(dataset$slope[i],5)))
  }
  
  # visualize results
  plot(dem_point)
  lines(coords_site, type="p", pch=19, cex=.6)
  lines(coords_site[13,], type="p", pch=19, cex=.6, col="red")
  
  
  