## Functions called in hypoxia project - spatial R scripts
## Last updated 19 March 2020




## ----- Function that calculates the UTM EPSG code associated with any point on earth ----- ##

lonlat2UTM = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }}


## ----- Function that determines whether a spatial point intersects the continents ----- ##

intersects.land <- function(location){  
  # note that "pt" input is a spatial point/location (class sf)
  pt <- st_as_sf(location, coords = c("Lon_WGS84", "Lat_WGS84"),crs=4326)
  
  # find the UTM projection associated with the spatial point/location (datum "WGS84"):
  pt.crs <- lonlat2UTM(st_coordinates(pt))
  
  # Adjust the size of the bounding box around the spatial point:
  new_bb = c(st_bbox(pt)$xmin-3,
             st_bbox(pt)$ymin-3, 
             st_bbox(pt)$xmax+3, 
             st_bbox(pt)$ymax+3)
  names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
  attr(new_bb, "class") = "bbox"
  attr(st_geometry(pt), "bbox") = new_bb
  
  # Convert the bounding box to a polygon:
  bb_polygon <- st_as_sfc(st_bbox(pt))
  
  # Clip the natural earth data "land" polygon to the extent of the adjusted bounding box:
  world.clipped <-st_intersection(st_geometry(land) %>% st_transform(.,3857),
                                  st_geometry(bb_polygon) %>% st_transform(.,3857))
  
  # Transform both the land and point features into the proper UTM projection:
  pt.transform <- pt %>% st_transform(.,crs=pt.crs)
  world.transform <- world.clipped %>% st_transform(.,crs=pt.crs)  
  
  # Does the point intersect land on earth? (land if "TRUE", ocean/large waterbody if "FALSE)
  land.intersect <- try(st_intersects(st_union(st_buffer(world.transform,0)),pt.transform,sparse=FALSE))
  
  return(land.intersect)    
}



## ----- Function to prepare NHDPlus data sets and write flowline data to output/data_spatial folder ----- ##

Read_NHDPlus <- function(read.path,output.file){
  
  # Define local path:
  nhdplusTools::nhdplus_path(read.path)
  
  # Stage the national data so that it is easier to access using nhdplusTools. Note that this takes a lot of memory so once this line has been run, can then call the saved RDS file:
  staged_data <- nhdplusTools::stage_national_data(include = "flowline",output_path = "./output/data_spatial")  
  
  # Read in flowline data and filter out coastlines:
  flowline <- readRDS(staged_data$flowline)
  flowline <- flowline[-which(flowline$FTYPE=="Coastline"),]
  saveRDS(flowline,file=output.file)
  rm(staged_data)                                         # remove staged_data to save memory
  unlink("./output/data_spatial/nhdplus_flowline.rds")    # remove nhdplus_flowlines data that contains coastlines
  
}



## ----- Function to match COMID with flowline data by VPU hydroregion ----- ##


##  1. this function subsets the national flowline data by vpu to speed up processing

  ## Function will return NA if distance to NHDPlusID is greater than maxDist ##
  
  link_comid <-  function (flowline.data,points, CRS, vpu, maxDist){
    
    #flowline.sub <- flowline.data[which(flowline.data$VPUID==vpu), ] %>% sf::st_transform(.,crs=CRS)
    
    if(is.na(points$nhd_vpu)){
      join.dat3 <- bind_cols(st_drop_geometry(hypoxia.dat.merged.usa[i,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag")]),
                       st_drop_geometry(flowline[1,c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")]))
      join.dat3[,c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")] <- NA  
      join.dat3$near_dist_m <- NA
    } else {
      
      flowline.sub <- flowline.data[which(flowline.data$VPUID==vpu), ] %>% sf::st_transform(.,crs=CRS)
      
      # find flowline reaches within maxDist of sample point and choose closest flowline
      join.dat <- sf::st_join(points,flowline.sub,sf::st_is_within_distance,dist=maxDist)
      
      if(length(join.dat$COMID)>1){
        # If discover_nhdplus_id returns more than one flowline segment, filter those where stream calc = stream order:
        join.dat <- join.dat[which(join.dat$STREAMCALC==join.dat$STREAMORDE),]
        dist <- st_distance(points,flowline.sub[which(flowline.sub$COMID %in% join.dat$COMID),])
        choose.dist <- which.min(dist) 
      } else {
        dist <- st_distance(points,flowline.sub[which(flowline.sub$COMID %in% join.dat$COMID),])
        choose.dist <- which.min(dist) 
      }
      
      if(!length(dist))  {
        join.dat <- sf::st_join(points,flowline.sub,sf::st_is_within_distance,dist=maxDist)
        join.dat <- join.dat[1,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag","COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")] 
        join.dat[,c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")] <- NA  
        join.dat2 <- join.dat} else {
          join.dat2 <- join.dat[choose.dist,]}
      
      join.dat2$near_dist_m <- if(!length(dist)) {NA} else {as.numeric(dist[choose.dist])}
      join.dat3 <- join.dat2[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag","COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID","near_dist_m")] %>%
        st_drop_geometry(.)
      
    }
    
    return(join.dat3)
  }


## 2. This function uses nearest neighbor calculation to do the spatial join rather than calculating a minimum distance (as in link_comid and link_comid2)

  link_comid2 <-  function (flowline.data,points, vpu){
    
    flowline.sub <- flowline.data[which(flowline.data$VPUID==vpu), ] 
    points.transform <- points %>% st_transform(.,st_crs(flowline.sub))
    
    # find the nearest neighbor flowline for each point, return a list
    nearest_flowline <- suppressMessages(nngeo::st_nn(points.transform,flowline.sub,returnDist=T,k=3,progress=FALSE)) 
      
    nearest_flowline_df <- data.frame(DB_ID = points$DB_ID,
                                      dist = t(do.call(rbind,nearest_flowline$dist)),
                                      flowline_row = t(do.call(rbind,nearest_flowline$nn))) %>%
                           distinct(.keep_all = T) %>%
                           mutate(COMID = flowline.sub$COMID[flowline_row])
  
    nearest_flowline_df2 <- flowline.sub[nearest_flowline_df$flowline_row[which.min(nearest_flowline_df$dist)],
                                         c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")]
    nearest_flowline_df2$near_dist_m <- if(!length(nearest_flowline_df$dist)) {NA} else {nearest_flowline_df$dist[which.min(nearest_flowline_df$dist)]}
    
    join.data <- cbind(st_drop_geometry(points),st_drop_geometry(nearest_flowline_df2))
    return(join.data)
    }
  

## 3. Function to match COMID with spatial points data using the nhdplusTools package

  ## Note that discover_nhdplusid finds which catchment (catchment ID) is closest to that point, and then returns the flowline COMID at the outlet of the catchment.
  ## This is different than asking - which flowline ID is closest to the sample point?
  ## get_flowline_index may be a better helper function here, is very similar to link_comid2 function above.
  
  ## Function will return NA if distance to NHDPlusID is greater than maxDist ##
  
  link_comid3 <- function(flowline.data,points,CRS,maxDist) {
  
    comid = suppressWarnings(nhdplusTools::discover_nhdplus_id(points))
  
     fline <- flowline.data[which(flowline.data$COMID %in% comid),c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")] %>%
             st_transform(.,CRS)
  
    if(length(fline$COMID)>0){
      dist <- st_distance(points,fline)
      
        if(length(dist)>1){
          # If discover_nhdplus_id returns more than one flowline segment, filter those where stream calc = stream order:
          fline <- fline[which(fline$STREAMCALC==fline$STREAMORDE),]
          
          if(length(fline$COMID)>1){
          # Find which comid minimizes the distance between the flowlines and the sample location:
          distances <- st_distance(points,fline)
          fline <- fline[which.min(distances),]
          dist <- st_distance(points,fline) } else { dist <- 9999 }
        }
    } else {dist <- 9999}
    
    if(as.numeric(dist) < maxDist){
      join.data <- cbind(st_drop_geometry(points),st_drop_geometry(fline))
      join.data$near_dist_m <- as.numeric(dist)} else{
      fline[1,c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")] <- NA   
      join.data <- cbind(st_drop_geometry(points),st_drop_geometry(fline))
      join.data[,c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")] <- NA  
      join.data$near_dist_m <- NA
      }
      
    return(join.data)
  }


## 4. This function calculates geodesic distances between GRDO points and flowlines to find nearest NHDPlusV2 flowline

  link_comid4 <-  function(flowline.data,points, vpu){
    
    flowline.sub <- flowline.data[which(flowline.data$VPUID==vpu), ] 
    #flowline.sub <- flowline.data
    points.transform <- points %>% st_transform(.,st_crs(flowline.sub))
  
    c <- st_distance(x=points.transform,y=flowline.sub,which="Great Circle",by_element = FALSE)
    join.dat <- bind_cols(points,st_drop_geometry(flowline.sub[which.min(c),c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")]))
    join.dat$near_dist_m <- as.numeric(c[which.min(c)])
    return(st_drop_geometry(join.dat))
  }
    
  
## ----- Function to filter GRDO database by the bounding box of a HydroATLAS shp file ----- ##

filter_data_HydroATLAS_shp <- function(GRDO_data,bbox_table,filename){
  data.subset <- GRDO_data[,c("X","V1",".id","SiteID","DB_Source","DB_ID","Lat_WGS84","Lon_WGS84","flag")] %>%
    filter(.,Lat_WGS84 > bbox_table$ymin[which(bbox_table$file==filename)] & 
             Lat_WGS84 < bbox_table$ymax[which(bbox_table$file==filename)] & 
             Lon_WGS84 > bbox_table$xmin[which(bbox_table$file==filename)] & 
             Lon_WGS84 < bbox_table$xmax[which(bbox_table$file==filename)])
  return(data.subset)
}


## ----- Function to create a polygon from HydroATLAS shp file bounding box ----- ##

create_bounding_box <- function(table,group_name,crs.code){
  
  bbox.polygon <- st_as_sfc(st_bbox(c(xmin = table$xmin[which(table$file==as.character(group_name))], 
                                      xmax = table$xmax[which(table$file==as.character(group_name))], 
                                      ymax = table$ymax[which(table$file==as.character(group_name))], 
                                      ymin = table$ymin[which(table$file==as.character(group_name))]), 
                                    crs = st_crs(crs.code)))
  return(bbox.polygon)
  
}


## ----- Function to join individual sites to nearest flowline within riverATLAS (using nearest geodesic distance) ----- ##

Join.RiverAtlas <- function(location,rvrAtlas.data){
  c<-st_distance(x = location,rvrAtlas.data,which="Great Circle",by_element = FALSE)
  join.dat <- bind_cols(st_drop_geometry(location),st_drop_geometry(rvrAtlas.data[which.min(c),]))
  join.dat$HydroATLAS_near_dist_m <- as.numeric(c[which.min(c)])
  
  return(join.dat)
}

## ----- Function to parallelize sf geospatial functions ----- ##

st_par <- function(sf_df, sf_func, n_cores,...){
  
  # Create a vector to split the data set up by.
  split_vector <- rep(1:n_cores, each = nrow(sf_df) / n_cores, length.out = nrow(sf_df))
  
  # Perform GIS analysis
  split_results <- split(sf_df, split_vector) %>%
    mclapply(function(x) sf_func(x, ...), mc.cores = n_cores)
  
  # Combine results back together. Method of combining depends on the output from the function.
  if ( class(split_results) == 'list' ){
    result <- do.call("c", split_results)
    names(result) <- NULL
  } else {
    result <- do.call("rbind", split_results)
  }
  
  # Return result
  return(result)
}

  



st_par2 <- function(sf_df, sf_func, n_cores,...){
  
  require(parallel)
  # Create a vector to split the data set up by.
  split_vector <- rep(1:n_cores, each = nrow(sf_df) / n_cores, length.out = nrow(sf_df))
  
  # Perform GIS analysis
  split_results <- split(sf_df, split_vector) %>%
    mclapply(function(x) sf_func(x, ...), mc.cores = n_cores)
  
  split_results2 <- lapply(split_results,function(x) apply(x,1,which.min))
  
  # Combine results back together. Method of combining depends on the output from the function.
  if ( class(split_results2) == 'list' ){
    result <- do.call("c", split_results2)
    names(result) <- NULL
  } else {
    "error: check output"
  }
  
  # Return result
  return(result)
}


## ----- Function to estimate streambed slope for GRDO sites ----- ##

Est_streambed_slope_DEM <- function(location,slope.dist,prj_dat){
  
  require(elevatr)
  require(tibble)
  require(dplyr)
  
  # 1. Obtain the coordinates
  coords_site <-  data.frame(x=location$Lon_WGS84, y=location$Lat_WGS84)
  
  # 2. Retrieve the coordinates of the points located within a defined range around the site:
  coords_corners <- destPoint(coords_site[1,], seq(0, 350, by=10), slope.dist)
  coords_site <- coords_site %>% add_row(x=coords_corners[,1], y=coords_corners[,2])
  
  # 3. Download raster DEM data for the site:
  dem_point <- elevatr::get_elev_raster(coords_site, prj = prj_dat,z = 12, src = "aws")
  
  # 4. Extract the elevation of the site as well as the points along the defined radius:
  coords_site <- tibble::add_column(coords_site, z=extract(dem_point, coords_site))
  
  # 5. Find the minimum elevation within the radius that corresponds to the downstream end of the reach:
  # 6. Calculate the slope between the site and the downstream end of the site reach
  slope <- abs(min(coords_site$z[-1])-coords_site$z[1])/slope.dist
  
  return(slope)
  
}
