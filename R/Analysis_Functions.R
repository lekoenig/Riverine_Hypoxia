## Functions called in hypoxia project - spatial R scripts
## Last updated 16 August 2019


## ----- Function to match COMID with flowline data by VPU hydroregion ----- ##

## 1. this function requires that the NHDV2 data is stored locally by VPU hydroregion

link_comid <-  function (points, CRS, path, vpu, maxDist){
  
  # set directory to appropriate vpu  
  directory <- grep(paste(vpu, "/NHDSnapshot/Hydrography",
    sep = ""), list.dirs(path, full.names = T), value = T)
   
  # read in flowlines and transform to albers equal area conic
  filename <- paste("./output/spatial/intermediate_NHDFlowline_data/","NHDFlowline_",vpu,sep="")
  if(file.exists(filename)){
    NHDFlowline <- readRDS(filename)
      if("Coastline" %in% NHDFlowline$FTYPE){
        NHDFlowline <- NHDFlowline[-which(NHDFlowline$FTYPE=="Coastline"),]
      }
    NHDFlowline <- sf::st_transform(NHDFlowline,crs=5070)} else {
      
    NHDFlowline <- sf::st_read(directory, layer = "NHDFlowline",quiet = T)
    saveRDS(NHDFlowline,filename)
    if("Coastline" %in% NHDFlowline$FTYPE){
      NHDFlowline <- NHDFlowline[-which(NHDFlowline$FTYPE=="Coastline"),]
    }
    NHDFlowline <- sf::st_transform(NHDFlowline,crs=5070)  
    }
    
  # drop Z and M dimensions from NHDFlowline$geometry (elevation + linear measurements)  
  NHDFlowline <- sf::st_zm(NHDFlowline, drop=T, what="ZM")  
  
  # replace geom labels with uppercase
  geom <- sf::st_geometry(NHDFlowline)
  sf::st_geometry(NHDFlowline) <- NULL
  names(NHDFlowline) <- toupper(names(NHDFlowline))
  NHDFlowline <- sf::st_sf(NHDFlowline, geom = geom)
   
  # read in vaa table  
  vaa.path <- grep(paste(vpu, "/NHDPlusAttributes", sep = ""),
                   list.dirs(path, full.names = T), value = T)    
  
  vaa.table <- grep("PlusFlowlineVAA.dbf", list.files(vaa.path[1],full.names=T),value=T)  
  vaa.table <- foreign::read.dbf(vaa.table)  
  names(vaa.table) <- toupper(names(vaa.table))
  
  # read in sample points
  points <- sf::st_as_sf(points,coords=c(2,3),crs=CRS)
  points <- sf::st_transform(points,crs=5070)
  
  # find flowline reaches within maxDist of sample point and choose closest flowline
  join.dat <- sf::st_join(points,NHDFlowline,sf::st_is_within_distance,dist=maxDist)
  dist <- st_distance(points,NHDFlowline[which(NHDFlowline$COMID %in% join.dat$COMID),])
  choose.dist <- which.min(dist)
  
  if(!length(dist))  {join.dat2 <- join.dat} else {
                      join.dat2 <- join.dat[choose.dist,]}
  join.dat2$near_dist_m <- if(!length(dist)) {NA} else {as.numeric(dist[choose.dist])}
  join.dat2 <- join.dat2[,c("db_ID","COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","near_dist_m")]
  
  # link vaa
  vaa.info <-vaa.table[which(vaa.table$COMID %in% join.dat2$COMID), c("TOTDASQKM","STREAMORDE")]
  if (length(vaa.info[, 1]) == 0) {
    vaa.info <- data.frame(TOTDASQKM = NA, STREAMORDE = NA)}
  
  join.dat3 <- cbind(as.data.frame(join.dat2),vaa.info)
  
  return(join.dat3)
}
 

##  2. this function subsets the national flowline data by vpu to speed up processing

link_comid2 <-  function (points, CRS, vpu, maxDist){
  
  flowline.sub <- flowline[which(flowline$VPUID==vpu), ]
  flowline.sub <- sf::st_transform(flowline.sub,crs=5070)
  
  # read in sample points
  points <- sf::st_as_sf(points,coords=c(2,3),crs=CRS)
  points <- sf::st_transform(points,crs=5070)
  
  # find flowline reaches within maxDist of sample point and choose closest flowline
  join.dat <- sf::st_join(points,flowline.sub,sf::st_is_within_distance,dist=maxDist)
  dist <- st_distance(points,flowline.sub[which(flowline.sub$COMID %in% join.dat$COMID),])
  choose.dist <- which.min(dist)
  
  if(!length(dist))  {join.dat2 <- join.dat} else {
    join.dat2 <- join.dat[choose.dist,]}
  join.dat2$near_dist_m <- if(!length(dist)) {NA} else {as.numeric(dist[choose.dist])}
  join.dat2 <- join.dat2[,c("db_ID","COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","VPUID","near_dist_m")]
  
  return(join.dat2)
}


## 3. this function uses nearest neighbor calculation to do the spatial join rather than calculating a minimum distance (as in link_comid and link_comid2)

link_comid3 <-  function (points, CRS, path, vpu, maxDist){
  
  flowline.sub <- flowline[which(flowline$VPUID==vpu), ]
  flowline.sub <- sf::st_transform(flowline.sub,crs=5070)
  
  # read in sample points
  points <- sf::st_as_sf(points,coords=c(2,3),crs=CRS)
  points <- sf::st_transform(points,crs=5070)
  
  # find the nearest neighbor flowline for each point, return a list
  
  nearest_flowline <- nngeo::st_nn(points,flowline.sub,returnDist=T,k=1,progress=FALSE)
  nearest_flowline_df <- data.frame(db_ID = points$db_ID,
                                dist=nearest_flowline$dist,
                                flowline_row=do.call(rbind,nearest_flowline$nn)) %>%
                         distinct(.keep_all = T)  

  nearest_flowline_df2 <- flowline.sub[nearest_flowline_df$flowline_row,c("COMID","GNIS_NAME",
                                                                          "FTYPE","FCODE",
                                                                          "StreamOrde","TotDASqKM")]
  nearest_flowline_df2$db_ID <- nearest_flowline_df$db_ID
  nearest_flowline_df2$nearest_dist_m <- nearest_flowline_df$dist
  st_geometry(nearest_flowline_df2) <- NULL
  return(nearest_flowline_df2)
  }
  

  
  




## ----- Function to match COMID with spatial points data using the nhdplusTools package ----- ##

comid_from_point = function(ID,lat, long, CRS) {
  pt = sf::st_point(c(long, lat))
  ptc = sf::st_sfc(pt, crs=CRS)
  COMID = nhdplusTools::discover_nhdplus_id(ptc)
  flag = length(COMID)
  if(! length(COMID)) COMID = NA
  if(length(COMID)>1) COMID = COMID[1]
  dat <- data.frame("db_ID" = ID, "NHDPlusV2_COMID" = COMID,"message" = flag)
  return(dat)
}



## ----- Modify link_to_flowlines function from hydrolinks package to match COMID with spatial points data ----- ##

link_to_flowlines2 <- function (lats, lons, ids, buffer, latlonCRS, dataset = c("nhdh", "nhdplusv2")){

  dataset = match.arg(dataset,c("nhdh","nhdplusv2"))
  dinfo = hydrolinks::dataset_info(dataset, "flowline")
  bbdf = NULL
  load(dinfo$bb_cache_path)
  sites = data.frame(lats, lons, ids)
  sites = sites[complete.cases(sites), ]
  pts = st_as_sf(sites, coords = c("lons", "lats"), crs = latlonCRS) 
  pts = st_transform(pts, st_crs(5070)) 
  bbdf = st_transform(bbdf, st_crs(5070))

  res = list()
  xmin = xmax = ymin = ymax = NULL
  for (i in 1:nrow(pts)) {
    vpu.assign <- bbdf[unlist(st_intersects(pts[i, ], bbdf)), "file", drop = TRUE]
    res[[i]] <- vpu.assign
  }
  
  # which vpu datasets to download
  to_check = as.data.frame(unique(unlist(res)), stringsAsFactors = FALSE)
  colnames(to_check)[1] = "file"
  
  if (nrow(to_check) == 0) {
    warning("hydrolinks::Supplied geopoints do not overlap ", 
            dataset, " dataset")
    ret = data.frame(MATCH_ID = rep(NA, 0))
    ret[, dinfo$id_column] = rep(NA, 0)
    return(ret)
  }
  
  match_res = list()
  for (i in 1:nrow(to_check)) {
    hydrolinks::check_dl_file(dinfo$file_index_path, to_check[i, "file"])
    shape = st_read(file.path(hydrolinks::cache_get_dir(), "unzip", to_check[i,"file"], dinfo$shapefile_name), 
                              stringsAsFactors = FALSE, quiet = TRUE)
    if("Coastline" %in% unique(shape$ftype)){
      shape = shape[-which(shape$ftype=="Coastline"),]
    }
    shape = st_transform(shape, 5070)

    shape_buffer = st_buffer(shape, buffer)
    matches = st_intersects(pts, shape_buffer)
   
    matches_multiple = which(lengths(matches) > 1)
    
    if (length(matches_multiple) > 0) {
      for (j in matches_multiple[seq_along(matches_multiple)]) {
        shape_rows = shape[matches[[j]],]
        distance = st_distance(pts[j, ], shape_rows)
        matches[[j]] <- matches[[j]][which.min(distance[1,])]
      }
    }
    
    shape_matched = shape[unlist(matches), ]
    shape_matched$MATCH_ID = pts[which(lengths(matches) > 0), ]$ids
    st_geometry(shape_matched) = NULL
    match_res[[i]] = data.frame(shape_matched, stringsAsFactors = FALSE)
  }
    
  unique_matches = unique(bind_rows(match_res))
  if (nrow(unique_matches) > 0) {
    return(unique_matches[!is.na(unique_matches[, dinfo$id_column]), ])
  } else {
    return(unique_matches)
  }
}
  