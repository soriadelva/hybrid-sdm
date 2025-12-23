#--------------------------------------------------
#--------------- Source scripts  ----------------
#--------------------------------------------------
source(file.path("src", "03_sdm_configurations.R"))
source(file.path("src", "helpers", "Main.Functions.R"))
source(file.path("src", "helpers", "niche_modelling_cross_validation_helpers.R"))


# --------------------------------------------------
#------------- Define files and folders-------------
#---------------------------------------------------
occs_all_folder<-file.path("data", "Input", "Occurrences")
meow<-file.path("data","Input", "Dependencies","Shapefiles","marine_ecoregions.shp")


# ---------------------------------------------------------------
#--------- Import processed occurrence records -------------
#----------------------------------------------------------------
#Import presences
occurrence.records <- read.csv2(file.path(occs_all_folder,"Processed_occurrences.csv"))
plot.google(occurrence.records,3,"black")


predict.times <- c("2100_RCP26", "2100_RCP45", "2100_RCP85")
Model_types <- c("Correlative", "Growth", "Germination", "FvFm")

for(Model_type in Model_types){
  for (predict.time  in predict.times){
    
    if(Model_type =="Correlative"){
      raster_folder<-file.path("data", "Input", "Rasters", paste0("Layers_future",predict.time))
    }else{
      raster_folder<-file.path("data", "Input", "Hybrid_rasters", paste0("Layers_future",predict.time))
    }
    
    #-------------------------------------
    #------Prepare rasters----------------
    #-------------------------------------
    year<-as.numeric(sub("_.*", "", predict.time))
    rcp<-sub(".*_", "", predict.time)
    rasters.names.future.old <- switch(Model_type,
                                       "Correlative"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                                        paste0(year,"AOGCM.",rcp,".Benthic.Min.Depth.Salinity.Lt.min.tif.BOv2_1"),
                                                        paste0(year,"AOGCM.",rcp,".Benthic.Min.Depth.Temperature.Lt.max.tif.BOv2_1"),
                                                        paste0(year,"AOGCM.",rcp,".Benthic.Min.Depth.Temperature.Lt.min.tif.BOv2_1")),
                                       "Germination"= c(paste0("Germination_LTmax_", predict.time),
                                                        paste0("Germination_LTmin_", predict.time)),
                                       "Growth"= c(paste0("Growth_LTmin_", predict.time),
                                                   paste0("Growth_LTmax_", predict.time)),
                                       "FvFm"= c(paste0("FvFm_LTmax_", predict.time),
                                                 paste0("FvFm_LTmin_", predict.time)))
    
    
    rasters.names.future <- switch(Model_type,
                                   "Correlative"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                                    "Present.Benthic.Min.Depth.Salinity.Lt.min",
                                                    "Present.Benthic.Min.Depth.Temperature.Lt.max",
                                                    "Present.Benthic.Min.Depth.Temperature.Lt.min"),
                                   "Germination"= c("Germination_LTmax_present",
                                                    "Germination_LTmin_present"),
                                   "Growth"= c("Growth_LTmin_present",
                                               "Growth_LTmax_present"),
                                   "FvFm"= c("FvFm_LTmax_present",
                                             "FvFm_LTmin_present"))
    
    
    future.rasters.files <- list.files(raster_folder, pattern="tif", full.names = TRUE)
    rasters <- terra::rast(future.rasters.files[as.vector(sapply(rasters.names.future.old,function(x) { which( grepl(x,future.rasters.files)) } ))])
    
    #Give the layers the same names as the present layer names: needed for modelling
    names(rasters)<-rasters.names.future
    
    #Overview
    overview<-data.frame(New_names=names(rasters),Old_Names= rasters.names.future.old)
    overview
    
    #--------------------------------------------------------------
    # -----------------Process rasters ----------------------------
    #--------------------------------------------------------------
    #- Crop to ext. occurrences + region.buffer (12°) -
    final.extent.rasters <- terra::ext( min(occurrence.records[,"Lon"] , na.rm=T) - region.buffer[1],
                                        max(occurrence.records[,"Lon"] , na.rm=T) + region.buffer[2],
                                        min(occurrence.records[,"Lat"] , na.rm=T) - region.buffer[3],
                                        max(occurrence.records[,"Lat"] , na.rm=T) + region.buffer[4] )
    
    rasters <- terra::crop(rasters, final.extent.rasters) # -29.2 42.2 15.8 73.6 
    
    
    #- Mask on ecoregions with occs + neighbours -
    ecoregions <- sf::st_read(meow, quiet=T) %>%
      group_by(PROV_CODE) %>%
      summarize(geometry = st_union(geometry), .groups = "drop")
    
    #Convert occurrences to sf dataframe
    occ_sf <- st_as_sf(occurrence.records,
                       coords = c("Lon", "Lat"),
                       crs = st_crs(ecoregions))
    
    # Get unique ecoregions with occurrences
    prov_with_occ  <- st_join(occ_sf, ecoregions) %>%
      pull(PROV_CODE) %>%
      unique()
    
    # compute neighbors using st_touches
    neighbors_list <- st_touches(ecoregions)
    
    #get neighbors of polygons that have occurrences
    neighbor_indices <- unique(unlist(neighbors_list[prov_with_occ]))
    
    # combine original polys and their neighbors
    final_ecoregions_idx <- sort(unique(c(prov_with_occ, neighbor_indices)))
    
    #Filter ecoregions
    ecoregions <- ecoregions %>%
      filter(PROV_CODE %in% final_ecoregions_idx)
    
    #Mask based on ecoregions
    rasters <- terra::mask(rasters, ecoregions)
    
    # Load bathymetry raster and crop to extent
    bathymetry <- terra::rast(bathymetry.file)%>%
      terra::crop(final.extent.rasters)
    
    # Create depth mask based on depth distribution of species + depth buffer
    depth_mask <- (bathymetry >= -(max.depth + depth.buffer)) & (bathymetry <= min.depth)
    
    #Mask with depth mask
    rasters <- terra::mask(rasters, depth_mask, maskvalues = FALSE)
    
    
    #-------------------------------------------------------
    #---------------Export rasters--------------------------
    #-------------------------------------------------------
    #Create future file names
    rasters.files.names.future <- switch(Model_type,
                                         "Correlative"= c("light.at.bottom.Benthic.Min.Var.Mean.tif",
                                                          paste0(predict.time,".Benthic.Min.Depth.Salinity.Lt.min.tif"),
                                                          paste0(predict.time,".Benthic.Min.Depth.Temperature.Lt.max.tif"),
                                                          paste0(predict.time,".Benthic.Min.Depth.Temperature.Lt.min.tif")),
                                         "Germination"= c(paste0("Germination_LTmax_", predict.time,".tif"),
                                                          paste0("Germination_LTmin_", predict.time,".tif")),
                                         "Growth"= c(paste0("Growth_LTmin_", predict.time,".tif"),
                                                     paste0("Growth_LTmax_", predict.time,".tif")),
                                         "FvFm"= c(paste0("FvFm_LTmax_", predict.time,".tif"),
                                                   paste0("FvFm_LTmin_", predict.time,".tif")))
    
    for(rasterlayer in 1:nlyr(rasters)){
      single_raster<-rasters[[rasterlayer]]
      writeRaster(single_raster, file.path(raster_folder,  rasters.files.names.future[[rasterlayer]]), overwrite=T)
      
    }
    write.csv2(overview, file.path(raster_folder, paste0(Model_type,"_renaming_scheme.csv")))
    
  }
}