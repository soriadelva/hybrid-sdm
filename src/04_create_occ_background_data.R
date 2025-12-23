# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
# Script for creating datasets of occurrences and background points for each run
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------


#--------------------------------------------------
#--------------- Start clean  ----------------
#--------------------------------------------------
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)


#--------------------------------------------------
#--------------- Source scripts  ----------------
#--------------------------------------------------
source(file.path("src", "03_sdm_configurations.R"))
source(file.path("src", "helpers", "Main.Functions.R"))
  
  
#--------------------------------------------------
#-------Define and Create folders  ----------------
#--------------------------------------------------
occs_all_folder<-file.path("data", "Input", "Occurrences")
data.folder.rasters  <- file.path("data", "Input", "Rasters", "Layers_present")
if(!dir.exists(occs_all_folder)) dir.create(occs_all_folder, recursive=TRUE)


#--------------------------------------------------
#--------------- Load Files ----------------
#--------------------------------------------------
occs_file <- file.path(occs_all_folder,"Dictyota_full_occurrences.csv")

if(!file.exists(occs_file)){
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.17865936", 
                       path= occs_all_folder,
                       files="Dictyota_full_occurrences.csv",
                       quiet=FALSE)
}

all.occurrences <- read.csv2(occs_file)



# -------------------------------------------------------------------------
#----------------------Import Occurrence records---------------------------
#--------------------------------------------------------------------------
occurrence.records <- all.occurrences %>%
  dplyr::rename(Lon = Longitude,
                Lat = Latitude) %>%
  dplyr::filter(!is.na(Lon) & !is.na(Lat),
                !duplicated(dplyr::select(., Lon, Lat),fromLast = TRUE),
                Lon <= 38 & Lon >= -25 & Lat <= 73 &  Lat >= 27)%>% #Removes 12 records
  dplyr::select(Lon, Lat)


# -----------------------------------------------------------
#--------------------Import Rasters--------------------------
#------------------------------------------------------------
rasters.files <- list.files(data.folder.rasters, pattern="tif", full.names = TRUE)
rasters <- terra::rast(rasters.files[as.vector(sapply(rasters.names.s,function(x) { which( grepl(x,rasters.files)) } ))])

data.frame(Rasters=names(rasters),Names=rasters.names.s,Resp=variable.monotonic.response)
names(rasters) <- rasters.names.s # [!]


#--------------------------------------------------------------
# -----------------Process rasters ----------------------------
#--------------------------------------------------------------
#----------------------------------
#- Crop to ext. occurrences + 12° -
#----------------------------------
final.extent.rasters <- terra::ext( min(occurrence.records[,"Lon"] , na.rm=T) - region.buffer[1],
                                       max(occurrence.records[,"Lon"] , na.rm=T) + region.buffer[2],
                                       min(occurrence.records[,"Lat"] , na.rm=T) - region.buffer[3],
                                       max(occurrence.records[,"Lat"] , na.rm=T) + region.buffer[4] )
    
  rasters <- terra::crop(rasters, final.extent.rasters) # -29.2 42.2 15.8 73.6 

  
  #---------------------------------------------
  #- Mask on ecoregions with occs + neighbours -
  #---------------------------------------------
  #Load marine ecoregions and group by PROV_CODE
  ecoregions <- sf::st_read(file.path("data","Input", "Dependencies","Shapefiles","marine_ecoregions.shp"), quiet=T)%>%
    group_by(PROV_CODE) %>%
    summarize(geometry = st_union(geometry), .groups = "drop")
  
  #Convert occurrences to sf dataframe
  occ_sf <- st_as_sf(occurrence.records,
                     coords = c("Lon", "Lat"),
                     crs = st_crs(ecoregions))
  
  # Get unique ecoregions with occurrences
  prov_with_occ  <- st_join(occ_sf, ecoregions)%>%
    pull(PROV_CODE)%>%
    unique()
  
    # compute adjacency (neighbors) using st_touches -> list of integer vectors
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
    
    #---------------------------------------------
    #- Mask based on depth -
    #---------------------------------------------
    # Load bathymetry raster and crop to extent
    bathymetry <- terra::rast(bathymetry.file)%>%
    terra::crop(final.extent.rasters)
    
    # Create depth mask based on depth distribution of species
    depth_mask_relocation <- (bathymetry >= -max.depth) & (bathymetry <= min.depth)
    
    #Mask with depth mask
    rasters.relocation <- terra::mask(rasters, depth_mask_relocation, maskvalues = FALSE)

    # Create depth mask based on depth distribution of species + depth buffer
    depth_mask <- (bathymetry >= -(max.depth + depth.buffer)) & (bathymetry <= min.depth)
    
    #Mask with depth mask
    rasters <- terra::mask(rasters, depth_mask, maskvalues = FALSE)


# ------------------------------------------------------------
#-----Test for correlations / collinearity between rasters-----
#-------------------------------------------------------------
pairs(rasters)
#In the study region, LTmax temperature and LTminTemperature are correlated with 0.7 --> at the edge but okay


# --------------------------------------------------------------------
# -----Relocate records to closest cell (those falling on land)-------
#---------------------------------------------------------------------
occurrence.records <- relocate.coordinates.na(occurrence.records,rasters.relocation,maximum.distance=relocate.occ.distance,use.species.depth=relocate.occ.species.depth,min.depth=min.depth,max.depth=max.depth)
#Relocate 57 points, lose 2

#-----------------------------------------------------------------------
#-------------Remove deep records outside known depth range-------------
#-----------------------------------------------------------------------
#Remove deep records
presences.depths <- abs( terra::extract(bathymetry, occurrence.records, ID=F) )%>%
      pull()
to.remove <- which( presences.depths > max.depth | is.na(presences.depths) )
if( sum(to.remove) != 0 ) { occurrence.records <- occurrence.records[ -to.remove,] }

# Visualize with histogram
hist(presences.depths,breaks=50,xlab="Depth (m)", freq = FALSE)
max(presences.depths) #49m


#------------------------------------------------------------------
#---------Remove points within a distance of 9.2 km -------
#------------------------------------------------------------------
# Use 9.2 (resolution bio oracle at equator)
distance.uncorr <- 9.2
occurrence.records <- spatial.autocorrelation.thinning(occurrence.records,distance.uncorr)
nrow(occurrence.records)

# Plot processed records
plot.google(occurrence.records,3,"black")

#Export processed records
write.csv2(occurrence.records,file.path(occs_all_folder,paste0("Processed_occurrences.csv")), row.names=F)


# ------------------------------------------------------------------------------------
#---Load polygon for exclusion of region for pseudoabsence selection-----
# ------------------------------------------------------------------------------------
polygon.pa <- st_read(file.path(data.folder, "Dependencies", "Spatial", paste0(project.name,".polygon.pa.shp")), quiet=T)


# --------------------------------------------------
#------Prepare dataset with all occurrences-----------
#---------------------------------------------------
#Generate dataset with all occurrences
all.occurrences  <- all.occurrences %>%
  dplyr::rename(Lon = Longitude,
                Lat = Latitude) %>%
  dplyr::select(Lon, Lat)


# --------------------------------------------------
#-------------------Specifiy specific run-----------
#---------------------------------------------------
Background.sizes<-c(640, 1000, 2000, 10000, 20000, 40000)


#---------------------------------------------------
#--------Start loop and set different seed----------
#---------------------------------------------------
i<-1
for(Background.size in Background.sizes){
  seed <- i
  set.seed(seed)
  cat("Using seed for background size", Background.size ,":", seed, "\n")
  
  
  #--------------------------------------------------
  #--------------- Create folders  ----------------
  #--------------------------------------------------
  background_folder<- file.path("data", "Input", "Background_points", paste0("Background_",Background.size))
  if(!dir.exists(background_folder)) dir.create(background_folder, recursive=TRUE)


#----------------------------------------------------------------
#----------------------Create absences and cv data---------------
#----------------------------------------------------------------
  #Set distance.uncorr to 5km for background selection or we won't be able to select 40000 points
  distance.uncorr<-5
  bias.surface.path <-  file.path("data", "Input", "Biasgrid", 
                                  "Seaweed_HybridSDM_log_counts_biasgrid.tif")
suppressMessages(
  background.points <- sample.background(occurrence.records= all.occurrences,
                                         rasters,
                                         polygon.pa,
                                         polygon.pa.type = "Exclude",
                                         pa.dist.min = 30, #km
                                         pa.type = "random",
                                         bias.surface.kernel = TRUE,
                                         n.background = Background.size,
                                         pa.environment.stratification = TRUE,
                                         plot.figure=FALSE,
                                         bias.surface.path = bias.surface.path)
  )


#Sample 50000 background points
# background.points.large <- sample.background.large(occurrence.records = all.occurrences, 
#                                                      raster_template = rasters[[1]],
#                                                      n_background = 50000,
#                                                      seed = seed)



# Plot background.points
plot.google(background.points,3,"black")
# plot.google(background.points.large,3,"black")

#Export absences
write.csv2(background.points,file.path(background_folder,paste0("Background_points_",Background.size,".csv")), row.names=F)
# write.csv2(background.points.large,file.path(background_folder,paste0("Background_points_50000_",Run,".csv")), row.names=F)

cat("Created processed occurrence and background datasets for background size:\n", Background.size, "\n")
assign(paste("absences_",Background.size),background.points)
# assign(paste("absenceslarge_",Run),background.points.large)
i<-i+1
}
