# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#                  Modelling the distribution of Dictyota dichotoma
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------


#--------------------------------------------------
#--------------- Start clean  ----------------
#--------------------------------------------------
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)


#--------------------------------------------------
#----TO DO: indicate background sizes for loop 
#--------------------------------------------------
#Determine background sizes
#Background_sizes<-c(640, 1000, 2000, 10000, 20000, 40000)

#If it fails run this script size per size!
Background_sizes<-c(40000) #640 1000, 2000, 10000, 20000, 40000
Model_types <- c("Correlative", "Growth", "Germination", "FvFm")


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
cv_folder<- file.path("data", "Input", "Cross_validation_folds")
if(!dir.exists(cv_folder)) dir.create(cv_folder, recursive=TRUE)


# ---------------------------------------------------------------
#--------- Import processed occurrence records -------------
#----------------------------------------------------------------
#Import presences
occurrence.records <- read.csv2(file.path(occs_all_folder,"Processed_occurrences.csv"))
plot.google(occurrence.records,3,"black")


# ---------------------------------------------------------------
#--------- Export cross validation folds -------------
#----------------------------------------------------------------
#Import background data from largest background size (40000)
background_folder <- file.path("data", "Input", "Background_points", "Background_40000")
background <- read.csv2(file.path(background_folder,paste0("Background_points_40000.csv")))

#Generate 7 (cv.k) cross validation blocks
cross.validation.data <- data.partitioning(occurrence.records,background,cv.type, cv.k)

#Export
write.csv2(cross.validation.data, file.path(cv_folder, "crossvalidation_blocks.csv"), row.names = F)
rm(background, background_folder)


#----------------------------------------------------------------
#--------Prepare cross validation combinations dataframe---------
#----------------------------------------------------------------
#Indicate number of cross validation blocks
cv.k.span <- 1:cv.k

#Remove outer two cross validation blocks as test dataset
cv.k.span <- cv.k.span[c(-1,-length(cv.k.span))] 

#make a df with all possible combinations of cv blocks, learning complex values, and tree depth
comb = expand.grid(cv.k = cv.k.span, learning.complex=brt.learning.complex.span,tree.depth=brt.max.tree.depth ) 


# ---------------------------------------------------------------
#-------------------- Start model type loop ---------------------
#----------------------------------------------------------------
for(Model_type in Model_types){
  
  # -----------------------------------------------------------
  #---Define raster names and monotonic responses flexibly-----
  #------------------------------------------------------------
  rasters.names.s <- switch(Model_type,
                            "Correlative"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                             "Present.Benthic.Min.Depth.Salinity.Lt.min",
                                             "Present.Benthic.Min.Depth.Temperature.Lt.max",
                                             "Present.Benthic.Min.Depth.Temperature.Lt.min"),
                            "Germination"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                             "Present.Benthic.Min.Depth.Salinity.Lt.min",
                                             "Germination_LTmax_present",
                                             "Germination_LTmin_present"),
                            "Growth"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                        "Present.Benthic.Min.Depth.Salinity.Lt.min",
                                        "Growth_LTmin_present",
                                        "Growth_LTmax_present"),
                            "FvFm"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                      "Present.Benthic.Min.Depth.Salinity.Lt.min",
                                      "FvFm_LTmax_present",
                                      "FvFm_LTmin_present")
  )
  
  if (Model_type == "Correlative") {
    variable.monotonic.response <- c(+1, +1, -1, +1)
  } else {
    variable.monotonic.response <- c(+1, +1, +1, +1)
  }
  
  
  # -----------------------------------------------------------
  #--------------------Import Rasters--------------------------
  #------------------------------------------------------------
  data.folder.rasters  <- file.path("data", "Input", "Rasters", "Layers_present")
  rasters.files <- list.files(data.folder.rasters, pattern="tif", full.names = TRUE)
  
  if(Model_type!="Correlative"){
    hybrid.rasters <- list.files(file.path("data", "Input", "Hybrid_rasters", "Layers_present"), pattern="tif", full.names = TRUE)
    rasters.files<-c(rasters.files, hybrid.rasters)
  }
  
  rasters <- terra::rast(rasters.files[as.vector(sapply(rasters.names.s,function(x) { which( grepl(x,rasters.files)) } ))])
  data.frame(Rasters=names(rasters),Names=rasters.names.s,Resp=variable.monotonic.response)
  names(rasters) <- rasters.names.s # [!]
  
  
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
  
  
  # ------------------------------------------------------------
  #-----Test for correlations / collinearity between rasters-----
  #-------------------------------------------------------------
  #In the study region, LTmax temperature and LTminTemperature are correlated with 0.7 --> at the edge but okay
  pairs(rasters)
  
  
  #----------------------------------------------------------------
  #------Write rasters to disk for niche modelling ----------------
  #----------------------------------------------------------------
  
  # # Clip rasters for performance
  # extent.t <- terra::ext( min(c(background[,1],occurrence.records[,1])) - 1 , max(c(background[,1],occurrence.records[,1])) + 1 , min(c(background[,2],occurrence.records[,2])) - 1 , max(c(background[,2],occurrence.records[,2])) + 1 )
  # rasters1 <- terra::crop(rasters,extent.t)
  
  #Export rasters so they can be read inside parallel loop
  interim.folder<-file.path("data", "Interim", Model_type, paste0(Model_type, "rasters"))
  if(!dir.exists(interim.folder)) dir.create(interim.folder, recursive=TRUE)
  outfile <- file.path(interim.folder, "predictors_stack.tif")
  terra::writeRaster(rasters, outfile, overwrite = TRUE)
  
  
  # ------------------------------------------------------------
  #-----Start loop for background size
  #-------------------------------------------------------------
  
  for (Background_size in Background_sizes){
    
    cat("Performing cross validation for",Model_type , "model using", 
        Background_size, "background points")
    start.time<-Sys.time()
    
    
    #----------------------------------------------------------------
    #---------------- Define pseudoabsence folders ------------------
    #----------------------------------------------------------------
    background_folder <- file.path("data", "Input", "Background_points", paste0("Background_",Background_size))
    
    
    #----------------------------------------------------------------
    #-------------------- Import background points ------------------
    #----------------------------------------------------------------
    #Import for each run
    background <- read.csv2(file.path(background_folder, paste0("Background_points_",Background_size,".csv")))
    
    
    #----------------------------------------------------------------
    #--Create overview of n presences and n background for each fold-
    #----------------------------------------------------------------
    CV.overview.file<-file.path(cv_folder, paste0("CV_folds_overview",Background_size,".csv"))
                                
    if(!file.exists(CV.overview.file)){
      n.presences<-data.frame()
      presences<-occurrence.records
      background.points<-background
      
      for (fold in 1:length(cv.k.span)){
        cv <- comb[fold,1]
        presences.train <- presences[ presences$Lat < cross.validation.data[cv,1] | presences$Lat > cross.validation.data[cv,2] , ]
        background.points.train <- background.points[ background.points$Lat < cross.validation.data[cv,1] | background.points$Lat > cross.validation.data[cv,2] , ]
        presences.test <- presences[ presences$Lat >= cross.validation.data[cv,1] & presences$Lat <= cross.validation.data[cv,2] , ]
        background.points.test <- background.points[ background.points$Lat >= cross.validation.data[cv,1] & background.points$Lat <= cross.validation.data[cv,2] , ]
        
        #Remove NAs  
        presences.train <- presences.train[complete.cases(presences.train), ]
        background.points.train <- background.points.train[complete.cases(background.points.train), ]
        presences.test <- presences.test[complete.cases(presences.test), ]
        background.points.test <- background.points.test[complete.cases(background.points.test), ]
        
        results<-data.frame(test_cv = cv,
                            n_train_presences = nrow(presences.train),
                            n.test.presences = nrow(presences.test),
                            n_train_background = nrow(background.points.train),
                            n.test.background = nrow(background.points.test))
        if(nrow(n.presences)==0){
          n.presences<-results
        }else{
          
          n.presences<-bind_rows(n.presences, results)
        }
        
      }
      
      #Export overview
      write.csv2(n.presences, CV.overview.file, row.names = F)
      rm(n.presences, background.points, results, presences.train, presences.test, background.points.train, background.points.test, cv, fold)
    }
    
    
    #----------------------------------------------------------------
    #--------Prepare environment for parallel processing --------
    #----------------------------------------------------------------
    plan(multisession, workers = 25)
    options(future.globals.maxSize = 4.5 * 1024^3)  # increase if needed
    
    
    #----------------------------------------------------------------
    #-------------Perform cross validated niche modelling -----------
    #----------------------------------------------------------------
    if( exists("model") ) { rm(model) }
    cv.accuracy <- future_lapply(1:nrow(comb), function(i){
      source("src/helpers/niche_modelling_cross_validation_helpers.R")
      run_brt_crossvalidation(i,  
                              comb = comb, 
                              cross.validation.data = cross.validation.data, 
                              presences = occurrence.records, 
                              background.points = background,
                              raster.file = outfile,
                              variable.monotonic.response.i = variable.monotonic.response)
     
    },  future.seed = TRUE
    )
    
    
    #-----------------------------------------------
    #------------Settings back to normal -----------
    #-----------------------------------------------
    plan(strategy = "sequential")   #Close parallel processing
    gc() # Clean up memory after processing
    options(future.globals.maxSize = 500 * 1024^2)  # Reset to 500 MB
    
    
    #---------------------------------
    #---------Combine results---------
    #---------------------------------
    # (i.e. each combo of lr, tc, and cross validation block)
    cv.accuracy <- do.call(rbind,cv.accuracy)
    
    
    #---------------------------------
    #---------Export results ---------
    #---------------------------------
    results_folder<-file.path("results", "cross_validation", Model_type, Background_size)
    dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
    write.csv2(cv.accuracy,file.path(results_folder,paste0("cv.accuracy",Background_size,".csv")), row.names=F)
    
    
    #---------------------------------
    #---------Clean up ---------
    #---------------------------------
    rm(background, variable.monotonic.response, results_folder)
    terra::tmpFiles(remove = TRUE)
  
    
    #---------------------------------
    #---------Mark time ---------
    #---------------------------------
    elapsed <- round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), 2)
    message(sprintf("[Model: %s | Background size: %d] Completed cross-validation in %0.2f minutes", 
                    Model_type, Background_size, elapsed))
    
  }
}


#
# #------------------------------------------------------------------
# #--------Visualize effect of background size on model metrics------
# #------------------------------------------------------------------
best_metrics<-data.frame()
Background_sizes <- c(640,1000, 2000,  10000, 20000, 40000)
Model_types <- c("Correlative", "Growth", "Germination", "FvFm")
for (Model_type in Model_types){
  for(Background_size in Background_sizes){

    # Read cv results
    cv_folder<-file.path("results", "cross_validation", Model_type, Background_size)
    cv_metrics<-read.csv2(file.path(cv_folder, paste0("cv.accuracy",Background_size, ".csv")))

    #Process cv results: calculate boyce/auc mean and se per combo of tree complexity and learning rate
    cv_metrics_summary <- cv_metrics %>%
      dplyr::group_by(tree.c, l.rate) %>%
      dplyr::summarise(
        mean_boyce = mean(boyce, na.rm = TRUE),
        se_boyce   = sd(boyce, na.rm = TRUE) / sqrt(dplyr::n()),
        mean_auc   = mean(auc, na.rm = TRUE),
        se_auc     = sd(auc, na.rm = TRUE) / sqrt(dplyr::n()),
        .groups = "drop"
      )

    #Select results for best parameter combination (max boyce)
    cv_best_metrics <- cv_metrics_summary%>%
      arrange(desc(mean_boyce), desc(mean_auc)) %>%  # sort by Boyce, then AUC
      slice(1) %>%# take top row
      mutate(Model = Model_type,
             Background_size = Background_size)

    #Export files
    write.csv2(cv_metrics_summary,
               file.path(cv_folder,paste0("cv_metrics_summary",Background_size,".csv")),
               row.names=F)
    write.csv2(cv_best_metrics,
               file.path(cv_folder,paste0("cv_best_metrics",Background_size,".csv")),
               row.names=F)


    #Create dataframe of best metrics per model and background size
    if(nrow(best_metrics)==0){
      best_metrics<-cv_best_metrics
    }else{
      best_metrics <- dplyr::bind_rows(best_metrics, cv_best_metrics)
    }
  }
}

best_metrics
write.csv2(best_metrics, file.path("results", "cross_validation", "best_metrics.csv"))

#------------------------------------
#---Visualize best metrics-----------
#------------------------------------
ggplot(best_metrics, aes(x = Background_size, y = mean_boyce)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = mean_boyce - se_boyce, ymax = mean_boyce + se_boyce), width = 0.2) +
  geom_line(color = "steelblue") +
  scale_x_continuous(breaks = unique(best_metrics$Background_size)) +
  labs(
    x = "Background sample size",
    y = "Boyce index"
  ) +
  facet_wrap(~Model,
             labeller = labeller(
               Model = c(Correlative = "Correlative",
                         Growth = "Growth",
                         Germination = "Germination",
                         FvFm = "F[v]/F[m]"),
               .default = label_parsed))+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=9)) #Export width 900 690

ggplot(best_metrics, aes(x = Background_size, y = mean_auc)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc), width = 0.2) +
  geom_line(color = "steelblue") +
  scale_x_continuous(breaks = unique(best_metrics$Background_size)) +
  labs(
    title = paste("Mean AUC for each model"),
    x = "Background sample size",
    y = "Mean AUC"
  ) +
  facet_wrap(~Model)+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=9))



