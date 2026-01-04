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
#------- Spefify background size and models -------
#--------------------------------------------------
Background_sizes<-1000
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
#--------- Import background records -------------
#----------------------------------------------------------------
#Import background data from 
background_folder <- file.path("data", "Input", "Background_points", paste0("Background_",Background_sizes))
background <- read.csv2(file.path(background_folder,paste0("Background_points_",Background_sizes,".csv")))
plot.google(background, 3, "black")


# --------------------------------------------------------------------------------
#-- Extract LTmin_temperature and LTmax temperature at training data locations ---
#---------------------------------------------------------------------------------
data.folder.rasters  <- file.path("data", "Input", "Rasters", "Layers_present")
LTmintemp_r<-terra::rast(file.path(data.folder.rasters, "Present.Benthic.Min.Depth.Temperature.Lt.min.tif"))
LTmaxtemp_r<-terra::rast(file.path(data.folder.rasters, "Present.Benthic.Min.Depth.Temperature.Lt.max.tif"))

LTmintemp_occ<-terra::extract(LTmintemp_r, rbind(occurrence.records, background), ID=F)
LTmaxtemp_occ<-terra::extract(LTmaxtemp_r, rbind(occurrence.records, background), ID=F)

LTmintemp_seq <- data.frame(Temperature=seq(min(LTmintemp_occ), max(LTmintemp_occ), length.out=100))
LTmaxtemp_seq <- data.frame(Temperature=seq(min(LTmaxtemp_occ), max(LTmaxtemp_occ), length.out=100))


# --------------------------------------------------------------------------------
#---------------------Fit physiological models on the data -----------------------
#---------------------------------------------------------------------------------
input_folder<-file.path("data", "Input", "Physiology")
Growth <- read.csv2(file.path(input_folder,"Dictyota_growth.csv"))%>%
  dplyr::rename(Value = Growth_rate,
                Temperature = Temperature_C) %>%
  dplyr::mutate(Trait = "Growth")%>%
  dplyr::select(Value, Trait, Temperature)

FvFm <- read.csv2(file.path(input_folder,"Dictyota_FvFm.csv"))%>%
  dplyr::rename(Value = FvFm,
                Temperature = Temperature_C) %>%
  dplyr::mutate(Trait = "FvFm")%>%
  dplyr::select(Value, Trait, Temperature)

Germination <- read.csv2(file.path(input_folder,"Dictyota_germination.csv"))%>%
  dplyr::rename(Value = Germination_rate_percent,
                Temperature = Temperature_C) %>%
  dplyr::mutate(Trait = "Germination")%>%
  dplyr::select(Value, Trait, Temperature)

Alldata<-bind_rows(FvFm, Growth, Germination)%>%
  mutate(TemperatureK = Temperature + 273.15)

#------------------------------------------------------------
#--------------- Define best TPC models----------------------
#------------------------------------------------------------
#For growth
formula5 <- as.formula( Value ~ ifelse(Temperature <= Topt,
                                       Gmax * exp(-((Temperature - Topt) / (2 * a))^2),
                                       Gmax - Gmax * (((Temperature - Topt) / (Topt - Tmax))^2)))#Deutsch

#For FvFm and Germination
formula7 <- as.formula(Value~ a*exp(-0.5*(abs(Temperature-tref)/b)^c)) # modification of a gaussian function


#------------------------------------------------------------
#--------------- Prepare trait datasets----------------------
#------------------------------------------------------------
#Initialize trait dataframes
growth_data<-filter(Alldata, Trait=="Growth")
pam_data<-filter(Alldata, Trait=="FvFm")
germination_data<-filter(Alldata, Trait=="Germination")


#--------------------------------------------------
# Create fits using the best model
#--------------------------------------------------
growth_fit<- nlsLM(formula5, 
                   data=growth_data, 
                   start = list(Gmax = 8, Tmax = 32, Topt = 22, a=4), 
                   control = list(maxiter = 100))

pam_fit<- nlsLM(formula7, 
                data=pam_data, 
                start = list(a = 10, b = 10, tref=20, c=10))

germination_fit<- nlsLM(formula7, 
                        data=germination_data, 
                        start = list(a = 100, b = 10, tref=20, c=10))


# --------------------------------------------------------------------------------
#--------Create hybrid temperature vectors for partial effect plots---------------
#---------------------------------------------------------------------------------

LTmintemp_seq$Growth<-predict(growth_fit, LTmintemp_seq)
LTmintemp_seq$FvFm <- predict(pam_fit, LTmintemp_seq)
LTmintemp_seq$Germination <- predict(germination_fit, LTmintemp_seq)

LTmaxtemp_seq$Growth<-predict(growth_fit, LTmaxtemp_seq)
LTmaxtemp_seq$FvFm <- predict(pam_fit, LTmaxtemp_seq)
LTmaxtemp_seq$Germination <- predict(germination_fit, LTmaxtemp_seq)

LTmintemp_seq<-LTmintemp_seq%>%
  mutate(Growth = ifelse(Growth < 0, 0, Growth))
LTmaxtemp_seq<-LTmaxtemp_seq%>%
  mutate(Growth = ifelse(Growth < 0, 0, Growth))


# --------------------------------------------------------------------------------
#----------Create other vectors for partial effect plots---------------------
#---------------------------------------------------------------------------------
Mean_light_r<-terra::rast(file.path(data.folder.rasters, "light.at.bottom.Benthic.Min.Var.Mean.tif"))
Min_salinity_r<-terra::rast(file.path(data.folder.rasters, "Present.Benthic.Min.Depth.Salinity.Lt.min.tif"))

Meanlight_occ<-terra::extract(Mean_light_r, rbind(occurrence.records, background), ID=F)
LTminsalinity_occ<-terra::extract(Min_salinity_r, rbind(occurrence.records, background), ID=F)

Meanlight_seq <- data.frame(Light=seq(min(Meanlight_occ), max(Meanlight_occ), length.out=100))
LTminsalinity_seq <- data.frame(Salinity=seq(min(LTminsalinity_occ), max(LTminsalinity_occ), length.out=100))


# ---------------------------------------------------------------
#--------- Import MedSea shape -------------
#----------------------------------------------------------------
MedSea<- st_read(file.path("Data", "Input", "MedSea", "MedSea.shp"))
MedSea<-terra::vect(MedSea)


# ---------------------------------------------------------------
#-------------------- Start model type loop ---------------------
#----------------------------------------------------------------

for(Model_type in Model_types){
  
  # ---------------------------------------------------------------
  #-------------Initiate dataframes ---------------------
  #----------------------------------------------------------------
  varimp_bootstrap_df<-data.frame()
  response_bootstrap_df<-data.frame()
  area_depth_bootstrap_df <- data.frame()
  logs_bootstrap_df<-data.frame()
  area_depth_future_bootstrap_df<-data.frame()
  
  
  # ---------------------------------------------------------------
  #-------------------- Start bootstrap loop ---------------------
  #----------------------------------------------------------------
  for(bootstrap_run in 1:500){
    
    start.time<-Sys.time()
    message("Predicting for ", Model_type,"model: bootstrap run", bootstrap_run)
    
    
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
    names(rasters) <- rasters.names.s # [!]
    data.frame(Rasters=names(rasters),Names=rasters.names.s,Resp=variable.monotonic.response)
    
    
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
      dplyr::group_by(PROV_CODE) %>%
      dplyr::summarize(geometry = sf::st_union(geometry), .groups = "drop")
    
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
    
    # Create depth mask based on depth distribution of species + depth buffer
    depth_mask <- (bathymetry >= -(max.depth + depth.buffer)) & (bathymetry <= min.depth)
    
    #Mask with depth mask
    rasters <- terra::mask(rasters, depth_mask, maskvalues = FALSE)
    
    
    #----------------------------------------------------------------
    #------Write rasters to disk for niche modelling ----------------
    #----------------------------------------------------------------
    
    # #Export rasters so they can be read inside parallel loop
    # interim.folder<-file.path("data", "Interim", Model_type, paste0(Model_type, "rasters"))
    # if(!dir.exists(interim.folder)) dir.create(interim.folder, recursive=TRUE)
    # outfile <- file.path(interim.folder, "predictors_stack.tif")
    # terra::writeRaster(rasters, outfile, overwrite = TRUE)
    
    
    #----------------------------------------------------------------
    #------Run models using a resampling of pseudoabsences and presences --------
    #----------------------------------------------------------------
    #Prepare rasters
    #rasters<-terra::rast(outfile)
    predictors<-names(rasters)
    
    #Resample
    set.seed(bootstrap_run)
    occurrence.resampled <- occurrence.records[sample(seq_len(nrow(occurrence.records)),
                                                      size = nrow(occurrence.records),
                                                      replace = TRUE),]
    set.seed(bootstrap_run)
    background.resampled <- background[sample(seq_len(nrow(background)),
                                              size = nrow(background),
                                              replace = TRUE),]
    
    #Use a resampling of pseudoabsences and presences to train the final model
    train.dataset <- data.frame( PA = c( rep(1,nrow(occurrence.resampled)) , rep(0,nrow(background.resampled)) ) , terra::extract( rasters, rbind( occurrence.resampled, background.resampled), ID=F ) )
    train.dataset[train.dataset == "NaN"] <- NA
    
    #Define tree complexity and learning rate of best model
    cv_folder<-file.path("results", "cross_validation", Model_type,Background_sizes)
    best_metrics <-read.csv2(file.path(cv_folder,paste0("cv_best_metrics",Background_sizes,".csv")))
    best.model.tc <- best_metrics$tree.c
    best.model.lr <-best_metrics$l.rate
    
    #Run final model
    set.seed(1234)
    model <- gbm.step( data=train.dataset, 
                       gbm.x = which(colnames(train.dataset) %in% predictors),
                       gbm.y = 1, 
                       family = "bernoulli", 
                       plot.main = FALSE,
                       tree.complexity = best.model.tc, # you run the model with the tree complexity and learning rate that gave the best TSS in previous step (i.e. with the parameters of the model that was identified as best)
                       learning.rate = best.model.lr, 
                       bag.fraction = 0.5, 
                       n.folds=10,
                       step.size=50,
                       max.trees=2000,
                       silent=TRUE,
                       var.monotone = variable.monotonic.response ,
                       verbose=FALSE)
    
    
    #---------------------------------
    #---------Mark time ---------
    #---------------------------------
    elapsed <- round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), 2)
    message(sprintf("[Model: %s] Completed model fitting in %0.2f minutes", 
                    Model_type,  elapsed))
    
    
    #----------------------------------------------------------------------------
    #----------------------------------Make predictions -------------------------
    #----------------------------------------------------------------------------
    final.prediction.brt <- predict.distribution(rasters, model , reclass.to.one=TRUE)
    
    
    #----------------------------------------------------------------------------
    #----------------------- Reclassify predictions-------------------------
    #----------------------------------------------------------------------------
    num.trees <- model$gbm.call$best.trees
    observed <- train.dataset$PA
    predicted <- predict(model , train.dataset[,-1] , n.trees = num.trees,type="response")
    predictions <- data.frame(observed= observed,
                              predicted= predicted)    
    pred_min <- min(predictions$predicted, na.rm = TRUE)
    pred_max <- max(predictions$predicted, na.rm = TRUE)
    
    # Rescale to 0-1
    predictions$predicted_rescaled <- (predictions$predicted - pred_min) / (pred_max - pred_min)
    
    #Define 10% mtp threshold (lowest probability of a presence)%
    mtp_prob<-0.1
    mtp_threshold<-predictions %>%
      filter(observed == 1) %>%
      summarise(threshold = quantile(predicted_rescaled, probs = mtp_prob, na.rm = TRUE)) %>%
      pull()
    
    #Reclassify
    final.prediction.brt.reclassified<- reclassify.predicted(final.prediction.brt,occurrence.resampled,background.resampled,method="direct.reclass",reclass.threshold=mtp_threshold)
    
    
    #--------------------------------------------------------------------
    #--------------------------Save predictions -------------------------
    #--------------------------------------------------------------------
    results_folder<-file.path("results", "bootstrap_resampling", "final_suitability_maps","Present")
    if(!dir.exists(results_folder))dir.create(results_folder,recursive=TRUE)
    terra::writeRaster(final.prediction.brt, file.path(results_folder,paste0("Final_prediction_brt_",Model_type,"_present_",bootstrap_run,".tif")), overwrite=T)
    terra::writeRaster(final.prediction.brt.reclassified, file.path(results_folder, paste0("Final_prediction_brt_reclassified_",Model_type,"_present_",bootstrap_run,".tif")), overwrite=T)
    
    
    # ------------------------------------------------------------------------------------
    # ----------------------------Variable contribution --------------------------------------
    # ------------------------------------------------------------------------------------
    suppressMessages({
      #Create folder for results
      varimp_folder<-file.path("results", "bootstrap_resampling", "variable_importance")
      if(!dir.exists(varimp_folder))dir.create(varimp_folder,recursive=TRUE)
      
      #Get variable importance dataframe
      varimp_df<-summary.model(model, print.data=T)
      varimp_df<-varimp_df%>%
        mutate(Variable = case_when(grepl("Salinity.Lt.min", Variable) ~ "Min. salinity",
                                    grepl("light", Variable, ignore.case = TRUE) ~ "Mean light",
                                    TRUE ~ Variable))%>%
        mutate(Variable = case_when(grepl (paste0("Temperature.Lt.max|", Model_type, "_LTmax"), Variable) ~ "Max. temperature",
                                    grepl(paste0("Temperature.Lt.min|", Model_type, "_LTmin"), Variable) ~ "Min. temperature",
                                    TRUE ~ Variable))%>%
        mutate(Variable = factor(Variable, levels = c("Min. temperature",
                                                      "Min. salinity",
                                                      "Max. temperature",
                                                      "Mean light")),
               Boostrap_run = bootstrap_run)
      
      if(nrow(varimp_bootstrap_df)==0){
        varimp_bootstrap_df <- varimp_df
      }else{
        varimp_bootstrap_df <- bind_rows(varimp_bootstrap_df, varimp_df)
      }
      
      if(bootstrap_run==500){
        #Export
        write.csv2(varimp_bootstrap_df, file.path(varimp_folder, paste0(Model_type, "_variable_importance.csv")), row.names=FALSE)
      }
      
    })
    # ggplot(varimp_df, aes(x = Percentage, y =Variable)) +
    #  geom_col(fill = "steelblue") +
    #  scale_x_continuous(limits=c(0,100))+
    #  labs( x = "Relative influence (%)",
    #        y = "",
    #        title = Model_type) +
    #  theme_bw()
    # 
  
    
    #--------------------------------------------------------------------
    #----------------------Create and save MESS maps ---------------
    #--------------------------------------------------------------------
    #Create mess map
    mess_present <- dismo::mess(x = raster::stack(rasters),
                                v = train.dataset[,-1])
    mess_present<-rast(mess_present)
    
    #Create folder
    present_mess_folder<-file.path("results", "bootstrap_resampling", "MESS_maps","Present")
    if(!dir.exists(present_mess_folder))dir.create(present_mess_folder,recursive=TRUE)
    
    #Write rasters
    terra::writeRaster(mess_present,file.path(present_mess_folder,paste0("Final_mess_map_",Model_type,"_present_",bootstrap_run,".tif")),overwrite=T)
    
    
    # ------------------------------------------------------------------------------------
    #--Calculate mean suitability, predicted area and bathymetry of the prediction---------
    # ----------------------------------------------------------------------------
    #Create folder for results
    deptharea_folder<-file.path("results", "bootstrap_resampling", "depth_area_overview", "Present")
    if(!dir.exists(deptharea_folder))dir.create(deptharea_folder,recursive=TRUE)
    
    for(Region in c("Full_range", "MedSea")){
      
      #Mask predictions to extent of MedSea
      if(Region=="MedSea"){
        final.prediction.brt<-terra::crop(final.prediction.brt, MedSea)
        final.prediction.brt<-terra::mask(final.prediction.brt, MedSea)
        final.prediction.brt.reclassified <- terra:: crop(final.prediction.brt.reclassified, MedSea)
        final.prediction.brt.reclassified <- terra:: mask(final.prediction.brt.reclassified, MedSea)
      }
      
      # Calculate mean probability of occurrence
      Mean_probability <- global(final.prediction.brt, fun = "mean", na.rm = TRUE)[1,1]
      
      #-----------------------------
      # Calculate predicted area in km2
      cell_areas<- cellSize(final.prediction.brt.reclassified, unit = "km", mask=TRUE) 
      area_raster <- cell_areas * final.prediction.brt.reclassified
      suitable_area <- global(area_raster, fun = "sum", na.rm = TRUE)[1,1]  # km2
      total_area <- global(cell_areas, fun = "sum", na.rm = TRUE)[1,1]  # km2
      
      #-----------------------------
      # Look at bathymetry of your final predictions
      bathy_raster <- rast(bathymetry.file)
      bathy_predicted <- crop(bathy_raster, final.prediction.brt.reclassified) *
        final.prediction.brt.reclassified * (-1)
      
      # Histogram of bathymetry values (excluding zeros)
      bathy_vals <- values(bathy_predicted)
      bathy_vals <- bathy_vals[!is.na(bathy_vals)]
      
      # Summary statistics
      min_depth <- min(bathy_vals)
      max_depth <- max(bathy_vals)
      quantile_90 <- quantile(bathy_vals, probs = 0.9)
      
      #-----------------------------
      # Combine results into a data.frame
      area_depth <- data.frame(
        Mean_occurrence_prob = Mean_probability,
        Total_area = total_area,
        Suitable_area = suitable_area,
        Min_depth = min_depth,
        Max_depth = max_depth,
        Depth_90quantile = quantile_90,
        MTP_prob = mtp_prob,
        MTP_threshold = mtp_threshold,
        learning_rate = best.model.lr,
        tree_complexity = best.model.tc,
        n_trees = num.trees ,
        Time = "Present",
        Model = Model_type,
        Study_area = Region,
        Bootstrap_run = bootstrap_run)
      
      #------Store in bootstrap data frame----
      if(nrow(area_depth_bootstrap_df) == 0){
        area_depth_bootstrap_df <- area_depth
      } else {
        area_depth_bootstrap_df <- bind_rows(area_depth_bootstrap_df, area_depth)
      }
      
      #--------------------
      if(bootstrap_run==500 & Region =="MedSea"){
        #Export
        write.csv(area_depth_bootstrap_df, file.path(deptharea_folder, paste0("depth_area_info_all_",Model_type,"_present.csv")), row.names = F)
      }
    }
    
    
    # -----------------------------------------------------------------
    #--------------------Predict future layers-------------------------
    # -----------------------------------------------------------------
    predict.times <- c("2100_RCP26", "2100_RCP45", "2100_RCP85")
    
    for (predict.time in predict.times){
      
      message("Predicting for ", predict.time)
      
      #-------------------------------------
      #------Prepare rasters----------------
      #-------------------------------------
      rasters.names.future <- switch(Model_type,
                                     "Correlative"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                                      paste0(predict.time,".Benthic.Min.Depth.Salinity.Lt.min"),
                                                      paste0(predict.time,".Benthic.Min.Depth.Temperature.Lt.max"),
                                                      paste0(predict.time,".Benthic.Min.Depth.Temperature.Lt.min")),
                                     "Germination"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                                      paste0(predict.time,".Benthic.Min.Depth.Salinity.Lt.min"),
                                                      paste0("Germination_LTmax_", predict.time),
                                                      paste0("Germination_LTmin_", predict.time)),
                                     "Growth"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                                 paste0(predict.time,".Benthic.Min.Depth.Salinity.Lt.min"),
                                                 paste0("Growth_LTmin_", predict.time),
                                                 paste0("Growth_LTmax_", predict.time)),
                                     "FvFm"= c("light.at.bottom.Benthic.Min.Var.Mean",
                                               paste0(predict.time,".Benthic.Min.Depth.Salinity.Lt.min"),
                                               paste0("FvFm_LTmax_", predict.time),
                                               paste0("FvFm_LTmin_", predict.time))
      )
      
      data.folder.future.rasters  <- file.path("data", "Input", "Rasters", paste0("Layers_future",predict.time))
      future.rasters.files <- list.files(data.folder.future.rasters, pattern="tif", full.names = TRUE)
      
      if(Model_type!="Correlative"){
        future.hybrid.rasters <- list.files(file.path("data", "Input", "Hybrid_rasters", paste0("Layers_future",predict.time)), pattern="tif", full.names = TRUE)
        future.rasters.files<-c(future.rasters.files, future.hybrid.rasters)
      }
      
      future.rasters <- terra::rast(future.rasters.files[as.vector(sapply(rasters.names.future,function(x) { which( grepl(x,future.rasters.files)) } ))])
      
      #Create and export raster overview
      fut_log_folder<-file.path("results", "bootstrap_resampling", "final_suitability_maps", predict.time, "Logs")
      if(!dir.exists(fut_log_folder))dir.create(fut_log_folder,recursive=TRUE)
      
      Overview<-data.frame(Rasters=names(future.rasters),Raster_file_names=rasters.names.future,Resp=variable.monotonic.response,  Time=predict.time, Bootstrap_run=bootstrap_run)
      
      if(nrow(logs_bootstrap_df)==0){
        logs_bootstrap_df <- Overview
      }else{
        logs_bootstrap_df <- bind_rows(logs_bootstrap_df, Overview)
      }
      
      if(bootstrap_run==500){
        write.csv2(logs_bootstrap_df, file.path(fut_log_folder, paste0(Model_type,"_raster_input_log_",predict.time,".csv")), row.names = FALSE)
      }
      
      #--------------------------------------------------------------------
      #----------------------Predict for future layers---------------
      #--------------------------------------------------------------------
      Predicted_future_distribution<- predict.distribution(future.rasters,model,reclass.to.one=TRUE)
      Reclassified_future_distribution<- reclassify.predicted(Predicted_future_distribution,occurrence.resampled,background.resampled,method="direct.reclass",reclass.threshold=mtp_threshold)
      
      
      #--------------------------------------------------------------------
      #----------------------Save future prediction layers ---------------
      #--------------------------------------------------------------------
      fut_pred_folder<-file.path("results", "bootstrap_resampling", "final_suitability_maps", predict.time)
      if(!dir.exists(fut_pred_folder))dir.create(fut_pred_folder,recursive=TRUE)
      
      terra::writeRaster(Predicted_future_distribution,file.path(fut_pred_folder,paste0("Final_prediction_brt_",Model_type,"_",predict.time,"_",bootstrap_run,".tif")),overwrite=T)
      terra::writeRaster(Reclassified_future_distribution,file.path(fut_pred_folder,paste0("Final_prediction_brt_reclassified_",Model_type,"_",predict.time,"_", bootstrap_run,".tif")),overwrite=T)
      
      
      #--------------------------------------------------------------------
      #----------------------Create and save MESS maps ---------------
      #--------------------------------------------------------------------
      #Create mess maps for future and present
      mess_future <- mess(x = raster::stack(future.rasters),
                          v = train.dataset[,-1])
      mess_future<-rast(mess_future)
      
      #Create folder
      fut_mess_folder<-file.path("results", "bootstrap_resampling", "MESS_maps", predict.time)
      if(!dir.exists(fut_mess_folder))dir.create(fut_mess_folder,recursive=TRUE)
      
      #Write rasters
      terra::writeRaster(mess_future,file.path(fut_mess_folder,paste0("Final_mess_map_",Model_type,"_",predict.time,"_",bootstrap_run,".tif")),overwrite=T)
      
      
      # ------------------------------------------------------------------------------------
      #--Calculate mean suitability, predicted area and bathymetry of the prediction---------
      # ----------------------------------------------------------------------------
      #Create folder for results
      future_deptharea_folder<-file.path("results", "bootstrap_resampling", "depth_area_overview", predict.time)
      if(!dir.exists(future_deptharea_folder))dir.create(future_deptharea_folder,recursive=TRUE)
      
      #Loop over region
      for(Region in c("Full_range", "MedSea")){
        
        #Mask predictions to extent of MedSea
        if(Region=="MedSea"){
          Predicted_future_distribution<-terra::crop(Predicted_future_distribution, MedSea)
          Predicted_future_distribution<-terra::mask(Predicted_future_distribution, MedSea)
          Reclassified_future_distribution <- terra:: crop(Reclassified_future_distribution, MedSea)
          Reclassified_future_distribution <- terra:: mask(Reclassified_future_distribution, MedSea)
        }
        
        # Calculate mean probability of occurrence
        Mean_probability <- global(Predicted_future_distribution, fun = "mean", na.rm = TRUE)[1,1]
        
        #-----------------------------
        # Calculate predicted area in km2
        cell_areas<- cellSize(Reclassified_future_distribution, unit = "km", mask=TRUE) 
        area_raster <- cell_areas * Reclassified_future_distribution
        suitable_area <- global(area_raster, fun = "sum", na.rm = TRUE)[1,1]  # km2
        total_area <- global(cell_areas, fun = "sum", na.rm = TRUE)[1,1]  # km2
        
        
        #-----------------------------
        # Look at bathymetry of your final predictions
        bathy_raster <- rast(bathymetry.file)
        bathy_predicted <- crop(bathy_raster, Reclassified_future_distribution) *
          Reclassified_future_distribution * (-1)
        
        # Histogram of bathymetry values (excluding zeros)
        bathy_vals <- values(bathy_predicted)
        bathy_vals <- bathy_vals[!is.na(bathy_vals)]
        
        # Summary statistics
        min_depth <- min(bathy_vals)
        max_depth <- max(bathy_vals)
        quantile_90 <- quantile(bathy_vals, probs = 0.9)
        
        #-----------------------------
        # Combine results into a data.frame
        area_depth_future <- data.frame(
          Mean_occurrence_prob = Mean_probability,
          Total_area = total_area,
          Suitable_area = suitable_area,
          Min_depth = min_depth,
          Max_depth = max_depth,
          Depth_90quantile = quantile_90,
          MTP_prob = mtp_prob,
          MTP_threshold = mtp_threshold,
          Time = predict.time,
          Model = Model_type,
          Study_area = Region,
          Bootstrap_run = bootstrap_run)
        
        #--------------------
        #Store bootstrap data in same dataframe
        if(nrow(area_depth_future_bootstrap_df)==0){
          area_depth_future_bootstrap_df <- area_depth_future
        }else{
          area_depth_future_bootstrap_df<- bind_rows(area_depth_future_bootstrap_df, area_depth_future)
        }
        
        
      }
      
    }
    
    #---------------------------------
    #---------Mark time ---------
    #---------------------------------
    elapsed <- round(as.numeric(difftime(Sys.time(), start.time, units = "mins")), 2)
    message(sprintf("[Model: %s | Bootstrap run: %s] Completed bootstrap round in %0.2f minutes", 
                    Model_type, bootstrap_run, elapsed))
    
  }
  
  #---------------------------------------
  # Export once per Model_type
  #---------------------------------------
  
  final_future_folder <- file.path("results",
                                   "bootstrap_resampling",
                                   "depth_area_overview",
                                   "Future_all")
  
  if(!dir.exists(final_future_folder))
    dir.create(final_future_folder, recursive = TRUE)
  
  write.csv(
    area_depth_future_bootstrap_df,
    file.path(final_future_folder,
              paste0("depth_area_info_all_", Model_type, "_future.csv")),
    row.names = FALSE
  )
  
  
}
