#-------------------------------
#-----   model.metrics   ------
#-------------------------------
model.metrics <- function(predictions, predicted.distribution) {
  
  #Get observed and predicted values
  observed <- predictions$observed
  predicted <- predictions$predicted
  
  #Get predicted values for occurrences and background
  predicted_occs <- predictions%>%
    dplyr::filter(observed == 1)%>%
    dplyr::pull(predicted)
  
  predicted_bg <- predictions%>%
    filter(observed == 0)%>%
    pull(predicted)
  
  
  #----------Calculate Boyce Index--------
  
  #Extract values in each cell of the raster
  all_suit_vals <- terra::values(predicted.distribution)
  all_suit_vals <- all_suit_vals[!is.na(all_suit_vals)]
  
  #Compute Boyce only if there are enough occurrences
  if (length(predicted_occs) > 0) {
    boyce_result <- ecospat.boyce(fit = all_suit_vals,
                                  obs = predicted_occs,
                                  nclass = 0) #Moving window boyce index
    boyce_val <- round(boyce_result$cor, 3)
  }else{
    boyce_val <-NA
    warning("No occurrence points available to calculate Boyce index")
  }
  
  
  #----------Calculate AUC----------
  
  if (length(predicted_occs) > 0 & length(predicted_bg) > 0) {
    auc_eval <- dismo::evaluate(p = predicted_occs, a = predicted_bg)
    auc_val <- ifelse(is.null(auc_eval@auc), NA, auc_eval@auc)
  } else {
    auc_val <- NA
    warning("Not enough points to calculate AUC")
  }
  
  
  #---------return---------
  predicted.accuracy <- data.frame( auc = auc_val,
                                    boyce = boyce_val)
  return(predicted.accuracy)
  
}


#-------------------------------
#---run_brt_crossvalidation-----
#-------------------------------
run_brt_crossvalidation <- function(c, comb, cross.validation.data, presences, background.points, raster.file, variable.monotonic.response.i){
  
  
  gc() #Free up memory
  
  cv <- comb[c,1]
  l.rate <- comb[c,2]
  tree.c <- comb[c,3]
  
  if( TRUE %in% grepl("lat", colnames(cross.validation.data)) ) {
    presences.train <- presences[ presences$Lat < cross.validation.data[cv,1] | presences$Lat > cross.validation.data[cv,2] , ]
    background.points.train <- background.points[ background.points$Lat < cross.validation.data[cv,1] | background.points$Lat > cross.validation.data[cv,2] , ]
    presences.test <- presences[ presences$Lat >= cross.validation.data[cv,1] & presences$Lat <= cross.validation.data[cv,2] , ]
    background.points.test <- background.points[ background.points$Lat >= cross.validation.data[cv,1] & background.points$Lat <= cross.validation.data[cv,2] , ]
  }
  
  if( TRUE %in% grepl("lon", colnames(cross.validation.data)) ) {
    presences.train <- presences[ presences$Lon <= cross.validation.data[cv,1] | presences$Lon >= cross.validation.data[cv,2] , ]
    background.points.train <- background.points[ background.points$Lon <= cross.validation.data[cv,1] | background.points$Lon >= cross.validation.data[cv,2] , ]
    presences.test <- presences[ presences$Lon >= cross.validation.data[cv,1] & presences$Lon <= cross.validation.data[cv,2] , ]
    background.points.test <- background.points[ background.points$Lon >= cross.validation.data[cv,1] & background.points$Lon <= cross.validation.data[cv,2] , ]
  }  
  
  #Remove NAs  
  presences.train <- presences.train[complete.cases(presences.train), ]
  background.points.train <- background.points.train[complete.cases(background.points.train), ]
  presences.test <- presences.test[complete.cases(presences.test), ]
  background.points.test <- background.points.test[complete.cases(background.points.test), ]
  
  if( nrow(presences.train) > 0 &nrow(background.points.train) > 0 & nrow(presences.test) > 0 & nrow(background.points.test) > 0 ) { 
    
    #Load rasters1
    rasters1<-terra::rast(raster.file)
    predictors <- names(rasters1)
    
    #Prepare datasets    
    train.dataset <- data.frame( 
      PA = c( rep(1,nrow(presences.train)) , rep(0,nrow(background.points.train)) ) , 
      terra::extract( rasters1 , rbind( presences.train, background.points.train), ID=F ) )
    train.dataset[train.dataset == "NaN"] <- NA
    colnames(train.dataset) <- c("PA",names(rasters1))
    
    test.dataset <- data.frame( 
      PA = c( rep(1,nrow(presences.test)) , rep(0,nrow(background.points.test)) ) , 
      terra::extract( rasters1 , rbind( presences.test, background.points.test), ID=F ) )
    test.dataset[test.dataset == "NaN"] <- NA
    colnames(test.dataset) <- c("PA",names(rasters1))
    
    #Fit brt model
    set.seed(200)
    model <- tryCatch({
      dismo::gbm.step( 
        data=train.dataset, 
        gbm.x = which(colnames(train.dataset) %in% predictors),
        gbm.y = 1, 
        family = "bernoulli", 
        plot.main = FALSE,
        tree.complexity = tree.c, 
        learning.rate = l.rate, 
        bag.fraction = 0.5, 
        n.folds=10,
        step.size=50,
        max.trees=2000,
        silent=TRUE,
        n.cores = 1,
        var.monotone = variable.monotonic.response.i ,
        verbose=FALSE)
    }, error = function(e) NULL)
    
    
    if( ! is.null(model)) {
      
      num.trees <- model$gbm.call$best.trees
      observed <- test.dataset$PA
      predicted <- predict(model , test.dataset[,-1] , n.trees = num.trees,type="response")

      predictions <- data.frame(observed= observed,
                                predicted= predicted)  
      
      #Get predicted distribution for the test region (used for boyce index calculation)
      # Define extent of test region
      test_extent <- terra::ext(ext(rasters1)[1], ext(rasters1)[2],
                                cross.validation.data[cv,1], cross.validation.data[cv,2])
      
      # Crop raster to test region
      rasters_cropped <- terra::crop(rasters1, test_extent)
      predicted.distribution <- predict(rasters_cropped, model , n.trees=num.trees,type="response", na.rm=TRUE)  
      
      predicted.accuracy <- model.metrics(predictions,
                                          predicted.distribution)
      predicted.accuracy <- data.frame( cv.round=cv,
                                        tree.c=tree.c,
                                        l.rate=l.rate, 
                                        n.trees= num.trees,
                                        predicted.accuracy)
      gc()
      rm(predicted.distribution, rasters_cropped)
      message(sprintf("[CV %d] Fitted combination: lr=%0.4f, tc=%d", cv, l.rate, tree.c))
      return(predicted.accuracy)
      
    }
  }
  # If anything fails, return NA row
  message(sprintf("[CV %d] Failed combination: lr=%0.4f, tc=%d", cv, l.rate, tree.c))
  return(data.frame(cv.round = cv, 
                    tree.c = tree.c, 
                    l.rate = l.rate, 
                    n.trees = NA,
                    auc = NA, 
                    boyce = NA))
  
}
