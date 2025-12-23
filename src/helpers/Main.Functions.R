# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
# One Pipeline for Modelling the distribution of marine Species
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

packages.to.use <- c(   "ecospat",
                        "future",
                        "future.apply",
                        "dplyr",
                        "terra",
                        "sf",
                        "ggplot2",
                        "spatstat",
                        "mboost",
                       "ecodist"
                       , "raster"
                       , "gdata"
                       , "leaflet"
                       , "dismo"
                       , "gbm"
                       , "sm" 
                       , "sp"
                       , "biomod2"
                       , "adehabitatHS"
                       #, "SDMTools"
                       , "parallel"
                       , "doParallel"
                       , "biganalytics"
                       , "nicheROVER"
                       , "vegan"
                       , "parallel"
                       #, "rgeos"
                       #, "rgdal"
                       , "sdmpredictors"
                       , "usdm"
                       , "doSNOW"
                       , "ENMeval"
                       #, "maptools"
                       , "rgl"
                       #, "rpanel"
                       , "rJava"
                       , "gstat"
                       ,"spThin")

#

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## -----------------------

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

## -----------------------

plot.google <- function(coordinates,radius,color) {
  
  set.seed(42)
  
  m <- leaflet()
  m <- addTiles(m)
  # m <- addMarkers(m, lng=coordinates[,1], lat=coordinates[,2], popup=paste0( "Species record ") , icon = greenLeafIcon)
  
  m <- addCircleMarkers(m, lng=coordinates[,1], lat=coordinates[,2], 
                        popup=paste0( "Species record ") , 
                        radius = radius, color = color , 
                        stroke = FALSE, fillOpacity = 0.5 )
  
  return(m)
  
}


## -----------------------

list.sdm.predictors <- function() {
  
  raster.predictors <- list_layers()
  raster.predictors <- raster.predictors[raster.predictors$dataset_code == "Bio-ORACLE",2:3]
  return(raster.predictors)
}

## -----------------------

correlated.pairs <- function(rasters,threhold) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasters))-1)) {
    
    for( j in (i+1):length(names(rasters))) {
      
      corr.val <- abs(cor(getValues(subset(rasters,i)),getValues(subset(rasters,j)),use = "pairwise.complete.obs"))
      
      if( corr.val >= threhold) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasters,i)),
                                                                               Var.2=names(subset(rasters,j)),
                                                                               Cor=corr.val,stringsAsFactors = FALSE)) }
      
    }
    
  }
  
  return(list.of.cor)
  
  
}

## -----------------------

relocate.coordinates.na <- function(occurrence.records,rasters,maximum.distance,use.species.depth,min.depth,max.depth) {
  
  set.seed(42)
  
  if(missing(use.species.depth)) { use.species.depth <- FALSE }
  
  if( use.species.depth ) {
    
    min.depth <- min.depth * (-1)
    rclmat <- matrix( c(max.depth*(-1), min.depth, 1, -9999999999 , max.depth*(-1) , NA , min.depth , 9999999999 , NA) , ncol=3, byrow=TRUE)
    
    polygon.to.intersect <- as(extent(c(min(occurrence.records$Lon)-2,max(occurrence.records$Lon)+2,min(occurrence.records$Lat)-2,max(occurrence.records$Lat)+2)), "SpatialPolygons")
    
    bathymetry <- raster(bathymetry.file)
    bathymetry <- crop(bathymetry, polygon.to.intersect)
    bathymetry <- reclassify(bathymetry, rclmat)
    
    to.relocate <- which(is.na(raster::extract(bathymetry,occurrence.records[,c("Lon","Lat")])))
    coordinates.to.relocate <- occurrence.records[to.relocate,]
    correct.points <- as.data.frame(bathymetry,xy=TRUE,na.rm=TRUE)[,1:2]
    
  } 
  
  if( ! use.species.depth ) {  
    
    to.relocate <- unique(which(is.na(raster::extract(rasters,occurrence.records[,c("Lon","Lat")])) , arr.ind = TRUE)[,1])
    coordinates.to.relocate <- occurrence.records[to.relocate,]
    correct.points <- as.data.frame(subset(rasters,1),xy=TRUE,na.rm=TRUE)[,1:2]
    
  }
  
  if( nrow(coordinates.to.relocate) > 0 ) { 
    
    cat( paste0("Relocating ",length(to.relocate)," Points that were falling out of range"))
    cat( paste0("\n"))
    
    near.cells <- numeric(nrow(coordinates.to.relocate))
    
    for(p in 1:nrow(coordinates.to.relocate)) {

      near.cell.p <- spDistsN1( as.matrix(correct.points), as.matrix(coordinates.to.relocate[p,c("Lon","Lat")]),longlat=TRUE)
      
      if( near.cell.p[which.min(near.cell.p)] <= maximum.distance ) {
        
        near.cell.p <- which.min(near.cell.p)
        
      } else {   near.cell.p <- NA }
      
      
      near.cells[p] <- near.cell.p
      
    }
    
    relocated <- which(!is.na(near.cells))
    
    if( length(relocated) > 0) {
      
      near.cells <- correct.points[near.cells[relocated],]
      old.presences <- occurrence.records[-to.relocate[relocated],]
      
      colnames(old.presences) <- c("Lon","Lat")
      colnames(near.cells) <- c("Lon","Lat")
      
      occurrence.records <- rbind(old.presences[,c("Lon","Lat")],near.cells)
      
    }
    
  }
  
  if( nrow(coordinates.to.relocate) == 0) { 
    
    cat( paste0("None to Relocate"))
    cat( paste0("\n"))
    
  }
  
  if( use.species.depth ) {
    
    to.remove <- which(is.na(raster::extract(bathymetry,occurrence.records)))
    
  }
  
  if( ! use.species.depth ) {
    
    to.remove <- unique(which(is.na(raster::extract(rasters,occurrence.records[,c("Lon","Lat")])), arr.ind = TRUE)[,1])
    
  }
  
  if( length(to.remove) > 0) { 
    
    occurrence.records <- occurrence.records[-to.remove,]
    
  }
  
  ## -----------------------
  
  return( occurrence.records )
  
}


## -----------------------

spatial.autocorrelation.thinning <- function(occurrence.records,min.distance) {
  
  coordinates.t <- occurrence.records
  
  space <- spDists(as.matrix(coordinates.t),as.matrix(coordinates.t),longlat=TRUE)
  diag(space) <- NA
  
  reclass <- space <= min.distance
  reclass[lower.tri(reclass, diag=TRUE)] <- NA
  
  v <- colSums(reclass, na.rm=TRUE) == 0
  coordinates.t <- coordinates.t[v,]
  
  # Number of All occurrences and number to keep
  
  cat( paste0("\n"))
  cat( paste0("\n"))
  
  cat( paste0("Input Records: ",nrow(occurrence.records)))
  cat( paste0("\n"))
  cat( paste0("Final Records: ",nrow(coordinates.t)))
  
  # Remove from main dataset of occurrences
  
  return(coordinates.t)
  
}

## -----------------------

sample.background <- function(occurrence.records,rasters,polygon.pa,polygon.pa.type,pa.dist.min,pa.type,bias.surface.kernel,n.background,pa.environment.stratification,plot.figure, bias.surface.path) {
  
  old_warn <- getOption("warn")
  on.exit(options(warn = old_warn), add = TRUE)
  options(warn = -1)
  
  ## ------------
  
  if( ! is.null(polygon.pa) ) { 
    if (inherits(polygon.pa, "sf")) polygon.pa <- terra::vect(polygon.pa)
    if(polygon.pa.type == "Include") { rasters_masked <- terra::mask(rasters,polygon.pa)  } 
    
    if(polygon.pa.type == "Exclude") { 
      #Mask raster with polygon
      rasters_masked <-terra::mask(rasters, polygon.pa, inverse=TRUE) 
    }
  }
 
  
  # -----------------
  # Create a start set of sink points: all pixels that have values in our raster layer (i.e. not too deep and inside range that we cropped them) and without counting pixels that fall in the polygon 
  sink.points <- as.data.frame(rasters_masked[[1]],xy=TRUE)
  sink.points <- sink.points[!is.na(sink.points[,3]),1:2]  #63888
  
  ## ------------
  
  # Removes points closer to a presence than pa.dist.min 
  if( ! is.null(pa.dist.min) ) {
    
    # Convert sink points and occurrence records to sf
    sink_sf <- sf::st_as_sf(sink.points, coords = c("x", "y"), crs = 4326)  # WGS84
    occ_sf <- sf::st_as_sf(occurrence.records, coords = c("Lon", "Lat"), crs = 4326)
    
    # st_is_within_distance uses meters when CRS is lon/lat 
    to_remove_idx_list <- sf::st_is_within_distance(sink_sf, occ_sf, dist = pa.dist.min * 1000)
    to_remove <- lengths(to_remove_idx_list) > 0
    
    # Keep only points that are farther than pa.dist.min: 55112
    sink.points <- sink.points[!to_remove, , drop = FALSE]
  }
  
  
  ## ------------
  
  if( exists("distance.uncorr") ) {
    sink_sf <- st_as_sf(sink.points, coords = c("x", "y"), crs = 4326)
    
    # randomize order to avoid bias from original ordering 
    sink_sf <- sink_sf[sample(nrow(sink_sf)), , drop = FALSE]
    
    # prepare
    keep_idx <- logical(nrow(sink_sf))
    remaining <- seq_len(nrow(sink_sf))
    
    # distance threshold in meters as a units object
    d_thresh <- units::set_units(distance.uncorr * 1000, "m")  
    
    # distances from chosen point to remaining candidates
    neighbors <- st_is_within_distance(sink_sf, sink_sf, dist = d_thresh)
    
    while (length(remaining) > 0) {
      i <- remaining[1]
      keep_idx[i] <- TRUE
      # remove all that are within distance of i
      too_close <- intersect(neighbors[[i]], remaining)
      remaining <- setdiff(remaining, too_close)
    }
    
    sink.points <- as.data.frame(st_coordinates(sink_sf[keep_idx, , drop = FALSE])) #41907
    colnames(sink.points) <- c("x", "y")
  }
  
  
  ## ------------
  sink.points.p <- NULL   # probabilities for sink points if computed
  
  # Generates bias surface kernel
  if( bias.surface.kernel ) {
  
    bias.surface <- terra::rast(bias.surface.path)
    
    #Mask with raster
    bias.surface <- terra::crop(bias.surface, rasters_masked[[1]] )
    bias.surface <- terra::mask(bias.surface, rasters_masked[[1]] )
    
    #Rescale so bias raster ranges from 1-20
    min_val <- global(bias.surface, fun = "min", na.rm = TRUE)[[1]]
    max_val <- global(bias.surface, fun = "max", na.rm = TRUE)[[1]]
    bias.surface <- ((bias.surface- min_val) / (max_val - min_val)) * 19 + 1
    
    sink.points.p <- terra::extract(bias.surface,sink.points, ID=F)%>% #extract occ. prob for each sink point based on the probability surface (Bias.surface)
      pull()
    keeppos <- !is.na(sink.points.p) 
    sink.points <- sink.points[keeppos, , drop = FALSE]
    sink.points.p <- sink.points.p[keeppos]
    
  }
  
  ## ------------
  pa.environment <- NULL
  
  # Generates environmental surface for stratification
  if( pa.environment.stratification ) {
    
    #Extract environmental data
    pa.environment <- terra::extract( rasters_masked,sink.points[,1:2], ID=F ) %>%
      as.data.frame()
    
    #Define number of clusters
    to.generate <- min(n.background , nrow(unique(pa.environment))) #lowest value of the two
    
    #Perform kmeans clustering
      set.seed(seed)
      pa.environment <- kmeans(pa.environment,centers=to.generate,iter.max = 100, nstart = 30)$cluster # we want 10000 clusters. The data given by x are clustered by the k-means method, which aims to partition the points into 10000 groups such that the sum of squares from points to the assigned cluster centres is minimized. At the minimum, all cluster centres are at the mean of their Voronoi sets (the set of data points which are nearest to the cluster centre). you get a vector of sink points which get assigned a number for the cluster they are partitioned in 
      pa.clusters <- unique( pa.environment )
      

  }
  
  ## ------------
  
  colnames(occurrence.records) <- c("Lon","Lat")
  
  ## ------------
  
  n.background <- min(n.background,nrow(sink.points))
  
  
  if(  bias.surface.kernel ) { 
    
    if( ! pa.environment.stratification ) {
      r.sample <- sample(1:nrow(sink.points),n.background,replace = FALSE, prob = sink.points.p)
    }
    
    if( pa.environment.stratification ) {
      set.seed(seed)
      r.sample <- sapply(1:n.background, function(x) { id.x <- which(pa.environment == x) ; id.x[sample(1:length(id.x),1,replace = FALSE, prob = sink.points.p[id.x])]  } )
    }
  }
  
  absences <- sink.points[r.sample,]%>%
    dplyr::rename(Lon = x, 
                  Lat = y)
  # ------------
  
  if( sum(is.na(absences$Lon)) > 0) {
    absences <- absences[!is.na(absences$Lon),]
  }
  
  if(plot.figure) {
    
    plot(absences[,1:2],col="grey", pch=21 , cex=1 , main="Dataset : Presences and Pseudo-absences")
    points(occurrence.records,col="black", pch=21 , cex=1)
  }
  
  #Clean up
  rm(r.sample, sink.points, pa.environment, occ_sf, sink_sf, to_remove_idx_list)
  
  return(absences)
  
}

## -----------------------
sample.background.large<-function(occurrence.records, raster_template, n_background, seed){
  
  occs <- terra::vect(occurrence.records, geom = c("Lon", "Lat"), crs = crs(raster_template))
  cells_with_occurrences <- terra::cellFromXY(raster_template, terra::crds(occs))
  raster_template[cells_with_occurrences] <- NA #64477
  
  set.seed(seed)
  global_points <- terra::spatSample(
    raster_template,
    size = n_background, #three times the number we need
    method = "random",     # weighted random sampling
    as.points = TRUE,       # return SpatVector of points
    na.rm = TRUE            # ignore NA pixels
  )
  global_points <- as.data.frame(crds(global_points))%>%
    rename(Lon = x, 
           Lat = y)
  
  return(global_points)
}


## -----------------------

data.partitioning <- function(presences,absences,type,k) {
  
  if( type=="blocks.longitudinal" ) {
    
    min.Lon <- min( presences$Lon )
    min.Lon.abs <- min(absences$Lon)
    max.Lon <- max( presences$Lon )
    max.Lon.abs <- max(absences$Lon)
    
    min.Lon <- min(min.Lon, min.Lon.abs)
    max.Lon <- max(max.Lon, max.Lon.abs)
    
    bands <- seq(from=min.Lon,to=max.Lon,length.out = k+1)
    bands <- data.frame(lon.from=bands[-length(bands)],lon.to=bands[-1])
    
  }
  
  if( type=="blocks.latitudinal" ) {
    
    min.lat <- min( presences$Lat )
    min.lat.abs <- min(absences$Lat)
    max.lat <- max( presences$Lat )
    max.lat.abs <- max(absences$Lat)
    
    min.lat <- min(min.lat, min.lat.abs)
    max.lat <- max(max.lat, max.lat.abs)

    bands <- seq(from=min.lat,to=max.lat,length.out = k+1)
    bands <- data.frame(lat.from=bands[-length(bands)],lat.to=bands[-1])
    
    
  }
  
  return( bands )
  
}


## -----------------------

summary.model <- function(model,print.data) {

  
  if( class(model) == "gbm" ) { 
    
    data.to.plot <- summary(model)
    colnames(data.to.plot) <- c("Variable","Percentage")
    rownames(data.to.plot) <- NULL
    
    if(print.data) { print(data.to.plot) }
    
  }else{
    data.to.plot <-NULL
  }
  # Plot variable importance 
  return(data.to.plot)
  
}


## -----------------------

model.plot <- function(model_type,
                       model,
                       training.data,
                       predictor.variable,
                       mintemp.seq,
                       maxtemp.seq){
  varname <- model$var.names[predictor.variable]
  nt <- model$gbm.call$best.trees
  
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  # Check if predictor exists in data
  if(!varname %in% names(training.data)) stop("Predictor not found in provided data.")
  
  # Get predictor values
  if(varname == paste0(Model_type,"_LTmin_present")){
    predictor.values<-mintemp.seq[Model_type][[1]]
  }else if (varname == paste0(Model_type,"_LTmax_present")){
    predictor.values<-maxtemp.seq[Model_type][[1]]
  }else{
    predictor.values<-seq(min(training.data[[varname]]), max(training.data[[varname]]),length.out =  100)
  }
  
  ## ---- compute partial dependence explicitly ----
  X <- data.frame(predictor.values)
  names(X) <- varname
  
  # Call gbm internal C function to compute partial dependence
  effect <- .Call("gbm_plot",
                  X = as.double(data.matrix(X)),
                  cRows = as.integer(nrow(X)),
                  cCols = as.integer(ncol(X)),
                  n.class = as.integer(model$num.classes),
                  i.var = as.integer(predictor.variable - 1),
                  n.trees = as.integer(model$gbm.call$best.trees),
                  initF = as.double(model$initF),
                  trees = model$trees,
                  c.splits = model$c.splits,
                  var.type = as.integer(model$var.type),
                  PACKAGE = "gbm")
  
  # ---- center effect (matches dismo::gbm.plot) ----
  effect <- scale(effect, scale=FALSE)
  
  # ---- convert to probability  ----
  effect <- logit2prob(effect)
  
  #---Store temperature data instead of trait values at temperature
  if(varname == paste0(Model_type,"_LTmin_present")){
    temp.values<-mintemp.seq["Temperature"][[1]]
  }else if (varname == paste0(Model_type,"_LTmax_present")){
    temp.values<-maxtemp.seq["Temperature"][[1]]
  }else{
    temp.values<-predictor.values
  }
  
  # ---- return data frame ----
  effect_df<-data.frame(Variable = predictor.values, Effect = effect, Variable_for_plotting = temp.values)
  return(effect_df)
}


## -----------------------

predict.distribution <- function(rasters,model,reclass.to.one) {

  
  if( class(model)[1] == "gbm" ) {
    
    num.tress <- model$gbm.call$best.trees
    if(is.null(num.tress)) { num.tress <- length(model$trees) }
    predicted.distribution <- predict( rasters , model , n.trees=num.tress,type="response", na.rm=TRUE)
  
  }
  
  if(reclass.to.one) {
    predicted.distribution <- predicted.distribution + (global(predicted.distribution, "min", na.rm = TRUE)[1,1] * (-1)) #add the minimum value (negative probability) to the predicted distribution so you don't have a probability of occurrence lower than 0
    predicted.distribution <- predicted.distribution / global(predicted.distribution, "max", na.rm = TRUE)[1,1] #Divide the predicted distribution by the maximum observed values so all values will be between 0 and 1 
  }
  
  return(predicted.distribution)
  
}

## -----------------------

accuracy.predicted <- function(predicted.distribution,presences,absences,type) {
  
  observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
  predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
  
  predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )
  
  if(type == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
  if(type == "tss") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
  if(type == "Kappa") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$Kappa),] }
  
  
  predicted.accuracy <- data.frame( threshold =  predicted.accuracy$threshold ,
                                    auc = predicted.accuracy$AUC ,
                                    specificity = predicted.accuracy$specificity ,
                                    sensitivity = predicted.accuracy$sensitivity ,
                                    tss = predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 )
  
  return(predicted.accuracy)
  
}

## -----------------------

reclassify.predicted <- function(predicted.distribution,presences,absences,method,reclass.threshold) {
  
  if( missing(presences) ) { presences <- NULL  }
  if( missing(absences) ) { absences <- NULL  }
  
  if(method == "direct.reclass") {
    
    
    predicted.distribution[predicted.distribution > reclass.threshold] <- 1
    predicted.distribution[predicted.distribution <= reclass.threshold] <- 0
    
  }
  
  if(method == "max.tss") {
    
    observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
    predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
    
    predicted.accuracy <- accuracy( observed , predicted , threshold =100 )
    predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  
  return(predicted.distribution)
}

## -----------------------
get.fitted.temperature<-function(Trait, Tmin, Topt, Tmax, target, trait_fit, Model_type){
  
trait_fun <- function(T) {
  stopifnot(length(T) == 1)
  as.numeric(predict(trait_fit, newdata = data.frame(Temperature = T)))
}

root_fun <- function(T) trait_fun(T) - target

# low side of Topt
T_low <- tryCatch(
  {uniroot(root_fun, lower = Tmin, upper = Topt)$root},
  error = function(e) {
    NA
  }
)

# high side of Topt
T_high <- tryCatch(
  {uniroot(root_fun, lower = Topt, upper = Tmax)$root},
  error = function(e) {
    NA
  }
)

round(c(T_low, T_high), 6)
output_df<-data.frame(Temperature = c(T_low, T_high), Value = c(target, target), Predictor_name=c(Trait, Trait), Model_type=c(Model_type, Model_type) )
output_df<-output_df[!is.na(output_df$Temperature),]
return(output_df)

}
