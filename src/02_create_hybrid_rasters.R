#--------------------------------------------------
#--------------- Clean up workspace----------------
#--------------------------------------------------
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)


#--------------------------------------------------
#--------------- Load packages ----------------
#--------------------------------------------------
library(dplyr)
library(nlstools)
library(minpack.lm)
library(zen4R)
library(readr)
library(terra)


#--------------------------------------------------
#--------------- Create folders  ----------------
#--------------------------------------------------
input_folder<-file.path("data", "Input", "Physiology")
output_folder<- file.path("data", "Output")

for(folder in c(input_folder, output_folder)){
  if(!dir.exists(folder)) dir.create(folder, recursive=TRUE)
}


#--------------------------------------------------
#-------Load physiological data Files -------------
#--------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.17865936", 
                       path=input_folder,
                       files=c("Dictyota_growth.csv",
                               "Dictyota_FvFm.csv",
                               "Dictyota_germination.csv"),
                       quiet=FALSE)

Growth <- read.csv2(file.path(input_folder,"Dictyota_growth.csv"))
FvFm <- read.csv2(file.path(input_folder,"Dictyota_FvFm.csv"))
Germination <- read.csv2(file.path(input_folder,"Dictyota_germination.csv"))


#--------------------------------------------------------------------------------
#----- Create hybrid rasters for LTmaxTemperature and LTminTemperature-----------
#--------------------------------------------------------------------------------
raster_folder<-(file.path("data", "Input", "Rasters"))

LTmaxT<-rast(file.path(raster_folder, "Layers_present","Present.Benthic.Min.Depth.Temperature.Lt.max.tif"))
LTminT<-rast(file.path(raster_folder, "Layers_present","Present.Benthic.Min.Depth.Temperature.Lt.min.tif"))

LTmaxT_2100_RCP26<-rast(file.path(raster_folder, "Layers_future2100_RCP26", "2100AOGCM.RCP26.Benthic.Min.Depth.Temperature.Lt.max.tif.BOv2_1.tif"))
LTminT_2100_RCP26<-rast(file.path(raster_folder, "Layers_future2100_RCP26", "2100AOGCM.RCP26.Benthic.Min.Depth.Temperature.Lt.min.tif.BOv2_1.tif"))

LTmaxT_2100_RCP45<-rast(file.path(raster_folder,"Layers_future2100_RCP45", "2100AOGCM.RCP45.Benthic.Min.Depth.Temperature.Lt.max.tif.BOv2_1.tif"))
LTminT_2100_RCP45<-rast(file.path(raster_folder,"Layers_future2100_RCP45", "2100AOGCM.RCP45.Benthic.Min.Depth.Temperature.Lt.min.tif.BOv2_1.tif"))

LTmaxT_2100_RCP85<-rast(file.path(raster_folder,"Layers_future2100_RCP85", "2100AOGCM.RCP85.Benthic.Min.Depth.Temperature.Lt.max.tif.BOv2_1.tif"))
LTminT_2100_RCP85<-rast(file.path(raster_folder,"Layers_future2100_RCP85", "2100AOGCM.RCP85.Benthic.Min.Depth.Temperature.Lt.min.tif.BOv2_1.tif"))

LTmaxT_2050_RCP26<-rast(file.path(raster_folder,"Layers_future2050_RCP26", "2050AOGCM.RCP26.Benthic.Min.Depth.Temperature.Lt.max.tif.BOv2_1.tif"))
LTminT_2050_RCP26<-rast(file.path(raster_folder,"Layers_future2050_RCP26", "2050AOGCM.RCP26.Benthic.Min.Depth.Temperature.Lt.min.tif.BOv2_1.tif"))

LTmaxT_2050_RCP45<-rast(file.path(raster_folder,"Layers_future2050_RCP45", "2050AOGCM.RCP45.Benthic.Min.Depth.Temperature.Lt.max.tif.BOv2_1.tif"))
LTminT_2050_RCP45<-rast(file.path(raster_folder,"Layers_future2050_RCP45", "2050AOGCM.RCP45.Benthic.Min.Depth.Temperature.Lt.min.tif.BOv2_1.tif"))

LTmaxT_2050_RCP85<-rast(file.path(raster_folder,"Layers_future2050_RCP85", "2050AOGCM.RCP85.Benthic.Min.Depth.Temperature.Lt.max.tif.BOv2_1.tif"))
LTminT_2050_RCP85<-rast(file.path(raster_folder, "Layers_future2050_RCP85", "2050AOGCM.RCP85.Benthic.Min.Depth.Temperature.Lt.min.tif.BOv2_1.tif"))

# Convert names to Temperature
names(LTmaxT)<-"Temperature"
names(LTminT)<-"Temperature"
names(LTmaxT_2100_RCP26)<-"Temperature"
names(LTminT_2100_RCP26)<-"Temperature"
names(LTmaxT_2100_RCP45)<-"Temperature"
names(LTminT_2100_RCP45)<-"Temperature"
names(LTmaxT_2100_RCP85)<-"Temperature"
names(LTminT_2100_RCP85)<-"Temperature"
names(LTmaxT_2050_RCP26)<-"Temperature"
names(LTminT_2050_RCP26)<-"Temperature"
names(LTmaxT_2050_RCP45)<-"Temperature"
names(LTminT_2050_RCP45)<-"Temperature"
names(LTmaxT_2050_RCP85)<-"Temperature"
names(LTminT_2050_RCP85)<-"Temperature"


#--------------------------------------------------
#-------Prepare physiological files ---------------
#--------------------------------------------------
FvFm <- FvFm %>%
  rename(Value = FvFm,
         Temperature = Temperature_C) %>%
  mutate(Trait = "FvFm")%>%
  select(Value, Trait, Temperature)

Growth <- Growth %>%
  rename(Value = Growth_rate,
         Temperature = Temperature_C) %>%
  mutate(Trait = "Growth")%>%
  select(Value, Trait, Temperature)

Germination <- Germination %>%
  rename(Value = Germination_rate_percent,
         Temperature = Temperature_C) %>%
  mutate(Trait = "Germination")%>%
  select(Value, Trait, Temperature)

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


#--------------------------------------------------
# Create fits using the best model
#--------------------------------------------------
# Prepare Hybrid_rasters folder
hybrid_folder <- file.path("data", "Input", "Hybrid_rasters")
if (!dir.exists(hybrid_folder)) dir.create(hybrid_folder)

periods <- c("present","2050","2100")
scenarios <- c("RCP26", "RCP45", "RCP85")

for (period in periods) {
  
  if (period == "present") {
    # No scenarios for present
    LTmax_r <- LTmaxT
    LTmin_r <- LTminT
    
    hybrid_raster_folder <- file.path(hybrid_folder, paste0("Layers_", period))
    if (!dir.exists(hybrid_raster_folder)) dir.create(hybrid_raster_folder)
    
    scenario_loop <- "present"  #for filename consistency
  } else {
    scenario_loop <- scenarios
  }
  
  for (scenario in scenario_loop) {
    
    if (period != "present") {
      LTmax_r <- get(paste0("LTmaxT_", period, "_", scenario))
      LTmin_r <- get(paste0("LTminT_", period, "_", scenario))
      
      hybrid_raster_folder <- file.path(hybrid_folder, paste0("Layers_future", period, "_", scenario))
      if (!dir.exists(hybrid_raster_folder)) dir.create(hybrid_raster_folder)
    }
    
    # Predict Growth, set negative values to 0
    Growth_LTmax <- predict(LTmax_r, growth_fit)
    Growth_LTmax[Growth_LTmax < 0] <- 0
    Growth_LTmin <- predict(LTmin_r, growth_fit)
    Growth_LTmin[Growth_LTmin < 0] <- 0
    
    # Predict FvFm
    FvFm_LTmax <- predict(LTmax_r, pam_fit)
    FvFm_LTmin <- predict(LTmin_r, pam_fit)
    
    # Predict Germination
    Germination_LTmax <- predict(LTmax_r, germination_fit)
    Germination_LTmin <- predict(LTmin_r, germination_fit)
    
    # Export layers
    extension <- ifelse(period == "present", period, paste0(period, "_", scenario))
    
    writeRaster(Growth_LTmax, 
                filename=file.path(hybrid_raster_folder, paste0("Growth_LTmax_", extension, ".tif")),
                overwrite=TRUE)
    writeRaster(Growth_LTmin, 
                filename=file.path(hybrid_raster_folder, paste0("Growth_LTmin_", extension, ".tif")),
                overwrite=TRUE)
    
    writeRaster(FvFm_LTmax, 
                filename=file.path(hybrid_raster_folder, paste0("FvFm_LTmax_", extension, ".tif")),
                overwrite=TRUE)
    writeRaster(FvFm_LTmin, 
                filename=file.path(hybrid_raster_folder, paste0("FvFm_LTmin_", extension, ".tif")),
                overwrite=TRUE)
    
    writeRaster(Germination_LTmax, 
                filename=file.path(hybrid_raster_folder, paste0("Germination_LTmax_", extension, ".tif")),
                overwrite=TRUE)
    writeRaster(Germination_LTmin, 
                filename=file.path(hybrid_raster_folder, paste0("Germination_LTmin_", extension, ".tif")),
                overwrite=TRUE)
  }
}
