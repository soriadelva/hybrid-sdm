#--------------------------------------------------
#--------------- Start clean  ----------------
#--------------------------------------------------
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)


#--------------------------------------------------
#--------------- Load packages ----------------
#--------------------------------------------------
library(dplyr)
library(ggplot2)

#---------------------------------------------
#-------Define folders and model types--------
#---------------------------------------------
varimp_folder<-file.path("results", "bootstrap_resampling", "variable_importance")
Model_types<-c("Correlative", "FvFm", "Germination", "Growth")

#---------------------------------------------
#-------Load data and merge in one df---------
#---------------------------------------------
varimp_df<-data.frame()
for(Model_type in Model_types){
  
  varimp_real<-read.csv2(file.path(varimp_folder, paste0(Model_type, "_variable_importance.csv")))%>%
    mutate(Model_type = Model_type)
  
  if(nrow(varimp_df)==0){
    varimp_df<-varimp_real
  }else{
    varimp_df<-bind_rows(varimp_df, varimp_real) 
  }
}


varimp_ci <- varimp_df %>%
  dplyr::group_by(Variable, Model_type) %>%
  dplyr::summarise(
    ci_lower = quantile(Percentage, probs = 0.025, na.rm = TRUE),
    ci_upper = quantile(Percentage, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )
  
  

#--------------------------------------------
#-------------  Create bargraphs -----------------
#--------------------------------------------
#Define predict times and model types
Predict.times <- c("Present","2100_RCP26", "2100_RCP45", "2100_RCP85")
Model_types <- c("Growth", "Germination", "FvFm", "Correlative")

#Start loop
deptharea_info<-data.frame()
ci_deptharea_info<-data.frame()

for(Predict.time in Predict.times){
  for(Model_type in Model_types){
    
    #Load data
    deptharea_folder<-file.path("results", "final_brt_predictions", "depth_area_overview", Predict.time)
    deptharea_all<-read.csv(file.path(deptharea_folder, paste0("depth_area_info_Full_range_", Model_type, "_", Predict.time, ".csv") ))
    deptharea_medsea<-read.csv(file.path(deptharea_folder, paste0("depth_area_info_MedSea_", Model_type, "_", Predict.time, ".csv") ))
    deptharea<-bind_rows(deptharea_all, deptharea_medsea)
    
   
    
    #Assign names based on model type
    if(nrow(deptharea_info)==0){
      deptharea_info<-deptharea
    }else{
      deptharea_info<-bind_rows(deptharea_info, deptharea)
    }
    
   
  }
}

#Define predict times and model types
Predict.times <- c("Present","2100_RCP85")
for(Predict.time in Predict.times){
  for(Model_type in Model_types){
#Load bootstrap data
ci_deptharea_folder<-file.path("results", "bootstrap_resampling", "depth_area_overview", Predict.time)
ci_deptharea_all<-read.csv(file.path(ci_deptharea_folder, paste0("depth_area_info_all_", Model_type, "_", Predict.time, ".csv") ))

#Assign names based on model type for ci
if(nrow(ci_deptharea_info)==0){
  ci_deptharea_info<-ci_deptharea_all
}else{
  ci_deptharea_info<-bind_rows(ci_deptharea_info, ci_deptharea_all)
}
  }
}

#--------------------------------------------
#---------Prepare datasets for plotting------
#--------------------------------------------
deptharea_info<- deptharea_info %>%
  dplyr::mutate(Time = factor(Time,
                              levels = c("Present", "2100_RCP26", "2100_RCP45", "2100_RCP85"),
                              labels = c("Present",
                                         "RCP 2.6",
                                         "RCP 4.5",
                                         "RCP 8.5")),
                Study_area = factor(Study_area,
                                    levels = c("Full_range", "MedSea"),
                                    labels = c("Total study area", "Mediterranean Sea")))

ci_deptharea<- ci_deptharea_info %>%
  dplyr::group_by(Model, Study_area, Bootstrap_run) %>%
  dplyr::mutate(baseline_area = Suitable_area[Time == "Present"]) %>%
  dplyr::ungroup() %>%
  
  # Compute % change
  dplyr::mutate(area_change = 100 * (Suitable_area - baseline_area) / baseline_area) %>%
  dplyr::mutate(Time = factor(Time,
                              levels = c("Present", "2100_RCP26", "2100_RCP45", "2100_RCP85"),
                              labels = c("Present",
                                         "RCP 2.6",
                                         "RCP 4.5",
                                         "RCP 8.5")),
                Study_area = factor(Study_area,
                                    levels = c("Full_range", "MedSea"),
                                    labels = c("Total study area", "Mediterranean Sea")))%>%

  
  # Keep only future scenarios
  dplyr::filter(Time %in% c("RCP 2.6",
                            "RCP 4.5",
                            "RCP 8.5"))%>%
  dplyr::group_by(Model,Time, Study_area) %>%
  dplyr::summarise(
    ci_lower_area_change = quantile(area_change, probs = 0.025, na.rm = TRUE),
    ci_upper_area_change = quantile(area_change, probs = 0.975, na.rm = TRUE),
    ci_lower_suitable_area = quantile(Suitable_area, probs = 0.025, na.rm = TRUE),
    ci_upper_suitable_area = quantile(Suitable_area, probs = 0.975, na.rm = TRUE),
    ci_lower_occ_prob = quantile(Mean_occurrence_prob, probs = 0.025, na.rm = TRUE),
    ci_upper_occ_prob = quantile(Mean_occurrence_prob, probs = 0.975, na.rm = TRUE),
    ci_lower_mtp = quantile(MTP_threshold, probs = 0.025, na.rm = TRUE),
    ci_upper_mtp = quantile(MTP_threshold, probs = 0.975, na.rm = TRUE))


#_-------------------------------------------
#-------------Plot bargraphs-----------------
#--------------------------------------------
suitable_area<-ggplot(deptharea_info,aes(x = Time,
                                         y = Suitable_area,
                                         fill = Model)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  facet_wrap(~ Study_area, nrow = 1) +
  # geom_hline(yintercept = 0,
  #            linetype = "dashed",
  #            color = "black") +
  scale_fill_manual(values = c(
    "Correlative" = "grey70",
    "FvFm"       =  "#f0a144ff",
    "Germination" = "#48a8b3ff",
    "Growth"      = "#163f6bff"),
    labels = c(
      "Correlative" = "Correlative",
      "FvFm" = expression(F[v]/F[m]),
      "Germination" = "Germination",
      "Growth" = "Growth"
    )
  ) +
  
  labs(x = NULL,
       y = expression("Suitable area (" * km^2 * ")"),
       fill = NULL ) +
  
  theme_test() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10, hjust=0.5),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "top"
  )


relative_suitability<-ggplot(deptharea_info,aes(x = Time,
                                                y = Mean_occurrence_prob,
                                                fill = Model)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  facet_wrap(~ Study_area, nrow = 1) +
  # geom_hline(yintercept = 0,
  #            linetype = "dashed",
  #            color = "black") +
  scale_fill_manual(values = c(
    "Correlative" = "grey70",
    "FvFm"       =  "#f0a144ff",
    "Germination" = "#48a8b3ff",
    "Growth"      = "#163f6bff")) +
  
  labs(x = NULL,
       y = "Avg. habitat suitability",
       fill = NULL ) +
  
  theme_test() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10,  hjust=0.5),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "none"
  )


deptharea_change<-deptharea_info %>%
  
  # Get baseline suitable area
  dplyr::group_by(Model, Study_area) %>%
  dplyr::mutate(baseline_area = Suitable_area[Time == "Present"]) %>%
  ungroup() %>%
  
  # Compute % change
  dplyr::mutate(area_change = 100 * (Suitable_area - baseline_area) / baseline_area) %>%
  
  # Keep only future scenarios
  dplyr::filter(Time %in% c("RCP 2.6",
                            "RCP 4.5",
                            "RCP 8.5"))


area_change <- ggplot() +
  geom_col(
    data = deptharea_change,
    aes(x = Time, y = area_change, fill = Model),
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_errorbar(
    data = ci_deptharea,
    aes(
      x = Time,
      ymin = ci_lower_area_change,
      ymax = ci_upper_area_change,
      group = Model
    ),
    color = "black",
    position = position_dodge(width = 0.8),
    width = 0.2,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(breaks = c(0, -25, -50, -75, 100), limits = c(15, -100)) +
  facet_wrap(~ Study_area, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Correlative" = "grey70",
      "FvFm"       = "#f0a144ff",
      "Germination" = "#48a8b3ff",
      "Growth"      = "#163f6bff"
    ),
    labels = c(
      "Correlative" = "Correlative",
      "FvFm" = expression(F[v]/F[m]),
      "Germination" = "Germination",
      "Growth" = "Growth"
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = NULL, y = "Area change (%)", fill = NULL) +
  theme_test() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "top"
  )

area_change



#suitable_area/
area_change #1000 - 533



#----------------------------------------------------------
#---------Plot prediction uncertainties--------------------
#----------------------------------------------------------
#Define predict times and model types
Predict.times <- c("Present","2100_RCP26", "2100_RCP45", "2100_RCP85")
Model_types <- c("Growth", "Germination", "FvFm", "Correlative")

for(Predict.time in Predict.times){
  for(Model_type in Model_types){
    
  raster_folder<-file.path("results", "bootstrap_resampling", "final_suitability_maps", Predict.time)
  # List all files in folder
  raster_files_all <- list.files(raster_folder, full.names = TRUE)
  
  # Pattern to match: Final_predictions_brt_"Model_type"_"Predict.time"_<anything>.tif
  if(Predict.time=="Present"){
  pattern <- paste0("Final_prediction_brt_reclassified_", Model_type, "_", tolower(Predict.time))
  }else{
  pattern<-paste0("Final_prediction_brt_reclassified_", Model_type, "_", Predict.time)
  }
  # Filter files using grep
  rast_files_to_load <- raster_files_all[grep(pattern, raster_files_all)]
  rast_files_to_load
  
  # Load rasters
  if(length(rast_files_to_load) > 0){
    rast_stack <- terra::rast(rast_files_to_load)
  } else {
    warning(paste("No rasters found for", Model_type, "at", Predict.time))
  }
  true_pred_folder<-file.path("results", "final_brt_predictions", "final_suitability_maps", Predict.time)
  if(Predict.time=="Present"){
  ref<-terra::rast(file.path(true_pred_folder, paste0("Final_prediction_brt_reclassified_",Model_type,"_", tolower(Predict.time),".tif")))
  }else{
  ref<-terra::rast(file.path(true_pred_folder, paste0("Final_prediction_brt_reclassified_",Model_type,"_", Predict.time,".tif")))
  }
  # Count 1s and 0s across the stack
  count_ones <- terra::app(rast_stack, function(x) {
  if(all(is.na(x))) {
    return(NA)
  } else {
    sum(x == 1, na.rm = TRUE)
  }
  })

count_zeros <- terra::app(rast_stack, function(x) {
  if(all(is.na(x))) {
    return(NA)
  } else {
    sum(x == 0, na.rm = TRUE)
  }
})
  # count_ones  <- terra::app(rast_stack, function(x) sum(x == 1, na.rm = TRUE))
  # count_zeros <- terra::app(rast_stack, function(x) sum(x == 0, na.rm = TRUE))

  # Create output raster based on ref
  output <- terra::ifel(ref == 1, count_ones, count_zeros)
  
  #output<-terra::app(rast_stack, sd, na.rm = TRUE)
  
  assign(paste0("Uncertainty_",Model_type, "_",Predict.time), output)
  #terra::writeRaster(output, file.path(output_folder, paste0("uncertainty_rast_",Model_type,"_",Predict.time)))
  
  }
}

terra::plot(Uncertainty_Growth_2100_RCP45)


#--------------------------------------------------
#--------------- Load packages  ----------------
#--------------------------------------------------
library(terra)
library(sf)
library(tidyterra)
library(rnaturalearth)
library(dplyr)
library(ggpubr)
library(patchwork)
library(viridis)


#-----------------------------------------------------
#-----------Load predicted distributions--------------
#-----------------------------------------------------
# Create an empty list to store plots
certainty_plots <- list()

#Define predict times and model types
Predict.times <- c("Present","2100_RCP26", "2100_RCP45", "2100_RCP85")
Model_types <- c("Correlative", "Growth", "Germination", "FvFm")

#Start loop
for(Predict.time in Predict.times){
  for(Model_type in Model_types){
  
    
    #--------------------------------------------
    #----- Download shape of Europe -------------
    #--------------------------------------------
    sf::sf_use_s2(FALSE)
    europe<-ne_countries(scale=10)%>%
      sf::st_make_valid()
    
    europe <- sf::st_crop(europe, 
                          c(xmin = -29.17,
                            xmax = 42.167,
                            ymin = 15.75,
                            ymax = 73.58333 ))%>%
      sf::st_union()
    
    
    #--------------------------------------------
    #-------------  Create plots  ---------------
    #--------------------------------------------
    
    #--------------------
    #--Suitability_plot--
    #--------------------
    brks <- seq(0, 50, by=1)
    nb <- length(brks) - 1
    viridis_palette <- viridis(nb)
    
    raster_file<-get(paste0("Uncertainty_",Model_type, "_",Predict.time))
    
    
    template<-raster_file
    values(template) <- rep(brks, length.out = ncell(template))
    template<-mask(template, raster_file)
    #Get extent
    exten<-as.vector(terra::ext(raster_file))
    
    text_label<-switch(Model_type,
                       "Correlative"="bold(Correlative)",
                       "FvFm"="bold(F[v]/F[m])",
                       "Germination"="bold(Germination)",
                       "Growth"= "bold(Growth)")
    
    certainty<-ggplot() +
      geom_sf(data=europe, color= "#f2f2f2", fill = "#f2f2f2")+
      geom_spatraster(data=template)+
      geom_spatraster(data = raster_file) +
      scale_fill_gradientn(colors = viridis_palette, 
                           values=scales::rescale(brks),
                           breaks = c(0,10,20,30,40,50), 
                           labels = c(0,10,20,30,40,50), 
                           na.value = NA) +
      labs(fill = "Bootstrap agreement (number of models)")+
      scale_y_continuous(breaks=c(29,44,59))+
      scale_x_continuous(position = "top", breaks=c(-12,7,26)) + 
      theme_bw() +
      theme(axis.title = element_blank(),
            #plot.margin = margin(0, 0, 0, 0),
            axis.text = element_text(size = 10))+
      #theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
      coord_sf(xlim = c(exten[1], exten[2]), 
               ylim = c(exten[3], exten[4]), 
               expand=FALSE)+
      theme(aspect.ratio=0.9)+ #0.85
      annotate(geom="label",
               x= -12,
               y=70,
               label=text_label, 
               fontface="bold",
               fill="white", 
               parse=TRUE,
               hjust=0, 
               size=4)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
   
    plot_type<-certainty

      if (Model_type == "Correlative") {
        plot_type <- plot_type +
          theme(
            axis.text.x = element_text(color = "transparent"),
            axis.ticks.x = element_line(color = "transparent"),
            axis.text.y = element_text(color = "black"),
            axis.ticks.y = element_line(color = "black")
          )
        
      } else if (Model_type == "FvFm") {
        plot_type <- plot_type +
          theme(
            axis.text.y = element_text(color = "transparent", angle=90),
            axis.ticks.y = element_line(color = "transparent"),
            axis.text.x = element_text(color = "transparent"),
            axis.ticks.x = element_line(color = "transparent")
          )
      } else if (Model_type == "Germination") {
        plot_type <-plot_type +
          theme(
            axis.text.y = element_text(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.text.x = element_text(color = "black"),
            axis.ticks.x = element_line(color = "black")
          )
      } else if (Model_type == "Growth") {
        plot_type <-plot_type +
          theme(
            axis.text.y = element_text(color = "transparent", angle=90),
            axis.ticks.y = element_line(color = "transparent"),
            axis.text.x = element_text(color = "black"),
            axis.ticks.x = element_line(color = "black")
          )}
      
      plot_name<-paste0("Certainty_",Model_type, "_", Predict.time)
      certainty_plots[[plot_name]]<-plot_type
      
    }
  }


#-------------------
#certainty_graphs
#-------------------

#--------------------------Create legends------------------------

legend_horizontal <- ggpubr::get_legend(
  certainty_plots$Certainty_Correlative_Present + 
    guides(
      fill = guide_colorbar(
        title.position = "bottom",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = 10,
        barheight = 1
      )
    ) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(2, "cm")
    )
)

legend_horizontal <- ggpubr::as_ggplot(legend_horizontal)

#------------------Create plots----------------------------- 
final_plot_present<-wrap_plots(
  certainty_plots$Certainty_Correlative_Present,
  certainty_plots$Certainty_FvFm_Present,
  certainty_plots$Certainty_Germination_Present,
  certainty_plots$Certainty_Growth_Present,
  ncol = 2
)/legend_horizontal &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP26<- wrap_plots(
  certainty_plots$Certainty_Correlative_2100_RCP26,
  certainty_plots$Certainty_FvFm_2100_RCP26,
  certainty_plots$Certainty_Germination_2100_RCP26,
  certainty_plots$Certainty_Growth_2100_RCP26,
  ncol = 2
)/legend_horizontal &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")


final_plot_2100_RCP45<- wrap_plots(
  certainty_plots$Certainty_Correlative_2100_RCP45,
  certainty_plots$Certainty_FvFm_2100_RCP45,
  certainty_plots$Certainty_Germination_2100_RCP45,
  certainty_plots$Certainty_Growth_2100_RCP45,
  ncol = 2
)/legend_horizontal &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP85 <- wrap_plots(
  certainty_plots$Certainty_Correlative_2100_RCP85,
  certainty_plots$Certainty_FvFm_2100_RCP85,
  certainty_plots$Certainty_Germination_2100_RCP85,
  certainty_plots$Certainty_Growth_2100_RCP85,
  ncol = 2
)/legend_horizontal &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")


plot_present <- final_plot_present + plot_layout(heights = c( 1, 0.1))
plot_2100_RCP26 <- final_plot_2100_RCP26 + plot_layout(heights = c( 1, 0.1))
plot_2100_RCP45 <- final_plot_2100_RCP45 + plot_layout(heights = c( 1, 0.1))
plot_2100_RCP85 <- final_plot_2100_RCP85 + plot_layout(heights = c( 1, 0.1))

plot_present #1000 width, 863 height
plot_2100_RCP45
plot_2100_RCP85

plot_path<-file.path("figures")
ggsave(filename="Certainty_present_Appendix.jpeg",plot=plot_present, width = 20, height = 20, units = "cm", path=plot_path)
ggsave(filename="Certainty_2100_RCP26_Appendix.jpeg",plot=plot_2100_RCP26, width = 20, height = 20, units = "cm", path=plot_path)
ggsave(filename="Certainty_2100_RCP45_Appendix.jpeg",plot=plot_2100_RCP45, width = 20, height = 20, units = "cm", path=plot_path)
ggsave(filename="Certainty_2100_RCP85_Appendix.jpeg",plot=plot_2100_RCP85, width = 20, height = 20, units = "cm", path=plot_path)

