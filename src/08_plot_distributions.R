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
suitability_plots <- list()

#Define predict times and model types
Predict.times <- c("Present","2100_RCP26", "2100_RCP45", "2100_RCP85")
Model_types <- c("Correlative", "Growth", "Germination", "FvFm")

#Start loop
for(Predict.time in Predict.times){
  for(Model_type in Model_types){
    
    extens<-ifelse(Predict.time=="Present", paste0("_",Predict.time,".tif"), paste0(Predict.time,".tif") )
    
    #Load data
    results_folder<-file.path("results", "final_brt_predictions", "final_suitability_maps", Predict.time)
    suit_prediction<-terra::rast(file.path(results_folder,paste0("Final_prediction_brt_",Model_type, extens)))
    reclas_prediction<-terra::rast( file.path(results_folder, paste0("Final_prediction_brt_reclassified_",Model_type,extens)))
    
    #Assign names based on model type
    assign(paste0(Model_type,"_suitability_", Predict.time), suit_prediction)
    assign(paste0(Model_type,"_reclas_suitability_", Predict.time), reclas_prediction)
    
    #Clean up
    rm(suit_prediction, reclas_prediction)
    
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
      st_union()
    
    
    #--------------------------------------------
    #-------------  Create plots  ---------------
    #--------------------------------------------
    
    #--------------------
    #--Suitability_plot--
    #--------------------
    brks <- seq(0, 1, by=0.2)
    nb <- length(brks) - 1
    viridis_palette <- viridis(nb)
    
    
    raster_file<-get(paste0(Model_type,"_suitability_", Predict.time))
    
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
    
    suitability<-ggplot() +
      geom_sf(data=europe, color= "#f2f2f2", fill = "#f2f2f2")+
      geom_spatraster(data=template)+
      geom_spatraster(data = raster_file) +
      scale_fill_gradientn(colors = viridis_palette, 
                           breaks = brks, 
                           labels = brks, 
                           na.value = NA) +
      labs(fill = "Suitability")+
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
    
    
    #--------------------
    #--Reclassified_plot--
    #--------------------
    #Get raster file
    raster_file<-get(paste0(Model_type,"_reclas_suitability_", Predict.time))
    raster_file <- as.factor( raster_file*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
    levels( raster_file) <- data.frame(ID = c(0, 1),
                                       class = c("Absent", "Present"))
    #Get extent
    exten<-as.vector(terra::ext(raster_file))
    
    #Create label
    text_label<-switch(Model_type,
                       "Correlative"="bold(Correlative)",
                       "FvFm"="bold(F[v]/F[m])",
                       "Germination"="bold(Germination)",
                       "Growth"= "bold(Growth)")
    
    reclas_pred <- ggplot() + 
      geom_sf(data=europe, color= "#e8e7e7ff", fill = "#e8e7e7ff")+
      geom_spatraster(data = raster_file) +
      scale_fill_manual(values = c("Absent" ="#fff9b7" , "Present" = "#085099"),#"##fff7aa"
                        na.value = "transparent",
                        na.translate=FALSE)+
      scale_y_continuous(breaks=c(29,44,59))+
      scale_x_continuous(position = "top", breaks=c(-12,7,26)) + 
      theme_bw() +
      theme(axis.title = element_blank(),
            legend.title= element_blank(),
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
    
    reclas_pred
    
    for(plot_typ in c("suitability", "reclas_pred")){
      if (plot_typ == "suitability") {
        plot_type <- suitability
      } else {
        plot_type <- reclas_pred
      }
      
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
      
      plot_name<-ifelse(plot_typ=="suitability",paste0(Model_type, "_", Predict.time), paste0(Model_type, "_", Predict.time, "_reclassified") )
      suitability_plots[[plot_name]]<-plot_type
      
    }
  }
}


#--------------------------------------------
#-------------  Plot graphs -----------------
#--------------------------------------------
#-------------------
#Suitability_graphs
#-------------------

#--------------------------Create legends------------------------

legend_horizontal <- ggpubr::get_legend(
  suitability_plots$Correlative_Present + 
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
  suitability_plots$Correlative_Present,
  suitability_plots$FvFm_Present,
  suitability_plots$Germination_Present,
  suitability_plots$Growth_Present,
  ncol = 2
)/legend_horizontal &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP26<- wrap_plots(
  suitability_plots$Correlative_2100_RCP26,
  suitability_plots$FvFm_2100_RCP26,
  suitability_plots$Germination_2100_RCP26,
  suitability_plots$Growth_2100_RCP26,
  ncol = 2
)/legend_horizontal &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP45<- wrap_plots(
  suitability_plots$Correlative_2100_RCP45,
  suitability_plots$FvFm_2100_RCP45,
  suitability_plots$Germination_2100_RCP45,
  suitability_plots$Growth_2100_RCP45,
  ncol = 2
)/legend_horizontal &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP85 <- wrap_plots(
  suitability_plots$Correlative_2100_RCP85,
  suitability_plots$FvFm_2100_RCP85,
  suitability_plots$Germination_2100_RCP85,
  suitability_plots$Growth_2100_RCP85,
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


#-------------------
#Reclassified_graphs
#-------------------

#--------------------------Create legends------------------------

legend_horizontal_reclas <- ggpubr::get_legend(
  suitability_plots$Correlative_Present_reclassified +
    guides(
      fill = guide_legend(
        title = NULL,          # position of the title
            # horizontal justification
        direction = "horizontal",   # make legend horizontal
        keywidth = unit(0.5 ,"cm"),   # width of each key
        keyheight = unit(0.5, "cm") # height of each key
      )
    ) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12)
    )
)

legend_horizontal_reclas <- ggpubr::as_ggplot(legend_horizontal_reclas)

#------------------Create plots----------------------------- 
final_plot_present_reclassified<-wrap_plots(
  suitability_plots$Correlative_Present_reclassified,
  suitability_plots$FvFm_Present_reclassified,
  suitability_plots$Germination_Present_reclassified,
  suitability_plots$Growth_Present_reclassified,
  ncol = 2
)/legend_horizontal_reclas&
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP26_reclassified<- wrap_plots(
  suitability_plots$Correlative_2100_RCP26_reclassified,
  suitability_plots$FvFm_2100_RCP26_reclassified,
  suitability_plots$Germination_2100_RCP26_reclassified,
  suitability_plots$Growth_2100_RCP26_reclassified,
  ncol = 2
)/legend_horizontal_reclas &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP45_reclassified<- wrap_plots(
  suitability_plots$Correlative_2100_RCP45_reclassified,
  suitability_plots$FvFm_2100_RCP45_reclassified,
  suitability_plots$Germination_2100_RCP45_reclassified,
  suitability_plots$Growth_2100_RCP45_reclassified,
  ncol = 2
)/legend_horizontal_reclas &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_2100_RCP85_reclassified <- wrap_plots(
  suitability_plots$Correlative_2100_RCP85_reclassified,
  suitability_plots$FvFm_2100_RCP85_reclassified,
  suitability_plots$Germination_2100_RCP85_reclassified,
  suitability_plots$Growth_2100_RCP85_reclassified,
  ncol = 2
)/legend_horizontal_reclas &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")


plot_present_reclassified <- final_plot_present_reclassified + plot_layout(heights = c( 1, 0.1))
plot_2100_RCP26_reclassified <- final_plot_2100_RCP26_reclassified + plot_layout(heights = c( 1, 0.1))
plot_2100_RCP45_reclassified <- final_plot_2100_RCP45_reclassified + plot_layout(heights = c( 1, 0.1))
plot_2100_RCP85_reclassified <- final_plot_2100_RCP85_reclassified + plot_layout(heights = c( 1, 0.1))

plot_present_reclassified #1000 width, 863 height
plot_2100_RCP45_reclassified
plot_2100_RCP85_reclassified


#---------------------------------------------------------
#-----Create combined prediction plots for hybrid models--
#---------------------------------------------------------
# Create an empty list to store plots
combined_plots <- list()

#Define predict times and model types
Predict.times <- c("Present","2100_RCP26", "2100_RCP45", "2100_RCP85")
Model_types <- c("Growth", "Germination", "FvFm", "Correlative")

#Start loop
for(Predict.time in Predict.times){
    for(Model_type in Model_types){
    extens<-ifelse(Predict.time=="Present", paste0("_",Predict.time,".tif"), paste0(Predict.time,".tif") )
    
    #Load data
    results_folder<-file.path("results", "final_brt_predictions", "final_suitability_maps", Predict.time)

    reclas_prediction<-terra::rast( file.path(results_folder, paste0("Final_prediction_brt_reclassified_",Model_type,extens)))
    
    #Assign names based on model type
    assign(paste0(Model_type,"_reclas_suitability_", Predict.time), reclas_prediction)
    
    #Clean up
    rm(suit_prediction, reclas_prediction)
    }
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
      st_union()
    
    sf::sf_use_s2(TRUE)
    #--------------------------------------------
    #-------------  Create plots  ---------------
    #--------------------------------------------
    
    #--------------------
    #--COmbined_plot--
    #--------------------
    #Load rasters
    correlative_raster <- get(paste0("Correlative_reclas_suitability_", Predict.time))
    growth_raster<-get(paste0("Growth_reclas_suitability_", Predict.time))
    FvFm_raster<-get(paste0("FvFm_reclas_suitability_", Predict.time))
    Germination_raster<-get(paste0("Germination_reclas_suitability_", Predict.time))
    
    # Process hybrid rasters
    rasters_stack <- c(growth_raster, FvFm_raster, Germination_raster, correlative_raster)
    
    # Initialize an empty raster for the combined categories
    combined_raster <- growth_raster
    combined_raster[] <- NA  # start with NA
    
    # Assign categories based on rules
    # 1 = Consensus presence
    combined_raster[growth_raster == 1 & FvFm_raster == 1 & Germination_raster == 1 & correlative_raster == 1] <- 1
    
    # 2 = Consensus absence
    combined_raster[growth_raster == 0 & FvFm_raster == 0 & Germination_raster == 0 & correlative_raster == 0] <- 2
    
    # 3 = Hybrid presence: only correlative is 0, others are 1
    combined_raster[growth_raster == 1 & FvFm_raster == 1 & Germination_raster == 1 & correlative_raster == 0] <- 3
    
    # 4 = Hybrid absence: only correlative is 1, others are 0
    combined_raster[growth_raster == 0 & FvFm_raster == 0 & Germination_raster == 0 & correlative_raster == 1] <- 4
    
    # 5 = Ambiguous: everything else
    original_na <- is.na(growth_raster) | is.na(FvFm_raster) | is.na(Germination_raster) | is.na(correlative_raster) #original NAs
    
    # Assign Ambiguous only to remaining cells that are NOT original NA
    combined_raster[is.na(combined_raster) & !original_na] <- 5
    
    # Optional: assign labels
    categories <- c("Consensus presence", "Consensus absence", "Hybrid presence",
                    "Hybrid absence", "Ambiguous")
    rat <- data.frame(ID = 1:5, Category = categories)
    levels(combined_raster) <- rat
    #Get extent
    exten<-as.vector(terra::ext(raster_file))
    
    #Create label
    text_label<-switch(Predict.time,
                       "Present"="Present",
                       "2100_RCP26"="2100 - RCP 2.6",
                       "2100_RCP45"="2100 - RCP 4.5",
                       "2100_RCP85"="2100 - RCP 8.5")
  
    
    comb_pred <- ggplot() + 
      geom_sf(data=europe, color= "#e8e7e7ff", fill = "#e8e7e7ff")+
      geom_spatraster(data = combined_raster, interpolate=F,  maxcell = Inf) +
      scale_fill_manual(values = c("Consensus absence" ="#fff9b7" ,
                                   "Consensus presence" = "#085099",
                                   "Ambiguous"= "#b21908",
                                   "Hybrid presence"="#65cfd8ff",
                                   "Hybrid absence"="#f68e3f"),#"##fff7aa"
                        na.value = "transparent",
                        na.translate=FALSE)+
      scale_y_continuous(breaks=c(29,44,59))+
      scale_x_continuous(position = "top", breaks=c(-12,7,26)) + 
      annotate(geom="label",
               x= ifelse(Predict.time=="Present",-12,-16),
               y=70,
               label=text_label,
               fontface="bold",
               fill="white",
               hjust=0,
               size=4)+
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 10))+
      coord_sf(xlim = c(exten[1], exten[2]), 
               ylim = c(exten[3], exten[4]), 
               expand=FALSE)+
      theme(aspect.ratio=0.9)+ #0.85
  
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
         
        comb_pred

        if(Predict.time=="Present"){
          comb_pred <- comb_pred +
            theme(
              axis.text.x = element_text(color = "transparent"),
              axis.ticks.x = element_line(color = "transparent"),
              axis.text.y = element_text(color = "black"),
              axis.ticks.y = element_line(color = "black")
            )
        }
        
        if(Predict.time=="2100_RCP26"){
          comb_pred <- comb_pred +
            theme(
              axis.text.x = element_text(color = "transparent"),
              axis.ticks.x = element_line(color = "transparent"),
              axis.text.y = element_text(color = "transparent", angle=90),
              axis.ticks.y = element_line(color = "transparent")
            )
        }
        
        if(Predict.time=="2100_RCP45"){
          comb_pred <- comb_pred +
            theme(
              axis.text.x = element_text(color = "black"),
              axis.ticks.x = element_line(color = "black"),
              axis.text.y = element_text(color = "black"),
              axis.ticks.y = element_line(color = "black")
            )
        }
        
        if(Predict.time=="2100_RCP85"){
          comb_pred <- comb_pred +
            theme(
              axis.text.x = element_text(color = "black"),
              axis.ticks.x = element_line(color = "black"),
              axis.text.y = element_text(color = "transparent", angle=90),
              axis.ticks.y = element_line(color = "transparent")
            )
        }
        
      plot_name<-paste0("Combined_reclassified_", Predict.time)
      combined_plots [[plot_name]]<-comb_pred
  
      
    }



#-------------------------------
#Reclassified combined graphs
#------------------------------

#--------------------------Create legends------------------------
legend_horizontal_comb <- ggpubr::get_legend(
  combined_plots$Combined_reclassified_Present +
    guides(
      fill = guide_legend(
        title = NULL,          # position of the title
        # horizontal justification
        direction = "horizontal",   # make legend horizontal
        keywidth = unit(0.4 ,"cm"),   # width of each key
        keyheight = unit(0.4, "cm") # height of each key
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.justification = "right")
)

legend_horizontal_comb <- ggpubr::as_ggplot(legend_horizontal_comb)


#--------------------------Create plots------------------------
final_plot_combined<-legend_horizontal_comb/wrap_plots(
  combined_plots$Combined_reclassified_Present,
  combined_plots$Combined_reclassified_2100_RCP26,
  combined_plots$Combined_reclassified_2100_RCP45,
  combined_plots$Combined_reclassified_2100_RCP85,
  ncol = 2
)&
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

final_plot_combined<- final_plot_combined+ plot_layout(heights = c( 0.1, 1))
final_plot_combined


#--------------------------------------------
#-------------  Save graphs -----------------
#--------------------------------------------

 plot_path<-file.path("figures")
 ggsave(filename="Suitability_present_Appendix.jpeg",plot=plot_present, width = 20, height = 20, units = "cm", path=plot_path)
 ggsave(filename="Suitability_2100_RCP26_Appendix.jpeg",plot=plot_2100_RCP26, width = 20, height = 20, units = "cm", path=plot_path)
 ggsave(filename="Suitability_2100_RCP45_Appendix.jpeg",plot=plot_2100_RCP45, width = 20, height = 20, units = "cm", path=plot_path)
 ggsave(filename="Suitability_2100_RCP85_Appendix.jpeg",plot=plot_2100_RCP85, width = 20, height = 20, units = "cm", path=plot_path)
 
 
 ggsave(filename="Reclassified_suitability_present_Appendix.jpeg",plot=plot_present_reclassified, width = 20, height = 20, units = "cm", path=plot_path)
 ggsave(filename="Reclassified_suitability_2100_RCP26_Appendix.jpeg",plot=plot_2100_RCP26_reclassified, width = 20, height = 20, units = "cm", path=plot_path)
 ggsave(filename="Reclassified_suitability_2100_RCP45_Appendix.jpeg",plot=plot_2100_RCP45_reclassified, width = 20, height = 20, units = "cm", path=plot_path)
 ggsave(filename="Reclassified_suitability_2100_RCP85_Appendix.jpeg",plot=plot_2100_RCP85_reclassified, width = 20, height = 20, units = "cm", path=plot_path)

 ggsave(filename="Fina_figure_5.jpeg",plot=final_plot_combined, width = 21, height = 21, units = "cm", path=plot_path,   dpi = 400)
# sf::sf_use_s2(TRUE)
