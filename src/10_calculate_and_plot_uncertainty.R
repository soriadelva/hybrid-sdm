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
area_change #600-350

