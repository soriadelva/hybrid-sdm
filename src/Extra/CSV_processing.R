#--------------------------------------------------
#--------------- Clean up workspace----------------
#--------------------------------------------------
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)


#--------------------------------------------------
#--------------- Load packages ----------------
#--------------------------------------------------
library(dplyr)


#--------------------------------------------------
#--------------- Rename columns ----------------
#--------------------------------------------------
Dictyota_germination <- Dictyota_germination %>%
  dplyr::rename(Population_id = Population,
                Strain_id = Replicate,
                Temperature_c = Temperature,
                Germinated_apical_rhizoidal = Germinated.with.rhizoid,
                Germinated_no_rhizoidal= Germinated.without.rhizoid,
                Not_germinated  = Not.germinated,
                Germination_rate_percent = Germination.rate....,
                Days_post_collection = Days.in.lab.when.isolated) %>%
  mutate(strain_id = paste0("Goes_", strain_id)) %>%
  relocate(days_post_collection, .after = temperature_c)

Dictyota_growth <- Dictyota_growth %>%
  dplyr::rename(Population_id = Population,
                Strain_id = Replicate,
                Temperature_C = Temperature,
                Surface_t0 = Initial_surface,
                Surface_t1 = End_surface,
                Growth_rate = Growth,
                Growth_old = Growth_old) %>%
  dplyr::select(-Growth_old)

Dictyota_FvFm <- Dictyota_FvFm %>%
  dplyr::mutate(Population_id = "Goes")%>%
  dplyr::rename(Strain_id = Strain,
                Temperature_C = Temperature,
                Corrected_Fm = corrected_Fm,
                FvFm = YII) %>%
  dplyr::select(-c(label, standaard_F0, standaard_Fm)) %>%
  dplyr::relocate(Population_id, .before = Strain_id)

Dictyota_full_occurrences <- Dictyota_full_occurrences %>%
  dplyr::mutate(species = "Dictyota dichotoma")%>%
  dplyr::rename(Species = species,
                Latitude = Lat,
                Longitude = Lon)
  
  
#--------------------------------------------------
#--------------- Export CSVs----------------
#--------------------------------------------------
write.csv2(Dictyota_FvFm,"C:/Users/soria_delva/Desktop/Dictyota_Zenodo/Dictyota_FvFm.csv", row.names = FALSE)
write.csv2(Dictyota_germination, "C:/Users/soria_delva/Desktop/Dictyota_Zenodo/Dictyota_germination.csv", row.names = FALSE)
write.csv2(Dictyota_growth, "C:/Users/soria_delva/Desktop/Dictyota_Zenodo/Dictyota_growth.csv", row.names = FALSE)
write.csv2(Dictyota_full_occurrences, "C:/Users/soria_delva/Desktop/Dictyota_Zenodo/Dictyota_full_occurrences.csv", row.names = FALSE)


